#!/usr/bin/python3

import argparse
import sys
from astropy.io import fits
from astropy.stats import mad_std
import numpy as np
import pyregion
from astropy.wcs import WCS

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Calculate noise statistics from RM cube using region file')
    parser.add_argument('input_cube', help='Input RM cube FITS file')
    parser.add_argument('region_file', help='PyRegion file defining noise regions')
    return parser.parse_args()

def load_data(cube_file):
    """Load RM cube data and WCS information."""
    try:
        with fits.open(cube_file) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            wcs = WCS(header)
        return data, wcs
    except Exception as e:
        print(f"Error loading cube file {cube_file}: {e}")
        sys.exit(1)

def load_regions(region_file):
    """Load regions from pyregion file."""
    try:
        regions = pyregion.open(region_file)
        return regions
    except Exception as e:
        print(f"Error loading region file {region_file}: {e}")
        sys.exit(1)

def extract_noise_data(data, wcs, regions):
    """Extract noise data from specified regions."""
    noise_arrays = []
    
    # Get the 2D spatial dimensions for masking
    spatial_shape = data.shape[-2:]  # (ny, nx)
    
    # Create a dummy HDU with the spatial shape and WCS
    dummy_data = np.zeros(spatial_shape)
    dummy_hdu = fits.PrimaryHDU(data=dummy_data, header=wcs.to_header())
        
    # Get mask from all regions
    mask = regions.get_mask(dummy_hdu)
        
    # Extract low and high RM data for all regions combined
    # Assuming RM axis is axis=1 (modify if different)
    low_rm_data = data[0, :250, mask]
    high_rm_data = data[0, -250:, mask]
        
    noise_arrays.extend([low_rm_data.flatten(), high_rm_data.flatten()])
    
    if len(noise_arrays) == 0:
        print("Error: No valid regions found or processed")
        sys.exit(1)
        
    return np.concatenate(noise_arrays)

def calculate_noise_stats(noise_data):
    """Calculate noise statistics."""
    # Remove NaN values
    clean_data = noise_data[~np.isnan(noise_data)]
    
    if len(clean_data) == 0:
        print("Error: No valid data found in noise regions")
        sys.exit(1)
    
    noise_median = np.median(clean_data)
    noise_std = mad_std(clean_data, ignore_nan=True)
    
    return noise_median, noise_std

def main():
    """Main function."""
    args = parse_arguments()
    
    if not "FDF_tot_dirty.fits" in args.input_cube:
        print("WARNING: Input cube filename should end with: 'FDF_tot_dirty.fits'")

    # Load data
    print(f"Loading RM cube: {args.input_cube}")
    data, wcs = load_data(args.input_cube)
    
    print(f"Loading regions: {args.region_file}")
    regions = load_regions(args.region_file)
    
    # Extract noise data from regions
    print(f"Extracting noise data from {len(regions)} regions...")
    noise_data = extract_noise_data(data, wcs, regions)
    
    # Calculate statistics
    noise_median, noise_std = calculate_noise_stats(noise_data)
    
    # Print results
    print(f"Number of noise samples: {len(noise_data)}")
    print(f"RM cube noise=({noise_median:.9f}±{noise_std:.9f})")
    print(f"In rmclean3d one can use 6x sigma = {6 * noise_std:.9f} in the -c option")

if __name__ == "__main__":
    main()