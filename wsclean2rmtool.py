#!/usr/bin/env python3

import argparse
import re, glob, sys
from astropy.io import fits
import numpy as np

def create_cube(inputfiles, outputfile):
    """
    Create a FITS cube from a list of input FITS files.
    Stacks 4D cubes along the frequency axis (axis 3 in FITS).
    Preserves the degenerate Stokes axis (axis 4 in FITS).
    """
    data_list = []; freqs = []
    header = None
    
    for file in sorted(inputfiles):
        with fits.open(file) as hdul:
            data = hdul[0].data
            header = hdul[0].header.copy()
            # Assume 4D data: (1, 1, ny, nx) -> extract frequency plane
            if data.ndim == 4:
                data = data[0, 0, :, :]  # Extract spatial plane (ny, nx)
            data_list.append(data)
            freqs.append(header.get('CRVAL3', 0.0))  # Frequency value from header
    
    # Stack along frequency axis
    cube_data = np.stack(data_list, axis=0)
    
    # Reshape to 4D: (1, 1, ny, nx) -> (1, nfreq, ny, nx)
    # Then transpose to FITS order: (1, 1, ny, nx, nfreq) requires shape (1, 1, ny, nx) repeated
    # FITS axis order is (naxis4, naxis3, naxis2, naxis1) = (freq, stokes, dec, ra)
    # NumPy order is reversed: (freq, stokes, ny, nx)
    cube_data = cube_data[np.newaxis, :, :, :]  # (1, nfreq, ny, nx)
    
    # Update header for frequency axis (NAXIS4) and preserve Stokes axis (NAXIS3)
    if header:
        nfreq = cube_data.shape[0]
        header['NAXIS3'] = nfreq
        header['NAXIS4'] = 1  # Degenerate Stokes axis
    
    hdu = fits.PrimaryHDU(cube_data, header=header)
    hdu.writeto(outputfile, overwrite=True)
    print(f"Created cube: {outputfile} with shape {cube_data.shape}")

    return freqs

def main():
    parser = argparse.ArgumentParser(
        description="Organize WSClean Stokes Q and U images, optionally normalize by Stokes I"
    )
    parser.add_argument("name_qu", type=str, help="Base name for Q/U files (e.g., 'img/test')")
    parser.add_argument("--name-i", "-i", help="Name for Stokes I files (e.g. 'img/testI')")
    parser.add_argument("--output", "-o", default="", help="Base name for output files or use StokesQ.fits and StokesU.fits")
    parser.add_argument("--freq-file", "-f", default="", help="Output file for frequency list")
    parser.add_argument("--convolve", "-c", action="store_true", help="Convolve Q/U images to minimum resolution before stacking")

    args = parser.parse_args()
    
    if args.output != "":
        output_q = args.output + "_StokesQ.fits"
        output_u = args.output + "_StokesU.fits"
    else:     
        output_q = "StokesQ.fits"
        output_u = "StokesU.fits"

    # Collect and move Q files
    pattern_q = re.compile(f"{re.escape(args.name_qu)}-(\d+)-Q-image\.fits")
    pattern_u = re.compile(f"{re.escape(args.name_qu)}-(\d+)-U-image\.fits")
    files_q = []; files_u = []

    for filename in glob.glob(f"{args.name_qu}*.fits"):
        match = pattern_q.match(filename)
        if match:
            files_q.append(filename)
        match = pattern_u.match(filename)
        if match:
            files_u.append(filename)
    print(f"Found {len(files_q)} Stokes Q files: {files_q}")
    if args.convolve:
        print("Convolution to minimum resolution is not implemented in this version.")
        sys.exit(1)
    freqs = create_cube(files_q, output_q)
    print(f"Found {len(files_u)} Stokes U files: {files_u}")
    freqs = create_cube(files_u, output_u)
    
    if args.freq_file:
        # Extract frequencies from one of the Q files
        np.savetxt(args.freq_file, freqs, delimiter="")

    print("\nDone!")


if __name__ == "__main__":
    main()