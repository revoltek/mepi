#!/usr/bin/env python3

import sys
import numpy as np
from astropy.io import fits
import argparse

def check_beam_equality(header1, header2, tolerance=1e-6):
    """Check if two FITS headers have equal beam parameters."""
    beam_keys = ['BMAJ', 'BMIN', 'BPA']
    
    for key in beam_keys:
        if key not in header1 or key not in header2:
            print(f"Warning: Beam parameter {key} not found in one or both headers")
            return False
        
        val1 = header1[key]
        val2 = header2[key]
        
        if abs(val1 - val2) > tolerance:
            print(f"Beam parameter {key} mismatch: {val1} vs {val2}")
            return False
    
    return True

def main():
    parser = argparse.ArgumentParser(description='Calculate fractional polarization P/I from radio maps')
    parser.add_argument('i_map', help='Input I (total intensity) FITS file')
    parser.add_argument('p_map', help='Input P (polarized intensity) FITS file')
    parser.add_argument('-o', '--output', default='fractpol.fits', help='Output fractional polarization map')
    
    args = parser.parse_args()
    
    # Load FITS files
    try:
        with fits.open(args.i_map) as hdul_i:
            i_data = hdul_i[0].data
            i_header = hdul_i[0].header
        
        with fits.open(args.p_map) as hdul_p:
            p_data = hdul_p[0].data
            p_header = hdul_p[0].header
    
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Check beam equality
    if not check_beam_equality(i_header, p_header):
        print("Error: Beam shapes are not equal between I and P maps")
        sys.exit(1)
    
    # Check data shapes match
    if i_data.shape != p_data.shape:
        print(f"Error: Data shapes don't match - I: {i_data.shape}, P: {p_data.shape}")
        sys.exit(1)
    
    # Calculate fractional polarization P/I
    # Avoid division by zero
    mask = i_data != 0
    fract_pol = np.zeros_like(i_data)
    fract_pol[mask] = p_data[mask] / i_data[mask]
    
    # Set invalid values to NaN
    fract_pol[~mask] = np.nan
    
    # Create output header based on I map
    output_header = i_header.copy()
    output_header['BUNIT'] = 'Fractional'
    output_header['HISTORY'] = f'Fractional polarization P/I created from {args.p_map} / {args.i_map}'
    
    # Save result
    fits.writeto(args.output, fract_pol, header=output_header, overwrite=True)
    print(f"Fractional polarization map saved to {args.output}")

if __name__ == "__main__":
    main()