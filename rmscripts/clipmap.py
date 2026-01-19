#!/usr/bin/env python3

import argparse
import numpy as np
from astropy.io import fits
import sys
import os

#!/usr/bin/env python3

def main():
    parser = argparse.ArgumentParser(description='Blank pixels in a map based on P-map threshold')
    parser.add_argument('-pmap_file', '-p', help='Polarisation intensity map FITS file')
    parser.add_argument('-noise', '-n', type=float, help='Noise threshold value')
    parser.add_argument('-sigma', '-s', type=int, help='Sigma multiplier')
    parser.add_argument('map_to_mask', nargs='+', help='FITS file(s) to mask')
    
    args = parser.parse_args()
    
    try:
        # Load P-map
        with fits.open(args.pmap_file) as pmap_hdul:
            pmap_data = pmap_hdul[0].data
        
        # Calculate threshold once
        threshold = args.noise * args.sigma
        mask = pmap_data < threshold
        
        # Process each map
        for map_file in args.map_to_mask:
            # Load map to mask
            with fits.open(map_file) as map_hdul:
                map_data = map_hdul[0].data.copy()
                map_header = map_hdul[0].header
            
            # Blank pixels below threshold to NaN
            map_data[mask] = np.nan
            
            # Create output filename
            base_name = os.path.splitext(map_file)[0]
            output_file = f"{base_name}-blanked.fits"
            
            # Save the result
            fits.writeto(output_file, map_data, map_header, overwrite=True)
            print(f"Blanked map saved as: {output_file}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()