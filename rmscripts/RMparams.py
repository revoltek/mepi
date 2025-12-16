#!/usr/bin/python3

import argparse
import sys
import numpy as np
from casacore.tables import table

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Derive RM synthesis parameters from measurement set')
    parser.add_argument('ms_path', help='Path to measurement set')
    return parser.parse_args()

def read_ms_spectral_info(ms_path):
    """Read frequency information from measurement set."""
    try:
        # Open the spectral window table
        with table(f'{ms_path}/SPECTRAL_WINDOW') as spw_table:
            chan_freq = spw_table.getcol('CHAN_FREQ')  # Shape: (nspws, nchans)
            
        # Flatten frequencies if multiple SPWs
        if chan_freq.ndim > 1:
            frequencies = chan_freq.flatten()
        else:
            frequencies = chan_freq
        
        # Remove any duplicate frequencies and sort
        frequencies = np.unique(frequencies)
        frequencies = np.sort(frequencies)
        
        return frequencies
        
    except Exception as e:
        print(f"Error reading measurement set {ms_path}: {e}")
        sys.exit(1)

def calculate_rm_parameters(frequencies):
    """Calculate RM synthesis parameters from frequency array."""
    # Convert frequency to wavelength
    c = 2.99792458e8  # Speed of light in m/s
    wavelengths = c / frequencies  # wavelengths in meters
    lambda_sq = wavelengths**2  # lambda squared
    
    # Sort by wavelength (descending frequency order)
    lambda_sq = np.sort(lambda_sq)
    
    # Calculate parameters
    lambda_min = np.min(wavelengths)
    lambda_max = np.max(wavelengths)
    lambda_sq_min = np.min(lambda_sq)
    lambda_sq_max = np.max(lambda_sq)
    
    # Channel width in lambda^2 space
    delta_lambda_sq = np.median(np.diff(lambda_sq))
    
    # Total bandwidth in lambda^2
    total_lambda_sq_range = lambda_sq_max - lambda_sq_min
    
    # RM synthesis parameters
    # FWHM of RMTF (rotation measure transfer function)
    fwhm_rmtf = 2 * np.sqrt(3) / total_lambda_sq_range
    
    # Maximum Faraday depth (to first null) - ?
    #phi_max = np.pi / total_lambda_sq_range
    
    # Absolute maximum Faraday depth (Nyquist sampling limit)
    phi_abs_max = np.sqrt(3) / delta_lambda_sq
    
    # Largest scale in Faraday space (maximum sensitivity)
    max_scale = np.pi / lambda_sq_min
    
    return {
        'frequencies': frequencies,
        'wavelengths': wavelengths,
        'lambda_min': lambda_min,
        'lambda_max': lambda_max,
        'lambda_sq_min': lambda_sq_min,
        'lambda_sq_max': lambda_sq_max,
        'delta_lambda_sq': delta_lambda_sq,
        'total_lambda_sq_range': total_lambda_sq_range,
    #    'phi_max': phi_max,
        'phi_abs_max': phi_abs_max,
        'fwhm_rmtf': fwhm_rmtf,
        'max_scale': max_scale,
        'n_channels': len(frequencies)
    }

def print_results(params):
    """Print the calculated RM parameters."""
    print("=" * 60)
    print("RM SYNTHESIS PARAMETERS")
    print("=" * 60)
    print(f"Number of channels: {params['n_channels']}")
    print(f"Frequency range: {params['frequencies'].min()/1e9:.3f} - {params['frequencies'].max()/1e9:.3f} GHz")
    print(f"Wavelength range: {params['lambda_max']:.4f} - {params['lambda_min']:.4f} m")
    print(f"λ² range: {params['lambda_sq_min']:.6f} - {params['lambda_sq_max']:.6f} m²")
    print(f"Δλ² (median): {params['delta_lambda_sq']:.8f} m²")
    print(f"Total λ² range: {params['total_lambda_sq_range']:.6f} m²")
    print("-" * 60)
    print("FARADAY DEPTH PARAMETERS:")
    #print(f"Maximum φ (half-width): ±{params['phi_max']:.2f} rad/m²")
    print(f"Absolute maximum φ (to null): ±{params['phi_abs_max']:.2f} rad/m²")
    print(f"RMTF FWHM: {params['fwhm_rmtf']:.2f} rad/m²")
    print(f"Largest detectable scale: {params['max_scale']:.2f} rad/m²")
    print("-" * 60)
    print("RECOMMENDED RM SYNTHESIS SETTINGS:")
    print(f"Up to a φ range: -{params['phi_abs_max']:.1f} to +{params['phi_abs_max']:.1f} rad/m²")
    print(f"Δφ step < {params['fwhm_rmtf']/3:.2f} rad/m² (FWHM/3)")

def main():
    """Main function."""
    args = parse_arguments()
    
    print(f"Reading measurement set: {args.ms_path}")
    
    # Read frequency information
    frequencies = read_ms_spectral_info(args.ms_path)
    
    # Calculate RM parameters
    params = calculate_rm_parameters(frequencies)
    
    # Print results
    print_results(params)

if __name__ == "__main__":
    main()