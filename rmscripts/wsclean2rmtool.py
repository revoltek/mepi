#!/usr/bin/env python3

import argparse
import re, glob, sys
from astropy.io import fits
import numpy as np

class Channel:
    def __init__(self, filename_q, filename_u, label=""):
        self.label = label
        self.filename_q = filename_q
        self.filename_u = filename_u
        with fits.open(filename_q) as hdul:
            self.data_q = hdul[0].data
            self.header_q = hdul[0].header.copy()
        with fits.open(filename_u) as hdul:
            self.data_u = hdul[0].data
            self.header_u = hdul[0].header.copy()
        self.frequency = self.header_q.get('CRVAL3', 0.0)
        self.beam = (self.header_q.get('BMAJ', 0.0),
                     self.header_q.get('BMIN', 0.0),
                     self.header_q.get('BPA', 0.0))
        print(f"{self.label}: Loading channel: Q={filename_q}, U={filename_u} - freq: {self.frequency} - \
              beam: {self.beam[0]*3600:.2f}\" x {self.beam[1]*3600:.2f}\" @ {self.beam[2]:.2f}°")

    def get_rms(self, max_iter=10, sigma_clip=5.0, convergence_tol=0.01):
        """
        Estimate the RMS noise of the data using the median absolute deviation.
        Iteratively removes pixels above sigma_clip sigma until convergence.
        Calculates RMS for both Q and U data and returns the average.
        
        Parameters:
        -----------
        max_iter : int
            Maximum number of iterations (default: 10)
        sigma_clip : float
            Sigma threshold for clipping (default: 5.0)
        convergence_tol : float
            Fractional convergence tolerance for RMS change (default: 0.01 = 1%)
        
        Returns:
        --------
        rms : float
            Average of estimated RMS noise from Q and U data
        """
        rms_q = self._estimate_rms_single(self.data_q, max_iter, sigma_clip, convergence_tol, "Q")
        rms_u = self._estimate_rms_single(self.data_u, max_iter, sigma_clip, convergence_tol, "U")
        rms_avg = (rms_q + rms_u) / 2.0
        print(f"{self.label}: Average RMS: {rms_avg:.6e} (Q: {rms_q:.6e}, U: {rms_u:.6e})")
        return rms_avg
    
    def _estimate_rms_single(self, data, max_iter, sigma_clip, convergence_tol, stokes):
        """
        Helper function to estimate RMS for a single data array.
        """
        mask = np.ones(data.shape, dtype=bool)
        rms_prev = None
        
        for iteration in range(max_iter):
            masked_data = data[mask]
            median = np.median(masked_data)
            mad = np.median(np.abs(masked_data - median))
            rms = 1.4826 * mad  # Convert MAD to RMS
            
            # Check convergence using fractional change
            if rms_prev is not None:
                rms_change = abs(rms - rms_prev) / rms_prev
                if rms_change < convergence_tol:
                    print(f"{self.label} ({stokes}): RMS converged after {iteration + 1} iterations: {rms:.6e}")
                    break
            
            # Update mask: remove pixels above sigma_clip * rms
            mask = np.abs(data - median) < sigma_clip * rms
            rms_prev = rms
        else:
            print(f"{self.label} ({stokes}): RMS estimation reached max iterations ({max_iter}): {rms:.6e}")
        
        return rms
    
    def convolve_to_resolution(self, target_beam):
        """
        Convolve the Q and U images to the target beam resolution.
        Placeholder function - actual convolution implementation is not provided.
        
        Parameters:
        -----------
        target_beam : tuple
            Target beam parameters (major, minor, pa)
        """
        print(f"Convolving to target beam: {target_beam} NOT IMPLEMENTED")
        # Placeholder: actual convolution code would go here
        pass

def create_cube(data_list, outputfile, header):
    """
    Create a FITS cube from a list of input FITS files.
    Stacks 4D cubes along the frequency axis (axis 3 in FITS).
    Preserves the degenerate Stokes axis (axis 4 in FITS).
    """
    # Stack along frequency axis
    cube_data = np.stack(data_list, axis=0)
    
    # Reshape to 4D: (1, 1, ny, nx) -> (1, nfreq, ny, nx)
    # Then transpose to FITS order: (1, 1, ny, nx, nfreq) requires shape (1, 1, ny, nx) repeated
    # FITS axis order is (naxis4, naxis3, naxis2, naxis1) = (freq, stokes, dec, ra)
    # NumPy order is reversed: (freq, stokes, ny, nx)
    cube_data = cube_data[np.newaxis, :, :, :]  # (1, nfreq, ny, nx)
    
    hdu = fits.PrimaryHDU(cube_data, header=header)
    hdu.writeto(outputfile, overwrite=True)
    print(f"Created cube: {outputfile} with shape {cube_data.shape}")

def main():
    parser = argparse.ArgumentParser(
        description="Organize WSClean Stokes Q and U images, optionally normalize by Stokes I"
    )
    parser.add_argument("name_qu", type=str, help="Base name for Q/U files (e.g., 'img/test')")
    parser.add_argument("--name-i", "-i", help="Name for Stokes I files (e.g. 'img/testI')")
    parser.add_argument("--output", "-o", default="", help="Base name for output files or use StokesQ.fits and StokesU.fits")
    parser.add_argument("--freq-file", "-f", default="", help="Output file for frequency list")
    parser.add_argument("--noise-file", "-n", default="", help="Estimate noise per channel and create noise files")
    parser.add_argument("--convolve", "-c", action="store_true", help="Convolve Q/U images to minimum resolution before stacking")
    parser.add_argument("--flag", type=float, default=0, help="Flag channels with noise above this threshold (in sigma) - default:0 (no flagging)")

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
    channels = []

    for filename in sorted(glob.glob(f"{args.name_qu}*.fits")):
        match_q = pattern_q.match(filename)
        match_u = pattern_u.match(filename)
        if match_q and match_u:
            channels.append(Channel(match_q.group(0), match_u.group(0), label=filename))
     
    if args.convolve:
        print("Convolving channels to common resolution...")
        target_beam = max([ch.beam for ch in channels], key=lambda b: b[0])  # Max major axis
        for channel in channels:
            channel.convolve_to_resolution(target_beam)

    if args.noise_file or args.flag > 0:
        print("Estimating noise per channel...")
        noises = []
        for channel in channels:
            rms = channel.get_rms()
            noises.append(rms)

        if args.noise_file:
            print(f"Writing noise values to {args.noise_file}")
            np.savetxt(args.noise_file, noises, delimiter="")

        if args.flag > 0:
            mean_noise = np.mean(noises)
            std_noise = np.std(noises)
            threshold = mean_noise + args.flag * std_noise
            print(f"Flagging channels with noise above {threshold:.6e} (mean: {mean_noise:.6e}, std: {std_noise:.6e})")
            for i, noise in enumerate(noises):
                if noise > threshold:
                    print(f"Flagging channel {i} with noise {noise:.6e}")

    print("Creating Stokes Q and U cubes...")
    create_cube([channel.data_q for channel in channels], output_q, channels[0].header_q)
    create_cube([channel.data_u for channel in channels], output_u, channels[0].header_u)
    
    freqs = []
    if args.freq_file:
        print(f"Writing frequency list to {args.freq_file}")
        for channel in channels:
            freqs.append(channel.frequency)
        # Extract frequencies from one of the Q files
        np.savetxt(args.freq_file, freqs, delimiter="")

    print("\nDone!")


if __name__ == "__main__":
    main()