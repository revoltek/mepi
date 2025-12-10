#!/usr/bin/env python3

import argparse
import re, glob, sys
from astropy.io import fits
import numpy as np
from astropy import units as u
from astropy.convolution import convolve_fft

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
        print(f"{self.label}: Loading channel: Q={filename_q}, U={filename_u} - freq: {self.frequency*1e-6:.0e} MHz - \
              beam: {self.beam[0]*3600:.0f}\" x {self.beam[1]*3600:.0f}\" @ {self.beam[2]:.1f}°")

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
        
        Parameters:
        -----------
        target_beam : tuple
            Target beam parameters (major, minor, pa) in degrees
        """
        
        # Create beam objects
        current_beam = Beam(major=self.beam[0]*u.deg, 
                           minor=self.beam[1]*u.deg, 
                           pa=self.beam[2]*u.deg)
        target = Beam(major=target_beam[0]*u.deg, 
                     minor=target_beam[1]*u.deg, 
                     pa=target_beam[2]*u.deg)
        
        # Calculate convolution kernel needed
        try:
            convolve_beam = target.deconvolve(current_beam)
        except ValueError:
            print(f"{self.label}: Target beam smaller than current beam, skipping convolution")
            return
        
        # Get pixel scale from header
        pixscale = abs(self.header_q.get('CDELT1', self.header_q.get('CD1_1', 1.0))) * u.deg
        
        # Create 2D Gaussian kernel
        kernel = convolve_beam.as_kernel(pixscale)
        
        # Convolve Q and U data
        
        # Handle 4D data (squeeze out degenerate axes for convolution)
        q_2d = np.squeeze(self.data_q)
        u_2d = np.squeeze(self.data_u)
        
        print(f"{self.label}: Convolving from {current_beam.major.to(u.arcsec):.2f} x {current_beam.minor.to(u.arcsec):.2f} to {target.major.to(u.arcsec):.2f} x {target.minor.to(u.arcsec):.2f}")
        
        q_convolved = convolve_fft(q_2d, kernel, preserve_nan=True, allow_huge=True)
        u_convolved = convolve_fft(u_2d, kernel, preserve_nan=True, allow_huge=True)
        
        # Restore original shape
        self.data_q = q_convolved.reshape(self.data_q.shape)
        self.data_u = u_convolved.reshape(self.data_u.shape)
        
        # Update beam parameters in headers
        self.beam = target_beam
        for header in [self.header_q, self.header_u]:
            header['BMAJ'] = target_beam[0]
            header['BMIN'] = target_beam[1]
            header['BPA'] = target_beam[2]

def create_cube(data_list, outputfile, header):
    """
    Create a FITS cube from a list of input FITS files.
    Stacks 4D cubes along the frequency axis (axis 3 in FITS).
    Preserves the degenerate Stokes axis (axis 4 in FITS).
    """
    for i, data in enumerate(data_list):
        data_list[i] = data[0, 0, :, :]  # Extract spatial plane (ny, nx)
    # Stack along frequency axis
    cube_data = np.stack(data_list, axis=0)
    
    # Reshape to 4D: (1, 1, ny, nx) -> (1, nfreq, ny, nx)
    # Then transpose to FITS order: (1, 1, ny, nx, nfreq) requires shape (1, 1, ny, nx) repeated
    # FITS axis order is (naxis4, naxis3, naxis2, naxis1) = (stokes, freq, dec, ra)
    # NumPy order is reversed: (stokes, freq, ny, nx)
    cube_data = cube_data[np.newaxis, :, :, :]  # (1, nfreq, ny, nx)
    
    nfreq = cube_data.shape[1]
    header['NAXIS3'] = nfreq

    hdu = fits.PrimaryHDU(cube_data, header=header)
    hdu.writeto(outputfile, overwrite=True)
    print(f"Created cube: {outputfile} with shape {cube_data.shape}")

def main():
    parser = argparse.ArgumentParser(
        description="Organize WSClean Stokes Q and U images, optionally normalize by Stokes I"
    )
    parser.add_argument("name_qu", type=str, help="Base name for Q/U files (e.g., 'img/test')")
    #parser.add_argument("--name-i", "-i", help="Name for Stokes I files (e.g. 'img/testI')")
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
    channels = []
    for filename in sorted(glob.glob(f"{args.name_qu}*.fits")):
        match_q = pattern_q.match(filename)
        if match_q:
            filename_q = filename
            filename_u = re.sub(r"-Q-image\.fits$", "-U-image.fits", filename)
            channels.append(Channel(filename_q, filename_u, label=filename))
     
    if len(channels) == 0:
        print("No Q/U files found matching the pattern.")
        sys.exit(1)

    if args.convolve:
        from radio_beam import Beam
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
            flagged_indices = []
            for i, noise in enumerate(noises):
                if noise > threshold:
                    print(f"Flagging channel {i} with noise {noise:.6e}")
                    flagged_indices.append(i)
            
            # Remove flagged channels in reverse order to preserve indices
            for i in reversed(flagged_indices):
                channels.pop(i)
                noises.pop(i)

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