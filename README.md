# mepi

A container file with all required softwares and packages can be found at docker://tpasini/pol_meerkat:latest

To run the singularity one can use (for instance):

`singularity build ~/MeerKATpol.simg docker://tpasini/pol_meerkat:latest`

`singularity shell --pid --writable-tmpfs --cleanenv -B/homes/fdg,/local/work/fdg,/iranet/groups/ulu/fdg/,/iranet/groups/lofar/containers/ ~/MeerKATpol.simg`

For L and S band the data must have X and Y flipped using the script "correct_parang.py" by Ben Hugo
https://github.com/bennahugo/LunaticPolarimetry/blob/master/correct_parang.py

`python correct_parang.py -f {field_id} --noparang  --applyantidiag MSFILE`

the script writes the output in the CORRECTED_DATA column and needs to be run for each field, so the corrected data needs to be split after the correction.

plots are done with shadems https://github.com/ratt-ru/shadeMS

### RM

## Wsclean

Use -fit-rm
Use -squared-channel-join
Do not use -multiscale (it cannot work with square channel join)

**todo**: find the right number of channels, larger than needed so flagging is better, important in Lband, ideally a factor 3-4 larger than the -s later
rm_max = sqrt(3)/deltaLmabda**2 (a freq + basse)

Example:
`wsclean -j 64 -name img/xxx -no-update-model-required -pol QU -reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 -baseline-averaging 12 -size 1500 1500 -scale 2arcsec -weight briggs -0.5 -minuv-l 80.0 -beam-size XXX -taper-gaussian XXXarcsec -niter 25000 -mgain 0.75 -nmiter 12 -join-channels -channels-out XXX -join-polarizations -squared-channel-joining -fit-rm XXX.MS > wsclean_QU.log`

**note**: best would be not to force a common restoring beam (and no tapering) but to convolve in the preprocess to the same beam

## Preprocess
create cubes, freq and noise files and flag bad chans:
`wsclean2rmtool.py -f freqs.dat -n noise.dat --flag 5 ../IMG/xxx`

TODO: rebin of images
TODO: add histogram of noise and beam size
TODO: add a check on the beam size to exclude those that are too large

## Do synthesis + clean
input: FITS files, and an ASCII file containing a list of channel frequencies

todo: add noise, add spidx, define -l and -s properly

https://github.com/CIRADA-Tools/RM-Tools/wiki/RMsynth3D
`rmsynth3d -v -l 200 -s 10 -o rmtoolsynth StokesQ.fits StokesU.fits freqs.dat`

-r? Optimise the resolution of the RMSF (as per Rudnick & Cotton).

Find the noise in P to limit the clean, the region should be away from sources:
`findPnoise.py rmtoolsynthFDF_tot_dirty.fits region.reg`

Do clean, set the noise to 6x the output of the previous script
`rmclean3d -v --ncores 64 -o rmtoolclean -c xxx rmtoolsynthFDF_tot_dirty.fits rmtoolsynthRMSF_tot.fits`

## Postprocess
create the smoothed RM map + angle map
`rmtools_peakfitcube -v -p rmtoolsynthFDF_tot_dirty.fits freqs.dat rmtoolpfc`

Output:
**ampPeakPIfitEff** and **dAmpPeakPIfit**: Polarized intensity found by fitting the peak. Corrected for polarization bias if snrPIfit > 5. Correction formula is PI_eff = sqrt(PI^2 - 2.3*dFDFth^2). And error in polarized intensity.

**polAngle0Fit_deg** and **dPolAngle0Fit_deg**: Derotated polarization angle (i.e. the angle at the point of emission), calculated by subtracting phiPeakPIfit_rm2*lam0Sq_m2 from the polarization angle above and error.

**phiPeakPIfit_rm2** and **dphiPeakPIfit_rm2**: Faraday depth found by fitting the peak and error.

**snrPIfit**: Signal to noise ratio in the (fitted) polarized intensity. Calculated as the ratio of polarized intensity (without bias correction) to theoretical noise.

**dFDFcorMAD**: Estimated noise in the Faraday depth spectrum, computed using the median deviation from the median (MADFM) of polarized intensity in the spectrum as described above, and then correcting to be equivalent to a Gaussian sigma.

QU fitting canbe run only in 1D

TODO: tagli, I>3sigma, spettroFaraday>6sigma (calcolato sul cubo dopo il clean usando i canali inizio e fine) - hales2012
