import os, glob
import casatasks as casa
from mepi import lib_mepi, lib_ms, lib_log, lib_runcode, lib_cfg, lib_walker, lib_sol

log = lib_log.log
cfg = lib_cfg.cfg
w = lib_walker.Walker("selfcal")

def run():

    ms_full_file = cfg['ms_full']
    ms_tgt_file = os.path.join(cfg['path_ms'], os.path.basename(ms_full_file.rstrip('/').rstrip('\\')).replace('.MS', '_tgt_avg.MS'))
    ms_tgt = lib_ms.MS(ms_tgt_file)

    if cfg['ref_ant'] is not None:
        ref_ant = cfg['ref_ant']
        log.info(f"Using user-specified reference antenna: {ref_ant}")
    else:
        ref_ant = ms_tgt.find_reference_antenna()
        log.info(f"Using automatically selected reference antenna: {ref_ant}")
    pixelscale = round(0.7 * (2.4/ms_tgt.freq_center*1e-9), 1) # arcsec

    # selfcal only on scalar amp and possibly diag phase.
    # If diag phase needed, only for stokes I and consider parang is amp rot matrix and doesn't commute
    casa.flagmanager(vis=ms_tgt_file, mode='save', versionname='PreSelfcal')
    casa.applycal(vis=ms_tgt_file, flagbackup=False, parang=True)
    for cc in range(10):
        log.info(f'Self-calibration cycle {cc+1}')
        # ok for m87 sband
        lib_runcode.run_command(f'{lib_runcode.wsclean_command} -name IMG/{Targets}-selfcal-c{cc}  -update-model-required -pol I \
            -reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 \
            -size 5000 5000 -scale {pixelscale}arcsec -weight briggs -0.2  \
            -niter 1000000 -mgain 0.8 \
            -join-channels -channels-out 12 -deconvolution-channels 3 -fit-spectral-pol 3 \
            -multiscale -multiscale-scales 1,4,8,16,32,64,128,256 \
            -auto-mask 5 -auto-threshold 3 \
            {ms_tgtavg_file}', cfg, logname='wsclean')
        
        if cc == 1 or cc == 3:
            log.info(f'Flagging data for self-calibration cycle {cc+1}')
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --corr XX,YY --png './PLOTS/Tgt-c{cc}.png' {ms_tgtavg_file}", cfg, logname='shadems')
            casa.flagdata(vis=ms_tgtavg_file, mode="rflag", datacolumn="residual", quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False)
            casa.flagdata(vis=ms_tgtavg_file, mode='extend', datacolumn='residual', growtime=80, growfreq=80, flagbackup=False, growaround=True, flagnearfreq=True)
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --corr XX,YY --png './PLOTS/Tgt-c{cc}-flag.png' {ms_tgtavg_file}", cfg, logname='shadems')

        log.info(f'Calibrating data for self-calibration cycle {cc+1}')
        casa.gaincal(vis=ms_tgtavg_file, caltable='CASA_Tables/selfcal%02i.K' %cc, gaintype='K', solint='32s', refant=ref_ant, parang=False)
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table CASA_Tables/selfcal{cc:02d}.K --plotname ./PLOTS/target-K-i{cc:02d}.png >> ragavi.log', cfg, logname='ragavi')
        # plotms(vis='CASA_Tables/selfcal%02i.K' %cc, coloraxis='antenna1', xaxis='time', yaxis='delay')
        casa.gaincal(vis=ms_tgtavg_file, caltable='CASA_Tables/selfcal%02i.Gp' %cc,  gaintype='G', calmode='p', solint='32s', refant=ref_ant, parang=False,
                    gaintable=['CASA_Tables/selfcal%02i.K' %cc])
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table CASA_Tables/selfcal{cc:02d}.Gp --yaxis phase --plotname ./PLOTS/target-Gp-i{cc:02d}.png >> ragavi.log', cfg, logname='ragavi')
        # plotms(vis='CASA_Tables/selfcal%02i.Gp' %cc, coloraxis='antenna1', xaxis='time', yaxis='phase', xconnector='line')
        casa.gaincal(vis=ms_tgtavg_file, caltable='CASA_Tables/selfcal%02i.Ga' %cc, gaintype='T', calmode='a', solint='128s', refant=ref_ant, solnorm=True, parang=True,
                    gaintable=['CASA_Tables/selfcal%02i.K' %cc, 'CASA_Tables/selfcal%02i.Gp' %cc])
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table CASA_Tables/selfcal{cc:02d}.Ga --yaxis amplitude --plotname ./PLOTS/target-Ga-i{cc:02d}.png >> ragavi.log', cfg, logname='ragavi')
        # plotms(vis='CASA_Tables/selfcal%02i.Ga' %cc, coloraxis='antenna1', xaxis='time', yaxis='amp', xconnector='line')
        #casa.bandpass(vis=ms_tgtavg_file, caltable='selfcal%02i.B' %cc, combine='', solint='300s', gaintable=['selfcal%02i.G' %cc, 'selfcal%02i.K' %cc], refant='m002', parang=False)
        casa.applycal(vis=ms_tgtavg_file, flagbackup=False, parang=True,
                    gaintable=['CASA_Tables/selfcal%02i.K' %cc, 'CASA_Tables/selfcal%02i.Gp' %cc, 'CASA_Tables/selfcal%02i.Ga' %cc])
        lib_mepi.print_flags(ms_tgtavg_file)

    # pol cleaning - possible problem with -squared-channel-joining when using -multiscale
    #os.system(f'{wsclean_command} -name IMG/{Targets}-selfcal-pol -update-model-required -pol IQUV '
    #          f'-reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 -baseline-averaging 12 '
    #          f'-size 2500 2500 -scale {pixelscale}arcsec -weight briggs -0.2 -minuv-l 80.0 '
    #          f'-niter 1000000 -mgain 0.7 '
    #         f'-join-channels -channels-out 32 -deconvolution-channels 6 -fit-spectral-pol 3 -squared-channel-joining '
    #          f'-multiscale -multiscale-scales 1,4,8,16,32,64,128,256 '
    #          f'-auto-threshold 3 -fits-mask m87-07asec-2500.fits '
    #          f'{ms_tgtavg_file} > wsclean_{Targets}-selfcal.log')

    # wsclean with rm
    restoring_beam = 6.0 # arcsec - this is ok for S1 band
    lib_runcode.run_command(f'{lib_runcode.wsclean_command} -name IMG/{Targets}-rm -no-update-model-required -pol QU '
            f'-reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 -baseline-averaging 12 '
            f'-size 5000 5000 -scale 2arcsec -weight briggs -0.2 -minuv-l 80.0 -beam-size {restoring_beam} -taper-gaussian {restoring_beam}arcsec '
            f'-niter 25000 -mgain 0.75 -nmiter 12 '
            f'-join-channels -channels-out 125 -join-polarizations -squared-channel-joining -fit-rm '
            f'{ms_tgtavg_file}', cfg, logname=f'wsclean_{Targets}-selfcal')
    # relative I stokes for fractional Pol
    lib_runcode.run_command(f'{lib_runcode.wsclean_command} -name IMG/{Targets}-rm -no-update-model-required -pol I '
            f'-reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 -baseline-averaging 12 '
            f'-size 5000 5000 -scale 2arcsec -weight briggs -0.2 -minuv-l 80.0 -beam-size {restoring_beam} -taper-gaussian {restoring_beam}arcsec '
            f'-niter 25000 -mgain 0.75 -nmiter 12 '
            f'-join-channels -channels-out 125 '
            f'{ms_tgtavg_file}', cfg, logname=f'wsclean_{Targets}-selfcal')