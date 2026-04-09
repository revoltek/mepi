import os, glob
import casatasks as casa
from mepi import lib_mepi, lib_ms, lib_log, lib_runcode

log = lib_log.logger

def run(cfg):

    ms_full_file = cfg['ms_full']
    ms_tgt_file = os.path.join(cfg['path_ms'], os.path.basename(ms_full_file.rstrip('/').rstrip('\\')).replace('.MS', '_tgt.MS'))
    ms_tgtavg_file = ms_tgt_file.replace('.MS', '_avg.MS')
    ms_full = lib_ms.MS(ms_full_file)

    Targets = ','.join(ms_full.find_targets())

    #########################################################################
    # Split the target
    log.info('Splitting target...')
    if not os.path.exists(ms_tgt_file):
        spw_selection = ms_full.get_spw_noedges() # select all channels but the first and last 5%
        casa.split(vis = ms_full_file, outputvis = ms_tgt_file, field = f"{Targets}", datacolumn = 'data', spw = spw_selection)
        log.info(f'Target split and saved in {ms_tgt_file}')
    else:
        log.info('Target has already been split previously')
    lib_mepi.print_flags(ms_tgt_file)
    ms_tgt = lib_ms.MS(ms_tgt_file)

    if cfg['ref_ant'] is not None:
        ref_ant = cfg['ref_ant']
        log.info(f"Using user-specified reference antenna: {ref_ant}")
    else:
        ref_ant = ms_tgt.find_reference_antenna()
        log.info(f"Using automatically selected reference antenna: {ref_ant}")
    pixelscale = round(0.7 * (2.4/ms_tgt.freq_center*1e-9), 1) # arcsec

    #############################################################################
    # First apply the calibration from crosscal to the target, flag, and split the target for self-calibration
    log.info('Applying calibration to target...')
    casa.applycal(vis=ms_tgt_file, parang=False, flagbackup=False, interp=['linear,linearflag', 'linear,linearflag', 'linear,linearflag'], \
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Xf_tab_ambcorr'], tab['Ga_tab'], tab['Ksec_tab'], tab['Gpsec_tab'], tab['Tsec_tab']])
    lib_mepi.print_flags(ms_tgt_file)

    # Standard flagging for shadowing, zero-clip, and auto-correlation
    log.info('Flagging target...')
    casa.flagdata(vis=ms_tgt_file, flagbackup=False, mode='shadow')
    casa.flagdata(vis=ms_tgt_file, flagbackup=False, mode='manual', autocorr=True)
    casa.flagdata(vis=ms_tgt_file, flagbackup=False, mode='clip', clipzeros=True, clipminmax=[0.0, 1000.0]) # high for virgo A, 100 is ok for others
    if ms_tgt.band == "UHF" or ms_tgt.band == "L": 
        lib_runcode.run_command(f"{lib_runcode.mask_ms_command} --mask {lib_runcode.rfimask} --accumulation_mode or --memory 4096 --uvrange 0~1000 --statistics {ms_tgt_file}", cfg, logname='rfimask')
    lib_mepi.print_flags(ms_tgt_file)
    casa.flagmanager(vis=ms_tgt_file, mode='save', versionname='PreAoflagger')
    lib_runcode.run_command(f"{lib_runcode.aoflagger_command} -strategy {lib_runcode.aoflagger_strategy1} -column CORRECTED_DATA {ms_tgt_file}", cfg, logname='aoflagger')
    lib_runcode.run_command(f"{lib_runcode.aoflagger_command} -strategy {lib_runcode.aoflagger_strategy1} -column CORRECTED_DATA {ms_tgt_file}", cfg, logname='aoflagger') # twice
    lib_mepi.print_flags(ms_tgt_file)

    # Split the target averaged in freq and time
    log.info('Splitting target avg...')
    if not os.path.exists(ms_tgtavg_file):
        casa.split(vis = ms_tgt_file, outputvis = ms_tgtavg_file, field = f"{Targets}", datacolumn = 'corrected',
                    width=cfg['freqbin'], timebin=cfg['timebin'])
        log.info(f'Target split and saved in {ms_tgtavg_file}')
    else:
        log.info('Target has already been split previously')
    lib_mepi.print_flags(ms_tgtavg_file)

    ##########################################################################################
    # selfcal only on scalar amp and possibly diag phase.
    # If diag phase needed, only for stokes I and consider parang is amp rot matrix and doesn't commute
    casa.flagmanager(vis=ms_tgtavg_file, mode='save', versionname='PreSelfcal')
    casa.applycal(vis=ms_tgtavg_file, flagbackup=False, parang=True)
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

    #######################################################################################################################################
    #######################################################################################################################################
    #######################################################################################################################################
    # Selfcal with DP3 - as an alternative method

    losoto_parset = os.path.join(cfg['mepi_dir'], 'parsets/losoto-plot.parset')
    dp3_sol_parset = os.path.join(cfg['mepi_dir'], 'parsets/DP3-sol.parset')
    dp3_cor_parset = os.path.join(cfg['mepi_dir'], 'parsets/DP3-cor.parset')

    # NOTE: DP3 still not work on multi-scan obs
    casa.applycal(vis=ms_tgtavg_file, flagbackup=False, parang=True) # apply parang and split corrected data
    scans = []
    for scan in casa.listobs(ms_tgtavg_file):
        if 'scan' in scan: scans.append(int(scan.split('_')[1]))
    for scan in scans:
        casa.split(vis = ms_tgtavg_file, outputvis = f'{ms_tgtavg_file.replace(".MS", f"-scan{scan}.MS")}', field = f"{Targets}", datacolumn = 'corrected', scan=str(scan))
    mss = sorted(glob.glob(f'{ms_tgtavg_file.replace(".MS", "")}-scan*.MS'))
    #for ms in mss:
    #    os.system(f'wsclean  -predict -padding 1.8 -j 64 -name {imgname} -channels-out 32 {ms} >> wsclean.log')
    for i in range(10):
        print(f'Cycle {i}: imaging...')
        os.system(f'{wsclean_command} -name IMG/{Targets}-selfcal-c{i}  -update-model-required -pol I \
            -reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 \
            -size 2500 2500 -scale {pixelscale}arcsec -weight briggs -0.2  -minuv-l 80.0 \
            -niter 1000000 -mgain 0.6 \
            -join-channels -channels-out 32 -deconvolution-channels 6 -fit-spectral-pol 3 \
            -multiscale -multiscale-scales 1,3,9,27,90,270 \
            -auto-threshold 3 -fits-mask m87-07asec-2500.fits \
            {' '.join(mss)} > wsclean_{Targets}-selfcal.log')
        for ms in mss:
            print(f'Cycle {i}: Solving+correcting {ms}...')
            # solve
            os.system(f'DP3 {dp3_sol_parset} msin={ms} msout=. sol.h5parm={ms}/ph-{i}.h5 sol.mode=diagonalphase sol.solint=1 sol.nchan=1 sol.smoothnessconstraint=10e6 >> DP3.log')
            os.system(f'DP3 {dp3_sol_parset} msin={ms} msout=. sol.h5parm={ms}/amp-{i}.h5 sol.mode=scalaramplitude sol.solint=20 sol.nchan=1 sol.smoothnessconstraint=30e6 >> DP3.log')
            # correct
            os.system(f'DP3 {dp3_cor_parset} msin={ms} msout=. cor1.parmdb={ms}/ph-{i}.h5 cor2.parmdb={ms}/amp-{i}.h5 >> DP3.log')
        # losoto
        os.system(f'H5parm_collector.py -V -s sol000 -o cal-ph-{i}.h5 '+' '.join([f'{ms}/ph-{i}.h5' for ms in mss]))
        os.system(f'losoto cal-ph-{i}.h5 {losoto_parset} >> losoto.log && mv plots plots-{i}')
        os.system(f'H5parm_collector.py -V -s sol000 -o cal-amp-{i}.h5 '+' '.join([f'{ms}/amp-{i}.h5' for ms in mss]))
        os.system(f'losoto cal-amp-{i}.h5 {losoto_parset} >> losoto.log && mv plots plots-{i}')
