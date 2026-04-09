import os, sys
import casatasks as casa
from casatools import msmetadata
import numpy as np
from mepi import lib_mepi, lib_ms, lib_log, lib_runcode

log = lib_log.logger

def run(cfg):
    """
    Run the cross-calibration part of the mepi pipeline on a given config file.
    """

    # set filenames
    ms_full_file = cfg['ms_full']
    ms_cal_file = os.path.join(cfg['path_ms'], os.path.basename(ms_full_file.rstrip('/').rstrip('\\')).replace('.MS', '_cal.MS'))
    ms_full = lib_ms.MS(ms_full_file)

    # Name your gain tables
    tab = {'Ginit_tab' : 'gain_init_bp.cal',
        'K_tab' : 'delay_bp.cal',
        'B_tab' : 'bandpass.cal',
        'Gp_tab': 'gain_p_bp.cal',
        'Ga_tab': 'gain_a_bp.cal',
        'Ksec_tab' : 'delay_sec.cal',
        'Tsec_tab' : 'T_sec.cal',
        'Gpsec_tab' : 'gain_p_sec.cal',
        'Kpol_tab' : 'delay_pol.cal',
        'Gppol_tab' : 'gain_p_pol.cal',
        # Pol cal tables
        'Xf_tab': 'Xf.cal',
        'Xf_tab_ambcorr': 'Xf_ambcorr.cal',
        'Df_tab': 'Df.cal'}

    # Fix some variables
    for name in tab:
        tab[name] = os.path.join(cfg['path_sols'], tab[name])

    PhaseCal = ','.join(ms_full.find_phase_calibrator())
    PolCal = ','.join(ms_full.find_pol_calibrator())
    BandPassCal = ','.join(ms_full.find_bandpass_calibrator())
    if PhaseCal is None or PolCal is None or BandPassCal is None:
        log.error("Could not find all calibrators in the MS.")
        sys.exit(1)

    ###########################
    # Split the calibrators
    log.info('Splitting calibrators...')
    if not os.path.exists(ms_cal_file):
        spw_selection = ms_full.get_spw_noedges() # select all channels but the first and last 5%
        casa.split(vis = ms_full_file, outputvis = ms_cal_file, field = f"{BandPassCal},{PolCal},{PhaseCal}", 
                   datacolumn = 'data', spw = spw_selection)
        log.info(f'Calibrators split and saved in {ms_cal_file}')
        lib_mepi.print_flags(ms_cal_file)
    else:
        log.info('Calibrators have already been split previously')

    ms_cal = lib_ms.MS(ms_cal_file)
    ms_cal.swap_feeds()
    if cfg['ref_ant'] is not None:
        ref_ant = cfg['ref_ant']
        log.info(f"Using user-specified reference antenna: {ref_ant}")
    else:
        ref_ant = ms_cal.find_reference_antenna()
        log.info(f"Using automatically selected reference antenna: {ref_ant}")

    # field ids for imaging
    msmd = msmetadata()
    msmd.open(ms_cal_file)
    PhaseCal_id = msmd.fieldsforname(PhaseCal)[0]
    PolCal_id = msmd.fieldsforname(PolCal)[0]
    msmd.close()
    # imaging parameters scaled from S1 band
    pixelscale = round(0.7 * (2.4/ms_cal.freq_center*1e-9), 1) # arcsec

    #######################################################################################################################################
    #######################################################################################################################################
    #######################################################################################################################################

    #################################
    # Standard flagging for shadowing, zero-clip, and auto-correlation
    casa.flagdata(vis=ms_cal_file, flagbackup=False, mode='shadow')
    casa.flagdata(vis=ms_cal_file, flagbackup=False, mode='manual', autocorr=True)
    casa.flagdata(vis=ms_cal_file, flagbackup=False, mode='clip', clipzeros=True)#, clipminmax=[0.0, 100.0])

    # clip known bad regions
    if ms_cal.band == "UHF": 
        casa.flagdata(vis=ms_cal_file, flagbackup=False, mode='manual', spw='*:925~945MHz, *:950~960MHz') # UHF bad data
    elif ms_cal.band == "L": 
        casa.flagdata(vis=ms_cal_file, flagbackup=False, mode='manual', spw='*:856~880MHz, *:1658~1800MHz, *:1419.8~1421.3MHz') # suggested by SARAO for Lband
    elif ms_cal.band == "S1": 
        casa.flagdata(vis=ms_cal_file, flagbackup=False, mode='manual', spw='*:2150~2161MHz, *:2312~2324MHz, *:2198~2203MHz, *:2251.5~2253.5MHz, *:2361~2366MHz, *:2491.5~2492.5MHz') # resonances+RFI S1 band
    else:
        log.warning(f"No specific flagging range defined for band {ms_cal.band}, only generic flags applied.")

    if ms_cal.band == "UHF" or ms_cal.band == "L": 
        lib_runcode.run_command(f"{lib_runcode.mask_ms_command} --mask {lib_runcode.rfimask} --accumulation_mode or --memory 4096 --uvrange 0~1000 --statistics {ms_cal_file}", cfg, logname='rfimask')
    lib_mepi.print_flags(ms_cal_file)

    # Set flux density scale
    for cal in set(BandPassCal.split(',')+PolCal.split(',')):
        log.info('Setting model for calibrator %s' % cal)
        if cal == 'J1939-6342':
            if ms_cal.band == "UHF":
                log.info('-- setting UHF-band model for flux calibrator J1939-6342')
                casa.clearcal(vis=ms_cal_file, addmodel=True)
                lib_runcode.run_command(f"{lib_runcode.crystalball_command} {ms_cal_file} -sm {os.path.join(cfg['mepi_dir'], 'parsets/J1939-6342_UHF.txt')} -f {cal} -j 2", cfg, logname='crystalball')
            else:
                casa.setjy(vis = ms_cal_file, field = cal, standard = 'Stevens-Reynolds 2016', usescratch = True)
        elif cal == 'J0408-6545':
            if ms_cal.band == "UHF":
                log.info('-- setting UHF-band model for flux calibrator J0408-6545')
                casa.clearcal(vis=ms_cal_file, addmodel=True)
                lib_runcode.run_command(f"{lib_runcode.crystalball_command} {ms_cal_file} -sm {os.path.join(cfg['mepi_dir'], 'parsets/J0408-6545_UHF.txt')} -f {cal} -j 2", cfg, logname='crystalball')
            else:
                a=-0.9790; b=3.3662; c=-1.1216; d=0.0861
                reffreq,fluxdensity,spix0,spix1,spix2 = lib_mepi.convert_flux_model(np.linspace(0.9,2,200)*1e9,a,b,c,d)
                casa.setjy(vis = ms_cal_file, field = cal, usescratch = True, standard = 'manual', \
                    spix = [spix0, spix1, spix2, 0], fluxdensity = fluxdensity, reffreq = '%f Hz'%(reffreq))
        elif cal == 'J1331+3030':
            if ms_cal.band == "UHF":
                I = 19.27475
                alpha = [-0.42727, -0.14583]
                reffreq = '0.816GHz'
                polfrac = 0.06398396324
                polangle = 0.41598756411
                rm = 0.12
            else:
                I= 14.7172
                alpha= [-0.4507, -0.1798, 0.0357]
                reffreq= '1.47GHz'
                polfrac= 0.098
                polangle = 0.575959
                rm=0.
            casa.setjy(vis=ms_cal_file, field=cal, usescratch = True, standard = 'manual', \
                    fluxdensity=[I,0,0,0], spix=alpha, reffreq=reffreq, polindex=polfrac, polangle=polangle, rotmeas=rm)
            # TODO: add iono corruption to the model
        else:
            print("Unknown calibrator ", cal)
            sys.exit()

    ################################################################
    ### Bandpass calibration
    log.info('Bandpass calibration...')
    casa.flagmanager(vis=ms_cal_file, mode='save', versionname='PreBP')

    lib_runcode.run_command(f"{lib_runcode.aoflagger_command} -strategy {lib_runcode.aoflagger_strategy1} -column DATA {ms_cal_file}", cfg, logname='aoflagger')
    lib_mepi.print_flags(ms_cal_file)
    lib_runcode.run_command(f"{lib_runcode.aoflagger_command} -strategy {lib_runcode.aoflagger_strategy1} -column DATA {ms_cal_file}", cfg, logname='aoflagger') # twice
    lib_mepi.print_flags(ms_cal_file)

    for cc in range(2):
        log.info(f'Calibration cycle {cc+1}')
        # Pre gaincal
        casa.gaincal(vis=ms_cal_file, field=BandPassCal, caltable=tab['Ginit_tab'], gaintype='G', calmode='ap', 
                    refant=ref_ant, solint='int', spw=ms_cal.get_good_spw(), solnorm=True)
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab['Ginit_tab']} --plotname ./PLOTS/cal_BP-Ginit-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['Ginit_tab'], coloraxis='antenna1', xaxis='time', yaxis='phase')
        # Delay calibration (fast to track the ionosphere)
        casa.gaincal(vis=ms_cal_file, field=BandPassCal, caltable=tab['K_tab'], gaintype='K', 
                    gaintable=[tab['Ginit_tab']], refant=ref_ant, solint='32s')
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab['K_tab']} --plotname ./PLOTS/cal_BP-K-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['K_tab'], coloraxis='antenna1', xaxis='time', yaxis='delay') # all plotms should be run separately in a casa session
        # one can now combine the scans and use different B as diagnostics
        casa.bandpass(vis=ms_cal_file, field=BandPassCal, caltable=tab['B_tab'], bandtype='B', 
                    gaintable=[tab['Ginit_tab'], tab['K_tab']], combine='scan', solint='inf', refant=ref_ant)
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab['B_tab']} --plotname ./PLOTS/cal_BP-B-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['B_tab'], coloraxis='antenna1', xaxis='freq', yaxis='amp')
        # plotms(vis=tab['B_tab'], coloraxis='antenna1', xaxis='freq', yaxis='phase')
        casa.gaincal(vis=ms_cal_file, field=BandPassCal, caltable=tab['Ga_tab'], gaintype='G', calmode='a', 
                    gaintable=[tab['B_tab'],tab['K_tab']], interp=['linear,linearflag'], refant=ref_ant)
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab['Ga_tab']} --yaxis amplitude --plotname ./PLOTS/cal_BP-Ga-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['Ga_tab'], coloraxis='antenna1', xaxis='time', yaxis='amp', xconnector='line')
        casa.gaincal(vis=ms_cal_file, field=BandPassCal, caltable=tab['Gp_tab'], gaintype='G', calmode='p', 
                    gaintable=[tab['B_tab'],tab['K_tab']], interp=['linear,linearflag'], refant=ref_ant, solint='32s')
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab['Gp_tab']} --yaxis phase --plotname ./PLOTS/cal_BP-Gp-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['Gp_tab'], coloraxis='antenna1', xaxis='time', yaxis='phase', xconnector='line')

        if cc == 0:
            casa.applycal(vis=ms_cal_file, interp=['linear,linearflag', 'nearest', 'nearest'], flagbackup=False, 
                        gaintable=[tab['B_tab'],tab['K_tab'],tab['Gp_tab'],tab['Ga_tab']]) # apply to all fields for better rfi flagging
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-amp.png' {ms_cal_file}", cfg, logname='shadems')
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:phase --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-ph.png' {ms_cal_file}", cfg, logname='shadems')
            
            lib_runcode.run_command(f"{lib_runcode.aoflagger_command} -strategy {lib_runcode.aoflagger_strategy2} -column CORRECTED_DATA {ms_cal_file}", cfg, logname='aoflagger')
            casa.flagdata(vis=ms_cal_file, mode="rflag", datacolumn="residual", field=BandPassCal, quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False)
            casa.flagdata(vis=ms_cal_file, mode='extend', datacolumn='residual', field=BandPassCal, growtime=80, growfreq=80, flagbackup=False)
            lib_mepi.print_flags(ms_cal_file)

            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-amp-flag.png' {ms_cal_file}", cfg, logname='shadems')
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:phase --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-ph-flag.png' {ms_cal_file}", cfg, logname='shadems')
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x CORRECTED_DATA:phase -y CORRECTED_DATA:amp --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass_ampph-flag.png' {ms_cal_file}", cfg, logname='shadems')

    #############################################
    ### Leackage calibration
    log.info('Leackage calibration...')
    casa.flagmanager(vis = ms_cal_file, mode = 'save', versionname = f'PreLeak')

    # plot and flagging
    casa.applycal(vis=ms_cal_file, interp=['linear,linearflag'], flagbackup=False, 
                field=BandPassCal, gaintable=[tab['B_tab'],tab['K_tab'],tab['Gp_tab'],tab['Ga_tab']])
    lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XY,YX --png './PLOTS/Bandpass-cross-preleak.png' {ms_cal_file}", cfg, logname='shadems')
    #casa.flagdata(vis=ms_cal_file, mode="rflag", datacolumn="corrected", field=BandPassCal, quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False, correlation='XY,YX')
    #casa.flagdata(vis=ms_cal_file, mode='extend', datacolumn="corrected", field=BandPassCal, growtime=80, growfreq=80, flagbackup=False, growaround=True, flagnearfreq=True, correlation='XY,YX')
    #lib_mepi.print_flags(ms_cal_file)
    #os.system(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XY,YX --png './PLOTS/Bandpass-cross-preleak-flag.png' {ms_cal_file} >> shadems.log")

    for cc in range(2):
        log.info(f'Calibration cycle {cc+1}')
        # Leackage on unpol calib
        casa.polcal(vis=ms_cal_file,
            caltable=tab['Df_tab'],field=BandPassCal, poltype='Dflls', solint='inf', combine='scan',
            gaintable=[tab['B_tab'], tab['K_tab'], tab['Gp_tab'], tab['Ga_tab']], interp=['linear,linearflag'], refant=ref_ant)
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab['Df_tab']} --xaxis channel --yaxis amplitude --plotname ./PLOTS/cal_BP-Df-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['Df_tab'], xaxis='frequency', yaxis='amplitude', coloraxis='antenna1')
        casa.flagdata(vis=tab['Df_tab'], mode='tfcrop', datacolumn="CPARAM", quackinterval=0.0, ntime="60s", combinescans=True, timecutoff=4.0, freqcutoff=3.0, usewindowstats="both", flagbackup=False)
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab['Df_tab']} --xaxis channel --yaxis amplitude --plotname ./PLOTS/cal_BP-Df-cc{cc}-flag.png', cfg, logname='ragavi')
        # plotms(vis=tab['Df_tab'], xaxis='frequency', yaxis='amplitude', coloraxis='antenna1')

        if cc == 0:
            casa.applycal(vis=ms_cal_file,field=BandPassCal, interp=['linear,linearflag','linear,linearflag'], flagbackup=False, 
                        gaintable=[tab['B_tab'],tab['Df_tab'],tab['K_tab'],tab['Gp_tab'],tab['Ga_tab']])
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XY,YX --png './PLOTS/Bandpass-cross-postleak.png' {ms_cal_file}", cfg, logname='shadems') # check if the amp are reduced and no big waves/spikes should be there
            casa.flagdata(vis=ms_cal_file, mode="rflag", datacolumn="residual", field=BandPassCal, quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False, correlation='XY,YX')
            casa.flagdata(vis=ms_cal_file, mode='extend', datacolumn="residual", field=BandPassCal, growtime=80, growfreq=80, flagbackup=False, growaround=True, flagnearfreq=True, correlation='XY,YX')
            lib_mepi.print_flags(ms_cal_file)
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XY,YX --png './PLOTS/Bandpass-cross-postleak-flag.png' {ms_cal_file}", cfg, logname='shadems') # check if the amp are reduced and no big waves/spikes should be there
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x ANTENNA1 -y CORRECTED_DATA:real --field {BandPassCal} --corr XY,YX --png './PLOTS/Bandpass-cross-ant.png' {ms_cal_file}", cfg, logname='shadems') #important check for chosing the reference antenna, make sure that no antenna with extreme leakage is chosen  

    # re-do BP to use better flags - usually it gives worse results...
    #log.info(f'Re-doing bandpass calibration with better flags...')
    #casa.bandpass(vis=ms_cal_file, field=BandPassCal, caltable=tab['B_tab'], bandtype='B', combine='scan', solint='inf',
    #              gaintable=[tab['Df_tab'],tab['K_tab'],tab['Ginit_tab']], interp=['linear,linearflag'], refant=ref_ant)
    #os.system(f'{ragavi_command} --table {tab['B_tab']} --plotname ./PLOTS/cal_BP-B-ccFINAL.png >> ragavi.log')
    # plotms(vis=tab['B_tab'], coloraxis='antenna1', xaxis='freq', yaxis='amp')
    # plotms(vis=tab['B_tab'], coloraxis='antenna1', xaxis='freq', yaxis='phase')
    #os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-amp-FINAL.png' {ms_cal_file} >> shadems.log")
    #os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:phase --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-ph-FINAL.png' {ms_cal_file} >> shadems.log")

    ############################################################################
    # Bootrap secondary calibrator
    log.info('Bootstrapping secondary calibrator...')
    casa.flagmanager(vis = ms_cal_file, mode = 'save', versionname = f'PreSec')

    for cc in range(3):
        log.info(f'Calibration cycle {cc+1}')
        casa.gaincal(vis=ms_cal_file, caltable=tab['Ksec_tab'], field=PhaseCal, gaintype='K', \
                    gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab']], interp=['linear,linearflag', 'linear,linearflag'], refant=ref_ant)
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab['Ksec_tab']} --plotname ./PLOTS/cal_Sec-K-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['Ksec_tab'], coloraxis='antenna1', xaxis='time', yaxis='delay', xconnector='line')
        casa.gaincal(vis=ms_cal_file, caltable=tab['Gpsec_tab'], field=PhaseCal, gaintype='G', calmode='p', \
                    gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Ksec_tab']], interp=['linear,linearflag', 'linear,linearflag'], refant=ref_ant)
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab['Gpsec_tab']} --yaxis phase --plotname ./PLOTS/cal_Sec-Gp-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['Gpsec_tab'], coloraxis='antenna1', xaxis='time', yaxis='phase', xconnector='line')
        casa.gaincal(vis=ms_cal_file, caltable=tab['Tsec_tab'], field=PhaseCal, gaintype='T', calmode='a', solnorm=True, \
                    gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Ksec_tab'], tab['Gpsec_tab']], interp=['linear,linearflag', 'linear,linearflag'], refant=ref_ant) # scalar as it can be polarised
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab['Tsec_tab']} --yaxis amplitude --plotname ./PLOTS/cal_Sec-Ta-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['Tsec_tab'], coloraxis='antenna1', xaxis='time', yaxis='amp', xconnector='line')

        if cc == 0:
            # image the secondary and selfcal to improve the local model
            casa.applycal(vis=ms_cal_file, field=PhaseCal, parang=True, flagbackup=False, interp=['linear,linearflag','linear,linearflag'], \
                        gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Ksec_tab'], tab['Gpsec_tab'], tab['Tsec_tab']])
            lib_runcode.run_command(f'{lib_runcode.wsclean_command} -name IMG/{PhaseCal}-selfcal -reorder -parallel-deconvolution 1024 -parallel-gridding 64 \
                    -update-model-required -weight briggs -0.2 -size 8000 8000 \
                    -scale {pixelscale}arcsec -channels-out 12 -pol I -data-column CORRECTED_DATA -niter 1000000 -mgain 0.8 -join-channels \
                    -fit-spectral-pol 3 -deconvolution-channels 3 -auto-mask 5 -auto-threshold 3 -field {PhaseCal_id} {ms_cal_file}', cfg, logname='wsclean')

        if cc == 1:
            # flagging on residuals
            casa.applycal(vis=ms_cal_file, field=PhaseCal, parang=True, flagbackup=False, interp=['linear,linearflag','linear,linearflag'], \
                        gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Ksec_tab'], tab['Gpsec_tab'], tab['Tsec_tab']])
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {PhaseCal} --corr XX,YY --png './PLOTS/Phasecal.png' {ms_cal_file}", cfg, logname='shadems') # check if the amp are reduced and no big waves/spikes should be there
            casa.flagdata(vis=ms_cal_file, mode="rflag", field=PhaseCal, datacolumn="residual", quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False)
            casa.flagdata(vis=ms_cal_file, mode='extend', field=PhaseCal, datacolumn='residual', growtime=80, growfreq=80, flagbackup=False, growaround=True, flagnearfreq=True)
            lib_mepi.print_flags(ms_cal_file)
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {PhaseCal} --corr XX,YY --png './PLOTS/Phasecal-flag.png' {ms_cal_file}", cfg, logname='shadems') # check if the amp are reduced and no big waves/spikes should be there
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x uv -y CORRECTED_DATA:phase -c ANTENNA1 --corr XX,YY --field {PhaseCal} --png './PLOTS/Phasecal_uvphase_XXYY.png' {ms_cal_file}", cfg, logname='shadems') # phase against uv-dsitance plot for the gain calibrator. Check both: scatter of phase per baseline (should be managable ~20deg) and if all long baselines are flagged (if the case need to adjust flagging again) 

    ##############################################################################
    # Solve for polarization alignment
    log.info('Solving for polarization alignment...')
    casa.flagmanager(vis=ms_cal_file, mode='save', versionname='PrePol')

    for cc in range(2):
        log.info(f'Calibration cycle {cc+1}')
        casa.gaincal(vis=ms_cal_file, caltable=tab['Kpol_tab'], field=PolCal, gaintype='K', interp=['linear,linearflag','linear,linearflag'], \
                    gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Gpsec_tab'], tab['Tsec_tab']], refant=ref_ant, solint='32s')
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab["Kpol_tab"]} --plotname ./PLOTS/cal_Pol-K-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['Kpol_tab'], coloraxis='antenna1', xaxis='time', yaxis='delay')
        # here we can use also secT to trace slow variations in the amp
        casa.gaincal(vis=ms_cal_file, caltable=tab['Gppol_tab'], field=PolCal, gaintype='G', calmode='p', interp=['linear,linearflag','linear,linearflag'], 
                    gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Tsec_tab'], tab['Kpol_tab']], refant=ref_ant, solint='32s')
        lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab["Gppol_tab"]} --yaxis phase --plotname ./PLOTS/cal_Pol-Gp-cc{cc}.png', cfg, logname='ragavi')
        # plotms(vis=tab['Gppol_tab'], coloraxis='antenna1', xaxis='time', yaxis='phase', xconnector='line')

        casa.applycal(vis=ms_cal_file, field=PolCal, parang=True, flagbackup=False, interp=['linear,linearflag','linear,linearflag'], \
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Tsec_tab'], tab['Kpol_tab'], tab['Gppol_tab']])

        if cc == 0:
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {PolCal} --corr XY,YX --png './PLOTS/PolCal-cross-preXf.png' {ms_cal_file}", cfg, logname='shadems')
            casa.flagdata(vis=ms_cal_file, mode="rflag", field=PolCal, datacolumn="residual", quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False)
            casa.flagdata(vis=ms_cal_file, mode='extend', field=PolCal, datacolumn='residual', growtime=80, growfreq=80, flagbackup=False, growaround=True, flagnearfreq=True)
            lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {PolCal} --corr XY,YX --png './PLOTS/PolCal-cross-preXf-flag.png' {ms_cal_file}", cfg, logname='shadems')

    # Xf that is constant within a scan, but it drift slowly with time, try not combining scans
    casa.polcal(vis=ms_cal_file, caltable=tab['Xf_tab'], field=PolCal, poltype='Xf', solint='inf,10MHz', refant=ref_ant, interp=['linear,linearflag','linear,linearflag'],
    combine='', preavg=-1., gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Tsec_tab'], tab['Kpol_tab'], tab['Gppol_tab']])
    lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab["Xf_tab"]} --yaxis phase --plotname ./PLOTS/cal_Pol-Xf.png', cfg, logname='ragavi')
    # plotms(vis=tab['Xf_tab'], xaxis='freq', yaxis='phase')
    log.info('Crosscal: Correcting for phase ambiguity')
    lib_mepi.xyamb(log, xytab=tab['Xf_tab'] ,xyout=tab['Xf_tab_ambcorr'])
    lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab["Xf_tab_ambcorr"]} --yaxis phase --plotname ./PLOTS/cal_Pol-Xfambcorr.png', cfg, logname='ragavi')
    # plotms(vis=tab['Xf_tab_ambcorr'], xaxis='freq', yaxis='phase')

    log.info('Applying calibration to PolCal and test imaging...')
    # Final applycal to PolCal to check pol quality
    casa.applycal(vis=ms_cal_file, field=PolCal, parang=True, flagbackup=False, interp=['linear,linearflag','linear,linearflag', 'linear,linearflag'], \
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Xf_tab_ambcorr'], tab['Ga_tab'], tab['Tsec_tab'], tab['Kpol_tab'], tab['Gppol_tab']])
    lib_runcode.run_command(f"{lib_runcode.shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {PolCal} --corr XY,YX --png './PLOTS/PolCal-cross-postXf.png' {ms_cal_file}", cfg, logname='shadems') # check if the amp are reduced and no big waves/spikes should be there

    # test image of the polcal - no update model!
    lib_runcode.run_command(f'{lib_runcode.wsclean_command} -name IMG/{PolCal}-selfcal -reorder -parallel-deconvolution 512 -parallel-gridding 64 \
            -no-update-model-required -weight briggs -0.2 -size 1000 1000 \
            -scale {pixelscale}arcsec -channels-out 12 -pol IQUV -data-column CORRECTED_DATA -niter 1000000 -mgain 0.8 -join-channels \
            -multiscale -fit-spectral-pol 3 -auto-mask 5 -auto-threshold 3 -field {PolCal_id} {ms_cal_file}', cfg, logname='wsclean')

    #####################################################################################
    # If the secondary is polarised we should re-do its calibration including Xf
    log.info('Re-doing secondary calibrator calibration including Xf...')
    casa.gaincal(vis=ms_cal_file, caltable=tab['Ksec_tab'], field=PhaseCal, gaintype='K', refant=ref_ant, interp=['linear,linearflag','linear,linearflag', 'linear,linearflag'], \
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Xf_tab_ambcorr'], tab['Ga_tab']])
    lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab["Ksec_tab"]} --plotname ./PLOTS/cal_Sec-K-FINAL.png', cfg, logname='ragavi')
    casa.gaincal(vis=ms_cal_file, caltable=tab['Gpsec_tab'], field=PhaseCal, gaintype='G', calmode='p', refant=ref_ant, interp=['linear,linearflag','linear,linearflag', 'linear,linearflag'], \
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Xf_tab_ambcorr'], tab['Ksec_tab'], tab['Ga_tab']])
    lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab["Gpsec_tab"]} --yaxis phase --plotname ./PLOTS/cal_Sec-Gp-FINAL.png', cfg, logname='ragavi')
    # parang=True for polarised sources and Xf should also be applied, otherwise it absorbs part of the effect
    casa.gaincal(vis=ms_cal_file, caltable=tab['Tsec_tab'], field=PhaseCal, gaintype='T', calmode='a', solnorm=True, refant=ref_ant, parang=True, interp=['linear,linearflag','linear,linearflag', 'linear,linearflag'], \
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Xf_tab_ambcorr'], tab['Ga_tab'], tab['Ksec_tab'], tab['Gpsec_tab']])
    lib_runcode.run_command(f'{lib_runcode.ragavi_command} --table {tab["Tsec_tab"]} --yaxis amplitude --plotname ./PLOTS/cal_Sec-Ta-FINAL.png', cfg, logname='ragavi')    