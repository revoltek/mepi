#!/usr/bin/python3
"""
best is to run this script by copypasting in ipython3
"""

import os, sys, logging, glob
import casatasks as casa
from casatools import msmetadata
from casatools import table
import numpy as np

############################################
# input variables:
invis   = 'RawData/m87sband-flipped.MS'
calms   = 'MS_Files/m87sband-cal.MS'
tgtms   = 'MS_Files/m87sband-tgt.MS'
tgtavgms   = 'MS_Files/m87sband-tgt-avg.MS'
FluxCal = 'J1939-6342' # one of the BandPassCals
BandPassCal = 'J1939-6342,J0408-6545'
PolCal = 'J1331+3030'
PhaseTargetDic = {'J1150-0023':'M87'} # PhaseCal <--> Target pairs
ref_ant = 'm063' # longest BL

invis   = 'RawData/a2163-flipped.MS/'
calms   = 'MS_Files/a2163-cal.MS'
tgtms   = 'MS_Files/a2163-tgt.MS'
tgtavgms   = 'MS_Files/a2163-tgt-avg.MS'
FluxCal = 'J1939-6342' # one of the BandPassCals
BandPassCal = 'J1939-6342' # J0408-6545
PolCal = 'J1331+3030'
PhaseTargetDic = {'J1550+0527':'A2163'} # PhaseCal <--> Target pairs
ref_ant = 'm063'

#script_dir = os.path.dirname(os.path.abspath(__file__))
script_dir = '/home/baq1889/opt/src/mepi/'
aoflagger_strategy1 = os.path.join(script_dir, 'parsets/aoflagger_StokesI.lua')
aoflagger_strategy2 = os.path.join(script_dir, 'parsets/aoflagger_StokesQUV.lua')
rfimask = os.path.join(script_dir, 'parsets/meerkat.rfimask.npy') # ok for UHF and L
losoto_parset = os.path.join(script_dir, 'parsets/losoto-plot.parset')
dp3_sol_parset = os.path.join(script_dir, 'parsets/DP3-sol.parset')
dp3_cor_parset = os.path.join(script_dir, 'parsets/DP3-cor.parset')
spw_selection = '' # channel selection - here is what we keep in the split command - '0:210~3841' range is for band=S1
freqbin = 1 # number of channel to average for the target split
timebin = '0s' # time binning for the target split
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)
from lib_mepi import *

############################################

############################################################
# external commands:
# tricolour_command = f'singularity run --bind $PWD -B /local/work/fdg ~/storage/tricolour.simg tricolour'
#shadems_command = f'SINGULARITY_TMPDIR=$PWD singularity run --bind $PWD -B /localwork/fdg /lofar/baq1889/flocs-latest.simg shadems --no-lim-save'
shadems_command = f'shadems --no-lim-save'
aoflagger_command = f'aoflagger -v -j 64'
wsclean_command = f'wsclean -j 64'
mask_ms_command = os.path.join(script_dir, 'mask_ms.py')
ragavi_command = 'ragavi-gains'
correct_parang_command = os.path.join(script_dir, 'correct_parang.py')

# Name your gain tables
tab = {'K_tab' : 'delay_bp.cal',
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
inpath = 'CASA_Tables/'

# Fix some variables
for name in tab:
    tab[name] = os.path.join(inpath, tab[name])
PhaseCal = ','.join(PhaseTargetDic.keys())
Targets = ','.join(PhaseTargetDic.values())
CalibFields = ','.join([BandPassCal, PolCal, PhaseCal])
# create dirs if missing
for d in ['IMG', 'PLOTS', 'MS_Files', inpath]:
    os.makedirs(d, exist_ok=True)

#############################
### Logs Setting up and functions
log_file = os.path.join('mepi.log')
casa_log = os.path.join('mepi_casa.log')

log_level = logging.DEBUG
logging.basicConfig(filename=log_file,
        format='%(asctime)s %(name)s:%(funcName)s\t%(message)s',
        datefmt='%Y-%m-%d %H:%M', filemode='w',
        level=log_level)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
logger = logging.getLogger(__name__)
old_log_filename = casa.casalog.logfile()
# Point the casa logsink to the new file
casa.casalog.setlogfile(filename=casa_log)
# Delete the old file
os.remove(old_log_filename)
# remove annoying warnings
logging.getLogger("asyncio").setLevel(logging.WARNING)

##############################
# Change RECEPTOR_ANGLE : DEFAULT IS -90DEG but should be fixed with the initial swap and set to 0 for the correct interpretation of the polarisation.
tb = table()
tb.open(invis+'/FEED', nomodify=False)
feed_angle = tb.getcol('RECEPTOR_ANGLE')
new_feed_angle = np.zeros(feed_angle.shape)
tb.putcol('RECEPTOR_ANGLE', new_feed_angle)
tb.close()

###########################
# Split the calibrators
logger.info('Splitting calibrators...')
if not os.path.exists(calms):
       casa.split(vis = invis, outputvis = calms, field = f"{BandPassCal},{PolCal},{PhaseCal}", datacolumn = 'data', spw = spw_selection)
       logger.info(f'Calibrators split and saved in {calms}')
       print_flags(calms)
else:
       logger.info('Calibrators have already been split previously')

# field ids for imaging
msmd = msmetadata()
msmd.open(calms)
PhaseCal_id = msmd.fieldsforname(PhaseCal)[0]
PolCal_id = msmd.fieldsforname(PolCal)[0]
central_freq = msmd.chanfreqs(0).mean()/1e9 # central freq in GHz
if central_freq < 1: band = "UHF"
elif central_freq > 2: band = "S"
else: band = "L"
msmd.close()

# imaging parameters scaled from S1 band
pixelscale = round(0.7 * (2.4/central_freq), 1) # arcsec

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

#################################
# Standard flagging for shadowing, zero-clip, and auto-correlation
casa.flagdata(vis=calms, flagbackup=False, mode='shadow')
casa.flagdata(vis=calms, flagbackup=False, mode='manual', autocorr=True)
casa.flagdata(vis=calms, flagbackup=False, mode='clip', clipzeros=True)#, clipminmax=[0.0, 100.0])
if band == "UHF": casa.flagdata(vis=calms, flagbackup=False, mode='manual', spw='*:925~945MHz, *:950~960MHz, *:1077~1090MHz') # UHF bad data
if band == "L": casa.flagdata(vis=calms, flagbackup=False, mode='manual', spw='*:856~880MHz, *:1658~1800MHz, *:1419.8~1421.3MHz') # suggested by SARAO for Lband
if band == "S": casa.flagdata(vis=calms, flagbackup=False, mode='manual', spw='0:850~900,0:1610~1660') # resonances S1 band

if band == "UHF" or band == "L": os.system(f"{mask_ms_command} --mask {rfimask} --accumulation_mode or --memory 4096 --uvrange 0~1000 --statistics {calms}")
print_flags(calms)

# Set flux density scale
for cal in set(FluxCal.split(',')+BandPassCal.split(',')+PolCal.split(',')):
    logger.info('Setting model for calibrator %s' % cal)
    if cal == 'J1939-6342':
        casa.setjy(vis = calms, field = cal, standard = 'Stevens-Reynolds 2016', usescratch = True)
    elif cal == 'J0408-6545':
        a=-0.9790; b=3.3662; c=-1.1216; d=0.0861
        reffreq,fluxdensity,spix0,spix1,spix2 =  convert_flux_model(np.linspace(0.9,2,200)*1e9,a,b,c,d)
        casa.setjy(vis = calms, field = cal, usescratch = True, standard = 'manual', \
            spix = [spix0, spix1, spix2, 0], fluxdensity = fluxdensity, reffreq = '%f Hz'%(reffreq))
    elif cal == 'J1331+3030':
        if central_freq > 1.:
            I= 14.7172
            alpha= [-0.4507, -0.1798, 0.0357]
            reffreq= '1.47GHz'
            polfrac= 0.098
            polangle = 0.575959
            rm=0.
        else:
            I = 19.27475
            alpha = [-0.42727, -0.14583]
            reffreq = '0.816GHz'
            polfrac = 0.06398396324
            polangle = 0.41598756411
            rm = 0.12
        casa.setjy(vis=calms, field=cal, usescratch = True, standard = 'manual', \
                   fluxdensity=[I,0,0,0], spix=alpha, reffreq=reffreq, polindex=polfrac, polangle=polangle, rotmeas=rm)
        # TODO: add iono corruption to the model
    else:
        print("Unknown calibrator ", cal)
        sys.exit()

################################################################
### Bandpass calibration
logger.info('Bandpass calibration...')
casa.flagmanager(vis=calms, mode='save', versionname='PreBP')

os.system(f"{aoflagger_command} -strategy {aoflagger_strategy1} -column DATA {calms} >> aoflagger.log")
print_flags(calms)
os.system(f"{aoflagger_command} -strategy {aoflagger_strategy1} -column DATA {calms} >> aoflagger.log") # twice
print_flags(calms)

for cc in range(2):
    logger.info(f'Calibration cycle {cc+1}')
    # Delay calibration (fast to track the ionosphere)
    casa.gaincal(vis=calms, field=BandPassCal, caltable=tab['K_tab'], gaintype='K', refant=ref_ant, solint='32s')
    os.system(f'{ragavi_command} --table {tab['K_tab']} --plotname ./PLOTS/cal_BP-K-cc{cc}.png >> ragavi.log')
    # plotms(vis=tab['K_tab'], coloraxis='antenna1', xaxis='time', yaxis='delay') # all plotms should be run separately in a casa session
    # Gain calibration (fast to track the ionosphere)
    casa.gaincal(vis=calms, field=BandPassCal, caltable=tab['Gp_tab'], gaintype='G', calmode='p', 
                 gaintable=[tab['K_tab']], refant=ref_ant, solint='32s')
    os.system(f'{ragavi_command} --table {tab['Gp_tab']} --yaxis phase --plotname ./PLOTS/cal_BP-Gp-cc{cc}.png >> ragavi.log')
    # plotms(vis=tab['Gp_tab'], coloraxis='antenna1', xaxis='time', yaxis='phase')
    # one can now combine the scans and use different B as diagnostics
    casa.bandpass(vis=calms, field=BandPassCal, caltable=tab['B_tab'], bandtype='B', 
                  gaintable=[tab['K_tab'],tab['Gp_tab']], combine='scan', solint='inf', refant=ref_ant)
    os.system(f'{ragavi_command} --table {tab['B_tab']} --plotname ./PLOTS/cal_BP-B-cc{cc}.png >> ragavi.log')
    # plotms(vis=tab['B_tab'], coloraxis='antenna1', xaxis='freq', yaxis='amp')
    # plotms(vis=tab['B_tab'], coloraxis='antenna1', xaxis='freq', yaxis='phase')
    casa.gaincal(vis=calms, field=BandPassCal, caltable=tab['Ga_tab'], gaintype='G', calmode='a', 
                 gaintable=[tab['B_tab'],tab['K_tab'],tab['Gp_tab']], interp=['linear,linearflag'], refant=ref_ant)
    os.system(f'{ragavi_command} --table {tab['Ga_tab']} --yaxis amplitude --plotname ./PLOTS/cal_BP-Ga-cc{cc}.png >> ragavi.log')
    # plotms(vis=tab['Ga_tab'], coloraxis='antenna1', xaxis='time', yaxis='amp', xconnector='line')

    if cc == 0:
        casa.applycal(vis=calms, interp=['linear,linearflag', 'nearest', 'nearest'], flagbackup=False, 
                      gaintable=[tab['B_tab'],tab['K_tab'],tab['Gp_tab'],tab['Ga_tab']]) # apply to all fields for better rfi flagging
        os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-amp.png' {calms} >> shadems.log")
        os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:phase --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-ph.png' {calms} >> shadems.log")
        
        os.system(f"{aoflagger_command} -strategy {aoflagger_strategy2} -column CORRECTED_DATA {calms} >> aoflagger.log")
        casa.flagdata(vis=calms, mode="rflag", datacolumn="corrected", quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False)
        casa.flagdata(vis=calms, mode='extend', datacolumn='corrected', growtime=80, growfreq=80, flagbackup=False)
        print_flags(calms)

        os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-amp-flag.png' {calms} >> shadems.log")
        os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:phase --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-ph-flag.png' {calms} >> shadems.log")
        os.system(f"{shadems_command} -x CORRECTED_DATA:phase -y CORRECTED_DATA:amp --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass_ampph-flag.png' {calms} >> shadems.log")

#############################################
### Leackage calibration
logger.info('Leackage calibration...')
casa.flagmanager(vis = calms, mode = 'save', versionname = f'PreLeak')

# plot and flagging
os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XY,YX --png './PLOTS/Bandpass-cross-preleak.png' {calms} >> shadems.log")
casa.applycal(vis=calms, interp=['linear,linearflag'], flagbackup=False, 
              field=BandPassCal, gaintable=[tab['B_tab'],tab['K_tab'],tab['Gp_tab'],tab['Ga_tab']])
casa.flagdata(vis=calms, mode="rflag", datacolumn="corrected", field=BandPassCal, quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False, correlation='XY,YX')
casa.flagdata(vis=calms, mode='extend', datacolumn="corrected", field=BandPassCal, growtime=80, growfreq=80, flagbackup=False, growaround=True, flagnearfreq=True, correlation='XY,YX')
print_flags(calms)
os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XY,YX --png './PLOTS/Bandpass-cross-preleak-flag.png' {calms} >> shadems.log")

for cc in range(2):
    logger.info(f'Calibration cycle {cc+1}')
    # Leackage on unpol calib
    casa.polcal(vis=calms,
        caltable=tab['Df_tab'],field=BandPassCal, poltype='Dflls', solint='inf', combine='scan',
        gaintable=[tab['B_tab'], tab['K_tab'], tab['Gp_tab'], tab['Ga_tab']], interp=['linear,linearflag'], refant=ref_ant)
    os.system(f'{ragavi_command} --table {tab['Df_tab']} --xaxis channel --yaxis amplitude --plotname ./PLOTS/cal_BP-Df-cc{cc}.png >> ragavi.log')
    # plotms(vis=tab['Df_tab'], xaxis='frequency', yaxis='amplitude', coloraxis='antenna1')
    casa.flagdata(vis=tab['Df_tab'], mode='tfcrop', datacolumn="CPARAM", quackinterval=0.0, ntime="60s", combinescans=True, timecutoff=4.0, freqcutoff=3.0, usewindowstats="both", flagbackup=False)
    os.system(f'{ragavi_command} --table {tab['Df_tab']} --xaxis channel --yaxis amplitude --plotname ./PLOTS/cal_BP-Df-cc{cc}-flag.png >> ragavi.log')
    # plotms(vis=tab['Df_tab'], xaxis='frequency', yaxis='amplitude', coloraxis='antenna1')

    if cc == 0:
        casa.applycal(vis=calms,field=BandPassCal, interp=['linear,linearflag'], flagbackup=False, 
                      gaintable=[tab['B_tab'],tab['K_tab'],tab['Gp_tab'],tab['Ga_tab'],tab['Df_tab']])
        os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XY,YX --png './PLOTS/Bandpass-cross-postleak.png' {calms} >> shadems.log") # check if the amp are reduced and no big waves/spikes should be there
        casa.flagdata(vis=calms, mode="rflag", datacolumn="residual", field=BandPassCal, quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False, correlation='XY,YX')
        casa.flagdata(vis=calms, mode='extend', datacolumn="residual", field=BandPassCal, growtime=80, growfreq=80, flagbackup=False, growaround=True, flagnearfreq=True, correlation='XY,YX')
        print_flags(calms)
        os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XY,YX --png './PLOTS/Bandpass-cross-postleak-flag.png' {calms} >> shadems.log") # check if the amp are reduced and no big waves/spikes should be there
        os.system(f"{shadems_command} -x ANTENNA1 -y CORRECTED_DATA:real --field {BandPassCal} --corr XY,YX --png './PLOTS/Bandpass-cross-ant.png' {calms} >> shadems.log") #important check for chosing the reference antenna, make sure that no antenna with extreme leakage is chosen  

# re-do BP to use better flags
logger.info(f'Re-doing bandpass calibration with better flags...')
casa.bandpass(vis=calms, field=BandPassCal, caltable=tab['B_tab'], bandtype='B', combine='scan', solint='inf',
              gaintable=[tab['Df_tab'],tab['K_tab'],tab['Gp_tab']], interp=['linear,linearflag'], refant=ref_ant)
os.system(f'{ragavi_command} --table {tab['B_tab']} --plotname ./PLOTS/cal_BP-B-ccFINAL.png >> ragavi.log')
# plotms(vis=tab['B_tab'], coloraxis='antenna1', xaxis='freq', yaxis='amp')
# plotms(vis=tab['B_tab'], coloraxis='antenna1', xaxis='freq', yaxis='phase')
os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-amp-FINAL.png' {calms} >> shadems.log")
os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:phase --field {BandPassCal} --corr XX,YY --png './PLOTS/Bandpass-ph-FINAL.png' {calms} >> shadems.log")

############################################################################
# Bootrap secondary calibrator
logger.info('Bootstrapping secondary calibrator...')
casa.flagmanager(vis = calms, mode = 'save', versionname = f'PreSec')

for cc in range(3):
    logger.info(f'Calibration cycle {cc+1}')
    casa.gaincal(vis=calms, caltable=tab['Ksec_tab'], field=PhaseCal, gaintype='K', \
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab']], interp=['linear,linearflag', 'linear,linearflag'], refant=ref_ant)
    os.system(f'{ragavi_command} --table {tab['Ksec_tab']} --plotname ./PLOTS/cal_Sec-K-cc{cc}.png >> ragavi.log')
    # plotms(vis=tab['Ksec_tab'], coloraxis='antenna1', xaxis='time', yaxis='delay', xconnector='line')
    casa.gaincal(vis=calms, caltable=tab['Gpsec_tab'], field=PhaseCal, gaintype='G', calmode='p', \
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Ksec_tab']], interp=['linear,linearflag', 'linear,linearflag'], refant=ref_ant)
    os.system(f'{ragavi_command} --table {tab['Gpsec_tab']} --yaxis phase --plotname ./PLOTS/cal_Sec-Gp-cc{cc}.png >> ragavi.log')
    # plotms(vis=tab['Gpsec_tab'], coloraxis='antenna1', xaxis='time', yaxis='phase', xconnector='line')
    casa.gaincal(vis=calms, caltable=tab['Tsec_tab'], field=PhaseCal, gaintype='T', calmode='a', solnorm=True, \
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Ksec_tab'], tab['Gpsec_tab']], interp=['linear,linearflag', 'linear,linearflag'], refant=ref_ant) # scalar as it can be polarised
    os.system(f'{ragavi_command} --table {tab['Tsec_tab']} --yaxis amplitude --plotname ./PLOTS/cal_Sec-Ta-cc{cc}.png >> ragavi.log')
    # plotms(vis=tab['Tsec_tab'], coloraxis='antenna1', xaxis='time', yaxis='amp', xconnector='line')

    if cc == 0:
        # image the secondary and selfcal to improve the local model
        casa.applycal(vis=calms, field=PhaseCal, parang=True, flagbackup=False, interp=['linear,linearflag','linear,linearflag'], \
                    gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Ksec_tab'], tab['Gpsec_tab'], tab['Tsec_tab']])
        os.system(f'{wsclean_command} -name IMG/{PhaseCal}-selfcal -reorder -parallel-deconvolution 1024 -parallel-gridding 64 \
                -update-model-required -weight briggs -0.2 -size 8000 8000 \
                -scale {pixelscale}arcsec -channels-out 12 -pol I -data-column CORRECTED_DATA -niter 1000000 -mgain 0.8 -join-channels \
                -fit-spectral-pol 3 -deconvolution-channels 3 -auto-mask 5 -auto-threshold 3 -field {PhaseCal_id} {calms} > wsclean_{PhaseCal}-selfcal.log')

    if cc == 1:
        # flagging on residuals
        casa.applycal(vis=calms, field=PhaseCal, parang=True, flagbackup=False, interp=['linear,linearflag','linear,linearflag'], \
                    gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Ksec_tab'], tab['Gpsec_tab'], tab['Tsec_tab']])
        os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {PhaseCal} --corr XX,YY --png './PLOTS/Phasecal.png' {calms} >> shadems.log") # check if the amp are reduced and no big waves/spikes should be there
        casa.flagdata(vis=calms, mode="rflag", field=PhaseCal, datacolumn="residual", quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False)
        casa.flagdata(vis=calms, mode='extend', field=PhaseCal, datacolumn='residual', growtime=80, growfreq=80, flagbackup=False, growaround=True, flagnearfreq=True)
        print_flags(calms)
        os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {PhaseCal} --corr XX,YY --png './PLOTS/Phasecal-flag.png' {calms} >> shadems.log") # check if the amp are reduced and no big waves/spikes should be there
        os.system(f"{shadems_command} -x uv -y CORRECTED_DATA:phase -c ANTENNA1 --corr XX,YY --field {PhaseCal} --png './PLOTS/Phasecal_uvphase_XXYY.png' {calms} >> shadems.log") # phase against uv-dsitance plot for the gain calibrator. Check both: scatter of phase per baseline (should be managable ~20deg) and if all long baselines are flagged (if the case need to adjust flagging again) 

##############################################################################
# Solve for polarization alignment
logger.info('Solving for polarization alignment...')
casa.flagmanager(vis=calms, mode='save', versionname='PrePol')

for cc in range(2):
    logger.info(f'Calibration cycle {cc+1}')
    casa.gaincal(vis=calms, caltable=tab['Kpol_tab'], field=PolCal, gaintype='K', interp=['linear,linearflag','linear,linearflag'], \
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Ga_tab'], tab['Gpsec_tab'], tab['Tsec_tab']], refant=ref_ant, solint='32s')
    os.system(f'{ragavi_command} --table {tab["Kpol_tab"]} --plotname ./PLOTS/cal_Pol-K-cc{cc}.png >> ragavi.log')
    # plotms(vis=tab['Kpol_tab'], coloraxis='antenna1', xaxis='time', yaxis='delay')
    # here we can use also secT to trace slow variations in the amp
    casa.gaincal(vis=calms, caltable=tab['Gppol_tab'], field=PolCal, gaintype='G', calmode='p', interp=['linear,linearflag','linear,linearflag'], 
                gaintable=[tab['B_tab'], tab['Df_tab'], tab['Kpol_tab'], tab['Ga_tab'], tab['Tsec_tab']], refant=ref_ant, solint='32s')
    os.system(f'{ragavi_command} --table {tab["Gppol_tab"]} --yaxis phase --plotname ./PLOTS/cal_Pol-Gp-cc{cc}.png >> ragavi.log')
    # plotms(vis=tab['Gppol_tab'], coloraxis='antenna1', xaxis='time', yaxis='phase', xconnector='line')

    casa.applycal(vis=calms, field=PolCal, parang=True, flagbackup=False, interp=['linear,linearflag','linear,linearflag'], \
              gaintable=[tab['B_tab'], tab['Df_tab'], tab['Kpol_tab'], tab['Ga_tab'], tab['Tsec_tab'], tab['Gppol_tab']])

    if cc == 0:
        os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {PolCal} --corr XY,YX --png './PLOTS/PolCal-cross-preXf.png' {calms} >> shadems.log")
        casa.flagdata(vis=calms, mode="rflag", field=PolCal, datacolumn="corrected", quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False)
        casa.flagdata(vis=calms, mode='extend', field=PolCal, datacolumn='corrected', growtime=80, growfreq=80, flagbackup=False, growaround=True, flagnearfreq=True)
        os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {PolCal} --corr XY,YX --png './PLOTS/PolCal-cross-preXf-flag.png' {calms} >> shadems.log")

# Xf that is constant within a scan, but it drift slowly with time, try not combining scans
casa.polcal(vis=calms, caltable=tab['Xf_tab'], field=PolCal, poltype='Xf', solint='inf,10MHz', refant=ref_ant, interp=['linear,linearflag','linear,linearflag'],
   combine='', preavg=-1., gaintable=[tab['B_tab'], tab['Df_tab'], tab['Kpol_tab'], tab['Ga_tab'], tab['Tsec_tab'], tab['Gppol_tab']])
os.system(f'{ragavi_command} --table {tab["Xf_tab"]} --yaxis phase --plotname ./PLOTS/cal_Pol-Xf.png >> ragavi.log')
# plotms(vis=tab['Xf_tab'], xaxis='freq', yaxis='phase')
logger.info('Crosscal: Correcting for phase ambiguity')
xyamb(logger, xytab=tab['Xf_tab'] ,xyout=tab['Xf_tab_ambcorr'])
os.system(f'{ragavi_command} --table {tab["Xf_tab_ambcorr"]} --yaxis phase --plotname ./PLOTS/cal_Pol-Xfambcorr.png >> ragavi.log')
# plotms(vis=tab['Xf_tab_ambcorr'], xaxis='freq', yaxis='phase')

logger.info('Applying calibration to PolCal and test imaging...')
# Final applycal to PolCal to check pol quality
casa.applycal(vis=calms, field=PolCal, parang=True, flagbackup=False, interp=['linear,linearflag','linear,linearflag', 'linear,linearflag'], \
              gaintable=[tab['B_tab'], tab['Df_tab'], tab['Xf_tab_ambcorr'], tab['Kpol_tab'], tab['Ga_tab'], tab['Tsec_tab'], tab['Gppol_tab']])
os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --field {PolCal} --corr XY,YX --png './PLOTS/PolCal-cross-postXf.png' {calms} >> shadems.log") # check if the amp are reduced and no big waves/spikes should be there

# test image of the polcal - no update model!
os.system(f'{wsclean_command} -name IMG/{PolCal}-selfcal -reorder -parallel-deconvolution 512 -parallel-gridding 64 \
          -no-update-model-required -weight briggs -0.2 -size 1000 1000 \
          -scale {pixelscale}arcsec -channels-out 12 -pol IQUV -data-column CORRECTED_DATA -niter 1000000 -mgain 0.8 -join-channels \
          -multiscale -fit-spectral-pol 3 -auto-mask 5 -auto-threshold 3 -field {PolCal_id} {calms} > wsclean_{PolCal}-selfcal.log')

#####################################################################################
# If the secondary is polarised we should re-do its calibration including Xf
logger.info('Re-doing secondary calibrator calibration including Xf...')
casa.gaincal(vis=calms, caltable=tab['Ksec_tab'], field=PhaseCal, gaintype='K', refant=ref_ant, interp=['linear,linearflag','linear,linearflag', 'linear,linearflag'], \
             gaintable=[tab['B_tab'], tab['Df_tab'], tab['Xf_tab_ambcorr'], tab['Ga_tab']])
os.system(f'{ragavi_command} --table {tab['Ksec_tab']} --plotname ./PLOTS/cal_Sec-K-FINAL.png >> ragavi.log')
casa.gaincal(vis=calms, caltable=tab['Gpsec_tab'], field=PhaseCal, gaintype='G', calmode='p', refant=ref_ant, interp=['linear,linearflag','linear,linearflag', 'linear,linearflag'], \
             gaintable=[tab['B_tab'], tab['Df_tab'], tab['Xf_tab_ambcorr'], tab['Ksec_tab'], tab['Ga_tab']])
os.system(f'{ragavi_command} --table {tab['Gpsec_tab']} --yaxis phase --plotname ./PLOTS/cal_Sec-Gp-FINAL.png >> ragavi.log')
# parang=True for polarised sources and Xf should also be applied, otherwise it absorbs part of the effect
casa.gaincal(vis=calms, caltable=tab['Tsec_tab'], field=PhaseCal, gaintype='T', calmode='a', solnorm=True, refant=ref_ant, parang=True, interp=['linear,linearflag','linear,linearflag', 'linear,linearflag'], \
             gaintable=[tab['B_tab'], tab['Df_tab'], tab['Xf_tab_ambcorr'], tab['Ksec_tab'], tab['Ga_tab'], tab['Gpsec_tab']])
os.system(f'{ragavi_command} --table {tab['Tsec_tab']} --yaxis amplitude --plotname ./PLOTS/cal_Sec-Ta-FINAL.png >> ragavi.log')

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# Target

# Split the target
logger.info('Splitting target...')
if not os.path.exists(tgtms):
       casa.split(vis = invis, outputvis = tgtms, field = f"{Targets}", datacolumn = 'data', spw = spw_selection)
       logger.info(f'Target split and saved in {tgtms}')
else:
       logger.info('Target has already been split previously')
print_flags(tgtms)

casa.applycal(vis=tgtms, parang=False, flagbackup=False, interp=['linear,linearflag', 'linear,linearflag', 'linear,linearflag'], \
              gaintable=[tab['B_tab'], tab['Df_tab'], tab['Xf_tab_ambcorr'], tab['Ksec_tab'], tab['Ga_tab'], tab['Gpsec_tab'], tab['Tsec_tab']])
print_flags(tgtms)

# Standard flagging for shadowing, zero-clip, and auto-correlation
casa.flagdata(vis=tgtms, flagbackup=False, mode='shadow')
casa.flagdata(vis=tgtms, flagbackup=False, mode='manual', autocorr=True)
casa.flagdata(vis=tgtms, flagbackup=False, mode='clip', clipzeros=True, clipminmax=[0.0, 1000.0]) # high for virgo A, 100 is ok for others
if central_freq < 2: os.system(f"{mask_ms_command} --mask {rfimask} --accumulation_mode or --memory 4096 --uvrange 0~1000 --statistics {tgtms}")
print_flags(tgtms)
casa.flagmanager(vis=tgtms, mode='save', versionname='PreAoflagger')
os.system(f"{aoflagger_command} -strategy {aoflagger_strategy1} -column CORRECTED_DATA {tgtms} >> aoflagger.log")
os.system(f"{aoflagger_command} -strategy {aoflagger_strategy1} -column CORRECTED_DATA {tgtms} >> aoflagger.log") # twice
print_flags(tgtms)

# Split the target averaged in freq and time
logger.info('Splitting target avg...')
if not os.path.exists(tgtavgms):
       casa.split(vis = tgtms, outputvis = tgtavgms, field = f"{Targets}", datacolumn = 'corrected',
                  width=freqbin, timebin=timebin)
       logger.info(f'Target split and saved in {tgtavgms}')
else:
       logger.info('Target has already been split previously')
print_flags(tgtavgms)

##########################################################################################
# selfcal only on scalar amp and possibly diag phase.
# If diag phase needed, only for stokes I and consider parang is amp rot matrix and doesn't commute
casa.flagmanager(vis=tgtavgms, mode='save', versionname='PreSelfcal')
casa.applycal(vis=tgtavgms, flagbackup=False, parang=True)
for cc in range(30):
    # ok for m87 sband
    os.system(f'{wsclean_command} -name IMG/{Targets}-selfcal-c{cc}  -update-model-required -pol I \
          -reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 \
          -size 2500 2500 -scale {pixelscale}arcsec -weight briggs -0.2  \
          -niter 1000000 -mgain 0.7 \
          -join-channels -channels-out 32 -deconvolution-channels 6 -fit-spectral-pol 3 \
          -multiscale -multiscale-scales 1,4,8,16,32,64,128,256 \
          -auto-threshold 3 \
          {tgtavgms} > wsclean_{Targets}-selfcal.log')
    
    os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --corr XX,YY --png './PLOTS/Tgt-c{cc}.png' {tgtavgms} >> shadems.log")
    casa.flagdata(vis=tgtavgms, mode="rflag", datacolumn="residual", quackinterval=0.0, timecutoff=4.0, freqcutoff=3.0, extendpols=False, flagbackup=False, outfile="",overwrite=True, extendflags=False)
    casa.flagdata(vis=tgtavgms, mode='extend', datacolumn='residual', growtime=80, growfreq=80, flagbackup=False, growaround=True, flagnearfreq=True)
    os.system(f"{shadems_command} -x FREQ -y CORRECTED_DATA:amp --corr XX,YY --png './PLOTS/Tgt-c{cc}-flag.png' {tgtavgms} >> shadems.log")

    casa.gaincal(vis=tgtavgms, caltable='CASA_Tables/selfcal%02i.K' %cc, gaintype='K', solint='32s', refant=ref_ant, parang=False)
    os.system(f'{ragavi_command} --table CASA_Tables/selfcal{cc:02d}.K --plotname ./PLOTS/target-K-i{cc:02d}.png >> ragavi.log')
    # plotms(vis='CASA_Tables/selfcal%02i.K' %cc, coloraxis='antenna1', xaxis='time', yaxis='delay')
    casa.gaincal(vis=tgtavgms, caltable='CASA_Tables/selfcal%02i.Gp' %cc,  gaintype='G', calmode='p', solint='8s', refant=ref_ant, parang=False,
                 gaintable=['CASA_Tables/selfcal%02i.K' %cc])
    os.system(f'{ragavi_command} --table CASA_Tables/selfcal{cc:02d}.Gp --plotname ./PLOTS/target-Gp-i{cc:02d}.png >> ragavi.log')
    # plotms(vis='CASA_Tables/selfcal%02i.Gp' %cc, coloraxis='antenna1', xaxis='time', yaxis='phase', xconnector='line')
    casa.gaincal(vis=tgtavgms, caltable='CASA_Tables/selfcal%02i.Ga' %cc, gaintype='T', calmode='a', solint='80s', refant=ref_ant, solnorm=True, parang=True,
                 gaintable=['CASA_Tables/selfcal%02i.K' %cc, 'CASA_Tables/selfcal%02i.Gp' %cc])
    os.system(f'{ragavi_command} --table CASA_Tables/selfcal{cc:02d}.Ga --plotname ./PLOTS/target-Ga-i{cc:02d}.png >> ragavi.log')
    # plotms(vis='CASA_Tables/selfcal%02i.Ga' %cc, coloraxis='antenna1', xaxis='time', yaxis='amp', xconnector='line')
    #casa.bandpass(vis=tgtavgms, caltable='selfcal%02i.B' %cc, combine='', solint='300s', gaintable=['selfcal%02i.G' %cc, 'selfcal%02i.K' %cc], refant='m002', parang=False)
    casa.applycal(vis=tgtavgms, flagbackup=False, parang=True,
                  gaintable=['CASA_Tables/selfcal%02i.K' %cc, 'CASA_Tables/selfcal%02i.Gp' %cc, 'CASA_Tables/selfcal%02i.Ga' %cc])
    print_flags(tgtavgms)

# pol cleaning - possible problem with -squared-channel-joining when using -multiscale
#os.system(f'{wsclean_command} -name IMG/{Targets}-selfcal-pol -update-model-required -pol IQUV '
#          f'-reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 -baseline-averaging 12 '
#          f'-size 2500 2500 -scale {pixelscale}arcsec -weight briggs -0.2 -minuv-l 80.0 '
#          f'-niter 1000000 -mgain 0.7 '
#         f'-join-channels -channels-out 32 -deconvolution-channels 6 -fit-spectral-pol 3 -squared-channel-joining '
#          f'-multiscale -multiscale-scales 1,4,8,16,32,64,128,256 '
#          f'-auto-threshold 3 -fits-mask m87-07asec-2500.fits '
#          f'{tgtavgms} > wsclean_{Targets}-selfcal.log')

# wsclean with rm
restoring_beam = 6.0 # arcsec - this is ok for S1 band
os.system(f'{wsclean_command} -name img/m87-rm -no-update-model-required -pol QU '
          f'-reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 -baseline-averaging 12 '
          f'-size 1500 1500 -scale 2arcsec -weight briggs -0.5 -minuv-l 80.0 -beam-size {restoring_beam} -taper-gaussian {restoring_beam}arcsec '
          f'-niter 25000 -mgain 0.75 -nmiter 12 '
          f'-join-channels -channels-out 125 -join-polarizations -squared-channel-joining -fit-rm '
          f'{tgtavgms} > wsclean_{Targets}-selfcal.log')
# relative I stokes for fractional Pol
os.system(f'{wsclean_command} -name img/m87-rm -no-update-model-required -pol I '
          f'-reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 -baseline-averaging 12 '
          f'-size 1500 1500 -scale 2arcsec -weight briggs -0.5 -minuv-l 80.0 -beam-size {restoring_beam} -taper-gaussian {restoring_beam}arcsec '
          f'-niter 25000 -mgain 0.75 -nmiter 12 '
          f'-join-channels -channels-out 125 '
          f'{tgtavgms} > wsclean_{Targets}-selfcal.log')

sys.exit(0)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# Selfcal with DP3 - as an alternative method
# NOTE: DP3 still not work on multi-scan obs
casa.applycal(vis=tgtavgms, flagbackup=False, parang=True) # apply parang and split corrected data
scans = []
for scan in casa.listobs(tgtavgms):
    if 'scan' in scan: scans.append(int(scan.split('_')[1]))
for scan in scans:
    casa.split(vis = tgtavgms, outputvis = f'{tgtavgms.replace(".MS", f"-scan{scan}.MS")}', field = f"{Targets}", datacolumn = 'corrected', scan=str(scan))
mss = sorted(glob.glob(f'{tgtavgms.replace(".MS", "")}-scan*.MS'))
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
