import os, glob
import casatasks as casa
from mepi import lib_mepi, lib_ms, lib_log, lib_runcode, lib_cfg, lib_walker, lib_sol

log = lib_log.log
cfg = lib_cfg.cfg
w = lib_walker.Walker("splittgt")

def run():

    ms_full_file = cfg['ms_full']
    ms_tgt_file = os.path.join(cfg['path_ms'], os.path.basename(ms_full_file.rstrip('/').rstrip('\\')).replace('.MS', '_tgt.MS').replace('.ms', '_tgt.ms'))
    ms_full = lib_ms.MS(ms_full_file)
    # solution tables
    tab = lib_sol.tab

    Targets = ','.join(ms_full.find_targets())

    #########################################################################
    # Split the target
    if not os.path.exists(ms_tgt_file):
        log.info('Splitting target...')
        spw_selection = ms_full.get_spw_noedges() # select all channels but the first and last 5%
        casa.split(vis = ms_full_file, outputvis = ms_tgt_file, field = f"{Targets}", datacolumn = 'data', spw = spw_selection)
        log.info(f'Target split and saved in {ms_tgt_file}')
    else:
        log.info('Target has already been split previously')
    lib_mepi.print_flags(ms_tgt_file)
    ms_tgt = lib_ms.MS(ms_tgt_file)
    ms_tgt.swap_feeds()

    if cfg['ref_ant'] is not None:
        ref_ant = cfg['ref_ant']
        log.info(f"Using user-specified reference antenna: {ref_ant}")
    else:
        ref_ant = ms_tgt.find_reference_antenna()
        log.info(f"Using automatically selected reference antenna: {ref_ant}")

    #############################################################################
    # First apply the calibration from crosscal to the target, flag, and split the target for self-calibration
    with w.if_todo("applycal"):
        log.info('Applying calibration to target...')
        casa.applycal(vis=ms_tgt_file, parang=False, flagbackup=False, interp=['linear,linearflag', 'linear,linearflag', 'linear,linearflag'], \
                    gaintable=[tab['B'], tab['Df'], tab['Xf_ambcorr'], tab['Ga'], tab['Ksec'], tab['Gpsec'], tab['Tsec']])
        lib_mepi.print_flags(ms_tgt_file)

    with w.if_todo("flag"):
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
    for target in ms_tgt.find_targets():
        log.info(f"Target {target}")
        ms_tgtavg_file = ms_tgt_file.replace('.MS', f"_{target}.MS").replace('.ms', f"_{target}.ms")
        if not os.path.exists(ms_tgtavg_file):
            casa.split(vis = ms_tgt_file, outputvis = ms_tgtavg_file, field = f"{target}", datacolumn = 'corrected',
                    width=cfg['freqbin'], timebin=cfg['timebin'])
            log.info(f'Target split and saved in {ms_tgtavg_file}')
            lib_mepi.print_flags(ms_tgtavg_file)
        else:
            log.info(f'Target {target} has already been split previously')
