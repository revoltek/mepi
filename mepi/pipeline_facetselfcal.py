import os, glob
import casatasks as casa
from mepi import lib_mepi, lib_ms, lib_log, lib_runcode, lib_cfg, lib_walker, lib_sol

log = lib_log.log
cfg = lib_cfg.cfg
w = lib_walker.Walker("facetselfcal")

def run():

    ms_full_file = cfg['ms_full']

    for ms_tgt_file in glob.glob(os.path.join(cfg['path_ms'], os.path.basename(ms_full_file.rstrip('/').rstrip('\\')).replace('.MS', '_tgt_*.MS').replace('.ms', '_tgt_*.ms'))):
        ms_tgt = lib_ms.MS(ms_tgt_file)
        target = ms_tgt.find_targets()[0] # todo: handle multiple targets
        w = lib_walker.Walker(f"facetselfcal-{target}")
        log.info(f"Found existing target MS: {ms_tgt_file} ({target})")

        with w.if_todo("parang"):
            # selfcal only on scalar amp and possibly diag phase.
            # If diag phase needed, only for stokes I and consider parang is amp rot matrix and doesn't commute
            casa.flagmanager(vis=ms_tgt_file, mode='save', versionname='PreSelfcal')
            casa.applycal(vis=ms_tgt_file, flagbackup=False, parang=True)

        # copy the data in a new dir
        facetselfcal_dir = os.path.abspath(f'{target}_facetselfcal')
        ms_fsc_file = os.path.join(facetselfcal_dir, os.path.basename(ms_tgt_file.rstrip('/\\')))

        with w.if_todo("copy-ms"):
            os.makedirs(facetselfcal_dir, exist_ok=True)
            log.info(f"Splitting corrected data: {ms_tgt_file} -> {ms_fsc_file}")
            casa.split(vis=ms_tgt_file, outputvis=ms_fsc_file, datacolumn='corrected')

        with w.if_todo("facetselfcal"):
            imsize             = cfg.get('facetselfcal_imsize', 12000)
            channelsout        = cfg.get('facetselfcal_channelsout', 12)
            niter              = cfg.get('facetselfcal_niter', 45000)
            stop               = cfg.get('facetselfcal_stop', 10)
            parallelgridding   = cfg.get('facetselfcal_parallelgridding', 2)
            solint_list        = cfg.get('facetselfcal_solint_list', "['1min','0.5min','0.5min']")
            soltype_list       = cfg.get('facetselfcal_soltype_list', "['scalarphase','scalarphase','scalarcomplexgain']")
            soltypecycles_list = cfg.get('facetselfcal_soltypecycles_list', '[0,0,2]')
            smoothness_list    = cfg.get('facetselfcal_smoothnessconstraint_list', '[100.]')
            aoflagger_strategy = os.path.join(cfg['mepi_dir'], 'parsets/aoflagger_StokesQUV.lua')

            parms  = f'-i complex --forwidefield --noarchive'
            parms += f' --fitspectralpol=9'
            parms += f' --solint-list="{solint_list}"'
            parms += f' --soltype-list="{soltype_list}"'
            parms += f' --soltypecycles-list={soltypecycles_list}'
            parms += f' --smoothnessconstraint-list={smoothness_list}'
            parms += f' --imsize={imsize} --channelsout={channelsout}'
            parms += f' --niter={niter} --stop={stop}'
            parms += f' --multiscale --multiscale-start=0'
            parms += f' --useaoflagger --aoflagger-strategy={aoflagger_strategy}'
            parms += f' --parallelgridding={parallelgridding}'
            if cfg.get('facetselfcal_msinnchan') is not None:
                parms += f' --msinnchan={cfg["facetselfcal_msinnchan"]}'
            if cfg.get('facetselfcal_msinstartchan') is not None:
                parms += f' --msinstartchan={cfg["facetselfcal_msinstartchan"]}'
            if cfg.get('facetselfcal_avgfreqstep') is not None:
                parms += f' --avgfreqstep={cfg["facetselfcal_avgfreqstep"]}'
            parms += f' {ms_fsc_file}'

            lib_runcode.run_facetselfcal.run(parms, logname=f'facetselfcal-{target}', cwd=facetselfcal_dir)

        # run facetselfcal as an external call (use absolute path)
        # facetselfcal -i complex --forwidefield --noarchive --fitspectralpol=9 --solint-list="['1min','0.5min','0.5min']" --soltype-list="['scalarphase','scalarphase','scalarcomplexgain']" --soltypecycles-list=[0,0,2] --smoothnessconstraint-list=[100.] --imsize=12000 --channelsout=12 --niter=45000 --stop=10 --multiscale --useaoflagger --aoflagger-strategy=default_StokesQUV.lua --multiscale-start=0 --parallelgridding=2 --msinnchan=3576 --msinstartchan=80 --avgfreqstep=2 /iranet/groups/ulu/t.pasini/MeerKAT/pilot/NGC_paper/NGC3411/LBAND/polarisation/uncalibrated/backup/target_avg.ms/

