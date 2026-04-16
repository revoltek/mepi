import os, glob
import casatasks as casa
from mepi import lib_mepi, lib_ms, lib_log, lib_runcode, lib_cfg, lib_walker, lib_sol

log = lib_log.log
cfg = lib_cfg.cfg
w = lib_walker.Walker("selfcaldp3")

def run():

    losoto_parset = os.path.join(cfg['mepi_dir'], 'parsets/losoto-plot.parset')
    dp3_sol_parset = os.path.join(cfg['mepi_dir'], 'parsets/DP3-sol.parset')
    dp3_cor_parset = os.path.join(cfg['mepi_dir'], 'parsets/DP3-cor.parset')

    ms_full_file = cfg['ms_full']

    for ms_tgt_file in glob.glob(os.path.join(cfg['path_ms'], os.path.basename(ms_full_file.rstrip('/').rstrip('\\')).replace('.MS', '_tgt_*.MS').replace('.ms', '_tgt_*.ms'))):
        ms_tgt = lib_ms.MS(ms_tgt_file)
        target = ms_tgt.find_targets()[0] # todo: handle multiple targets
        w = lib_walker.Walker(f"selfcal-{target}")
        log.info(f"Found existing target MS: {ms_tgt_file} ({target})")

        if cfg['ref_ant'] is not None:
            ref_ant = cfg['ref_ant']
            log.info(f"Using user-specified reference antenna: {ref_ant}")
        else:
            ref_ant = ms_tgt.find_reference_antenna()
            log.info(f"Using automatically selected reference antenna: {ref_ant}")
        pixelscale = round(0.7 * (2.4/(ms_tgt.freq_center*1e-9)), 1) # arcsec

        with w.if_todo("parang"):
            # selfcal only on scalar amp and possibly diag phase.
            # If diag phase needed, only for stokes I and consider parang is amp rot matrix and doesn't commute
            casa.flagmanager(vis=ms_tgt_file, mode='save', versionname='PreSelfcal')
            casa.applycal(vis=ms_tgt_file, flagbackup=False, parang=True)

        # NOTE: DP3 still doesn't work on multi-scan data
        scans = []
        for scan in casa.listobs(ms_tgt_file):
            if 'scan' in scan: scans.append(int(scan.split('_')[1]))
        mss = []
        for scan in scans:
            ms_tgt_scan_file = ms_tgt_file.replace('.MS', f'-scan{scan}.MS').replace('.ms', f'-scan{scan}.ms')
            log.info(f"Splitting scan {scan} for DP3 self-calibration {ms_tgt_file} -> {ms_tgt_scan_file}")
            casa.split(vis = ms_tgt_file, outputvis = ms_tgt_scan_file, datacolumn = 'corrected', scan=str(scan))
            mss.append(ms_tgt_scan_file)

        #for ms in mss:
        #    os.system(f'wsclean  -predict -padding 1.8 -j 64 -name {imgname} -channels-out 32 {ms} >> wsclean.log')
        for i in range(10):
            print(f'Cycle {i}: imaging...')
            lib_runcode.run_wsclean.run(f'-name IMG/{target}-selfcal-c{i}  -update-model-required -pol I \
                -reorder -parallel-reordering 5 -parallel-gridding 64 -parallel-deconvolution 1024 \
                -size 2500 2500 -scale {pixelscale}arcsec -weight briggs -0.2  -minuv-l 80.0 \
                -niter 1000000 -mgain 0.6 \
                -join-channels -channels-out 32 -deconvolution-channels 6 -fit-spectral-pol 3 \
                -multiscale -multiscale-scales 1,3,9,27,90,270 \
                -auto-threshold 3 {' '.join(mss)}')
            for ms in mss:
                print(f'Cycle {i}: Solving+correcting {ms}...')
                # solve
                lib_runcode.run_dp3.run(f'{dp3_sol_parset} msin={ms} msout=. sol.h5parm={ms}/ph-{i}.h5 sol.mode=diagonalphase sol.solint=1 sol.nchan=1 sol.smoothnessconstraint=10e6')
                lib_runcode.run_dp3.run(f'{dp3_sol_parset} msin={ms} msout=. sol.h5parm={ms}/amp-{i}.h5 sol.mode=scalaramplitude sol.solint=20 sol.nchan=1 sol.smoothnessconstraint=30e6')
                # correct
                lib_runcode.run_dp3.run(f'{dp3_cor_parset} msin={ms} msout=. cor1.parmdb={ms}/ph-{i}.h5 cor2.parmdb={ms}/amp-{i}.h5')
            # losoto
            os.system(f'H5parm_collector.py -V -s sol000 -o cal-ph-{i}.h5 '+' '.join([f'{ms}/ph-{i}.h5' for ms in mss]))
            os.system(f'losoto cal-ph-{i}.h5 {losoto_parset} >> losoto.log && mv plots plots-{i}')
            os.system(f'H5parm_collector.py -V -s sol000 -o cal-amp-{i}.h5 '+' '.join([f'{ms}/amp-{i}.h5' for ms in mss]))
            os.system(f'losoto cal-amp-{i}.h5 {losoto_parset} >> losoto.log && mv plots plots-{i}')