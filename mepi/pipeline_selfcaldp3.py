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