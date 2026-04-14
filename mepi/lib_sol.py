import os, time
from casatools import table
import numpy as np
from mepi import lib_log, lib_cfg

log = lib_log.log
cfg = lib_cfg.cfg

# Name your gain tables
tab = {'Ginit' : 'gain_init_bp.cal',
        'K' : 'delay_bp.cal',
        'B' : 'bandpass.cal',
        'Gp' : 'gain_p_bp.cal',
        'Ga' : 'gain_a_bp.cal',
        'Ksec' : 'delay_sec.cal',
        'Tsec' : 'T_sec.cal',
        'Gpsec' : 'gain_p_sec.cal',
        'Kpol' : 'delay_pol.cal',
        'Gppol' : 'gain_p_pol.cal',
        # Pol cal tables
        'Xf' : 'Xf.cal',
        'Xf_ambcorr' : 'Xf_ambcorr.cal',
        'Df' : 'Df.cal'}

# Fix some variables
for name in tab:
    tab[name] = os.path.join(cfg['path_sols'], tab[name])


def xyamb(logger, xytab, xyout=''):
    """
    Resolve the 180-degree cross-hand phase ambiguity in a CASA calibration table.
    Calculates the mean phase and shifts every point deviating more then 90 degrees from the mean phase by 180 degrees.

    Parameters:
    xytab : str
        Path to the input calibration table.
    xyout : str, optional
        Path to the output calibration table. If not specified, the input table is modified in place.
    """
    tb=table()

    if xyout == '':
        xyout = xytab
    if xyout != xytab:
        tb.open(xytab)
        tb.copy(xyout)
        tb.clearlocks()
        tb.close()

    tb.open(xyout, nomodify=False)

    spw_ids = np.unique(tb.getcol('SPECTRAL_WINDOW_ID'))

    for spw in spw_ids:
        st = tb.query('SPECTRAL_WINDOW_ID=='+str(spw))
        if st.nrows() > 0:
            c = st.getcol('CPARAM')
            fl = st.getcol('FLAG')
            num_channels = c.shape[1]
            num_scans = c.shape[2]
            flipped_channels = 0
            avg_phase = np.angle(np.mean(c[0, :, :][~fl[0,:,:]]), True)
            logger.info('xyamb: Average phase = '+str(avg_phase))
            for ch in range(num_channels):
                for scan in range(num_scans):
                    valid_data = c[0, ch, scan][~fl[0, ch, scan]]
                    if np.size(valid_data) > 0:
                        xyph0 = np.angle(np.mean(valid_data), True)
                        phase_diff = np.abs(xyph0 - avg_phase)
                        if phase_diff >= 100.0:
                            flipped_channels += 1
                            c[0, ch, scan] *= -1.0
                            st.putcol('CPARAM', c)
            logger.info('xyamb: Flipped '+str(flipped_channels)+' channels in SPW '+str(spw))
            st.close()
            time.sleep(1)
    
    tb.clearlocks()
    tb.flush()
    tb.close()