import numpy as np
import casatasks as casa
from casatools import table
import time 
from mepi import lib_log

log = lib_log.logger

#### Functions needed for J0408-6545 from https://skaafrica.atlassian.net/wiki/spaces/ESDKB/pages/1481408634/Flux+and+bandpass+calibration
def casa_flux_model(lnunu0, iref, *args):
    """
    Compute model:
    iref * 10**lnunu0 ** (args[0] + args[1] * lnunu0 + args[1] * lnunu0 ** 2 + args[0] * lnunu0 ** 3)
    """
    exponent = np.sum([arg * (lnunu0 ** (power))
                       for power, arg in enumerate(args)], axis=0)
    return iref * (10**lnunu0) **(exponent)

def fit_flux_model(nu, s, nu0, sigma, sref, order=5):
    """
    Fit a flux model of given order from :
    S = fluxdensity *(freq/reffreq)**(spix[0]+spix[1]*log(freq/reffreq)+..)
    Very rarely, the requested fit fails, in which case fall
    back to a lower order, iterating until zeroth order. If all
    else fails return the weighted mean of the components.
    Finally convert the fitted parameters to a
    katpoint FluxDensityModel:
    log10(S) = a + b*log10(nu) + c*log10(nu)**2 + ...
    Parameters
    ----------
    nu : np.ndarray
        Frequencies to fit in Hz
    s : np.ndarray
        Flux densities to fit in Jy
    nu0 : float
        Reference frequency in Hz
    sigma : np.ndarray
        Errors of s
    sref : float
        Initial guess for the value of s at nu0
    order : int (optional)
        The desired order of the fitted flux model (1: SI, 2: SI + Curvature ...)
    """
    from scipy.optimize import curve_fit

    init = [sref, -0.7] + [0] * (order - 1)
    lnunu0 = np.log10(nu/nu0)
    for fitorder in range(order, -1, -1):
        try:
            popt, _ = curve_fit(casa_flux_model, lnunu0, s, p0=init[:fitorder + 1], sigma=sigma)
        except RuntimeError:
            print("Fitting flux model of order %d to CC failed. Trying lower order fit." %
                           (fitorder,))
        else:
            coeffs = np.pad(popt, ((0, order - fitorder),), "constant")
            return [nu0] +  coeffs.tolist()
    # Give up and return the weighted mean
    coeffs = [np.average(s, weights=1./(sigma**2))] + [0] * order
    return [nu0]+  coeffs.tolist()

def convert_flux_model(nu=None, a=1, b=0, c=0, d=0, Reffreq=1.0e9):
    """
    Convert a flux model from the form:
    log10(S) = a + b*log10(nu) + c*log10(nu)**2 + ...
    to an ASA style flux model in the form:
    S = fluxdensity *(freq/reffreq)**(spix[0]+spix[1]*log(freq/reffreq)+..)
    Parameters
    ----------
    nu : np.ndarray
        Frequencies to fit in Hz
    a,b,c,d : float
        parameters of a log flux model.
    Reffreq : float
        Reference frequency in Hz
    returns :
    reffreq,fluxdensity,spix[0],spix[1],spix[2]
    """
    if nu is None:
        nu = np.linspace(0.9, 2, 200) * 1e9
    MHz = 1e6
    S = 10**(a + b*np.log10(nu/MHz) + c*np.log10(nu/MHz)**2 + d*np.log10(nu/MHz)**3)
    return fit_flux_model(nu, S, Reffreq, np.ones_like(nu), sref=1, order=3)

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

def print_flags(vis):
    ##############
    # Print flagging summary
    ##############
    s = casa.flagdata(vis=vis, mode='summary')
    # Print per-antenna flags on one line
    ant_flags = ', '.join([f"{ant}: {100.0*info.get('flagged', 0)/info.get('total', 0):.1f}%" 
                           for ant, info in s['antenna'].items()])
    print(f"Antenna flags: {ant_flags}")
    
    # Print per-scan flags on one line
    scan_flags = ', '.join([f"{scan}: {100.0*info.get('flagged', 0)/info.get('total', 0):.1f}%" 
                           for scan, info in sorted(s['scan'].items(), key=lambda x: int(x[0]))])
    print(f"Scan flags: {scan_flags}")
    
    # Print total flags percentage
    total_flagged = s.get('flagged', 0)
    total_points = s.get('total', 0)
    total_pct = 100.0 * total_flagged / total_points if total_points > 0 else 0
    print(f"Total flags: {total_flagged}/{total_points} ({total_pct:.2f}%)")

def ionosphere_rm(logger, pol_ms, obs_id, path, caltype='polcal'):
    from pathlib import Path
    from spinifex import h5parm_tools
    from spinifex.vis_tools import ms_tools
    import os
    """
    Calculate the ionospheric RM for the specified calibrator and create a h5param file.
    """

    ms_path = Path(pol_ms)
    ms_metadata = ms_tools.get_metadata_from_ms(ms_path)

    rms = ms_tools.get_rm_from_ms(ms_path, use_stations=ms_metadata.station_names)
    h5parm_name = f"{path}/CAL_TABLES/{obs_id}_{caltype}.h5parm"
    if os.path.exists(h5parm_name):
        logger.info(f"ionosphere_rm: Found {h5parm_name}. Deleting it.")
        os.system(f"rm -rf {h5parm_name}")
    h5parm_tools.write_rm_to_h5parm(rms=rms, h5parm_name=h5parm_name)
    logger.info(f"ionosphere_rm: Created h5parm file {h5parm_name} with ionospheric RM.")
    return h5parm_name