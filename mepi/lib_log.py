import os, logging
import casatasks as casa

import warnings
from astropy.wcs import FITSFixedWarning
# Suppress FITSFixedWarning
warnings.simplefilter('ignore', FITSFixedWarning)

def setup_logging(log_level, log_file, casa_log):
    """
    Configure the root logger with a console handler and a file handler,
    and redirect CASA logs to a separate file.

    Parameters
    ----------
    log_level : str
        Console verbosity ('DEBUG', 'INFO', 'WARNING', 'ERROR').
    log_file : str
        Path for the mepi log file.
    casa_log : str
        Path for the CASA log file.

    Returns
    -------
    logging.Logger
        The 'mepi' logger.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)  # handlers filter their own level

    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(logging.Formatter(
        "%(asctime)s %(levelname)s - %(message)s", datefmt='%Y-%m-%d %H:%M:%S'))

    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s %(levelname)s - %(name)s:%(funcName)s\t%(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'))

    root_logger.addHandler(console_handler)
    root_logger.addHandler(file_handler)

    log = logging.getLogger("mepi")

    # Redirect CASA logs
    old_log_filename = casa.casalog.logfile()
    log.debug(f"Redirecting casa logs from {old_log_filename} to {casa_log}")
    casa.casalog.setlogfile(filename=casa_log)
    os.remove(old_log_filename)

    # Suppress noisy third-party loggers
    logging.getLogger("asyncio").setLevel(logging.WARNING)

log = logging.getLogger("mepi")