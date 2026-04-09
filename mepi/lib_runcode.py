from datetime import datetime
import os, sys
import subprocess
from mepi import lib_log

log = lib_log.logger

# external commands:
# tricolour_command = f'singularity run --bind $PWD -B /local/work/fdg ~/storage/tricolour.simg tricolour'
shadems_command = f'shadems --no-lim-save'
aoflagger_command = f'aoflagger -v -j 64'
wsclean_command = f'wsclean -j 64'
mask_ms_command = os.path.join(cfg['mepi_dir'], 'mask_ms.py')
ragavi_command = 'ragavi-gains'
crystalball_command = 'crystalball'

# some parsets
aoflagger_strategy1 = os.path.join(cfg['mepi_dir'], 'parsets/aoflagger_StokesI.lua')
aoflagger_strategy2 = os.path.join(cfg['mepi_dir'], 'parsets/aoflagger_StokesQUV.lua')
rfimask = os.path.join(cfg['mepi_dir'], 'parsets/meerkat.rfimask.npy') # ok for UHF and L

def run_command(cmd, cfg, logname=None):
    """
    Run an external command, writing stdout and stderr to a log file under cfg['path_logs'].

    Parameters
    ----------
    cmd : str or list
        Command to run (passed to subprocess with shell=True if str).
    cfg : dict
        Pipeline config dict; must contain 'path_logs'.
    logname : str, optional
        Base name for the log file (without extension). Defaults to the
        first word of the command.

    Raises
    ------
    RuntimeError
        If the command returns a non-zero exit code.
    """
    if logname is None:
        first = cmd if isinstance(cmd, str) else cmd[0]
        logname = os.path.basename(first.split()[0])

    logdir = cfg['path_logs']
    logfile = os.path.join(logdir, f'{logname}.log')

    shell = isinstance(cmd, str)
    log.info(f'Running: {cmd if shell else " ".join(cmd)}')

    with open(logfile, 'a') as fh:
        result = subprocess.run(cmd, shell=shell, stdout=fh, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        msg = f'Command failed (exit {result.returncode}): {cmd if shell else " ".join(cmd)}'
        log.error(msg)
        raise RuntimeError(msg)

class Walker():
    """
    An object of this class may be used to re-run a pipeline without repeating steps that were completed previously.
    Use like:
    w = Walker("filename.walker")
    with w.if_todo("stepname"):
        Do whatever...

    Adopted from https://stackoverflow.com/questions/12594148/skipping-execution-of-with-block
    """
    def __init__(self, filename):
        open(filename, 'a').close() # create the file if doesn't exists
        self.filename = os.path.abspath(filename)
        self.__skip__ = False
        self.__step__ = None
        self.__inittime__ = None
        self.__globaltimeinit__ = datetime.datetime.now()

    def if_todo(self, stepname):
        """
        This is basically a way to get a context manager to accept an argument. Will return "self" as context manager
        if called as context manager.
        """
        self.__skip__ = False
        self.__step__ = stepname
        with open(self.filename, "r") as f:
            for stepname_done in f:
                if stepname == stepname_done.split('#')[0].rstrip():
                    self.__skip__ = True
        return self

    def __enter__(self):
        """
        Skips body of with-statement if __skip__.
        This uses some kind of dirty hack that might only work in CPython.
        """
        if self.__skip__:
            sys.settrace(lambda *args, **keys: None)
            frame = sys._getframe(1)
            frame.f_trace = self.trace
        else:
            log.info(20, '>> start >> {}'.format(self.__step__))
            self.__timeinit__ = datetime.datetime.now()

    def trace(self, frame, event, arg):
        raise Skip()

    def __exit__(self, type, value, traceback):
        """
        Catch "Skip" errors, if not skipped, write to file after exited without exceptions.
        """
        if type is None:
            with open(self.filename, "a") as f:
                delta_td = datetime.datetime.now() - self.__timeinit__
                total_seconds = int(delta_td.total_seconds())
                days, rem = divmod(total_seconds, 86400)
                hours, rem = divmod(rem, 3600)
                minutes, seconds = divmod(rem, 60)
                delta = f"{days}d {hours}h {minutes}m {seconds}s"
                f.write(f"{self.__step__} # {delta}\n")
            log.info('<< done << {}'.format(self.__step__))
            return  # No exception
        if issubclass(type, Skip):
            log.warning('>> skip << {}'.format(self.__step__))
            return True  # Suppress special SkipWithBlock exception
        if issubclass(type, Exit):
            log.error('<< exit << {}'.format(self.__step__))
            return True

    def alldone(self):
        delta = 'h '.join(str(datetime.datetime.now() - self.__globaltimeinit__).split(':')[:-1])+'m'
        log.info('Done. Total time: '+delta)
