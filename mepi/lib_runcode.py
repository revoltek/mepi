from datetime import datetime
import os, sys
import subprocess
from mepi import lib_log, lib_cfg

log = lib_log.log
cfg = lib_cfg.cfg

class Run_command():
    def __init__(self, type):

        # external commands:
        if type == 'tricolour':
            self.command = 'tricolour'
            self.logname = 'tricolour'
        elif type == 'shadems':
            self.command = 'shadems --no-lim-save'
            self.logname = 'shadems'
        elif type == 'aoflagger':
            self.command = 'aoflagger -v -j 64'
            self.logname = 'aoflagger'
        elif type == 'wsclean':
            self.command = 'wsclean -j 64'
            self.logname = 'wsclean'
        elif type == 'mask_ms':
            self.command = os.path.join(cfg['mepi_dir'], 'bin/mask_ms.py')
            self.logname = 'mask_ms'
        elif type == 'ragavi':
            self.command = 'ragavi-gains'
            self.logname = 'ragavi'
        elif type == 'crystalball':
            self.command = 'crystalball'
            self.logname = 'crystalball'

    def run(self, params, logname=None):
        """
        Run an external command, writing stdout and stderr to a log file under cfg['path_logs'].

        Parameters
        ----------
        params : str
            Parameters to run (passed to subprocess with shell=True if str).
        logname : str, optional
            Base name for the log file (without extension). Defaults to the
            first word of the command.

        Raises
        ------
        RuntimeError
            If the command returns a non-zero exit code.
        """

        params = ' '.join(params.split()) # normalize whitespace
        cmd = self.command + ' ' + params

        if logname is None:
            logname = self.logname

        logdir = cfg['path_logs']
        logfile = os.path.join(logdir, f'{logname}.log')

        shell = isinstance(cmd, str)
        log.info(f'Running: {cmd if shell else " ".join(cmd)}')

        with open(logfile, 'a') as fh:
            result = subprocess.run(cmd, shell=shell, stdout=fh, stderr=subprocess.STDOUT)

        if result.returncode != 0:
            msg = f'Command failed (exit {result.returncode}): {cmd if shell else " ".join(cmd)}'
            msg += f'\nSee log file: {logfile}'
            log.error(msg)
            raise RuntimeError(msg)

run_wsclean = Run_command('wsclean')
run_tricolour = Run_command('tricolour')
run_shadems = Run_command('shadems')
run_aoflagger = Run_command('aoflagger')
run_mask_ms = Run_command('mask_ms')
run_ragavi = Run_command('ragavi')
run_crystalball = Run_command('crystalball')

# some parsets
aoflagger_strategy1 = os.path.join(cfg['mepi_dir'], 'parsets/aoflagger_StokesI.lua')
aoflagger_strategy2 = os.path.join(cfg['mepi_dir'], 'parsets/aoflagger_StokesQUV.lua')
rfimask = os.path.join(cfg['mepi_dir'], 'parsets/meerkat.rfimask.npy') # ok for UHF and L bands