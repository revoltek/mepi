import os, sys, datetime

from mepi import lib_log, lib_cfg

log = lib_log.log
cfg = lib_cfg.cfg

class Skip(Exception):
    pass

class Exit(Exception):
    pass

class Walker():
    """
    An object of this class may be used to re-run a pipeline without repeating steps that were completed previously.
    Use like:
    w = Walker("filename.walker")
    with w.if_todo("stepname"):
        Do whatever...

    Adopted from https://stackoverflow.com/questions/12594148/skipping-execution-of-with-block
    """
    def __init__(self, filename='mepi.walker', pipeline=None):
        open(filename, 'a').close() # create the file if doesn't exists
        self.filename = os.path.abspath(filename)
        self.pipeline = pipeline
        self.__skip__ = False
        self.__step__ = None
        self.__inittime__ = None
        self.__globaltimeinit__ = datetime.datetime.now()

    def if_todo(self, stepname):
        """
        This is basically a way to get a context manager to accept an argument. Will return "self" as context manager
        if called as context manager.
        """
        if self.pipeline is not None:
            stepname = f"{self.pipeline}:{stepname}"
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
            log.info('\033[92m>> start >> {}\033[0m'.format(self.__step__))
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
                hours, rem = divmod(total_seconds, 3600)
                minutes, seconds = divmod(rem, 60)
                delta = f"{hours}h {minutes}m {seconds}s"
                f.write(f"{self.__step__} # {delta}\n")
            log.info('\033[92m<< done << {}\033[0m - {}'.format(self.__step__, delta))
            return  # No exception
        if issubclass(type, Skip):
            log.warning('\033[33m>> skip << {}\033[0m'.format(self.__step__))
            return True  # Suppress special SkipWithBlock exception
        if issubclass(type, Exit):
            log.error('<< exit << {}'.format(self.__step__))
            return True

    def alldone(self):
        delta = 'h '.join(str(datetime.datetime.now() - self.__globaltimeinit__).split(':')[:-1])+'m'
        log.info('Done. Total time: '+delta)
