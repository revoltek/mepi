import os, sys, ast, argparse

_REQUIRED_KEYS = [
    "ms_full",
    "path_ms",
    "path_plots",
    "path_imgs",
    "path_logs",
    "path_sols",
    "freqbin",
    "timebin",
    "pipeline_steps",
]

class _CfgSingleton(dict):
    """A dict that can only be populated once via setup_cfg()."""
    _configured = False

    def setup_cfg(self, config_file):
        if self._configured:
            raise RuntimeError("setup_cfg() called more than once")
        if not os.path.isfile(config_file):
            print(f"Config file not found: {config_file}")
            sys.exit(1)

        with open(config_file) as fh:
            for lineno, raw_line in enumerate(fh, 1):
                line = raw_line.split("#")[0].strip()
                if not line:
                    continue
                if "=" not in line:
                    print(f"{config_file}:{lineno}: skipping unrecognised line: {raw_line.rstrip()}")
                    continue
                key, _, value_str = line.partition("=")
                key = key.strip()
                value_str = value_str.strip()
                try:
                    self[key] = ast.literal_eval(value_str)
                except (ValueError, SyntaxError):
                    self[key] = value_str

        missing = [k for k in _REQUIRED_KEYS if k not in self]
        if missing:
            print(f"Missing required parameter(s) in {config_file}: {', '.join(missing)}")
            sys.exit(1)

        _CfgSingleton._configured = True
        return self


# Single instance — import and use as:
#   from mepi import lib_cfg
#   cfg = lib_cfg.cfg          # dict-like access
#   lib_cfg.cfg.setup_cfg(...)  # called once at startup
cfg = _CfgSingleton()
def setup_cfg(config_file):
    """Parse *config_file* and populate the singleton ``cfg`` dict. Call once at startup."""
    return cfg.setup_cfg(config_file)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="mepi",
        description="MeerKAT pipeline launcher",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "config",
        metavar="CONFIG_FILE",
        help="Pipeline config file (e.g. mepi.conf or mepi.config)",
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Parse config and print parameters without running the pipeline",
    )
    parser.add_argument(
        "--log-level", "-l",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )
    return parser