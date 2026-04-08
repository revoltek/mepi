#!/usr/bin/env python3
"""
mepi - MeerKAT pipeline launcher
Usage: mepi.py [options] <config_file>
"""

import os, sys, argparse, ast, importlib

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

def parse_config(config_file):
    """
    Parse a mepi config file.

    The config file uses simple Python assignment syntax, e.g.::

        invis = 'RawData/a2163-flipped.MS/'
        freqbin = 4
        PhaseTargetDic = {'J1550+0527': 'A2163'}

    Comments (``#``) and blank lines are ignored.
    Values are evaluated with ast.literal_eval so that strings, numbers,
    dicts and lists are returned as native Python objects.

    Returns
    -------
    dict
        Mapping of parameter name -> value.
    """
    if not os.path.isfile(config_file):
        print(f"Config file not found: {config_file}")
        sys.exit(1)

    cfg = {}
    with open(config_file) as fh:
        for lineno, raw_line in enumerate(fh, 1):
            # Strip trailing comments and whitespace
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
                cfg[key] = ast.literal_eval(value_str)
            except (ValueError, SyntaxError):
                # Fall back to raw string if literal_eval fails
                cfg[key] = value_str

    # Check for required keys
    missing = [k for k in _REQUIRED_KEYS if k not in cfg]
    if missing:
        print(f"Missing required parameter(s) in {config_file}: {', '.join(missing)}")
        sys.exit(1)

    return cfg


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


def main():
    mepi_dir = os.path.dirname(os.path.abspath(__file__))

    parser = build_parser()
    args = parser.parse_args()
    
    from mepi import lib_log
    lib_log.setup_logging(
        log_level=args.log_level,
        log_file='mepi.log',
        casa_log='mepi_casa.log',
    )
    log = lib_log.logger

    # Add the mepi package directory to sys.path so pipeline modules are importable
    cfg['script_dir'] = mepi_dir
    pkg_dir = os.path.join(mepi_dir, "mepi")
    if pkg_dir not in sys.path:
        sys.path.insert(0, pkg_dir)

    log.info(f"Reading config: {args.config}")
    cfg = parse_config(args.config)

    log.info("Pipeline parameters:")
    for key in _REQUIRED_KEYS:
        log.info(f"  {key} = {cfg[key]!r}")
    extra = {k: v for k, v in cfg.items() if k not in _REQUIRED_KEYS}
    for key, val in extra.items():
        log.info(f"  {key} = {val!r}  [extra]")

    if args.dry_run:
        log.info("Dry run — exiting before pipeline execution.")
        return

    # ------------------------------------------------------------------
    # Pipeline execution starts here
    # ------------------------------------------------------------------

    # create dirs if missing
    for d in [cfg["path_ms"], cfg["path_imgs"], cfg["path_plots"], cfg["path_logs"], cfg["path_sols"]]:
        os.makedirs(d, exist_ok=True)

    for step in cfg.get("pipeline_steps", []):
        log.info(f"Running pipeline step: {step}")
        module_name = f"pipeline_{step}"
        try:
            module = importlib.import_module(module_name)
            module.run(cfg)
        except ImportError as e:
            log.error(f"Failed to import module for step '{step}': {e}")
            sys.exit(1)
        except Exception as e:
            log.error(f"Error during pipeline step '{step}': {e}", exc_info=True)
            sys.exit(1)

if __name__ == "__main__":
    main()
