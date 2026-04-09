#!/usr/bin/env python3
"""
mepi - MeerKAT pipeline launcher
Usage: mepi.py [options] <config_file>
"""

import os, sys, importlib, glob, shutil

try:
    from mepi import lib_cfg, lib_log, lib_walker
except ImportError as e:
    # Add the mepi package directory to sys.path so pipeline modules are importable
    mepi_dir = os.path.dirname(os.path.abspath(__file__))
    pkg_dir = os.path.join(mepi_dir, "mepi")
    if pkg_dir not in sys.path:
        sys.path.insert(0, pkg_dir)
    from mepi import lib_cfg, lib_log, lib_walker

def main():
    
    parser = lib_cfg.build_parser()
    args = parser.parse_args()
    
    lib_log.setup_logging(
        log_level=args.log_level,
        log_file='mepi.log',
        casa_log='mepi_casa.log',
    )
    log = lib_log.log
    log.info(f"Reading config: {args.config}")
    cfg = lib_cfg.setup_cfg(args.config)
    w = lib_walker.Walker()

    log.info("Pipeline parameters:")
    for key in lib_cfg._REQUIRED_KEYS:
        log.info(f"  {key} = {cfg[key]!r}")
    extra = {k: v for k, v in cfg.items() if k not in lib_cfg._REQUIRED_KEYS}
    for key, val in extra.items():
        log.info(f"  {key} = {val!r}  [extra]")

    if args.dry_run:
        log.info("Dry run — exiting before pipeline execution.")
        return

    # ------------------------------------------------------------------
    # Pipeline execution starts here
    # ------------------------------------------------------------------

    # todo: add cleanup
    with w.if_todo("cleanup"):
        log.info("Cleaning up intermediate files...")
        for pattern in [cfg["path_ms"], cfg["path_imgs"], cfg["path_plots"], cfg["path_logs"], cfg["path_sols"]]:
            log.info(f"Removing {pattern}")
            shutil.rmtree(pattern, ignore_errors=True)

    # create dirs if missing
    for d in [cfg["path_ms"], cfg["path_imgs"], cfg["path_plots"], cfg["path_logs"], cfg["path_sols"]]:
        os.makedirs(d, exist_ok=True)
    # add mepi dir to cfg
    cfg['mepi_dir'] = os.path.dirname(os.path.abspath(__file__))

    for step in cfg.get("pipeline_steps", []):
        log.info(f"Running pipeline step: {step}")
        module_name = f"mepi.pipeline_{step}"
        try:
            module = importlib.import_module(module_name)
            module.run()
        except ImportError as e:
            log.error(f"Failed to import module for step '{step}': {e}")
            sys.exit(1)
        except Exception as e:
            log.error(f"Error during pipeline step '{step}': {e}", exc_info=True)
            sys.exit(1)

if __name__ == "__main__":
    main()
