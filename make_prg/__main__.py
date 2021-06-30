import argparse
import sys

from loguru import logger

from make_prg import __version__
from make_prg.subcommands import from_msa, update


def main():
    parser = argparse.ArgumentParser(
        prog="make_prg",
        usage="make_prg <subcommand> <options>",
        description="Subcommand entrypoint",
    )

    parser.add_argument("-V", "--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )

    msa_parser = from_msa.register_parser(subparsers)
    update_parser = update.register_parser(subparsers)

    for par in [msa_parser, update_parser]:
        par.add_argument(
            "-v",
            "--verbose",
            action="store_true",
            default=False,
            help="Increase output verbosity",
        )

    args = parser.parse_args()

    log_lvl = "DEBUG" if args.verbose else "INFO"
    handlers = [
        dict(sink=sys.stderr, enqueue=True, level=log_lvl),
    ]
    logger.configure(handlers=handlers)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
