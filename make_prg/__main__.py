import argparse

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

    from_msa.register_parser(subparsers)
    update.register_parser(subparsers)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
