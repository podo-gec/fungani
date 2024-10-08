import argparse
import sys

from fungani.core import main

ARG_MAX_WIDTH = 60


def parse_args(args):
    parser = argparse.ArgumentParser(
        prog="fungani",
        description="Compute ANI using Blast",
        formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(
            prog, max_help_position=ARG_MAX_WIDTH
        ),
    )
    parser.add_argument("reference", type=str, help="Reference genome")
    parser.add_argument("test", type=str, help="Test genome")
    parser.add_argument(
        "-t",
        "--threshold",
        dest="threshold",
        type=int,
        default=80,
        help="ANI threshold",
    )
    parser.add_argument(
        "-p",
        "--percent",
        dest="percent",
        type=int,
        default=10,
        help="Genome fraction",
    )
    parser.add_argument(
        "-w", "--size", dest="size", type=int, default=1000, help="Window size"
    )
    parser.add_argument(
        "-g", "--overlap", dest="overlap", type=int, default=500, help="Step size"
    )
    parser.add_argument(
        "-j", "--cpus", dest="cpus", type=int, default=4, help="Number of CPU cores"
    )
    parser.add_argument(
        "-o", "--output", dest="outdir", type=str, help="Output directory"
    )
    parser.add_argument(
        "-u",
        "--onepass",
        action="store_true",
        help="Only in one direction only",
    )
    parser.add_argument(
        "-c", "--clean", action="store_true", help="Clean intermediate files"
    )

    return parser.parse_args(args)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    if args.onepass:
        args.mode = "fwd"
        main(args)
    else:
        args.mode = "fwd"
        main(args)
        args.mode = "rev"
        args.test, args.reference = args.reference, args.test
        main(args)
