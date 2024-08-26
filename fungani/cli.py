import argparse

ARG_MAX_WIDTH = 60

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
    "-p", "--percent", dest="percent", type=int, default=100, help="Genome sample size"
)
parser.add_argument(
    "-w", "--size", dest="size", type=int, default=1000, help="Window size"
)
parser.add_argument(
    "-g", "--overlap", dest="overlap", type=int, default=500, help="Window overlap"
)
parser.add_argument(
    "-j", "--cpus", dest="cpus", type=int, default=4, help="Number of CPU cores"
)
parser.add_argument("-o", "--output", dest="outdir", type=str, help="Output directory")

args = parser.parse_args()


if __name__ == "__main__":
    pass
