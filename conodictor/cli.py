"""Command line interface."""

import argparse
import sys
from pathlib import Path

VERSION = "2.3.6"

parser = argparse.ArgumentParser(
    prog="conodictor",
    usage="conodictor [options] <FILE>",
    add_help=False,
    epilog="Citation: Koua et al., 2021, Bioinformatics Advances",
)

parser.add_argument(
    "file",
    nargs="?",
    type=argparse.FileType("r"),
    default=sys.stdin,
    help=argparse.SUPPRESS,
)
parser.add_argument(
    "-d",
    "--dir",
    type=str,
    metavar="STR",
    help="specify database path",
)
parser.add_argument(
    "-o",
    "--out",
    nargs="?",
    type=Path,
    metavar="DIR",
    default="ConoDictor",
    help="output result to DIR [ConoDictor]",
)
parser.add_argument(
    "--mlen",
    type=int,
    metavar="INT",
    help="minimum length of sequences to be considered [off]",
)
parser.add_argument(
    "--ndup",
    type=int,
    metavar="INT",
    help="minimum occurence sequences to be considered [off]",
)
parser.add_argument(
    "--faa",
    action="store_true",
    help="dump a fasta file of matched sequences [false]",
)
parser.add_argument(
    "--filter",
    action="store_true",
    help="only keep sequences matching sig, pro and mat regions [false]",
)
parser.add_argument(
    "-a",
    "--all",
    action="store_true",
    help="add unclassified sequences in result [false]",
)
parser.add_argument(
    "-j",
    "--cpus",
    type=int,
    metavar="INT",
    default=1,
    help="number of threads [1]",
)
parser.add_argument(
    "--force",
    action="store_true",
    help="re-use output directory [false]",
)
parser.add_argument(
    "-q",
    "--quiet",
    action="store_true",
    help="decrease program verbosity",
)
parser.add_argument(
    "-v",
    "--version",
    action="version",
    version="%(prog)s " + f"{VERSION}",
)
parser.add_argument(
    "-h",
    "--help",
    action="help",
    help="show this help message and exit",
)
parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

args = parser.parse_args()
