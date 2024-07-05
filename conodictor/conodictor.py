#!/usr/bin/env python3
"""Main code of conodictor."""

# ConoDictor: Prediction and classification of conopeptides
# Copyright (C) 2019-2024  Koualab
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
from __future__ import annotations

import datetime
import importlib.resources
import logging
import os
import platform
import sys
from pathlib import Path

import psutil
import pyfastx
import pyhmmer

from . import cli, conolib

AUTHOR = "Anicet Ebou and Dominique Koua"
URL = "https://github.com/koualab/conodictor.git"
VERSION = "2.4"

# Some global variables
UNKNOWN_FAM = "UNKNOWN"
CONFLICT_FAM = "CONFLICT"
PSSM_SEQ_ID = 3

# Define command-line arguments----------------------------------------------
args = cli.args

QUIETNESS_LEVEL = logging.CRITICAL if args.quiet else logging.INFO

logging.basicConfig(
    format="[%(asctime)s][%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    level=QUIETNESS_LEVEL,
)


def main() -> None:
    """Contains main program of conodictor."""
    # Record start time
    startime = datetime.datetime.now(tz=datetime.timezone.utc)
    # Start program ---------------------------------------------------------
    logging.info("This is conodictor %s", VERSION)
    logging.info("Written by %s", AUTHOR)
    logging.info("Available at %s", URL)
    logging.info(
        "Localtime is %s",
        datetime.datetime.now(tz=datetime.timezone.utc).strftime("%H:%M:%S"),
    )
    # Get current user name
    try:
        user = os.environ["USER"]
    except KeyError:
        user = "not telling me who you are"
    logging.info("You are %s", user)
    logging.info("Operating system is %s", platform.system())

    # Decompressing input file if needed
    uncompressed_input = conolib.decompress_file(str(args.file.name))

    # Check input for various potential problems
    check = conolib.check_input(Path(uncompressed_input))
    if check != "":
        logging.error(check)
        sys.exit(1)

    # Handling number of threads -----------------------------------------------
    available_cpus = os.cpu_count()
    logging.info("System has %s cores", available_cpus)
    cpus = conolib.handle_cpus(args.cpus, available_cpus)
    logging.info("We will use maximum of %s cores", cpus)

    # Translating sequences if needed
    fa = pyfastx.Fasta(uncompressed_input)
    fatype = conolib.isdnaorproteins(fa[0].seq)
    if fatype == "DNA":
        logging.info("Translating sequences into 6 frames")
        # translate sequences
        conolib.do_translation(
            uncompressed_input,
            Path(args.out, "conodictor_translate.fa"),
        )
        seqdata = Path(args.out, "conodictor_translate.fa")
    else:
        seqdata = Path(uncompressed_input)

    hres = {}

    # HMM search pipeline
    logging.info("Classifying sequences into conotoxins superfamilies")
    logging.info("Step 1. using hidden Markov models")
    available_memory = psutil.virtual_memory().available
    target_size = Path.stat(seqdata).st_size
    hmm = importlib.resources.files("conodictor").joinpath("conodictor.hmm").open("rb")

    with pyhmmer.plan7.HMMFile(hmm) as hmm_file, pyhmmer.easel.SequenceFile(  # type: ignore[type-supported]
        seqdata,
        digital=True,
    ) as seq_file:
        if target_size < available_memory * 0.1:
            logging.info("Pre-fetching hidden Markov models database into memory")
            targets = seq_file.read_block()
            mem = sys.getsizeof(targets) + sum(
                sys.getsizeof(target) for target in targets
            )
            mem = mem / 1024
            logging.info("Database in-memory size: %s KiB", f"{mem:.1f}")
        else:
            targets = seq_file
        logging.info(
            "Searching hidden markov models profiles against the sequences",
        )
        conolib.search_hmm_profiles(hres, hmm_file, cpus, targets)

    logging.info("Parsing result")
    new_dict = conolib.transform_hmm_result(hres)

    # Compute the combined evalue for each fam and find the maximum combined evalue
    logging.info("Compute combined evalue for each superfamily")
    results = conolib.compute_combined_evalue(new_dict)

    # Compare combined_evalues
    logging.info("Predicting sequence superfamily")
    seqfam = conolib.get_hmm_superfamily(results)
    print(seqfam)

    logging.info("Walltime used (hh:mm:ss.ms): %s", conolib.elapsed_since(startime))


if __name__ == "__main__":
    main()
