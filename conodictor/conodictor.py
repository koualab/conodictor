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

import contextlib
import datetime
import logging
import os
import platform
import sys
from pathlib import Path

import pyfastx

from . import cli, conolib

AUTHOR = "Anicet Ebou and Dominique Koua"
URL = "https://github.com/koualab/conodictor.git"
VERSION = "2.4"

# Define command-line arguments----------------------------------------------
args = cli.args

QUIETNESS_LEVEL = logging.CRITICAL if args.quiet else logging.INFO

logging.basicConfig(
    format="[%(asctime)s][%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    level=QUIETNESS_LEVEL,
)


def main() -> None:
    """Predict and classify conopeptides in conotoxins superfamilies."""
    # Record start time
    startime = datetime.datetime.now(tz=datetime.timezone.utc)
    seqdata = None
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
    if not Path.exists(args.out):
        Path.mkdir(args.out)
    uncompressed_input = conolib.decompress_file(str(args.file.name), args.out)

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
    first_fasta_sequence = conolib.read_first_fasta_record(Path(uncompressed_input))
    fatype = conolib.isdnaorproteins(first_fasta_sequence)
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
    infile = pyfastx.Fasta(str(seqdata))
    seqids = infile.keys()
    # If --ndup is specified, get sequence ids of duplicate sequence
    dupdata = conolib.is_min_occurence_activated(
        args.ndup,
        args.mlen,
        seqdata,
        args.out,
    )

    # HMM search pipeline
    logging.info("Classifying sequences into conotoxins superfamilies")
    logging.info(
        "Searching hidden markov models profiles against the sequences",
    )
    conolib.run_hmm(seqdata, hres, cpus)

    logging.info("Parsing result")
    new_dict = conolib.transform_hmm_result(hres)
    logging.info("Compute combined evalue for each superfamily")
    results = conolib.compute_combined_evalue(new_dict)
    logging.info("Predicting sequence superfamily base on HMM")
    hmmfam = conolib.get_hmm_superfamily(results)

    # PSSM search pipeline
    logging.info("Running PSSM prediction")
    pssm_out = conolib.run_pssm(seqdata, args.out, cpus)
    logging.info("Parsing PSSM result")
    pssmdict = conolib.parse_pssm_result(pssm_out)
    pssmdict = conolib.clear_dict(pssmdict)
    logging.info("Predicting sequence superfamily base on PSSM")
    pssmfam = conolib.get_pssm_fam(pssmdict)

    # Writing output---------------------------------------------------------
    logging.info("Writing output")

    # Final families dict to store both predicted families
    finalfam = conolib.get_final_superfamilies(hmmfam, pssmfam, seqids)

    # Enter "reads" mode
    if args.mlen:
        conolib.write_result_read_mode(
            Path(args.out, "summary.csv"),
            seqdata,
            {"dupdata": dupdata, "finalfam": finalfam, "program": args.filter},
            report_all_seqs=args.all,
        )
    # "Transcriptome mode"
    else:
        conolib.write_result_transcriptome_mode(
            Path(args.out, "summary.csv"),
            finalfam,
            args.filter,
            report_all_seqs=args.all,
        )

    # Finishing -------------------------------------------------------------
    if not args.debug:
        logging.info("Cleaning around")
        with contextlib.suppress(OSError):
            Path.unlink(seqdata)
            Path.unlink(Path(f"{seqdata}.fxi"))
            Path.unlink(Path(args.out, "out.pssm"))

    logging.info("Creating donut plot")
    if args.mlen:
        conolib.donut_graph(
            6,
            Path(args.out, "summary.csv"),
            Path(args.out, "superfamilies_distribution.png"),
        )
    else:
        conolib.donut_graph(
            3,
            Path(args.out, "summary.csv"),
            Path(args.out, "superfamilies_distribution.png"),
        )
    logging.info("Done creating donut plot")
    logging.info("Classification finished successfully")
    logging.info("Check %s folder for results", args.out)
    logging.info("Walltime used (hh:mm:ss.ms): %s", conolib.elapsed_since(startime))
    if len(seqids) % 2:
        logging.info("Nice to have you. Share, enjoy and come back!")
    else:
        logging.info("Thanks you, come again.")


if __name__ == "__main__":
    main()
