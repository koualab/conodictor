"""Sub-program for conodictor."""

from __future__ import annotations

import contextlib
import csv
import datetime
import importlib.resources
import logging
import operator
import os
import shutil
import subprocess
import sys
import warnings
from collections import Counter, defaultdict
from functools import reduce
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
import psutil
import pyfastx
import pyhmmer
import xphyle
from Bio.Seq import reverse_complement, translate
from matplotlib import pyplot as plt
from pyhmmer.hmmer import hmmsearch
from xphyle import paths

if TYPE_CHECKING:
    from pyhmmer.easel import Sequence, SequenceBlock, SequenceFile
    from pyhmmer.plan7 import HMMFile


VERSION = "2.4.1"
PSSM_SEQ_ID = 3
CONOPEP_FAMILY = 0
CONOPEP_FAMILY_NAME = 1
PRO_REGION = 2


def write_fasta(fampos: dict, hmmfam: dict, seqdata: Path, out: str) -> None:
    """Write a fasta file of matched sequences."""
    infile = pyfastx.Fasta(str(seqdata))
    output_path = Path(out) / f"{seqdata.stem}_predicted.fa"

    logging.info("Writing out fasta file of matched sequences")

    with output_path.open("w") as faah:
        for seqid, pos in fampos.items():
            start, end = map(int, pos.split("-"))
            sequence = infile[seqid].seq[start:end]
            faah.write(f">{seqid} conodictor={hmmfam[seqid]}\n{sequence}\n")


def is_min_occurence_activated(
    min_occ: int,
    min_len: int,
    seqdata: Path,
    out: str,
) -> dict:
    """Return duplicated sequence stats."""
    infile = pyfastx.Fasta(str(seqdata))
    seqids = infile.keys()
    # If --ndup is specified, get sequence ids of duplicate sequence
    dupdata = {}
    if min_occ:
        dupdata = get_dup_seqs(seqdata, seqids, min_occ)
        ldu = len(dupdata)
        if min_len is None:
            logging.info(
                "Input file contains %s sequences with at least %s occurences."
                " Only these sequences will be used for prediction.",
                ldu,
                min_occ,
            )
        elif ldu == 0:
            logging.info("We have 0 sequences with at least %s occurences.", min_occ)
            logging.info("No prediction will therefore be made. Stopping...")
            sys.exit(1)
        else:
            logging.info(
                "And from them we have %s sequences with at least %s occurences",
                len(dupdata),
                min_occ,
            )
            logging.info("Only these sequences will be used for prediction")

    # Create a fasta file of sequence after filtering
    if min_occ is not None:
        with Path.open(Path(out, "filtfa.fa"), "w") as fih:
            for kid in dupdata:
                fih.write(f">{infile[kid].description}\n{infile[kid].seq}\n")
        fih.close()

        # Use the filtered file as input of further commands
        seqdata = Path(out, "filtfa.fa")
    elif min_occ is None and min_len is not None:
        with Path.open(Path(out, "filtfa.fa"), "w") as fih:
            for kid in seqids:
                fih.write(f">{infile[kid].description}\n{infile[kid].seq}\n")
        fih.close()

        # Use the filtered file as input of further commands
        seqdata = Path(out, "filtfa.fa")

    return dupdata


def get_fam_or_unknown(seqid: str, hmmfam: dict, pssmfam: dict) -> list[str]:
    """Helper function to get family or unknown."""  # noqa: D401
    return [
        hmmfam.get(seqid, "UNKNOWN"),
        pssmfam.get(seqid, "UNKNOWN"),
        definitive_prediction(
            hmmfam.get(seqid, "UNKNOWN"),
            pssmfam.get(seqid, "UNKNOWN"),
        ),
    ]


def get_final_superfamilies(hmmfam: dict, pssmfam: dict, seqids: list) -> defaultdict:
    """Get final superfamilies."""
    finalfam = defaultdict(list)

    # Collect all known sequence IDs
    known_seqs = set(hmmfam) | set(pssmfam)

    # Assign families to known sequences
    for seqid in known_seqs:
        finalfam[seqid] = get_fam_or_unknown(seqid, hmmfam, pssmfam)

    # Assign "UNKNOWN" to unknown sequences
    for sid in seqids:
        if sid not in finalfam:
            finalfam[sid] = ["UNKNOWN", "UNKNOWN", "UNKNOWN"]

    return finalfam


def run_hmm(seqdata: Path, hres: dict, cpus: int) -> None:
    """Run hmm search pipeline."""
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

        search_hmm_profiles(hres, hmm_file, cpus, targets)


def parse_pssm_result(pssm_result: Path) -> defaultdict:
    """Parse pssm result."""
    pssmdict = defaultdict(lambda: defaultdict(list))

    with Path.open(pssm_result) as pssmfile:
        rd = csv.reader(pssmfile, delimiter="\t")
        for row in rd:
            pssmdict[row[PSSM_SEQ_ID]][
                (row[CONOPEP_FAMILY].split("|")[CONOPEP_FAMILY]).split("_")[
                    CONOPEP_FAMILY_NAME
                ]
            ].append(
                (row[CONOPEP_FAMILY].split("|")[CONOPEP_FAMILY]).split("_")[PRO_REGION],
            )
    return pssmdict


def write_result_read_mode(
    out: str,
    seqdata: Path,
    data: dict,
    *,
    report_all_seqs: bool,
) -> None:
    """Write result to output in read mode.

    Args:
    ----
    out: Output directory path.
    seqdata: Path to sequence data.
    data: Dictionary containing `dupdata`, `finalfam`, and `program`.
    report_all_seqs: Boolean flag to report all sequences.

    """
    outfile_path = Path(out, "summary.csv")

    uniq_final = select_filter(data["program"], data["finalfam"])

    with Path.open(outfile_path, "w") as outfile:
        outfile.write(
            "sequence,length,num_cysteines,occurence,"
            "hmm_pred,pssm_pred,definitive_pred\n",
        )

        def write_sequence(sequence_id: str, sequence_data: list) -> None:
            length, num_cysteines = get_stats(sequence_id, str(seqdata))
            occurence = data["dupdata"].get(sequence_id, 1)
            outfile.write(
                f"{sequence_id},{length},{num_cysteines},{occurence},"
                f"{sequence_data[0]},{sequence_data[1]},{sequence_data[2]}\n",
            )

        if report_all_seqs:
            for seq_id, seq_data in data["finalfam"].items():
                write_sequence(seq_id, seq_data)
        else:
            for seq_id, seq_data in uniq_final.items():
                write_sequence(seq_id, seq_data)


def write_result_transcriptome_mode(
    out: Path,
    finalfam: defaultdict,
    program: str,
    *,
    report_all_seqs: bool,
) -> None:
    """Write result to output in transcriptome mode."""
    uniq_final = select_filter(program, finalfam)

    with Path.open(out, "w") as outfile:
        outfile.write("sequence,hmm_pred,pssm_pred,definitive_pred\n")
        data_to_write = finalfam.items() if report_all_seqs else uniq_final.items()

        for k, v in data_to_write:
            outfile.write(f"{k},{v[0]},{v[1]},{v[2]}\n")


def select_filter(profile: str, mydict: dict) -> dict:
    """Filter result by program.

    Args:
    ----
    profile: type of profile. Either HMM or PSSM
    mydict: dictionary of classifications

    Returns:
    -------
    A dict with the unknown superfamilies removed

    """
    if profile not in {"pssm", "hmm"}:
        return {
            k: v for k, v in mydict.items() if v != ["UNKNOWN", "UNKNOWN", "UNKNOWN"]
        }

    index = 0 if profile == "pssm" else 1
    return {
        k: v
        for k, v in mydict.items()
        if v != ["UNKNOWN", "UNKNOWN", "UNKNOWN"] and v[index] != "UNKNOWN"
    }


def donut_graph(ncol: int, stat_file: Path, outfile: Path) -> None:
    """Create a donut graph from the outputted stats of predicted sequences.

    Args:
    ----
    ncol: The column number in the CSV file to use for the graph.
    stat_file: Path to the CSV file containing the statistics.
    outfile: Path to save the output graph.

    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Read the CSV file and extract the relevant column data
        data = pd.read_csv(stat_file)
        plot_data = data.iloc[:, ncol].tolist()

        # Count occurrences of each unique value in the column
        dtc = Counter(plot_data)

        # Filter out "CONFLICT" and "UNKNOWN" entries from labels and values
        filtered_items = [
            (k, v)
            for k, v in sorted(dtc.items())
            if not k.startswith(("CONFLICT", "UNKNOWN"))
        ]
        labels = [f"{k}: {v}" for k, v in filtered_items]
        values = [v for _, v in filtered_items]

        # Create the donut plot
        fig, ax = plt.subplots(figsize=(13, 10), subplot_kw={"aspect": "equal"})
        wedges, _ = ax.pie(  # type: ignore[type-supported]
            values,
            wedgeprops={"width": 0.5},
            startangle=-40,
            shadow=False,
        )

        # Add legend and title
        ax.legend(wedges, labels, loc="lower center", ncol=6)
        ax.set_title("ConoDictor Predictions")

        # Add version text
        plt.text(-2, -1.5, f"Made with ConoDictor v{VERSION}")

        # Save the figure
        plt.savefig(outfile, dpi=300)
        plt.close(fig)  # Close the figure to avoid resource warning


def get_stats(seqid: str, file_path: str) -> list[int]:
    """get_stats return sequence length and number of cysteines in a sequences for a sequence id.

    Args:
    ----
    seqid: Input sequence id list.
    file_path: Fasta file to use to retrieve sequence.

    """  # noqa: E501
    stats = []
    infile = pyfastx.Fasta(file_path)

    # Sequence length
    stats.append(len(infile[seqid]))
    # Number of cysteines in sequence
    stats.append(infile[seqid].seq.count("C"))

    return stats


def get_dup_seqs(file_path: Path, idslist: list[str], mnoc: int) -> dict:
    """Search provided fasta file for duplicate sequences and return sequence ids of duplicate sequences.

    Args:
    ----
    file_path: Input fasta file to use for search.
    idslist: Sequence ids list to consider.
    mnoc: Minimum number of occurrences wanted.

    Returns:
    -------
    A dictionary of sequence ids with their duplicate count.

    """  # noqa: E501
    dupid = {}
    infile = pyfastx.Fasta(str(file_path))
    flipped = defaultdict(list)

    # Create flipped dictionary directly from infile
    for sid in idslist:
        seq = infile[sid].seq
        flipped[seq].append(sid)

    # Only consider sequences with enough duplicates
    for v in flipped.values():
        if len(v) >= mnoc:
            dupid[v[0]] = len(v)

    return dupid


def definitive_prediction(hmmclass: str, pssmclass: str) -> str:
    """Combine HMM and PSSM predictions to give a definitive classification.

    Args:
    ----
    hmmclass: HMM predicted superfamily, required (string)
    pssmclass: PSSM predicted superfamily, required (string)

    Returns:
    -------
    The definitve classification string.

    """
    if hmmclass == pssmclass:
        return hmmclass

    if "CONFLICT" in pssmclass and "CONFLICT" in hmmclass:
        fams_pssm = pssmclass.replace("CONFLICT ", "").split(" and ")
        fams_hmm = hmmclass.replace("CONFLICT ", "").split(" and ")
        return f"CONFLICT {', '.join(fams_pssm + fams_hmm)}"

    if "CONFLICT" in pssmclass and "CONFLICT" not in hmmclass:
        return hmmclass

    if ("CONFLICT" in hmmclass and "CONFLICT" not in pssmclass) or (
        "UNKNOWN" in hmmclass and "UNKNOWN" not in pssmclass
    ):
        return pssmclass

    if "UNKNOWN" not in hmmclass and "UNKNOWN" in pssmclass:
        return hmmclass

    return f"CONFLICT {hmmclass} and {pssmclass}"


def get_pssm_fam(mdict: dict) -> dict:
    """Return the superfamily based on PSSM with the highest occurrence.

    Args:
    ----
    mdict: Dictionary containing sequences and their associated families with matches.

    Returns:
    -------
    Dictionary mapping sequence IDs to the most common family based on PSSM matches.

    """
    pssmfam = {}
    for key, subdict in mdict.items():
        # Flatten the list of families and count occurences
        families = [fam for sublist in subdict.values() for fam in sublist]
        count = Counter(families)

        # Find the most common families
        most_common_fams = count.most_common(2)

        if len(most_common_fams) == 1:
            fam = most_common_fams[0][0]
        else:
            top_fam, second_fam = most_common_fams
            if top_fam[1] == second_fam[1]:
                fam = f"CONFLICT {top_fam[0]} and {second_fam[0]}"
            else:
                fam = top_fam[0]
        pssmfam[key] = fam
    return pssmfam


def clear_dict(
    hdict: defaultdict[str, defaultdict[str, list[str]]],
) -> dict:
    """Filter out sequences without a MATURE PSSM profile."""
    # Using a list to collect empty top-level keys for removal
    empty_keys = []
    for k in list(hdict.keys()):
        subdict = hdict[k]
        for a in list(subdict.keys()):
            if "MAT" not in subdict[a]:
                del subdict[a]
        # Check if the subdictionary is empty after deletion
        if not subdict:
            empty_keys.append(k)
    # Remove top-level keys with empty subdictionaries
    for k in empty_keys:
        del hdict[k]
    return hdict


def run_pssm(infile: Path, outdir: str, cpus: int) -> Path:
    """Run PSSM pipeline."""
    pssm_executable = shutil.which("pfscanV3")
    pssm = importlib.resources.files("conodictor").joinpath("conodictor.pssm")
    if is_num_records_greater_than(infile, 100000):
        split_dir = split_fasta_file(infile)
        subfiles = os.listdir(Path(split_dir))
        with Path.open(Path(outdir, "out.pssm"), "a") as pssm_out:
            for file in subfiles:
                pssm_run = subprocess.run(  # noqa: S603
                    [
                        str(pssm_executable),
                        "--nthreads",
                        str(cpus),
                        "-o",
                        "7",
                        str(pssm),
                        "-f",
                        file,
                    ],
                    capture_output=True,
                    check=False,
                )
                pssm_out.write(pssm_run.stdout.decode("utf-8"))
    else:
        pssm_run = subprocess.run(  # noqa: S603
            [
                str(pssm_executable),
                "--nthreads",
                str(cpus),
                "-o",
                "7",
                str(pssm),
                "-f",
                infile,
            ],
            capture_output=True,
            check=False,
        )
        with Path.open(Path(outdir, "out.pssm"), "w") as pssm_out:
            pssm_out.write(pssm_run.stdout.decode("utf-8"))

    return Path(outdir, "out.pssm")


def split_fasta_file(input_file: Path, max_records: int = 25000) -> Path:
    """Splits a FASTA file into multiple files with a maximum number of records in each.

    Args:
    ----
    input_file: Path to the input FASTA file.
    max_records: Maximum number of records per split file.

    Returns:
    -------
    Path to the directory containing the split files.

    """  # noqa: D401
    # Create a directory to store the split files
    base_name = Path(input_file).name
    split_dir = Path(input_file).parent / f"split_{base_name}"
    Path(split_dir).mkdir(parents=True)

    file_count = 0
    record_count = 0
    output_file = None

    with Path.open(input_file) as infile:
        for line in infile:
            if line.startswith(">"):
                if record_count >= max_records:
                    if output_file:
                        output_file.close()
                    record_count = 0
                    file_count += 1

                if record_count == 0:
                    output_file_path = (
                        Path(split_dir) / f"{base_name}_part{file_count + 1}.fa"
                    )

                    output_file = Path.open(output_file_path, "w")

                record_count += 1

            if output_file:
                output_file.write(line)

    if output_file:
        output_file.close()

    return split_dir


def is_num_records_greater_than(file_path: Path, limit: int) -> bool:
    """Check number of records of fasta file."""
    i = 0
    with Path.open(file_path, "r") as file:
        for line in file:
            if line.startswith(">"):
                i += 1
    return i >= limit


def read_first_fasta_record(file_path: Path) -> str:
    """Return the first record of a FASTA file.

    :param file_path: Path to the FASTA file.
    :return: A tuple containing the header and sequence of the first record.
    """
    header = None
    sequence = []

    with Path.open(file_path) as file:
        for line in file:
            if line.startswith(">"):
                if header is not None:
                    # If we already have a header, we've read the first record
                    break
                header = line.strip()
            else:
                sequence.append(line.strip())

    # Combine the sequence list into a single string
    return "".join(sequence)


def isdnaorproteins(s: str) -> str:
    """Determine if the input sequence is DNA or proteins.

    Args:
    ----
    s: input sequence

    Returns:
    -------
    'DNA' if the sequence is DNA, 'protein' if it is a protein, otherwise 'unknown'

    """
    dna_set = {"A", "T", "C", "G"}
    protein_set = {
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
        "Z",
        "*",
        "X",
    }

    if set(s).issubset(dna_set):
        return "DNA"
    elif set(s).issubset(protein_set):
        return "protein"
    else:
        return "unknown"


def decompress_file(input_file: str, outdir: str) -> str:
    """Decompress a file and return the path to the uncompressed file."""
    output_path = Path(outdir, "conodictor_input.fa")
    # Handling input file supply
    if input_file == "<stdin>":
        with output_path.open("w+", encoding="utf-8") as binput, xphyle.xopen(
            paths.STDIN,
            context_wrapper=True,
        ) as infile:
            for line in infile:
                binput.write(line)  # type: ignore[type-supported]
    else:
        with output_path.open("w+", encoding="utf-8") as binput, xphyle.xopen(
            input_file,
            context_wrapper=True,
        ) as infile:
            for line in infile:
                binput.write(line)  # type: ignore[type-supported]
    return str(output_path)


def check_input(infile: Path) -> str:
    """Check input fasta file."""
    error = ""
    fasta = pyfastx.Fasta(str(infile))
    # Check if fasta file does not contain duplicate sequences
    # which would break hmmsearch
    ids = fasta.keys()
    if len(ids) != len(set(ids)):
        error = "Supplied FASTA file contains duplicate sequences."
        with contextlib.suppress(OSError):
            Path.unlink(Path(f"{Path}.fxi"))
    return error


def handle_cpus(asked_cpus: int, available_cpus: int | None) -> int:
    """Allocate the good number of CPUs based on asked cpus vs available cpus."""
    cpus = 1
    if asked_cpus == 0 or asked_cpus > int(available_cpus or 1):
        cpus = available_cpus
    else:
        cpus = asked_cpus
    return int(cpus or 1)


def search_hmm_profiles(
    result: dict,
    hmmfile: HMMFile,
    cpus: int,
    seqs: SequenceBlock[Sequence] | SequenceFile,
) -> None:
    """Perform hmmsearch of `hmmfile` agains `seqs` using `cpus` CPUs."""
    for hits in hmmsearch(hmmfile, seqs, cpus=cpus):  # type: ignore[type-supported]
        hmm_id = hits.query_name
        for hit in hits:
            if hit.included:
                hit_name = str(hit.name, encoding="utf-8")
                parts = str(hmm_id, encoding="utf-8").split("_")  # type: ignore[str-conversion-from-bytes]
                for domain in hit.domains:
                    hmmresult = {
                        "fam": parts[1],
                        "part": parts[2],
                        "evalue": domain.i_evalue,
                        "score": hit.score,
                        "start": domain.env_from,
                        "end": domain.env_to,
                    }

                if hit_name in result:
                    result[hit_name].append(hmmresult)
                else:
                    result[hit_name] = [hmmresult]


def get_hmm_superfamily(
    result_with_combined_evalue: dict,
) -> tuple[dict, dict]:
    """Compare combined_evalues to get superfamilies."""
    seqfam = {}
    fampose = {}
    for key, (
        fam_list,
        max_combined_evalue,
        max_fam,
    ) in result_with_combined_evalue.items():
        res = ""
        for fam_dict in fam_list:
            if (
                max_combined_evalue >= 100 * fam_dict["combined_evalue"]
                and fam_dict["combined_evalue"] != max_combined_evalue
            ):
                res = fam_dict["fam"]
        conflicts = [
            fam_dict["fam"]
            for fam_dict in fam_list
            if fam_dict["combined_evalue"] == max_combined_evalue
        ]
        if len(conflicts) > 1:
            res = f"Conflict between {', '.join(conflicts)}"
        else:
            min_start = min(info["start"] for info in max_fam["info"])
            max_end = max(info["end"] for info in max_fam["info"])
            res = conflicts[0]
        seqfam[key] = res
        fampose[key] = f"{min_start}-{max_end}"
    return seqfam, fampose


def compute_combined_evalue(transformed_result: dict) -> dict:
    """Compute the combined evalue for each fam and find the maximum combined evalue."""
    results = {}
    for key, fam_list in transformed_result.items():
        max_combined_evalue = None
        max_fam = None
        for fam_dict in fam_list:
            combined_evalue = reduce(
                operator.mul,
                (info["evalue"] for info in fam_dict["info"]),
                1,
            )
            fam_dict["combined_evalue"] = combined_evalue
            # Save start and end of the fam with the max combined evalue
            if max_combined_evalue is None or combined_evalue > max_combined_evalue:
                max_combined_evalue = combined_evalue
                max_fam = fam_dict
        results[key] = fam_list, max_combined_evalue, max_fam
    return results


def transform_hmm_result(hmm_result: dict) -> dict:
    """Transform HMM result into a new dictionnary."""
    new_dict = {}
    for key, value in hmm_result.items():
        fam_dict = defaultdict(list)
        for item in value:
            fam = item.pop("fam")
            fam_dict[fam].append(item)
        # Filter out 'fam' that do not contain 'part' equal to 'MAT'
        filtered_fam_dict = {
            fam: info
            for fam, info in fam_dict.items()
            if any(i["part"] == "MAT" for i in info)
        }
        if filtered_fam_dict:
            new_dict[key] = [
                {"fam": fam, "info": info} for fam, info in filtered_fam_dict.items()
            ]

    return new_dict


def elapsed_since(start: datetime.datetime) -> datetime.timedelta:
    """Compute time between two datetime."""
    return datetime.datetime.now(tz=datetime.timezone.utc) - start


def do_translation(infile: str, outfile: Path, sw: int = 60) -> None:
    """Translate a DNA fasta file into a proteins fasta file.

    Args:
    ----
    infile: input file path as string.
    outfile: Output file as path.
    sw: Sequence width. Default: 60.

    Returns:
    -------
    Translated sequences

    """
    fa = pyfastx.Fasta(infile)
    output_path = Path(outfile)
    # Use buffer to reduce I/O operations
    buffer_size = 1024 * 1024  # 1 MB buffer
    buffer = []
    with output_path.open("w") as protfile:
        for sequence in fa:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                protseq = _translate_seq(sequence.seq)
                for idx, frame in enumerate(protseq):
                    seq_letters = [frame[i : i + sw] for i in range(0, len(frame), sw)]
                    nl = "\n"
                    buffer.append(
                        f">{sequence.name}_frame={idx + 1}\n{nl.join(seq_letters)}\n",
                    )
                    if len(buffer) > buffer_size:
                        protfile.write("".join(buffer))
                        buffer = []
        # Write remaining buffer to file
        if buffer:
            protfile.write("".join(buffer))


def _translate_seq(seq: str) -> list[str]:
    """Translate DNA sequence to proteins in the six frames.

    Args:
    ----
    seq: DNA sequence to translate.

    Returns:
    -------
    A list of amino-acids in six frames

    """
    rc_seq = reverse_complement(seq)
    return [
        str(translate(seq)),
        str(translate(seq[1:])),
        str(translate(seq[2:])),
        str(translate(rc_seq)),
        str(translate(rc_seq[1:])),
        str(translate(rc_seq[2:])),
    ]
