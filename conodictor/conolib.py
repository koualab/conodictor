"""Sub-program for conodictor."""

from __future__ import annotations

import contextlib
import csv
import datetime
import importlib.resources
import logging
import operator
import os
import re
import shutil
import subprocess
import sys
import warnings
from collections import Counter, defaultdict
from functools import reduce
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
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


VERSION = "2.4"
PSSM_SEQ_ID = 3
CONOPEP_FAMILY = 0
CONOPEP_FAMILY_NAME = 1
PRO_REGION = 2


def get_final_superfamilies(hmmfam: dict, pssmfam: dict, seqids: list) -> defaultdict:
    """Get final superfamilies."""
    # Final families dict to store both predicted families
    finalfam = defaultdict(list)

    # Itterate through all submitted sequence to assign families
    known_seqs = []
    known_seqs.extend([*hmmfam])
    known_seqs.extend([*pssmfam])

    for seqid in known_seqs:
        finalfam[seqid].extend(get_fam_or_unknown(seqid, hmmfam, pssmfam))

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
    out: Path,
    seqdata: Path,
    dupdata: dict,
    finalfam: defaultdict,
    program: str,
    *,
    report_all_seqs: bool,
) -> None:
    """Write result to output in read mode."""
    infile = pyfastx.Fasta(str(seqdata))
    outfile = Path.open(Path(out, "summary.csv"), "w")
    uniq_final = select_filter(program, finalfam)
    # write summary.txt file with sequence stats
    outfile.write(
        "sequence,length,num_cysteines,occurence,"
        "hmm_pred,pssm_pred,definitive_pred\n",
    )

    if not report_all_seqs:
        for uk, uv in uniq_final.items():
            outfile.write(
                (
                    f"{uk},{get_stats(uk, infile)[0]},"
                    f"{get_stats(uk, infile)[1]},"
                    f"{dupdata[uk]},{uv[0]},{uv[1]},{uv[2]}\n"
                ),
            )
        outfile.close()
    else:
        for k, v in finalfam.items():
            outfile.write(
                (
                    f"{k},{get_stats(k, infile)[0]},"
                    f"{get_stats(k, infile)[1]},"
                    f"{dupdata[k]},"
                    f"{v[0]},{v[1]},{v[2]}\n"
                ),
            )
        outfile.close()


def write_result_transcriptome_mode(
    out: Path,
    finalfam: defaultdict,
    program: str,
    *,
    report_all_seqs: bool,
) -> None:
    """Write result to output in transcriptome mode."""
    outfile = Path.open(Path(out, "summary.csv"), "w")
    uniq_final = select_filter(program, finalfam)

    # Open output file for writing
    outfile.write("sequence,hmm_pred,pssm_pred,definitive_pred\n")
    # Make reporting unclassified sequences optional
    if not report_all_seqs:
        # Write output
        for uk, uv in uniq_final.items():
            outfile.write(
                f"{uk},{uv[0]},{uv[1]},{uv[2]}\n",
            )
        outfile.close()
    else:
        for k, v in finalfam.items():
            outfile.write(f"{k},{v[0]},{v[1]},{v[2]}\n")
        outfile.close()


def select_filter(program: str, mydict: dict) -> dict:
    """Filter result by program."""
    uniq_final = {}
    # Get final sequences which has been classified
    if program == "pssm":
        uniq_final = {
            k: v
            for k, v in mydict.items()
            if v != ["UNKNOWN", "UNKNOWN", "UNKNOWN"]
            if v[0] != "UNKNOWN"
        }
    elif program == "hmm":
        uniq_final = {
            k: v
            for k, v in mydict.items()
            if v != ["UNKNOWN", "UNKNOWN", "UNKNOWN"]
            if v[1] != "UNKNOWN"
        }
    else:
        uniq_final = {
            k: v for k, v in mydict.items() if v != ["UNKNOWN", "UNKNOWN", "UNKNOWN"]
        }
    return uniq_final


def donut_graph(ncol: int, stat_file: Path, outfile: Path) -> None:
    """donut_graph make a donut graph from outputed stats of predicted sequences."""
    data = pd.read_csv(stat_file)
    plot_data = data[data.columns[ncol]].tolist()
    dtc = Counter(plot_data)
    labels = [
        f"{k1}: {v1}"
        for k1, v1 in sorted(dtc.items())
        if not k1.startswith(("CONFLICT", "UNKNOWN"))
    ]
    values = [
        x for k2, x in sorted(dtc.items()) if not k2.startswith(("CONFLICT", "UNKNOWN"))
    ]

    # White circle
    _, ax = plt.subplots(figsize=(13, 10), subplot_kw={"aspect": "equal"})
    wedges, _ = ax.pie(  # type: ignore[type-supported]
        np.array(values).ravel(),
        wedgeprops={"width": 0.5},
        startangle=-40,
        shadow=False,
    )
    # bbox: x, y, width, height
    ax.legend(wedges, labels, loc="lower center", ncol=6)
    ax.set_title("ConoDictor Predictions")
    plt.text(-2, -1.5, f"Made with ConoDictor v{VERSION}")
    plt.savefig(outfile, dpi=300)


def get_stats(seqid: str, file_path: Path) -> list[int]:
    """get_stats return sequence length and number of cysteines in a sequences for a sequence id.

    :id: Input sequence id list.
    :infile: Fasta file to use to retrieve sequence.
    """  # noqa: E501
    stats = []
    infile = pyfastx.Fasta(file_path)

    # Sequence length
    stats.append(len(infile[seqid]))
    # Number of cysteines in sequence
    stats.append(infile[seqid].seq.count("C"))

    return stats


def get_dup_seqs(file_path: Path, idslist: list[str], mnoc: int) -> dict:
    """get_dup_seqs search provided fasta file for duplicate sequence.

    Return sequence ids of duplicate sequences.

    :infile: Input fasta file to use for search
    :idslist: Sequence ids list to consider
    :mnoc: Minimum number of occurence wanted
    """
    dupid = {}
    infile = pyfastx.Fasta(file_path)
    flipped = defaultdict(set)
    seqdict = defaultdict()
    for sid in idslist:
        seqdict[sid] = infile[sid].seq

    flipped = _flip_dict(seqdict)

    # The flipped dict is a dict of list like: dict = {"ATCT": [id1, id2], "GCTA": [id4, id5]}  # noqa: E501
    # returning only the first element of the value
    # let us consider only one occurence of the sequence
    # Therefore we will really only predict one sequence
    # which can have multiple occurence
    for v in flipped.values():
        if len(v) >= mnoc:
            dupid[v[0]] = len(v)

    return dupid


def _flip_dict(mydict: dict) -> defaultdict:
    """Return a flipped dict of input dict."""
    flipped_ = defaultdict(list)
    for k, v in mydict.items():
        if v not in flipped_:
            flipped_[v] = [k]
        else:
            flipped_[v].append(k)

    return flipped_


def definitive_prediction(hmmclass: str, pssmclass: str) -> str:
    """definitive_prediction gives definitive classification by combining HMM and PSSM.

    :hmmclass: HMM predicted family, required (string)
    :pssmclass: PSSM predicted family, required (string)
    """
    deffam = None

    if hmmclass == pssmclass:
        deffam = hmmclass
    elif "CONFLICT" in pssmclass and "CONFLICT" in hmmclass:
        fams_pssm = re.search("(?<=CONFLICT)(.*)and(.*)", pssmclass)
        fams_hmm = re.search("(?<=CONFLICT)(.*)and(.*)", hmmclass)
        deffam = f"CONFLICT {fams_pssm.group(1)}, {fams_pssm.group(2)}, {fams_hmm.group(1)}, and {fams_hmm.group(2)}"  # type: ignore[method-supported]  # noqa: E501
    elif "CONFLICT" in pssmclass and "CONFLICT" not in hmmclass:
        deffam = hmmclass
    elif (
        "CONFLICT" in hmmclass
        and "CONFLICT" not in pssmclass
        or "UNKNOWN" in hmmclass
        and "UNKNOWN" not in pssmclass
    ):
        deffam = pssmclass
    elif "UNKNOWN" not in hmmclass and "UNKNOWN" in pssmclass:
        deffam = hmmclass
    elif pssmclass != hmmclass:
        deffam = f"CONFLICT {hmmclass} and {pssmclass}"

    return deffam or ""


def get_fam_or_unknown(seqid: str, hmmfam: dict, pssmfam: dict) -> list[str]:
    """Get fam from both dict."""
    fam = []
    if seqid in hmmfam:
        fam.append(hmmfam[seqid])
    else:
        fam.append("UNKNOWN")

    if seqid in pssmfam:
        fam.append(pssmfam[seqid])
    else:
        fam.append("UNKNOWN")

    fam.append(definitive_prediction(fam[0], fam[1]))

    return fam


def get_pssm_fam(mdict: dict) -> dict:
    """get_pssm_fam return the family based on pssm.

    With the highest number of
    occurence in PSSM profile match recorded as list for each
    sequence id.

    >>> my_dict = {ID1: defaultdict(<class 'list'>,
                                    { 'A' : ['SIG', 'MAT']},
                                    {'B': ['MAT']}
                                    ),
                   ID2: defaultdict(<class 'list'>,
                                    {'M': ['MAT']},
                                    {'P': ['MAT']},
                                    {'O1': ['PRO', 'MAT']}
                                    )
                   }
    >>> get_pssm_fam(my_dict)
    {ID1: 'A', ID2: 'O1'}

    :mdict: Dictionnary, required (dict)
    """
    fam = ""
    pssmfam = {}
    for key in mdict:
        x = Counter(mdict[key])
        # Take the top 2 item with highest count in list
        possible_fam = x.most_common(2)

        if len(possible_fam) == 1:
            fam = possible_fam[0][0]
        elif len(possible_fam) > 1:
            if len(possible_fam[0][1]) == len(possible_fam[1][1]):  # type: ignore[type-supported]
                fam = f"CONFLICT {possible_fam[0][0]}" + f" and {possible_fam[1][0]}"
            elif len(possible_fam[0][1]) > len(possible_fam[1][1]):  # type: ignore[type-supported]
                fam = possible_fam[0][0]
            else:
                fam = possible_fam[1][0]

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

    :param input_file: Path to the input FASTA file.
    :param max_records: Maximum number of records per split file.
    :return: Path to the directory containing the split files.
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
    """Isdnaorproteins test if input sequence is DNA or proteins.

    s input sequence
    """
    dna = "ATCG"
    prot = "ABCDEFGHIKLMNPQRSTVWYZ*X"
    stype = ""

    if all(i in dna for i in s):
        stype = "DNA"
    elif all(i in prot for i in s):
        stype = "protein"
    else:
        stype = "unknown"

    return stype


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


def get_hmm_superfamily(result_with_combined_evalue: dict) -> dict:
    """Compare combined_evalues to get superfamilies."""
    seqfam = {}
    for key, (fam_list, max_combined_evalue) in result_with_combined_evalue.items():
        res = ""
        conflicts = [
            fam_dict["fam"]
            for fam_dict in fam_list
            if fam_dict["combined_evalue"] == max_combined_evalue
        ]
        if len(conflicts) > 1:
            res = f"Conflict between {', '.join(conflicts)}"
        else:
            res = conflicts[0]
        seqfam[key] = res
    return seqfam


def compute_combined_evalue(transformed_result: dict) -> dict:
    """Compute the combined evalue for each fam and find the maximum combined evalue."""
    results = {}
    for key, fam_list in transformed_result.items():
        max_combined_evalue = None
        for fam_dict in fam_list:
            combined_evalue = reduce(
                operator.mul,
                (info["evalue"] for info in fam_dict["info"]),
                1,
            )
            fam_dict["combined_evalue"] = combined_evalue
            if max_combined_evalue is None or combined_evalue > max_combined_evalue:
                max_combined_evalue = combined_evalue
        results[key] = fam_list, max_combined_evalue
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

    `infile` Pyfasta object.
    `outfile` Output file.
    `sw` Sequence width. Default: 60.
    """
    fa = pyfastx.Fasta(infile)
    output_path = Path(f"{outfile}_allpep.fa")
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

    :seq: DNA sequence to translate.
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
