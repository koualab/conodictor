"""Sub-program for conodictor."""

from __future__ import annotations

import contextlib
import datetime
import operator
import warnings
from collections import defaultdict
from functools import reduce
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

import pyfastx
import xphyle
from Bio.Seq import reverse_complement, translate
from pyhmmer.hmmer import hmmsearch
from xphyle import paths

if TYPE_CHECKING:
    from pyhmmer.easel import Sequence, SequenceBlock, SequenceFile
    from pyhmmer.plan7 import HMMFile


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


def decompress_file(input_file: str) -> str:
    """Decompress a file and return the path to the uncompressed file."""
    with TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir, "conodictor_input.fa")

        # Handling input file supply
        if input_file == "<stdin>":
            with output_path.open("w", encoding="utf-8") as binput, xphyle.xopen(
                paths.STDIN,
                context_wrapper=True,
            ) as infile:
                for line in infile:
                    binput.write(line)  # type: ignore[type-supported]
        else:
            with output_path.open("w", encoding="utf-8") as binput, xphyle.xopen(
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
    """`do_translation` translate a DNA fasta file into proteins fasta file.

    `infile` Pyfasta object.
    `outfile` Output file.
    `sw` Sequence width. Default: 60.
    """
    fa = pyfastx.Fasta(infile)
    with Path.open(Path(f"{outfile}_allpep.fa"), "w") as protfile:
        for sequence in fa:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                protseq = _translate_seq(sequence.seq)
                for idx, frame in enumerate(protseq):
                    seq_letters = [frame[i : i + sw] for i in range(0, len(frame), sw)]
                    nl = "\n"
                    protfile.write(
                        f">{sequence.name}_frame={idx + 1}\n"
                        f"{nl.join(map(str, seq_letters))}\n",
                    )


def _translate_seq(seq: str) -> list[str]:
    """_translate_seq translate DNA sequence to proteins in the six frames.

    :seq: DNA sequence to translate.
    """
    seqlist = []
    # frame 1
    seqlist.append(translate(seq))
    # frame 2
    seqlist.append(translate(seq[1:]))
    # frame 3
    seqlist.append(translate(seq[2:]))
    # frame 4
    seqlist.append(translate(reverse_complement(seq)))
    # frame 5
    seqlist.append(translate(reverse_complement(seq)[1:]))
    # frame 6
    seqlist.append(translate(reverse_complement(seq)[2:]))

    return seqlist
