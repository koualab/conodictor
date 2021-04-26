from Bio.Seq import reverse_complement, translate
from collections import Counter, defaultdict
from decimal import Decimal
from functools import reduce
from heapq import nsmallest
from operator import mul
import pathlib
import pyfastx
import re
import warnings


def clear_dict(hdict):
    """
    clear_dict filter out sequences without a MATURE HMM or PSSM profile
    matched by hmmsearch or pfscan. It return the input dict with sequences
    without MATURE profile match filtered out.

    :hdict: A dictionnary containing matching profiles names by family by
            sequence id.

    example: {'sp|P0C640|CT55_CONPL':
              defaultdict(<class 'list'>, {'A': ['MAT', 'SIG']}),
              'sp|Q1A3Q6|CT57_CONLT':
               defaultdict(<class 'list'>, {'T': ['MAT', 'PRO', 'SIG']})}

    Such dict is created with defaultdict(lambda: defaultdict(list))
    """
    remove = defaultdict(list)

    # Get list of seq without mature sequence match
    for k in hdict.keys():
        for a, b in hdict[k].items():
            if "MAT" not in b:
                remove[k].append(a)

    # Remove families without mature sequence match
    for k, v in remove.items():
        for x in v:
            del hdict[k][x]

    # Remove sequence with no match with mature sequence
    for k in remove.keys():
        if not hdict[k]:
            del hdict[k]

    return hdict


def definitive_prediction(hmmclass, pssmclass):
    """
    definitive_prediction gives definitive classification by
    combining HMM and PSSM classification.

    :hmmclass: HMM predicted family, required (string)
    :pssmclass: PSSM predicted family, required (string)
    """

    deffam = None

    if hmmclass == pssmclass:
        deffam = hmmclass
    elif "CONFLICT" in pssmclass and "CONFLICT" in hmmclass:
        fams_pssm = re.search("(?<=CONFLICT)(.*)and(.*)", pssmclass)
        fams_hmm = re.search("(?<=CONFLICT)(.*)and(.*)", hmmclass)
        deffam = f"CONFLICT {fams_pssm.group(1)}, {fams_pssm.group(2)},"
        +f" {fams_hmm.group(1)}, and {fams_hmm.group(2)}"
    elif "CONFLICT" in pssmclass and "CONFLICT" not in hmmclass:
        deffam = hmmclass
    elif "CONFLICT" in hmmclass and "CONFLICT" not in pssmclass:
        deffam = pssmclass
    elif pssmclass != hmmclass:
        deffam = f"CONFLICT {hmmclass} and {pssmclass}"

    return deffam


def isdnaorproteins(s):
    """
    isdnaorproteins test if input sequence is DNA or proteins.

    :s: input sequence
    """

    dna = "ATCG"
    prot = "ABCDEFGHIKLMNPQRSTVWYZ*"
    stype = ""

    if all(i in dna for i in s):
        stype = "DNA"
    elif all(i in prot for i in s):
        stype = "protein"
    else:
        stype = "unknown"

    return stype


def get_pssm_fam(mdict):
    """
    get_pssm_fam return the family with the highest number of
    occurence in PSSM profile match recorded as list for each
    sequence id.

    >>> my_dict = {ID1: ['A', 'A', 'B', 'M'], ID2: ['M', 'P', 'O1', 'O1']}
    >>> get_pssm_fam(my_dict)
    {ID1: 'A', ID2: 'O1'}

    :mdict: Dictionnary, required (dict)
    """

    fam = ""
    pssmfam = {}
    for key in mdict.keys():
        x = Counter(mdict[key])
        # Take the top 2 item with highest count in list
        possible_fam = x.most_common(2)

        if len(possible_fam) == 1:
            fam = possible_fam[0][0]
        elif len(possible_fam) > 1:
            if possible_fam[0][1] == possible_fam[1][1]:
                fam = (
                    f"CONFLICT {possible_fam[0][0]} and {possible_fam[1][0]}"
                )
            elif possible_fam[0][1] > possible_fam[1][1]:
                fam = possible_fam[0][0]
            else:
                fam = possible_fam[1][0]

        pssmfam[key] = fam

    return pssmfam


def hmm_threshold(mdict):
    """
    hmm_threshold calculate evalue by family for each sequence
    and return a dict with the evalue for each family.

    :mdict: Dictionnary, required (dict)
    """

    score = defaultdict(dict)
    for key in mdict.keys():
        for k, v in mdict[key].items():
            # v has the format evalue|hsp_start#hsp_end
            score[key][k] = reduce(
                mul, [Decimal(x.split("#")[0]) for x in v], 1
            )

    return score


def pssm_sorter_func(x):
    seq, reg = x.split("#")
    return reg


def get_pssm_seq(pdict, id, pclass):
    """
    get_pssm_seq returns sequence matched by PSSM profile.
    It take as input the pssmdict created from pfscan output and
    the id of the sequence wanted along with the final PSSM classification.

    :pdict: PSSM output file created from pfscan output.
    :id: sequence id to get matched sequence.
    :pclass: PSSM classification of desired sequence.
    """

    pseq = []

    for k in pdict.keys():
        for x, j in pdict[k].items():
            if k == id and x == pclass:
                v = sorted(
                    j,
                    key=pssm_sorter_func,
                    reverse=True,
                )
                pseq.extend([x.split("#")[0] for x in v])

    return "".join(pseq)


def get_hmm_seq(mydict, mykey, file):
    """
    get_hmm_seq returns sequence matched by HMM profile.
    It take as input the hmmdict created from hmmsearch output and
    the id of the sequence to retrieve sequence from. It takes also
    the original fasta file to retrieve sequence from using sequence index.

    :mydict: Hmmdict created from hmmsearch output
    :mykey: Sequence id
    :file: input original file
    """

    hsp_id = {}
    hspseq = []
    ml = []

    for k in mydict.keys():
        for mv in mydict[k].values():
            if k == mykey:
                for x in mv:
                    a = x.split("#")[1]
                    ml.extend(a.split("|"))
                    fi = [int(b) for b in ml]
                hsp_id[k] = fi

    for k, v in hsp_id.items():
        v = sorted(v)
        hspseq = [(file[k].seq)[v[0] : v[-1]]]  # noqa

    return "".join(hspseq)


def get_hmm_fam(mdict):
    """
    get_hmm_fam get sequence family from hmm dictionnary.

    :mdict: Dictionnary of evalues by families.
    """

    conofam = ""
    seqfam = {}
    for key in mdict.keys():
        two_smallest = nsmallest(2, mdict[key].values())

        if len(two_smallest) == 1:
            conofam = next(iter(mdict[key]))
        elif two_smallest[0] * 100 != two_smallest[1]:
            conofam = list(mdict[key].keys())[
                list(mdict[key].values()).index(two_smallest[0])
            ]
        elif two_smallest[0] * 100 == two_smallest[1]:
            fam1 = list(mdict[key].keys())[
                list(mdict[key].values()).index(two_smallest[0])
            ]
            fam2 = list(mdict[key].keys())[
                list(mdict[key].values()).index(two_smallest[1])
            ]
            conofam = f"CONFLICT {fam1} and {fam2}"

        seqfam[key] = conofam

    return seqfam


def _translate_seq(seq):
    """
    _translate_seq translate DNA sequence to proteins in the six frames.

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


def do_translation(infile, outfile, sw=60):
    """
    do_translation translate a DNA fasta file into proteins
    fasta file.

    :infile: Input DNA fasta file.
    :outfile: Output file.
    :sw: Sequence width. Default: 60.
    """

    seqin = pyfastx.Fasta(infile)
    with open(pathlib.Path(f"{outfile}_allpep.fa"), "w") as protfile:
        for sequence in seqin:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                protseq = _translate_seq(sequence.seq)
                for idx, frame in enumerate(protseq):
                    # Rule E203 from flacke8 check for extraneous whitespace
                    # before a colon. But black follow PEP8 rules.
                    # A PR is open to resolve this issue:
                    # https://github.com/PyCQA/pycodestyle/pull/914
                    seq_letters = [
                        frame[i : i + sw]  # noqa: E203
                        for i in range(0, len(frame), sw)
                    ]
                    nl = "\n"
                    protfile.write(
                        f">{sequence.name}_frame={idx + 1}\n"
                        + f"{nl.join(map(str, seq_letters))}\n"
                    )
