[![PyPI](https://img.shields.io/pypi/v/conodictor.svg)](https://pypi.org/project/conodictor)
[![Wheel](https://img.shields.io/pypi/wheel/conodictor.svg)](https://pypi.org/project/conodictor)
[![Language](https://img.shields.io/pypi/implementation/conodictor)](https://pypi.org/project/conodictor)
[![Pyver](https://img.shields.io/pypi/pyversions/conodictor.svg)](https://pypi.org/project/conodictor)
[![Downloads](https://img.shields.io/pypi/dm/conodictor)](https://pypi.org/project/conodictor)
[![Docker](https://img.shields.io/docker/pulls/ebedthan/conodictor.svg)]()
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## ConoDictor: A fast and accurate prediction and classification tool for conopeptides


### Important
If using conodictor and have issue like [CONODB issue](https://github.com/koualab/conodictor/issues/18), please update to v2.3.6 that provide a fix.

conodictor v2.3.6 introduce the -d option to specify the path to the db folder containing HMM and PSSM files for classification.

This is a temporary solution while I am working on the next big release. Thanks.

### Introduction

Cone snails are among the richest sources of natural peptides with promising pharmacological and therapeutic applications. With the reduced costs of RNAseq, scientists now heavily rely on venom gland transcriptomes for the mining of novel bioactive conopeptides, but the bioinformatic analyses often hamper the discovery process.

ConoDictor 2 is a standalone and user-friendly command-line program. We have updated the program originally published as a web server 10 years ago using novel and updated tools and algorithms and improved our classification models with new and higher quality sequences. ConoDictor 2 is now more accurate, faster, multiplatform, and able to deal with a whole cone snail venom gland transcriptome (raw reads or contigs) in a very short time.

The only input ConoDictor 2 requires is the assembled transcriptome or the raw reads file either in DNA or amino acid: the used alphabet is automatically recognized. ConoDictor 2 runs predictions directly on the proteins file (submitted or dynamically generated) and tries to report the longest conopeptide precursor-like sequence.

### Installation

#### Install from Pip

You will have first to install [HMMER 3](https://hmmer.org) and [Pftools](https://github.com/sib-swiss/pftools3) to be able to run conodictor.

```bash
pip install conodictor
```

#### Using containers

#### Docker

Accessible at https://hub.docker.com/u/ebedthan or on [BioContainers](https://github.com/BioContainers/containers/tree/master/conodictor/2.2.2).


```bash
docker pull ebedthan/conodictor:latest
docker run ebedthan/conodictor:latest conodictor -h
```

Example of a run

```bash
docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) ebedthan/conodictor:latest conodictor --out /data/outdir /data/input.fa.gz
```

See https://staph-b.github.io/docker-builds/run_containers/ for more informations on how to properly run a docker container.


#### Singularity

The singularity container does not need admin privileges making it
suitable for university clusters and HPC.

```bash
singularity build conodictor.sif docker://ebedthan/conodictor:latest
singularity exec conodictor.sif conodictor -h
```


#### Install from source

```bash
# Download ConoDictor development version
git clone https://github.com/koualab/conodictor.git conodictor

# Navigate to directory
cd conodictor

# Install with poetry: see https://python-poetry.org
poetry install --no-dev

# Enter the Python virtual environment with
poetry shell

# Test conodictor is correctly installed
conodictor -h
```

If you do not want to go into the virtual environment just do:

```bash
poetry run conodictor -h
```

## Test

* Type `conodictor -h` and it should output something like:

```
usage: conodictor [options] <FILE>

optional arguments:
  -o DIR, --out DIR   output result to DIR [ConoDictor]
  --mlen INT          minimum length of sequences to be considered [off]
  --ndup INT          minimum occurence sequences to be considered [off]
  --faa               dump a fasta file of matched sequences [false]
  --filter            only keep sequences matching sig, pro and mat regions [false]
  -a, --all           add unclassified sequences in result [false]
  -j INT, --cpus INT  number of threads [1]
  --force             re-use output directory [false]
  -q, --quiet         decrease program verbosity
  -v, --version       show program's version number and exit
  -h, --help          show this help message and exit

Citation: Koua et al., 2021, Bioinformatics Advances
```


## Invoking conodictor

```bash
conodictor file.fa.gz
conodictor --out outfolder --cpus 4 --mlen 51 file.fa
```
  

## Output files

The comma separeted-values file summary.csv can be easily viewed with any office suite,
or text editor.

```csv
sequence,hmm_pred,pssm_pred definitive_pred
SEQ_ID_1,A,A,A
SEQ_ID_2,B,D,CONFLICT B and D
SEQ_ID_3,O1,O1,O1
...

```

## Citation

When using ConoDictor2 in your work, you should cite:

Dominique Koua, Anicet Ebou, Sébastien Dutertre, Improved prediction of conopeptide superfamilies with ConoDictor 2.0, Bioinformatics Advances, Volume 1, Issue 1, 2021, vbab011, https://doi.org/10.1093/bioadv/vbab011.
  
## Bugs

Submit problems or requests to the [Issue Tracker](https://github.com/koualab/conodictor/issues).


## Dependencies

### Mandatory

* [**HMMER 3**](https://hmmer.org)  
  Used for HMM profile prediction.   
  *Eddy SR, Accelerated Profile HMM Searches. PLOS Computational Biology 2011, 10.1371/journal.pcbi.1002195*

* [**Pftools**](https://github.com/sib-swiss/pftools3)  
  Used for PSSM prediction.    
  *Schuepbach P et al. pfsearchV3: a code acceleration and heuristic to search PROSITE profiles. Bioinformatics 2013, 10.1093/bioinformatics/btt129*


## Licence

[GPL v3](https://github.com/koualab/conodictor/blob/main/LICENSE).

For commercial uses please contact Dominique Koua at dominique.koua@inphb.ci.

## Authors

* [Anicet Ebou](https://orcid.org/0000-0003-4005-177X)
* [Dominique Koua](https://www.researchgate.net/profile/Dominique_Koua)