[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# ConoDictor: improved prediction and classification of conopeptides

## Introduction

Conopeptides are the main component of cone snails venom. It has proven to have diverse pharmaceutical, physiological and therapeutic application. In order to accelerate the drug
discovery process, ConoDictor predict and classifiy amino acid precursors of conopeptides using
hidden Markov models and position-specific scoring matrix. 

## Installation

### Docker

Accessible at https://hub.docker.com/u/ebedthan

```
docker pull ebedthan/conodictor:latest
docker run ebedthan/conodictor:latest conodictor -h

# Example of a run
docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) ebedthan/conodictor:latest conodictor --out /data/outdir /data/input.fa.gz

```

See https://staph-b.github.io/docker-builds/run_containers/ for more informations about the run.

### Singularity

```
singularity build conodictor.sif docker://ebedthan/conodictor:latest
singularity exec conodictor.sif conodictor -h
```


### Unix-like operating systems 

```
# Download ConoDictor latest version
wget https://github.com/koualab/conodictor/releases/tag/v2.2.0

# or Download ConoDictor development version
git clone https://github.com/koualab/conodictor.git conodictor

cd conodictor

# Copy conodictor bin folder: Requires admin privileges
sudo cp conodictor /usr/local/bin

# Copy conodictor databases: Requires admin privileges
sudo mkdir /usr/share/conodictor/db
sudo cp db/* /usr/share/conodictor/db

# Add path to databases in .bashrc and/or .bash_profile
export CONODB=/usr/share/conodictor/db

# Create new python environment conoenv
cd ..
python3 -m venv conoenv
python3 -m pip install -r conodictor/requirements.txt

# Restart terminal or source .bashrc or .bash_profile
# Test conodictor is correctly installed
conodictor -h
```

## Test

* Type `conodictor -h` and it should output something like:

```
usage: conodictor [options] seqs.fa.gz

positional arguments:
  seqs         Specify input fasta file.

optional arguments:
  -h, --help   show this help message and exit
  --out OUT    Specify output folder.
  --mlen MLEN  Minimum sequence length to consider for prediction.
  --ndup NDUP  Minimum number of sequence occurence for a sequence to be
               considered.
  --all        Display sequence without hits in output. Default: False.
  --cpus CPUS  Specify the number of threads. Default: 1.
  --force      Force re-use output directory. Default: Off.
  --quiet      Decrease program verbosity
  --debug      Activate debug mode

Version:   2.2.0
Licence:   GPL-3
Homepage:  https://github.com/koualab/conodictor.git
Author:    Anicet Ebou <anicet.ebou@gmail.com>
Last Run:  Sat, 06 Mar 2021 13:26:59.
```


## Invoking conodictor

```
conodictor seqs.fa.gz
```
  

## Output files

```
summary.txt


sequence  length  num_cysteines occurence hmm_pred  pssm_pred definitive_pred
SEQ_ID_1  56  4 2 A A A
SEQ_ID_2  60  0 1 B D CONFLICT B and D
SEQ_ID_3  145 8 1 O1  O1  O1
...

```

## Command line options

```
General:
         seqs          Specify input fasta file [required]

Outputs:
         --out         Specify output folder.
         --mlen        Minimum sequence length to consider for prediction.
         --ndup        Minimum number of sequence occurence for a sequence to be considered.
         --all         Display sequence without hits in output. Default: False.
         --force       Force re-use output directory. Default: Off.
Computation:
         --cpus        Specify number of threads. Default: 1.
         
Setup:
         --quiet       Decrease verbosity
         --debug       Activate debug mode

Standard meta-options:
         --help, -h    Print help and exit

```
  
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

[GPL v3](https://github.com/koualab/conodictor/blob/main/LICENSE)

## Authors

* [Anicet Ebou](https://orcid.org/0000-0003-4005-177X)
* [Dominique Koua](https://www.researchgate.net/profile/Dominique_Koua)