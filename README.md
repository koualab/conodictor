[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# ConoDictor: improved prediction and classification of conopeptides

## Introduction

Conopeptides are the main component of cone snails venom. It has proven to have diverse pharmaceutical, physiological and therapeutic application. In order to accelerate the drug
discovery process, ConoDictor predict and classifiy amino acid precursors of conopeptides using
hidden Markov models and position-specific scoring matrix. 

## Installation

### Unix-like operating systems 

```
# Download ConoDictor repository and change directory
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
pip3 install -r conodictor/requirements.txt

# Restart terminal or source .bashrc or .bash_profile
# Test conodictor is correctly installed
conodictor -h
```

## Test

* Type `conodictor -h` and it should output something like:

```
usage: conodictor [options] seqs.fa.gz

positional arguments:
  seqs         Specify input sequences fasta file

optional arguments:
  -h, --help   show this help message and exit
  --out OUT    Specify Output directory
  --all        Display unpredicted sequence in output
  --cpus CPUS  Specify the number of threads
  --force      Force re-use output directory
  --quiet      Decrease program verbosity
  --debug      Activate debug mode

Version:   2.1.3
Licence:   GPL-3
Homepage:  https://github.com/koualab/conodictor.git
Author:    Anicet Ebou <anicet.ebou@gmail.com>
Last Run: Tue, 02 Feb 2021 23:33:00.
```


## Invoking conodictor

```
conodictor seqs.fa.gz
```
  

## Output files

```
summary.txt


sequence  hmm_pred  pssm_pred definitive_pred
SEQ_ID_1  A A A
SEQ_ID_2  B D CONFLICT B and D
SEQ_ID_3  O1  O1  O1
...

```

## Command line options

```
General:
         seqs          Specify input fasta file [required]

Outputs:
         --out         Specify output folder name
         --all         Add unpredicted sequence in output.
         --force       Force reuse of output folder

Computation:
         --cpus        Specify number of threads
         
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