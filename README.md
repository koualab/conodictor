[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# ConoDictor: improved prediction and classification of conopeptides

## Introduction

Conopeptides are the main component of cone snails venom. It has proven to have diverse pharmaceutical, physiological and therapeutic application. In order to accelerate the drug
discovery process, ConoDictor predict and classifiy amino acid precursors of conopeptides using
hidden Markov models and position-specific scoring matrix. 

## Installation

### Bioconda


### Ubuntu/Debian/Mint

```
sudo apt-get update
pip3 install datetime
git clone https://github.com/koualab/conodictor.git
cd conodictor-main
sudo cp conodictor /usr/local/bin
sudo mkdir /usr/share/conodictor/db
sudo cp db/* /usr/share/conodictor/db
# Edit .bashrc and/or .bash_profile
export CONODB=/usr/share/conodictor/db
conodictor -h
```

### Centos/Fedora/RHEL

```
sudo yum update
pip3 install datetime
git clone https://github.com/koualab/conodictor.git
cd conodictor-main
sudo cp conodictor /usr/local/bin
sudo mkdir /usr/share/conodictor/db
sudo cp db/* /usr/share/conodictor/db
# Edit .bashrc and/or .bash_profile
export CONODB=/usr/share/conodictor/db
conodictor -h
```

### MacOS

```
pip3 install datetime
git clone https://github.com/koualab/conodictor.git
cd conodictor-main
sudo cp conodictor /usr/local/bin
sudo mkdir /usr/share/conodictor/db
sudo cp db/* /usr/share/conodictor/db
# Edit .bashrc and/or .bash_profile
export CONODB=/usr/share/conodictor/db
conodictor -h
```
  

## Test

* Type `itap -h` and it should output its help screen.
* Type `itap --man` and you should see a man page at the screen.
  


## Invoking conodictor

```
conodictor seqs.fa.gz
```
  

## Output files

```
conodictor_summary.txt


sequence id HMM prediction  PSSM prediction conoDictor definitive prediction
SEQ_ID_1    A   A   A
SEQ_ID_2    B   D   CONFLICT B and D
SEQ_ID_3    O1  O1  O1
...

```

## Command line options

```
General:
    seqs          Specify input fasta file [required]

Outputs:
	--out         Specify output folder name
	--force       Force reuse of output folder

Setup:
	--quiet       Decrease verbosity
    --debug       Activate debug mode

Standard meta-options:
	--usage, -u     Print program usage and exit
	--man           Print man page
	--help, -h      Print help and exit

```
  
  
## Bugs

Submit problems or requests to the [Issue Tracker](https://github.com/conodictor/issues).

  


## Dependencies

### Mandatory

* **HMMER 3**  
  Used for HMM profile prediction 
  *Eddy SR, Accelerated Profile HMM Searches. PLOS Computational Biology 2011, 10.1371/journal.pcbi.1002195*

* **Pftools**  
  Used for PSSM prediction.  
  *Schuepbach P et al. pfsearchV3: a code acceleration and heuristic to search PROSITE profiles. Bioinformatics 2013, 10.1093/bioinformatics/btt129*


## Licence

GPL v3

## Authors

* Anicet Ebou
* Dominique Koua