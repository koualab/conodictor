[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


<h1 align="center">
  ConoDictor
</h1>
<h3 align="center">
  Improved prediction and classification of conopeptides
</h3>

## Introduction

Cone snails are among the richest sources of natural peptides with promising pharmacological and therapeutic applications. With the reduced costs of RNAseq, scientists now heavily rely on venom gland transcriptomes for the mining of novel bioactive conopeptides, but the bioinformatic analyses often hamper the discovery process.

ConoDictor 2 is a standalone and user-friendly command-line program. We have updated the program originally published as a web server 10 years ago using novel and updated tools and algorithms and improved our classification models with new and higher quality sequences. ConoDictor 2 is now more accurate, faster, multiplatform, and able to deal with a whole cone snail venom gland transcriptome (raw reads or contigs) in a very short time. The new version of Conodictor also improves the identification and subsequent classification for entirely novel or relatively distant conopeptides.

The only input ConoDictor 2 requires is the assembled transcriptome or the raw reads file either in DNA or amino acid: used alphabet is automatically recognized. ConoDictor 2 run predictions directly on the proteins file (submitted or dynamically generated) and tries to report the longest conopeptide  precursor-like sequence.

## Installation

### Docker

Accessible at https://hub.docker.com/u/ebedthan or on [BioContainers](https://github.com/BioContainers/containers/tree/master/conodictor/2.2.2).

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
wget https://github.com/koualab/conodictor/releases/tag/v2.3.0

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
  _____                  _____  _      _             
 / ____|                |  __ \(_)    | |            
| |     ___  _ __   ___ | |  | |_  ___| |_ ___  _ __ 
| |    / _ \| '_ \ / _ \| |  | | |/ __| __/ _ \| '__|
| |___| (_) | | | | (_) | |__| | | (__| || (_) | |   
 \_____\___/|_| |_|\___/|_____/|_|\___|\__\___/|_|    v2.3.0

conodictor [FLAGS/OPTIONS] <file>
Examples:
	conodictor file.fa.gz
	conodictor --out outfolder --cpus 4 --mlen 51 file.fa

positional arguments:
  file                  Specifies input file.

optional arguments:
  -h, --help            show this help message and exit
  -o OUT, --out OUT     Specify output folder.
  --mlen MLEN           Set the minimum length of the sequence to be
                        considered as a match
  --ndup NDUP           Minimum sequence occurence of a sequence to be
                        considered
  --faa                 Create a fasta file of matched sequences. Default:
                        False.
  --filter              Activate the removal of sequences that matches only
                        the signal and/or proregions. Default: False.
  -a, --all             Display sequence without hits in output. Default:
                        False.
  -j CPUS, --cpus CPUS  Specify the number of threads. Default: 1.
  --force               Force re-use output directory. Default: Off.
  -q, --quiet           Decrease program verbosity
  --debug               Activate debug mode
```


## Invoking conodictor

```
conodictor file.fa.gz
conodictor --out outfolder --cpus 4 --mlen 51 file.fa
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
         file              Specify input fasta file [required]

Outputs:
         -o, --out         Specify output folder.
         --faa             Create a fasta file of matched sequences. Default: False.
         -a, --all         Display sequence without hits in output. Default: False.
         --force           Force re-use output directory. Default: Off.
Computation:
         -j, --cpus        Specify number of threads. Default: 1.
         
Setup:
         -q, --quiet       Decrease verbosity
         --debug           Activate debug mode

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