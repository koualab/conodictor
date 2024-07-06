# ConoDictor

*A fast and accurate prediction and classification tool for conopeptides*

[![PyPI](https://img.shields.io/pypi/v/conodictor.svg)](https://pypi.org/project/conodictor)
[![Wheel](https://img.shields.io/pypi/wheel/conodictor.svg)](https://pypi.org/project/conodictor)
[![Language](https://img.shields.io/pypi/implementation/conodictor)](https://pypi.org/project/conodictor)
[![Pyver](https://img.shields.io/pypi/pyversions/conodictor.svg)](https://pypi.org/project/conodictor)
[![Downloads](https://img.shields.io/pypi/dm/conodictor)](https://pypi.org/project/conodictor)
[![Docker](https://img.shields.io/docker/pulls/ebedthan/conodictor.svg)]()
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## 🗺️ Overview
### Unlocking the Potential of Cone Snail Venom
Cone snails are a treasure trove of natural peptides with immense pharmacological and therapeutic potential. The advent of affordable RNA sequencing (RNAseq) has revolutionized the mining of novel bioactive conopeptides from venom gland transcriptomes. However, the complexity of bioinformatic analyses often impedes the discovery process.

### Introducing ConoDictor 2
ConoDictor 2 is a standalone, user-friendly command-line tool designed to streamline the discovery of conopeptides. Building on a decade-old web server, we have significantly upgraded ConoDictor with modern tools and algorithms, and enhanced our classification models using new, high-quality sequences. The result is a program that is more accurate, faster, and compatible across multiple platforms.

### Key Features
* **Enhanced Accuracy and Speed**: ConoDictor 2 processes entire venom gland transcriptomes, whether from raw reads or assembled contigs, in record time.
* **Ease of Use**: The program requires only the assembled transcriptome or raw reads file, in either DNA or amino acid format. ConoDictor 2 automatically recognizes the alphabet used.
* **Advanced Prediction Capabilities**: It runs predictions directly on the submitted or dynamically generated proteins file, aiming to identify the longest conopeptide precursor-like sequences.

### Simplified Bioinformatics for Breakthrough Discoveries
With ConoDictor 2, researchers can bypass the intricate bioinformatic challenges and focus on uncovering the next generation of bioactive peptides from cone snail venom. Its robust performance and user-centric design make it an indispensable tool in venom research and drug discovery.

## Installing

### Install from Pip

You will first have to install ~~[HMMER 3](https://hmmer.org) and~~ [Pftools](https://github.com/sib-swiss/pftools3) to be able to run conodictor (**as of version 2.4, conodictor does not need hmmer anymore as it use the wonderful [pyhmmer](https://github.com/althonos/pyhmmer) library**).

```bash
pip install conodictor
```

### Using containers

### Docker

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


### Singularity

The singularity container does not need admin privileges making it
suitable for university clusters and HPC.

```bash
singularity build conodictor.sif docker://ebedthan/conodictor:latest
singularity exec conodictor.sif conodictor -h
```


### Install from source

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


## 💡 Example

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

## 💭 Feedback

### Issue tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/koualab/conodictor/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

## ⚖️ License

[GPL v3](https://github.com/koualab/conodictor/blob/main/LICENSE).

For commercial uses please contact Dominique Koua at dominique.koua@inphb.ci.

## 🔖 Citation

ConoDictor is a scientifc software, with a [published paper](https://doi.org/10.1093/bioadv/vbab011) in the [Bioinformatics Advances](https://academic.oup.com/bioinformaticsadvances) journal. Please cite this article if you are using it in an academic work, for instance as: 
Koua, D., Ebou, A., & Dutertre, S. (2021). Improved prediction of conopeptide superfamilies with ConoDictor 2.0. Bioinformatics Advances, 1(1), vbab011. https://doi.org/10.1093/bioadv/vbab011


## Dependencies

* [**Pftools**](https://github.com/sib-swiss/pftools3)  
  Used for PSSM prediction.    
  *Schuepbach P et al. pfsearchV3: a code acceleration and heuristic to search PROSITE profiles. Bioinformatics 2013, 10.1093/bioinformatics/btt129*


## 📚 References

* [**HMMER 3**](https://hmmer.org)  
  Used for HMM profile prediction.   
  *Eddy SR, Accelerated Profile HMM Searches. PLOS Computational Biology 2011, 10.1371/journal.pcbi.1002195*

* [**Pftools**](https://github.com/sib-swiss/pftools3)  
  Used for PSSM prediction.    
  *Schuepbach P et al. pfsearchV3: a code acceleration and heuristic to search PROSITE profiles. Bioinformatics 2013, 10.1093/bioinformatics/btt129*


## Authors

* [Anicet Ebou](https://orcid.org/0000-0003-4005-177X)
* [Dominique Koua](https://www.researchgate.net/profile/Dominique_Koua)