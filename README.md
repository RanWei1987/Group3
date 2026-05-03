# Group3
Repository for APHL Hackathon Group3
## Project Outline
This project is a bioinformatics pipeline that assembles and analyzes the the genome of *Legionella pneumophila* using Illumina or ONT reads. The core of the analysis module is **el_gato**, a bioinformatics tool that derives Sequence Type (ST) information.  
The workflow of this project:  
Illumina: el_gato  
ONT: NanoPlot -> Seqtk (Subsample) -> RAVEN (De novo Assemble) -> Medaka (Polishing) -> Prokka (Annotation) -> AmrFinder -> el_gato
## Installation
### Prerequisite
Hardware: 8GB RAM, 2 CPUs
Software: Nextflow, docker

## Codespace setup
See ```.devcontainer``` files, ```devcontainer.json``` and ```setup.sh```

## Local Setup
### Prerequisites
- Python 3.11+
- Java 11+ (`sudo apt install default-jre` / `brew install java`)
- Docker ([install guide](https://docs.docker.com/engine/install/)) - Docker Desktop or Docker Engine
- Graphviz (`sudo apt install graphviz` / `brew install graphviz`) - system binary needs to be installed before the pip installation
- Nextflow (`curl -s https://get.nextflow.io | bash`)
```
# Obtain the git project using either git clone or downloading the code manually. 
git clone https://github.com/APHL-Infectious-Disease/Group3.git
cd Group3

# Install Python dependencies for the dashboard
pip install -r scripts/requirements.txt

## Install Nextflow if needed
# curl -s https://get.nextflow.io | bash
#sudo mv nextflow /usr/local/bin/

# Example command to run the pipeline
nextflow run scripts/el_gato_run.nf --read1 /path/to/your/ont_reads.fastq.gz --sample_id your_sample
```

If using the ONT option, the workflow will automatically build and launch a Streamlit app to view the results interactively. Users can upload their own files to the app interface if any are missing.

The devcontainer setup will install Python and necessary dependencies for the dashboard. Note: if running in a GitHub Codespace, ```streamlit run``` will work BUT the URL printed to terminal isn't clickable and won't let you navigate to localhost. Codespaces rewrites the app to the dashboard URL ```*.app.github.dev``` (which will appear in the Ports tab) and the "onAutoForward": "openBrowser" line in setup enables port forwarding to your browser. 


### Run examples
Illumina:
```
nextflow scripts/el_gato_run.nf --read1 data/illumina/sample01/ERR10879100_1.fastq.gz --read2 data/illumina/sample01/ERR10879100_2.fastq.gz --sample_id ERR10879100
```
ONT:
```
nextflow scripts/el_gato_run.nf --read1 data/ont/sample01/ERR10879072.fastq.gz --sample_id ERR10879072
```
