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
### Procedures
Obtain the git project using either git clone or downloading the code manually. Then run scripts/download_database.sh to download the database required by CheckM2
### Run example
Illumina:
```
nextflow scripts/el_gato_run.nf --read1 data/illumina/sample01/ERR10879100_1.fastq.gz --read2 data/illumina/sample01/ERR10879100_2.fastq.gz --sample_id ERR10879100
```
ONT:
```
nextflow scripts/el_gato_run.nf --read1 data/ont/sample01/ERR10879072.fastq.gz --sample_id ERR10879072
```
