#!/bin/bash
# Made by Guido Cordoni guido.cordoni@apha.gov.uk
# Create conda environment
mamba create -n apha_meta python=3.6

# Activate the environment
mamba activate apha_meta
#pip install
pip install pyfiglet
pip install biopython
pip install pysam
pip install pandas

# Install required software and dependencies
mamba install -c bioconda fastp multiqc trimmomatic megahit spades metabat2 pplacer checkm-genome kraken kraken2 krakentools krona prokka

#sudo checkm data setRoot /home/guidocordoni/fsx/ranch-44/Other_projects/guido/db/checkm_db
# Deactivate the environment
mamba deactivate

echo "Remember to activate apha-meta environment before to run the analysis (mamba activate mamba_meta)

