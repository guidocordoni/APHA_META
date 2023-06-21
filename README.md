# APHA_META V.1.3
Metagenomic analysis of Illumina sequences.

1) Use the apha_met_gc_mamba_installer.sh to install all software needed by the pypeline to run
2) activate apha_meta mamba environment.
3) create a folder called RAW_READS and copy all your metagenomic sequences there.
4) Copy the APHA_META_GC_V1.2_working.py file in the folder that contains the folder RAW_READS thay tou have just created and run it (python APHA_META_GC_V1.2_working.py)

What this pipeline will do?
1) Cleans the sequences using fastP 
2) perform the metagenomic assembly using Megahit
3) Runs Kraken2 on assemblies and generate a Krona graph
4) Bins the sequences in your assemblies according to Kraken2 taxonomy. You will get multi-fasta files containing all the sequences of a specific organism and the files will be named after the Kraken2 taxonomy)

The script intentionally stops here leaving you the possibilty to use the results for further metagenomic or comparative genomics analysis.

It's suggested to use at least a 16 cores machine with 128GB memory. It will not stop running if you are using a 4 core machine but it can take a very long time to get the results.
REMEMBER TO CHANGE THE NUMBER OF THREADS VALUE ACCORDING TO YOUR CORE AVAILABILITY (default is set to 16 cores)

SCRIPT UNDER DEVELOPMENT
KNOWN ISSUES:
No known issue
