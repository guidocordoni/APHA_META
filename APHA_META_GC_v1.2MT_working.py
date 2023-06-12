import os
import glob
import multiprocessing
import pandas as pd
import sys
from Bio import SeqIO
import pyfiglet
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

graphic = pyfiglet.figlet_format("APHA META GC V 1 . 2")
print(graphic)
print("APHA_META_GC is a pipeline to work on metagenomic sequences. This version can handle Illumina generated sequences and can detect viruses and bacteria.\n\nBefore starting you MUST create a folder called RAW_READS and put the sequences that you want to analyse there.\n\nThe pipeline will perform sequences cleaning and quality check using fastp, metagenomic assembly using megahit, Kraken 2 to assess the taxonomy of the organisms present in your samples.\n\nThe last part of this pipeline will bin the contigs obtained with megahit according to Kraken2 taxonomy found (each file will be named according to the taxonomy.\n\nIt's a long process (depending on how many sequences you have can go from hours to days), please be patient and let the computer work for you.\n\nThis pipeline intentionally stops after the binning to give you the possibility to choose what to do with the data generated (i.e. abundance stats, comparative genomics of some selected groups, intra/cross-domain associations). \n\nBefore starting, install the necessary mamba environment using the bash script apha_met_gc_mamba_installer.sh provided. ***REMEMBER TO COPY THIS SCRIPT IN THE FOLDER WHERE YOU HAVE CREATED YOUR RAW_READS FOLDER***\n\n")

# Check if mamba environment is installed
mamba_installed = os.system("mamba --version >/dev/null 2>&1") == 0

if not mamba_installed:
    print("Mamba environment is not installed. Please install it before running this script.")
    sys.exit(1)

# Check if mamba environment apha-meta is activated
env_activated = "CONDA_DEFAULT_ENV" in os.environ and os.environ["CONDA_DEFAULT_ENV"] == "apha_meta"

if not env_activated:
    print("Mamba environment 'apha_meta' is not activated. Please activate it before running this script.")
    sys.exit(1)

# Check if RAW_READS folder exists
raw_reads_folder_exists = os.path.isdir("RAW_READS")

if not raw_reads_folder_exists:
    print("RAW_READS folder does not exist. Please create the RAW_READS folder and put your sequences there.")
    sys.exit(1)

# Ask the user if they want to proceed with the metagenomic analysis
user_input = input("Do you want to proceed with the metagenomic analysis? (Y/N) ")

while user_input not in ("Y", "y", "N", "n"):
    user_input = input("Please enter Y or N. ")

if user_input in ("Y", "y"):
    # Variable to store the user's input
    proceed = user_input

# Check if the user wants to proceed
if proceed == "Y" or proceed == "y":
    print("Running the metagenomic analysis...")
    # Folders creation. This bit will create all the folders needed for the software to run. Please put your samples in the RAW_READS folder, all the other folders will be populated automatically, and the subfolders in them will have the name of your file.
    os.makedirs("READ_QC", exist_ok=True)
    os.makedirs("CLEAN_READS", exist_ok=True)
    os.makedirs("ASSEMBLY", exist_ok=True)
    os.makedirs("KRAKEN", exist_ok=True)

    # Cleaning and Quality Check
    def run_fastp(f):
        R = f.replace("_1.fq", "_2.fq")
        BASE = os.path.basename(f)
        SAMPLE = BASE[:-16]  # Remove the "_1.fq" extension
        os.system(
            f'fastp -i "RAW_READS/{f}" -o "CLEAN_READS/{SAMPLE}_1.fq" -I "RAW_READS/{R}" -O "CLEAN_READS/{SAMPLE}_2.fq" -w 16 -V -j "READ_QC/{SAMPLE}.json" -h "READ_QC/{SAMPLE}.html"'
        )

    with ThreadPoolExecutor() as executor:
        for f in os.listdir("RAW_READS"):
            if f.endswith("_R1_001.fastq.gz"):
                executor.submit(run_fastp, f)

    # ASSEMBLY + KRAKEN2
    clean_reads_folder = "CLEAN_READS"

    def run_megahit_kraken(f):
        R = f.replace("_1.fq", "_2.fq")
        BASE = os.path.basename(f)
        SAMPLE = BASE[:-5]  # Remove the "_1.fq" extension
        assembly_output_dir = os.path.join("ASSEMBLY", SAMPLE)
        KRAKEN_output_dir = os.path.join("KRAKEN", SAMPLE)
        KRAKEN_report_path = os.path.join(KRAKEN_output_dir, "output_report.txt")

        os.makedirs(assembly_output_dir, exist_ok=True)
        os.makedirs(KRAKEN_output_dir, exist_ok=True)

        os.system(f"megahit -1 {f} -2 {R} -t 16 -o {assembly_output_dir} -f")
        os.system(
            f'kraken2 --db /home/guidocordoni/fsx/ranch-44/Bactipipes/BactiPipes_SCE3/Programs/Kraken2/standard_db_2021-05-17 --threads 16 --memory-mapping --use-names --output {KRAKEN_report_path} "{assembly_output_dir}/final.contigs.fa"'
        )

    with ThreadPoolExecutor() as executor:
        for f in glob.glob(os.path.join(clean_reads_folder, "*_1.fq")):
            executor.submit(run_megahit_kraken, f)

    # Taxonomic binning based on Kraken2 output
    # Input Kraken2 output files
    kraken_output_files = glob.glob("KRAKEN/*/output_report.txt")
    # Input fasta file (contigs)
    fasta_file = glob.glob("ASSEMBLY/*/final.contigs.fa")

    # Output directory for bins
    bins_output_dir = "BINS"

    # Number of processors for parallel processing (set to maximum CPU count)
    num_processors = multiprocessing.cpu_count()

    # Create output directory if it doesn't exist
    os.makedirs(bins_output_dir, exist_ok=True)

    def process_kraken_output(kraken_output_file):
        # Extract the sample name from the Kraken2 output file path
        sample_name = os.path.basename(os.path.dirname(kraken_output_file))

        # Get the corresponding fasta file
        fasta_file = os.path.join("ASSEMBLY", sample_name, "final.contigs.fa")

        # Read the Kraken2 output file using Pandas
        df = pd.read_csv(kraken_output_file, sep="\t", header=None)

        # Check if column 2 exists in the DataFrame
        if 2 in df.columns:
            # Group sequences by taxonomy
            grouped = df.groupby(2)

            # Iterate over each taxonomy group and write sequences to the bin files
            for taxonomy, group in grouped:
                # Convert taxonomy to string and replace problematic characters in the taxonomy name
                taxonomy_str = str(taxonomy)
                taxonomy_filename = taxonomy_str.replace("/", "_").replace(":", "_")
                bin_dir = os.path.join(bins_output_dir, sample_name)
                os.makedirs(bin_dir, exist_ok=True)
                bin_file = os.path.join(bin_dir, f"{taxonomy_filename}.fa")

                # Get the sequence IDs belonging to the taxonomy group
                sequence_ids = group[1]

                # Open the fasta file and write sequences to the bin file
                with open(fasta_file, "r") as fa:
                    with open(bin_file, "w") as bin_fa:
                        for record in SeqIO.parse(fa, "fasta"):
                            record_id = record.id
                            if any(seq_id in record_id for seq_id in sequence_ids):
                                SeqIO.write(record, bin_fa, "fasta")

    with ThreadPoolExecutor() as executor:
        executor.map(process_kraken_output, kraken_output_files)

    # Function to process a single output_report.txt file and generate a Krona graph
    def process_output_report(file_path):
        # Directory containing the output_report.txt file
        output_dir = os.path.dirname(file_path)

        # Dictionary to store the counts for each taxon
        taxon_counts = defaultdict(int)

        # Dictionary to store the taxonomic hierarchy
        taxon_hierarchy = {}

        # Read the input file
        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith("C"):
                    fields = line.split("\t")
                    taxon = fields[2].split("(")[0].strip()  # Exclude specific taxonomic information
                    count = int(fields[3])
                    taxon_counts[taxon] += count
                    # Add taxonomic hierarchy information
                    hierarchy = fields[4]
                    if hierarchy:
                        taxon_hierarchy[taxon] = hierarchy

        # Write the Krona output to a file
        krona_output_path = os.path.join(output_dir, "krona_output.txt")
        with open(krona_output_path, "w") as file:
            for taxon, count in taxon_counts.items():
                file.write(f"{count}\t{taxon}\n")

        # Generate the Krona graph
        krona_graph_path = os.path.join(output_dir, "krona_graph.html")
        subprocess.run(["ktImportText", krona_output_path, "-o", krona_graph_path])

        # Modify the Krona graph file to include taxonomic hierarchy
        with open(krona_graph_path, "r") as file:
            krona_html = file.read()

        # Modify the HTML to include taxonomic hierarchy information
        for taxon, hierarchy in taxon_hierarchy.items():
            krona_html = krona_html.replace(f">{taxon}<", f">{taxon}<br><span class='subtaxa'>{hierarchy}</span><")

        # Add custom grouping for viruses and bacteria
        krona_html = krona_html.replace("Viruses", "Viruses<br><span class='subtaxa'>Virus Taxon</span>")
        krona_html = krona_html.replace("Bacteria", "Bacteria<br><span class='subtaxa'>Bacteria Taxon</span>")

        # Write the modified Krona graph HTML
        with open(krona_graph_path, "w") as file:
            file.write(krona_html)

    # Set the path to the directory containing the output_report.txt files
    directory = "KRAKEN/"

    # Find all output_report.txt files recursively in the directory
    file_paths = glob.glob(os.path.join(directory, "*", "output_report.txt"))

    # Process each output_report.txt file and generate a Krona graph
    with ThreadPoolExecutor() as executor:
        executor.map(process_output_report, file_paths)

    print("Metagenomic analysis completed.")

