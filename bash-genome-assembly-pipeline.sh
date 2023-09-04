#!/bin/bash -l

# Assembly pipeline information:
# phyloFlash is used for contamination analysis of the raw reads.
# fastp is used for preprocessing of the raw reads prior to assembly.
# FastQC is used for quality control of the processed reads.
# SPAdes is used to assembly the preprocessed reads into contigs/scaffolds.
# QUAST is used to generate general assembly statistics like N50.
# BUSCO is used for genome completeness assessment with a chosen lineage.
# Whokaryote is used for classification of the scaffolds.
# MetaEuk is used for gene discovery.

# This pipeline was only tested with gzipped paired-end short reads in FASTQ format.
# It runs offline, the data required for phyloFlash and BUSCO must be provided.

# Author: Justin Teixeira Pereira Bassiaridis
# Date: 2023-08-07
# License: MIT

# Setup following conda environments before running:
# Environment named "pf" containing phyloFlash 3.4.1.
# Environment named "assembly" containing fastp 0.23.4, FastQC 0.12.1 and SPAdes 3.15.5.
# Environment named "analysis" containing QUAST 5.2.0, BUSCO 5.4.5 and MetaEuk 6.
# Environment named "whokaryote" containing Whokaryote 1.1.2.
# May work with other versions of the software.

####################################
# Basic setup, change values here! #
####################################
# Working directory
work_dir=""
# phyloFlash data directory
pf_data=""
# BUSCO data directory
busco_data=""
# MetaEuk protein reference file
protein_data=""
# Sample name
sample=""
# Number of threads to use
threads=""
# Amount of memory to use
memory=""
# Raw short-read files (must end with fastq.gz)
reads1=""
reads2=""
####################################

# Work directory setup
log_dir="$work_dir/Logs/$sample/"
mkdir -p "$log_dir"
builtin cd "$work_dir"

# Text function for displaying software name and version at start
text_function() {
    # Short repeat function from:
    # https://stackoverflow.com/questions/5349718/how-can-i-repeat-a-character-in-bash
    repeat() { printf "$1""%.s" $(eval "echo {1.."$(($2))"}"); }
    printf "\n\n"
    printf "$(repeat = $((${#1}*3)))" # Repeat "=" 3x length of version output
    printf "\n\n"
    printf "$(repeat " " $((${#1}*1)))$1" # Repeat " " up to length of version output
    printf "\n\n"
    printf "$(repeat = $((${#1}*3)))"
    printf "\n\n"
}

# phyloFlash for contamination analysis
phyloflash_function() {
    conda activate pf
    # phyloFlash setup
    sample_fixed="${sample//./_}" # Because phyloFlash does not accept library names with dots
    pf_out_dir="$work_dir/Analyses/phyloFlash/$sample/"
    mkdir -p "$pf_out_dir"
    builtin cd "$pf_out_dir" # Because phyloFlash does not take output parameter
    # phyloFlash version display
    version=$(phyloFlash.pl -v 2>&1) # Redirect version output to stdout
    version="${version: -17}" # Show only name and version
    text_function "${version/v/}" # Remove "v" from version
    # phyloFlash task
    phyloFlash.pl \
    -read1 "$reads1" \
    -read2 "$reads2" \
    -lib "$sample_fixed" \
    -dbhome "$pf_data" \
    -sc \
    -CPUs "$threads"
    builtin cd "$work_dir"
    conda deactivate
}

# fastp setup
fastp_out_dir="$work_dir/Processed_reads/"
fastp_out1="$fastp_out_dir/$(basename "$reads1" .fastq.gz)_processed.fastq.gz"
fastp_out2="$fastp_out_dir/$(basename "$reads2" .fastq.gz)_processed.fastq.gz"   

# fastp for preprocessing of raw reads
fastp_function() {
    conda activate assembly
    # fastp setup
    fastp_report_dir="$fastp_out_dir/Reports/"
    mkdir -p "$fastp_report_dir" # Because fastp does not create directories
    # fastp version display
    version=$(fastp -v 2>&1) # Redirect version output to stdout
    text_function "$version"
    # fastp task
    fastp \
    --in1 "$reads1" \
    --in2 "$reads2" \
    --out1 "$fastp_out1" \
    --out2 "$fastp_out2" \
    --html "$fastp_report_dir/$sample.html" \
    --json "$fastp_report_dir/$sample.json" \
    --correction \
    --detect_adapter_for_pe \
    --low_complexity_filter \
    --overrepresentation_analysis \
    --thread "$threads"
    conda deactivate
}

# FastQC for quality control of processed reads
fastqc_function() {
    conda activate assembly
    # FastQC setup
    fastqc_out_dir="$work_dir/Processed_reads/FastQC/"
    mkdir -p "$fastqc_out_dir" # Because FastQC does not create directories
    # FastQC version display
    version=$(fastqc -v 2> /dev/null) # Send stderr to null
    text_function "${version::6} ${version: -6}"  # Remove "v" from version
    # FastQC task
    fastqc \
    "$fastp_out1" \
    "$fastp_out2" \
    --outdir "$fastqc_out_dir" \
    --threads "$threads"

    # Copy archives to respective folder
    mkdir -p "$fastqc_out_dir/Archives/"
    mv $fastqc_out_dir/*.zip "$fastqc_out_dir/Archives/"
    conda deactivate
}

# SPAdes setup
spades_out_dir="$work_dir/Assemblies/SPAdes/$sample/"
scaffolds="$spades_out_dir/scaffolds.fasta"
scaffolds_out_dir="$work_dir/Assemblies/Scaffolds/$sample/"

# SPAdes for assembly of processed reads
spades_function() {
    conda activate assembly
    # SPAdes version display
    version=$(spades.py -v 2>&1)
    text_function "${version::6} ${version: -6}" # Show only name and remove "v"
    # SPAdes task
    spades.py \
    -1 "$fastp_out1" \
    -2 "$fastp_out2" \
    -o "$spades_out_dir" \
    --sc \
    --memory "$memory" \
    --threads "$threads"

    # Copy scaffolds.fasta to collective folder
    mkdir -p "$scaffolds_out_dir"
    cp "$scaffolds" "$scaffolds_out_dir/$sample.fasta"
    cp "$spades_out_dir/assembly_graph.fastg" "$scaffolds_out_dir/$sample.fastg"
    conda deactivate
}

# QUAST for general assembly statistics like N50 of raw genome
quast_raw_function() {
    conda activate analysis
    # QUAST setup
    quast_out_dir="$work_dir/Analyses/QUAST/raw/$sample/"
    # QUAST version display
    version=$(quast.py -v)
    text_function "${version::5} ${version:7:5}" # Remove "v" and build number
    # QUAST task
    quast.py \
    "$scaffolds" \
    --output "$quast_out_dir" \
    --label "$sample" \
    --split-scaffolds \
    --threads "$threads"
    conda deactivate
}

# BUSCO for raw genome completeness assessment
busco_raw_function() {
    conda activate analysis
    # BUSCO setup
    busco_lineages=("eukaryota_odb10" "alveolata_odb10" "bacteria_odb10")
    busco_out_dir="$work_dir/Analyses/BUSCO/raw/"
    # BUSCO version display
    version=$(busco -v)
    text_function "$version"
    # BUSCO task
    for busco_lineage in "${busco_lineages[@]}"; do
        busco \
        --in "$scaffolds" \
        --out "${sample}_${busco_lineage%_odb10}" \
        --out_path "$busco_out_dir" \
        --lineage_dataset "$busco_lineage" \
        --mode genome \
        --offline \
        --download_path "$busco_data" \
        --cpu "$threads" \
        --force # Overwrites data
    done
    conda deactivate
}

# Whokaryote setup
whokaryote_out_dir="$work_dir/Assemblies/Whokaryote/$sample/"
eukaryotic_fasta="$whokaryote_out_dir/eukaryotes.fasta"
prokaryotic_fasta="$whokaryote_out_dir/prokaryotes.fasta"
contig_size="2000"

whokaryote_function() {
    conda activate whokaryote
    # Whokaryote setup
    mkdir -p "$whokaryote_out_dir"
    # Whokaryote version display
    text_function "Whokaryote 1.1.2" # Has no version parameter
    # Whokaryote task
    whokaryote.py \
    --contigs "$scaffolds" \
    --outdir "$whokaryote_out_dir" \
    --minsize "$contig_size" \
    --f \
    --threads "$threads"

    # Copy eukaryotic and prokaryotic contigs to collective folder
    cp "$eukaryotic_fasta" "$scaffolds_out_dir/$sample.euk$contig_size.fasta"
    cp "$prokaryotic_fasta" "$scaffolds_out_dir/$sample.prok$contig_size.fasta"
    conda deactivate
}

# QUAST for general assembly statistics like N50 of raw genome
quast_euk_function() {
    conda activate analysis
    # QUAST setup
    quast_out_dir="$work_dir/Analyses/QUAST/eukaryotic_$contig_size/$sample.euk$contig_size/"
    # QUAST version display
    version=$(quast.py -v)
    text_function "${version::5} ${version:7:5}" # Remove "v" and build number
    # QUAST task
    quast.py \
    "$eukaryotic_fasta" \
    --output "$quast_out_dir" \
    --pe1 "$fastp_out1" \
    --pe2 "$fastp_out2" \
    --label "$sample.euk$contig_size" \
    --min-contig "$contig_size" \
    --contig-thresholds "5000,10000,25000,50000" \
    --eukaryote \
    --rna-finding \
    --split-scaffolds \
    --threads "$threads"
    conda deactivate
}

# BUSCO for eukaryotic genome completeness assessment
busco_euk_function() {
    conda activate analysis
    # BUSCO setup
    busco_lineages=("eukaryota_odb10" "alveolata_odb10" "bacteria_odb10")
    busco_out_dir="$work_dir/Analyses/BUSCO/eukaryotic_$contig_size/"
    # BUSCO version display
    version=$(busco -v)
    text_function "$version"
    # BUSCO task
    for busco_lineage in "${busco_lineages[@]}"; do
        busco \
        --in "$eukaryotic_fasta" \
        --out "$sample.euk${contig_size}_${busco_lineage%_odb10}" \
        --out_path "$busco_out_dir" \
        --lineage_dataset "$busco_lineage" \
        --mode genome \
        --offline \
        --download_path "$busco_data" \
        --cpu "$threads" \
        --force # Overwrites data
    done
    conda deactivate
}

metaeuk_function() {
    conda activate analysis
    # MetaEuk setup
    metaeuk_out_dir="$work_dir/Predictions/$sample.euk$contig_size/"
    genes_out_dir="$work_dir/Predictions/Genes/"
    proteins_out_dir="$work_dir/Predictions/Proteins/"
    mkdir -p "$metaeuk_out_dir" "$genes_out_dir" "$proteins_out_dir"
    builtin cd "$metaeuk_out_dir"
    # MetaEuk version display
    text_function "MetaEuk 6" # Has no version parameter
    # MetaEuk task
    metaeuk easy-predict \
    "$eukaryotic_fasta" \
    "$protein_data" \
    "$sample.euk$contig_size" \
    "tmp"
    # Remove temporary folder
    rm -r "tmp"

    # Copy eukaryotic and prokaryotic contigs to collective folder
    cp "$sample.euk$contig_size.codon.fas" "$genes_out_dir/$sample.euk$contig_size.cds"
    cp "$sample.euk$contig_size.fas" "$proteins_out_dir/$sample.euk$contig_size.pep"
    conda deactivate
}

# Call functions and print stdout and stderr to both terminal and log file
phyloflash_function |& tee "$log_dir/$sample.phyloflash.log"
fastp_function |& tee "$log_dir/$sample.fastp.log"
fastqc_function |& tee "$log_dir/$sample.fastqc.log"
spades_function |& tee "$log_dir/$sample.spades.log"
quast_raw_function |& tee "$log_dir/$sample.quast.raw.log"
busco_raw_function |& tee "$log_dir/$sample.busco.raw.log"
whokaryote_function |& tee "$log_dir/$sample.whokaryote.log"
quast_euk_function |& tee "$log_dir/$sample.quast.euk$contig_size.log"
busco_euk_function |& tee "$log_dir/$sample.busco.euk$contig_size.log"
metaeuk_function |& tee "$log_dir/$sample.metaeuk.log"

# Delete empty folders and temporary files in case of error
builtin cd "$work_dir"
find . -type d -empty -delete
rm busco*
rm tmp*
