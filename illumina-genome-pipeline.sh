#!/bin/bash -l

# Author: Justin Teixeira Pereira Bassiaridis
# Date: 2024-05-06
# License: MIT

# This pipeline assembles paired-end Illumina short reads into a genome.
# It then decontaminates it, analyzes its quality, and predicts gene-coding proteins.
    # phyloFlash is used for contamination analysis of the raw reads.
    # fastp is used for preprocessing of the raw reads prior to assembly.
    # FastQC is used for quality control of the processed reads.
    # SPAdes is used to assemble the preprocessed reads into contigs/scaffolds.
    # QUAST is used to generate general assembly statistics like N50.
    # BUSCO is used for genome/proteome completeness assessment with several lineages.
    # Whokaryote is used for classification/decontamination of the scaffolds.
    # RepeatModeler and RepeatMasker are used to softmask the genome.
    # BRAKER is used to predict genes using a protein reference file.
# It runs offline, the data required for phyloFlash, BUSCO and BRAKER must be provided.

# Setup following conda environments before running:
    # Environment named "pf" containing phyloFlash 3.4.2.
    # Environment named "assembly" containing fastp 0.23.4, FastQC 0.12.1 and SPAdes 3.15.5.
    # Environment named "analysis" containing QUAST 5.2.0 and BUSCO 5.7.1.
    # Environment named "whokaryote" containing Whokaryote 1.1.2.
    # Environment named "softmasking" containing RepeatModeler 2.0.5 and RepeatMasker 4.1.5
    # May work with other versions of the listed software.
# Singularity/apptainer and a singularity/apptainer image of BRAKER 3.0.7 is needed.

###################################
# Basic setup, enter values here! #
###################################
# Working directory
work_dir=""

# Number of threads to use
threads=

# Amount of memory to use
memory=

# Whokaryote minimum contig size
contig_size=2000

# phyloFlash data directory
pf_data=""

# BRAKER 3.0.7 singularity/apptainer image file path and protein reference file path
braker_sif=""
ref_proteins=""

# BUSCO data directory and lineages to use
busco_data=""
busco_lineages=("eukaryota_odb10" "alveolata_odb10")  # Must be a list

# Sample name
sample=""

# Raw Illumina read file paths (must end with fastq.gz)
reads1=""
reads2=""
###################################
####### End of basic setup! #######
###################################

# Automatic setup
log_dir="$work_dir/Logs/$sample/"
out_dir="$work_dir/Output/$sample/"
assembly_out_dir="$out_dir/Assemblies/"
analysis_out_dir="$out_dir/Analyses/"
predict_out_dir="$out_dir/Predictions/"
mkdir -p "$log_dir" "$assembly_out_dir" "$analysis_out_dir" "$predict_out_dir"
cd "$work_dir"

# Text function for displaying software name and version
text_function() {
    # Short repeat function
    repeat() { for ((i=1; i<=$2; i++)); do printf "$1"; done; }
    printf "\n\n"
    printf "$(repeat = $((${#1} * 3)))"  # Repeat "=" 3x length of version output
    printf "\n\n"
    printf "$(repeat ' ' $((${#1} * 1)))$1"  # Repeat whitespace up to length of version output
    printf "\n\n"
    printf "$(repeat = $((${#1} * 3)))"
    printf "\n\n"
}

# phyloFlash for contamination analysis
phyloflash_function() {
    conda activate pf

    # phyloFlash setup
    sample_fixed="${sample//./_}"  # Because phyloFlash does not accept library names with dots
    pf_work_dir="$work_dir/Analyses/phyloFlash/$sample/"
    mkdir -p "$pf_work_dir"
    cd "$pf_work_dir"  # Because phyloFlash does not have an output parameter

    # phyloFlash version display
    version=$(phyloFlash.pl --version 2>&1)  # Redirect version output to stdout
    text_function "${version: -17}"  # Show only name and version

    # phyloFlash task
    phyloFlash.pl \
    -read1 "$reads1" \
    -read2 "$reads2" \
    -lib "$sample_fixed" \
    -dbhome "$pf_data" \
    -CPUs "$threads" \
    -sc  # For single-cell data
    cd "$work_dir"

    # Copy phyloFlash result file to Output folder
    cp "$pf_work_dir/$sample.phyloFlash.html" "$analysis_out_dir/$sample.reads.phyloflash.html"
    cp "$pf_work_dir/$sample.all.final.fasta" "$assembly_out_dir/$sample.reads.phyloflash_18S.fasta"

    conda deactivate
}

# fastp setup
fastp_work_dir="$work_dir/Processed_reads/"
fastp_reads1="$fastp_work_dir/$(basename "$reads1" .fastq.gz).processed.fastq.gz"
fastp_reads2="$fastp_work_dir/$(basename "$reads2" .fastq.gz).processed.fastq.gz"   

# fastp for preprocessing of raw reads
fastp_function() {
    conda activate assembly

    # fastp setup
    fastp_report_dir="$fastp_work_dir/Reports/"
    fastp_report="$fastp_report_dir/$sample.fastp.html"
    mkdir -p "$fastp_report_dir"  # Because fastp does not create directories

    # fastp version display
    version=$(fastp -v 2>&1)  # Redirect version output to stdout
    text_function "$version"

    # fastp task
    fastp \
    --in1 "$reads1" \
    --in2 "$reads2" \
    --out1 "$fastp_reads1" \
    --out2 "$fastp_reads2" \
    --html "$fastp_report" \
    --json "$fastp_report_dir/$sample.fastp.json" \
    --thread "$threads" \
    --correction \
    --detect_adapter_for_pe \
    --overrepresentation_analysis
    conda deactivate

    # Copy fastp result file to Output folder
    cp "$fastp_report" "$analysis_out_dir/$sample.reads.fastp.html"
}

# FastQC for quality control of processed reads
fastqc_function() {
    conda activate assembly

    # FastQC setup
    fastqc_work_dir="$fastp_work_dir/FastQC/$sample/"
    mkdir -p "$fastqc_work_dir"  # Because FastQC does not create directories

    # FastQC version display
    version=$(fastqc --version 2> /dev/null)  # Send stderr to null
    text_function "$version"

    # FastQC task
    fastqc \
    "$fastp_reads1" \
    "$fastp_reads2" \
    --outdir "$fastqc_work_dir" \
    --threads "$threads"

    # Copy archives to respective folder
    mkdir -p "$fastqc_work_dir/Archives/"
    mv "$fastqc_work_dir"/*.zip "$fastqc_work_dir/Archives/"

    # Copy FastQC result files to Output folder
    cp "$fastqc_work_dir"/*.html "$analysis_out_dir"

    conda deactivate
}

# SPAdes setup
assembly_work_dir="$work_dir/Assemblies/$sample/"
scaffolds="$assembly_work_dir/scaffolds.fasta"

# SPAdes for assembly of processed reads
spades_function() {
    conda activate assembly

    # SPAdes version display
    version=$(spades.py --version 2>&1)  # Redirect version output to stdout
    text_function "${version::6} ${version: -7}"  # Show only name and version

    # SPAdes task
    spades.py \
    -1 "$fastp_reads1" \
    -2 "$fastp_reads2" \
    -o "$assembly_work_dir" \
    --memory "$memory" \
    --threads "$threads" \
    --meta  # For metagenomes

    # Copy scaffolds and assembly graph to Output folder
    cp "$scaffolds" "$assembly_out_dir/$sample.assembly.fasta"
    cp "$assembly_work_dir/assembly_graph.fastg" "$assembly_out_dir/$sample.assembly_graph.fastg"
    cp "$assembly_work_dir/scaffolds.paths" "$assembly_out_dir/$sample.assembly.paths"

    conda deactivate
}

# Whokaryote setup
decont_work_dir="$work_dir/Decontamination/$sample/"
eukaryotic_scaffolds="$decont_work_dir/eukaryotes.fasta"
eukaryotic_sample="$sample.euk$contig_size"

# Whokaryote for decontamination of the genome
whokaryote_function() {
    conda activate whokaryote

    # Whokaryote setup
    mkdir -p "$decont_work_dir"

    # Whokaryote version display
    text_function "Whokaryote 1.1.2"  # Has no version parameter

    # Whokaryote task
    whokaryote.py \
    --contigs "$scaffolds" \
    --outdir "$decont_work_dir" \
    --minsize "$contig_size" \
    --threads "$threads" \
    --f  # Create filtered FASTA files

    # Copy eukaryotic and prokaryotic contigs to Output folder
    cp "$eukaryotic_scaffolds" "$assembly_out_dir/$eukaryotic_sample.fasta"
    cp "$decont_work_dir/prokaryotes.fasta" "$assembly_out_dir/$sample.prok$contig_size.fasta"
    cp "$decont_work_dir/unclassified.fasta" "$assembly_out_dir/$sample.unclassified$contig_size.fasta"

    conda deactivate
}

# Softmasking setup
softmask_work_dir="$work_dir/Softmasking/$sample/"
masked_scaffolds="$softmask_work_dir/eukaryotes.fasta.masked"

# RepeatModeler and RepeatMasker for softmasking decontaminated genome
softmasking_function() {
    conda activate softmasking

    # Softmasking setup
    export BLAST_USAGE_REPORT=false  # Do not send BLAST usage report over network
    mkdir -p "$softmask_work_dir"
    cd "$softmask_work_dir"

    # RepeatModeler version display
    version=$(RepeatModeler -version)
    text_function "$version"

    # Build database for RepeatModeler
    BuildDatabase \
    -name "$eukaryotic_sample" \
    "$eukaryotic_scaffolds"

    # Model repeats using database
    RepeatModeler \
    -threads="$threads" \
    -database "$eukaryotic_sample"

    # RepeatMasker version display
    version=$(RepeatMasker -v)
    text_function "$version"

    # Softmask genome based on model
    RepeatMasker \
    -engine "ncbi" \
    -pa "$(($threads / 4))" \
    -lib "$eukaryotic_sample-families.fa" \
    -dir "." \
    -xsmall \
    -gff \
    -alignments \
    -poly \
    -html \
    -source \
    "$eukaryotic_scaffolds"

    # Remove intermittent files
    rm -r RM_* *.euk2000.n*
    cd "$work_dir"

    # Copy masked genome, masking information and repeat families to Output folder
    cp "$masked_scaffolds" "$assembly_out_dir/$eukaryotic_sample.masked.fasta"
    cp "$softmask_work_dir/eukaryotes.fasta.out.html" "$predict_out_dir/$eukaryotic_sample.masked.html"
    cp "$softmask_work_dir/eukaryotes.fasta.out.gff" "$predict_out_dir/$eukaryotic_sample.masked.gff"
    cp "$softmask_work_dir/$eukaryotic_sample-families.fa" "$predict_out_dir/$eukaryotic_sample.families.fasta"

    conda deactivate
}

# BRAKER setup
predict_work_dir="$work_dir/Prediction/$sample/"
predicted_proteins="$predict_work_dir/braker.aa"

# BRAKER in protein-only mode to predict protein-coding genes
braker_function() {
    # BRAKER setup
    export APPTAINER_BIND="$work_dir"
    augustus_config="$work_dir/.augustus/config/"
    if [[ ! -d "$augustus_config" ]]; then
        mkdir --parents "$work_dir/.augustus/"
        apptainer exec "$braker_sif" cp --recursive /usr/share/augustus/config/ "$work_dir/.augustus/"
    elif [[ -d "$augustus_config/species/$eukaryotic_sample/" ]]; then
        rm --recursive "$augustus_config/species/$eukaryotic_sample/"
    fi
    mkdir --parents "$predict_work_dir"
    prot_seq="$predict_work_dir/prot_seq.fasta"
    cp "$ref_proteins" "$prot_seq"

    # BRAKER version display
    version=$(apptainer exec "$braker_sif" braker.pl -version)
    text_function "$version"

    # Run BRAKER via Singularity/Apptainer
    apptainer exec --cleanenv "$braker_sif" braker.pl \
    --species="$eukaryotic_sample" \
    --genome="$masked_scaffolds" \
    --prot_seq="$prot_seq" \
    --min_contig=10000 \
    --workingdir="$predict_work_dir" \
    --threads "$threads" \
    --gff3 \
    --AUGUSTUS_ab_initio \
    --AUGUSTUS_CONFIG_PATH="$augustus_config"

    # Copy predictions to Output folder and remove copy of protein reference file
    cp "$predicted_proteins" "$predict_out_dir/$eukaryotic_sample.braker.pep"
    cp "$predict_work_dir/braker.codingseq" "$predict_out_dir/$eukaryotic_sample.braker.cds"
    cp "$predict_work_dir/braker.gff3" "$predict_out_dir/$eukaryotic_sample.braker.gff"
    rm "$prot_seq"
}

# QUAST for general assembly statistics like N50 of decontaminated genome
quast_function() {
    conda activate analysis

    # QUAST version display
    version=$(quast.py --version)
    text_function "$version"

    # QUAST for raw assembly
        # QUAST setup
        quast_work_dir="$work_dir/Analyses/QUAST/raw/$sample/"

        # QUAST task
        quast.py \
        "$scaffolds" \
        --output "$quast_work_dir" \
        --label "$sample" \
        --min-contig "$contig_size" \
        --contig-thresholds "0,5000,10000,25000,50000,100000,250000,500000,1000000" \
        --threads "$threads" \
        --split-scaffolds

        # Copy QUAST result file to Output folder
        cp "$quast_work_dir/report.html" "$analysis_out_dir/$sample.assembly.quast.html"

    # QUAST for decontaminated assembly
        # QUAST setup
        quast_work_dir="$work_dir/Analyses/QUAST/euk$contig_size/$eukaryotic_sample/"

        # QUAST task
        quast.py \
        "$eukaryotic_scaffolds" \
        --output "$quast_work_dir" \
        --pe1 "$fastp_reads1" \
        --pe2 "$fastp_reads2" \
        --label "$eukaryotic_sample" \
        --min-contig "$contig_size" \
        --contig-thresholds "5000,10000,25000,50000,100000,250000,500000,1000000" \
        --threads "$threads" \
        --split-scaffolds

        # Copy QUAST result file to Output folder
        cp "$quast_work_dir/report.html" "$analysis_out_dir/$eukaryotic_sample.quast.html"

    conda deactivate
}

# BUSCO for completeness assessment
busco_function() {
    conda activate analysis

    # BUSCO version display
    version=$(busco --version)
    text_function "$version"

    # BUSCO loop for raw assembly
    for busco_lineage in "${busco_lineages[@]}"; do
        # BUSCO setup
        busco_work_dir="$work_dir/Analyses/BUSCO/raw/${busco_lineage%_odb10}/"

        # BUSCO task
        busco \
        --in "$scaffolds" \
        --out "$sample" \
        --out_path "$busco_work_dir" \
        --lineage_dataset "$busco_lineage" \
        --mode genome \
        --cpu "$threads" \
        --offline \
        --download_path "$busco_data"

        # Copy BUSCO result file to Output folder
        cp "$busco_work_dir/$sample/short_summary.specific.$busco_lineage.$sample.txt" \
        "$analysis_out_dir/$sample.assembly.busco_${busco_lineage%_odb10}.txt"
    done

    # BUSCO loop for decontaminated assembly
    for busco_lineage in "${busco_lineages[@]}"; do
        # BUSCO setup
        busco_work_dir="$work_dir/Analyses/BUSCO/euk$contig_size/${busco_lineage%_odb10}/"

        # BUSCO task
        busco \
        --in "$eukaryotic_scaffolds" \
        --out "$eukaryotic_sample" \
        --out_path "$busco_work_dir" \
        --lineage_dataset "$busco_lineage" \
        --mode genome \
        --cpu "$threads" \
        --offline \
        --download_path "$busco_data"
        
        # Copy BUSCO result file to Output folder
        cp "$busco_work_dir/$eukaryotic_sample/short_summary.specific.$busco_lineage.$eukaryotic_sample.txt" \
        "$analysis_out_dir/$eukaryotic_sample.busco_${busco_lineage%_odb10}.txt"
    done

    # BUSCO loop for predicted proteins
    for busco_lineage in "${busco_lineages[@]}"; do
        # BUSCO setup
        busco_work_dir="$work_dir/Analyses/BUSCO/euk${contig_size}_proteins/${busco_lineage%_odb10}/"
        busco_sample="$sample.euk${contig_size}_proteins"

        # BUSCO task
        busco \
        --in "$predicted_proteins" \
        --out "$busco_sample" \
        --out_path "$busco_work_dir" \
        --lineage_dataset "$busco_lineage" \
        --mode proteins \
        --cpu "$threads" \
        --offline \
        --download_path "$busco_data"
        
        # Copy BUSCO result file to Output folder
        cp "$busco_work_dir/$busco_sample/short_summary.specific.$busco_lineage.$busco_sample.txt" \
        "$analysis_out_dir/$busco_sample.busco_${busco_lineage%_odb10}.txt"
    done
    conda deactivate
}

# Call functions and print stdout and stderr to both terminal and log file
phyloflash_function |& tee "$log_dir/0_$sample.phyloflash.log"
fastp_function |& tee "$log_dir/1_$sample.fastp.log"
fastqc_function |& tee "$log_dir/2_$sample.fastqc.log"
spades_function |& tee "$log_dir/3_$sample.spades.log"
whokaryote_function |& tee "$log_dir/4_$sample.whokaryote.log"
softmasking_function |& tee "$log_dir/5_$sample.softmasking.log"
braker_function |& tee "$log_dir/6_$sample.braker.log"
quast_function |& tee "$log_dir/7_$sample.quast.log"
busco_function |& tee "$log_dir/8_$sample.busco.log"

# Copy logs to output folder
cp --recursive "$log_dir" "$out_dir/Logs/"
