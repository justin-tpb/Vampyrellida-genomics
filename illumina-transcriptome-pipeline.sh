#!/bin/bash -l

# Author: Justin Teixeira Pereira Bassiaridis
# Date: 2024-05-09
# License: MIT

# This pipeline assembles paired-end Illumina RNA reads into a transcriptome.
# It then analyzes its quality, and performs RNAseq quantification.
    # phyloFlash is used for contamination analysis of the raw reads.
    # RiboDetector is used to remove remaining rRNA reads.
    # Rcorrector is used to correct low quality bases based on k-mers.
    # fastp is used for preprocessing of the raw reads prior to assembly.
    # FastQC is used for quality control of the processed reads.
    # Trinity is used to assembly the preprocessed reads into transcripts.
    # Salmon is used to quantify the mRNA expression.
    # TransDecoder is used to translate transcripts into proteins.
    # QUAST is used to generate general assembly statistics like N50.
    # BUSCO is used for transcriptome completeness assessment with various lineages.
# It runs offline, the data required for phyloFlash and BUSCO must be provided.

# Setup following conda environments before running:
    # Environment named "pf" containing phyloFlash 3.4.2.
    # Environment named "transcriptomics" containing RiboDetector 0.3.1, Rcorrector 1.0.7
        # fastp 0.23.4, FastQC 0.12.1 and TransDecoder 5.7.1.
    # Environment named "analysis" containing QUAST 5.2.0 and BUSCO 5.7.1.
    # May work with other versions of the listed software.
# Singularity/apptainer and a singularity/apptainer image of Trinity 2.15.1 is needed.

###################################
# Basic setup, enter values here! #
###################################
# Working directory
work_dir=""

# Number of threads to use
threads=

# Amount of memory to use
memory=

# phyloFlash data directory
pf_data=""

# Trinity singularity/apptainer image file path and protein reference file path
trinity_img=""

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

# Work directory setup
log_dir="$work_dir/Logs/$sample/"
out_dir="$work_dir/Output/$sample/"
analysis_out_dir="$work_dir/Output/$sample/Analyses"
mkdir --parents "$log_dir" "$analysis_out_dir"
cd "$work_dir"
export APPTAINER_BIND="$work_dir"

# Read output directory setup
reads_work_dir="$work_dir/Processed_reads/"
mkdir --parents "$reads_work_dir"

# Text function for displaying software name and version at start
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
    pf_work_dir="$work_dir/Analysis/phyloFlash/$sample/"
    mkdir --parents "$pf_work_dir"
    cd "$pf_work_dir"  # Because phyloFlash does not take output parameter

    # phyloFlash version display
    version=$(phyloFlash.pl --version 2>&1)  # Redirect version output to stdout
    text_function "${version: -17}"  # Show only name and version

    # phyloFlash task
    phyloFlash.pl \
    -read1 "$reads1" \
    -read2 "$reads2" \
    -lib "$sample" \
    -dbhome "$pf_data" \
    -CPUs "$threads"
    cd "$work_dir"

    # Copy phyloFlash overview file and extracted SSU sequences to Output folder
    cp "$pf_work_dir/$sample.phyloFlash.html" "$analysis_out_dir/$sample.reads.phyloflash.html"
    cp "$pf_work_dir/$sample.all.final.fasta" "$analysis_out_dir/$sample.reads.phyloflash_18S.fasta"

    conda deactivate
}

# RiboDetector setup
ribodetector_out1="$reads_work_dir/$(basename "$reads1" .fastq.gz).non-rrna.fq"
ribodetector_out2="$reads_work_dir/$(basename "$reads2" .fastq.gz).non-rrna.fq"
ribodetector_rrna1="$reads_work_dir/$(basename "$reads1" .fastq.gz).rrna.fastq.gz"
ribodetector_rrna2="$reads_work_dir/$(basename "$reads2" .fastq.gz).rrna.fastq.gz"

# RiboDetector to remove rRNA reads
ribodetector_function() {
    conda activate transcriptomics

    # RiboDetector version display
    text_function "$(ribodetector --version)"

    # RiboDetector task
    ribodetector_cpu \
    --threads "$threads" \
    --len 101 \
    --input "$reads1" "$reads2" \
    --ensure rrna \
    --chunk_size 256 \
    --output "$ribodetector_out1" "$ribodetector_out2" \
    --rrna "$ribodetector_rrna1" "$ribodetector_rrna2"

    conda deactivate
}

# Rcorrector to correct RNA reads based on k-mers
rcorrector_function() {
    conda activate transcriptomics

    # Rcorrector version display
    text_function "Rcorrector 1.0.7"  # Has no version parameter

    # Rcorrector task
    run_rcorrector.pl \
    -1 "$ribodetector_out1" \
    -2 "$ribodetector_out2" \
    -od "$reads_work_dir" \
    -t "$threads"

    # Remove uncorrected reads
    rm "$ribodetector_out1"
    rm "$ribodetector_out2"

    conda deactivate
}

# fastp setup
fastp_in1="$reads_work_dir/$(basename "$reads1" .fastq.gz).non-rrna.cor.fq"
fastp_in2="$reads_work_dir/$(basename "$reads2" .fastq.gz).non-rrna.cor.fq"  
fastp_out1="$reads_work_dir/$(basename "$reads1" .fastq.gz).non-rrna.cor.processed.fastq.gz" 
fastp_out2="$reads_work_dir/$(basename "$reads2" .fastq.gz).non-rrna.cor.processed.fastq.gz" 

# fastp for preprocessing of raw reads
fastp_function() {
    conda activate transcriptomics

    # fastp setup
    fastp_report_dir="$reads_work_dir/fastp/"
    fastp_report="$fastp_report_dir/$sample.fastp.html"
    mkdir --parents "$fastp_report_dir"  # Because fastp does not create directories

    # fastp version display
    version=$(fastp -v 2>&1)  # Redirect version output to stdout
    text_function "$version"

    # fastp task
    fastp \
    --in1 "$fastp_in1" \
    --in2 "$fastp_in2" \
    --out1 "$fastp_out1" \
    --out2 "$fastp_out2" \
    --html "$fastp_report" \
    --json "$fastp_report_dir/$sample.fastp.json" \
    --thread "$threads" \
    --correction \
    --detect_adapter_for_pe \
    --length_required 36


    # Remove corrected but untrimmed reads
    rm "$fastp_in1"
    rm "$fastp_in2"

    # Copy fastp result file to Output folder
    cp "$fastp_report" "$analysis_out_dir/$sample.reads.fastp.html"

    conda deactivate
}

# FastQC for quality control of processed reads
fastqc_function() {
    conda activate transcriptomics

    # FastQC setup
    fastqc_work_dir="$reads_work_dir/FastQC/"
    mkdir --parents "$fastqc_work_dir"  # Because FastQC does not create directories

    # FastQC version display
    version=$(fastqc -v 2> /dev/null)  # Send stderr to null
    text_function "$version"

    # FastQC task
    fastqc \
    "$fastp_out1" \
    "$fastp_out2" \
    --outdir "$fastqc_work_dir" \
    --threads "$threads"

    # Copy FastQC result files to Output folder
    cp "$fastqc_work_dir"/$(basename "$fastp_out1" .fastq.gz)_fastqc.html "$analysis_out_dir/$sample.reads1.fastqc.html"
    cp "$fastqc_work_dir"/$(basename "$fastp_out2" .fastq.gz)_fastqc.html "$analysis_out_dir/$sample.reads2.fastqc.html"

    conda deactivate
}

# Trinity setup
trinity_work_dir="$work_dir/Assembly/Trinity/$sample-trinity/"
transcripts="${trinity_work_dir%/}.Trinity.fasta"

# Trinity for assembly of processed reads
trinity_function() {
    # Trinity version display
    version=$(apptainer exec "$trinity_img" Trinity --version | grep "Trinity version")
    text_function "${version: -15}"

    # Trinity task
    apptainer exec \
    --cleanenv \
    "$trinity_img" \
    Trinity \
    --seqType fq \
    --SS_lib_type RF \
    --left "$fastp_out1" \
    --right "$fastp_out2"  \
    --output "$trinity_work_dir" \
    --CPU "$threads" \
    --max_memory "${memory}G"

    # Copy transcripts and gene/transript mappings to Output folder
    cp "$transcripts" "$out_dir/$sample.trinity.fasta"
    cp "${trinity_work_dir%/}.Trinity.fasta.gene_trans_map" "$out_dir/$sample.trinity.gene_trans_map"

    # Generate assembly statistics
    apptainer exec \
    --cleanenv \
    "$trinity_img" \
    /usr/local/bin/util/TrinityStats.pl "$transcripts" \
    > "$analysis_out_dir/$sample.trinity_stats.txt"

    # Remove Trinity work directory
    rm --recursive "$trinity_work_dir"
}

salmon_function() {
    # Salmon setup
    salmon_work_dir="$work_dir/Quantification/$sample/"
    mkdir --parents "$salmon_work_dir"
    cd "$salmon_work_dir"

    # Salmon version display
    text_function "$(apptainer exec "$trinity_img" salmon --version)"

    # Salmon task
    apptainer exec \
    --cleanenv \
    "$trinity_img" \
    /usr/local/bin/util/align_and_estimate_abundance.pl \
    --seqType fq \
    --SS_lib_type RF \
    --est_method salmon \
    --transcripts "$transcripts" \
    --left "$fastp_out1" \
    --right "$fastp_out2" \
    --output_dir "$salmon_work_dir" \
    --thread_count "$threads" \
    --trinity_mode \
    --prep_reference

    # Copy Salmon quantifications on transcript and gene level to Output folder
    cp "$salmon_work_dir/quant.sf" "$out_dir/$sample.quant.sf"
    cp "$salmon_work_dir/quant.sf.genes" "$out_dir/$sample.quant.sf.genes"
}

transdecoder_function(){
    conda activate transcriptomics

    # TransDecoder setup
    transdecoder_work_dir="$work_dir/Assembly/TransDecoder/$sample/"
    mkdir --parents "$transdecoder_work_dir"
    cd "$transdecoder_work_dir"

    # TransDecoder version display
    text_function "TransDecoder 5.7.1"  # Has no version parameter

    # TransDecoder task
    TransDecoder.LongOrfs -t "$transcripts"
    TransDecoder.Predict -t "$transcripts"

    # Copy TransDecoder results to Output folder
    cp "$transdecoder_work_dir/$sample-trinity.Trinity.fasta.transdecoder.bed" \
    "$out_dir/$sample.transdecoder.bed"
    cp "$transdecoder_work_dir/$sample-trinity.Trinity.fasta.transdecoder.cds" \
    "$out_dir/$sample.transdecoder.cds"
    cp "$transdecoder_work_dir/$sample-trinity.Trinity.fasta.transdecoder.gff3" \
    "$out_dir/$sample.transdecoder.gff3"
    cp "$transdecoder_work_dir/$sample-trinity.Trinity.fasta.transdecoder.pep" \
    "$out_dir/$sample.transdecoder.pep"

    conda deactivate
}

# QUAST for general assembly statistics
quast_function() {
    conda activate analysis

    # QUAST setup
    quast_work_dir="$work_dir/Analysis/QUAST/$sample/"

    # QUAST version display
    version=$(quast.py --version)
    text_function "$version"

    # QUAST task
    quast.py "$transcripts" \
    --output "$quast_work_dir" \
    --pe1 "$fastp_out1" \
    --pe2 "$fastp_out2" \
    --label "$sample" \
    --threads "$threads"

    conda deactivate

    # Copy QUAST report to Output folder
    cp "$quast_work_dir/report.html" "$analysis_out_dir/$sample.quast.html"
}

# BUSCO for genome completeness assessment
busco_function() {
    conda activate analysis

    # BUSCO version display
    version=$(busco --version)
    text_function "$version"

    # BUSCO loop
    for busco_lineage in "${busco_lineages[@]}"; do
        # BUSCO setup
        busco_work_dir="$work_dir/Analysis/BUSCO/${busco_lineage%_odb10}"

        # BUSCO task
        busco \
        --in "$transcripts" \
        --out "$sample" \
        --out_path "$busco_work_dir" \
        --lineage_dataset "$busco_lineage" \
        --mode transcriptome \
        --cpu "$threads" \
        --offline \
        --download_path "$busco_data"

        # Copy BUSCO result file to Output folder
        cp "$busco_work_dir/$sample/short_summary.specific.$busco_lineage.$sample.txt" \
        "$analysis_out_dir/$sample.busco_${busco_lineage%_odb10}.html"
    done
    conda deactivate
}

# Call functions and print stdout and stderr to both terminal and log file
phyloflash_function |& tee "$log_dir/0_$sample.phyloflash.log"
ribodetector_function |& tee "$log_dir/1_$sample.ribodetector.log"
rcorrector_function |& tee "$log_dir/2_$sample.rcorrector.log"
fastp_function |& tee "$log_dir/3_$sample.fastp.log"
fastqc_function |& tee "$log_dir/4_$sample.fastqc.log"
trinity_function |& tee "$log_dir/5_$sample.trinity.log"
transdecoder_function |& tee "$log_dir/6_$sample.transdecoder.log"
salmon_function |& tee "$log_dir/7_$sample.salmon.log"
quast_function |& tee "$log_dir/8_$sample.quast.log"
busco_function |& tee "$log_dir/9_$sample.busco.log"

# Copy logs to output folder
cp --recursive "$log_dir" "$out_dir/Logs/"
