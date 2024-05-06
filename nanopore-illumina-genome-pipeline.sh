#!/bin/bash -l

# Author: Justin Teixeira Pereira Bassiaridis
# Date: 2024-05-06
# License: MIT

# This pipeline assembles Nanopore long reads reads into a genome and polishes it with paired-end Illumina short reads.
# It then decontaminates it, analyzes its quality, and predicts gene-coding proteins.
    # Fastp is used for preprocessing of the raw Illumina reads before polishing.
    # Filtlong is used for preprocessing of the raw Nanopore reads before assembly.
    # Flye is used to assemble the preprocessed Nanopore reads into contigs/scaffolds.
    # NextPolish is used to polish the assembly with the preprocessed Illumina reads.
    # Pypolca is used to further polish the assembly with the preprocessed Illumina reads.
    # Whokaryote is used for classification/decontamination of the scaffolds.
    # RepeatModeler and RepeatMasker are used to softmask the genome.
    # BRAKER is used to predict genes within the softmasked genome using a protein reference file.
    # QUAST is used to generate general assembly statistics like N50.
    # BUSCO is used for genome/proteome completeness assessment with several lineages.
# It runs offline, the data required for BUSCO and BRAKER must be provided.

# Setup following conda environments before running:
    # Environment named "nanopore" containing fastp 0.23.4, Filtlong v0.2.1, Flye 2.9.3,
        # bwa-mem2 2.2.1, nextPolish 1.4.1, SeqKit 2.8.0 and pypolca 0.3.1.
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

# Polishing iterations
flye_iterations=3
nextpolish_iterations=3
pypolca_iterations=3

# Whokaryote minimum contig size
contig_size=2000

# BRAKER 3.0.7 singularity/apptainer image file path and protein reference file path
braker_sif=""
ref_proteins=""

# BUSCO data directory and lineages to use
busco_data=""
busco_lineages=("eukaryota_odb10" "alveolata_odb10")  # Must be a list

# Sample name
sample=""

# Raw Nanopore read file path (must end with fastq.gz)
nanopore_reads=""

# Raw Illumina read file paths (must end with fastq.gz)
illumina_reads1=""
illumina_reads2=""
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
builtin cd "$work_dir"

# Text function for displaying software name and version at start
text_function() {
    # Short repeat function from:
    # https://stackoverflow.com/questions/5349718/how-can-i-repeat-a-character-in-bash
    repeat() { printf "$1""%.s" $(eval "echo {1.."$(($2))"}"); }
    printf "\n\n"
    printf "$(repeat = $((${#1}*3)))"  # Repeat "=" 3x length of version output
    printf "\n\n"
    printf "$(repeat " " $((${#1}*1)))$1"  # Repeat " " up to length of version output
    printf "\n\n"
    printf "$(repeat = $((${#1}*3)))"
    printf "\n\n"
}

# fastp setup
read_work_dir="$work_dir/Processed_reads/"
fastp_reads1="$read_work_dir/$(basename "$illumina_reads1" .fastq.gz).processed.fastq.gz"
fastp_reads2="$read_work_dir/$(basename "$illumina_reads2" .fastq.gz).processed.fastq.gz"   

# fastp for preprocessing of raw Illumina short reads
fastp_function() {
    conda activate nanopore

    # fastp setup
    fastp_report_dir="$read_work_dir/Reports/"
    fastp_report="$fastp_report_dir/$sample.fastp.html"
    mkdir -p "$fastp_report_dir"  # Because fastp does not create directories

    # fastp version display
    version=$(fastp --version 2>&1)  # Redirect version output to stdout
    text_function "$version"

    # fastp task
    fastp \
    --in1 "$illumina_reads1" \
    --in2 "$illumina_reads2" \
    --out1 "$fastp_reads1" \
    --out2 "$fastp_reads2" \
    --html "$fastp_report" \
    --json "$fastp_report_dir/$sample.fastp.json" \
    --thread "$threads" \
    --cut_tail \
    --cut_mean_quality 25 \
    --average_qual 25 \
    --length_required 100 \
    --correction \
    --detect_adapter_for_pe \
    --overrepresentation_analysis
    conda deactivate

    # Copy fastp result file to Output folder
    cp "$fastp_report" "$analysis_out_dir/$sample.reads.fastp.html"
}

# Filtlong setup
filtlong_reads="$read_work_dir/$(basename "$nanopore_reads" .fastq.gz).processed.fastq.gz"

# Filtlong for preprocessing of the Nanopore long reads using the Illumina short reads as control
filtlong_function() {
    conda activate nanopore

    # Filtlong version display
    version=$(filtlong --version)
    text_function "$version"

    # Filtlong task
    filtlong \
    -1 "$fastp_reads1" \
    -2 "$fastp_reads2" \
    --keep_percent 90 \
    --min_length 1000 \
    --trim \
    --split 500 \
    "$nanopore_reads" |
    gzip > "$filtlong_reads"

    # SeqKit for basic information about the raw and processed Nanopore long reads
    seqkit stats \
    --all \
    --basename \
    --out-file "$analysis_out_dir/$sample.reads.seqkit.txt" \
    "$nanopore_reads" "$filtlong_reads"

    conda deactivate
}

# Flye setup
assembly_work_dir="$work_dir/Assemblies/$sample/"
scaffolds="$assembly_work_dir/assembly.fasta"
    
# Flye for assembly of Nanopore long reads
flye_function() {
    conda activate nanopore

    # Flye setup
    mkdir -p "$assembly_work_dir"

    # Flye version display
    version=$(flye --version)
    version="flye $version"  # Show name and version
    text_function "$version"

    # Flye task
    flye \
    --nano-hq "$nanopore_reads" \
    --read-error 0.03 \
    --out-dir "$assembly_work_dir" \
    --iterations "$flye_iterations" \
    --threads "$threads" \
    --scaffold \
    --meta  # For metagenomes

    # Copy scaffolds and assembly graphs/info to Output folder
    cp "$scaffolds" "$assembly_out_dir/$sample.assembly.fasta"
    cp "$assembly_work_dir/assembly_graph.gfa" "$assembly_out_dir/$sample.assembly_graph.gfa"
    cp "$assembly_work_dir/assembly_graph.gv" "$assembly_out_dir/$sample.assembly_graph.gv"
    cp "$assembly_work_dir/assembly_info.txt" "$assembly_out_dir/$sample.assembly_info.txt"
    cp "$assembly_work_dir/flye.log" "$assembly_out_dir/$sample.assembly.log"

    conda deactivate
}

# NextPolish setup
nextpolish_work_dir="$work_dir/Polishing/NextPolish/$sample/"
nextpolish_scaffolds="$nextpolish_work_dir/$sample.nextpolish.final.fasta"

# NextPolish for polishing the assembly with Illumina short reads
nextpolish_function() {
    conda activate nanopore

	# NextPolish setup
	mkdir -p "$nextpolish_work_dir"
    builtin cd "$nextpolish_work_dir"
    cp "$scaffolds" .
    input="assembly.fasta"

    # NextPolish version display
    version=$(nextPolish --version)
    text_function "$version"

    # NextPolish task, one loop for each round
    for ((i=1; i<=nextpolish_iterations; i++)); do
        # Step 1:
        # Index the genome file and do alignment
        bwa-mem2 index "$input"
        bwa-mem2 mem -t "$threads" "$input" "$fastp_reads1" "$fastp_reads2" |
        samtools view --threads 3 --excl-flags 0x4 --bam - |
        samtools fixmate -m --threads 3 - - |
        samtools sort -m 2G --threads 5 - |
        samtools markdup --threads 5 -r - sgs.sort.bam
        
        # Index bam and genome files
        samtools index -@ "$threads" sgs.sort.bam
        samtools faidx "$input"

        # Polish genome file
        python "$CONDA_PREFIX/share/nextpolish-1.4.1/lib/nextpolish1.py" \
        -g "$input" \
        -t 1 \
        -p "$threads" \
        -s sgs.sort.bam \
        > "$sample.polishtemp.fasta"
        input="$sample.polishtemp.fasta"

        # Step 2:
        # Index genome file and do alignment
        bwa-mem2 index "$input"
        bwa-mem2 mem -t "$threads" "$input" "$fastp_reads1" "$fastp_reads2" |
        samtools view --threads 3 --excl-flags 0x4 --bam - |
        samtools fixmate -m --threads 3 - - |
        samtools sort -m 2G --threads 5 - |
        samtools markdup --threads 5 -r - sgs.sort.bam
        
        # Index bam and genome files
        samtools index -@ "$threads" sgs.sort.bam
        samtools faidx "$input"

        # Polish genome file
        python "$CONDA_PREFIX/share/nextpolish-1.4.1/lib/nextpolish1.py" \
        -g "$input" \
        -t 2 \
        -p "$threads" \
        -s sgs.sort.bam \
        > "$sample.nextpolish.fasta"
        input="$sample.nextpolish.fasta"

        # Create iteration folder and copy corresponding assembly
	    mkdir -p "$nextpolish_work_dir/Iteration_$i/"
        cp "$sample.nextpolish.fasta" "$nextpolish_work_dir/Iteration_$i/"
    done
    # Remove intermittent files
    rm sgs.sort.* assembly.* *.polishtemp.* *.nextpolish.*

    # Simplify contig names, force uppercase letters and sort contigs
    sed "s/_np12.*//" "$nextpolish_work_dir/Iteration_$nextpolish_iterations/$sample.nextpolish.fasta" |
    seqkit seq \
    --upper-case |
    seqkit sort \
    --by-name \
    --natural-order \
    --out-file "$nextpolish_scaffolds"

    builtin cd "$work_dir"
    conda deactivate
}

# pypolca setup
pypolca_work_dir="$work_dir/Polishing/pypolca/$sample/"
polished_scaffolds="$pypolca_work_dir/Iteration_$pypolca_iterations/${sample}_corrected.fasta"

# pypolca for polishing the assembly with Illumina short reads
pypolca_function() {
    conda activate nanopore

    # pypolca setup
    input="$nextpolish_scaffolds"
    
    # pypolca version display
    version=$(pypolca --version)
    text_function "$version"

    # pypolca task, one loop for each round
    for ((i=1; i<=pypolca_iterations; i++)); do
        pypolca run \
        --assembly "$input" \
        --reads1 "$fastp_reads1" \
        --reads2 "$fastp_reads2" \
        --output "$pypolca_work_dir/Iteration_$i/" \
        --prefix "$sample" \
        --threads "$threads"
        input="$pypolca_work_dir/Iteration_$i/${sample}_corrected.fasta"
    done

    # Copy polished assembly to Output folder
    cp "$polished_scaffolds" "$assembly_out_dir/$sample.polished.fasta"

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
    cp "$decont_work_dir/featuretable_predictions_T.tsv" "$predict_out_dir/$sample.whokaryote.tsv"

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
    builtin cd "$softmask_work_dir"

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
    -engine ncbi \
    -pa "$(($threads / 4))" \
    -lib "$eukaryotic_sample-families.fa" \
    -dir . \
    -xsmall \
    -gff \
    -alignments \
    -poly \
    -html \
    -source \
    "$eukaryotic_scaffolds"

    # Remove intermittent files
    rm -r RM_* *.euk2000.n*
    builtin cd "$work_dir"

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
        "$polished_scaffolds" \
        --output "$quast_work_dir" \
        --label "$sample" \
        --min-contig "$contig_size" \
        --contig-thresholds "0,5000,10000,25000,50000,100000,250000,500000,1000000" \
        --threads "$threads" \
        --split-scaffolds

        # Copy QUAST result file to Output folder
        cp "$quast_work_dir/report.html" "$analysis_out_dir/$sample.polished.quast.html"

    # QUAST for decontaminated assembly
        # QUAST setup
        quast_work_dir="$work_dir/Analyses/QUAST/euk$contig_size/$eukaryotic_sample/"

        # QUAST task
        quast.py \
        "$eukaryotic_scaffolds" \
        --output "$quast_work_dir" \
        --nanopore "$nanopore_reads" \
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
        --in "$polished_scaffolds" \
        --out "$sample" \
        --out_path "$busco_work_dir" \
        --lineage_dataset "$busco_lineage" \
        --mode genome \
        --cpu "$threads" \
        --offline \
        --download_path "$busco_data"

        # Copy BUSCO result file to Output folder
        cp "$busco_work_dir/$sample/short_summary.specific.$busco_lineage.$busco_sample.txt" \
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
        cp "$busco_work_dir/$eukaryotic_sample/short_summary.specific.$busco_lineage.$busco_sample.txt" \
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
fastp_function |& tee "$log_dir/0_$sample.fastp.log"
filtlong_function |& tee "$log_dir/1_$sample.filtlong.log"
flye_function |& tee "$log_dir/2_$sample.flye.log"
nextpolish_function |& tee "$log_dir/3_$sample.nextpolish.log"
pypolca_function |& tee "$log_dir/4_$sample.pypolca.log"
whokaryote_function |& tee "$log_dir/5_$sample.whokaryote.log"
softmasking_function |& tee "$log_dir/6_$sample.softmasking.log"
braker_function |& tee "$log_dir/7_$sample.braker.log"
quast_function |& tee "$log_dir/8_$sample.quast.log"
busco_function |& tee "$log_dir/9_$sample.busco.log"

# Copy logs to output folder
cp --recursive "$log_dir" "$out_dir/Logs/"
