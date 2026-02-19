#!/bin/bash -l

# Author: Justin Teixeira Pereira Bassiaridis
# Date: 2025-12-22
# License: MIT

# This pipeline assembles Nanopore long reads into a genome and polishes it with paired-end Illumina short reads.
# It then decontaminates it, analyzes its quality, and predicts gene-coding proteins.
    # Fastp is used for preprocessing of the raw Illumina reads before polishing.
    # Filtlong is used for preprocessing of the raw Nanopore reads before assembly.
    # Flye is used to assemble the preprocessed Nanopore reads into contigs/scaffolds.
    # Polypolish is used to polish the assembly with the preprocessed Illumina reads.
    # Pypolca is used to further polish the assembly with the preprocessed Illumina reads.
    # Whokaryote is used for classification/decontamination of the contigs/scaffolds.
    # RepeatModeler and RepeatMasker are used to softmask the genome.
    # BRAKER is used to predict genes within the softmasked genome using a protein reference file.
    # QUAST is used to generate general assembly statistics like N50.
    # BUSCO is used for genome/proteome completeness assessment with several lineages.

# Use the parameter -h when running this script to display the usage information.
# It runs offline, the lineage data required for BUSCO must be manually downloaded first.

# Mamba is required with following environments:
    # Environment named "hybrid-assembly" containing fastp 1.0.1, Filtlong 0.3.1, Flye 2.9.6, bwa-mem2 2.3,
        # Polypolish 0.6.1, pypolca 0.4.0, RepeatModeler 2.0.7, RepeatMasker 4.2.2, QUAST 5.3.0 and BUSCO 6.0.0.
    # Environment named "whokaryote" containing Whokaryote 1.1.2.
    # Use the -m parameter to automatically install the environments.
    # May work with other versions of the listed software.
# Apptainer and an apptainer image of BRAKER 3.0.8 is needed.

set -eo pipefail  # Exit on error and fail on pipe errors
IFS=$'\n\t' # Split fields on newline and tab only
eval "$(mamba shell hook --shell bash)" # Initialize mamba in case of a non-interactive shell

# Default values
threads=12
contig_size=2000
busco_lineages=("eukaryota_odb12")
mamba_setup=false

# Help function to display usage
usage() {
    printf "Usage: %s [OPTIONS]\n" "$0"
    printf "Options:\n"
    printf "  -s, --sample TEXT            Sample name (required).\n"
    printf "  -d, --pipeline_dir PATH      Directory for the pipeline working and output data (required).\n"
    printf "  -n, --nanopore_reads PATH    Raw Nanopore read file path (fastq.gz) (required).\n"
    printf "  -1, --illumina_reads1 PATH   Raw Illumina R1 read file path (fastq.gz) (required).\n"
    printf "  -2, --illumina_reads2 PATH   Raw Illumina R2 read file path (fastq.gz) (required).\n"
    printf "  -i, --braker_sif PATH        BRAKER 3.0.8 apptainer image file path (required).\n"
    printf "  -r, --ref_proteins PATH      Protein reference file path for BRAKER (required).\n"
    printf "  -b, --busco_data PATH        BUSCO lineage data directory (required).\n"
    printf "  -l, --busco_lineages LIST    Comma-separated list of lineages (default: eukaryota_odb12).\n"
    printf "  -t, --threads INT            Number of threads to use (minimum: 4) (default: 12).\n"
    printf "  -c, --contig_size INT        Whokaryote minimum contig size (default: 2000).\n"
    printf "  -m, --mamba                  Automatically set up required mamba environments.\n"
    printf "  -h, --help                   Show this help message and exit.\n"
    exit 0
}

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -s|--sample) sample="$2"; shift 2 ;;
        -d|--pipeline_dir) pipeline_dir="$2"; shift 2 ;;
        -n|--nanopore_reads) nanopore_reads="$2"; shift 2 ;;
        -1|--illumina_reads1) illumina_reads1="$2"; shift 2 ;;
        -2|--illumina_reads2) illumina_reads2="$2"; shift 2 ;;
        -i|--braker_sif) braker_sif="$2"; shift 2 ;;
        -r|--ref_proteins) ref_proteins="$2"; shift 2 ;;
        -b|--busco_data) busco_data="$2"; shift 2 ;;
        -l|--busco_lineages) 
            # Convert comma-separated string into bash array
            IFS="," read -r -a busco_lineages <<< "$2"; shift 2 ;;
        -t|--threads) threads="$2"; shift 2 ;;
        -c|--contig_size) contig_size="$2"; shift 2 ;;
        -m|--mamba) mamba_setup=true; shift ;;
        -h|--help) usage ;;
        *) printf "Unknown option: %s\n" "$1"; usage ;;
    esac
done

# Validate threads parameter
if (( threads < 4 )); then
    printf "Error: --threads must be >= 4\n"
    exit 1
fi

# Function to handle automatic mamba environment setup
setup_mamba_envs() {
    printf "\n"
    printf "WARNING: Mamba environment setup requested!\n"
    printf "This will attempt to create the following environments:\n"
    printf "  - hybrid-assembly\n"
    printf "  - whokaryote\n"
    printf "\nIf these environments already exist, they will be overwritten.\n"
    printf "Do you want to proceed? [y/n] "
    read -r response

    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
        # Create mamba environments
        printf "\nStarting environment creation...\n"
        printf "\n>>> Creating 'hybrid-assembly' environment...\n"
        mamba create --name hybrid-assembly --channel bioconda --yes \
            fastp=1.0.1 \
            filtlong=0.3.1 \
            seqkit=2.12.0 \
            flye=2.9.6 \
            bwa-mem2=2.3 \
            polypolish=0.6.1 \
            pypolca=0.4.0 \
            repeatmodeler=2.0.7 \
            repeatmasker=4.2.2 \
            quast=5.3.0 \
            busco=6.0.0
        printf "\n>>> Creating 'whokaryote' environment...\n"
        mamba create --name whokaryote --channel bioconda --yes \
            whokaryote=1.1.2
        printf "\nAll environments created successfully.\n"
        printf "You can now run the pipeline without the -m flag.\n"
        exit 0
    else
        printf "Setup cancelled by user.\n"
        exit 0
    fi
}

# Run mamba setup if requested
if [ "$mamba_setup" = true ]; then
    setup_mamba_envs
fi

# Check if required arguments are provided
if [[ -z "$sample" || -z "$pipeline_dir" || -z "$nanopore_reads" || \
      -z "$illumina_reads1" || -z "$illumina_reads2" || -z "$braker_sif" || \
      -z "$ref_proteins" || -z "$busco_data" || ${#busco_lineages[@]} -eq 0 ]]; then
    printf "Error: Missing required arguments.\n"
    usage
fi

# Set directory variables
work_dir="$pipeline_dir/Work_dir/"
sample_dir="$work_dir/$sample/"
log_dir="$sample_dir/Logs/"
out_dir="$pipeline_dir/Output/$sample/"
assembly_out_dir="$out_dir/Assemblies/"
analysis_out_dir="$out_dir/Analyses/"
predict_out_dir="$out_dir/Predictions/"

# Print config, wait for user validation, create directories, and write log file
initialize_pipeline() {
    # Store configuration summary in a variable
    config_text=$(cat <<EOF

################################################################################
#######################  Pipeline Configuration Summary  #######################
################################################################################
Date:                 $(date)
Sample Name:          $sample
Pipeline Directory:   $pipeline_dir
Working Directory:    ${sample_dir//\/\///}
Output Directory:     ${out_dir//\/\///}
Threads:              $threads
Min Contig Size:      $contig_size bp

--- Read Data ---
Nanopore Reads:       $nanopore_reads
Illumina R1:          $illumina_reads1
Illumina R2:          $illumina_reads2

--- Other Data ---
BRAKER Image:         $braker_sif
Protein Reference:    $ref_proteins
BUSCO Data Dir:       ${busco_data//\/\///}
BUSCO Lineages:       $(IFS=", "; echo "${busco_lineages[*]}")
################################################################################

EOF
    )
    # Display configuration in terminal
    printf "%s\n\n" "$config_text"
    printf "Pipeline will start in 30 seconds! Press Ctrl+C to abort.\n\n"
    sleep 30

    # Create directories and change to sample directory
    mkdir -p "$log_dir" "$assembly_out_dir" "$analysis_out_dir" "$predict_out_dir"
    cd "$sample_dir"

    # Write configuration to log file
    printf "%s\n" "$config_text" > "$log_dir/00_$sample.config.log"
}

# Text function for displaying software name and version
print_version() {
    # Short repeat function
    repeat() { for ((i=1; i<=$2; i++)); do printf "%s" "$1"; done; }
    printf "\n\n"
    repeat "=" $((${#1} * 3))  # Repeat "=" 3x length of version output
    printf "\n\n"
    repeat " " $((${#1} * 1))  # Repeat whitespace up to length of version output
    printf "%s\n\n" "$1"
    repeat "=" $((${#1} * 3))
    printf "\n\n"
}

# fastp setup
read_sample_dir="$sample_dir/Processed_reads/"
fastp_reads1="$read_sample_dir/$(basename "$illumina_reads1" .fastq.gz).processed.fastq.gz"
fastp_reads2="$read_sample_dir/$(basename "$illumina_reads2" .fastq.gz).processed.fastq.gz"   

# fastp for preprocessing of raw Illumina short reads
run_fastp() {
    mamba activate hybrid-assembly

    # fastp setup
    fastp_report_dir="$read_sample_dir/Reports/"
    fastp_report="$fastp_report_dir/$sample.fastp.html"
    mkdir -p "$fastp_report_dir"  # Because fastp does not create directories

    # fastp version display
    version=$(fastp --version 2>&1)  # Redirect version output to stdout
    print_version "$version"

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
    mamba deactivate

    # Copy fastp result file to Output folder
    cp "$fastp_report" "$analysis_out_dir/$sample.reads.fastp.html"
}

# Filtlong setup
filtlong_reads="$read_sample_dir/$(basename "$nanopore_reads" .fastq.gz).processed.fastq.gz"

# Filtlong for preprocessing of the Nanopore long reads using the Illumina short reads as reference
run_filtlong() {
    mamba activate hybrid-assembly

    # Filtlong version display
    version=$(filtlong --version)
    print_version "$version"

    # Filtlong task
    filtlong \
    --short_1 "$fastp_reads1" \
    --short_2 "$fastp_reads2" \
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

    mamba deactivate
}

# Flye setup
assembly_sample_dir="$sample_dir/Assembly/"
assembly="$assembly_sample_dir/assembly.sorted.fasta"
    
# Flye for assembly of Nanopore long reads
run_flye() {
    mamba activate hybrid-assembly

    # Flye setup
    mkdir -p "$assembly_sample_dir"

    # Flye version display
    version=$(flye --version)
    version="flye $version"  # Show name and version
    print_version "$version"

    # Flye task
    flye \
    --nano-hq "$filtlong_reads" \
    --out-dir "$assembly_sample_dir" \
    --iterations 3 \
    --threads "$threads" \
    --scaffold \
    --meta  # For uneven coverage

    # Reorder assembly and copy output and assembly graphs/info/logs to Output folder
    seqkit sort \
    --natural-order \
    --out-file "$assembly" \
    "$assembly_sample_dir/assembly.fasta" 
    cp "$assembly" "$assembly_out_dir/$sample.assembly.fasta"
    cp "$assembly_sample_dir/assembly_graph.gfa" "$assembly_out_dir/$sample.assembly_graph.gfa"
    cp "$assembly_sample_dir/assembly_graph.gv" "$assembly_out_dir/$sample.assembly_graph.gv"
    cp "$assembly_sample_dir/assembly_info.txt" "$assembly_out_dir/$sample.assembly_info.txt"
    cp "$assembly_sample_dir/flye.log" "$assembly_out_dir/$sample.assembly.log"

    mamba deactivate
}

# Polypolish setup
polypolish_sample_dir="$sample_dir/Polishing/Polypolish/"
polypolish_assembly="$polypolish_sample_dir/$sample.polypolish.fasta"

# Polypolish for polishing the assembly with Illumina short reads
run_polypolish() {
    mamba activate hybrid-assembly

    # Polypolish setup
    mkdir -p "$polypolish_sample_dir"
    cd "$polypolish_sample_dir"
    cp "$assembly" "$polypolish_sample_dir"
    input=$(basename "$assembly")

    # Polypolish version display
    version=$(polypolish --version)
    print_version "$version"

    # Polypolish task
    bwa-mem2 index "$input"
    bwa-mem2 mem -t "$threads" -a "$input" "$fastp_reads1" > "$polypolish_sample_dir/alignments_1.sam"
    bwa-mem2 mem -t "$threads" -a "$input" "$fastp_reads2" > "$polypolish_sample_dir/alignments_2.sam"
    polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
    polypolish polish "$input" filtered_1.sam filtered_2.sam > "$polypolish_assembly"

    # Remove " polypolish" from polished headers and remove intermittent files
    sed --in-place "/^>/s/ polypolish$//" "$polypolish_assembly"
    rm --force "$input" *.amb *.ann *.0123 *.64 *.pac *.sam

    cd "$sample_dir"
    mamba deactivate
}

# pypolca setup
pypolca_sample_dir="$sample_dir/Polishing/pypolca/"
polished_assembly="$pypolca_sample_dir/${sample}_corrected.fasta"

# pypolca for polishing the assembly with Illumina short reads
run_pypolca() {
    mamba activate hybrid-assembly

    # pypolca setup
    input="$polypolish_assembly"
    
    # pypolca version display
    version=$(pypolca --version)
    print_version "$version"

    # pypolca task
    pypolca run \
    --assembly "$input" \
    --reads1 "$fastp_reads1" \
    --reads2 "$fastp_reads2" \
    --output "$pypolca_sample_dir" \
    --prefix "$sample" \
    --threads "$threads" \
    --careful  # Equivalent to --min_alt 4 --min_ratio 3

    # Copy polished assembly to Output folder
    cp "$polished_assembly" "$assembly_out_dir/$sample.polished.fasta"

    mamba deactivate
}

# Whokaryote setup
decont_sample_dir="$sample_dir/Decontamination/"
eukaryotic_assembly="$decont_sample_dir/eukaryotes.fasta"
eukaryotic_sample="$sample.euk$contig_size"

# Whokaryote for decontamination of the genome
run_whokaryote() {
    mamba activate whokaryote

    # Whokaryote setup
    mkdir -p "$decont_sample_dir"

    # Whokaryote version display
    print_version "Whokaryote 1.1.2"  # Has no version parameter

    # Whokaryote task
    whokaryote.py \
    --contigs "$polished_assembly" \
    --outdir "$decont_sample_dir" \
    --minsize "$contig_size" \
    --threads "$threads" \
    --f  # Create filtered FASTA files

    # Copy eukaryotic, prokaryotic and unclassified contigs to Output folder
    cp "$eukaryotic_assembly" "$assembly_out_dir/$eukaryotic_sample.fasta"
    cp "$decont_sample_dir/prokaryotes.fasta" "$assembly_out_dir/$sample.prok$contig_size.fasta"
    cp "$decont_sample_dir/unclassified.fasta" "$assembly_out_dir/$sample.unclassified$contig_size.fasta"
    cp "$decont_sample_dir/featuretable_predictions_T.tsv" "$predict_out_dir/$sample.whokaryote.tsv"

    mamba deactivate
}

# Softmasking setup
softmask_sample_dir="$sample_dir/Softmasking/"
masked_assembly="$softmask_sample_dir/eukaryotes.fasta.masked"

# RepeatModeler and RepeatMasker for softmasking decontaminated genome
run_softmasking() {
    mamba activate hybrid-assembly

    # Softmasking setup
    mkdir -p "$softmask_sample_dir"
    cd "$softmask_sample_dir"

    # RepeatModeler version display
    version=$(RepeatModeler -version)
    print_version "$version"

    # Build database for RepeatModeler
    BuildDatabase \
    -name "$eukaryotic_sample" \
    "$eukaryotic_assembly"

    # Model repeats using database
    RepeatModeler \
    -threads="$threads" \
    -database "$eukaryotic_sample"

    # RepeatMasker version display
    version=$(RepeatMasker -v)
    print_version "$version"

    # Softmask genome based on model
    RepeatMasker \
    -engine ncbi \
    -pa "$((threads / 4))" \
    -lib "$eukaryotic_sample-families.fa" \
    -dir . \
    -xsmall \
    -gff \
    -html \
    "$eukaryotic_assembly"

    # Remove intermittent files
    rm --force --recursive RM_* "$eukaryotic_sample".n*
    cd "$sample_dir"

    # Copy masked genome, masking information and repeat families to Output folder
    cp "$masked_assembly" "$assembly_out_dir/$eukaryotic_sample.masked.fasta"
    cp "$softmask_sample_dir/eukaryotes.fasta.out.html" "$predict_out_dir/$eukaryotic_sample.masked.html"
    cp "$softmask_sample_dir/eukaryotes.fasta.out.gff" "$predict_out_dir/$eukaryotic_sample.masked.gff"
    cp "$softmask_sample_dir/$eukaryotic_sample-families.fa" "$predict_out_dir/$eukaryotic_sample.families.fasta"

    mamba deactivate
}

# BRAKER setup
predict_sample_dir="$sample_dir/Prediction/"
predicted_proteins="$predict_sample_dir/braker.aa"

# BRAKER in protein-only mode to predict protein-coding genes
run_braker() {
    # BRAKER setup
    export APPTAINER_BIND="$work_dir"
    augustus_config="$work_dir/.augustus/config/"
    if [[ ! -d "$augustus_config" ]]; then
        mkdir --parents $(dirname "$augustus_config")
        apptainer exec "$braker_sif" cp --recursive /opt/Augustus/config/ "$augustus_config"
    elif [[ -d "$augustus_config/species/$eukaryotic_sample/" ]]; then
        rm --force --recursive "$augustus_config/species/$eukaryotic_sample/"
    fi
    mkdir --parents "$predict_sample_dir"
    prot_seq="$predict_sample_dir/prot_seq.fasta"
    cp "$ref_proteins" "$prot_seq"

    # BRAKER version display
    version=$(apptainer exec "$braker_sif" braker.pl -version)
    print_version "$version"

    # Run BRAKER via Apptainer
    apptainer exec --cleanenv "$braker_sif" braker.pl \
    --species="$eukaryotic_sample" \
    --genome="$masked_assembly" \
    --prot_seq="$prot_seq" \
    --min_contig=10000 \
    --workingdir="$predict_sample_dir" \
    --threads "$threads" \
    --gff3 \
    --AUGUSTUS_ab_initio \
    --AUGUSTUS_CONFIG_PATH="$augustus_config"

    # Copy predictions to Output folder and remove copy of protein reference file
    cp "$predicted_proteins" "$predict_out_dir/$eukaryotic_sample.braker.pep"
    cp "$predict_sample_dir/braker.codingseq" "$predict_out_dir/$eukaryotic_sample.braker.cds"
    cp "$predict_sample_dir/braker.gff3" "$predict_out_dir/$eukaryotic_sample.braker.gff3"
    cp "$predict_sample_dir/braker.gtf" "$predict_out_dir/$eukaryotic_sample.braker.gtf"
    rm --force "$prot_seq"
}

# QUAST for general assembly statistics like N50 of decontaminated genome
run_quast() {
    mamba activate hybrid-assembly

    # QUAST version display
    version=$(quast.py --version)
    print_version "$version"

    # QUAST setup
    quast_sample_dir="$sample_dir/Analysis/QUAST/$sample/"
    quast_eukaryotic_sample_dir="$sample_dir/Analysis/QUAST/$eukaryotic_sample/"

    # QUAST task for polished assembly
    quast.py \
    "$polished_assembly" \
    --output "$quast_sample_dir" \
    --label "$sample" \
    --min-contig "$contig_size" \
    --contig-thresholds "0,5000,10000,25000,50000,100000,250000,500000,1000000" \
    --threads "$threads" \
    --split-scaffolds

    # QUAST task for polished and decontaminated assembly
    quast.py \
    "$eukaryotic_assembly" \
    --output "$quast_eukaryotic_sample_dir" \
    --nanopore "$nanopore_reads" \
    --label "$eukaryotic_sample" \
    --min-contig "$contig_size" \
    --contig-thresholds "5000,10000,25000,50000,100000,250000,500000,1000000" \
    --threads "$threads" \
    --split-scaffolds

    # Copy QUAST reports to Output folder
    cp "$quast_sample_dir/report.html" "$analysis_out_dir/$sample.polished.quast.html"
    cp "$quast_eukaryotic_sample_dir/report.html" "$analysis_out_dir/$eukaryotic_sample.quast.html"

    mamba deactivate
}

# BUSCO for completeness assessment
run_busco() {
    mamba activate hybrid-assembly

    # BUSCO version display
    version=$(busco --version)
    print_version "$version"

    # Run BUSCO for each lineage
    for busco_lineage in "${busco_lineages[@]}"; do
        # BUSCO setup
        busco_sample_dir="$sample_dir/Analysis/BUSCO/$busco_lineage/"
        proteins_sample="${eukaryotic_sample}_proteins"

        # BUSCO task for polished assembly
        busco \
        --in "$polished_assembly" \
        --out "$sample" \
        --out_path "$busco_sample_dir" \
        --lineage_dataset "$busco_lineage" \
        --mode genome \
        --metaeuk \
        --cpu "$threads" \
        --offline \
        --download_path "$busco_data"

        # BUSCO task for polished and decontaminated assembly
        busco \
        --in "$eukaryotic_assembly" \
        --out "$eukaryotic_sample" \
        --out_path "$busco_sample_dir" \
        --lineage_dataset "$busco_lineage" \
        --mode genome \
        --metaeuk \
        --cpu "$threads" \
        --offline \
        --download_path "$busco_data"

        # BUSCO task for predicted proteins
        busco \
        --in "$predicted_proteins" \
        --out "$proteins_sample" \
        --out_path "$busco_sample_dir" \
        --lineage_dataset "$busco_lineage" \
        --mode proteins \
        --cpu "$threads" \
        --offline \
        --download_path "$busco_data"

        # Remove temporary BUSCO files 
        rm --force --recursive "$busco_sample_dir/$sample/tmp/"
        rm --force --recursive "$busco_sample_dir/$eukaryotic_sample/tmp/"
        rm --force --recursive "$busco_sample_dir/$proteins_sample/tmp/"

        # Copy BUSCO result files to Output folder
        cp "$busco_sample_dir/$sample/short_summary.specific.$busco_lineage.$sample.txt" \
        "$analysis_out_dir/$sample.polished.busco_$busco_lineage.txt"
        cp "$busco_sample_dir/$eukaryotic_sample/short_summary.specific.$busco_lineage.$eukaryotic_sample.txt" \
        "$analysis_out_dir/$eukaryotic_sample.busco_$busco_lineage.txt"
        cp "$busco_sample_dir/$proteins_sample/short_summary.specific.$busco_lineage.$proteins_sample.txt" \
        "$analysis_out_dir/$proteins_sample.busco_$busco_lineage.txt"
    done

    mamba deactivate
}

# Call functions and print stdout and stderr to both terminal and log file
# Comment out functions if you need to skip steps during a restart
initialize_pipeline
run_fastp |& tee "$log_dir/01_$sample.fastp.log"
run_filtlong |& tee "$log_dir/02_$sample.filtlong.log"
run_flye |& tee "$log_dir/03_$sample.flye.log"
run_polypolish |& tee "$log_dir/04_$sample.polypolish.log"
run_pypolca |& tee "$log_dir/05_$sample.pypolca.log"
run_whokaryote |& tee "$log_dir/06_$sample.whokaryote.log"
run_softmasking |& tee "$log_dir/07_$sample.softmasking.log"
run_braker |& tee "$log_dir/08_$sample.braker.log"
run_quast |& tee "$log_dir/09_$sample.quast.log"
run_busco |& tee "$log_dir/10_$sample.busco.log"

# Copy logs to output folder
cp --recursive "$log_dir" "$out_dir/Logs/"

# Final message
printf "\n\nPipeline finished. Output files are in: $out_dir\n"
