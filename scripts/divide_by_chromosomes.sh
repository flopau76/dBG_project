### THIS SCRIPTS EXPECTS TO RECEIVE AS INPUT A TXT FILE WITH A LIST OF GENOMES IN FASTA FILES
### IT EXPECTS THE FIRST GENOME TO BE A REFERENCE, I.E. BEING ASSEMBLED AT CHROMOSOME LEVEL
### SO ALL THE SEQUENCES OF THE OTHER GENOMES WILL BE COMPARED TO IT TO DECIDE IN WHICH CLUSTER
### THEY WILL BE PLACED IN

### IT OUTPUTS A LIST OF FASTAS CONTAINING THE READS CLUSTERED BY REFERENCE CHROMOSOME

#!/bin/bash

set -e 

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# cd "$SCRIPT_DIR"


CORES_NUMBER=$(getconf _NPROCESSORS_ONLN)
# Default values
DEFAULT_INPUT_FILES_LIST="./input_fastas.txt"
DEFAULT_OUTPUT_DIRECTORY="../data"
DEFAULT_KMER_SIZE=31
## MODE 1: SAVE THE FASTAS (FOR REFERENCE)
## MODE 2: DELTE THE FASTAS AFTER APPENDING (FOR THE REST)
DEFAULT_MODE=1

# Parse command-line options
while getopts "f:k:o:r:" opt; do
  case $opt in
    f)
      INPUT_FILE="$OPTARG"
      ;;
    k)
      KMER_SIZE="$OPTARG"
      ;;
    o)
      OUTPUT_DIRECTORY="$OPTARG"
      ;;
    r)
      MODE="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

INPUT_FILE="${INPUT_FILE:-$DEFAULT_INPUT_FILES_LIST}"
KMER_SIZE="${KMER_SIZE:-$DEFAULT_KMER_SIZE}"
OUTPUT_DIRECTORY="${OUTPUT_DIRECTORY:-$DEFAULT_OUTPUT_DIRECTORY}"
MODE="${MODE:-$DEFAULT_MODE}"

### CHECKING FILES EXIST
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: INPUT FILE (-f) '$INPUT_FILE' not found."
    exit 1
fi

### TEMPORARY FILES AND FOLDERS
TEMP_DIR=$(mktemp -d -t "divide_XXXXX")
# trap "rm -rf $TEMP_DIR" EXIT
echo "Generated temp dir $TEMP_DIR"

# TEMP FILES FOR REFERENCE
REF_GENOME_FILE=$TEMP_DIR/ref_file.txt

CLUSTERS_DIR=$TEMP_DIR/clusters
CLUSTERS_LIST=$TEMP_DIR/clusters.txt
CLUSTER_SKETCHES_DIR=$TEMP_DIR/clusters_sketches
CLUSTER_SKETCHES=$TEMP_DIR/clusters_sketch.txt

# TEMP FILES FOR OTHER GENOMES
CURRENT_GENOME_DIR=$TEMP_DIR/temp_genomes
CURRENT_GENOME_SKETCH_DIR=$TEMP_DIR/temp_sketches
CURRENT_GENOME_LIST=$TEMP_DIR/temp_genomes.txt
CURRENT_GENOME_SKETCH_LIST=$TEMP_DIR/temp_genomes.txt

DISTANCE_MATRIX=$TEMP_DIR/distance.tsv 

### 1. DIVIDE THE GENOME BY ITS FASTA SEQUENCES (WILL USE THE FASTA TO APPEND THE OTHER SEQUENCES IN THE SAME COMMUNITY)
mkdir -p $CLUSTERS_DIR
REF_GENOME_FILE=$(head -n1 $INPUT_FILE)
echo "REF GENOME: $REF_GENOME_FILE"
./scripts/divide_by_sequences.sh $CLUSTERS_DIR > $CLUSTERS_LIST < "$REF_GENOME_FILE"
echo "CLUSTERS LIST : $CLUSTERS_LIST"
### 1b. RUN DASHING SKETCHING TO SKETCH THE SEQUENCE
# Sketching the chromosomes and creating a combined archived
mkdir -p $CLUSTER_SKETCHES_DIR
dashing sketch -k$KMER_SIZE -p$CORES_NUMBER --use-bb-minhash -F $CLUSTERS_LIST -P $CLUSTER_SKETCHES_DIR/
find $CLUSTER_SKETCHES_DIR/ -maxdepth 1 -type f > $CLUSTER_SKETCHES

### 2 FOR EACH FILE, DIVIDE IT BY SEQUENCES, ESTIMATE DISTANCE FOR THE SEQUENCES, ASSIGN BY CHROMOSOME
while IFS= read -r CURRENT_FASTA; do
    echo "Processing file: '$CURRENT_FASTA'"

    if [[ -f "$CURRENT_FASTA" ]]; then
        mkdir $CURRENT_GENOME_DIR $CURRENT_GENOME_SKETCH_DIR
    
        # dividing sequences into single files
        ./scripts/divide_by_sequences.sh $CURRENT_GENOME_DIR > $CURRENT_GENOME_LIST < $CURRENT_FASTA
        
        # sketching sequences
        echo "RUNNING dashing sketch -k$KMER_SIZE -p$CORES_NUMBER --use-bb-minhash -F $CURRENT_GENOME_LIST -P $CURRENT_GENOME_SKETCH_DIR/"
        dashing sketch -k$KMER_SIZE -p$CORES_NUMBER --use-bb-minhash -F $CURRENT_GENOME_LIST -P $CURRENT_GENOME_SKETCH_DIR/
        find $CURRENT_GENOME_SKETCH_DIR/ -maxdepth 1 -type f > $CURRENT_GENOME_SKETCH_LIST
        
        #calculating distance
        echo "RUNNING dashing dist --presketched -p$CORES_NUMBER -F $CLUSTER_SKETCHES -Q $CURRENT_GENOME_SKETCH_LIST -O $DISTANCE_MATRIX --full-tsv --containment-index"
        dashing dist --presketched -k$KMER_SIZE -p$CORES_NUMBER --use-bb-minhash -F $CLUSTER_SKETCHES -Q $CURRENT_GENOME_SKETCH_LIST -O $DISTANCE_MATRIX --containment-index

        # moving each file to its fasta
        while IFS=$'\t' read -r sketch_file values; do
            # Extract base filename
            filename=$(basename "$sketch_file" | sed 's/\.w\..*//')
            tmpdir=$(dirname "$(dirname "$sketch_file")")
            
            # Find max value and position
            max_val=0
            max_pos=1
            pos=1
            for val in $values; do
                if (( $(echo "$val > $max_val" | bc -l) )); then
                    max_val=$val
                    max_pos=$pos
                fi
                ((pos++))
            done
            
            # Get target file from clusters.txt
            target_sketch=$(sed -n "${max_pos}p" "$CLUSTER_SKETCHES")
            target_fasta=$(basename "$target_sketch" | sed 's/\.w\..*//')
            # echo ">>>> MY FILENAME: $target_sketch"
            # echo ">>>> MY PATH: $CURRENT_GENOME_DIR/$target_fasta"
            # Append genome to target
            cat "$CURRENT_GENOME_DIR/$filename" >> "$CLUSTERS_DIR/$target_fasta"
        done < $DISTANCE_MATRIX

        # ereasing the folder with the fastas
        rm -r $CURRENT_GENOME_DIR $CURRENT_GENOME_SKETCH_DIR
    else
        echo "Warning: File '$CURRENT_FASTA' not found. Skipping."
    fi
done < <(tail -n +2 "$INPUT_FILE")


mkdir -p $OUTPUT_DIRECTORY
mv $CLUSTERS_DIR $OUTPUT_DIRECTORY