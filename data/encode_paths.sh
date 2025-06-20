#!/bin/bash
#
#SBATCH --job-name=encode_paths
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=seqbio
#SBATCH --qos=seqbio
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
# #SBATCH --array=1-8

module purge
module load SeqKit
module load ggcat


n=8
k=23


# Define the paths to the files and executable
path_data="./scerevisiae8"
export PATH=$PATH:"../rust_dbg/target/release"

path_fasta="$path_data/fasta.fa.gz"
path_split="$path_data/splits"
path_graphs="$path_data/graphs"
path_encoding="$path_data/encoding"

mkdir -p $path_encoding

#____________________________________________________
# Encode paths
#____________________________________________________

input_files=("$path_split"/*.fa.gz)
for path_input in "${input_files[@]}"; do
    fasta_id=$(basename "$path_input" | sed -E 's/.*\.part_([^.]+)\.fa\.gz/\1/')
    path_output="$path_encoding/n${n}_k${k}_${fasta_id}.bin"

    echo "Encoding fasta: $fasta_id" >&2
    srun rust_dbg -k $k -g $path_graphs/n${n}_k${k}.bin encode -i "$path_input" -o "$path_output"
done

#____________________________________________________
# Get stats about paths
#____________________________________________________
# n=8
# k=23

# for file in $path_encoding/*; do
#     if [[ -f $file ]]; then
#         srun rust_dbg -k $k -g $path_graphs/n${n}_k${k}.bin stats-p -i $file
#     fi
# done