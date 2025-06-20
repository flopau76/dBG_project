#!/bin/bash
#
#SBATCH --job-name=test_graphs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=seqbio
#SBATCH --qos=seqbio
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-8

module purge
module load SeqKit
module load ggcat


# Define the paths to the files and executable
path_data=".data/scerevisiae8"
export PATH=$PATH:"$PWD/rust_dbg/target/release"

path_fasta="$path_data/fasta.fa.gz"
path_split="$path_data/splits"
path_graphs="$path_data/graphs"
path_encoding="$path_data/encoding"

#____________________________________________________
# Encode paths
#____________________________________________________
n=8
k=23

mkdir -p $path_encoding

input_files=("$path_split"/*.fa.gz)
path_input="${input_files[$((SLURM_ARRAY_TASK_ID - 1))]}"
fasta_id=$(basename "$path_input" | sed -E 's/.*\.part_([^.]+)\.fa\.gz/\1/')
path_output="$path_encoding/n${n}_k${k}_${fasta_id}.bin"

echo "Encoding fasta: $file_name" >&2
srun rust_dbg -k $k -g $path_graphs/n${n}_k${k}.bin encode -i "$path_input" -o "$path_output"

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