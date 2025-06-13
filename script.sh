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
export PATH=$PATH:"$PWD/rust_dbg/target/release"
path_fasta="./data/scerevisiae8.fa.gz"
path_split="./data/scerevisiae8_splits"
path_graphs="./data/scerevisiae8_graphs"
path_encoding="./data/scerevisiae8_encoding"

#____________________________________________________
# Split fasta by individuals
#____________________________________________________

# srun seqkit split $path_fasta -i --id-regexp "^([\w]+)#" -O $path_split

#____________________________________________________
# Create graphs containing different number of individuals
#____________________________________________________
# k=31

# fasta_names=("$path_split"/*.fa.gz)
# IFS=$'\n' sorted_fasta_names=($(sort <<<"${fasta_names[*]}"))
# unset IFS

# mkdir -p $path_graphs
# for i in $(seq 1 ${#sorted_fasta_names[@]}); do
#     srun ggcat build -k $k "${sorted_fasta_names[@]:0:i}" -o $path_graphs/n${i}_k${k}.fna --min-multiplicity 1 -j $SLURM_CPUS_PER_TASK
#     srun rust_dbg -k $k -g $path_graphs/n${i}_k${k}.bin build -i $path_graphs/n${i}_k${k}.fna
# done
# rm $path_graphs/*.stats.log

#____________________________________________________
# Get stats about graphs
#____________________________________________________

# k=31

# for n in {1..8}; do
#     srun rust_dbg -k $k -g $path_graphs/n${i}_k${k}.bin stats-g
# done

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