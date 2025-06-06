#!/bin/bash
#
#SBATCH --job-name=test_graphs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --partition=seqbio
#SBATCH --qos=seqbio
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-8

module purge
module load SeqKit
module load ggcat

# define the parameter k
# k_values=(15 19 23 27 31)
# k=${k_values[$SLURM_ARRAY_TASK_ID]}
# k=$SLURM_ARRAY_TASK_ID

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
# Get the name of the newly created files
#____________________________________________________

# Initialize an empty array to store file names
file_names=()

# Iterate over the items in the folder
for item in "$path_split"/*; do
    # Check if the item is a file
    if [ -f "$item" ]; then
        # Append the file name to the array
        file_names+=("$(basename "$item")")
    fi
done

# Sort the file names
IFS=$'\n' sorted_file_names=($(sort <<<"${file_names[*]}"))
unset IFS

#____________________________________________________
# Create graphs containing different number of individuals
#____________________________________________________

# mkdir -p $path_graphs
# for n in $(seq 1 ${#sorted_file_names[@]}); do
#     temp_file_list=()
#     for j in $(seq 0 $((i - 1))); do
#         temp_file_list+=("$path_split/${sorted_file_names[j]}")
#     done
#     srun ggcat build -k $k ${temp_file_list[@]} -o $path_graphs/n${n}_k${k}.fna --min-multiplicity 1 -j 4
#     srun rust_dbg -k $k -g $path_graphs/n${n}_k${k}.bin build -i $path_graphs/n${n}_k${k}.fna
# done
# rm $path_graphs/*.stats.log

#____________________________________________________
# Encode paths
#____________________________________________________
n=8
k=13

fasta_name=${sorted_file_names[$((SLURM_ARRAY_TASK_ID-1))]}
fasta_id=$(echo "$fasta_name" | sed -E 's/.*\.part_([^.]+)\.fa\.gz/\1/')
path_input=$path_split/$fasta_name
path_output=$path_encoding/n${n}_k${k}_$fasta_id.bin
path_bin_graph=$path_graphs/n${n}_k${k}.bin

echo "Encoding path $path_input in graph: $path_bin_graph" >&2
srun rust_dbg -k $k -g $path_bin_graph encode -i $path_input -o $path_output