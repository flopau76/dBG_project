#!/bin/bash
#
#SBATCH --job-name=test_graphs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --partition=seqbio
#SBATCH --qos=seqbio
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module purge
module load SeqKit
module load ggcat


# Define the paths to the files and executable
export PATH=$PATH:"$PWD/rust_dbg/target/release"
path_fasta="./data/scerevisiae8.fa.gz"
path_split="./data/scerevisiae8_splits"
path_graphs="./data/scerevisiae8_graphs"
k=31

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

mkdir -p $path_graphs
for i in $(seq 1 ${#sorted_file_names[@]}); do
    temp_file_list=()
    for j in $(seq 0 $((i - 1))); do
        temp_file_list+=("$path_split/${sorted_file_names[j]}")
    done
    srun ggcat build -k $k ${temp_file_list[@]} -o $path_graphs/n${i}_k${k}.fna --min-multiplicity 1 -j 4
    srun rust_dbg -k $k -g $path_graphs/n${i}_k${k}.bin build -i $path_graphs/n${i}_k${k}.fna
done
rm $path_graphs/*.stats.log