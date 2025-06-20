#!/bin/bash
#
#SBATCH --job-name=make_graphs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --partition=seqbio
#SBATCH --qos=seqbio
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

module purge
module load SeqKit
module load ggcat


k=23


# Define the paths to the files and executable
path_data="./scerevisiae8"
export PATH=$PATH:"../rust_dbg/target/release"

path_fasta="$path_data/fasta.fa.gz"
path_split="$path_data/splits"
path_graphs="$path_data/graphs"

mkdir -p $path_graphs

#____________________________________________________
# Split fasta by individuals
#____________________________________________________

# srun seqkit split $path_fasta -i --id-regexp "^([\w]+)#" -O $path_split

#____________________________________________________
# Create graphs containing different number of individuals
#____________________________________________________

fasta_names=("$path_split"/*.fa.gz)
IFS=$'\n' sorted_fasta_names=($(sort <<<"${fasta_names[*]}"))
unset IFS

for i in $(seq 1 ${#sorted_fasta_names[@]}); do
    srun ggcat build -k $k "${sorted_fasta_names[@]:0:i}" -o $path_graphs/n${i}_k${k}.fna --min-multiplicity 1 -j $SLURM_CPUS_PER_TASK
    srun rust_dbg -k $k -g $path_graphs/n${i}_k${k}.bin build -i $path_graphs/n${i}_k${k}.fna
done
rm $path_graphs/*.stats.log

#____________________________________________________
# Get stats about graphs
#____________________________________________________

# for n in {1..8}; do
#     srun rust_dbg -k $k -g $path_graphs/n${i}_k${k}.bin stats-g
# done