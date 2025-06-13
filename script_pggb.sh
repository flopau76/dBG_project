#!/bin/bash
#
#SBATCH --job-name=build_pggb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=5G
#SBATCH --partition=seqbio
#SBATCH --qos=seqbio
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
# #SBATCH --array=1-8

input_fasta="./data/scerevisiae8.fa.gz"
output_dir="./data/pggb_output"
mkdir -p $output_dir

module load wfmash seqwish smoothxg GFAffix vg samtools odgi MultiQC pggb

pggb -i $input_fasta -o $output_dir -t 16