#!/bin/bash

#_______________________________________________________________________________
# Download sequences with edirect from NCBI
#_______________________________________________________________________________

# cat accession_ids.txt | efetch -db nucleotide -format fasta > genomes.fasta

#_______________________________________________________________________________
# Download sequences with seqdd
#_______________________________________________________________________________

seqdd init -r aedes.reg
seqdd download -d input

# clean up
rm -rf logs
# rm -rf .register

#_______________________________________________________________________________
# Rename sequences with seqkit, to follow PAN-Spec convention
#_______________________________________________________________________________

echo "Renaming sequences with seqkit..."

#AalbF3
seqkit replace -p '.*chr([\d\.]+).*' -r 'AalbF3#1#chr$1' input/GCA_018104305.1/GCA_018104305.1_AalbF3_genomic.fna.gz -o input/AalbF3.fna.gz

#AalbF5
seqkit replace -p '.*chromosome ([\d\.]+).*' -r 'AalbF5#1#chr$1' input/GCF_035046485.1/GCF_035046485.1_AalbF5_genomic.fna.gz | \
seqkit replace -p '.*scaffold_([\d\.]+).*' -r 'AalbF5#1#scaf$1' | \
seqkit replace -p '.*mitochondrion.*' -r 'AalbF5#1#mit' \
    -o input/AalbF5.fna.gz

# remove original files
rm -rf input/GCA_018104305.1
rm -rf input/GCF_035046485.1

#_______________________________________________________________________________
# Split the sequences by chromosome
#_______________________________________________________________________________

echo "Splitting sequences by chromosome..."

seq 1 3 | while read i; do
    mkdir -p input/chr$i
    seqkit grep input/AalbF3.fna.gz -r -p chr$i -o input/chr$i/AalbF3.fna.gz
    seqkit grep input/AalbF5.fna.gz -r -p chr$i -o input/chr$i/AalbF5.fna.gz
done

#_______________________________________________________________________________
# Split the sequences by hardmasked regions
#_______________________________________________________________________________

# echo "Splitting sequences by hardmasked regions..."

# seq 1 3 | while read i; do
#     seqkit seq input/chr$i/AalbF5.fna.gz -w 0 | \
#         perl -pe "s/N+/\n>AalbF5#1#chr$i\n/g" | \
#         seqkit rename -1 | \
#         seqkit seq -m 31  -o input/chr$i/AalbF5_splitN.fna.gz

#     seqkit seq input/chr$i/AalbF3.fna.gz -w 0 | \
#         perl -pe "s/N+/\n>AalbF5#1#chr$i\n/g" | \
#         seqkit rename -1 | \
#         seqkit seq -m 31  -o input/chr$i/AalbF3_splitN.fna.gz
# done

#_______________________________________________________________________________
# Run GGCAT
#_______________________________________________________________________________
seq 1 3 | while read i; do
    mkdir -p ggcat/chr$i
    ggcat build -k 31 -s 1 -j 10 input/chr$i/AalbF3.fna.gz input/chr$i/AalbF5.fna.gz -o output/chr$i/graph_k31.fna
    ggcat build -k 31 -s 1 -j 10 input/chr$i/AalbF3.fna.gz -o output/chr$i/AalbF3_k31.fna
    ggcat build -k 31 -s 1 -j 10 input/chr$i/AalbF5.fna.gz -o output/chr$i/AalbF5_k31.fna
done

# args: -k: kmer size
#       -s: minimal kmer multiplicity for filtering
#       -j: number of threads

#_______________________________________________________________________________
# Draft: align reconstruction with input
#_______________________________________________________________________________
# mkdir test_align
# seq 1 6 | while read i; do
#     seqkit grep -p "AalbF5#1#chr1_$i" output/chr1/decoding_AalbF5_in_1.txt | seqkit replace -w 0 -p $ -r "#reconstruct" -o test_align/AalbF5#1#chr1_$i.fna
#     seqkit grep -p "AalbF5#1#chr1_$i" input/chr1/AalbF5_splitN.fna | seqkit replace -w 0 -p $ -r "#input" >> test_align/AalbF5#1#chr1_$i.fna
#     clustalo -i test_align/AalbF5#1#chr1_$i.fna -o test_align/AalbF5#1#chr1_$i.aln --force
#     seqkit -w 0 seq test_align/AalbF5#1#chr1_$i.aln -o test_align/AalbF5#1#chr1_$i.aln
# done