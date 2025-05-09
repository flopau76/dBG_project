#!/bin/bash

#_______________________________________________________________________________
# Download sequences with seqdd
#_______________________________________________________________________________

seqdd init -r aedes.reg
seqdd download -d input

# clean up
rm -rf logs
rm -rf .register

#_______________________________________________________________________________
# Rename sequences with seqkit, to follow PAN-Spec convention
#_______________________________________________________________________________

#AalbF3
seqkit replace -p '.*chr([\d\.]+).*' -r 'AalbF3#1#chr$1' input/GCA_018104305.1/GCA_018104305.1_AalbF3_genomic.fna.gz -o input/AalbF3.fna

#AalbF5
seqkit replace -p '.*chromosome ([\d\.]+).*' -r 'AalbF5#1#chr$1' input/GCF_035046485.1/GCF_035046485.1_AalbF5_genomic.fna.gz | \
seqkit replace -p '.*scaffold_([\d\.]+).*' -r 'AalbF5#1#scaf$1' | \
seqkit replace -p '.*mitochondrion.*' -r 'AalbF5#1#mit' \
    -o input/AalbF5.fna

# remove original files
rm -rf input/GCA_018104305.1
rm -rf input/GCF_035046485.1

#_______________________________________________________________________________
# Split the sequences by chromosome
#_______________________________________________________________________________

seq 1 3 | while read i; do
    mkdir -p input/chr$i
    seqkit grep input/AalbF3.fna chr$i -r -o input/chr$i/AalbF3.fna
    seqkit grep input/AalbF5.fna chr$i -r -o input/chr$i/AalbF5.fna
done

#_______________________________________________________________________________
# Split the sequences into chunks separated by NNNs
#_______________________________________________________________________________

seqkit seq input/chr1/AalbF5.fna -w 0 | \
    perl -pe 's/N+/\n>AalbF5#1#chr1\n/g' | \
    seqkit rename -1 | \
    seqkit seq -m 31  -o AalbF5_splitN.fna

#_______________________________________________________________________________
# Run GGCAT
#_______________________________________________________________________________
seq 1 1 | while read i; do
    ggcat build -k 31 -s 1 -j 10 input/chr$i/AalbF3.fna input/chr$i/AalbF5.fna -o output/chr${i}/graph_k31.fna
    ggcat build -k 31 -s 1 -j 10 input/chr$i/AalbF3.fna -o output/chr${i}/AalbF3_k31.fna
    ggcat build -k 31 -s 1 -j 10 input/chr$i/AalbF5.fna -o output/chr${i}/AalbF5_k31.fna
done

# args: -k: kmer size
#       -s: minimal kmer multiplicity for filtering
#       -j: number of threads

#_______________________________________________________________________________
# Draft: align reconstruction with input
#_______________________________________________________________________________
mkdir test_align
seq 1 10 | while read i; do
    # seqkit grep -p "AalbF5#1#chr1_$i" output/chr1/reconstruct_AalbF5_in_1.fna | seqkit replace -w 0 -p $ -r "#reconstruct" -o test_align/AalbF5#1#chr1_$i.fna
    # seqkit grep -p "AalbF5#1#chr1_$i" input/chr1/AalbF5_splitN.fna | seqkit replace -w 0 -p $ -r "#input" >> test_align/AalbF5#1#chr1_$i.fna
    clustalo -i test_align/AalbF5#1#chr1_$i.fna -o test_align/AalbF5#1#chr1_$i.aln --force
    seqkit -w 0 seq test_align/AalbF5#1#chr1_$i.aln -o test_align/AalbF5#1#chr1_$i.aln2
done