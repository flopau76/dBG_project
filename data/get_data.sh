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
    seqkit faidx input/AalbF3.fna chr$i -r -o input/chr$i/AalbF3.fna
    seqkit faidx input/AalbF5.fna chr$i -r -o input/chr$i/AalbF5.fna
done

#_______________________________________________________________________________
# Split the sequences to remove Ns
#_______________________________________________________________________________

seqkit seq AalbF5.fna -w 0 | \
    perl -pe 's/N+/\n>AalbF5#1#chr1\n/g' | \
    seqkit rename -1 | \
    seqkit seq -m 31  -o AalbF5_splitN.fna
