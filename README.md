# De-Bruijn Experiment

## Getting started
This project is about the reconstruction of haplotypes in de Bruijn graphs. It aims to enhance de Bruijn graphs with minimal supplementary information, to enable lossles reconstruction of the input haplotypes stored in the graph.

### Dataset generation
An experimental dataset can be created using the script `data/get_data.sh`. Required dependencies are:
- [seqdd](https://github.com/yoann-dufresne/seqdd)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [ggcat](https://github.com/algbio/ggcat)

This script first uses seqdd to fetch data from the NCBI database, given the accession IDs in `data/aedes.reg`. The samples are then renamed with seqkit to follow the [panSN-spec](https://github.com/pangenome/PanSN-spec) naming convention, and partitionned by chromosomes. Finally, a de Bruijn Graph is created using GGCAT.

### Rust package
The core of the project is the rust package contained in `rust_dbg`. It is separated into three modules:
- `fasta_reader` parses fasta files and enumerates their kmer.
- `graph` implements the Graph structure, which allows to query the presence/absence of a kmer and to obtain its edges, if present.
- `path` performs breadth-first search in the graph.

A running example can be seen in `main.rs`. It takes as an input a set of unitigs as generated by ggcat and transforms them into a graph. The graph can be dumped into and reconstructed quickly from a binary file. Various functions are available to perform statistics or look for haplotypes in the graph.

For more details about the package, see the documentation generated by cargo doc. 
