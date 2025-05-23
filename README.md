# De-Bruijn Experiment

## Getting started
This project is about the reconstruction of haplotypes in de Bruijn graphs. It aims to enhance de Bruijn graphs with minimal supplementary information, to enable lossles reconstruction of the input haplotypes stored in the graph.

### Dataset generation
An experimental dataset can be created using the script `data/get_data.sh`. Required dependencies are:
- [seqdd](https://github.com/yoann-dufresne/seqdd)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [ggcat](https://github.com/algbio/ggcat)

This script first uses seqdd to fetch data from the NCBI database, given the accession IDs in `data/aedes.reg`. The samples are then renamed with seqkit to follow the [panSN-spec](https://github.com/pangenome/PanSN-spec) naming convention, and partitionned by chromosomes. Finally, a de Bruijn Graph is created using GGCAT.

## Rust package
### Generating a graph
This tool uses a binary graph format, wich must first be generated from a list of unitigs, with the following command:
```bash
cargo run -- -g $PATH_GRAPH -k $K_SIZE build $INPUT_PATH
```
By default, kmers are considered as canonical. If not, use the flag `--forward-only`. Howevere, note that the arguments used here must be the same as those used to generate the input file.

Some basic stats on the graph can be obtained with the command `stats-g`:
```bash
cargo run -- -g $PATH_GRAPH -k $K_SIZE stats-g
```
### En-/decoding a path
The core of the project consists in embedding continuous sequences as paths in the graph. Those paths are described in a custom text format. To generate one, use the command `encode`:
```bash
cargo run -- -g $PATH_GRAPH -k $K_SIZE encode -i $INPUT_FASTA -o $OUTPUT_TEXT
```
The encoding is composed of a start node, and a list of extensions of the following format:
- SP(target node): shortest path to the given node
- NN(next node): 2 bit encoding of the following node to add
- R(pattern size)x(nb of repetitions): coming soon

> [!NOTE]
> The encoding is not lossless. For now, when a shortest path is specified, it is garantied that the input path is one of the shortest paths to the target node, but it may not be the only one.

The command `decode` performs the reverse operation, by retrieving a sequence from a path:
```bash
cargo run -- -g $PATH_GRAPH -k $K_SIZE decode -i $INPUT_TEXT  -o $OUTPUT_FASTA
```
### Stats
For an allready encoded path, the command `stats-p` will print some information about the weight and the efficiency of the encoding.
```bash
cargo run -- -g $PATH_GRAPH -k $K_SIZE stats-p
```