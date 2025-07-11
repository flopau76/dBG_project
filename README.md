# De-Bruijn Experiment

This project is about the reconstruction of haplotypes in de Bruijn graphs. It aims to enhance de Bruijn graphs with minimal supplementary information, to enable lossless reconstruction of the input haplotypes stored in the graph.

## Compiling the code
This project is coded in rust. To install Rust, please follow the instruction on this [page](https://www.rust-lang.org/tools/install).

Once this is done, the recommended way to install the tool is to clone the git repository and compile it using cargo, the official rust package manager:

```bash
git clone https://github.com/flopau76/dBG_project.git
cd dBG_project
cargo build --release --manifest-path rust_dbg/Cargo.toml
export PATH=$PATH:"$PWD/rust_dbg/target/release"
```

## Getting input files
Several pangenomes in FASTA format can be downloaded from [Zenodo](https://zenodo.org/records/7937947). Those datasets are compressed with bgzip and sequences follow the [panSN-spec](https://github.com/pangenome/PanSN-spec) naming convention.

## Generating a graph
The tool takes as an input a list of unitigs, as generated by [ggcat](https://github.com/algbio/ggcat). The unitigs can be created from (possibly compressed) fasta files with the following command:
```bash
ggcat build -k $K_SIZE $INPUT_FASTA -o $OUTPUT_UNITIGS --min-multiplicity 1
```
This set of unitigs is then converted into a binary graph using the following command:
```bash
rust_dbg build -k $K_SIZE -i $INPUT_PATH -o $OUTPUT_PATH 
```
By default, kmers are considered as canonical. If not, use the flag `--forward-only`.

> [!NOTE]
> It is necessary to use the same arguments `-k` and `--forward-only` (if used) with ggcat and with rust_dbg.

## En-/decoding a path
The core of the project consists in embedding continuous sequences as paths in the graph. For this, it takes as an input the previously generated graph and a set of sequences, in (potentially compressed) fasta format. Those sequences are embedded in the graph using the command `encode`:
```bash
rust_dbg encode -i $INPUT_FASTA -o $OUTPUT_ENCODING -g $PATH_GRAPH 
```
The encoding is composed of a list of the following extensions:
- TN(target node): shortest path to the given **T**arget **N**ode
- NN(next node): 2 bit encoding of the **N**ext **N**ucleotide to add
- R(pattern size)-(pattern offset): **R**epetition

The command `decode` performs the reverse operation, by retrieving a sequence from a path.
```bash
rust_dbg decode -i $INPUT_ENCODING -o $OUTPUT_FASTA -g $PATH_GRAPH
```
## Stats
Some basic stats on the graph can be obtained with the command `stats-g`.
```bash
rust_dbg stats-g -g $PATH_GRAPH
```
For an already encoded path, the command `stats-p` will print some information the encoding.
```bash
rust_dbg stats-p -i $INPUT_ENCODING -g $PATH_GRAPH
```