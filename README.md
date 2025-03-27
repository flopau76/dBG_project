# De-Bruijn Experiment

## Getting started

### Data preparation
Experimental dataset can be created using the scrip `data/get_data.sh`. This script uses [seqdd](https://github.com/yoann-dufresne/seqdd) to download data given their NCBI accession ids specified in `data/aedes.reg`.  
The samples are then renamed with seqkit to follow the [panSN-spec](https://github.com/pangenome/PanSN-spec) naming convention, before partitioning them by chromosomes.

### GGCAT
<!-- ```bash
seq 1 3 | while read i; do
    ggcat build -k 31 -s 1 -j 10 data/input/chr$i/AalbF3.fna data/input/chr$i/AalbF5.fna -o data/output/chr${i}/graph_k31.fna
done
``` -->

```bash
ggcat build -k 31 -s 1 -j 10 <input files> -o <output_file>
```
Used options are `-k` the size of the kmers, `-s` the minimum multiplicity to keep a kmer and `-j` the number of threads.  
Other possible options are 