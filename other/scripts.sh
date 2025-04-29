INPUT_FOLDER=input/data/per_chromosome
OUTPUT_FOLDER=output

#_______________________________________________________________________________
# PGGB
#_______________________________________________________________________________

seq 1 3 | while read i; do
    echo "   Processing chromosome $i"
    pggb -n 2 -t 8 -i $INPUT_FOLDER/chr$i.fna.gz -o $OUTPUT_FOLDER/pggb/chr$i
done

#_______________________________________________________________________________
# Bifrost
#_______________________________________________________________________________
# option -c (for colors) overrides option -f (for output in fasta format)

seq 1 3 | while read i; do
    echo "   Processing chromosome $i"
    Bifrost build -k 31 -f -t 8 -r $INPUT_FOLDER/chr$i.fna.gz -o $OUTPUT_FOLDER/bifrost/chr$i
    Bifrost build -k 31 -c -t 8 -r $INPUT_FOLDER/chr$i.fna.gz -o $OUTPUT_FOLDER/bifrost/chr$i
done

#_______________________________________________________________________________
# GGCAT
#_______________________________________________________________________________
# opt-in feature when compiled -> add KC+km; recompilation necessary to activate/deactivate it
# option -e to add edges

seq 1 3 | while read i; do
    echo "   Processing chromosome $i"
    ggcat build -k 31 -s 1 -e -j 8 $INPUT_FOLDER/chr$i.fna.gz -o $OUTPUT_FOLDER/ggcat/chr$i.fna
done

#_______________________________________________________________________________
# Cactus (via docker image)
#_______________________________________________________________________________

docker run --rm --mount  type=bind,source=.,target=/data cactus \
    cactus-pangenome /temp /data/input/data/path_cactus.txt --outDir /data/output/cactus \
    --outName chr1 --reference aalbf5 --gfa full clip filter

# change ownership of output files
sudo chown -R $(whoami):$(whoami) output/cactus


#_______________________________________________________________________________
# ODGI: analyse graph properties
#_______________________________________________________________________________

GRAPH=$OUTPUT_FOLDER/pggb/chr1/chr1.fna.gz.6771046.7bdde5a.6abc6ee.smooth.fix.og
OUTPUT=$OUTPUT_FOLDER/odgi/pggb_chr1.w5kbps
odgi depth -i $GRAPH -r GCF#1#chr1 | \
    bedtools makewindows -b /dev/stdin -w 5000 > $OUTPUT.bed    # create bed file with 5k bp interval windows along a path

odgi depth -i $GRAPH -b $OUTPUT.bed --threads 8 | \
    bedtools sort > ${OUTPUT}_depth.bed                         # compute depth for each window

#_______________________________________________________________________________
# SSHASH: create dict to query dBG
#_______________________________________________________________________________

GRAPH=$OUTPUT_FOLDER/ggcat/chr1.fna
tools/sshash/build/sshash build --verbose -i $GRAPH -k 31 -m 13 -o ${GRAPH}.sshash_index   #  --weighted

tools/sshash/build/sshash query -i $GRAPH.sshash_index -q query_test.fasta    # --multiline


#_______________________________________________________________________________
# map AalbF3 (scaffold genome) to AalbF5 (assembled chromosomes)
#_______________________________________________________________________________

# input
QUERY=database/aedes/data/GCA_018104305.1/GCA_018104305.1_AalbF3_genomic.fna.gz
TARGET=database/aedes/data/GCF_035046485.1/GCF_035046485.1_AalbF5_genomic.fna.gz
MERGED=database/aedes/data.fna.gz

# output
PAF_FILE=outputs/wfmash/aedes.mapping.paf
OUTPUT_FASTA=database/aedes_community

# 1: create pairwise alignment file using wfmash
wfmash $TARGET $QUERY -t 8 > $PAF_FILE

# 2: decompose paf into network (edges, weights and vertices)
python3 tools/pggb/scripts/paf2net.py -p $PAF_FILE

# 3: cluster network into communities
python3 tools/pggb/scripts/net2communities.py \
    -e $PAF_FILE.edges.list.txt \
    -w $PAF_FILE.edges.weights.txt \
    -n $PAF_FILE.vertices.id2name.txt \
    --plot

# 4: split the communities into separate fasta
seq 0 2 | while read i; do
    echo "community $i"
    samtools faidx $MERGED $(cat $PAF_FILE.edges.weights.txt.community.$i.txt) -o $OUTPUT_FASTA/community.$i.fna.gz
    samtools faidx $OUTPUT_FASTA/community.$i.fna.gz
done

#_______________________________________________________________________________
# run pggb separately on each community
#_______________________________________________________________________________

# OUTPUT_GRAPH=outputs/pggb/aedes

# seq 0 2 | while read i; do
#     echo "     community $i"
#     pggb -i $OUTPUT_FASTA/community.$i.fna.gz -o $OUTPUT_GRAPH/community.$i -n 2 -t 8
# done