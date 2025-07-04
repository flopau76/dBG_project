process build_unitigs {
    conda 'bioconda::ggcat'
    publishDir('data/graphs', pattern: "*.unitigs")
    input: path (fasta_file, arity: 1)
    output: tuple path(fasta_file), path ("${fasta_file}_k${params.k}.unitigs")

    script:
    """
    ggcat build -k ${params.k} $fasta_file --min-multiplicity 1 -o "${fasta_file}_k${params.k}.unitigs"
    """
}
process build_graph {
    publishDir('data/graphs', pattern: "*.bin")
    input: tuple path(fasta_file), path (unitig_file)
    output: tuple path(fasta_file), path ("${fasta_file}_k${params.k}.bin")

    script:
    """
    ${params.rust} -k ${params.k} build -i $unitig_file -o "${fasta_file}_k${params.k}.bin"
    """
}
process encode_paths {
    publishDir 'data/encodings'
    input: tuple path(fasta_file), path(bin_graph)
    output:
        path "${fasta_file}.encoding"
        path "${fasta_file}.encoding.txt"

    script:
    """
    echo "Encoding paths for ${fasta_file} using graph ${bin_graph}"
    ${params.rust} -k ${params.k} encode -i $fasta_file -g $bin_graph  -o "${fasta_file}.encoding" > "${fasta_file}.encoding.txt"
    """
}

// This should be replaced by an awk script
process split_communities {
    input: val input_files
    output: path "split_per_ch/*"

    script:
    """
    echo "Splitting communities from input files: ${input_files}"
    seqkit split ${input_files} -i --by-id-prefix part_ --id-regexp '^[a-zA-Z0-9_]+#1#([a-zA-Z0-9_]+)' -O split_per_ch
    """
}

workflow {
    split_communities(params.input_files)
    split_files = split_communities.output.flatten()
    build_unitigs(split_files)
    build_graph(build_unitigs.output)
    encode_paths(build_graph.output)
}