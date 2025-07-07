process build_graph {
    conda 'bioconda::ggcat'
    publishDir('data/graphs/bin', pattern: "*.bin", mode: 'copy')
    publishDir('data/graphs/unitigs', pattern: "*.unitigs", mode: 'move')
    input: path (fasta_file, arity: 1)
    output:
        path ("${fasta_file}_k${params.k}.unitigs"), emit: unitig_file
        tuple path(fasta_file), path ("${fasta_file}_k${params.k}.bin"), emit: bin_graph

    script:
    """
    ggcat build -k ${params.k} $fasta_file --min-multiplicity 1 -o "${fasta_file}_k${params.k}.unitigs"
    ${params.rust} -k ${params.k} build -i "${fasta_file}_k${params.k}.unitigs" -o "${fasta_file}_k${params.k}.bin"
    """
}
process encode_paths {
    publishDir ('data/encodings', mode: 'move')
    input: tuple path(fasta_file), path(bin_graph)
    output:
        path "${fasta_file}.encoding"
        path "${fasta_file}.encoding.txt"

    script:
    """
    ${params.rust} -k ${params.k} encode -i $fasta_file -g $bin_graph  -o "${fasta_file}.encoding" > "${fasta_file}.encoding.txt"
    """
}

// Actual work of building the graph and encoding the paths
// Input is a single fasta file containing multiple sequences
workflow encoding {
    take:
    fasta_file

    main:
    build_graph(fasta_file)
    encode_paths(build_graph.output.bin_graph)
}

process split_communities {
    publishDir './'
    input: val input_files
    output: path "data/clusters/*.fa"

    script:
    """
    divide_by_chromosomes.sh \
        -f ${input_files} \
        -k ${params.k} \
        -o data
    """
}
process concat_fasta {
    input: val input_files
    output: path "concat.fa"

    script:
    """
    xargs cat < ${input_files} > concat.fa
    """
}

// Preprocessing step, which (if enabled) splits the input files into communities
workflow preprocess {
    take:
    samples_list

    main:
    if (params.preprocess)
    out = split_communities(samples_list)
    else
    out = concat_fasta(samples_list)

    emit:
    out
}

workflow {
    preprocess(params.input_files)
    encoding(preprocess.out.flatten())
}