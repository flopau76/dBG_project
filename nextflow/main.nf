process build_graph {
    conda 'bioconda::ggcat'
    publishDir("${params.outdir}/graphs/bin", pattern: "*.bin", mode: 'copy')
    publishDir("${params.outdir}/graphs/unitigs", pattern: "*.unitigs", mode: 'move')
    input: path (fasta_file, arity: 1)
    output:
        path ("${fasta_file}_k${params.k}.unitigs"), emit: unitig_file
        tuple path(fasta_file), path ("${fasta_file}_k${params.k}.bin"), emit: bin_graph

    script:
    """
    ggcat build -k ${params.k} $fasta_file --min-multiplicity 1 -o "${fasta_file}_k${params.k}.unitigs"
    ${params.rust} build -k ${params.k} -i "${fasta_file}_k${params.k}.unitigs" -o "${fasta_file}_k${params.k}.bin"
    """
}
process encode_paths {
    publishDir ("${params.outdir}/encodings", mode: 'move')
    input: tuple path(fasta_file), path(bin_graph)
    output:
        path "${fasta_file}.encoding"
        path "${fasta_file}.encoding.txt"

    script:
    """
    ${params.rust} encode -i $fasta_file -g $bin_graph  -o "${fasta_file}.encoding" > "${fasta_file}.encoding.txt"
    """
}

process split_communities {
    publishDir "${params.outdir}"
    input: val input_files
    output: path "fasta/clusters/*.fa"

    script:
    """
    divide_by_chromosomes.sh \
        -f ${input_files} \
        -k ${params.k} \
        -o fasta
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

workflow {
    samples_list = params.input_files
    if (params.preprocess)
    out = split_communities(samples_list)
    else
    out = concat_fasta(samples_list)

    build_graph(out.flatten())
    encode_paths(build_graph.output.bin_graph)
}