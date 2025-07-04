process make_unitigs {
    conda 'bioconda::ggcat'
    
    input:
    path fasta_files
    
    output:
    path "k${params.k}.unitigs"

    script:
    """
    ggcat build -k ${params.k} $fasta_files --min-multiplicity 1 -o "k${params.k}.unitigs"
    """
}

process build_graph {
    input:
    path unitigs

    output:
    path "k${params.k}.bin"

    script:
    """
    ${params.rust} -k ${params.k} build -i $unitigs -o "k${params.k}.bin"
    """
}

process encode_path {
    input:
    path bin_graph
    path fasta_file

    output:
    path "${fasta_file}.encoding"
    path "${fasta_file}.encoding.txt"

    script:
    """
    ${params.rust} -k ${params.k} encode -i $fasta_file -g $bin_graph  -o "${fasta_file}.encoding" > "${fasta_file}.encoding.txt"
    """
}

// Workflow block
workflow {
    input_files = Channel.fromPath(params.input_files)
    println "Input sequences: ${input_files}"
    make_unitigs(input_files)
    build_graph(make_unitigs.out)
    encode_path(build_graph.out, input_files)   // TODO: is this done separately for each file in input_files ?
}