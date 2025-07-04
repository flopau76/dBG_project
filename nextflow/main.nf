params.stranded = false
params.k = 23

params.input_files = "/home/florence/Documents/dbg_project/data/scerevisiae/*.fa.gz"

params.rust = "/home/florence/Documents/dbg_project/rust_dbg/target/release/rust_dbg"

process make_unitigs {
    conda 'bioconda::ggcat'
    
    input:
    path input_files
    
    output:
    path "k${params.k}.unitigs"

    script:
    """
    echo launching command: ggcat build -k ${params.k} $input_files --min-multiplicity 1 -o "k${params.k}.unitigs"
    ggcat build -k ${params.k} $input_files --min-multiplicity 1 -o "k${params.k}.unitigs"
    """
}

process build_graph {
    input:
    path unitigs

    output:
    path "k${params.k}.unitigs"

    script:
    """
    ${params.rust} -k ${params.k} build -i $unitigs -o "k${params.k}.unitigs"
    """
}

// Workflow block
workflow {
    input_files = Channel.fromPath(params.input_files)
    println "Input sequences: ${input_files}"
    make_unitigs(input_files)
    build_graph(make_unitigs.out)
}