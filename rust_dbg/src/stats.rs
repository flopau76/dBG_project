//! Minor stats about debruijn graphs and their paths

use std::io::Write;
use std::time::Instant;

use crate::graph::Graph;
use crate::fasta_reader::FastaReader;
use crate::path::NodeIterator;
use debruijn::Kmer;

/// Print some stats about the graph
pub fn stats_graph<K: Kmer>(graph: &Graph<K>) {
    eprint!("Iterating graph... ");
    std::io::stderr().flush().unwrap();
    let start = Instant::now();
    let nb_nodes = graph.len();
    let mut nb_edges: usize = 0;
    let mut degree_histo = vec![0; 9];
    for exts in graph.base.exts.iter() {
        let degree = (exts.num_exts_l() + exts.num_exts_r()) as usize;
        nb_edges += degree;
        degree_histo[degree] += 1;
    };
    nb_edges /= 2;
    let duration = start.elapsed();
    eprintln!("done in {:?}", duration);
    eprintln!("Graph contains:\n  - {} nodes\n  - {} edges\n  degree histogramm: {:?}", nb_nodes, nb_edges, degree_histo);
}

/// Count the number of breakpoints in a graph, for a given haplotype
#[deprecated]
pub fn stats_haplo<K: Kmer>(graph: &Graph<K>, haplo: FastaReader) {
    eprint!("Iterating fasta file... ");
    std::io::stderr().flush().unwrap();
    let start = Instant::now();

    let mut kmer_count: usize = 0;
    let mut node_count: usize = 0;

    let mut junctions: usize = 0;
    let mut dead_ends: usize = 0;

    for record in haplo {
        let seq = record.dna_string();
        let mut node_iter = NodeIterator::new(graph, &seq);
        while let Some((node_id, dir)) = node_iter.next().unwrap() {
            let node = graph.get_node(node_id);
            let node_length = node.len();
            kmer_count += node_length;
            node_count += 1;

            // get the edges of the current unitig
            match node.exts().num_ext_dir(dir.flip()) {
                0 => { dead_ends += 1; },
                1 => { },
                _ => { junctions += 1; },
            }
        }
    }
    let duration = start.elapsed();
    eprintln!("done in {:?}", duration);
    eprintln!("Haplo contains:\n  - nodes {}\n  - kmers: {}\n  - breakpoints: {}  (>1)\t\t{}  (<1)", node_count, kmer_count, junctions, dead_ends);
}