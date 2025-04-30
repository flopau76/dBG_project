//! Minor stats about debruijn graphs and their paths

use std::io::Write;
use std::time::Instant;

use crate::graph::Graph;
use crate::fasta_reader::FastaReader;
use crate::path::NodeIterator;
use debruijn::{Dir, Kmer};

use ahash::AHashSet;

use rand::Rng;

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
        kmer_count += seq.len() - K::k() + 1;
        let mut node_iter = NodeIterator::new(graph, &seq).unwrap();
        while let Some((node_id, dir)) = node_iter.next().unwrap() {
            let node = graph.get_node(node_id);
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
    eprintln!("Haplo contains:\n  - nodes: {}\n  - kmers: {}\n  - breakpoints: {}  (>1)\t\t{}  (<1)", node_count, kmer_count, junctions, dead_ends);
}

/// Get the number of nodes in the frontier set, depending on the depth in the BFS
pub fn stats_depth_bfs<K: Kmer>(graph: &Graph<K>, max_depth: usize, nb_repeats: usize) {
    eprintln!("Performing random BFS... ");

    let mut total_depth_histo = vec![0; max_depth];

    for _ in 0..nb_repeats {
        // generate a random start node
        let mut rng = rand::rng();
        let start_id = rng.random_range(0..graph.len());
        let start_dir = match rng.random_bool(0.5) {
            true => Dir::Left,
            false => Dir::Right,
        };

        // initialize the BFS
        let mut visited = AHashSet::new();
        let mut frontier_set = vec![(start_id, start_dir)];
        let mut depth = 0;
        let mut depth_histo = vec![0; max_depth];
        while depth < max_depth {
            let mut new_frontier_set = vec![];
            for (node_id, node_dir) in frontier_set {
                // get the edges of the current unitig
                for (neigh_id, neigh_dir, _) in graph.get_node(node_id).edges(node_dir.flip()) {
                    let new = visited.insert((neigh_id, neigh_dir));
                    if new {
                        new_frontier_set.push((neigh_id, neigh_dir));
                    }
                }
            }
            depth_histo[depth] = new_frontier_set.len();
            frontier_set = new_frontier_set;
            depth += 1;
        }

        println!("{:?}", depth_histo);
        // add the depth histogram to the total
        for (depth, count) in depth_histo.iter().enumerate() {
            total_depth_histo[depth] += count;
        }
    }
    for depth in 0..max_depth {
        total_depth_histo[depth] /= nb_repeats;
    }
    // print the depth histogram
    eprintln!("Mean histogram: {:?}", total_depth_histo);
}