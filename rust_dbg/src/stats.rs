//! Minor stats about debruijn graphs and their paths

use std::io::Write;
use std::time::Instant;

use crate::Graph;
use debruijn::Kmer;

/// Print some stats about the graph
pub fn stats_graph<K: Kmer>(graph: &Graph<K>) {
    print!("Iterating graph... ");
    std::io::stdout().flush().unwrap();
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
    println!("done in {:?}", duration);
    println!("Graph contains:\n  - {} nodes\n  - {} edges\n  degree histogramm: {:?}", nb_nodes, nb_edges, degree_histo);
}

// /// Count the number of breakpoints in a graph, for a given haplotype
// #[deprecated]
// pub fn stats_haplo<K: Kmer>(graph: &Graph<K>, haplo: FastaReader) {
//     print!("Iterating haplo... ");
//     std::io::stdout().flush().unwrap();
//     let start = Instant::now();

//     let mut kmer_count: usize = 0;

//     let mut junctions: usize = 0;
//     let mut dead_ends: usize = 0;

//     for record in haplo {
//         let mut kmer_iter = record.iter_kmers(false);
//         let mut current_kmer = kmer_iter.next();
//         while let Some(kmer) = current_kmer {
//             // check that the kmer is in the graph at the beginning of a unitig
//             let lookup = graph.find_link(kmer, Dir::Right);
//             if lookup.is_none() {
//                 println!("Kmer {} not found at the beginning of a unitig\n", kmer_count);
//                 continue;
//             }
//             let (node_id, dir, _) = lookup.unwrap();
//             let node = graph.get_node(node_id);
//             let node_length = node.len();
//             kmer_count += node_length;

//             // go to the end of the current unitig
//             // TODO: check that the unitig matches the input haplotype
//             if node_length >= K::k() +1 {
//                 current_kmer = kmer_iter.nth(node_length-K::k()-1);
//                 if current_kmer.is_none() {
//                     println!("More kmers in the unitig than expected");
//                     continue;
//                 }
//             }
//             // go to the start of the next unitig, if any
//             current_kmer = kmer_iter.next();

//             // get the edges of the current unitig
//             // TODO: if current kmer is none, we expect to find a dead end, as we are at the end of the sequence
//             match node.exts().num_ext_dir(dir.flip()) {
//                 0 => { dead_ends += 1; },
//                 1 => { },
//                 _ => { junctions += 1; },
//             }
//         }
//     }
//     let duration = start.elapsed();
//     println!("done in {:?}", duration);
//     println!("Haplo contains:\n  - kmers: {}\n  - breakpoints: {}  (>1)\t\t{}  (<1)", kmer_count, junctions, dead_ends);
// }