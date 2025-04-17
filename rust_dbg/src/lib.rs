/*! # rust_dbg: a crate for the manipulation of de Bruijn graphs
This crates provides methods for the construction and the manipulation of de Bruijn graphs. It builds heavily upon the crate `debruijn` developped by 10x Genomics.
*/

pub mod fasta_reader;
pub mod graph;
pub mod path;
pub mod stats;

// pub use ... for re-exports
pub use debruijn::{Kmer, Dir};
use debruijn::graph::DebruijnGraph;

pub type Graph<K> = DebruijnGraph<K, ()>;

/// Search a kmer at one extremity of the unitigs.
fn search_kmer<K: Kmer>(graph: &Graph<K>, kmer: K, side: Dir) -> Option<(usize, Dir)> {
    match side {
        Dir::Left => graph.find_link(kmer, Dir::Right).map(|(id, dir, _)| (id, dir)),
        Dir::Right => graph.find_link(kmer, Dir::Left).map(|(id, dir, _)| (id, dir.flip()))
    }
}

/// Search a kmer by iterating over all the graph.
/// Result = ((node_id, direction), offset) with:  
/// - node_id: id of the node where the kmer is found
/// - direction: direction of the node (left=normal or right=reverse complement)
/// - offset: offset of the kmer from the given side of the node
fn search_kmer_offset<K: Kmer>(graph: &Graph<K>, kmer: K, side: Dir) -> Option<((usize, Dir), usize)> {
    // search at the given extremity
    if let Some((node_id, dir)) = search_kmer(graph, kmer, side) {
        return Some(((node_id, dir), 0));
    }
    // if not found, iterate over all kmers
    let rc = kmer.rc();
    for node_id in 0..graph.len() {
        for (offset, k) in graph.get_node_kmer(node_id).into_iter().enumerate() {
            if k == kmer {
                 let offset = match side {
                    Dir::Left => offset,
                    Dir::Right => graph.get_node(node_id).len() - K::k() - offset,
                };
                return Some(((node_id, side), offset));
            }
            // if not stranded, look for the reverse complement
            else if !graph.base.stranded && k == rc {
                let offset = match side {
                    Dir::Left => graph.get_node(node_id).len() - K::k() - offset,
                    Dir::Right => offset,
                };
                return Some(((node_id, side.flip()), offset));
            }
        }
    }
    None
}