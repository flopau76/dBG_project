//! Path finding algorithms for de Bruijn graphs

use debruijn::{Dir, DnaBytes, Kmer};
use debruijn::graph::DebruijnGraph;
use debruijn::dna_string::DnaString;

use ahash::{AHashMap, AHashSet};
use std::collections::VecDeque;
use std::collections::hash_map::Entry::Vacant;
use std::error::Error;

/// Custom error type for pathway search operations
#[derive(Debug)]
pub enum PathwayError<K: Kmer> {
    KmerNotFound(K),
    NoPathExists,
}
impl<K: Kmer> Error for PathwayError<K> {}
impl<K: Kmer> std::fmt::Display for PathwayError<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PathwayError::KmerNotFound(kmer) => write!(f, "Kmer not found: {:?}", kmer),
            PathwayError::NoPathExists => write!(f, "No path found between the given k-mers"),
        }
    }
}

/// Search the shortest path between two k-mers in a graph using the Breadth-First Search algorithm.
pub fn get_shortest_path<K: Kmer>(graph: &DebruijnGraph<K,()>, start: K, end: K) -> Result<DnaString, PathwayError<K>> {
    let mut parents= AHashMap::default();
    let mut queue = VecDeque::new();

    // edge case: start and end are the same
    if start == end {
        let seq = start.iter().collect::<Vec<u8>>();
        return Ok(DnaString::from_bytes(&seq));
    }

    let (start_id, start_dir, _flip) = graph.find_link(start, Dir::Right)
        .ok_or(PathwayError::KmerNotFound(start))?;
    let (end_id, end_dir, _flip) = graph.find_link(end, Dir::Left)
        .ok_or(PathwayError::KmerNotFound(end))?;
    let end_dir = end_dir.flip();

    // mark the start kmer as visited and add it to the queue
    parents.insert((start_id, start_dir), (start_id, start_dir));
    queue.push_back((start_id, start_dir));

    // perform BFS
    'kmer: while let Some((node_id, dir)) = queue.pop_front() {
        let edges = graph.get_node(node_id).edges(dir.flip());
        for (neigh_id, neigh_dir, _flip) in edges {
            let entry = parents.entry((neigh_id, neigh_dir));
            if let Vacant(e) = entry {
                e.insert((node_id, dir));
                queue.push_back((neigh_id, neigh_dir));
            }
            if (neigh_id, neigh_dir) == (end_id, end_dir) {
                break 'kmer;
            }
        }
    }

    if queue.is_empty() {
        return Err(PathwayError::NoPathExists);
    }

    let mut path = Vec::new();
    let mut node = (end_id, end_dir);

    while node != (start_id, start_dir) {
        path.push(node);
        node = *parents.get(&node).unwrap();
    }
    path.push(node);
    path.reverse();

    Ok(graph.sequence_of_path(path.iter()))
}


#[cfg(test)]
mod unit_test {
    use super::*;

    use crate::graph;

    use debruijn::kmer::Kmer3;
    use debruijn::{base_to_bits, DnaBytes, Vmer};

    static PATH: &str = "../data/input/test.fna";
    static SEQ: &[u8] = b"GgctGAGCTGAGTT";
    static SHORTEST_NO_C: &[u8] = b"GGCTGAGTT";
    static SHORTEST_C_1: &[u8] = b"GGCTGAGTT";
    static SHORTEST_C_2: &[u8] = b"GGCTCAGTT";

    fn make_seq() -> DnaBytes {
        DnaBytes(SEQ.iter().map(|&b| base_to_bits(b)).collect())
    }

    #[test]
    fn test_shortest_path() {
        let seq = make_seq();
        let (start, end) = seq.both_term_kmer::<Kmer3>();
        let graph = graph::graph_from_unitigs_serial(PATH, false);
        println!("{:?}", graph.base.sequences);
        graph.print();
        println!();
        println!("Start: {:?}", start);
        println!("End: {:?}", end);
        let path = get_shortest_path(&graph, start, end).unwrap();
        println!("Shortest path: {:?}", path);
    }
}