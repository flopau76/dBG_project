//! Path finding algorithms for de Bruijn graphs

use debruijn::{Dir, DnaBytes, Kmer};
use debruijn::graph::DebruijnGraph;
use debruijn::dna_string::DnaString;

use ahash::{AHashMap, AHashSet};
use std::collections::VecDeque;
use std::collections::hash_map::Entry::Vacant;
use std::error::Error;

/// Search the shortest path between two k-mers in a graph using the Breadth-First Search algorithm.
pub fn get_shortest_path<K: Kmer>(graph: &DebruijnGraph<K,()>, start: K, end: K) -> Result<DnaString, Box<dyn Error>> {
    let mut parents= AHashMap::default();
    let mut queue = VecDeque::new();

    // edge case: start and end are the same
    if start == end {
        let seq = start.iter().collect::<Vec<u8>>();
        return Ok(DnaString::from_bytes(&seq));
    }

    let (start_id, _side, start_flip) = graph.find_link(start, Dir::Right).expect("start kmer does not correspond to the beginning of a unitig");
    let (end_id, _side, end_flip) = graph.find_link(end, Dir::Left).expect("end kmer does not correspond to the ending of a unitig");

    // mark the start kmer as visited and add it to the queue
    parents.insert((start_id, start_flip), (start_id, start_flip));
    queue.push_back((start_id, start_flip));

    // perform BFS
    'kmer: while let Some((node_id, flip)) = queue.pop_front() {
        let edges = graph.get_node(node_id).edges(Dir::Right.cond_flip(flip));
        for (neigh_id, _dir, neigh_flip) in edges {
            let entry = parents.entry((neigh_id, neigh_flip));
            if let Vacant(e) = entry {
                e.insert((node_id, flip));
                queue.push_back((neigh_id, neigh_flip));
            }
            if (neigh_id, neigh_flip) == (end_id, end_flip) {
                break 'kmer;
            }
        }
    }

    if queue.is_empty() {
        return Err("No path found".into());
    }

    let mut path = Vec::new();
    let (mut kmer_id, mut flip) = (end_id, end_flip);

    while (kmer_id, flip) != (start_id, start_flip) {
        path.push((kmer_id, Dir::Left.cond_flip(flip)));
        (kmer_id, flip) = *parents.get(&(kmer_id, flip)).expect("parent should be in graph");
    }
    path.push((kmer_id, Dir::Right.cond_flip(flip)));
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
        let graph = graph::graph_from_unitigs(PATH, false);
        graph.print();
        let path = get_shortest_path(&graph, start, end).unwrap();
        println!("Shortest path: {:?}", path);
    }
}