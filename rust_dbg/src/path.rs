//! Path finding algorithms for de Bruijn graphs
use crate::fasta_reader::DnaRecord;

use debruijn::{Dir, Kmer, Mer};
use debruijn::graph::DebruijnGraph;
use debruijn::dna_string::DnaString;

use ahash::{AHashMap, AHashSet};
use std::collections::{VecDeque, BinaryHeap};
use std::collections::hash_map::Entry;
use std::cmp::Reverse;
use std::error::Error;

/// Custom error type for pathway search operations
#[derive(Debug)]
pub enum PathwayError<K: Kmer> {
    FirstKmerNotFound(K),
    LastKmerNotFound(K),
    NoPathExists,
}
impl<K: Kmer> Error for PathwayError<K> {}
impl<K: Kmer> std::fmt::Display for PathwayError<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PathwayError::FirstKmerNotFound(kmer) => write!(f, "First Kmer not found in the graph: {:?}", kmer),
            PathwayError::LastKmerNotFound(kmer) => write!(f, "Last Kmer not found in the graph: {:?}", kmer),
            PathwayError::NoPathExists => write!(f, "No path found between the given k-mers"),
        }
    }
}

/// Search the shortest path (in number of nodes) in a compacted De Bruijn Graph using the Breadth-First Search algorithm.
/// Arguments:
/// - `graph`: the de Bruijn graph to search in
/// - `start_node`: the starting node in the graph as a tuple of (node_id, start position (left -> node read from left to right))
/// - `end_node`: the ending node in the graph as a tuple of (node_id, start position (left -> node read from left to right))
pub fn get_shortest_path_bfs<K: Kmer>(graph: &DebruijnGraph<K,()>, start_node: (usize, Dir), end_node: (usize, Dir)) -> Result<DnaString, PathwayError<K>> {
    let mut parents= AHashMap::default();
    let mut queue = VecDeque::new();

    // edge case: start and end are the same
    if start_node == end_node {
        let seq = graph.get_node(start_node.0).sequence();
        if start_node.1 == Dir::Left {
            seq.rc();
        }
        return Ok(seq.to_owned());
    }

    // mark the start kmer as visited and add it to the queue
    parents.insert(start_node, start_node);
    queue.push_back(start_node);

    // perform BFS
    'kmer: while let Some(current_node) = queue.pop_front() {
        let node = graph.get_node(current_node.0);
        let edges = node.edges(current_node.1.flip());
        for neigh_node in edges.into_iter().map(|(id, dir, _)| (id, dir)) {
            let entry = parents.entry(neigh_node);
            if let Entry::Vacant(e) = entry {
                e.insert( current_node);
                queue.push_back(neigh_node);
            }
            if neigh_node == end_node {
                break 'kmer;
            }
        }
    }

    if queue.is_empty() {
        return Err(PathwayError::NoPathExists);
    }

    let mut path = Vec::new();
    let mut node = end_node;

    while node != start_node {
        path.push(node);
        node = *parents.get(&node).unwrap();
    }
    path.push(node);
    path.reverse();

    Ok(graph.sequence_of_path(path.iter()))
}

/// Search the shortest path in a compacted De Bruijn Graph. Contrarily to [get_shortest_path_bfs], this algorithm works with a custom distance function `f_dist(x)`,
/// which indicates how much a path is elongated when adding a node of length `x`.
pub fn get_shortest_path_djk<K: Kmer, D>(graph: &DebruijnGraph<K,()>, start_node: (usize, Dir), end_node: (usize, Dir), f_dist: D) -> Result<DnaString, PathwayError<K>>
where D: Fn(usize) -> usize
{
    let mut parents= AHashMap::default();
    let mut queue = BinaryHeap::new();

    // edge case: start and end are the same
    if start_node == end_node {
        let seq = graph.get_node(start_node.0).sequence();
        if start_node.1 == Dir::Left {
            seq.rc();
        }
        return Ok(seq.to_owned());
    }

    // mark the start kmer as visited and add it to the queue
    parents.insert(start_node, (0, start_node));
    queue.push(Reverse((0, start_node)));

    // perform BFS
    while let Some(Reverse((distance, current_node))) = queue.pop() {
        if current_node == end_node {
            queue.clear();
            queue.push(Reverse((distance, current_node)));
            break;
        }
        let node = graph.get_node(current_node.0);
        let edges = node.edges(current_node.1.flip());
        let neigh_dist = distance + f_dist(node.len());
        for neigh_node in edges.into_iter().map(|(id, dir, _)| (id, dir)) {
            let entry = parents.entry(neigh_node);
            match entry {
                Entry::Vacant(e) => {
                    e.insert((neigh_dist, current_node));
                    queue.push(Reverse((neigh_dist, neigh_node)));
                }
                Entry::Occupied(mut e) => {
                    let (dist, _) = *e.get();
                    if dist > neigh_dist {
                        e.insert((neigh_dist, current_node));
                        queue.push(Reverse((neigh_dist, neigh_node)));
                    }
                }
            }
        }
    }

    if queue.is_empty() {
        return Err(PathwayError::NoPathExists);
    }

    let mut path = Vec::new();
    let mut node = end_node;

    while node != start_node {
        path.push(node);
        (_, node) = *parents.get(&node).unwrap();
    }
    path.push(node);
    path.reverse();

    Ok(graph.sequence_of_path(path.iter()))
}

// search a kmer by iterating over all the graph.
// Result = (node_id, direction) where (direction == left) => kmer found when reading from left to right
fn search_kmer<K: Kmer>(graph: &DebruijnGraph<K, ()>, kmer: K) -> Option<(usize, Dir)> {
    let rc = kmer.rc();
    for node_id in 0..graph.len() {
        for k in graph.get_node_kmer(node_id) {
            if k == kmer {
                return Some((node_id, Dir::Left));
            } else if k == rc {
                return Some((node_id, Dir::Right));
            }
        }
    }
    None
}
// search a kmer at the beginning of a unitig. if not found, search by iterating the whole graph.
fn get_start_node<K: Kmer>(graph: &DebruijnGraph<K, ()>, kmer: K) -> Result<(usize, Dir), PathwayError<K>> {
    graph.find_link(kmer, Dir::Right).map(|(id, dir, _)| (id, dir))
        .or_else(|| search_kmer(graph, kmer))
        .ok_or(PathwayError::FirstKmerNotFound(kmer))
}
// search a kmer at the end of a unitig. if not found, search by iterating the whole graph.
fn get_end_node<K: Kmer>(graph: &DebruijnGraph<K, ()>, kmer: K) -> Result<(usize, Dir), PathwayError<K>> {
    graph.find_link(kmer, Dir::Left).map(|(id, dir, _)| (id, dir.flip()))
        .or_else(|| search_kmer(graph, kmer))
        .ok_or(PathwayError::LastKmerNotFound(kmer))
}

/// Wrapper for finding the shortest path in terms of nodes between two kmers.
/// The kmers must mark the extremities of a unitig in the graph.
pub fn get_shortest_path_nodes<K: Kmer>(graph: &DebruijnGraph<K, ()>, start: K, end: K) -> Result<DnaString, PathwayError<K>> {
    let start = get_start_node(graph, start)?;
    let end = get_end_node(graph, end)?;
    get_shortest_path_bfs::<K>(graph, start, end)
}

/// Wrapper for finding the shortest path in terms of nucleotides between two kmers.
/// The kmers must mark the extremities of a unitig in the graph.
pub fn get_shortest_path_distance<K: Kmer>(graph: &DebruijnGraph<K, ()>, start: K, end: K) -> Result<DnaString, PathwayError<K>> {
    let start = get_start_node(graph, start)?;
    let end = get_end_node(graph, end)?;
    get_shortest_path_djk(graph, start, end, |node_len| node_len)
}

// Get the node corresponding to the next kmer in the haplotype
// and advances the kmer iterator to the end of the node
// returns None if the iterator is empty, and an error if the kmer is not found in the graph
fn get_next_node<K: Kmer>(graph: &DebruijnGraph<K, ()>, kmer_iter: &mut impl Iterator<Item=K>) -> Result<Option<(usize, Dir)>, PathwayError<K>> {
    let kmer = kmer_iter.next();
    if kmer.is_none() {
        return Ok(None);
    }
    let kmer = kmer.unwrap();
    let node = get_start_node(graph, kmer)?;
    let node_length = graph.get_node(node.0).len();
    if node_length > K::k() {
        kmer_iter.nth(node_length-K::k()-1);
        // TODO: check that the sequence of the current unitig corresponds to the skipped kmers in the haplotype ?
    }
    Ok(Some(node))
}

/// Breaks the input haplotype into segments corresponding to shortest paths (nb of nodes) in the graph.
pub fn get_checkpoints_bfs<K: Kmer>(graph: &DebruijnGraph<K,()>, haplo: &DnaRecord) -> Result<Vec<(usize, Dir)>, PathwayError<K>> {
    let mut chunks = Vec::new();
    let mut kmer_iter = haplo.iter_kmers::<K>(false);

    let mut current_node;
    let mut next_node = get_next_node(graph, &mut kmer_iter)?;
    chunks.push(next_node.expect("Haplotype is empty"));

    // start a new BFS
    loop {
        // initialise the BFS
        current_node = next_node.unwrap();
        let mut frontier_set: Vec<(usize, Dir)> = Vec::new();
        let mut next_set: Vec<(usize, Dir)> = vec![current_node];
        let mut visited: AHashSet<(usize, Dir)> = AHashSet::default();
        visited.insert(current_node);

        // continue the current BFS
        loop {
            // advance in the haplotype
            current_node = next_node.unwrap();
            next_node = get_next_node(graph, &mut kmer_iter)?;
            if next_node.is_none() {
                chunks.push(current_node);
                return Ok(chunks);
            }
            let next_node = next_node.unwrap();

            // if the next node was allready visited: end the current BFS and start a new one
            if visited.contains(&next_node) {
                chunks.push(current_node);
                break;
            }

            // explore the next level of the BFS
            frontier_set = next_set;
            next_set = Vec::new();
            while let Some(node) = frontier_set.pop() {
                let edges = graph.get_node(node.0).edges(node.1.flip());
                for neigh_node in edges.into_iter().map(|(id, dir, _)| (id, dir)) {
                    let new = visited.insert(neigh_node);
                    if new {
                        next_set.push(neigh_node);
                    }
                }
            }

            // should not happen: the next node was not visited in this level. Some edges are missing ?
            // assert!(visited.contains(&next_node), "next node not found in the frontier set");
        }
    }
}


#[cfg(test)]
mod unit_test {
    use super::*;

    use crate::graph;

    use debruijn::kmer::Kmer3;
    use debruijn::{base_to_bits, DnaBytes, Vmer};

    static PATH: &str = "../data/input/test.fna";
    static SEQ: &[u8] = b"GgctGAGCTGAGTT";
    static _SHORTEST_NO_C: &[u8] = b"GGCTGAGTT";
    static _SHORTEST_C_1: &[u8] = b"GGCTGAGTT";
    static _SHORTEST_C_2: &[u8] = b"GGCTCAGTT";

    fn make_seq() -> DnaBytes {
        DnaBytes(SEQ.iter().map(|&b| base_to_bits(b)).collect())
    }

    #[test]
    fn test_shortest_path() {
        let seq = make_seq();
        let (start, end) = seq.both_term_kmer::<Kmer3>();
        let graph = graph::graph_from_unitigs_serial(PATH, true);
        let path = get_shortest_path_distance(&graph, start, end).unwrap();

        println!("{:?}", graph.base.sequences);
        graph.print();
        println!();
        println!("Start: {:?}", start);
        println!("End: {:?}", end);
        println!("Shortest path: {:?}", path);

        assert!(path == DnaString::from_acgt_bytes(_SHORTEST_NO_C));
    }

    #[test]
    fn test_get_checkpoints() {
        let haplo = DnaRecord::new(String::from("test"), SEQ.to_vec());
        let graph = graph::graph_from_unitigs_serial::<Kmer3>(PATH, true);

        println!("{:?}", graph.base.sequences);
        graph.print();
        println!();
        let chunks = get_checkpoints_bfs(&graph, &haplo).unwrap();
        println!("Checkpoints: {:?}", chunks);
    }
}