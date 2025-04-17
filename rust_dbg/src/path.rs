//! Path finding algorithms for de Bruijn graphs
use crate::fasta_reader::DnaRecord;
use crate::{search_kmer, search_kmer_offset, Graph};

use debruijn::{Dir, Kmer, Mer};
use debruijn::dna_string::DnaString;

use ahash::{AHashMap, AHashSet};
use std::collections::{VecDeque, BinaryHeap};
use std::collections::hash_map::Entry;
use std::cmp::Reverse;
use std::error::Error;
use std::io::Write;

/// Custom error type for pathway search operations
#[derive(Debug)]
pub enum PathwayError<K> {
    NoPathExists,
    KmerNotFound(K),
    UnitigNotMatching,
}
impl<K: Kmer> Error for PathwayError<K> {}
impl<K: Kmer> std::fmt::Display for PathwayError<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PathwayError::NoPathExists => write!(f, "No path found between the given k-mers"),
            PathwayError::KmerNotFound(kmer) => write!(f, "Kmer not found: {:?}", kmer),
            PathwayError::UnitigNotMatching => write!(f, "Unitig does not match the input haplotype"),
        }
    }
}

//####################################################################################
//                       Find the shortest path                                     //
//####################################################################################

/// Search the shortest path (in number of nodes) in a compacted De Bruijn Graph using the Breadth-First Search algorithm.
/// Arguments:
/// - `graph`: the de Bruijn graph to search in
/// - `start_node`: the starting node in the graph as a tuple of (node_id, start position (left -> node read from left to right))
/// - `end_node`: the ending node in the graph as a tuple of (node_id, start position (left -> node read from left to right))

// TODO: Handle cases with several shortest paths
pub fn get_shortest_path_bfs<K: Kmer>(graph: &Graph<K>, start_node: (usize, Dir), end_node: (usize, Dir)) -> Result<DnaString, PathwayError<K>> {
    let mut parents= AHashMap::default();
    let mut queue = VecDeque::new();

    // edge case: start and end are the same
    if start_node == end_node {
        let mut seq = graph.get_node(start_node.0).sequence();
        if start_node.1 == Dir::Left {
            seq = seq.rc();
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

    // reconstruct the path by backtracking the visited nodes
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

// TODO: Handle cases with several shortest paths
pub fn get_shortest_path_djk<K: Kmer, D>(graph: &Graph<K>, start_node: (usize, Dir), end_node: (usize, Dir), f_dist: D) -> Result<DnaString, PathwayError<K>>
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

    // reconstruct the path by backtracking the visited nodes
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

/// Wrapper for finding the shortest path in terms of nodes between two kmers.
/// The kmers must mark the extremities of a unitig in the graph.
pub fn get_shortest_path_nodes<K: Kmer>(graph: &Graph<K>, start: K, end: K) -> Result<DnaString, PathwayError<K>> {
    let (start, start_offset) = search_kmer_offset(graph, start, Dir::Left)
        .ok_or(PathwayError::KmerNotFound(start))?;
    let (end, end_offset) = search_kmer_offset(graph, end, Dir::Right)
        .ok_or(PathwayError::KmerNotFound(end))?;
    let path = get_shortest_path_bfs::<K>(graph, start, end)?;
    let path = path.slice(start_offset, path.len()-end_offset);
    Ok(path.to_owned())
}

/// Wrapper for finding the shortest path in terms of nucleotides between two kmers.
/// The kmers must mark the extremities of a unitig in the graph.
pub fn get_shortest_path_distance<K: Kmer>(graph: &Graph<K>, start: K, end: K) -> Result<DnaString, PathwayError<K>> {
    let (start, start_offset) = search_kmer_offset(graph, start, Dir::Left)
        .ok_or(PathwayError::KmerNotFound(start))?;
    let (end, end_offset) = search_kmer_offset(graph, end, Dir::Right)
        .ok_or(PathwayError::KmerNotFound(end))?;
    let path = get_shortest_path_djk(graph, start, end, |node_len| node_len)?;
    let path = path.slice(start_offset, path.len()-end_offset);
    println!("End offset: {}", end_offset);
    Ok(path.to_owned())
}

//####################################################################################
//                       Find the checkpoints                                       //
//####################################################################################

// Advances the kmer iterator to the start of the next node, if any.
// Raises an error if the current node does not coincide with the kmer_iterator.
fn get_next_node<K: Kmer>(graph: &Graph<K>, kmer_iter: &mut impl Iterator<Item=K>, node: (usize, Dir), offset: usize) -> Result<Option<(usize, Dir)>, PathwayError<K>> {
    let (node_id, dir) = node;
    let mut node_seq = graph.get_node(node_id).sequence();
    if dir == Dir::Right { node_seq = node_seq.rc(); }
    let nb_kmers_in_node = node_seq.len() - K::k() + 1;

    // advance the iterator to the last kmer the current node
    for i in offset+1..nb_kmers_in_node {
        let expected_base = node_seq.get(i + K::k()-1);
        let kmer = kmer_iter.next();
        if kmer.is_none() {
            return Ok(None);
        }
        let kmer = kmer.unwrap();
        if expected_base != kmer.get(K::k()-1) {
            return Err(PathwayError::UnitigNotMatching);
        }
    }

    // get the next node
    let next_kmer = kmer_iter.next();
    if next_kmer.is_none() {
        return Ok(None);
    }
    let next_kmer = next_kmer.unwrap();
    let next_node = search_kmer(graph, next_kmer, Dir::Left);
    Ok(next_node)
}

/// Breaks the input haplotype into segments corresponding to shortest paths (nb of nodes) in the graph.

// TODO: save start and end offset ?
// TODO: this will loop forever if start and end are the same
pub fn get_checkpoints_bfs<K: Kmer>(graph: &Graph<K>, haplo: &DnaRecord) -> Result<Vec<(usize, Dir)>, PathwayError<K>> {
    let mut chunks = Vec::new();
    let mut kmer_iter = haplo.iter_kmers::<K>(false);

    let start_kmer = kmer_iter.next().expect("Empty haplotype");
    let (start_node, start_offset) = search_kmer_offset(graph, start_kmer, Dir::Left)
        .ok_or(PathwayError::KmerNotFound(start_kmer))?;


    let mut current_node = start_node;
    let mut next_node = get_next_node(graph, &mut kmer_iter, start_node, start_offset)?;

    chunks.push(start_node);


    // start a new BFS
    loop {
        // initialise the BFS
        let mut frontier_set: Vec<(usize, Dir)> = vec![current_node];
        let mut next_set: Vec<(usize, Dir)> = Vec::new();
        let mut visited: AHashSet<(usize, Dir)> = AHashSet::default();
        visited.insert(current_node);

        let mut depth: usize = 0;

        // continue the current BFS
        loop {
            if next_node.is_none() {
                chunks.push(current_node);
                return Ok(chunks);
            }

            let target = next_node.unwrap();
            assert!(target != current_node);
            if visited.contains(&target) {
                chunks.push(current_node);
                println!("\tFound a cycle of size {}", depth);
                break;
            }

            // explore the next level of the BFS
            depth += 1;
            while let Some(node) = frontier_set.pop() {
                let edges = graph.get_node(node.0).edges(node.1.flip());
                for neigh_node in edges.into_iter().map(|(id, dir, _)| (id, dir)) {
                    let new = visited.insert(neigh_node);
                    if new {
                        next_set.push(neigh_node);
                    }
                }
            }
            frontier_set = next_set;
            next_set = Vec::new();

            assert!(visited.contains(&target));

            // advance the kmer iterator to the next node
            current_node = target;
            next_node = get_next_node(graph, &mut kmer_iter, current_node, 0)?;
        }
    }
    Ok(chunks)
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

    fn expected_checkpoints() -> Vec<(usize, Dir)> {
        vec![
            (0, Dir::Left),
            (3, Dir::Left),
            (2, Dir::Left),
        ]
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
        let checkpoints = get_checkpoints_bfs(&graph, &haplo).unwrap();

        // println!("{:?}", graph.base.sequences);
        // graph.print();
        // println!();
        // println!("Checkpoints: {:?}", checkpoints);

        let expected = expected_checkpoints();
        assert_eq!(checkpoints, expected);
        
    }
}