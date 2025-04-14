//! Path finding algorithms for de Bruijn graphs
use crate::fasta_reader::DnaRecord;

use debruijn::{Dir, Kmer};
use debruijn::graph::DebruijnGraph;
use debruijn::dna_string::DnaString;

use ahash::{AHashMap, AHashSet};
use std::collections::VecDeque;
use std::collections::hash_map::Entry;
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
/// `f_dist` mesures the cost of adding a unitig to the path. it takes as an input the length of this node
// TODO: modify call to `f_dist` to enable more complex functions, depending on the directed edge and not only on its tail node
pub fn get_shortest_path<K: Kmer, D>(graph: &DebruijnGraph<K,()>, start: K, end: K, f_dist: D) -> Result<DnaString, PathwayError<K>>
where D: Fn(usize) -> usize
{
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
    parents.insert((start_id, start_dir), (0, (start_id, start_dir)));
    queue.push_back((0, start_id, start_dir));

    // perform BFS
    'kmer: while let Some((distance, node_id, dir)) = queue.pop_front() {
        let node = graph.get_node(node_id);
        let edges = node.edges(dir.flip());
        let neigh_dist = distance + f_dist(node.len());
        for (neigh_id, neigh_dir, _) in edges {
            let entry = parents.entry((neigh_id, neigh_dir));
            match entry {
                Entry::Vacant(e) => {
                    e.insert((neigh_dist, (node_id, dir)));
                    queue.push_back((neigh_dist, neigh_id, neigh_dir));
                }
                Entry::Occupied(mut e) => {
                    let (dist, _) = *e.get();
                    if dist > neigh_dist {
                        e.insert((neigh_dist, (node_id, dir)));
                        queue.push_back((neigh_dist, neigh_id, neigh_dir));
                    }
                }
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
        (_, node) = *parents.get(&node).unwrap();
    }
    path.push(node);
    path.reverse();

    Ok(graph.sequence_of_path(path.iter()))
}

/// Wrapper for finding the shortest path in terms of nodes
pub fn get_shortest_path_nodes<K: Kmer>(graph: &DebruijnGraph<K, ()>, start: K, end: K) -> Result<DnaString, PathwayError<K>> {
    get_shortest_path(graph, start, end, |_| 1)
}

/// Wrapper for finding the shortest path in terms of nucleotides
pub fn get_shortest_path_distance<K: Kmer>(graph: &DebruijnGraph<K, ()>, start: K, end: K) -> Result<DnaString, PathwayError<K>> {
    get_shortest_path(graph, start, end, |node_len| node_len)
}

/// Breaks the input haplotype into segments corresponding to shortest paths (nb of nodes) in the graph.
pub fn get_checkpoints<K: Kmer>(graph: &DebruijnGraph<K,()>, haplo: &DnaRecord) -> Result<Vec<( (usize, Dir), (usize, Dir) )>, PathwayError<K>> {
    let mut chunks = Vec::new();
    let mut kmer_iter = haplo.iter_kmers::<K>(false);

    let mut start = kmer_iter.next();

    // while not at the end of the haplotype: start a new BFS
    while start.is_some() {

        let start_kmer = start.unwrap();
        let (start_id, start_dir, _) = graph.find_link(start_kmer, Dir::Right)
        .ok_or(PathwayError::KmerNotFound(start_kmer))?;

        let mut current_node = (start_id, start_dir);
        let node_length = graph.get_node(start_id).len();
        let mut next_kmer = kmer_iter.nth(node_length-K::k());   // first kmer of the next unitig in the haplo

        let mut FS = Vec::new();        // frontier set
        let mut NS = vec![current_node];  // next set
        let mut visited = AHashSet::default();

        // while BFS coincides with the haplotype: explore next level of BFS
        while next_kmer.is_some() {

            let (next_id, next_dir, _) = graph.find_link(next_kmer.unwrap(), Dir::Right)
                .ok_or(PathwayError::KmerNotFound(next_kmer.unwrap()))?;
            // we allready visited the node so we are not on the shortest path anymore
            if visited.contains(&(next_id, next_dir)) {
                break;
            }

            current_node = (next_id, next_dir);
            let node_length = graph.get_node(next_id).len();
            next_kmer = kmer_iter.nth(node_length-K::k());

            FS = NS;
            NS = Vec::new();

            // while there are kmers in the frontier set: explore their neighbors
            while let Some((node_id, dir)) = FS.pop() {
                let node = graph.get_node(node_id);
                let edges = node.edges(dir.flip());
                for (neigh_id, neigh_dir, _) in edges {
                    let new = visited.insert((neigh_id, neigh_dir));
                    if new {
                        NS.push((neigh_id, neigh_dir));
                    }
                }
            }
        }

        println!("Found chunk");
        chunks.push(((start_id, start_dir), current_node));
        start = next_kmer;
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
        let path = get_shortest_path_nodes(&graph, start, end).unwrap();
        println!("Shortest path: {:?}", path);
    }
}