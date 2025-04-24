//! Path finding algorithms for de Bruijn graphs
use crate::fasta_reader::DnaRecord;
use crate::Graph;

use debruijn::{Dir, Kmer, Mer, Vmer, KmerIter};
use debruijn::dna_string::DnaString;

use ahash::{AHashMap, AHashSet};
use std::collections::{VecDeque, BinaryHeap};
use std::collections::hash_map::Entry;
use std::cmp::Reverse;
use std::error::Error;

/// Custom error type for pathway search operations
#[derive(Debug)]
pub enum PathwayError<K> {
    NoPathExists,
    KmerNotFound(K),
    UnitigNotMatching(DnaString, K, usize),
}
impl<K: Kmer> Error for PathwayError<K> {}
impl<K: Kmer> std::fmt::Display for PathwayError<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PathwayError::NoPathExists => write!(f, "No path found between the given k-mers"),
            PathwayError::KmerNotFound(kmer) => write!(f, "Kmer not found: {:?}", kmer),
            PathwayError::UnitigNotMatching(seq, kmer, i) => write!(f, "Expected kmer {:?} at position {} in unitig {}", kmer, i, seq),
        }
    }
}

//####################################################################################
//                            Utility functions                                     //
//####################################################################################

/// Search a kmer at one extremity of the unitigs.
fn search_kmer<K: Kmer>(graph: &Graph<K>, kmer: K, side: Dir) -> Option<(usize, Dir)> {
    match side {
        Dir::Left => graph.find_link(kmer, Dir::Right).map(|(id, dir, _)| (id, dir)),
        Dir::Right => graph.find_link(kmer, Dir::Left).map(|(id, dir, _)| (id, dir.flip()))
    }
}

/// Search a kmer by iterating over all the graph. Returns ((node_id, dir), offset) where offset is counted from the given side
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

fn get_left_kmer<K: Kmer>(graph: &Graph<K>, node: (usize, Dir)) -> K {
    let (node_id, dir) = node;
    match dir {
        Dir::Left => graph.get_node(node_id).sequence().first_kmer(),
        Dir::Right => graph.get_node(node_id).sequence().last_kmer::<K>().rc(),
    }
}
fn get_right_extensions<K: Kmer>(graph: &Graph<K>, node: (usize, Dir)) -> Vec<K> {
    let (node_id, dir) = node;
    let node = graph.get_node(node_id);
    let exts = match dir {
        Dir::Left => node.exts(),
        Dir::Right => node.exts().rc(),
    };
    let right_kmer = match dir {
        Dir::Left => node.sequence().last_kmer(),
        Dir::Right => node.sequence().first_kmer::<K>().rc(),
    };
    right_kmer.get_extensions(exts, Dir::Right)
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

    // edge case: start and end are the same
    if start_node == end_node {
        let mut seq = graph.get_node(start_node.0).sequence();
        if start_node.1 == Dir::Left {
            seq = seq.rc();
        }
        return Ok(seq.to_owned());
    }

    // mark the start kmer as visited and add it to the queue
    let start_left_kmer = get_left_kmer(graph, start_node);
    let end_left_kmer = get_left_kmer(graph, end_node);
    let mut parents= AHashMap::default();
    let mut queue = VecDeque::new();
    parents.insert(start_left_kmer, start_node);
    queue.push_back(start_left_kmer);

    // perform BFS
    'kmer: while let Some(current_left_kmer) = queue.pop_front() {
        let current_node = search_kmer(graph, current_left_kmer, Dir::Left)
            .ok_or(PathwayError::KmerNotFound(current_left_kmer))?;
        for neigh_left_kmer in get_right_extensions(graph, current_node) {
            let entry = parents.entry(neigh_left_kmer);
            if let Entry::Vacant(e) = entry {
                e.insert(current_node);
                queue.push_back(neigh_left_kmer);
            }
            if neigh_left_kmer == end_left_kmer {
                break 'kmer;
            }
        }
    }

    if queue.is_empty() {
        return Err(PathwayError::NoPathExists);
    }

    // reconstruct the path by backtracking the visited nodes
    let mut left_kmer = end_left_kmer;
    let mut path = Vec::new();
    path.push(end_node);

    while left_kmer != start_left_kmer {
        let node = *parents.get(&left_kmer).expect("Node not found in parents map");
        left_kmer = get_left_kmer(graph, node);
        path.push(node);
    }
    path.reverse();

    Ok(graph.sequence_of_path(path.iter()))
}

/// Search the shortest path in a compacted De Bruijn Graph. Contrarily to [get_shortest_path_bfs], this algorithm works with a custom distance function `f_dist(x)`,
/// which indicates how much a path is elongated when adding a node of length `x`.

// TODO: Handle cases with several shortest paths
pub fn get_shortest_path_djk<K: Kmer, D>(graph: &Graph<K>, start_node: (usize, Dir), end_node: (usize, Dir), f_dist: D) -> Result<DnaString, PathwayError<K>>
where D: Fn(usize) -> usize
{
    // edge case: start and end are the same
    if start_node == end_node {
        let seq = graph.get_node(start_node.0).sequence();
        if start_node.1 == Dir::Left {
            seq.rc();
        }
        return Ok(seq.to_owned());
    }

    // mark the start kmer as visited and add it to the queue
    let start_left_kmer = get_left_kmer(graph, start_node);
    let end_left_kmer = get_left_kmer(graph, end_node);
    let mut parents= AHashMap::default();
    let mut queue = BinaryHeap::new();
    parents.insert(start_left_kmer, (0, start_node));
    queue.push(Reverse((0, start_left_kmer)));

    // perform BFS
    while let Some(Reverse((distance, current_left_kmer))) = queue.pop() {
        if current_left_kmer == end_left_kmer {
            queue.clear();
            queue.push(Reverse((distance, current_left_kmer)));
            break;
        }
        let current_node = search_kmer(graph, current_left_kmer, Dir::Left)
            .ok_or(PathwayError::KmerNotFound(current_left_kmer))?;
        let current_node_len = graph.get_node(current_node.0).len();
        let neigh_dist = distance + f_dist(current_node_len);
        for neigh_left_kmer in get_right_extensions(graph, current_node) {
            let entry = parents.entry(neigh_left_kmer);
            match entry {
                Entry::Vacant(e) => {
                    e.insert((neigh_dist, current_node));
                    queue.push(Reverse((neigh_dist, neigh_left_kmer)));
                }
                Entry::Occupied(mut e) => {
                    let (dist, _) = *e.get();
                    if dist > neigh_dist {
                        e.insert((neigh_dist, current_node));
                        queue.push(Reverse((neigh_dist, neigh_left_kmer)));
                    }
                }
            }
        }
    }

    if queue.is_empty() {
        return Err(PathwayError::NoPathExists);
    }

    // reconstruct the path by backtracking the visited nodes
    let mut left_kmer = end_left_kmer;
    let mut path = Vec::new();
    path.push(end_node);

    while left_kmer != start_left_kmer {
        let (_, node) = *parents.get(&left_kmer).expect("Node not found in parents map");
        left_kmer = get_left_kmer(graph, node);
        path.push(node);
    }
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
    Ok(path.to_owned())
}

//####################################################################################
//                            Unitig iterator                                       //
//####################################################################################

/// Iterator over a dna sequence, following the nodes in the graph. It will raise an error if the sequence from kmer_iter can not be found in the graph.
struct UnitigIterator<'a, K: Kmer, D: Vmer> {
    graph: &'a Graph<K>,
    kmer_iter: KmerIter<'a, K, D>,
    current_node: Option<(usize, Dir)>,
    start_offset: Option<usize>,
    end_offset: Option<usize>,

    current_pos: usize,
    next_pos: usize,
}

impl<'a, K: Kmer, D: Vmer> UnitigIterator<'a, K, D> {
    fn new(graph: &'a Graph<K>, sequence: &'a D) -> Self {
        let kmer_iter = sequence.iter_kmers::<K>();
        UnitigIterator { graph, kmer_iter, current_node:None, start_offset:None, end_offset:None, current_pos:0, next_pos:0 }
    }

    /// get the start offset of the current node in the iterator
    fn position(&self) -> usize {
        self.current_pos
    }

    /// get the current node in the iterator, without advancing the iterator (except for initialisation)
    fn current(&mut self) -> Result<Option<(usize, Dir)>, PathwayError<K>> {
        if self.current_node.is_none() {
            self.current_node = self.next()?;
        }
        return Ok(self.current_node)
    }

    /// get the next node in the iterator, advancing the iterator
    fn next(&mut self) -> Result<Option<(usize, Dir)>, PathwayError<K>> {
        self.current_pos = self.next_pos;
        // get the next kmer in the kmer iterator
        let kmer = match self.kmer_iter.next() {
            Some(kmer) => kmer,
            None => {
                self.current_node = None;
                return Ok(None)
            },
        };
        self.next_pos += 1;

        // get the corresponding node
        let node: (usize, Dir);
        let offset: usize;
        // if it is the first kmer of the sequence, it may have an offset
        if self.start_offset.is_none() {
            (node, offset) = search_kmer_offset(self.graph, kmer, Dir::Left)
                .ok_or(PathwayError::KmerNotFound(kmer))?;
            self.start_offset = Some(offset);
        }
        // otherwise the kmer must correspond to the beginning of a unitig
        else {
            node = search_kmer(self.graph, kmer, Dir::Left)
                .ok_or(PathwayError::KmerNotFound(kmer))?;
            offset = 0;
        }

        // get the sequence of this node
        let node_seq = match node.1 {
            Dir::Left => self.graph.get_node(node.0).sequence(),
            Dir::Right => self.graph.get_node(node.0).sequence().rc(),
        };
        let expected_seq = node_seq.slice(K::k() + offset, node_seq.len());

        // skip the kmers in the iterator until the end of the unitig and check that the bases coincide
        for (idx, expected_base) in expected_seq.iter().enumerate() {
            let kmer = match self.kmer_iter.next() {
                Some(kmer) => kmer,
                None => {
                    self.end_offset = Some(expected_seq.len()-idx);
                    break
                },
            };
            self.next_pos += 1;
            if expected_base != kmer.get(K::k()-1) {
                return Err(PathwayError::UnitigNotMatching(node_seq.to_owned(), kmer, offset+1));
            }
        }
        self.current_node = Some(node);
        Ok(Some(node))
    }
}

//####################################################################################
//                       Find the checkpoints                                       //
//####################################################################################

// performs a BFS to find the next checkpoints in the graph
fn get_next_checkpoint<K: Kmer, D: Vmer>(graph: &Graph<K>, unitig_iter: &mut UnitigIterator<K,D>) -> Result<Option<((usize, Dir),(usize, Dir))>, PathwayError<K>> {
    let start_node = match unitig_iter.current()? {
        Some(node) => node,
        None => return Ok(None),
    };
    let start_position = unitig_iter.position();
    let start_left_kmer = get_left_kmer(graph, start_node);

    let mut current_node = start_node;
    let mut next_node = match unitig_iter.next()? {
        Some(node) => node,
        None => {
            println!("1 unitig, starting at {}:\t {:?}\t {:?}", start_position, start_node, start_node);
            return Ok(Some((start_node, start_node)))
        },
    };
    let mut next_left_kmer = get_left_kmer(graph, next_node);

    // mark the start kmer as visited and add it to the queue
    let mut frontier_set: Vec<K> = vec![start_left_kmer];
    let mut next_set: Vec<K> = Vec::new();
    let mut visited: AHashSet<K> = AHashSet::default();
    visited.insert(start_left_kmer);

    // perform the BFS
    let mut depth: usize = 1;
    while !visited.contains(&next_left_kmer) {
        // explore the next level of the BFS
        while let Some(left_kmer) = frontier_set.pop() {
            let node = search_kmer(graph, left_kmer, Dir::Left)
                .ok_or(PathwayError::KmerNotFound(left_kmer))?;
            for neigh_left_kmer in get_right_extensions(graph, node) {
                let new = visited.insert(neigh_left_kmer);
                if new {
                    next_set.push(neigh_left_kmer);
                }
            }
        }
        frontier_set = next_set;
        next_set = Vec::new();

        assert!(visited.contains(&next_left_kmer));

        // advance the kmer iterator to the next node
        current_node = next_node;
        next_node = match unitig_iter.next()? {
            Some(node) => node,
            None => break
        };
        next_left_kmer = get_left_kmer(graph, next_node);
        depth += 1;
    }
    println!("{} unitigs, starting at {}:\t {:?}\t {:?}", depth, start_position, start_node, current_node);
    Ok(Some((start_node, current_node)))
}

/// Breaks the input haplotype into segments corresponding to shortest paths (nb of nodes) in the graph.
pub fn get_checkpoints_bfs<K: Kmer>(graph: &Graph<K>, haplo: &DnaRecord) -> Result<Vec<((usize, Dir),(usize, Dir))>, PathwayError<K>> {
    let seq = haplo.dna_string();
    let mut unitig_iter = UnitigIterator::new(graph, &seq);
    let mut checkpoints = Vec::new();

    while let Some(node) = get_next_checkpoint(graph, &mut unitig_iter)? {
        checkpoints.push(node);
    }
    Ok(checkpoints)
}


#[cfg(test)]
mod unit_test {
    use super::*;

    use crate::graph;

    use debruijn::kmer::Kmer3;
    use debruijn::{base_to_bits, DnaBytes, Vmer};

    static PATH_S: &str = "../data/input/test_stranded.fna";
    static PATH_NS: &str = "../data/input/test_not_stranded.fna";
    static SEQ: &[u8] = b"gggccccgggaaaaa";
    static SHORTEST_NO_C: &[u8] = b"gggaaa";

    fn make_seq() -> DnaBytes {
        DnaBytes(SEQ.iter().map(|&b| base_to_bits(b)).collect())
    }

    #[test]
    fn test_shortest_path() {
        let seq = make_seq();
        let (start, end) = seq.both_term_kmer::<Kmer3>();
        let graph = graph::graph_from_unitigs_serial(PATH_NS, false);
        let path_bfs = get_shortest_path_nodes(&graph, start, end).unwrap();
        let path_djk = get_shortest_path_distance(&graph, start, end).unwrap();

        assert!(path_bfs == DnaString::from_acgt_bytes(SHORTEST_NO_C));
        assert!(path_djk == DnaString::from_acgt_bytes(SHORTEST_NO_C));

        // println!("{:?}", graph.base.sequences);
        // graph.print();
        // println!();
        // println!("Start: {:?}", start);
        // println!("End: {:?}", end);
        // println!("Shortest path BFS: {:?}", path_bfs);
    }

    #[test]
    fn test_unitig_iterator() {
        let haplo = DnaRecord::new(String::from("test"), SEQ.to_vec());
        let dna_sting = haplo.dna_string();
        let graph = graph::graph_from_unitigs_serial::<Kmer3>(PATH_NS, false);
        let mut unitig_iter = UnitigIterator::new(&graph, &dna_sting);

        println!("{:?}", graph.base.sequences);
        graph.print();
        println!();

        let mut path = Vec::new();
        while let Some(node) = unitig_iter.next().unwrap() {
            path.push(node);
            println!("position: {}, node: {:?}", unitig_iter.position(), graph.get_node(node.0).sequence());
        }
        println!("{:?}", path);
    }

    #[test]
    fn test_get_checkpoints() {
        let haplo = DnaRecord::new(String::from("test"), SEQ.to_vec());
        let graph = graph::graph_from_unitigs_serial::<Kmer3>(PATH_S, true);
        let checkpoints = get_checkpoints_bfs(&graph, &haplo).unwrap();

        // println!("{:?}", graph.base.sequences);
        // graph.print();
        // println!();
        println!("Checkpoints: {:?}", checkpoints);
        
    }
}