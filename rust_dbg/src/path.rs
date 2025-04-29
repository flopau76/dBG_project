//! Path finding algorithms for de Bruijn graphs
use crate::fasta_reader::DnaRecord;
use crate::graph::Graph;

use debruijn::{Dir, Kmer, KmerIter, Mer, Vmer};
use debruijn::dna_string::DnaString;

use ahash::{AHashMap, AHashSet};
use std::collections::hash_map::Entry;
use std::error::Error;

//####################################################################################
//                              Custom errors                                       //
//####################################################################################

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
//                            Unitig iterator                                       //
//####################################################################################

/// Iterator over a dna sequence, following the nodes in the graph. This will raise an error if the sequence from kmer_iter is not present in the graph.
pub struct NodeIterator<'a, K: Kmer, D: Vmer> {
    graph: &'a Graph<K>,
    kmer_iter: KmerIter<'a, K, D>,
    current_node: Option<(usize, Dir)>,
    current_start: usize,
    current_end: usize,
    pub start_offset: usize,
    pub end_offset: Option<usize>,
}

impl<'a, K: Kmer, D: Vmer> NodeIterator<'a, K, D> {
    pub fn new(graph: &'a Graph<K>, sequence: &'a D) -> Result<Self, PathwayError<K>> {
        let kmer_iter = sequence.iter_kmers::<K>();
        let mut node_iter = Self{
            graph,
            kmer_iter,
            current_node: None,
            current_start: 0,
            current_end: 1,
            start_offset: 0,
            end_offset: None,
        };

        // initialise the iterator by looking for the first kmer
        let mut kmer = node_iter.kmer_iter.next().expect("kmer iter is empty");
        let node = search_kmer(graph, kmer, Dir::Left);

        // first kmer corresponds to the start of a node
        if node.is_some() {
            node_iter.current_node = node;
            // get the sequence of this node
            let node = node.unwrap();
            let expected_seq = match node.1 {
                Dir::Left => graph.get_node(node.0).sequence(),
                Dir::Right => graph.get_node(node.0).sequence().rc(),
            };
            // advance the kmer iterator to the end of this node and verify that both coincide
            for (idx, expected_base) in expected_seq.slice(K::k(), expected_seq.len()).iter().enumerate() {
                let kmer = match node_iter.kmer_iter.next(){
                    Some(kmer) => kmer,
                    None => {
                        node_iter.end_offset = Some(expected_seq.len()-idx);
                        break
                    },
                };
                node_iter.current_end += 1;
                if expected_base != kmer.get(K::k()-1) {
                    return Err(PathwayError::UnitigNotMatching(expected_seq.to_owned(), kmer, 1));
                }
            }
        }
        // first kmer does not correspond to the start of a node
        else {
            let mut node = search_kmer(graph, kmer, Dir::Right);
            let mut kmer_seq = DnaString::new();
            // advance in kmer_iter until we find the end of a node
            while node.is_none() {
                kmer_seq.push(kmer.get(0));
                node_iter.current_end += 1;
                kmer = match node_iter.kmer_iter.next() {
                    Some(kmer) => kmer,
                    None => {
                        // the sequence does not even correspond to a full node => we have to iterate the whole graph to find it
                        let (node, start_offset) = search_kmer_offset(graph, kmer, Dir::Left).ok_or(PathwayError::KmerNotFound(kmer))?;
                        node_iter.current_node = Some(node);
                        node_iter.start_offset = start_offset;
                        node_iter.end_offset = Some(graph.get_node(node.0).len()-start_offset-node_iter.current_end);   // TODO: check value (edge case anyway)
                        return Ok(node_iter)
                    },
                };
                node = search_kmer(graph, kmer, Dir::Right);
            }
            node_iter.current_node = node;
            let node = node.unwrap();
            // verify that the skipped kmers coincide with the beginning of the node
            let expected_seq = match node.1 {
                Dir::Left => graph.get_node(node.0).sequence(),
                Dir::Right => graph.get_node(node.0).sequence().rc(),
            };
            node_iter.start_offset = expected_seq.len()-K::k()+1-node_iter.current_end;
            let expected_seq = expected_seq.slice(node_iter.start_offset, expected_seq.len()-K::k()).to_owned();
            if expected_seq != kmer_seq {
                return Err(PathwayError::UnitigNotMatching(expected_seq, kmer, 0));
            }
        }
        Ok(node_iter)
    }

    /// get the position of the current node in the kmer iterator
    pub fn start_position(&self) -> usize {
        self.current_start
    }
    pub fn end_position(&self) -> usize {
        self.current_end
    }

    /// get the next node in the iterator, without advancing
    pub fn peek(&self) -> Option<(usize, Dir)> {
        self.current_node
    }

    /// get the next node in the iterator, advancing the iterator
    pub fn next(&mut self) -> Result<Option<(usize, Dir)>, PathwayError<K>> {
        let current = self.peek();
        self.advance()?;
        Ok(current)
    }

    /// advance the iterator
    pub fn advance(&mut self) -> Result<(), PathwayError<K>> {
        // get the next kmer in the iterator, if any
        let kmer = match self.kmer_iter.next() {
            Some(kmer) => kmer,
            None => {
                self.current_node = None;
                return Ok(())
            },
        };
        // look for the kmer at the beginning of the unitigs
        let node = search_kmer(self.graph, kmer, Dir::Left).ok_or(PathwayError::KmerNotFound(kmer))?;   // TODO: change message: the kmer might exist but not at the start of a node, as expected
        self.current_node = Some(node);
        self.current_start = self.current_end;
        self.current_end += 1;
        // get the sequence of this node
        let expected_seq = match node.1 {
            Dir::Left => self.graph.get_node(node.0).sequence(),
            Dir::Right => self.graph.get_node(node.0).sequence().rc(),
        };
        let expected_seq = expected_seq.slice(K::k(), expected_seq.len());

        // advance the iterator to the beginning of the next node and
        // check that it coincides with the whole node, not only with its first kmer
        for (idx, expected_base) in expected_seq.iter().enumerate() {
            let kmer = match self.kmer_iter.next() {
                Some(kmer) => kmer,
                None => {
                    self.end_offset = Some(expected_seq.len()-idx);
                    break
                },
            };
            self.current_end += 1;
            if expected_base != kmer.get(K::k()-1) {
                return Err(PathwayError::UnitigNotMatching(expected_seq.to_owned(), kmer, 1));
            }
        }
        Ok(())
    }
}

//####################################################################################
//                            Utility functions                                     //
//####################################################################################

// Search a kmer at an extremity of the nodes.
fn search_kmer<K: Kmer>(graph: &Graph<K>, kmer: K, side: Dir) -> Option<(usize, Dir)> {
    match side {
        Dir::Left => graph.find_link(kmer, Dir::Right).map(|(id, dir, _)| (id, dir)),
        Dir::Right => graph.find_link(kmer, Dir::Left).map(|(id, dir, _)| (id, dir.flip()))
    }
}

// Search a kmer within the nodes, by iterating over all the graph. Returns ((node_id, dir), offset) where offset is counted from the given side.
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

// Get the leftmost kmer of a node.
fn get_left_kmer<K: Kmer>(graph: &Graph<K>, node: (usize, Dir)) -> K {
    let (node_id, dir) = node;
    match dir {
        Dir::Left => graph.get_node(node_id).sequence().first_kmer(),
        Dir::Right => graph.get_node(node_id).sequence().last_kmer::<K>().rc(),
    }
}
// Get the right kmer extensions of a node.
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

// Advances the bfs by one depth and returns the new frontier set. Direction indicates in wich direction the strand is elongated.
fn advance_bfs<K: Kmer>(graph: &Graph<K>, queue: &mut Vec<(usize, Dir)>, parents: &mut AHashMap<(usize, Dir), (usize, Dir)>, direction:Dir) -> Vec<(usize, Dir)> {
    let mut next_queue = Vec::new();
    while let Some(current_node) = queue.pop() {
        for (neigh_id, neigh_dir, _) in graph.get_node(current_node.0).edges(current_node.1.cond_flip(direction==Dir::Left)) {
            let entry = parents.entry((neigh_id, neigh_dir.cond_flip(direction!=Dir::Left)));
            if let Entry::Vacant(e) = entry {
                e.insert(current_node);
                next_queue.push((neigh_id, neigh_dir.cond_flip(direction!=Dir::Left)));
            }
        }
    }
    next_queue
}

/// Search the shortest path (in number of nodes) in a compacted De Bruijn Graph using a double ended breadth-first search.  
/// Note: when several paths exist, the first one found is returned.
pub fn get_shortest_path_double_bfs<K: Kmer>(graph: &Graph<K>, start_node: (usize, Dir), end_node: (usize, Dir)) -> Result<Vec<(usize, Dir)>, PathwayError<K>> {

    // edge case: start and end are the same
    if start_node == end_node {
        return Ok(vec![(start_node)]);
    }

    // initialize the two queues
    let mut parents_left= AHashMap::default();
    let mut parents_right= AHashMap::default();
    let mut queue_left = Vec::new();
    let mut queue_right = Vec::new();
    parents_left.insert(start_node, start_node);
    parents_right.insert(end_node, end_node);
    queue_left.push(start_node);
    queue_right.push(end_node);

    // perform BFS
    let mut middle_node = None;
    'kmer: while !queue_left.is_empty() && !queue_right.is_empty() {
        if queue_left.len() < queue_right.len() {
            // elongate from the left
            queue_left = advance_bfs(graph, &mut queue_left, &mut parents_left, Dir::Left);
            // check if the two queues have met
            for node in &queue_left {
                if parents_right.contains_key(node) {
                    middle_node = Some(*node);
                    break 'kmer;
                }
            }
        } else {
            // elongate from the right
            queue_right = advance_bfs(graph, &mut queue_right, &mut parents_right, Dir::Right);
            // check if the two queues have met
            for node in &queue_right {
                if parents_left.contains_key(node) {
                    middle_node = Some(*node);
                    break 'kmer;
                }
            }
        }
    }

    if middle_node.is_none() {
        return Err(PathwayError::NoPathExists);
    }
    let middle_node = middle_node.unwrap();

    // traceback left side by baktracking from the middle node
    let mut node = middle_node;
    let mut path = Vec::new();
    path.push(node);

    while node != start_node {
        node = *parents_left.get(&node).expect("Node not found in parents map");
        path.push(node);
    }
    path.reverse();

    // traceback right side by baktracking from the middle node
    node = middle_node;
    while node != end_node {
        node = *parents_right.get(&node).expect("Node not found in parents map");
        path.push(node);
    }

    Ok(path)
}

/// Wrapper to find the shortest path in terms of nodes between two kmers.
pub fn get_shortest_path_nodes<K: Kmer>(graph: &Graph<K>, start: K, end: K) -> Result<DnaString, PathwayError<K>> {
    let (start, start_offset) = search_kmer_offset(graph, start, Dir::Left)
        .ok_or(PathwayError::KmerNotFound(start))?;
    let (end, end_offset) = search_kmer_offset(graph, end, Dir::Right)
        .ok_or(PathwayError::KmerNotFound(end))?;
    let path = get_shortest_path_double_bfs::<K>(graph, start, end)?;
    let path = graph.sequence_of_path(path.iter());
    let path = path.slice(start_offset, path.len()-end_offset);
    Ok(path.to_owned())
}

//####################################################################################
//                       Find the checkpoints                                       //
//####################################################################################

// Performs a BFS to find the next checkpoints in the graph
fn get_next_checkpoint<K: Kmer, D: Vmer>(graph: &Graph<K>, unitig_iter: &mut NodeIterator<K,D>) -> Result<Option<((usize, Dir),(usize, Dir))>, PathwayError<K>> {
    let start_node = match unitig_iter.peek() {
        Some(node) => node,
        None => return Ok(None),
    };
    let start_position = unitig_iter.start_position();
    let start_left_kmer = get_left_kmer(graph, start_node);

    let mut current_node = start_node;
    let mut next_node = match unitig_iter.next()? {
        Some(node) => node,
        None => {
            println!("1 unitigs, starting at {}:\t {:?}\t {:?}", start_position, start_node, start_node);
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
                // second time we find the expected node -> several different paths exists
                if neigh_left_kmer == next_left_kmer && !new { break; }
            }
        }
        assert!(visited.contains(&next_left_kmer));

        frontier_set = next_set;
        next_set = Vec::new();


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
    let mut unitig_iter = NodeIterator::new(graph, &seq)?;
    let mut checkpoints = Vec::new();

    while let Some(node) = get_next_checkpoint(graph, &mut unitig_iter)? {
        checkpoints.push(node);
    }
    Ok(checkpoints)
}


#[cfg(test)]
mod unit_test {
    use super::*;

    use crate::graph::Graph;

    use debruijn::kmer::Kmer3;
    use debruijn::{DnaSlice, dna_string::DnaString, Vmer, bits_to_ascii};

    const STRANDED : bool = true;
    const SEQ: DnaSlice = DnaSlice(&[2,2,2,1,1,1,1,2,2,2,0,0,0,0,0,1]);    // gggccccgggaaaaac
    const SHORTEST: DnaSlice = DnaSlice(&[2,2,2,0,0,1]); // gggaac
    // TODO: example with different result if stranded or not + example with different paths of the same length

    #[test]
    fn test_shortest_path() {
        let (start, end) = SEQ.both_term_kmer::<Kmer3>();
        let graph = Graph::from_seq_serial(SEQ, STRANDED);
        let path_bfs = get_shortest_path_nodes(&graph, start, end).unwrap();

        let shortest = DnaString::from_bytes(SHORTEST.0);
        assert!(path_bfs == shortest);
    }

    #[test]
    fn test_unitig_iterator() {
        let graph = Graph::<Kmer3>::from_seq_serial(SEQ, STRANDED);
        let mut unitig_iter = NodeIterator::new(&graph, &SEQ).unwrap();

        let mut path = Vec::new();
        while let Some(node) = unitig_iter.next().unwrap() {
            path.push(node);
            // println!("position: {}, node: {:?}", unitig_iter.position(), graph.get_node(node.0).sequence());
        }
        let path = graph.sequence_of_path(path.iter());
        let seq = DnaString::from_bytes(SEQ.0);
        assert!(path == seq);
    }

    #[test]
    #[ignore]
    fn test_get_checkpoints() {
        let graph = Graph::<Kmer3>::from_seq_serial(SEQ, STRANDED);
        let haplo = DnaRecord::new(String::from("test"), SEQ.iter().map(|x| bits_to_ascii(x)).collect::<Vec<u8>>());
        let _checkpoints = get_checkpoints_bfs(&graph, &haplo).unwrap();        
    }
}