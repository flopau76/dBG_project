//! To create graphs from fasta/unitig files

use smallvec::SmallVec;
use std::{collections::HashSet, error::Error, path::Path};

use boomphf::hashmap::NoKeyBoomHashMap;
use needletail::parse_fastx_file;
use packed_seq::{PackedSeq, PackedSeqVec, Seq, SeqVec};

use crate::KmerStorage;

mod node_iterator;
pub use node_iterator::NodeIterator;

pub mod shortest_path;

use epserde::prelude::*;

//####################################################################################
//                              Custom errors                                       //
//####################################################################################

/// Custom error type for pathway search operations
#[derive(Debug)]
pub enum PathwayError {
    KmerNotFound(String),
    UnitigNotMatching,
}
impl Error for PathwayError {}
impl std::fmt::Display for PathwayError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PathwayError::KmerNotFound(kmer) => write!(f, "Kmer not found: {}", kmer),
            PathwayError::UnitigNotMatching => write!(f, "Unitig not matching"),
        }
    }
}

//####################################################################################
//                              SequenceSet                                         //
//####################################################################################

/// A 2-bit packed owned set of multiple DNA sequences.
#[derive(Epserde, Debug, Default)]
pub struct SequenceSet {
    sequences: PackedSeqVec,
    start_positions: Vec<usize>,
    lengths: Vec<u32>,
    // Note: if we want nodes to start at the begining of an u8, there might be a gap of 1-6 nucleotides between two sequences.
    // If it is not the case, there might be undefined behaviour or panic later on.
}

impl SequenceSet {
    /// Add a packed sequence to the set.
    pub fn push_seq(&mut self, seq: PackedSeq) -> () {
        let range = self.sequences.push_seq(seq);
        self.start_positions.push(range.start);
        self.lengths.push(seq.len() as u32);
    }

    /// Add an ascii sequence to the set.
    /// Undefined behavior if the sequence contains characters other than ACGTacgt.
    pub fn push_ascii(&mut self, seq: &[u8]) -> () {
        let range = self.sequences.push_ascii(seq);
        self.start_positions.push(range.start);
        self.lengths.push(seq.len() as u32);
    }

    /// Get a PackedSeq containing the 'i'th sequence in the set.
    pub fn get(&self, i: usize) -> PackedSeq {
        self.sequences
            .slice(self.start_positions[i]..self.start_positions[i] + self.lengths[i] as usize)
    }

    /// Get the number of sequences in the set.
    pub fn len(&self) -> usize {
        self.start_positions.len()
    }
}

/// Oriented node in a canonical de Bruijn graph
#[derive(Epserde, Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(C)]
#[zero_copy]
pub struct Node {
    pub id: usize,
    pub is_rc: bool,
}

impl Node {
    pub fn new(id: usize, is_rc: bool) -> Self {
        Self { id, is_rc }
    }
}

#[derive(Epserde, Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
#[zero_copy]
pub enum Side {
    Left,
    Right,
}

impl Side {
    /// Get the opposite side.
    pub fn opposite(self) -> Self {
        match self {
            Side::Left => Side::Right,
            Side::Right => Side::Left,
        }
    }
    /// Choose one of the two possibilities, depending on the orientation of the node.
    pub fn choose<D>(self, left: D, right: D) -> D {
        match self {
            Side::Left => left,
            Side::Right => right,
        }
    }
}

//####################################################################################
//                                  BaseGraph                                       //
//####################################################################################

/// A base graph, containing only unitigs. No edges nor hash-functions.
/// This is only an itermediate step used to build the final graph.
#[derive(Epserde, Debug)]
pub struct BaseGraph {
    pub k: usize,
    pub stranded: bool,
    pub sequences: SequenceSet,
}

impl BaseGraph {
    /// Get the number of nodes in the graph.
    #[inline]
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }
    #[inline]
    pub fn stranded(&self) -> bool {
        self.stranded
    }

    /// Get the sequence of a node.
    #[inline]
    #[deprecated = "PackedSeq.to_vec() is not working properly if the underlying data representation possesses an offset"]
    pub fn node_seq(&self, node: Node) -> PackedSeqVec {
        if node.is_rc {
            self.sequences.get(node.id).to_revcomp()
        } else {
            self.sequences.get(node.id).to_vec()
        }
    }

    /// Get the kmer at the given `side` of a `node`.
    pub fn node_kmer<KS: KmerStorage>(&self, node: Node, side: Side) -> KS {
        match node.is_rc {
            true => {
                let seq_vec = self.sequences.get(node.id).to_revcomp();
                let seq = seq_vec.as_slice();
                let pos = side.choose(0, seq.len() - self.k());
                KS::get_kmer(self.k(), seq, pos)
            }
            false => {
                let seq = self.sequences.get(node.id);
                let pos = side.choose(0, seq.len() - self.k());
                KS::get_kmer(self.k(), seq, pos)
            }
        }
    }

    /// Create a graph from a fasta file containing unitigs (as returned by ggcat for example).
    pub fn from_unitig_file(path: &Path, k: usize, stranded: bool) -> Self {
        // iterate over the unitigs and add them to the sequence set
        let mut sequences = SequenceSet::default();
        let mut fasta_reader = parse_fastx_file(path).unwrap();
        while let Some(record) = fasta_reader.next() {
            let record = record.unwrap();
            let seq = record.seq();
            sequences.push_ascii(&seq);
        }
        // create the graph
        Self {
            k,
            stranded,
            sequences,
        }
    }

    /// Create a (non-compacted) graph from a packed sequence. (For debugging mainly)
    pub fn from_seq<'a>(seq: PackedSeqVec, k: usize, stranded: bool) -> Self {
        if k > 32 {
            panic!("Unsuported kmer size: {}. Maximum is 32.", k);
        }
        // get first occurence of unique kmers
        let mut kmer_set = HashSet::new();
        let mut sequences = SequenceSet::default();
        for pos in 0..seq.len() - k + 1 {
            let kmer = seq.slice(pos..pos + k);
            let mut kmer_u64 = kmer.as_u64();
            if !stranded {
                let rc_kmer = kmer.to_revcomp();
                let rc_kmer_u64 = rc_kmer.as_slice().as_u64();
                kmer_u64 = kmer_u64.min(rc_kmer_u64);
            }
            if kmer_set.insert(kmer_u64) {
                sequences.push_seq(kmer);
            }
        }
        // create the graph
        Self {
            k,
            stranded,
            sequences,
        }
    }

    /// Create a (non-compacted) graph from an ascii sequence. (For debugging mainly)
    pub fn from_seq_ascii<'a>(seq: &[u8], k: usize, stranded: bool) -> Self {
        let seq = PackedSeqVec::from_ascii(seq);
        Self::from_seq(seq, k, stranded)
    }

    pub fn finish<KS: KmerStorage>(self) -> Graph<KS> {
        // compute the edges
        // TODO

        // compute the kmers
        let mut left_kmers = Vec::with_capacity(self.len());
        let mut right_kmers = Vec::with_capacity(self.len());
        for i in 0..self.len() {
            let node = Node::new(i, false);
            let left_kmer = self.node_kmer(node, Side::Left);
            let right_kmer = self.node_kmer(node, Side::Right);
            left_kmers.push(left_kmer);
            right_kmers.push(right_kmer);
        }
        let left_ids: Vec<usize> = (0..self.len()).collect();
        let right_ids: Vec<usize> = (0..self.len()).collect();

        // create the mphf for left and right kmers
        let left_kmers_mphf = NoKeyBoomHashMap::new(left_kmers, left_ids);
        let right_kmers_mphf = NoKeyBoomHashMap::new(right_kmers, right_ids);

        // create the graph
        Graph {
            base: self,
            edges: vec![],
            left_kmers: left_kmers_mphf,
            right_kmers: right_kmers_mphf,
        }
    }
}

//####################################################################################
//                                  Graph                                           //
//####################################################################################

/// A de Bruijn graph, containing unitigs and edges.
#[derive(Epserde, Debug)]
pub struct Graph<KS: KmerStorage> {
    base: BaseGraph,
    edges: Vec<usize>,
    left_kmers: NoKeyBoomHashMap<KS, usize>,
    right_kmers: NoKeyBoomHashMap<KS, usize>,
}

/// For Graph manipulation
impl<KS: KmerStorage> Graph<KS> {
    /// Return the kmer size of the graph.
    #[inline]
    pub fn k(&self) -> usize {
        self.base.k()
    }
    /// Return wether the graph is stranded or not.
    #[inline]
    pub fn stranded(&self) -> bool {
        self.base.stranded()
    }
    /// Get the sequence of a `node`.
    #[inline]
    pub fn node_seq(&self, node: Node) -> PackedSeqVec {
        self.base.node_seq(node)
    }
    /// Get the kmer of a `node` at a given `side`.
    #[inline]
    pub fn node_kmer(&self, node: Node, side: Side) -> KS {
        self.base.node_kmer(node, side)
    }

    /// Get the neighbors of a `node` on a given `side`.
    pub fn node_neigh(&self, node: Node, side: Side) -> SmallVec<[Node; 4]> {
        let mut neighbors = SmallVec::new();
        let kmer = self.node_kmer(node, side);

        // TODO: use self.edges to reduce the nb of calls to search_kmer
        for base in 0..4 {
            let mut new_kmer = kmer.clone();
            new_kmer.extend_side(self.k(), base, side);
            if let Some(neigh) = self.search_kmer(new_kmer, side.opposite()) {
                neighbors.push(neigh);
            }
        }
        neighbors
    }

    /// Get the sequence of a `path`.
    pub fn path_seq(&self, path: &Vec<Node>) -> PackedSeqVec {
        let mut seq = self.node_seq(*path.first().expect("path is empty"));
        for node in path {
            let node_seq = self.node_seq(*node);
            seq.push_seq(node_seq.slice(self.k() - 1..node_seq.len()));
        }
        seq
    }

    /// Search a `kmer` at one `side` of the nodes.
    pub fn search_kmer(&self, kmer: KS, side: Side) -> Option<Node> {
        println!("Searching kmer: {kmer:?} on side: {side:?}");
        let map = side.choose(&self.left_kmers, &self.right_kmers);
        if let Some(id) = map.get(&kmer) {
            println!("Found kmer in mphf: {id}");
            let node = Node::new(*id, false);
            if self.node_kmer(node, side) == kmer {
                return Some(node);
            }
        }
        if !self.stranded() {
            let rc_kmer = kmer.rc(self.k());
            let map = side.choose(&self.right_kmers, &self.left_kmers);
            if let Some(id) = map.get(&rc_kmer) {
                let node = Node::new(*id, true);
                if self.node_kmer(node, side) == rc_kmer {
                    return Some(node);
                }
            }
        }
        None
    }

    /// Search a `kmer` by iterating over all the nodes (this may take some time!).
    /// Returns the node and the offset of the kmer in the node (from the given `side`).
    pub fn search_kmer_offset(&self, kmer: KS, side: Side) -> Option<(Node, usize)> {
        // search first at the given extremity
        if let Some(node) = self.search_kmer(kmer, side) {
            return Some((node, 0));
        }
        // if not found, iterate over all kmers
        unimplemented!();
    }
}

// // Dump and load from binary
// impl<KS: KmerStorage + Serialize + Deserialize> Graph<KS> {
//     /// Write the graph to a binary file.
//     pub fn save_to_binary(&self, path: &Path) -> Result<(), Box<dyn Error>> {
//         let mut file = std::fs::File::create(path)?;
//         self.serialize(&mut file, self)?;
//         Ok(())
//     }
// }

#[cfg(test)]
mod test_packed_seq {
    use super::*;

    #[test]
    fn push_seq() {
        let seq1 = PackedSeqVec::from_ascii(b"ACGTACGTTT");
        let seq2 = PackedSeqVec::from_ascii(b"TTTT");
        let mut new_seq = PackedSeqVec::default();

        let range_1 = new_seq.push_seq(seq1.as_slice());
        let range_2 = new_seq.push_seq(seq2.as_slice());

        let seq1b = new_seq.slice(range_1).unpack();
        let seq2b = new_seq.slice(range_2).unpack();

        assert_eq!(seq1b, seq1.as_slice().unpack());
        assert_eq!(seq2b, seq2.as_slice().unpack());
    }

    #[test]
    fn push_ascii() {
        let seq1 = b"ACGTACGTTT";
        let seq2 = b"TTTT";
        let mut new_seq = PackedSeqVec::default();

        let range_1 = new_seq.push_ascii(seq1);
        let range_2 = new_seq.push_ascii(seq2);

        let seq1b = new_seq.slice(range_1).unpack();
        let seq2b = new_seq.slice(range_2).unpack();

        assert_eq!(seq1b, seq1.as_slice());
        assert_eq!(seq2b, seq2.as_slice());
    }
}

#[cfg(test)]
mod unit_test {
    use super::*;

    const STRANDED: bool = false;
    const SEQ: &[u8; 16] = b"gggccccgggaaaaac";

    const K: usize = 3;
    type KS = u8;

    #[ignore]
    #[test]
    fn graph() {
        let base = BaseGraph::from_seq_ascii(SEQ, K, STRANDED);
        let graph: Graph<KS> = base.finish();

        for i in 0..graph.base.len() {
            let node = Node::new(i, false);
            let l_kmer = graph.node_kmer(node, Side::Left);
            let r_kmer = graph.node_kmer(node, Side::Right);
            let l = graph.search_kmer(l_kmer, Side::Left);
            let r = graph.search_kmer(r_kmer, Side::Right);
            assert_eq!(l, Some(Node::new(i, false)));
            assert_eq!(r, Some(Node::new(i, false)));
        }
    }

    #[ignore]
    #[test]
    fn node_iterator_offset() {
        let base = BaseGraph::from_seq_ascii(SEQ, K, STRANDED);
        let graph: Graph<KS> = base.finish();
        println!("Graph: {:?}", graph);

        let seq_start = PackedSeqVec::from_ascii(b"aacc");
        let mut unitig_iter = NodeIterator::new(&graph, seq_start).unwrap();
        while let Some(node) = unitig_iter.next().unwrap() {
            println!("Node: {:?}", node);
        }
        assert_eq!(unitig_iter.start_offset, Some(0));
        assert_eq!(unitig_iter.end_offset, Some(4));

        let seq_end = PackedSeqVec::from_ascii(b"ggtt");
        let mut unitig_iter = NodeIterator::new(&graph, seq_end).unwrap();
        while let Some(node) = unitig_iter.next().unwrap() {
            println!("Node: {:?}", node);
        }
        assert_eq!(unitig_iter.start_offset, Some(4));
        assert_eq!(unitig_iter.end_offset, Some(0));

        // // unimplemented: if unitig does not correspond to the extremity of a node, it is not detected
        // let seq_middle = PackedSeqVec::from_ascii(b"ccgg");
        // let mut unitig_iter = NodeIterator::new(&graph, seq_middle).unwrap();
        // while let Some(node) = unitig_iter.next().unwrap() {
        //     println!("Node: {:?}", node);
        // }
        // assert_eq!(unitig_iter.start_offset, Some(2));
        // assert_eq!(unitig_iter.end_offset, Some(2));
    }

    // #[test]
    // fn test_node_iterator() {
    //     // seq -> node_list -> seq
    //     let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);
    //     let mut unitig_iter = NodeIterator::new(&graph, &SEQ).unwrap();

    //     let mut path = Vec::new();
    //     while let Some(node) = unitig_iter.next().unwrap() {
    //         path.push(node);
    //     }
    //     println!("{:?}", path);
    //     let path_seq = graph.sequence_of_path(path.iter());
    //     let seq = DnaString::from_bytes(SEQ.0);
    //     assert!(path_seq == seq);
    // }

    // #[ignore]
    // #[test]
    // fn print_encode_seq() {
    //     // seq -> (node_list) -> path
    //     let seq = DnaString::from_bytes(SEQ.0);
    //     let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);
    //     let path = ContinuousPath::encode_seq(&graph, &seq);
    //     println!("Start node: {:?}", path.start_node);
    //     for ext in path.extensions.iter() {
    //         println!("{}", ext);
    //     }
    // }

    // #[test]
    // fn test_code_seq() {
    //     // seq -> (node_list) -> path -> (node_list) -> seq
    //     let seq = DnaString::from_bytes(SEQ.0);
    //     let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);

    //     let decoded_seq = ContinuousPath::encode_seq(&graph, &seq).decode_seq(&graph);
    //     assert_eq!(seq, decoded_seq);
    // }

    // #[ignore]
    // #[test]
    // fn print_get_repetitions() {
    //     let test = vec![(0, Dir::Left); 5]
    //         .into_iter()
    //         .chain(vec![(5, Dir::Left); 5])
    //         .chain(vec![(0, Dir::Left); 5])
    //         .chain(vec![(5, Dir::Left); 5])
    //         .collect::<Vec<_>>();
    //     let rep = get_repetitions(&test);
    //     println!("Repetitions:");
    //     for (start, nb_repeats, offset) in rep.iter() {
    //         println!("{} {} {}", start, nb_repeats, offset);
    //     }
    // }
}
