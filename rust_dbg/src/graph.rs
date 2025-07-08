//! To create graphs from fasta/unitig files

use smallvec::SmallVec;
use std::fmt::Debug;
use std::ops::Range;
use std::{collections::HashSet, error::Error, path::Path};

use boomphf::hashmap::NoKeyBoomHashMap;
use epserde::prelude::*;
use needletail::parse_fastx_file;
use packed_seq::{PackedSeq, PackedSeqVec, Seq, SeqVec};

use crate::KmerStorage;

mod node_iterator;
pub use node_iterator::NodeIterator;
pub mod shortest_path;

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

/// 2-bit packed owned set of multiple DNA sequences.
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

//####################################################################################
//                                     Node                                         //
//####################################################################################

/// Oriented node in a canonical de Bruijn graph
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Node {
    pub id: usize,
    pub is_rc: bool,
}

impl Node {
    pub fn new(id: usize, is_rc: bool) -> Self {
        Self { id, is_rc }
    }
}

impl Debug for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_rc {
            write!(f, "Node(-{})", self.id)
        } else {
            write!(f, "Node(+{})", self.id)
        }
    }
}

/// Oriented node in a canonical de Bruijn graph. `NodeSeq` may or may not own the data, depending if it is canonical or reverse complement.
#[derive(Debug)]
pub enum NodeSeq<'a> {
    Owned(PackedSeqVec),
    Borrowed(PackedSeq<'a>),
}

impl<'a> NodeSeq<'a> {
    pub fn as_slice(&self) -> PackedSeq {
        match self {
            NodeSeq::Owned(seq) => seq.as_slice(),
            NodeSeq::Borrowed(seq) => *seq,
        }
    }
    pub fn slice(&self, range: Range<usize>) -> PackedSeq {
        self.as_slice().slice(range)
    }
    pub fn len(&self) -> usize {
        self.as_slice().len()
    }
}

//####################################################################################
//                                     Side                                         //
//####################################################################################

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

/// Base graph, containing only unitigs, without indexing.
/// This is only an itermediate step to build the final graph.
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
    pub fn node_seq(&self, node: Node) -> NodeSeq {
        if node.is_rc {
            NodeSeq::Owned(self.sequences.get(node.id).to_revcomp())
        } else {
            NodeSeq::Borrowed(self.sequences.get(node.id))
        }
    }

    /// Get the kmer at the given `side` of a `node`.
    pub fn node_kmer<KS: KmerStorage>(&self, node: Node, side: Side) -> KS {
        let seq = self.node_seq(node);
        let pos = side.choose(0, seq.len() - self.k());
        KS::get_kmer(self.k(), seq.as_slice(), pos)
    }

    /// Create a graph from a fasta file containing unitigs (as returned by ggcat for example).
    pub fn from_unitig_file(path: &Path, k: usize, stranded: bool) -> Self {
        // iterate over the unitigs and add them to the sequence set
        let mut sequences = SequenceSet::default();
        let mut fasta_reader = parse_fastx_file(path)
            .expect(format!("Could not read unitigs from {}", path.display()).as_str());
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
    pub fn from_seq<'a>(seq: &PackedSeqVec, k: usize, stranded: bool) -> Self {
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
                kmer_u64 = kmer_u64.min(rc_kmer_u64); // TODO: this is not really the canonical kmers, as ordering is ACTG
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
        Self::from_seq(&seq, k, stranded)
    }

    pub fn finish<KS: KmerStorage>(self) -> Graph<KS> {
        assert!(
            KS::capacity() >= self.k(),
            "KmerStorage {} cant hold k = {} nucleotides. Maximum is {}.",
            std::any::type_name::<KS>(),
            self.k(),
            KS::capacity()
        );

        // compute the kmers
        let mut left_kmers = Vec::with_capacity(self.len());
        let mut right_kmers = Vec::with_capacity(self.len());
        for i in 0..self.len() {
            let node = Node::new(i, false);
            let left_kmer: KS = self.node_kmer(node, Side::Left);
            let right_kmer: KS = self.node_kmer(node, Side::Right);
            left_kmers.push(left_kmer);
            right_kmers.push(right_kmer);
        }
        let left_ids: Vec<usize> = (0..self.len()).collect();
        let right_ids: Vec<usize> = (0..self.len()).collect();

        // check that the kmers are unique
        let mut left_kmers_set = HashSet::new();
        let mut right_kmers_set = HashSet::new();
        for kmer in &left_kmers {
            if !left_kmers_set.insert(kmer.clone()) {
                panic!(
                    "Left kmer {} is not unique in the graph.",
                    kmer.print(self.k())
                );
            }
        }
        for kmer in &right_kmers {
            if !right_kmers_set.insert(kmer.clone()) {
                panic!(
                    "Right kmer {} is not unique in the graph.",
                    kmer.print(self.k())
                );
            }
        }

        // create the mphf for left and right kmers
        let left_kmers_mphf = NoKeyBoomHashMap::new(left_kmers, left_ids);
        let right_kmers_mphf = NoKeyBoomHashMap::new(right_kmers, right_ids);

        // create the graph
        Graph {
            base: self,
            left_kmers: left_kmers_mphf,
            right_kmers: right_kmers_mphf,
        }
    }

    /// Print the base graph, with all its nodes and edges
    pub fn print(&self) {
        println!("BaseGraph: k = {}, stranded = {}", self.k, self.stranded);
        for i in 0..self.len() {
            let node = Node::new(i, false);
            let seq = self.node_seq(node);
            let seq_str = unsafe { String::from_utf8_unchecked(seq.as_slice().unpack()) };
            println!("Node {}: {}", i, seq_str); // TODO: add edges
        }
    }

    /// Print some summary stats
    pub fn print_stats(&self) {
        println!(
            "BaseGraph: k = {}, stranded = {}, nodes = {}",
            self.k,
            self.stranded,
            self.len(),
        );
    }

    /// Write the base graph to a binary file.
    pub fn save_to_binary(&self, path: &Path) -> Result<(), Box<dyn Error>> {
        self.store(path)?;
        Ok(())
    }

    /// Load the graph from a binary file.
    pub fn load_from_binary(path: &Path) -> Result<Self, Box<dyn Error>> {
        let b = std::fs::read(&path)?;
        let base = BaseGraph::deserialize_eps(&b)?;
        Ok(base)
    }
}

//####################################################################################
//                                  Graph                                           //
//####################################################################################

/// De Bruijn graph, containing unitigs and edges.
#[derive(Debug)]
pub struct Graph<KS: KmerStorage> {
    base: BaseGraph,
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
    pub fn node_seq(&self, node: Node) -> NodeSeq {
        self.base.node_seq(node)
    }
    /// Get the kmer of a `node` at a given `side`.
    #[inline]
    pub fn node_kmer(&self, node: Node, side: Side) -> KS {
        self.base.node_kmer(node, side)
    }
    /// Print the graph.
    #[inline]
    pub fn print(&self) {
        self.base.print();
    }
    /// Print some summary stats.
    #[inline]
    pub fn print_stats(&self) {
        self.base.print_stats();
        println!("Underlying kmer storage: {}", std::any::type_name::<KS>());
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

    /// Get the ascii sequence of a `path`.
    pub fn path_seq(&self, path: &Vec<Node>) -> Vec<u8> {
        let mut path_seq = self
            .node_seq(*path.first().expect("path is empty"))
            .as_slice()
            .unpack();
        for node in path.iter().skip(1) {
            let node_seq = self.node_seq(*node);
            let node_seq = node_seq.slice(self.k() - 1..node_seq.len()).unpack();
            path_seq.extend(node_seq);
        }
        path_seq
    }

    /// Search a `kmer` at one `side` of the nodes.
    pub fn search_kmer(&self, kmer: KS, side: Side) -> Option<Node> {
        let map = side.choose(&self.left_kmers, &self.right_kmers);
        if let Some(id) = map.get(&kmer) {
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
                if self.node_kmer(node, side) == kmer {
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

/// Graph construction
impl<KS: KmerStorage> Graph<KS> {
    /// Create a graph from a fasta file containing unitigs (as returned by ggcat for example).
    pub fn from_unitig_file(path: &Path, k: usize, stranded: bool) -> Self {
        let base = BaseGraph::from_unitig_file(path, k, stranded);
        base.finish()
    }

    /// Create a (non-compacted) graph from a packed sequence. (For debugging mainly)
    #[deprecated = "warning: the obtained graph will not be compacted"]
    pub fn from_seq<'a>(seq: &PackedSeqVec, k: usize, stranded: bool) -> Self {
        let base = BaseGraph::from_seq(seq, k, stranded);
        base.finish()
    }

    /// Create a (non-compacted) graph from an ascii sequence. (For debugging mainly)
    #[deprecated = "warning: the obtained graph will not be compacted"]
    pub fn from_seq_ascii<'a>(seq: &[u8], k: usize, stranded: bool) -> Self {
        let base = BaseGraph::from_seq_ascii(seq, k, stranded);
        base.finish()
    }
}
