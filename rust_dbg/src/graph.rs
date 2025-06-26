//! To create graphs from fasta/unitig files

use crate::kmer::KmerStorage;

use std::{collections::HashSet, path::Path};

use boomphf::Mphf;
use needletail::parse_fastx_file;
use packed_seq::{PackedSeq, PackedSeqVec, Seq, SeqVec};
use smallvec::SmallVec;

pub mod node_iterator;

/// A 2-bit packed owned set of multiple DNA sequences.
#[derive(Debug, Default)]
pub struct SequenceSet {
    sequences: PackedSeqVec,
    start_positions: Vec<usize>,
    length: Vec<u32>,
}

impl SequenceSet {
    /// Add a packed sequence to the set.
    pub fn push_seq(&mut self, seq: PackedSeq) -> () {
        let range = self.sequences.push_seq(seq);
        self.start_positions.push(range.start);
        self.length.push(seq.len() as u32);
        println!("Adding seq: {} of length {}", seq.as_u64(), seq.len());
        println!("range: {:?}", range.clone());
        println!("{:?}", self.sequences);
        assert_eq!(self.sequences.slice(range).unpack(), seq.unpack());
    }

    /// Add an ascii sequence to the set.
    /// Undefined behavior if the sequence contains characters other than ACGTacgt.
    pub fn push_ascii(&mut self, seq: &[u8]) -> () {
        let range = self.sequences.push_ascii(seq);
        self.start_positions.push(range.start);
        self.length.push(seq.len() as u32);
    }

    /// Get a PackedSeq containing the 'i'th sequence in the set.
    pub fn get(&self, i: usize) -> PackedSeq {
        self.sequences
            .slice(self.start_positions[i]..self.start_positions[i] + self.length[i] as usize)
    }

    /// Get the number of sequences in the set.
    pub fn len(&self) -> usize {
        self.start_positions.len()
    }
}

/// Oriented node in a canonical de Bruijn graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Node {
    pub id: usize,
    pub is_rc: bool,
}

impl Node {
    pub fn new(id: usize, is_rc: bool) -> Self {
        Self { id, is_rc }
    }
}

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

/// A base graph, containing only unitigs. No edges nor hash-functions.
pub struct BaseGraph {
    pub k: usize,
    pub stranded: bool,
    pub sequences: SequenceSet,
}

impl BaseGraph {
    /// Get the number of nodes in the graph.
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// Create a graph from a fasta file containing unitigs (as returned by ggcat for example).
    pub fn from_unitigs_serial(path: &Path, k: usize, stranded: bool) -> Self {
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

    /// Create a (non-compacted) graph from a sequence. (For debugging mainly)
    pub fn from_seq_serial<'a>(seq: &[u8], k: usize, stranded: bool) -> Self {
        if k > 32 {
            panic!("Unsuported kmer size: {}. Maximum is 32.", k);
        }
        // get unique kmers and add them to the sequence set
        let seq = PackedSeqVec::from_ascii(seq);
        let mut sequences = SequenceSet::default();
        let mut kmer_set = HashSet::new();
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

    pub fn finish<KS: KmerStorage>(self) -> Graph<KS> {
        // compute the edges
        // TODO

        // compute the kmers
        let mut left_kmers = Vec::with_capacity(self.sequences.len());
        let mut right_kmers = Vec::with_capacity(self.sequences.len());
        for i in 0..self.sequences.len() {
            let seq = self.sequences.get(i);
            let left_kmer = KS::get_kmer(self.k, seq, 0);
            let right_kmer = KS::get_kmer(self.k, seq, seq.len() - self.k);
            left_kmers.push(left_kmer);
            right_kmers.push(right_kmer);
        }
        let left_kmers_mphf = Mphf::new(2.0, &left_kmers);
        let right_kmers_mphf = Mphf::new(2.0, &right_kmers);

        // create the graph
        Graph {
            base: self,
            edges: vec![],
            left_kmers: left_kmers_mphf,
            right_kmers: right_kmers_mphf,
        }
    }
}

pub struct Graph<KS: KmerStorage> {
    base: BaseGraph,
    edges: Vec<usize>,
    left_kmers: Mphf<KS>,
    right_kmers: Mphf<KS>,
}

/// For Graph manipulation
impl<KS: KmerStorage> Graph<KS> {
    /// Return the kmer size of the graph.
    pub fn k(&self) -> usize {
        self.base.k
    }
    /// Return wether the graph is stranded or not.
    pub fn stranded(&self) -> bool {
        self.base.stranded
    }

    /// Get the sequence of a node.
    pub fn node_seq(&self, node: Node) -> PackedSeqVec {
        if node.is_rc {
            self.base.sequences.get(node.id).to_revcomp()
        } else {
            self.base.sequences.get(node.id).to_vec()
        }
    }

    /// Get the kmer at the given `side` of a `node`.
    pub fn node_kmer(&self, node: Node, side: Side) -> KS {
        match node.is_rc {
            true => {
                let seq_vec = self.base.sequences.get(node.id).to_revcomp();
                let seq = seq_vec.as_slice();
                let pos = side.choose(0, seq.len() - self.k());
                KS::get_kmer(self.k(), seq, pos)
            }
            false => {
                let seq = self.base.sequences.get(node.id);
                let pos = side.choose(0, seq.len() - self.k());
                KS::get_kmer(self.k(), seq, pos)
            }
        }
    }

    /// Get the neighbors of a `node` on a given `side`.
    pub fn node_neigh(&self, node: Node, side: Side) -> SmallVec<[Node; 4]> {
        let mut neighbors = SmallVec::new();
        let kmer = self.node_kmer(node, side);
        // TODO: use self.edges to reduce calls to search_kmer ?
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
        let mphf = side.choose(&self.left_kmers, &self.right_kmers);
        if let Some(id) = mphf.try_hash(&kmer) {
            let node = Node::new(id as usize, false);
            if self.node_kmer(node, side) == kmer {
                return Some(node);
            }
        }
        if !self.stranded() {
            let rc_kmer = kmer.rc(self.k());
            let mphf = side.choose(&self.right_kmers, &self.left_kmers);
            if let Some(id) = mphf.try_hash(&rc_kmer) {
                let node = Node::new(id as usize, true);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seq_set() {
        let mut seq_set = SequenceSet::default();
        let seq1 = PackedSeqVec::from_ascii(b"ACGTACGT");
        let seq2 = PackedSeqVec::from_ascii(b"TTTT");
        seq_set.push_seq(seq1.as_slice());
        seq_set.push_seq(seq2.as_slice());

        let seq1b = seq_set.get(0).unpack();
        let seq2b = seq_set.get(1).unpack();
        println!("Seq1: {seq1b:?}");
        println!("Seq2: {seq2b:?}");
        println!("{seq_set:?}");
        assert_eq!(seq1b, seq1.as_slice().unpack());
        assert_eq!(seq2b, seq2.as_slice().unpack());
    }

    #[ignore]
    #[test]
    fn test_from_seq() {
        let seq = b"AACCCCGGTTTTT";
        let graph = BaseGraph::from_seq_serial(seq, 3, true);
        println!("{:?}", graph.sequences.start_positions);
        for i in 0..graph.sequences.len() {
            let seq = graph.sequences.get(i);
            let kmer = seq.as_u64();
            println!("Kmer {i} of len {}: {kmer:b}", seq.len());
        }
        assert_eq!(graph.len(), 8); // 8 unique kmers of size 3
        println!("{:?}", graph.sequences);
    }
}
