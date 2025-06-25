//! To create graphs from fasta/unitig files

use crate::kmer::KmerStorage;

use std::ops::Range;
use std::path::Path;

use boomphf::hashmap::NoKeyBoomHashMap;
use needletail::parse_fastx_file;
use packed_seq::{PackedSeq, PackedSeqVec, Seq, SeqVec};

mod node_iterator;

/// A 2-bit packed owned set of multiple DNA sequences.
pub struct SequenceSet {
    sequences: PackedSeqVec,
    start_positions: Vec<usize>,
}

impl SequenceSet {
    /// Create an empty `SequenceSet`
    pub fn new() -> Self {
        Self {
            sequences: PackedSeqVec::default(),
            start_positions: vec![0],
        }
    }

    /// Add an ascii sequence to the set.
    /// Undefined behavior if the sequence contains characters other than ACGTacgt.
    pub fn push_ascii(&mut self, seq: &[u8]) -> () {
        self.sequences.push_ascii(seq);
        self.start_positions.push(self.sequences.len());
    }

    /// Get a PackedSeq containing the 'i'th sequence in the set.
    pub fn get(&self, i: usize) -> PackedSeq {
        self.sequences.slice(Range {
            start: self.start_positions[i],
            end: self.start_positions[i + 1],
        })
    }

    /// Get the number of sequences in the set.
    pub fn len(&self) -> usize {
        self.start_positions.len()
    }
}

/// Oriented node in a canonical de Bruijn graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Node {
    /// The id of the node
    pub id: usize,
    /// The orientation of the node (only needed in canonical graphs)
    pub is_rc: bool,
}

pub enum Side {
    Left,
    Right,
}

pub struct Graph<KS: KmerStorage> {
    k: usize,
    stranded: bool,
    sequences: SequenceSet,
    edges: Vec<usize>,
    hash_map: NoKeyBoomHashMap<KS, usize>,
}

impl<KS: KmerStorage> Graph<KS> {
    pub fn k(&self) -> usize {
        self.k
    }

    /// Create a graph from a fasta file containing unitigs (as returned by ggcat for example).
    pub fn from_unitigs_serial(path: &Path, k: usize, stranded: bool) -> Self {
        let mut sequences = SequenceSet::new();

        // iterate over the unitigs and add them to the sequence set
        let mut fasta_reader = parse_fastx_file(path).unwrap();
        while let Some(record) = fasta_reader.next() {
            let record = record.unwrap();
            let seq = record.seq();
            sequences.push_ascii(&seq);
        }

        // compute the edges

        todo!();
    }

    /// Search a `kmer` at one `side` of the nodes.
    pub fn search_kmer(&self, kmer: &KS, side: Side) -> Option<Node> {
        todo!();
    }

    /// Search a `kmer` by iterating over all the nodes.
    /// Returns the node and the offset of the kmer in the node (from the given `side`).
    pub fn search_kmer_offset(&self, kmer: &KS, side: Side) -> Option<(Node, usize)> {
        todo!();
    }

    /// Get the sequence of a node.
    pub fn node_seq(&self, node: Node) -> PackedSeqVec {
        if node.is_rc {
            self.sequences.get(node.id).to_revcomp()
        } else {
            self.sequences.get(node.id).to_vec()
        }
    }
}
