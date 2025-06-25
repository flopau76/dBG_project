//! Transform a dna sequence into a list of nodes in the graph

use super::{Graph, Node, Side};
use crate::kmer::KmerStorage;

use packed_seq::{PackedSeqVec, Seq, SeqVec};

use std::{cmp::min, error::Error, ops::Range};

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
//                              Node Iterator                                       //
//####################################################################################

/// Iterator over nodes of a contig present in the graph.
pub struct NodeIterator<'a, KS: KmerStorage> {
    graph: &'a Graph<KS>,
    seq: PackedSeqVec,
    current_kmer: KS,
    current_pos: usize,
    next_node: Option<Node>,
    pub start_offset: usize,
    pub end_offset: Option<usize>,
}

impl<'a, KS: KmerStorage> NodeIterator<'a, KS> {
    pub fn new(graph: &Graph<KS>, seq: PackedSeqVec) -> Result<Self, PathwayError> {
        let mut node_iter = Self {
            graph,
            seq,
            current_kmer: KS::new(),
            current_pos: 0,
            next_node: None,
            start_offset: 0,
            end_offset: None,
        };

        // compute the first kmer from the sequence
        while node_iter.current_pos < graph.k() {
            node_iter
                .current_kmer
                .extend_right(graph.k, seq.as_slice().get(node_iter.current_pos));
            node_iter.current_pos += 1;
        }
        // search for the kmer at the start of a node
        let node = graph.search_kmer(&node_iter.current_kmer, Side::Left);

        // first kmer corresponds to the start of a node
        if node.is_some() {
            node_iter.next_node = node;
            let node = node.unwrap();
            let node_seq = graph.node_seq(node);

            // seq contains only a single node
            if seq.len() <= node_seq.len() {
                node_iter.end_offset = Some(node_seq.len() - seq.len());
            }

            // check that the suite of the sequence matches the node sequence
            let node_end = min(node_seq.len(), seq.len());
            let expected_seq = node_seq.slice(Range {
                start: graph.k(),
                end: node_end,
            });
            let actual_seq = seq.slice(Range {
                start: graph.k(),
                end: node_end,
            });
            if expected_seq != actual_seq {
                return Err(PathwayError::UnitigNotMatching);
            }
        }
        // first kmer does not correspond to the start of a node
        // => advance in the sequence until we find the last kmer of a node
        else {
            let mut node = graph.search_kmer(kmer, Dir::Right);
            let mut skipped_bases = DnaString::new();
            while node.is_none() {
                skipped_bases.push(kmer.get(0));
                kmer = node_iter
                    .kmer_iter
                    .next()
                    .expect("sequence contains neither begining nor end of a node");
                node = graph.search_kmer(kmer, Dir::Right);
            }
            node_iter.next_node = node;
            let node = node.unwrap();
            let node_seq = match node.1 {
                Dir::Left => graph.get_node(node.0).sequence(),
                Dir::Right => graph.get_node(node.0).sequence().rc(),
            };
            node_iter.start_offset = node_seq.len() - skipped_bases.len() - K::k();
            let expected_seq = node_seq
                .slice(node_iter.start_offset, node_seq.len() - K::k())
                .to_owned();
            if expected_seq != skipped_bases {
                println!("skipped bases: {:?}", skipped_bases);
                println!("expected bases: {:?}", expected_seq);
                return Err(PathwayError::UnitigNotMatching(expected_seq, kmer, 0));
                // todo: proper values here
            }
        }
        Ok(node_iter)
    }

    /// get the start position (in nucleotides) of the next node in the iterator
    pub fn position(&self) -> usize {
        self.kmer_iter.pos - K::k()
    }

    /// get the next node in the iterator, without advancing
    pub fn peek(&self) -> Option<Node> {
        self.next_node
    }

    /// get the next node in the iterator, before advancing the iterator
    pub fn next(&mut self) -> Result<Option<Node>, PathwayError> {
        let next = self.next_node;
        self.advance()?;
        Ok(next)
    }

    /// advance the iterator
    fn advance(&mut self) -> Result<(), PathwayError> {
        // get the next kmer in the iterator, if any
        let kmer = match self.kmer_iter.next() {
            Some(kmer) => kmer,
            None => {
                self.next_node = None;
                return Ok(());
            }
        };
        // look for the kmer at the beginning of the unitigs
        let node = self
            .graph
            .search_kmer(kmer, Dir::Left)
            .ok_or(PathwayError::KmerNotFound(kmer))?;
        self.next_node = Some(node);
        // get the sequence of this node
        let expected_seq = match node.1 {
            Dir::Left => self.graph.get_node(node.0).sequence(),
            Dir::Right => self.graph.get_node(node.0).sequence().rc(),
        };
        let expected_seq = expected_seq.slice(K::k(), expected_seq.len());

        // advance the kmer iterator to the end of this node and
        // check that the sequence from kmer_iter coincides with the whole node, not only with its first kmer
        for (idx, expected_base) in expected_seq.iter().enumerate() {
            let kmer = match self.kmer_iter.next() {
                Some(kmer) => kmer,
                None => {
                    self.end_offset = Some(expected_seq.len() - idx);
                    break;
                }
            };
            if expected_base != kmer.get(K::k() - 1) {
                return Err(PathwayError::UnitigNotMatching(
                    expected_seq.to_owned(),
                    kmer,
                    1,
                ));
            }
        }
        Ok(())
    }
}

impl Iterator for NodeIterator<'_> {
    type Item = Node;
    fn next(&mut self) -> Option<Self::Item> {
        self.next().unwrap()
    }
}
