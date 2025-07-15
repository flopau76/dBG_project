//! Transform a dna sequence into a list of nodes in the graph

use crate::{Graph, KmerStorage, Node, PathwayError, Side};

use packed_seq::{PackedSeqVec, Seq, SeqVec};

use std::{cmp::min, ops::Range};

//####################################################################################
//                              Node Iterator                                       //
//####################################################################################

/// Iterator over nodes spelling a given sequence.
pub struct NodeIterator<'a, KS: KmerStorage> {
    graph: &'a Graph<KS>,
    seq: PackedSeqVec,
    next_pos: usize,
    current_node: Option<Node>,
    pub start_offset: Option<usize>,
    pub end_offset: Option<usize>,
}

impl<'a, KS: KmerStorage> NodeIterator<'a, KS> {
    pub fn new(graph: &'a Graph<KS>, seq: PackedSeqVec) -> Result<Self, PathwayError> {
        if seq.len() < graph.k() {
            panic!("sequence is shorter than kmer size");
        }
        let mut node_iter = Self {
            graph,
            seq,
            next_pos: 0,
            current_node: None,
            start_offset: None,
            end_offset: None,
        };
        node_iter.advance()?;
        Ok(node_iter)
    }

    /// get the start position (in nucleotides) of the next node in the iterator
    pub fn position(&self) -> usize {
        self.next_pos
    }

    /// get the current node in the iterator, without advancing
    pub fn peek(&self) -> Option<Node> {
        self.current_node
    }

    /// get the current node in the iterator, before advancing the iterator
    pub fn next(&mut self) -> Result<Option<Node>, PathwayError> {
        let next = self.current_node;
        self.advance()?;
        Ok(next)
    }

    /// advance the iterator
    fn advance(&mut self) -> Result<(), PathwayError> {
        // check if we reached the end of the sequence
        if self.next_pos + self.graph.k() > self.seq.len() {
            self.current_node = None;
            self.end_offset = Some(self.next_pos - 1 + self.graph.k() - self.seq.len());
            return Ok(());
        }

        // search for the next kmer at the start of a node
        let mut kmer = KS::get_kmer(self.graph.k(), self.seq.as_slice(), self.next_pos);
        let node = self.graph.search_kmer(kmer, Side::Left);

        // kmer corresponds to the start of a node
        if node.is_some() {
            self.current_node = node;
            let node = node.unwrap();
            let node_seq = self.graph.node_seq(node);

            // check that the rest of the node sequence matches the input sequence
            let nb_bases = min(node_seq.len(), self.seq.len() - self.next_pos);
            let expected_seq = node_seq.slice(Range {
                start: self.graph.k(),
                end: nb_bases,
            });
            let actual_seq = self.seq.slice(Range {
                start: self.next_pos + self.graph.k(),
                end: self.next_pos + nb_bases,
            });
            if expected_seq != actual_seq {
                return Err(PathwayError::UnitigNotMatching);
            }
            self.next_pos += node_seq.len() - self.graph.k() + 1;
            if self.start_offset.is_none() {
                self.start_offset = Some(0);
            }
        }
        // first kmer does not correspond to the start of a node
        else {
            if self.start_offset.is_some() {
                return Err(PathwayError::KmerNotFound(format!("{:?}", kmer)));
            }
            // the node_iterator has not been initialized yet. We might have an offset in the first node
            // advance in the sequence until we find the last kmer of a node
            let mut node = self.graph.search_kmer(kmer, Side::Right);
            while node.is_none() {
                if self.next_pos + self.graph.k() >= self.seq.len() {
                    panic!("sequence contains neither begining nor end of a node")
                }
                kmer.extend_right(
                    self.graph.k(),
                    self.seq.as_slice().get(self.next_pos + self.graph.k()),
                );
                self.next_pos += 1;
                node = self.graph.search_kmer(kmer, Side::Right);
            }
            // we found the end of a node
            self.current_node = node;
            let node = node.unwrap();
            let node_seq = self.graph.node_seq(node);
            self.start_offset = Some(node_seq.len() - self.next_pos - self.graph.k());

            // check that the beginning of the node sequence matches the input sequence
            let expected_seq = node_seq.slice(Range {
                start: self.start_offset.unwrap(),
                end: self.start_offset.unwrap() + self.next_pos,
            });
            let actual_seq = self.seq.slice(Range {
                start: 0,
                end: self.next_pos,
            });
            if expected_seq != actual_seq {
                return Err(PathwayError::UnitigNotMatching);
            }
            self.next_pos += 1;
        }
        Ok(())
    }
}

impl<KS: KmerStorage> Iterator for NodeIterator<'_, KS> {
    type Item = Node;
    fn next(&mut self) -> Option<Self::Item> {
        self.next().unwrap()
    }
}
