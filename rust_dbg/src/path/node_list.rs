//! Transform a dna sequence into a list of nodes in the graph

use crate::PathwayError;
use crate::graph::Graph;

use debruijn::{Dir, Mer, Vmer, Kmer, KmerIter};
use debruijn::dna_string::DnaString;

use super::ElementaryPath;

/// A list of consecutive nodes in the graph, representing a sequence
pub type NodeList = Vec<(usize, Dir)>;

impl<K: Kmer> ElementaryPath<K> for NodeList {
    fn decode_seq(&self, graph: &Graph<K>) -> DnaString {
        graph.sequence_of_path(self.iter())
    }
    fn encode_seq(seq: &DnaString, graph: &Graph<K>) -> Self {
        // TODO: handle start/end offset
        NodeIterator::new(graph, seq).unwrap().collect()
    }
}

/// Iterator over a dna sequence, following the nodes in the graph. This will raise an error if the sequence from kmer_iter is not present in the graph.
pub struct NodeIterator<'a, K: Kmer, D: Vmer> {
    graph: &'a Graph<K>,
    kmer_iter: KmerIter<'a, K, D>,
    next_node: Option<(usize, Dir)>,
    pub start_offset: usize,
    pub end_offset: Option<usize>,
}

impl<'a, K: Kmer, D: Vmer> NodeIterator<'a, K, D> {
    pub fn new(graph: &'a Graph<K>, sequence: &'a D) -> Result<Self, PathwayError<K>> {
        let kmer_iter = sequence.iter_kmers::<K>();
        let mut node_iter = Self{
            graph,
            kmer_iter,
            next_node: None,
            start_offset: 0,
            end_offset: None,
        };

        // initialise the iterator by looking for the first kmer
        let mut kmer = node_iter.kmer_iter.next().expect("kmer iter is empty");
        let node = graph.search_kmer(kmer, Dir::Left);

        // first kmer corresponds to the start of a node
        if node.is_some() {
            node_iter.next_node = node;
            let node = node.unwrap();
            // get the sequence of this node
            let expected_seq = match node.1 {
                Dir::Left => graph.get_node(node.0).sequence(),
                Dir::Right => graph.get_node(node.0).sequence().rc(),
            };
            let expected_seq = expected_seq.slice(K::k(), expected_seq.len());

            // advance the kmer iterator to the end of this node and
            // check that the sequence from kmer_iter coincides with the whole node, not only with its first kmer
            for (idx, expected_base) in expected_seq.iter().enumerate() {
                let kmer = match node_iter.kmer_iter.next(){
                    Some(kmer) => kmer,
                    None => {
                        node_iter.end_offset = Some(expected_seq.len()-idx);
                        break
                    },
                };
                if expected_base != kmer.get(K::k()-1) {
                    return Err(PathwayError::UnitigNotMatching(expected_seq.to_owned(), kmer, 1));
                }
            }
        }
        // first kmer does not correspond to the start of a node
        // => advance in kmer_iter until we find the end of a node
        else {
            let mut node = graph.search_kmer(kmer, Dir::Right);
            let mut skipped_bases = DnaString::new();
            while node.is_none() {
                skipped_bases.push(kmer.get(0));
                kmer = match node_iter.kmer_iter.next() {
                    Some(kmer) => kmer,
                    None => {
                        // edge case: the sequence from kmer_iter is contained within a single node
                        todo!("fix case when sequence is included within a single node");
                    },
                };
                node = graph.search_kmer(kmer, Dir::Right);
            }
            node_iter.next_node = node;
            let node = node.unwrap();
            node_iter.start_offset = graph.get_node(node.0).len() - skipped_bases.len();
            // verify that the skipped bases coincide with the beginning of the node
            let expected_seq = match node.1 {
                Dir::Left => graph.get_node(node.0).sequence(),
                Dir::Right => graph.get_node(node.0).sequence().rc(),
            };
            let expected_seq = expected_seq.slice(node_iter.start_offset, expected_seq.len()-K::k()).to_owned();
            if expected_seq != skipped_bases {
                return Err(PathwayError::UnitigNotMatching(expected_seq, kmer, 0));
            }
        }
        Ok(node_iter)
    }

    /// get the next node in the iterator, without advancing
    pub fn peek(&self) -> Option<(usize, Dir)> {
        self.next_node
    }

    /// get the next node in the iterator, before advancing the iterator
    pub fn next(&mut self) -> Result<Option<(usize, Dir)>, PathwayError<K>> {
        let next = self.next_node;
        self.advance()?;
        Ok(next)
    }

    /// advance the iterator
    pub fn advance(&mut self) -> Result<(), PathwayError<K>> {
        // get the next kmer in the iterator, if any
        let kmer = match self.kmer_iter.next() {
            Some(kmer) => kmer,
            None => {
                self.next_node = None;
                return Ok(())
            },
        };
        // look for the kmer at the beginning of the unitigs
        let node = self.graph.search_kmer(kmer, Dir::Left).ok_or(PathwayError::KmerNotFound(kmer))?;
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
                    self.end_offset = Some(expected_seq.len()-idx);
                    break
                },
            };
            if expected_base != kmer.get(K::k()-1) {
                return Err(PathwayError::UnitigNotMatching(expected_seq.to_owned(), kmer, 1));
            }
        }
        Ok(())
    }
}

impl<K: Kmer, D: Vmer> Iterator for NodeIterator<'_, K, D> {
    type Item = (usize, Dir);
    fn next(&mut self) -> Option<Self::Item> {
        self.next().unwrap()
    }
}