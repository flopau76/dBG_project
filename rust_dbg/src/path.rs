//! Defines how to encode a path in a debruijn Graph
 
use debruijn::Kmer;
use debruijn::dna_string::DnaString;

use crate::graph::Graph;

pub mod node_list;
pub mod shortest_path;

/// An elementary path represents a continuous sequence of nucleotides
pub trait ElementaryPath<K: Kmer> {
    /// Return the sequence of nucleotides encoded in the path
    fn decode_seq(&self, graph: &Graph<K>) -> DnaString;

    /// Encode a sequence of nucleotides into a given graph
    fn encode_seq(seq: &DnaString, graph: &Graph<K>) -> Self;
}

impl<K: Kmer> ElementaryPath<K> for DnaString {
    fn decode_seq(&self, _graph: &Graph<K>) -> DnaString {
        self.to_owned()
    }
    fn encode_seq(seq: &DnaString, _graph: &Graph<K>) -> Self {
        seq.to_owned()
    }
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use crate::graph::Graph;
    use node_list::NodeList;
    use shortest_path::ShortestPathList;

    use debruijn::kmer::Kmer3;
    use debruijn::{DnaSlice, dna_string::DnaString};

    const STRANDED : bool = true;
    const SEQ: DnaSlice = DnaSlice(&[2,2,2,1,1,1,1,2,2,2,0,0,0,0,0,1]);    // gggccccgggaaaaac
    const _SHORTEST: DnaSlice = DnaSlice(&[2,2,2,0,0,1]); // gggaac

    #[test]
    fn node_list() {
        let graph = Graph::<Kmer3>::from_seq_serial(SEQ, STRANDED);
        let seq = DnaString::from_bytes(SEQ.0);

        let node_list = NodeList::encode_seq(&seq, &graph);
        let reconstructed_seq = node_list.decode_seq(&graph);

        assert!(reconstructed_seq == seq);
    }

    #[test]
    fn shortest_path() {
        let graph = Graph::<Kmer3>::from_seq_serial(SEQ, STRANDED);
        let seq = DnaString::from_bytes(SEQ.0);

        let shortest_path = ShortestPathList::encode_seq(&seq, &graph);
        let reconstructed_seq = shortest_path.decode_seq(&graph);

        assert!(reconstructed_seq == seq);
    }
}