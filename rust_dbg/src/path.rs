//! Defines how to encode a path in a debruijn Graph. Pri
 
use debruijn::{Dir, Kmer};

use crate::graph::Graph;

pub mod node_iterator;
pub mod shortest_path;

/// Some information about how to extend a path in a given graph
pub trait Extension<K: Kmer> {
    /// Get the next extension, given the full path and the current position in the path
    fn next(graph: &Graph<K>, path: &Vec<(usize, Dir)>, position: usize) -> Self;
    /// Extend the current path given the extension
    fn extend_path(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>);
}

/// The easiest way to extend the current path is to give the following node
impl<K: Kmer> Extension<K> for (usize, Dir) {
    fn next(_graph: &Graph<K>, path: &Vec<(usize, Dir)>, position: usize) -> Self {
        *path.get(position+1).expect("We are allready at the end of the path")
    }
    fn extend_path(&self, _graph: &Graph<K>, path: &mut Vec<(usize, Dir)>) {
        path.push(*self);
    }
}

// #[cfg(test)]
// mod unit_test {
//     use super::*;

//     use crate::graph::Graph;
//     use node_list::NodeList;

//     use debruijn::kmer::Kmer3;
//     use debruijn::{DnaSlice, dna_string::DnaString};

//     const STRANDED : bool = true;
//     const SEQ: DnaSlice = DnaSlice(&[2,2,2,1,1,1,1,2,2,2,0,0,0,0,0,1]);    // gggccccgggaaaaac
//     const _SHORTEST: DnaSlice = DnaSlice(&[2,2,2,0,0,1]); // gggaac

// }