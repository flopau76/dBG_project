/*! # rust_dbg: a crate for the manipulation of de Bruijn graphs
This crates provides methods for the construction and the manipulation of de Bruijn graphs. It builds heavily upon the crate `debruijn` developped by 10x Genomics.
*/

pub mod fasta_reader;
pub mod graph;

pub mod path;

// pub use ... for re-exports

//####################################################################################
//                              Custom errors                                       //
//####################################################################################
use std::error::Error;
use debruijn::dna_string::DnaString;
use debruijn::Kmer;

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
