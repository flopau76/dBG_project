/*! # rust_dbg: a crate for the manipulation of de Bruijn graphs
This crates provides methods for the construction and the manipulation of de Bruijn graphs. It builds heavily upon the crate `debruijn` developped by 10x Genomics.
*/

pub mod fasta_reader;
pub mod graph;
pub mod path;

pub mod graph_old;
pub mod path_old;