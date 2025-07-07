/*! # rust_dbg: a crate for the manipulation of de Bruijn graphs
This crates provides methods for the construction and the manipulation of de Bruijn graphs. It builds heavily upon the crate `debruijn` developped by 10x Genomics.
*/

// pub mod encoder;
pub mod encoder;
pub mod graph;
pub mod kmer;

// pub use ... for re-exports;
use std::io::Write;

pub use graph::{BaseGraph, Graph, Node, NodeIterator, PathwayError, Side};
pub use kmer::KmerStorage;

//####################################################################################
//                             Utility functions                                    //
//####################################################################################

/// Removes the previous line and print a new progress bar
pub fn print_progress_bar(current: usize, total: usize) {
    let bar_width = 40;
    let progress = current as f32 / total as f32;
    let filled = (progress * bar_width as f32).round() as usize;
    let empty = bar_width - filled;
    let bar = format!(
        "[{}{}] {}/{} ({:.0}%)",
        "#".repeat(filled),
        " ".repeat(empty),
        current,
        total,
        progress * 100.0
    );
    eprint!("\r{}", bar); // carriage return to overwrite the current line
    std::io::stderr().flush().unwrap();
}

/// Format a long integer with commas
pub fn format_int(n: usize) -> String {
    let s = n.to_string();
    let mut result = String::new();
    let mut chars = s.chars().rev().peekable();

    while let Some(c) = chars.next() {
        result.push(c);
        if chars.peek().is_some() && result.len() % 4 == 3 {
            result.push(',');
        }
    }

    result.chars().rev().collect()
}
