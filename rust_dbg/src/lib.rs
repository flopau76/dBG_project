/*! # rust_dbg: a crate for the manipulation of de Bruijn graphs
This crates provides methods for the construction and the manipulation of de Bruijn graphs. It builds heavily upon the crate `debruijn` developped by 10x Genomics.
*/

pub mod graph;
pub mod path;

// pub use ... for re-exports
// pub use debruijn::Dir;

//####################################################################################
//                              Custom errors                                       //
//####################################################################################
use debruijn::dna_string::DnaString;
use debruijn::Kmer;
use std::error::Error;
use std::io::Write;

/// Custom error type for pathway search operations
#[derive(Debug)]
pub enum PathwayError<K> {
    NoPathExists,
    KmerNotFound(K),
    NodeNotMatching(DnaString, K, usize),
}
impl<K: Kmer> Error for PathwayError<K> {}
impl<K: Kmer> std::fmt::Display for PathwayError<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PathwayError::NoPathExists => write!(f, "No path found between the given k-mers"),
            PathwayError::KmerNotFound(kmer) => write!(f, "Kmer not found: {:?}", kmer),
            PathwayError::NodeNotMatching(seq, kmer, i) => write!(
                f,
                "Expected kmer {:?} at position {} in unitig {}",
                kmer, i, seq
            ),
        }
    }
}

//####################################################################################
//                             Utility functions                                    //
//####################################################################################

/// Parses a string representing a node in the format "(id, direction)"
pub fn parse_node(s: &str) -> Result<(usize, debruijn::Dir), Box<dyn std::error::Error>> {
    let (id, dir) = s
        .strip_prefix('(')
        .and_then(|s| s.strip_suffix(')'))
        .and_then(|s| s.split_once(','))
        .ok_or(format!("Invalid node format: {}", s))?;

    let id = id.parse::<usize>()?;

    let dir = match dir.trim() {
        "Left" => debruijn::Dir::Left,
        "Right" => debruijn::Dir::Right,
        other => return Err(format!("Invalid direction: {}", other).into()),
    };

    Ok((id, dir))
}

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
