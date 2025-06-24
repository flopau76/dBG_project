/*! # rust_dbg: a crate for the manipulation of de Bruijn graphs
This crates provides methods for the construction and the manipulation of de Bruijn graphs. It builds heavily upon the crate `debruijn` developped by 10x Genomics.
*/

pub mod encoder;
pub mod graph;
pub mod path;

// pub use ... for re-exports;
use debruijn::{dna_only_base_to_bits, Dir, DnaBytes, Kmer};
use needletail::Sequence;
use std::io::Write;

pub use crate::graph::Graph;
use crate::path::node_iterator::{NodeIterator, PathwayError};

//####################################################################################
//                             Basic structures                                     //
//####################################################################################

/// Oriented node in a canonical de Bruijn graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Node(pub usize, pub Dir);

impl Node {
    /// Returns the id of the node
    pub fn id(&self) -> usize {
        self.0
    }

    /// Returns the orientation of the node: Left=canonical, Right=non-canonical
    pub fn dir(&self) -> Dir {
        self.1
    }

    pub fn to_string(&self) -> String {
        format!("({}, {:?})", self.0, self.1)
    }

    pub fn from_string(s: &str) -> Result<Self, Box<dyn std::error::Error>> {
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

        Ok(Self(id, dir))
    }
}

/// Contig in a de Bruijn graph, represented as a suite of nodes
#[derive(Debug, Clone)]
pub struct ContigNodes {
    pub nodes: Vec<Node>,
    pub start_offset: usize,
    pub end_offset: usize,
}

impl ContigNodes {
    /// Creates a ContigNodes from an ascii encoded sequence containing only ACGT bases.
    pub fn from_sequence<'a, K: Kmer>(
        sequence: &'a [u8],
        graph: &graph::Graph<K>,
    ) -> Result<Self, PathwayError<K>> {
        let sequence_bits = DnaBytes(
            sequence
                .iter()
                .map(|b| dna_only_base_to_bits(*b).expect("Contig contains non-ACGT base"))
                .collect(),
        );
        let mut iterator = NodeIterator::new(graph, &sequence_bits)?;
        let mut nodes = Vec::new();
        while let Some(node) = iterator.next()? {
            nodes.push(node);
        }
        Ok(Self {
            nodes,
            start_offset: iterator.start_offset,
            end_offset: iterator.end_offset.expect("End offset should be set"),
        })
    }
}

/// Scaffold in a de Bruijn graph, represented as an alternance of NNNs and contigs
#[derive(Debug, Clone)]
pub struct ScaffoldNodes {
    pub header: String,
    pub contigs: Vec<(usize, ContigNodes)>,
    pub end_gap: usize,
}

impl ScaffoldNodes {
    /// Creates a ScaffoldNodes from an ascii encoded sequence.
    pub fn from_string<'a, K: Kmer>(
        header: String,
        sequence: &'a [u8],
        graph: &graph::Graph<K>,
    ) -> Result<Self, PathwayError<K>> {
        let mut contigs = Vec::new();
        let mut count_n = 0;
        for contig in sequence.normalize(false).split(|c| *c == b'N') {
            if contig.is_empty() {
                count_n += 1;
                continue;
            };
            let contig = ContigNodes::from_sequence(&contig, graph)?;
            contigs.push((count_n, contig));
            count_n = 0;
        }
        Ok(Self {
            header,
            contigs,
            end_gap: count_n,
        })
    }
}

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
