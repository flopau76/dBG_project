/*! # rust_dbg: a crate for the manipulation of de Bruijn graphs
This crates provides methods for the construction and the manipulation of de Bruijn graphs. It builds heavily upon the crate `debruijn` developped by 10x Genomics.
*/

pub mod encoder;
pub mod graph;

// pub use ... for re-exports;
use debruijn::Dir;
use std::io::Write;

use crate::graph::{
    node_iterator::{NodeIterator, PathwayError},
    Graph,
};

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

pub fn print_stats(nb_nn: usize, nb_r: usize, nb_sp: usize, nodes_r: usize, nodes_sp: usize) {
    const NN_COST: usize = 2;
    const SP_COST: usize = 32;
    const R_COST: usize = 24;

    let total_nodes = nb_nn + nodes_r + nodes_sp;
    let total_cost = NN_COST * nb_nn + SP_COST * nb_sp + R_COST * nb_r;

    eprintln!("\n       Method | Number of encoded nodes | Memory cost");
    eprintln!("--------------|-------------------------|-----------------");
    for (name, nodes, cost) in [
        ("Next node", nb_nn, NN_COST * nb_nn),
        ("Shortest path", nodes_sp, SP_COST * nb_sp),
        ("Repetition", nodes_r, R_COST * nb_r),
    ] {
        eprintln!(
            "{:>13} | {:>14}  ({:>4.1}%) | {:>11} bits ({:>4.1}%)",
            name,
            format_int(nodes),
            nodes as f64 / total_nodes as f64 * 100.0,
            format_int(cost),
            cost as f64 / total_cost as f64 * 100.0
        );
    }
    eprintln!(
        "        Total | {:>14}  ( 100%) | {:>11} bits ( 100%)",
        format_int(total_nodes),
        format_int(total_cost),
    );
}
