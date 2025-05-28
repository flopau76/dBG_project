//! Defines how to encode a path in a debruijn Graph. Pri

use debruijn::{dna_string::DnaString, Dir, Kmer, Vmer};
use needletail::Sequence;

use self::node_iterator::NodeIterator;
use crate::format_int;
use crate::graph::Graph;

use std::collections::VecDeque;
use std::fmt::Display;

use std::io::{Read, Write};

use bincode::{deserialize_from, serialize_into};
use serde::{Deserialize, Serialize};

pub mod node_iterator;
mod shortest_path;

// for shortest path
pub const MIN_PATH_LENGTH: usize = 17; // path encoded on 32 bits
pub const MAX_PATH_LENGTH: usize = 60;
// for repetitions
pub const MIN_NB_REPEATS: u16 = 13; // repetition encoded on 24 bits
pub const MAX_OFFSET: u8 = 255;

#[derive(Serialize, Deserialize, Copy, Clone)]
/// An enum containing all possible ways to encode a path extension.
enum MyExtension {
    ShortestPath((usize, Dir), usize), // target_node, length. Note: not necessary to encode length, but usefull for stats
    NextNode((usize, Dir)),
    Repetition((u16, u8)), // nb_repeats, offset (-1)
}

impl MyExtension {
    fn to_string(&self) -> String {
        match self {
            MyExtension::ShortestPath(target_node, length) => {
                format!("SP:{}:{:?}", length, target_node)
            }
            MyExtension::NextNode(next_node) => format!("NN{:?}", next_node),
            MyExtension::Repetition((nb_repeats, offset)) => format!("R:{}x{}", nb_repeats, offset),
        }
    }
    fn extend_path<K: Kmer>(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>) {
        match self {
            MyExtension::ShortestPath(target_node, _) => path.extend(
                shortest_path::get_shortest_path(graph, *path.last().unwrap(), *target_node)
                    .unwrap()
                    .iter(),
            ),
            MyExtension::NextNode(nn) => path.push(*nn),
            MyExtension::Repetition((nb_repeats, offset)) => {
                for _ in 0..*nb_repeats {
                    let prev = path[path.len() - 1 - *offset as usize];
                    path.push(prev);
                }
            }
        }
    }
}
impl Display for MyExtension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

/// A struct that represents a continuous path in a de Bruijn graph using a start node and a list of extensions.
#[derive(Serialize, Deserialize)]
struct ContinuousPath {
    pub start_node: (usize, Dir),
    pub extensions: Vec<MyExtension>,
}

// Return a Vector of (start_position, nb_repeats, offset),
// so that path[start_position .. start_position + nb_repeats] = path[start_position -1 - offset .. start_position + nb_repeats - 1- offset]
fn get_repetitions(path: &Vec<(usize, Dir)>) -> VecDeque<(usize, u16, u8)> {
    let mut repetitions = VecDeque::new();
    let mut current_pos = 1;

    while current_pos < path.len() {
        let mut best_nb_repeats = 0;
        let mut best_offset = 0;

        // find the longest repetition in the sliding window
        for offset in 0..=(MAX_OFFSET as usize).min(current_pos - 1) {
            let mut nb_repeats = 0;
            while path[current_pos - 1 - offset + nb_repeats] == path[current_pos + nb_repeats] {
                nb_repeats += 1;
                if current_pos + nb_repeats >= path.len() {
                    break;
                }
            }
            if nb_repeats > best_nb_repeats {
                best_nb_repeats = nb_repeats;
                best_offset = offset;
            }
        }
        // if the repetition is long enough, add it to the list
        if best_nb_repeats >= MIN_NB_REPEATS as usize {
            repetitions.push_back((current_pos, best_nb_repeats as u16, best_offset as u8));
            current_pos += best_nb_repeats;
        } else {
            current_pos += 1;
        }
    }
    repetitions
}

impl ContinuousPath {
    /// Encode the list of nodes from a sequence into a MixedPath.
    pub fn encode_seq<K: Kmer, V: Vmer>(graph: &Graph<K>, seq: &V) -> Self {
        let nodes = NodeIterator::new(graph, seq).unwrap().collect::<Vec<_>>();
        let start_node = nodes[0];
        let mut extensions = Vec::new();

        let mut current_position = 0;
        let mut repetitions = get_repetitions(&nodes);
        repetitions.push_back((nodes.len(), 0, 0)); // add a sentinel to avoid adding the last part separately

        while let Some(repetition) = repetitions.pop_front() {
            let target_pos = repetition.0 - 1;
            while current_position < target_pos {
                let (mut shortest_path, mut length) =
                    shortest_path::get_next_target_node(graph, &nodes, current_position)
                        .unwrap()
                        .unwrap();
                if current_position + length >= target_pos {
                    length = target_pos - current_position;
                    shortest_path = nodes[target_pos];
                }
                // if the shortest path is long enough, we use it to extend the path
                if length >= MIN_PATH_LENGTH {
                    extensions.push(MyExtension::ShortestPath(shortest_path, length));
                }
                // otherwise, we encode its nodes directly (2 bits per node)
                else {
                    extensions.extend(
                        nodes[current_position + 1..current_position + 1 + length]
                            .iter()
                            .map(|&node| MyExtension::NextNode(node)),
                    );
                }
                current_position += length as usize;
            }
            // add the repetition
            extensions.push(MyExtension::Repetition((repetition.1, repetition.2)));
            current_position += repetition.1 as usize;
        }
        extensions.pop(); // remove the sentinel
        ContinuousPath {
            start_node,
            extensions,
        }
    }

    /// Decode the path into a DnaString.
    pub fn decode_seq<K: Kmer>(&self, graph: &Graph<K>) -> DnaString {
        let mut path = vec![self.start_node];
        for ext in self.extensions.iter() {
            ext.extend_path(graph, &mut path);
        }
        graph.sequence_of_path(path.iter())
    }

    /// Transform the path into a string representation (to save to text file)
    fn to_string(&self) -> String {
        let mut result = String::new();
        result.push_str(&format!("Start node: {:?}", self.start_node));
        for ext in &self.extensions {
            result.push_str(&format!("\n{}", ext.to_string()));
        }
        result
    }
}

impl Display for ContinuousPath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

/// A struct encoding a sequence in a de Bruijn graph.
#[derive(Serialize, Deserialize)]
pub struct DiscontinuousPath {
    id: String,
    sequences: Vec<ContinuousPath>,
    masks: Vec<usize>,
}

impl DiscontinuousPath {
    pub fn header(&self) -> &str {
        &self.id
    }
    /// Encode a sequence into a Path in a graph.
    pub fn encode_record<'a, K: Kmer, S: Sequence<'a>>(
        graph: &'a Graph<K>,
        record: &'a S,
        id: String,
    ) -> Self {
        let mut sequences = Vec::new();
        let mut masks = Vec::new();

        let mut count_n = 0;
        for seq in record.normalize(false).split(|c| *c == b'N') {
            if seq.is_empty() {
                count_n += 1;
                continue;
            }
            masks.push(count_n);
            let seq = DnaString::from_acgt_bytes(seq);
            let path = ContinuousPath::encode_seq(graph, &seq);
            sequences.push(path);
            count_n = 0;
        }
        masks.push(count_n);

        Self {
            id,
            sequences,
            masks,
        }
    }

    /// Decode the encoded path into a string representation.
    pub fn decode_record<K: Kmer>(&self, graph: &Graph<K>) -> String {
        let mut record = String::new();
        assert_eq!(self.masks.len(), self.sequences.len() + 1);
        for i in 0..self.sequences.len() {
            record.push_str("N".repeat(self.masks[i]).as_str());
            let seq = self.sequences[i].decode_seq(graph);
            record.push_str(&seq.to_string());
        }
        record.push_str("N".repeat(*self.masks.last().unwrap()).as_str());
        record
    }

    /// Append a single path to a binary file
    pub fn append_to_binary<W: Write>(
        &self,
        file_writer: &mut W,
    ) -> Result<(), Box<dyn std::error::Error>> {
        serialize_into(file_writer, self)?;
        Ok(())
    }

    /// Load a single path from the given file reader
    pub fn load_from_binary<R: Read>(
        file_reader: &mut R,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let result = deserialize_from(file_reader)?;
        Ok(result)
    }

    /// Print some stats about the path.
    pub fn print_stats(&self) {
        let mut nb_nn = 0;
        let mut nb_r = 0;
        let mut nb_sp = 0;
        let mut nodes_r = 0;
        let mut nodes_sp = 0;

        for seq in self.sequences.iter() {
            for ext in seq.extensions.iter() {
                match ext {
                    MyExtension::NextNode(_) => {
                        nb_nn += 1;
                    }
                    MyExtension::Repetition((nb_repeats, _offset)) => {
                        nb_r += 1;
                        nodes_r += *nb_repeats as usize;
                    }
                    MyExtension::ShortestPath(_target_node, length) => {
                        nb_sp += 1;
                        nodes_sp += *length;
                    }
                }
            }
        }
        // print the stats
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
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use crate::graph::Graph;

    use debruijn::kmer::Kmer3;
    use debruijn::{dna_string::DnaString, DnaSlice};

    const STRANDED: bool = true;
    const SEQ: DnaSlice = DnaSlice(&[2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 1]); // gggccccgggaaaaac

    #[test]
    fn test_node_iterator_offset() {
        let seq = DnaString::from_dna_string("aaccggtt");
        let graph = Graph::<Kmer3>::from_seq_serial(&seq, true);

        let seq_start = DnaString::from_dna_string("aacc");
        let unitig_iter_start = NodeIterator::new(&graph, &seq_start).unwrap();
        assert_eq!(unitig_iter_start.start_offset, 0);
        assert_eq!(unitig_iter_start.end_offset, Some(4));

        let seq_end = DnaString::from_dna_string("ggtt");
        let unitig_iter_end = NodeIterator::new(&graph, &seq_end).unwrap();
        assert_eq!(unitig_iter_end.start_offset, 4);
        assert_eq!(unitig_iter_end.end_offset, None);

        // let seq_middle = DnaString::from_dna_string("ccgg");
        // let unitig_iter_middle = NodeIterator::new(&graph, &seq_middle).unwrap();
    }

    #[test]
    fn test_node_iterator() {
        // seq -> node_list -> seq
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);
        let mut unitig_iter = NodeIterator::new(&graph, &SEQ).unwrap();

        let mut path = Vec::new();
        while let Some(node) = unitig_iter.next().unwrap() {
            path.push(node);
        }
        let path_seq = graph.sequence_of_path(path.iter());
        let seq = DnaString::from_bytes(SEQ.0);
        assert!(path_seq == seq);
    }

    #[ignore]
    #[test]
    fn print_encode_seq() {
        // seq -> (node_list) -> path
        let seq = DnaString::from_bytes(SEQ.0);
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);
        let path = ContinuousPath::encode_seq(&graph, &seq);
        for ext in path.extensions.iter() {
            println!("{}", ext);
        }
    }

    #[test]
    fn test_code_seq() {
        // seq -> (node_list) -> path -> (node_list) -> seq
        let seq = DnaString::from_bytes(SEQ.0);
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);

        let decoded_seq = ContinuousPath::encode_seq(&graph, &seq).decode_seq(&graph);
        assert_eq!(seq, decoded_seq);
    }

    #[ignore]
    #[test]
    fn print_get_repetitions() {
        let test = vec![(0, Dir::Left); 5]
            .into_iter()
            .chain(vec![(5, Dir::Left); 5])
            .chain(vec![(0, Dir::Left); 5])
            .chain(vec![(5, Dir::Left); 5])
            .collect::<Vec<_>>();
        let rep = get_repetitions(&test);
        println!("Repetitions:");
        for (start, nb_repeats, offset) in rep.iter() {
            println!("{} {} {}", start, nb_repeats, offset);
        }
    }
}
