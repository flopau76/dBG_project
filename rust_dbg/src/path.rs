//! Defines how to encode a path in a debruijn Graph. Pri

use debruijn::{Dir, Kmer, dna_string::DnaString};

use self::node_iterator::NodeIterator;
use crate::graph::Graph;
use crate::parse_node;

use std::collections::VecDeque;
use std::fmt::Debug;

pub mod node_iterator;
mod shortest_path;

// for shortest path
pub const MAX_PATH_LENGTH: usize = 60;
pub const MIN_PATH_LENGTH: usize = 17; // path encoded on 32 bits
// for repetitions
pub const MAX_OFFSET: u8 = 255;
pub const MIN_NB_REPEATS: u16 = 13; // repetition encoded on 24 bits

#[derive(Debug, Copy, Clone)]
/// An enum containing all possible ways to encode a path extension.
enum MyExtension {
    ShortestPath((usize, Dir)),
    NextNode((usize, Dir)),
    Repetition((u16, u8)), // nb_repeats, offset (-1)
}

impl MyExtension {
    fn to_string(&self) -> String {
        match self {
            MyExtension::ShortestPath(target_node) => format!("SP{:?}", target_node),
            MyExtension::NextNode(next_node) => format!("NN{:?}", next_node),
            MyExtension::Repetition((size, position)) => format!("R:{}x{}", size, position),
        }
    }
    fn from_string(s: &str) -> Self {
        if s.starts_with("SP") {
            let target_node = parse_node(s[2..].trim()).unwrap();
            MyExtension::ShortestPath(target_node)
        } else if s.starts_with("NN") {
            let target_node = parse_node(s[2..].trim()).unwrap();
            MyExtension::NextNode(target_node)
        } else if s.starts_with("R:") {
            let (size, count) = s[2..].split_once("x").unwrap();
            let nb_repeats = size.parse::<u16>().unwrap();
            let offset = count.parse::<u8>().unwrap();
            MyExtension::Repetition((nb_repeats, offset))
        } else {
            panic!("Unknown extension type: {}", s);
        }
    }
    fn extend_path<K: Kmer>(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>) {
        match self {
            MyExtension::ShortestPath(target_node) => path.extend(
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

/// A struct that represents a path in a de Bruijn graph using a start node and a list of extensions.
pub struct MixedPath<'a, K: Kmer> {
    graph: &'a Graph<K>,
    start_node: (usize, Dir),
    extensions: Vec<MyExtension>,
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

impl<'a, K: Kmer> MixedPath<'a, K> {
    /// Encode the list of nodes from a sequence into a MixedPath.
    pub fn encode_seq(graph: &'a Graph<K>, seq: &DnaString) -> Self {
        let nodes = NodeIterator::new(graph, seq).unwrap().collect::<Vec<_>>();
        let start_node = nodes[0];
        let mut extensions = Vec::new();

        let mut current_position = 0;
        let mut repetitions = get_repetitions(&nodes);
        repetitions.push_back((nodes.len(), 0, 0)); // add a sentinel to avoid adding the last part separately

        while let Some(repetition) = repetitions.pop_front() {
            let target_pos = repetition.0 - 1;
            while current_position < target_pos {
                let (shortest_path, length) =
                    shortest_path::get_next_target_node(graph, &nodes, current_position)
                        .unwrap()
                        .unwrap();
                let length = length.min(target_pos - current_position);
                // if the shortest path is long enough, we use it to extend the path
                if length >= MIN_PATH_LENGTH {
                    extensions.push(MyExtension::ShortestPath(shortest_path));
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
        MixedPath {
            graph,
            start_node,
            extensions,
        }
    }

    /// Decode the path into a DnaString.
    pub fn decode_seq(&self) -> DnaString {
        let mut path = vec![self.start_node];
        for ext in self.extensions.iter() {
            ext.extend_path(self.graph, &mut path);
        }
        self.graph.sequence_of_path(path.iter())
    }

    /// Transform the path into a string representation (to save to text file)
    pub fn to_string(&self) -> String {
        let mut result = String::new();

        let git_hash = option_env!("GIT_COMMIT_HASH").unwrap_or("unknown");
        result.push_str(&format!("Code version: {}\n", git_hash));
        result.push_str("Constants:\n");
        result.push_str(&format!("\tMIN_PATH_LENGTH: {}\n", MIN_PATH_LENGTH));
        result.push_str(&format!("\tMAX_PATH_LENGTH: {}\n", MAX_PATH_LENGTH));
        result.push_str(&format!("\tMIN_NB_REPEATS: {}\n", MIN_NB_REPEATS));
        result.push_str(&format!("\tMAX_OFFSET: {}\n", MAX_OFFSET));

        result.push_str(&format!("Start node: {:?}", self.start_node));
        for ext in &self.extensions {
            result.push_str(&format!("\n{}", ext.to_string()));
        }
        result
    }

    /// Create a path from its string representation (to load from text file)
    pub fn from_string(s: &str, graph: &'a Graph<K>) -> Result<Self, Box<dyn std::error::Error>> {
        let mut lines = s.lines().skip(6);
        let start_node_line = lines.next().ok_or("Missing start node line")?;
        let start_node = parse_node(&start_node_line[12..])?;

        let mut extensions = Vec::new();
        for line in lines {
            let ext = MyExtension::from_string(line);
            extensions.push(ext);
        }

        Ok(MixedPath {
            graph,
            start_node,
            extensions,
        })
    }
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use crate::graph::Graph;

    use debruijn::kmer::Kmer3;
    use debruijn::{DnaSlice, dna_string::DnaString};

    const STRANDED: bool = true;
    const SEQ: DnaSlice = DnaSlice(&[2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 1]); // gggccccgggaaaaac
    const _SHORTEST: DnaSlice = DnaSlice(&[2, 2, 2, 0, 0, 1]); // gggaac

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
        let path = MixedPath::encode_seq(&graph, &seq);
        for ext in path.extensions.iter() {
            println!("{:?}", ext);
        }
    }

    #[test]
    fn test_code_seq() {
        // seq -> (node_list) -> path -> (node_list) -> seq
        let seq = DnaString::from_bytes(SEQ.0);
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);

        let decoded_seq = MixedPath::encode_seq(&graph, &seq).decode_seq();
        assert_eq!(seq, decoded_seq);
    }

    #[test]
    fn test_path_to_string() {
        // seq -> (node_list) -> path -> string -> path -> (node_list) -> seq
        let seq = DnaString::from_bytes(SEQ.0);
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);

        let path_str = MixedPath::encode_seq(&graph, &seq).to_string();
        let decoded_seq = MixedPath::from_string(&path_str, &graph)
            .unwrap()
            .decode_seq();
        assert_eq!(seq, decoded_seq);
    }

    #[ignore]
    #[test]
    fn print_get_repetitions() {
        let test = vec![
            (0, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (0, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
            (5, Dir::Left),
        ];
        let rep = get_repetitions(&test);
        println!("Repetitions:");
        for (start, nb_repeats, offset) in rep.iter() {
            println!("{} {} {}", start, nb_repeats, offset);
        }
    }
}
