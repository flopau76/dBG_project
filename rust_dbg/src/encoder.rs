use crate::{Dir, Graph, Node};
use debruijn::Kmer;
use std::collections::VecDeque;

mod shortest_path;

#[derive(Copy, Clone)]
/// An enum containing all possible ways to encode a path extension.
enum Extension {
    ShortestPath(Node),
    NextNode(Node),
    Repetition((u16, u8)), // nb_repeats, offset (-1)
}

pub struct EncoderParams {
    pub min_sp_length: usize,
    pub max_sp_length: usize,
    pub min_nb_repeats: u16,
    pub max_offset: u8,
}

impl Default for EncoderParams {
    fn default() -> Self {
        Self {
            min_sp_length: 1,
            max_sp_length: 100,
            min_nb_repeats: 10,
            max_offset: 255,
        }
    }
}

/// Same as ContigNodes but with a list of extensions instead of nodes.
pub struct ContigExtensions {
    pub start_node: Node,
    pub extensions: Vec<Extension>,
    pub start_offset: usize,
    pub end_offset: usize,
}

/// Same as ScaffoldNodes but with a list of extensions instead of contigs.
pub struct ScaffoldExtensions {
    pub header: String,
    pub contigs: Vec<(usize, ContigExtensions)>,
    pub end_gap: usize,
}

pub struct Encoder<'a, K: Kmer> {
    params: EncoderParams,
    graph: &'a Graph<K>,
}

impl<'a, K: Kmer> Encoder<'a, K> {
    fn encode_contig(&self, nodes: &Vec<Node>) -> ContigExtensions {
        let start_node = nodes[0];
        let mut extensions = Vec::new();

        let mut current_position = 0;
        let mut repetitions = self.get_repetitions(&nodes);
        repetitions.push_back((nodes.len(), 0, 0)); // add a sentinel to avoid adding the last part separately

        while let Some(repetition) = repetitions.pop_front() {
            let target_pos = repetition.0 - 1;
            while current_position < target_pos {
                let (mut shortest_path, mut length) =
                    shortest_path::get_next_target_node(self.graph, &nodes, current_position)
                        .unwrap();
                if current_position + length >= target_pos {
                    length = target_pos - current_position;
                    shortest_path = nodes[target_pos];
                }
                println!("SP:{}:{:?}", length, shortest_path);
                // if the shortest path is long enough, we use it to extend the path
                if length >= self.params.min_sp_length {
                    extensions.push(Extension::ShortestPath(shortest_path));
                }
                // otherwise, we encode its nodes directly (2 bits per node)
                else {
                    extensions.extend(
                        nodes[current_position + 1..current_position + 1 + length]
                            .iter()
                            .map(|&node| Extension::NextNode(node)),
                    );
                }
                current_position += length as usize;
            }
            // add the repetition
            extensions.push(Extension::Repetition((repetition.1, repetition.2)));
            current_position += repetition.1 as usize;
            println!("R:{}x{}", repetition.1, repetition.2);
        }
        extensions.pop(); // remove the sentinel
        (start_node, extensions)
    }

    pub fn decode(&self, extensions: Vec<Extension>) -> Vec<Node> {
        let nodes = vec![extensions[0].0];
        todo!()
                match self {
            Extension::ShortestPath(target_node, _) => path.extend(
                shortest_path::get_shortest_path(graph, *path.last().unwrap(), *target_node)
                    .unwrap()
                    .iter(),
            ),
            Extension::NextNode(nn) => path.push(*nn),
            Extension::Repetition((nb_repeats, offset)) => {
                for _ in 0..*nb_repeats {
                    let prev = path[path.len() - 1 - *offset as usize];
                    path.push(prev);
                }
            }
        }
    }

    // Return a Vector of (start_position, nb_repeats, offset),
    // so that path[start_position .. start_position + nb_repeats] = path[start_position -1 - offset .. start_position + nb_repeats - 1- offset]
    fn get_repetitions(&self, path: &Vec<Node>) -> VecDeque<(usize, u16, u8)> {
        let mut repetitions = VecDeque::new();
        let mut current_pos = 1;

        while current_pos < path.len() {
            let mut best_nb_repeats = 0;
            let mut best_offset = 0;

            // find the longest repetition in the sliding window
            for offset in 0..=(self.params.max_offset as usize).min(current_pos - 1) {
                let mut nb_repeats = 0;
                while path[current_pos - 1 - offset + nb_repeats] == path[current_pos + nb_repeats]
                {
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
            if best_nb_repeats >= self.params.min_nb_repeats as usize {
                repetitions.push_back((current_pos, best_nb_repeats as u16, best_offset as u8));
                current_pos += best_nb_repeats;
            } else {
                current_pos += 1;
            }
        }
        repetitions
    }
}
