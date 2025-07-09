use std::collections::VecDeque;

use packed_seq::Seq;

use crate::embeddings::{Extension, Pathway, VecExtensions, VecNodes};
use crate::graph::shortest_path;
use crate::{Graph, KmerStorage};

//####################################################################################
//                                 Encoder                                          //
//####################################################################################

pub trait Encoder<P: Pathway> {
    // transform a list of nodes into any type of encoding
    fn encode_path(&self, nodes: VecNodes, graph: &Graph<impl KmerStorage>) -> P;
}

/// Greedy encoder, transforming a list of nodes into some extensions based on specific cutoffs.
pub struct GreedyEncoder {
    pub min_sp_length: usize,
    pub max_sp_length: usize,
    pub min_nb_repeats: u16,
    pub max_offset: u8,
}

impl Default for GreedyEncoder {
    fn default() -> Self {
        GreedyEncoder {
            min_sp_length: 10,
            max_sp_length: 100,
            min_nb_repeats: 10,
            max_offset: u8::MAX,
        }
    }
}

impl GreedyEncoder {
    // Search for repetitions in a list using Lempel-Ziv
    // Return a Vector of (start, nb_repeats, offset) so that list[start .. start + nb_repeats] = list[start-1 - offset .. start-1- offset + nb_repeats]
    // Maximal offset and minimal number of repeats are defined in the parameters
    fn get_repetitions<D: Eq>(&self, list: &Vec<D>) -> VecDeque<(usize, (u16, u8))> {
        let mut repetitions = VecDeque::new();
        let mut current_pos = 1;

        while current_pos < list.len() {
            let mut best_nb_repeats = 0;
            let mut best_offset = 0;

            // find the longest repetition in the sliding window
            for offset in 0..=(self.max_offset as usize).min(current_pos - 1) {
                let mut nb_repeats = 0;
                while list[current_pos - 1 - offset + nb_repeats] == list[current_pos + nb_repeats]
                {
                    nb_repeats += 1;
                    if current_pos + nb_repeats >= list.len() {
                        break;
                    }
                }
                if nb_repeats > best_nb_repeats {
                    best_nb_repeats = nb_repeats;
                    best_offset = offset;
                }
            }
            // if the repetition is long enough, add it to the list
            if best_nb_repeats >= self.min_nb_repeats as usize {
                repetitions.push_back((current_pos, (best_nb_repeats as u16, best_offset as u8)));
                current_pos += best_nb_repeats;
            } else {
                current_pos += 1;
            }
        }
        repetitions
    }
}

impl Encoder<VecExtensions> for GreedyEncoder {
    fn encode_path(&self, nodes: VecNodes, graph: &Graph<impl KmerStorage>) -> VecExtensions {
        let mut extensions = vec![Extension::TargetNode(
            *nodes.first().expect("The node vector is empty"),
        )];

        let mut current_position = 0;
        let mut repetitions = self.get_repetitions(&nodes);
        repetitions.push_back((nodes.len(), (0, 0))); // add a sentinel to avoid adding the last part separately

        while let Some((rep_pos, rep)) = repetitions.pop_front() {
            let target_pos = rep_pos - 1;
            while current_position < target_pos {
                let (mut target_node, mut length) = shortest_path::get_next_target_node_naive(
                    graph,
                    &nodes,
                    current_position,
                    self.max_sp_length,
                )
                .unwrap();
                if current_position + length >= target_pos {
                    length = target_pos - current_position;
                    target_node = nodes[target_pos];
                }
                // if the shortest path is long enough, we use it to extend the path
                if length >= self.min_sp_length {
                    let ext = Extension::TargetNode(target_node);
                    extensions.push(ext);
                }
                // otherwise, we encode its nodes directly (2 bits per node)
                else {
                    extensions.extend(
                        nodes[current_position + 1..current_position + 1 + length]
                            .iter()
                            .map(|&node| {
                                Extension::NextNucleotide(
                                    graph.node_seq(node).as_slice().get(graph.k() - 1),
                                )
                            }),
                    );
                }
                current_position += length as usize;
            }
            // add the repetition
            let ext = Extension::Repetition(rep.0, rep.1);
            extensions.push(ext);
            current_position += rep.0 as usize;
        }
        extensions.pop(); // remove the sentinel
        VecExtensions(extensions)
    }
}
