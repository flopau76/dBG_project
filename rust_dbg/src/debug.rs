use std::collections::VecDeque;

use packed_seq::{PackedSeqVec, Seq, SeqVec};

use crate::graph::shortest_path;
use crate::pathway::{Contig, Extension, Scaffold};
use crate::{Graph, KmerStorage, Node, NodeIterator, PathwayError};
//####################################################################################
//                                  Encoder                                         //
//####################################################################################

/// Parameters for the encoder
pub struct EncoderParams {
    pub min_sp_length: usize,
    pub max_sp_length: usize,
    pub min_nb_repeats: u16,
    pub max_offset: u8,
}

/// Encoder for paths in a de Bruijn graph.
pub struct Encoder<'a, K: KmerStorage> {
    pub params: EncoderParams,
    pub graph: &'a Graph<K>,
}

impl<'a, K: KmerStorage> Encoder<'a, K> {
    /// Encode a list of nodes using the different extensions.
    pub fn encode_path(&self, nodes: &Vec<Node>) -> VecExtension {
        let mut extensions = vec![Extension::TargetNode(
            *nodes.first().expect("The node vector is empty"),
        )];
        println!("{:?}", extensions[0]);

        let mut current_position = 0;
        let mut repetitions = self.get_repetitions(&nodes);
        repetitions.push_back((nodes.len(), 0, 0)); // add a sentinel to avoid adding the last part separately

        while let Some(repetition) = repetitions.pop_front() {
            let target_pos = repetition.0 - 1;
            while current_position < target_pos {
                let (mut target_node, mut length) = shortest_path::get_next_target_node_naive(
                    self.graph,
                    &nodes,
                    current_position,
                    self.params.max_sp_length,
                )
                .unwrap();
                if current_position + length >= target_pos {
                    length = target_pos - current_position;
                    target_node = nodes[target_pos];
                }
                // if the shortest path is long enough, we use it to extend the path
                if length >= self.params.min_sp_length {
                    let ext = Extension::TargetNode(target_node);
                    println!("{:?}:{}", ext, length);
                    extensions.push(ext);
                }
                // otherwise, we encode its nodes directly (2 bits per node)
                else {
                    extensions.extend(
                        nodes[current_position + 1..current_position + 1 + length]
                            .iter()
                            .map(|&node| {
                                let ext = Extension::NextNucleotide(
                                    self.graph.node_seq(node).as_slice().get(self.graph.k() - 1),
                                );
                                println!("{:?}", ext);
                                ext
                            }),
                    );
                }
                current_position += length as usize;
            }
            // add the repetition
            let ext = Extension::Repetition {
                nb_repeats: repetition.1,
                offset: repetition.2,
            };
            println!("{:?}", ext);
            extensions.push(ext);
            current_position += repetition.1 as usize;
        }
        extensions.pop(); // remove the sentinel
        ExtensionVec(extensions)
    }

    // Search for repetitions in a list using Lempel-Ziv
    // Return a Vector of (start, nb_repeats, offset) so that list[start .. start + nb_repeats] = list[start-1 - offset .. start-1- offset + nb_repeats]
    // Maximal offset and minimal number of repeats are defined in the parameters
    fn get_repetitions<D: Eq>(&self, list: &Vec<D>) -> VecDeque<(usize, u16, u8)> {
        let mut repetitions = VecDeque::new();
        let mut current_pos = 1;

        while current_pos < list.len() {
            let mut best_nb_repeats = 0;
            let mut best_offset = 0;

            // find the longest repetition in the sliding window
            for offset in 0..=(self.params.max_offset as usize).min(current_pos - 1) {
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
