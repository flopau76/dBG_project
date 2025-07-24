use derive_more::with_trait::{Add, AddAssign, Deref, DerefMut, Sum};
use std::fmt::Display;

use std::collections::VecDeque;

use packed_seq::Seq;

use super::{Encoder, Encoding, EncodingStats};

use crate::{graph::shortest_path, Graph, KmerStorage, Node, Side};

//####################################################################################
//                                     Extensions                                   //
//####################################################################################

/// An enum containing different ways to encode a path extension.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Extension {
    TargetNode(Node),
    NextNucleotide(u8),
    Repetition(u16, u8), // nb_repeats, offset
}

impl Display for Extension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Extension::TargetNode(node) => write!(f, "TN({})", node),
            Extension::NextNucleotide(base) => write!(f, "NN({})", base),
            Extension::Repetition(nb_repeats, offset) => {
                write!(f, "R({}-{})", nb_repeats, offset)
            }
        }
    }
}

//####################################################################################
//                                   VecExtensions                                  //
//####################################################################################

#[derive(Default, Debug, Deref, DerefMut)]
pub struct VecExtensions(pub Vec<Extension>);

impl Encoding for VecExtensions {
    type Stats = ExtensionStats;
    fn decode(&self, graph: &Graph<impl KmerStorage>) -> Vec<Node> {
        let mut nodes = Vec::new();

        match self.first() {
            Some(Extension::TargetNode(start_node)) => {
                nodes.push(*start_node);
            }
            None => panic!("The extension vector is empty"),
            _ => panic!("Invalid encoding: the first extension should be a TargetNode"),
        }

        for extension in &self[1..] {
            match extension {
                Extension::TargetNode(target_node) => {
                    nodes.extend(
                        crate::graph::shortest_path::get_shortest_path(
                            graph,
                            *nodes.last().unwrap(),
                            *target_node,
                        )
                        .unwrap(),
                    );
                }
                Extension::Repetition(nb_repeats, offset) => {
                    for _ in 0..*nb_repeats {
                        let prev = nodes[nodes.len() - 1 - *offset as usize];
                        nodes.push(prev);
                    }
                }
                Extension::NextNucleotide(base) => {
                    let last_node = *nodes.last().unwrap();
                    let mut kmer = graph.node_kmer(last_node, Side::Right);
                    kmer.extend_right(graph.k(), *base);
                    let next_node = graph.search_kmer(kmer, Side::Left).expect(&format!(
                        "Next expected seq {} after node {:?}:{} not found in graph",
                        kmer.print(graph.k()),
                        last_node,
                        unsafe {
                            String::from_utf8_unchecked(
                                graph.node_seq(last_node).as_slice().unpack(),
                            )
                        }
                    ));
                    nodes.push(next_node);
                }
            }
        }
        nodes
    }

    fn get_encoding_stats(&self, graph: &Graph<impl KmerStorage>) -> Self::Stats {
        let mut nucleo_count = 0;
        let mut targets_count = 0;
        let mut reps_count = 0;
        let mut reps_nodes = 0;

        for ext in self.iter() {
            match ext {
                Extension::TargetNode(_) => {
                    targets_count += 1;
                }
                Extension::NextNucleotide(_) => nucleo_count += 1,
                Extension::Repetition(nb_repeats, _) => {
                    reps_count += 1;
                    reps_nodes += *nb_repeats as usize;
                }
            }
        }

        let total_nodes = self.decode(graph).len();

        ExtensionStats {
            nucleo_count,
            targets_count,
            reps_count,
            nucleo_nodes: nucleo_count,
            targets_nodes: total_nodes - nucleo_count - reps_nodes,
            reps_nodes,
        }
    }
}

//####################################################################################
//                                   Greedy Encoder                                 //
//####################################################################################

/// Greedy encoder, which transforms a list of nodes into some extensions based on specific cutoffs.
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

impl Encoder for GreedyEncoder {
    type Encoding = VecExtensions;
    fn encode_path(&self, nodes: Vec<Node>, graph: &Graph<impl KmerStorage>) -> VecExtensions {
        let mut extensions = vec![Extension::TargetNode(
            *nodes.first().expect("The node vector is empty"),
        )];

        let mut current_position = 0;
        let mut repetitions = self.get_repetitions(&nodes);
        repetitions.push_back((nodes.len(), (0, 0))); // add a sentinel to avoid adding the last part separately

        while let Some((rep_pos, rep)) = repetitions.pop_front() {
            let target_pos = rep_pos - 1;
            while current_position < target_pos {
                let (mut target_node, mut length) = shortest_path::get_next_target_node(
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
                    println!("{}", ext);
                }
                // otherwise, we encode its nodes directly (2 bits per node)
                else {
                    extensions.extend(
                        nodes[current_position + 1..current_position + 1 + length]
                            .iter()
                            .map(|&node| {
                                let ext = Extension::NextNucleotide(
                                    graph.node_seq(node).as_slice().get(graph.k() - 1),
                                );
                                println!("{}", ext);
                                ext
                            }),
                    );
                }
                current_position += length as usize;
            }
            // add the repetition
            let ext = Extension::Repetition(rep.0, rep.1);
            extensions.push(ext);
            println!("{}", ext);
            current_position += rep.0 as usize;
        }
        extensions.pop(); // remove the sentinel
        VecExtensions(extensions)
    }
}

//####################################################################################
//                                ExtensionStats                                    //
//####################################################################################

#[derive(Default, Debug, Add, AddAssign, Sum)]
pub struct ExtensionStats {
    nucleo_count: usize,
    targets_count: usize,
    reps_count: usize,
    nucleo_nodes: usize,
    targets_nodes: usize,
    reps_nodes: usize,
}
impl Display for ExtensionStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let total_count = self.nucleo_count + self.targets_count + self.reps_count;
        let total_nodes = self.nucleo_nodes + self.targets_nodes + self.reps_nodes;

        writeln!(
            f,
            "       Method | Number of extensions | Number of encoded nodes "
        )?;
        writeln!(
            f,
            "--------------|-------------------------|-----------------"
        )?;
        for (name, count, nodes) in [
            ("Next node", self.nucleo_count, self.nucleo_nodes),
            ("Target node", self.targets_count, self.targets_nodes),
            ("Repetition", self.reps_count, self.reps_nodes),
        ] {
            writeln!(
                f,
                "{:>13} | {:>14}  ({:>4.1}%) | {:>11} bits ({:>4.1}%)",
                name,
                crate::format_int(count),
                count as f64 / total_count as f64 * 100.0,
                crate::format_int(nodes),
                nodes as f64 / total_nodes as f64 * 100.0,
            )?;
        }
        writeln!(
            f,
            "        Total | {:>14}          | {:>11} bits",
            crate::format_int(total_count),
            crate::format_int(total_nodes),
        )?;
        Ok(())
    }
}

impl EncodingStats for ExtensionStats {}
