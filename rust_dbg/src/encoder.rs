use std::collections::VecDeque;
use std::fmt::Debug;
use std::ops::Add;

use bincode::{Decode, Encode};
use needletail::Sequence;
use packed_seq::{PackedSeqVec, Seq, SeqVec};

use crate::graph::shortest_path;
use crate::{Graph, KmerStorage, Node, NodeIterator, PathwayError, Side};

#[allow(non_snake_case)]
mod ExtensionVec_Encode;
#[allow(non_snake_case)]
mod Node_Encode;

//####################################################################################
//                        Extension  &  ExtensionVec                                //
//####################################################################################

/// An enum containing different ways to encode a path extension.
#[derive(Encode, Decode, Copy, Clone, PartialEq, Eq)]
pub enum Extension {
    TargetNode(Node),
    NextNucleotide(u8),
    Repetition { nb_repeats: u16, offset: u8 },
}

impl Debug for Extension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Extension::TargetNode(node) => write!(f, "TN({:?})", node),
            Extension::NextNucleotide(base) => write!(f, "NN({})", base),
            Extension::Repetition { nb_repeats, offset } => {
                write!(f, "R({}-{})", nb_repeats, offset)
            }
        }
    }
}

/// A vector of extensions, describing a path in a de Bruijn graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtensionVec(pub Vec<Extension>);

impl ExtensionVec {
    /// Converts the vector of extensions into a suite of nodes in the graph.
    pub fn decode(&self, graph: &Graph<impl KmerStorage>) -> Vec<Node> {
        let mut nodes = Vec::new();

        match self.0.first() {
            Some(Extension::TargetNode(start_node)) => {
                nodes.push(*start_node);
            }
            None => panic!("The extension vector is empty"),
            _ => panic!("Invalid encoding: the first extension should be a TargetNode"),
        }

        for extension in &self.0[1..] {
            match extension {
                Extension::TargetNode(target_node) => {
                    nodes.extend(
                        shortest_path::get_shortest_path(
                            graph,
                            *nodes.last().unwrap(),
                            *target_node,
                        )
                        .unwrap(),
                    );
                }
                Extension::Repetition { nb_repeats, offset } => {
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

    /// Returns the sequence of the path represented by the extensions.
    pub fn sequence(&self, graph: &Graph<impl KmerStorage>) -> String {
        let nodes = self.decode(graph);
        unsafe { String::from_utf8_unchecked(graph.path_seq(&nodes)) }
    }
}

//####################################################################################
//                             Contig  &  Scaffold                                  //
//####################################################################################

/// Contig in a de Bruijn graph, represented as a suite of nodes
#[derive(Encode, Decode, Debug, Clone)]
pub struct Contig {
    pub nodes_encoding: ExtensionVec,
    pub start_offset: usize,
    pub end_offset: usize,
}

impl Contig {
    /// Returns the sequence of the contig
    pub fn sequence(&self, graph: &Graph<impl KmerStorage>) -> String {
        let seq = self.nodes_encoding.sequence(graph);
        seq[self.start_offset..seq.len() - self.end_offset].to_string()
    }
}

/// Scaffold in a de Bruijn graph, represented as an alternance of contigs and gaps of a given size.
#[derive(Encode, Decode, Debug, Clone)]
pub struct Scaffold {
    pub id: String,
    pub contigs: Vec<Contig>,
    pub gaps: Vec<usize>,
}

impl Scaffold {
    /// Returns the sequence of the scaffold
    pub fn sequence(&self, graph: &Graph<impl KmerStorage>) -> String {
        assert_eq!(
            self.contigs.len() + 1,
            self.gaps.len(),
            "Expected one mor gap than contigs, but got {} contigs and {} gaps",
            self.contigs.len(),
            self.gaps.len()
        );
        let mut seq = "N".repeat(self.gaps[0]);
        for (contig, gap) in std::iter::zip(&self.contigs, &self.gaps[1..]) {
            seq.push_str(&contig.sequence(graph));
            seq.push_str(&"N".repeat(*gap));
        }
        seq
    }

    /// Prints statistics about the scaffold
    pub fn print_stats(&self, graph: &Graph<impl KmerStorage>) {
        println!("Scaffold ID: {}", self.id);
        println!("Number of contigs: {}", self.contigs.len());
        let sum = self
            .contigs
            .iter()
            .map(|c| c.nodes_encoding.0.len())
            .sum::<usize>();
        println!("Number of extensions in contigs: {}", sum);
        println!("Total length: {}", self.sequence(graph).len());
    }
}

//####################################################################################
//                                  Stats                                           //
//####################################################################################

struct Stats {
    nb_repeats: usize,
    nodes_repeats: usize,
    nb_target: usize,
    nodes_target: usize,
    nb_nucleo: usize,
    gap_count: usize,
    gap_size: usize,
}

impl Add for Stats {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Stats {
            nb_repeats: self.nb_repeats + other.nb_repeats,
            nodes_repeats: self.nodes_repeats + other.nodes_repeats,
            nb_target: self.nb_target + other.nb_target,
            nodes_target: self.nodes_target + other.nodes_target,
            nb_nucleo: self.nb_nucleo + other.nb_nucleo,
            gap_count: self.gap_count + other.gap_count,
            gap_size: self.gap_size + other.gap_size,
        }
    }
}

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
    /// Encode a scaffold from a sequence (containing non ACGT bases).
    pub fn encode_record(&self, id: String, seq: &[u8]) -> Result<Scaffold, PathwayError> {
        let mut contigs = Vec::new();
        let mut gaps = Vec::new();
        let mut gap_size = 0;
        for contig in seq.normalize(false).split(|c| *c == b'N') {
            if contig.is_empty() {
                gap_size += 1;
                continue;
            };
            gaps.push(gap_size);
            let contig = self.encode_contig(&contig)?;
            contigs.push(contig);
            gap_size = 0;
        }
        gaps.push(gap_size);

        Ok(Scaffold { id, contigs, gaps })
    }

    /// Encode a contig from a sequence (containing only ACGT bases).
    pub fn encode_contig(&self, seq_in: &[u8]) -> Result<Contig, PathwayError> {
        let seq = PackedSeqVec::from_ascii(seq_in);
        let mut iterator = NodeIterator::new(self.graph, seq)?;
        let mut nodes = Vec::new();
        while let Some(node) = iterator.next()? {
            nodes.push(node);
        }
        let nodes_encoding = self.encode_path(&nodes);
        Ok(Contig {
            nodes_encoding,
            start_offset: iterator.start_offset.expect("Start offset should be set"),
            end_offset: iterator.end_offset.expect("End offset should be set"),
        })
    }

    /// Encode a list of nodes using the different extensions.
    pub fn encode_path(&self, nodes: &Vec<Node>) -> ExtensionVec {
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
