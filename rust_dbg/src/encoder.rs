use crate::graph::shortest_path;
use crate::{Graph, KmerStorage, Node, NodeIterator, PathwayError, Side};

use needletail::Sequence;
use packed_seq::{PackedSeqVec, Seq, SeqVec};
use std::collections::VecDeque;
use std::fmt::Debug;

#[derive(Copy, Clone)]
/// An enum containing different ways to encode a path extension.
pub enum Extension {
    ShortestPath(Node),
    NextNucleotide(u8),
    Repetition((u16, u8)), // nb_repeats, offset (-1)
}

impl Debug for Extension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Extension::ShortestPath(node) => write!(f, "SP({:?})", node),
            Extension::NextNucleotide(base) => write!(f, "NN({})", base),
            Extension::Repetition((nb_repeats, offset)) => {
                write!(f, "R({}-{})", nb_repeats, offset)
            }
        }
    }
}

/// An encoded suite of nodes, described as a start node and a list of extensions
#[derive(Debug, Clone)]
pub struct NodesEncoding {
    start_node: Node,
    extensions: Vec<Extension>,
}

impl NodesEncoding {
    pub fn decode(&self, graph: &Graph<impl KmerStorage>) -> Vec<Node> {
        let mut nodes = vec![self.start_node];
        for extension in &self.extensions {
            match extension {
                Extension::ShortestPath(target_node) => {
                    nodes.extend(
                        shortest_path::get_shortest_path(
                            graph,
                            *nodes.last().unwrap(),
                            *target_node,
                        )
                        .unwrap(),
                    );
                }
                Extension::Repetition((nb_repeats, offset)) => {
                    for _ in 0..*nb_repeats {
                        let prev = nodes[nodes.len() - 1 - *offset as usize];
                        nodes.push(prev);
                    }
                }
                Extension::NextNucleotide(base) => {
                    let last_node = *nodes.last().unwrap();
                    let mut kmer = graph.node_kmer(last_node, Side::Right);
                    kmer.extend_right(graph.k(), *base);
                    let next_node = graph
                        .search_kmer(kmer, Side::Left)
                        .expect("Next node should exist");
                    nodes.push(next_node);
                }
            }
        }
        nodes
    }

    pub fn sequence(&self, graph: &Graph<impl KmerStorage>) -> String {
        let nodes = self.decode(graph);
        unsafe { String::from_utf8_unchecked(graph.path_seq(&nodes)) }
    }
}

/// Contig in a de Bruijn graph, represented as a suite of nodes
#[derive(Debug, Clone)]
pub struct Contig {
    pub nodes_encoding: NodesEncoding,
    pub start_offset: usize,
    pub end_offset: usize,
}

impl Contig {
    /// Returns the sequence of the contig
    pub fn sequence(&self, graph: &Graph<impl KmerStorage>) -> String {
        let seq = self.nodes_encoding.sequence(graph);
        seq[self.start_offset..seq.len() - self.end_offset].to_string() // TODO: more efficient way of slicing ?
    }
}

/// Scaffold in a de Bruijn graph, represented as an alternance of contigs and gaps of a given size.
#[derive(Debug, Clone)]
pub struct Scaffold {
    pub id: String,
    pub contigs: Vec<(usize, Contig)>,
    pub end_gap: usize,
}

impl Scaffold {
    /// Returns the sequence of the scaffold
    pub fn sequence(&self, graph: &Graph<impl KmerStorage>) -> String {
        let mut seq = String::new();
        for (gap_size, contig) in &self.contigs {
            seq.push_str(&"N".repeat(*gap_size));
            seq.push_str(&contig.sequence(graph));
        }
        seq.push_str(&"N".repeat(self.end_gap));
        seq
    }
}

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
        let mut count_n = 0;
        for contig in seq.normalize(false).split(|c| *c == b'N') {
            if contig.is_empty() {
                count_n += 1;
                continue;
            };
            let contig = self.encode_contig(&contig)?;
            contigs.push((count_n, contig));
            count_n = 0;
        }
        Ok(Scaffold {
            id,
            contigs,
            end_gap: count_n,
        })
    }

    /// Encode a contig from a sequence (containing only ACGT bases).
    pub fn encode_contig(&self, seq: &[u8]) -> Result<Contig, PathwayError> {
        let seq = PackedSeqVec::from_ascii(seq);
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
    pub fn encode_path(&self, nodes: &Vec<Node>) -> NodesEncoding {
        let start_node = nodes[0];
        let mut extensions = Vec::new();

        let mut current_position = 0;
        let mut repetitions = self.get_repetitions(&nodes);
        repetitions.push_back((nodes.len(), 0, 0)); // add a sentinel to avoid adding the last part separately

        while let Some(repetition) = repetitions.pop_front() {
            let target_pos = repetition.0 - 1;
            while current_position < target_pos {
                let (mut target_node, mut length) = shortest_path::get_next_target_node(
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
                    let ext = Extension::ShortestPath(target_node);
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
            let ext = Extension::Repetition((repetition.1, repetition.2));
            println!("{:?}", ext);
            extensions.push(ext);
            current_position += repetition.1 as usize;
        }
        extensions.pop(); // remove the sentinel
        NodesEncoding {
            start_node,
            extensions,
        }
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
