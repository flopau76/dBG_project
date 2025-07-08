use std::collections::VecDeque;
use std::error::Error;
use std::fmt::Debug;
use std::iter::Sum;
use std::ops::Add;
use std::path::Path;

use bincode::{Decode, Encode};
use needletail::Sequence;
use packed_seq::{PackedSeqVec, Seq, SeqVec};

use crate::format_int;
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

    /// Get statistics about the node encoding
    fn get_stats(&self, graph: &Graph<impl KmerStorage>) -> StatsEncoding {
        let mut nucleo_count = 0;
        let mut targets_count = 0;
        let mut reps_count = 0;
        let mut reps_nodes = 0;

        for ext in &self.0 {
            match ext {
                Extension::TargetNode(_) => {
                    targets_count += 1;
                }
                Extension::NextNucleotide(_) => nucleo_count += 1,
                Extension::Repetition {
                    nb_repeats,
                    offset: _,
                } => {
                    reps_count += 1;
                    reps_nodes += *nb_repeats as usize;
                }
            }
        }

        let total_nodes = self.decode(graph).len();

        StatsEncoding {
            nucleo_count,
            targets_count,
            reps_count,
            nucleo_nodes: nucleo_count,
            targets_nodes: total_nodes - nucleo_count - reps_nodes,
            reps_nodes,
        }
    }
}

struct StatsEncoding {
    nucleo_count: usize,
    targets_count: usize,
    reps_count: usize,
    nucleo_nodes: usize,
    targets_nodes: usize,
    reps_nodes: usize,
}

impl Add for StatsEncoding {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        StatsEncoding {
            nucleo_count: self.nucleo_count + other.nucleo_count,
            targets_count: self.targets_count + other.targets_count,
            reps_count: self.reps_count + other.reps_count,
            nucleo_nodes: self.nucleo_nodes + other.nucleo_nodes,
            targets_nodes: self.targets_nodes + other.targets_nodes,
            reps_nodes: self.reps_nodes + other.reps_nodes,
        }
    }
}

impl Sum for StatsEncoding {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(
            StatsEncoding {
                nucleo_count: 0,
                targets_count: 0,
                reps_count: 0,
                nucleo_nodes: 0,
                targets_nodes: 0,
                reps_nodes: 0,
            },
            |acc, x| acc + x,
        )
    }
}

impl StatsEncoding {
    fn print(&self) {
        let total_count = self.nucleo_count + self.targets_count + self.reps_count;
        let total_nodes = self.nucleo_nodes + self.targets_nodes + self.reps_nodes;

        eprintln!("       Method | Number of extensions | Number of encoded nodes ");
        eprintln!("--------------|-------------------------|-----------------");
        for (name, count, nodes) in [
            ("Next node", self.nucleo_count, self.nucleo_nodes),
            ("Target node", self.targets_count, self.targets_nodes),
            ("Repetition", self.reps_count, self.reps_nodes),
        ] {
            eprintln!(
                "{:>13} | {:>14}  ({:>4.1}%) | {:>11} bits ({:>4.1}%)",
                name,
                format_int(count),
                count as f64 / total_count as f64 * 100.0,
                format_int(nodes),
                nodes as f64 / total_nodes as f64 * 100.0,
            );
        }
        eprintln!(
            "        Total | {:>14}          | {:>11} bits",
            format_int(total_count),
            format_int(total_nodes),
        );
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
    /// Get statistics about the underlying node encoding
    fn get_stats(&self, graph: &Graph<impl KmerStorage>) -> StatsEncoding {
        self.nodes_encoding.get_stats(graph)
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
        println!(">{}", self.id);
        let stats = self
            .contigs
            .iter()
            .map(|c| c.get_stats(graph))
            .sum::<StatsEncoding>();
        stats.print();
        eprintln!();
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
    /// Encode all records present in a fasta file
    pub fn encode_file(&self, input: &Path, output: &Path) -> Result<(), Box<dyn Error>> {
        let mut input_reader = needletail::parse_fastx_file(input)?;
        let mut output_writer = std::io::BufWriter::new(std::fs::File::create(output)?);

        while let Some(record) = input_reader.next() {
            let record = record?;
            let id = unsafe { String::from_utf8_unchecked(record.id().to_owned()) };
            let seq = record.seq();
            eprintln!("Encoding record {}", id);
            println!(">{}", id);
            let scaffold = self.encode_record(id, &seq)?;
            bincode::encode_into_std_write(
                scaffold,
                &mut output_writer,
                bincode::config::standard(),
            )?;
        }
        Ok(())
    }

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
