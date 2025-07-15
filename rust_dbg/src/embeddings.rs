//! Embedds various dna sequence into a graph, from the simple path to the scaffold with gaps.

use std::fmt::Debug;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Deref, DerefMut};

use bincode::{Decode, Encode};
use needletail::Sequence;
use packed_seq::{PackedSeqVec, SeqVec};

use crate::encoder::Encoder;
use crate::graph::shortest_path;
use crate::{Graph, KmerStorage, Node, NodeIterator, PathwayError, Side};

//####################################################################################
//                           traits Pathway & Embedding                             //
//####################################################################################

/// A way to represent a continuous suite of nodes in a graph.
pub trait Pathway: Sized + Debug {
    /// Decode the pathway into a suite of nodes.
    fn decode(&self, graph: &Graph<impl KmerStorage>) -> VecNodes;

    /// Encode a pathway from a suite of nodes.
    fn encode(
        nodes: VecNodes,
        graph: &Graph<impl KmerStorage>,
        encoder: &impl Encoder<Self>,
    ) -> Self {
        encoder.encode_path(nodes, graph)
    }

    /// Get some statistics.
    fn get_stats_p(&self, graph: &Graph<impl KmerStorage>) -> PathwayStats;
}

/// A way to embedd any sequence into a graph, as a combination of pathways.
pub trait Embedding<P: Pathway>: Sized {
    /// Retrieve the embeddded sequence.
    fn get_seq(&self, graph: &Graph<impl KmerStorage>) -> Vec<u8>;

    /// Embedd a sequence into a graph.
    fn from_seq(
        seq: &[u8],
        graph: &Graph<impl KmerStorage>,
        encoder: &impl Encoder<P>,
    ) -> Result<Self, PathwayError>;

    // Get some statistics.
    fn get_stats(&self, graph: &Graph<impl KmerStorage>) -> PathwayStats;
}

impl<P: Pathway> Embedding<P> for P {
    fn get_seq(&self, graph: &Graph<impl KmerStorage>) -> Vec<u8> {
        let nodes = self.decode(graph);
        graph.path_seq(&nodes)
    }
    fn from_seq(
        seq: &[u8],
        graph: &Graph<impl KmerStorage>,
        encoder: &impl Encoder<P>,
    ) -> Result<Self, PathwayError> {
        let seq = PackedSeqVec::from_ascii(seq);
        let mut iterator = NodeIterator::new(graph, seq)?;
        let mut nodes = VecNodes::default();
        while let Some(node) = iterator.next()? {
            nodes.push(node);
        }
        Ok(encoder.encode_path(nodes, graph))
    }
    fn get_stats(&self, graph: &Graph<impl KmerStorage>) -> PathwayStats {
        self.get_stats_p(graph)
    }
}

//####################################################################################
//                                  VecNodes                                        //
//####################################################################################

#[derive(Default, Debug)]
pub struct VecNodes(pub Vec<Node>);
impl Deref for VecNodes {
    type Target = Vec<Node>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for VecNodes {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Pathway for VecNodes {
    fn decode(&self, _graph: &Graph<impl KmerStorage>) -> VecNodes {
        VecNodes(self.0.clone())
    }
    fn get_stats_p(&self, _graph: &Graph<impl KmerStorage>) -> PathwayStats {
        PathwayStats {
            nucleo_count: self.0.len() - 1,
            targets_count: 1,
            reps_count: 0,
            nucleo_nodes: self.0.len(),
            targets_nodes: 1,
            reps_nodes: 0,
        }
    }
}

//####################################################################################
//                                  VecExtensions                                   //
//####################################################################################

/// An enum containing different ways to encode a path extension.
#[derive(Encode, Decode, Copy, Clone, PartialEq, Eq)]
pub enum Extension {
    TargetNode(Node),
    NextNucleotide(u8),
    Repetition(u16, u8), // nb_repeats, offset
}

impl Debug for Extension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Extension::TargetNode(node) => write!(f, "TN({:?})", node),
            Extension::NextNucleotide(base) => write!(f, "NN({})", base),
            Extension::Repetition(nb_repeats, offset) => {
                write!(f, "R({}-{})", nb_repeats, offset)
            }
        }
    }
}

#[derive(Default, Debug)]
pub struct VecExtensions(pub Vec<Extension>);
impl Deref for VecExtensions {
    type Target = Vec<Extension>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for VecExtensions {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Pathway for VecExtensions {
    fn decode(&self, graph: &Graph<impl KmerStorage>) -> VecNodes {
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
                        shortest_path::get_shortest_path(
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
        VecNodes(nodes)
    }

    fn get_stats_p(&self, graph: &Graph<impl KmerStorage>) -> PathwayStats {
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

        PathwayStats {
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
//                                   Contig                                         //
//####################################################################################

/// Contig in a de Bruijn graph, represented as a pathway of nodes with start and end offsets.
#[derive(Encode, Decode, Debug, Clone)]
pub struct Contig<P: Pathway> {
    pub nodes: P,
    pub start_offset: usize,
    pub end_offset: usize,
}

impl<P: Pathway> Embedding<P> for Contig<P> {
    fn get_seq(&self, graph: &Graph<impl KmerStorage>) -> Vec<u8> {
        let seq = self.nodes.get_seq(graph);
        seq[self.start_offset..seq.len() - self.end_offset].to_owned()
    }

    fn from_seq(
        seq: &[u8],
        graph: &Graph<impl KmerStorage>,
        encoder: &impl Encoder<P>,
    ) -> Result<Self, PathwayError> {
        let seq = PackedSeqVec::from_ascii(seq);
        let mut iterator = NodeIterator::new(graph, seq)?;
        let mut nodes = VecNodes::default();
        while let Some(node) = iterator.next()? {
            nodes.push(node);
        }
        let nodes = encoder.encode_path(nodes, graph);
        Ok(Contig {
            nodes,
            start_offset: iterator.start_offset.expect("Start offset should be set"),
            end_offset: iterator.end_offset.expect("End offset should be set"),
        })
    }

    fn get_stats(&self, graph: &Graph<impl KmerStorage>) -> PathwayStats {
        self.nodes.get_stats(graph)
    }
}

//####################################################################################
//                                  Scaffold                                        //
//####################################################################################

/// Scaffold in a de Bruijn graph, represented as an alternance of contigs and gaps of a given size.
#[derive(Encode, Decode, Debug, Clone)]
pub struct Scaffold<P: Pathway> {
    pub id: String,
    pub contigs: Vec<Contig<P>>,
    pub gaps: Vec<usize>,
}

impl<P: Pathway> Embedding<P> for Scaffold<P> {
    fn get_seq(&self, graph: &Graph<impl KmerStorage>) -> Vec<u8> {
        assert_eq!(
            self.contigs.len() + 1,
            self.gaps.len(),
            "Expected one more gap than contigs, but got {} contigs and {} gaps",
            self.contigs.len(),
            self.gaps.len()
        );
        let mut seq = vec![b'N'; self.gaps[0]];
        for (contig, gap) in self.contigs.iter().zip(&self.gaps[1..]) {
            seq.extend(&contig.get_seq(graph));
            seq.extend(std::iter::repeat(b'N').take(*gap));
        }
        seq
    }

    fn from_seq(
        seq: &[u8],
        graph: &Graph<impl KmerStorage>,
        encoder: &impl Encoder<P>,
    ) -> Result<Self, PathwayError> {
        let mut contigs = Vec::new();
        let mut gaps = Vec::new();
        let mut gap_size = 0;
        for contig in seq.normalize(false).split(|c| *c == b'N') {
            if contig.is_empty() {
                gap_size += 1;
                continue;
            };
            if contig.len() < graph.k() {
                eprintln!(
                    "[warning] Contig of length {} < k has been skipped",
                    contig.len(),
                );
                gap_size += contig.len();
                continue;
            }
            gaps.push(gap_size);
            let contig = Contig::from_seq(contig, graph, encoder)?;
            contigs.push(contig);
            gap_size = 0;
        }
        gaps.push(gap_size);

        Ok(Scaffold {
            id: String::default(),
            contigs,
            gaps,
        })
    }

    fn get_stats(&self, graph: &Graph<impl KmerStorage>) -> PathwayStats {
        self.contigs
            .iter()
            .map(|c| c.get_stats(graph))
            .sum::<PathwayStats>()
    }
}

//####################################################################################
//                               PathwayStats                                       //
//####################################################################################

#[derive(Default, Copy, Clone)]
pub struct PathwayStats {
    nucleo_count: usize,
    targets_count: usize,
    reps_count: usize,
    nucleo_nodes: usize,
    targets_nodes: usize,
    reps_nodes: usize,
}

impl AddAssign for PathwayStats {
    fn add_assign(&mut self, other: Self) {
        self.nucleo_count += other.nucleo_count;
        self.targets_count += other.targets_count;
        self.reps_count += other.reps_count;
        self.nucleo_nodes += other.nucleo_nodes;
        self.targets_nodes += other.targets_nodes;
        self.reps_nodes += other.reps_nodes;
    }
}

impl Add for PathwayStats {
    type Output = Self;

    fn add(mut self, other: Self) -> Self {
        self += other;
        self
    }
}

impl Sum for PathwayStats {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(PathwayStats::default(), |acc, x| acc + x)
    }
}

/// Format a long integer with commas
fn format_int(n: usize) -> String {
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

impl PathwayStats {
    pub fn print(&self) {
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
