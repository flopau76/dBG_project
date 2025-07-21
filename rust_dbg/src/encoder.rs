//! Encodes various dna sequence into a graph, from the simple path to the scaffold with gaps.

use derive_more::with_trait::{Add, AddAssign, Sum};

use std::{
    fmt::{Debug, Display},
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

use bincode::{Decode, Encode};
use needletail::Sequence;
use packed_seq::{PackedSeqVec, SeqVec};

use crate::{Graph, KmerStorage, Node, NodeIterator, PathwayError};

mod basic_encoder;
mod gnome_encoder;
mod greedy_encoder;
mod serialize;

pub use basic_encoder::BasicEncoder;
pub use gnome_encoder::GnomeEncoder;
pub use greedy_encoder::{Extension, GreedyEncoder, VecExtensions};

//####################################################################################
//                                    Encoding                                      //
//####################################################################################

/// The encoding of a path in a de Bruijn graph.
pub trait Encoding: Sized + Debug {
    type Stats: EncodingStats;

    /// Decode the data into a list of nodes.
    fn decode(&self, graph: &Graph<impl KmerStorage>) -> Vec<Node>;

    /// Get some statistics.
    fn get_encoding_stats(&self, graph: &Graph<impl KmerStorage>) -> Self::Stats;
}

/// Some stats describing an encoding.
pub trait EncodingStats: Sized + AddAssign + Add<Self> + Sum<Self> + Display + Default {}

//####################################################################################
//                                    Encoder                                       //
//####################################################################################

/// A structure transforming a path in a graph into an encoding.
pub trait Encoder: Sized {
    type Encoding: Encoding + Encode + Decode<()>;

    /// Transform a list of nodes into any type of encoding
    fn encode_path(&self, nodes: Vec<Node>, graph: &Graph<impl KmerStorage>) -> Self::Encoding;

    /// Encode a fasta file into a binary format using the provided graph.
    fn encode_from_fasta(
        &self,
        input: &Path,
        output: &Path,
        graph: &Graph<impl KmerStorage>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut input_reader = needletail::parse_fastx_file(input)?;
        let mut output_writer = BufWriter::new(File::create(output)?);

        while let Some(record) = input_reader.next() {
            let record = record?;
            let id = unsafe { String::from_utf8_unchecked(record.id().to_owned()) };
            let seq = record.seq();
            eprintln!("Encoding record {}", id);
            println!(">{}", id);
            let mut scaffold: Scaffold<Self::Encoding> = Scaffold::from_seq(&seq, graph, self)?;
            scaffold.id = id;
            bincode::encode_into_std_write(
                scaffold,
                &mut output_writer,
                bincode::config::standard(),
            )?;
        }
        Ok(())
    }

    /// Decode a binary file into a fasta format using the provided graph.
    fn decode_to_fasta(
        &self,
        input: &Path,
        output: &Path,
        graph: &Graph<impl KmerStorage>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut input_reader = BufReader::new(File::open(input)?);
        let mut output_writer = BufWriter::new(File::create(output)?);

        while let Ok(scaffold) = bincode::decode_from_std_read::<Scaffold<Self::Encoding>, _, _>(
            &mut input_reader,
            bincode::config::standard(),
        ) {
            eprintln!("Decoding scaffold: {}", scaffold.id);
            let seq = scaffold.get_seq(graph);
            output_writer.write_all(b">")?;
            output_writer.write_all(scaffold.id.as_bytes())?;
            output_writer.write_all(b"\n")?;
            output_writer.write_all(&seq)?;
            output_writer.write_all(b"\n")?;
        }
        Ok(())
    }

    /// Print some stats about the encoding.
    fn print_stats(
        &self,
        input: &Path,
        graph: &Graph<impl KmerStorage>,
    ) -> Result<(), Box<dyn std::error::Error>>
    where
        Self::Encoding: Decode<()>,
    {
        let mut input_reader = BufReader::new(File::open(input)?);

        let mut stats = <Self::Encoding as Encoding>::Stats::default();
        while let Ok(scaffold) = bincode::decode_from_std_read::<Scaffold<Self::Encoding>, _, _>(
            &mut input_reader,
            bincode::config::standard(),
        ) {
            stats += scaffold.get_stats(&graph);
        }
        eprintln!("Stats for {}: ", input.display());
        eprintln!("{}", stats);

        Ok(())
    }
}

//####################################################################################
//                                   Embedding                                      //
//####################################################################################

/// Any sequence, which can be embedded into a graph as a combination of multiple encodings.
pub trait Embedding<P: Encoding>: Sized {
    /// Retrieve the embeddded sequence.
    fn get_seq(&self, graph: &Graph<impl KmerStorage>) -> Vec<u8>;

    /// Embed a sequence into a graph.
    fn from_seq(
        seq: &[u8],
        graph: &Graph<impl KmerStorage>,
        encoder: &impl Encoder<Encoding = P>,
    ) -> Result<Self, PathwayError>;

    // Get some statistics.
    fn get_stats(&self, graph: &Graph<impl KmerStorage>) -> P::Stats;
}

impl<P: Encoding> Embedding<P> for P {
    fn get_seq(&self, graph: &Graph<impl KmerStorage>) -> Vec<u8> {
        let nodes = self.decode(graph);
        graph.path_seq(&nodes)
    }
    fn from_seq(
        seq: &[u8],
        graph: &Graph<impl KmerStorage>,
        encoder: &impl Encoder<Encoding = P>,
    ) -> Result<Self, PathwayError> {
        let seq = PackedSeqVec::from_ascii(seq);
        let mut iterator = NodeIterator::new(graph, seq)?;
        let mut nodes = Vec::default();
        while let Some(node) = iterator.next()? {
            nodes.push(node);
        }
        Ok(encoder.encode_path(nodes, graph))
    }
    fn get_stats(&self, graph: &Graph<impl KmerStorage>) -> P::Stats {
        self.get_encoding_stats(graph)
    }
}

//####################################################################################
//                                   Contig                                         //
//####################################################################################

/// Contig in a de Bruijn graph, represented as a path, along with a start and end offsets.
#[derive(Encode, Decode, Debug, Clone)]
pub struct Contig<P: Encoding> {
    pub nodes: P,
    pub start_offset: usize,
    pub end_offset: usize,
}

impl<P: Encoding> Embedding<P> for Contig<P> {
    fn get_seq(&self, graph: &Graph<impl KmerStorage>) -> Vec<u8> {
        let seq = self.nodes.get_seq(graph);
        seq[self.start_offset..seq.len() - self.end_offset].to_owned()
    }

    fn from_seq(
        seq: &[u8],
        graph: &Graph<impl KmerStorage>,
        encoder: &impl Encoder<Encoding = P>,
    ) -> Result<Self, PathwayError> {
        let seq = PackedSeqVec::from_ascii(seq);
        let mut iterator = NodeIterator::new(graph, seq)?;
        let mut nodes = Vec::default();
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

    fn get_stats(&self, graph: &Graph<impl KmerStorage>) -> P::Stats {
        self.nodes.get_stats(graph)
    }
}

//####################################################################################
//                                  Scaffold                                        //
//####################################################################################

/// Scaffold in a de Bruijn graph, represented as an alternance of contigs and gaps of a given size.
#[derive(Encode, Decode, Debug, Clone)]
pub struct Scaffold<P: Encoding> {
    pub id: String,
    pub contigs: Vec<Contig<P>>,
    pub gaps: Vec<usize>,
}

impl<P: Encoding> Embedding<P> for Scaffold<P> {
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
        encoder: &impl Encoder<Encoding = P>,
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

    fn get_stats(&self, graph: &Graph<impl KmerStorage>) -> P::Stats {
        self.contigs.iter().map(|c| c.get_stats(graph)).sum()
    }
}
