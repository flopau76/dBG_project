use std::collections::{BTreeMap, HashMap, VecDeque};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

use bincode::{Decode, Encode};
use packed_seq::Seq;

use crate::embeddings::{
    Embedding, Extension, Pathway, PathwayStats, Scaffold, VecExtensions, VecNodes,
};
use crate::graph::shortest_path;
use crate::{Graph, KmerStorage, Node};

//####################################################################################
//                                 Greedy Encoder                                   //
//####################################################################################

pub trait Encoder<P: Pathway>: Sized {
    /// Transform a list of nodes into any type of encoding
    fn encode_path(&self, nodes: VecNodes, graph: &Graph<impl KmerStorage>) -> P;

    /// Encode a fasta file into a binary format using the provided graph.
    fn encode_from_fasta(
        &self,
        input: &Path,
        output: &Path,
        graph: &Graph<impl KmerStorage>,
    ) -> Result<(), Box<dyn std::error::Error>>
    where
        P: Encode,
    {
        let mut input_reader = needletail::parse_fastx_file(input)?;
        let mut output_writer = BufWriter::new(File::create(output)?);

        while let Some(record) = input_reader.next() {
            let record = record?;
            let id = unsafe { String::from_utf8_unchecked(record.id().to_owned()) };
            let seq = record.seq();
            eprintln!("Encoding record {}", id);
            println!(">{}", id);
            let mut scaffold = Scaffold::from_seq(&seq, graph, self)?;
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
    ) -> Result<(), Box<dyn std::error::Error>>
    where
        P: Decode<()>,
    {
        let mut input_reader = BufReader::new(File::open(input)?);
        let mut output_writer = BufWriter::new(File::create(output)?);

        while let Ok(scaffold) = bincode::decode_from_std_read::<Scaffold<VecExtensions>, _, _>(
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

    /// Get some stats about the binary encoding.
    fn print_stats(
        &self,
        input: &Path,
        graph: &Graph<impl KmerStorage>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut input_reader = BufReader::new(File::open(input)?);

        let mut stats = PathwayStats::default();
        while let Ok(scaffold) = bincode::decode_from_std_read::<Scaffold<VecExtensions>, _, _>(
            &mut input_reader,
            bincode::config::standard(),
        ) {
            stats += scaffold.get_stats(&graph);
        }
        stats.print();

        Ok(())
    }
}

/// Basic encoder, which transforms a list of nodes into a sequence of extensions.
#[derive(Default)]
pub struct BasicEncoder;
impl Encoder<VecNodes> for BasicEncoder {
    fn encode_path(&self, nodes: VecNodes, _graph: &Graph<impl KmerStorage>) -> VecNodes {
        nodes
    }
}

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
                    println!("{:?}", ext);
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
                                println!("{:?}", ext);
                                ext
                            }),
                    );
                }
                current_position += length as usize;
            }
            // add the repetition
            let ext = Extension::Repetition(rep.0, rep.1);
            extensions.push(ext);
            println!("{:?}", ext);
            current_position += rep.0 as usize;
        }
        extensions.pop(); // remove the sentinel
        VecExtensions(extensions)
    }
}

//####################################################################################
//                                   G-node-mer                                     //
//####################################################################################

use textplots::{Chart, Plot, Shape};

pub struct GnomeEncoder {
    g: usize,
    dict: HashMap<Vec<Node>, Vec<(usize, usize, usize)>>,
    is_init: bool,
}

impl Default for GnomeEncoder {
    fn default() -> Self {
        GnomeEncoder {
            g: 10,
            dict: HashMap::new(),
            is_init: false,
        }
    }
}

impl GnomeEncoder {
    pub fn encode_from_fasta(
        &mut self,
        input: &Path,
        output: &Path,
        graph: &Graph<impl KmerStorage>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut input_reader = needletail::parse_fastx_file(input)?;
        // let mut output_writer = BufWriter::new(File::create(output)?);

        // Create a dictionary containing all g-node-mers from the fasta
        let mut i = 0;
        while let Some(record) = input_reader.next() {
            let record = record?;
            let id = unsafe { String::from_utf8_unchecked(record.id().to_owned()) };
            let seq = record.seq();
            eprintln!("Encoding record {}", id);
            println!(">{}", id);
            let scaffold = Scaffold::from_seq(&seq, graph, &BasicEncoder::default())?;
            for (j, contig) in scaffold.contigs.iter().enumerate() {
                let gnome_iter = contig.nodes.windows(self.g);
                for (k, gnome) in gnome_iter.enumerate() {
                    self.dict
                        .entry(gnome.to_vec())
                        .or_insert_with(Vec::new)
                        .push((i, j, k));
                }
            }
            i += 1;
        }
        // inspect the dictionary
        let data: BTreeMap<usize, u32> = self.dict.values().fold(BTreeMap::new(), |mut acc, v| {
            *acc.entry(v.len()).or_insert(0) += 1;
            acc
        });
        for (k, v) in &data {
            println!("{}: {}", k, v);
        }
        let data = data
            .into_iter()
            .map(|(k, v)| (k as f32, v as f32))
            .collect::<Vec<_>>();
        let max = data.iter().map(|(_, v)| *v).fold(0.0, f32::max);
        let histo = textplots::utils::histogram(&data, 0.0, max, 20);
        println!("\ny = histogram bars");
        Chart::new(180, 60, 6.0, 14.0)
            .lineplot(&Shape::Bars(&histo))
            .nice();
        println!("\ny = bars");
        Chart::new(180, 60, 0.0, 100.0)
            .lineplot(&Shape::Bars(&data))
            .display();

        Ok(())
    }
}
