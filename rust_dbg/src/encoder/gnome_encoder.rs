use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use crate::{Graph, KmerStorage, Node};

use super::{BasicEncoder, Embedding, Scaffold};

pub struct GnomeEncoder {
    g: usize,
    dict: BTreeMap<Vec<Node>, Vec<(u8, u8, u32)>>,
    record_names: Vec<String>,
    is_init: bool,
}

impl Default for GnomeEncoder {
    fn default() -> Self {
        GnomeEncoder {
            g: 10,
            dict: BTreeMap::new(),
            record_names: Vec::new(),
            is_init: false,
        }
    }
}

impl GnomeEncoder {
    pub fn new(g: usize) -> Self {
        GnomeEncoder {
            g,
            dict: BTreeMap::new(),
            record_names: Vec::new(),
            is_init: false,
        }
    }

    pub fn encode_from_fasta(
        &mut self,
        input: &Path,
        output: &Path,
        graph: &Graph<impl KmerStorage>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // let mut input_reader = needletail::parse_fastx_file(input)?;
        let mut output_writer = BufWriter::new(File::create(output)?);

        if !self.is_init {
            self.init_from_fasta(input, graph)?;
        }

        // Write stats to the output file
        writeln!(output_writer, "# input fasta: {}", input.display())?;
        writeln!(output_writer, "# g={}", self.g)?;
        self.print_stats_gnomes(&mut output_writer)?;
        output_writer.flush()?;

        Ok(())
    }

    /// Initialize the encoder from a fasta file, creating a dictionary of g-node-mers.
    /// The dictionary maps g-node-mer sequences to their occurrences in the scaffolds.
    /// Each entry in the dictionary is a vector of tuples (scaffold_index, contig_index, position_in_contig).
    fn init_from_fasta(
        &mut self,
        input: &Path,
        graph: &Graph<impl KmerStorage>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        if self.is_init {
            return Ok(());
        }
        eprint!(
            "Initializing GnomeEncoder with g = {}, from fasta {}...",
            self.g,
            input.display()
        );
        let mut input_reader = needletail::parse_fastx_file(input)?;

        let mut i = 0;
        while let Some(record) = input_reader.next() {
            let record = record?;
            let record_name = unsafe { String::from_utf8_unchecked(record.id().to_owned()) };
            self.record_names.push(record_name);
            let record_seq = record.seq();
            let scaffold = Scaffold::from_seq(&record_seq, graph, &BasicEncoder::default())?;
            for (j, contig) in scaffold.contigs.iter().enumerate() {
                let gnome_iter = contig.nodes.windows(self.g);
                for (k, gnome) in gnome_iter.enumerate() {
                    self.dict
                        .entry(gnome.to_vec())
                        .or_insert_with(Vec::new)
                        .push((i, j as u8, k as u32));
                }
            }
            i += 1;
        }
        self.is_init = true;
        eprintln!("done",);
        Ok(())
    }

    // Creates a tsv with the number of occurrences of each gnome, per scaffold.
    fn print_stats_gnomes(
        &self,
        output_writer: &mut BufWriter<File>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // group occurences per haplotype
        let mut dict_haplo = self
            .dict
            .iter()
            .map(|(_k, v)| {
                let mut haplo_counts = vec![0 as u16; self.record_names.len()];
                for (i, _j, _k) in v {
                    // increment the count for the corresponding scaffold
                    haplo_counts[*i as usize] += 1;
                }
                haplo_counts
            })
            .collect::<Vec<_>>();

        // sort by decreasing number of occurrences
        dict_haplo.sort_by(|a, b| {
            let sum_a: u16 = a.iter().sum();
            let sum_b: u16 = b.iter().sum();
            match sum_b.cmp(&sum_a) {
                std::cmp::Ordering::Equal => a.cmp(b),
                ord => ord,
            }
        });

        writeln!(
            output_writer,
            "# Number of records: {}",
            self.record_names.len()
        )?;
        writeln!(
            output_writer,
            "# Number of different gnomes: {}",
            self.dict.len()
        )?;

        for (i, name) in self.record_names.iter().enumerate() {
            if i > 0 {
                write!(output_writer, "\t")?;
            }
            write!(output_writer, "{}", name)?;
        }
        writeln!(output_writer)?;

        for v in &dict_haplo {
            for (i, count) in v.into_iter().enumerate() {
                if i > 0 {
                    write!(output_writer, "\t")?;
                }
                write!(output_writer, "{}", count)?;
            }
            writeln!(output_writer)?;
        }

        Ok(())
    }

    /// Find all occurrences of a g-node-mer of size <=g in the dictionary.
    fn find_gnome(&self, gnome: &[Node]) -> Vec<(u8, u8, u32)> {
        if gnome.len() > self.g {
            panic!(
                "G-node-mer length {} exceeds the maximum allowed size {}",
                gnome.len(),
                self.g
            );
        }
        if gnome.len() == self.g {
            return self.dict.get(gnome).unwrap_or(&vec![]).to_vec();
        }
        let mut left = gnome.to_vec();
        left.extend_from_slice(&vec![Node::MIN; self.g - gnome.len()]);
        let mut right = gnome.to_vec();
        right.extend_from_slice(&vec![Node::MAX; self.g - gnome.len()]);
        self.dict
            .range(left..=right)
            .flat_map(|(_k, v)| v.to_vec())
            .collect::<Vec<_>>()
    }
}
