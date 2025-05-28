//! To create graphs from fasta/unitig files

use debruijn::compression;
use debruijn::graph::{BaseGraph, DebruijnGraph};
use debruijn::{Dir, Exts, Kmer, Vmer};

use boomphf::hashmap::BoomHashMap2;

use std::collections::HashSet;
use std::ops::{Deref, DerefMut};

use bincode::{deserialize_from, serialize_into};
use serde::{Deserialize, Serialize};

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use crate::fasta_reader::FastaReader;

/// Wrapper around the DebruijnGraph from the debruijn crate.
#[derive(Serialize, Deserialize, Debug)]
pub struct Graph<K: Kmer>(DebruijnGraph<K, ()>);

impl<K: Kmer> Deref for Graph<K> {
    type Target = DebruijnGraph<K, ()>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl<K: Kmer> DerefMut for Graph<K> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<K: Kmer> Graph<K> {
    /// Print some stats about the graph
    pub fn print_stats(&self) {
        let mut node_length: usize = 0;
        let mut nb_edges: usize = 0;
        let mut edges_histo: [[u32; 5]; 5] = [[0; 5]; 5];

        for node in self.iter_nodes() {
            node_length += node.len();
            let exts = node.exts();
            nb_edges += (exts.num_exts_l() + exts.num_exts_r()) as usize;
            edges_histo[exts.num_exts_l() as usize][exts.num_exts_r() as usize] += 1;
        }
        eprintln!(
            "Graph contains:\n   - {} nodes\n   - {} edges\nAverage node length: {}",
            self.len(),
            nb_edges / 2,
            node_length / self.len()
        );
        eprintln!("Edges histogram:");
        for i in 0..5 {
            for j in 0..5 {
                eprint!("{:>7} ", edges_histo[i][j]);
            }
            eprintln!();
        }
    }

    /// Search a kmer at one side of the nodes
    pub fn search_kmer(&self, kmer: K, side: Dir) -> Option<(usize, Dir)> {
        match side {
            Dir::Left => self
                .find_link(kmer, Dir::Right)
                .map(|(id, dir, _)| (id, dir)),
            Dir::Right => self
                .find_link(kmer, Dir::Left)
                .map(|(id, dir, _)| (id, dir.flip())),
        }
    }

    /// Search a kmer within the nodes, by iterating over all the graph. Returns ((node_id, dir), offset) where offset is counted from the given side.
    pub fn search_kmer_offset(&self, kmer: K, side: Dir) -> Option<((usize, Dir), usize)> {
        // search first at the given extremity
        if let Some((node_id, dir)) = self.search_kmer(kmer, side) {
            return Some(((node_id, dir), 0));
        }
        // if not found, iterate over all kmers
        let rc = kmer.rc();
        for node_id in 0..self.len() {
            for (offset, k) in self.get_node_kmer(node_id).into_iter().enumerate() {
                if k == kmer {
                    let offset = match side {
                        Dir::Left => offset,
                        Dir::Right => self.get_node(node_id).len() - K::k() - offset,
                    };
                    return Some(((node_id, side), offset));
                }
                // if not stranded, look for the reverse complement
                else if !self.base.stranded && k == rc {
                    let offset = match side {
                        Dir::Left => self.get_node(node_id).len() - K::k() - offset,
                        Dir::Right => offset,
                    };
                    return Some(((node_id, side.flip()), offset));
                }
            }
        }
        None
    }

    /// Utility function to get the correct extention of a node.
    /// Returns Some(exts) if the extensions are incorrect, None otherwise.
    fn get_new_exts(&self, node_id: usize) -> Option<(usize, Exts)> {
        let node = self.get_node(node_id);
        let l_kmer: K = node.sequence().first_kmer();
        let r_kmer: K = node.sequence().last_kmer();
        let mut new_exts = Exts::empty();

        for i in 0..4 {
            if self.find_link(l_kmer.extend_left(i), Dir::Left).is_some() {
                new_exts = new_exts.set(Dir::Left, i);
            }
            if self.find_link(r_kmer.extend_right(i), Dir::Right).is_some() {
                new_exts = new_exts.set(Dir::Right, i);
            }
        }
        if new_exts == node.exts() {
            return None;
        } else {
            return Some((node_id, new_exts));
        }
    }

    // Utility function to add all missing edges to the graph
    fn fix_exts_serial(&mut self) {
        // Get edges wich need update
        let updates: Vec<(usize, Exts)> = (0..self.len())
            .filter_map(|i| self.get_new_exts(i))
            .collect();

        // Apply all the updates to the original graph
        for (node_id, new_exts) in updates {
            self.base.exts[node_id] = new_exts;
        }
    }

    /// Create a graph from a sequence of kmers. (For debugging mainly)
    pub fn from_seq_serial(seq: &impl Vmer, stranded: bool) -> Self {
        let can = |k: K| {
            if stranded {
                k
            } else {
                k.min_rc()
            }
        };

        let kmer_set: HashSet<K> = seq.iter_kmers().map(|k| can(k)).collect();
        let mut keys = Vec::with_capacity(kmer_set.len());
        let mut exts = Vec::with_capacity(kmer_set.len());
        let data = vec![(); kmer_set.len()];

        for kmer in kmer_set.iter() {
            let mut e = Exts::empty();
            for base in 0..4 {
                let new = can(kmer.extend_left(base));
                if kmer_set.contains(&new) {
                    e = e.set(Dir::Left, base);
                }
                let new = can(kmer.extend_right(base));
                if kmer_set.contains(&new) {
                    e = e.set(Dir::Right, base);
                }
            }
            keys.push(*kmer);
            exts.push(e);
        }

        let spec = compression::ScmapCompress::<()>::new();
        let index = BoomHashMap2::new(keys, exts, data);
        let graph = compression::compress_kmers_with_hash(stranded, &spec, &index);
        Self(graph.finish_serial())
    }

    /// Create a graph from a fasta file containing unitigs (as returned by ggcat for example).
    pub fn from_unitigs_serial(path: &Path, stranded: bool) -> Self {
        let mut base_graph: BaseGraph<K, ()> = BaseGraph::new(stranded);

        // Iterate over unitigs and add them to the graph
        let fasta_reader = FastaReader::new(path).unwrap();
        fasta_reader.into_iter().for_each(|record| {
            let seq_bytes = record.dna_string().to_bytes();
            base_graph.add(seq_bytes, Exts::empty(), ());
        });

        // Finish the graph (computes boomphf to retrieve unitigs giving their edge kmers)
        let mut graph = Self(base_graph.finish_serial());

        // Update the links between the unitigs
        graph.fix_exts_serial();
        graph
    }
}

// Parallel construction
impl<K: Kmer + Send + Sync> Graph<K> {
    /// Create a graph from a fasta file containing unitigs (as returned by ggcat for example).
    pub fn from_unitigs(path: &Path, stranded: bool) -> Self {
        let mut base_graph: BaseGraph<K, ()> = BaseGraph::new(stranded);

        // Iterate over unitigs and add them to the graph
        let fasta_reader = FastaReader::new(path).unwrap(); // TODO: parallelise this
        fasta_reader.into_iter().for_each(|record| {
            let seq_bytes = record.dna_string().to_bytes();
            base_graph.add(seq_bytes, Exts::empty(), ());
        });

        // Finish the graph (computes boomphf to retrieve unitigs giving their edge kmers)
        let mut graph = Self(base_graph.finish());

        // Update the links between the unitigs
        graph.fix_exts_serial();
        graph
    }
}

// Dump and load from binary
impl<K: Kmer + Send + Sync + Serialize + for<'a> Deserialize<'a>> Graph<K> {
    /// Write the graph as a binary
    pub fn save_to_binary(&self, path_bin: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::create(path_bin)?;
        let mut writer = BufWriter::new(file);
        serialize_into(&mut writer, &self.base)?;
        Ok(())
    }
    /// Load the graph from a binary file.
    pub fn load_from_binary(path_bin: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(path_bin)?;
        let reader = BufReader::new(file);
        let base: BaseGraph<K, ()> = deserialize_from(reader)?;
        let graph = Self(base.finish());

        Ok(graph)
    }
}

#[cfg(test)]
mod unit_test {
    use std::path::PathBuf;

    use super::Graph;

    use debruijn::kmer::Kmer3;
    use debruijn::DnaSlice;

    const SEQ: DnaSlice = DnaSlice(&[2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 1]); // gggccccgggaaaaac
    const STRANDED: bool = true;
    const PATH_UNITIGS: &str = "../data/input/test.fna";

    #[test]
    #[ignore]
    fn test_from_seq() {
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);

        println!("{:?}", graph.base.sequences);
        graph.print();
    }

    #[test]
    #[ignore]
    fn test_from_unitigs() {
        let graph = Graph::<Kmer3>::from_unitigs_serial(&PathBuf::from(PATH_UNITIGS), STRANDED);

        println!("{:?}", graph.base.sequences);
        graph.print();
    }
}
