//! To create graphs from fasta/unitig files

use debruijn::{Kmer, Vmer, Exts, Dir};
use debruijn::graph::{BaseGraph, DebruijnGraph};
use debruijn::compression;

use ahash::AHashSet;
use std::ops::{Deref, DerefMut};
use serde::{Serialize, Deserialize};

use std::sync::{Arc, Mutex};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::fasta_reader::FastaReader;

/// Wrapper around the DebruijnGraph from the debruijn crate.
#[derive(Serialize, Deserialize, Debug)]
pub struct Graph<K: Kmer> (DebruijnGraph<K, ()>);

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

/// Serial construction
impl<K: Kmer> Graph<K> {

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
            .filter_map(|i| {self.get_new_exts(i)})
            .collect();

        // Apply all the updates to the original graph
        for (node_id, new_exts) in updates {
            self.base.exts[node_id] = new_exts;
        }
    }
    
    /// Creates a graph from a sequence of kmers. (For debugging mainly)
    pub fn from_seq_serial(seq: impl Vmer, stranded: bool) -> Self {
        let can = |k: K| {if stranded {k} else {k.min_rc()}};
        let unique_kmers = seq.iter_kmers().map(|k| can(k)).collect::<AHashSet<K>>().into_iter().map(|k| (k, ())).collect::<Vec<_>>();
        let compression = compression::ScmapCompress::<()>::new();
        let graph = compression::compress_kmers_no_exts(stranded, &compression,  &unique_kmers);
        Self(graph.finish_serial())
    }

    /// Creates a graph from a fasta file containing unitigs (as returned by ggcat for example).
    pub fn from_unitigs_serial(path:&str, stranded: bool) -> Graph<K> {
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

/// Parallel construction
impl<K: Kmer+Send+Sync> Graph<K> {

    /// Same as [Graph::fix_exts_serial] but with some parallelisation.
    fn fix_exts(&mut self) {
        // Get edges wich need update
        let updates: Vec<(usize, Exts)> = (0..self.len()).into_par_iter()
            .filter_map(|i| {self.get_new_exts(i)})
            .collect();

        // Apply all the updates to the original graph
        for (node_id, new_exts) in updates {
            self.base.exts[node_id] = new_exts;
        }
    }

    /// Same as [from_unitigs_serial](Graph::from_unitigs_serial) but with some parallelisation.
    pub fn from_unitigs(path:&str, stranded: bool) -> Self {
        let base_graph: BaseGraph<K, ()> = BaseGraph::new(stranded);
        let base_graph = Arc::new(Mutex::new(base_graph));

        // Iterate over unitigs and add them to the graph
        let fasta_reader = FastaReader::new(path).unwrap();
        fasta_reader.into_par_iter().for_each(|record| {
            let seq_bytes = record.dna_string().to_bytes();
            let mut graph = base_graph.lock().unwrap();
            graph.add(seq_bytes, Exts::empty(), ());
        });

        // Finish the graph (computes boomphf to retrieve unitigs giving their edge kmers)
        let base_graph = Arc::into_inner(base_graph).unwrap().into_inner().unwrap();
        let mut graph = Self(base_graph.finish());

        // Update the links between the unitigs
        graph.fix_exts();
        graph
    }
}

#[cfg(test)]
mod unit_test {
    use super::Graph;

    use debruijn::kmer::Kmer3;
    use debruijn::DnaSlice;

    const SEQ: DnaSlice = DnaSlice(&[2,2,2,1,1,1,1,2,2,2,0,0,0,0,0,1]);    // gggccccgggaaaaac
    const STRANDED : bool = true;
    const PATH_UNITIGS : &str = "../data/input/test.fna";

    #[test]
    #[ignore]
    fn test_from_seq() {
        let graph = Graph::<Kmer3>::from_seq_serial(SEQ, STRANDED);

        println!("{:?}", graph.base.sequences);
        graph.print();
    }

    #[test]
    #[ignore]
    fn test_from_unitigs() {
        let graph = Graph::<Kmer3>::from_unitigs_serial(PATH_UNITIGS, STRANDED);

        println!("{:?}", graph.base.sequences);
        graph.print();
    }
}