//! To create graphs from fasta/unitig files

use debruijn::{Kmer, Vmer, Exts, Dir};
use debruijn::base_to_bits;
use debruijn::graph::{BaseGraph, DebruijnGraph};
use debruijn::compression;

use ahash::AHashSet;

use std::sync::{Arc, Mutex};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::fasta_reader::FastaReader;


pub type Graph<K> = DebruijnGraph<K, ()>;

/// Creates a graph from a sequence of kmers. (For debugging mainly)
pub fn graph_from_seq_serial<K: Kmer>(seq: impl Vmer, stranded: bool) -> Graph<K> {
    let can = |k: K| {if stranded {k} else {k.min_rc()}};
    let unique_kmers = seq.iter_kmers().map(|k| can(k)).collect::<AHashSet<K>>().into_iter().map(|k| (k, ())).collect::<Vec<_>>();
    let compression = compression::ScmapCompress::<()>::new();
    let graph = compression::compress_kmers_no_exts(stranded, &compression,  &unique_kmers);
    graph.finish_serial()
}

/// Create a graph from a fasta file containing unitigs (as returned by ggcat for example).
pub fn graph_from_unitigs_serial<K: Kmer>(path:&str, stranded: bool) -> Graph<K> {
    let mut graph: BaseGraph<K, ()> = BaseGraph::new(stranded);

    // Iterate over unitigs and add them to the graph
    let fasta_reader = FastaReader::new(path).unwrap();
    fasta_reader.into_iter().for_each(|record| {
        let seq_ascii = record.sequence();
        let seq_bytes = seq_ascii.iter().map(|&b| base_to_bits(b)).collect::<Vec<u8>>();
        graph.add(seq_bytes, Exts::empty(), ());
    });

    // Finish the graph (computes boomphf to retrieve unitigs giving their edge kmers)
    let mut db_graph = graph.finish_serial();

    // Update the links between the unitigs
    fix_exts_serial(&mut db_graph);
    db_graph
}

/// Same as `graph_from_unitigs_serial` with some parallelisation.
pub fn graph_from_unitigs<K: Kmer+Send+Sync>(path:&str, stranded: bool) -> Graph<K> {
    let graph: BaseGraph<K, ()> = BaseGraph::new(stranded);
    let graph_arc = Arc::new(Mutex::new(graph));

    // Iterate over unitigs and add them to the graph
    let fasta_reader = FastaReader::new(path).unwrap();
    fasta_reader.into_par_iter().for_each(|record| {
        let seq_ascii = record.sequence();
        let seq_bytes = seq_ascii.iter().map(|&b| base_to_bits(b)).collect::<Vec<u8>>();
        let mut graph = graph_arc.lock().unwrap();
        graph.add(seq_bytes, Exts::empty(), ());
    });

    // Finish the graph (computes boomphf to retrieve unitigs giving their edge kmers)
    let graph = Arc::into_inner(graph_arc).unwrap().into_inner().unwrap();
    let mut db_graph = graph.finish();

    // Update the links between the unitigs
    fix_exts(&mut db_graph);
    db_graph
}

// utilitary function to add missing edges to the graph
fn fix_exts_serial<K: Kmer>(graph: &mut Graph<K>) {
    // Get edges wich need update
    let updates: Vec<(usize, Exts)> = (0..graph.len())
        .filter_map(|i| {get_new_exts(&graph, i)})
        .collect();

    // Apply all the updates to the original graph
    for (node_id, new_exts) in updates {
        graph.base.exts[node_id] = new_exts;
    }
}

fn fix_exts<K: Kmer+Send+Sync>(graph: &mut Graph<K>) {
    // Get edges wich need update
    let updates: Vec<(usize, Exts)> = (0..graph.len()).into_par_iter()
        .filter_map(|i| {get_new_exts(&graph, i)})
        .collect();

    // Apply all the updates to the original graph
    for (node_id, new_exts) in updates {
        graph.base.exts[node_id] = new_exts;
    }
}

/// Checks if the edges of a node are correct. If not, return the new edges.
fn get_new_exts<K: Kmer>(graph: &Graph<K>, node_id: usize) -> Option<(usize, Exts)> {
    let node = graph.get_node(node_id);
    let l_kmer: K = node.sequence().first_kmer();
    let r_kmer: K = node.sequence().last_kmer();
    let mut new_exts = Exts::empty();

    for i in 0..4 {
        if graph.find_link(l_kmer.extend_left(i), Dir::Left).is_some() {
            new_exts = new_exts.set(Dir::Left, i);
        }
        if graph.find_link(r_kmer.extend_right(i), Dir::Right).is_some() {
            new_exts = new_exts.set(Dir::Right, i);
        }
    }
    if new_exts == node.exts() {
        return None;
    } else {
        return Some((node_id, new_exts));
    }
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use debruijn::kmer::Kmer3;
    use debruijn::DnaSlice;

    const SEQ: DnaSlice = DnaSlice(&[2,2,2,1,1,1,1,2,2,2,0,0,0,0,0,1]);    // gggccccgggaaaaac
    const STRANDED : bool = true;
    const PATH_UNITIGS : &str = "../data/input/test.fna";

    #[test]
    #[ignore]
    fn test_from_seq() {
        let graph = graph_from_seq_serial::<Kmer3>(SEQ, STRANDED);

        println!("{:?}", graph.base.sequences);
        graph.print();
    }

    #[test]
    #[ignore]
    fn test_from_unitigs() {
        let graph = graph_from_unitigs_serial::<Kmer3>(PATH_UNITIGS, STRANDED);

        println!("{:?}", graph.base.sequences);
        graph.print();
    }
}