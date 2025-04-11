use debruijn::{Kmer, Vmer, Exts, Dir};
use debruijn::base_to_bits;
use debruijn::graph::{BaseGraph, DebruijnGraph};

use std::sync::{Arc, Mutex};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::fasta_reader::FastaReader;

pub fn graph_from_unitigs_serial<K: Kmer>(path:&str, stranded: bool) -> DebruijnGraph<K, ()> {
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

pub fn graph_from_unitigs<K: Kmer+Send+Sync>(path:&str, stranded: bool) -> DebruijnGraph<K, ()> {
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

fn fix_exts_serial<K: Kmer>(graph: &mut DebruijnGraph<K, ()>) {
    // Get edges wich need update
    let updates: Vec<(usize, Exts)> = (0..graph.len())
        .filter_map(|i| {get_new_exts(&graph, i)})
        .collect();

    // Apply all the updates to the original graph
    for (node_id, new_exts) in updates {
        graph.base.exts[node_id] = new_exts;
    }
}

fn fix_exts<K: Kmer+Send+Sync>(graph: &mut DebruijnGraph<K, ()>) {
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
fn get_new_exts<K: Kmer>(graph: &DebruijnGraph<K, ()>, node_id: usize) -> Option<(usize, Exts)> {
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

    const PATH : &str = "../data/input/test.fna";

    #[test]
    fn test_from_fasta() {
        let graph = graph_from_unitigs_serial::<Kmer3>(PATH, true);
        graph.print();
    }
}

#[cfg(test)]
mod test_time {
    use super::*;
    use std::time::Instant;

    use debruijn::kmer;
    use debruijn::compression::ScmapCompress;

    const PATH : &str = "../data/output/chr1/graph_k31.fna";
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    #[test]
    fn time_from_unitigs() {
        let start = Instant::now();
        let graph = graph_from_unitigs_serial::<Kmer31>(PATH, false);
        let duration = start.elapsed();
        println!("Time to create graph: {:?}", duration);

        // check that the graph is compressed (ie composed of maximal unitigs)
        let compression = ScmapCompress::new();
        println!("Graph size: {}", graph.len());

        let test = graph.is_compressed(&compression);
        match test {
            None => println!("Graph is compressed"),
            Some((id1, id2)) => {
                println!("Graph is not compressed");
                let (node1, node2) = (graph.get_node(id1), graph.get_node(id2));
                println!("{:?}\t{:?}", node1.exts(), node1.sequence());
                println!("{:?}\t{:?}", node2.exts(), node2.sequence());
            }
        }
    }

    #[test]
    fn time_from_unitigs_parallel() {
        let start = Instant::now();
        let graph = graph_from_unitigs::<Kmer31>(PATH, false);
        let duration = start.elapsed();
        println!("Time to create graph (parallel): {:?}", duration);

        // check that the graph is compressed (ie composed of maximal unitigs)
        let compression = ScmapCompress::new();
        println!("Graph size: {}", graph.len());

        let test = graph.is_compressed(&compression);
        match test {
            None => println!("Graph is compressed"),
            Some((id1, id2)) => {
                println!("Graph is not compressed");
                let (node1, node2) = (graph.get_node(id1), graph.get_node(id2));
                println!("{:?}\t{:?}", node1.exts(), node1.sequence());
                println!("{:?}\t{:?}", node2.exts(), node2.sequence());
            }
        }
    }
}