use debruijn::{Kmer, Vmer, Exts, Dir};
use debruijn::base_to_bits;
use debruijn::graph::{BaseGraph, DebruijnGraph};

use crate::fasta_reader::FastaReader;

pub fn graph_from_unitigs_serial<K: Kmer>(path:&str, stranded: bool) -> DebruijnGraph<K, ()> {
    let mut graph: BaseGraph<K, ()> = BaseGraph::new(stranded);

    let fasta_reader = FastaReader::new(path).unwrap();

    for record in fasta_reader {
        let seq_ascii = record.sequence();
        let seq_bytes = seq_ascii.iter().map(|&b| base_to_bits(b)).collect::<Vec<u8>>();
        graph.add(seq_bytes, Exts::empty(), ());
    }
    let mut db_graph = graph.finish_serial();
    fix_exts(&mut db_graph);
    db_graph
}

pub fn graph_from_unitigs<K: Kmer+Send+Sync>(path:&str, stranded: bool) -> DebruijnGraph<K, ()> {
    let mut graph: BaseGraph<K, ()> = BaseGraph::new(stranded);

    let fasta_reader = FastaReader::new(path).unwrap();

    for record in fasta_reader {
        let seq_ascii = record.sequence();
        let seq_bytes = seq_ascii.iter().map(|&b| base_to_bits(b)).collect::<Vec<u8>>();
        graph.add(seq_bytes, Exts::empty(), ());
    }
    let mut db_graph = graph.finish();
    fix_exts(&mut db_graph);
    db_graph
}

fn fix_exts<K: Kmer>(graph: &mut DebruijnGraph<K, ()>) {
    for i in 0..graph.len() {
        let exts = get_exts(&graph, i);
        graph.base.exts[i] = exts;
    }
}

/// Checks for the presence of edges in the graph
fn get_exts<K: Kmer>(graph: &DebruijnGraph<K, ()>, node_id: usize) -> Exts {
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
    new_exts
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