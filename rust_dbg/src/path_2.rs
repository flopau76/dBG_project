pub use crate::graph_2::Graph;

use debruijn::{Kmer, Dir, complement, DnaBytes};

use std::collections::{HashMap, VecDeque};
use std::collections::hash_map::Entry::Vacant;
use std::error::Error;

/// Search the shortest path between two k-mers in a graph using the Breadth-First Search algorithm.
pub fn get_shortest_path<K: Kmer>(graph: &Graph<K>, start: K, end: K) -> Result<DnaBytes, Box<dyn Error>> {
    let mut parents= HashMap::new();
    let mut queue = VecDeque::new();

    // edge case: start and end are the same
    if start == end {
        return Ok(DnaBytes(start.iter().collect()));
    }

    // mark the start kmer as visited and add it to the queue
    parents.insert(start, start);
    queue.push_back(start);

    // perform BFS
    'kmer: while let Some(kmer) = queue.pop_front() {
        // println!("   Processing kmer: {:?}", kmer);
        let exts = graph.get_exts(&kmer).expect("kmer was added to queue so should be present");
        for &base in exts.get(Dir::Right).iter() {
            let neigh = kmer.extend_right(base);
            // println!("Neighbor: {:?}", neigh);
            let entry = parents.entry(neigh);
            if let Vacant(e) = entry {
                e.insert(kmer);
                queue.push_back(neigh);
            }
            if neigh == end {
                // println!("Found end, exiting BFS");
                break 'kmer;
            }
        }
        // println!("Queue: {:?}", queue);
        // println!("Parents: {:?}", parents);
    }

    if queue.is_empty() {
        return Err("No path found".into());
    }
    // println!("Parents: {:?}", parents);

    // trace back to reconstruct the path
    let mut seq: Vec<u8> = end.iter().collect();
    seq.reverse();

    let mut kmer = end;
    while kmer != start {
        kmer = *parents.get(&kmer).expect("parent should be in graph");
        seq.push(kmer.get(0));
    }
    seq.reverse();
    Ok(DnaBytes(seq))
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use crate::fasta_reader::FastaReader;

    use debruijn::kmer::{Kmer3, K31, VarIntKmer};
    use debruijn::{base_to_bits, bits_to_ascii, DnaBytes, Kmer, Vmer};

    type Kmer31 = VarIntKmer<u64, K31>;

    static SEQ: &[u8] = b"GgctGAGCTGAGTT";
    static SHORTEST_NO_C: &[u8] = b"GGCTGAGTT";
    static SHORTEST_C_1: &[u8] = b"GGCTGAGTT";
    static SHORTEST_C_2: &[u8] = b"GGCTCAGTT";

    fn make_seq() -> DnaBytes {
        DnaBytes(SEQ.iter().map(|&b| base_to_bits(b)).collect())
    }

    fn make_graph(canon: bool) -> Graph<Kmer3> {
        let seq = make_seq();
        let kmers: Vec<Kmer3> = if canon { seq.iter_kmers::<Kmer3>().map(|k| k.min_rc()).collect() } else { seq.iter_kmers().collect() };
        Graph::from_kmers(kmers, canon)
    }

    #[test]
    fn test_shortest_path() {
        let seq = make_seq();
        let (start, end) = seq.both_term_kmer::<Kmer3>();
        for &canon in [false, true].iter() {
            let graph = make_graph(canon);
            let path = get_shortest_path(&graph, start, end).unwrap();
            let path_ascii: Vec<u8> = path.0.iter().map(|&b| bits_to_ascii(b)).collect();
            if canon {
                assert!(path_ascii == SHORTEST_C_1 || path_ascii == SHORTEST_C_2, "Path: {:?}, Expected: {:?} or {:?}", path_ascii, SHORTEST_C_1, SHORTEST_C_2);
            } else {
                assert_eq!(path_ascii, SHORTEST_NO_C, "Path: {:?}, Expected: {:?}", path_ascii, SHORTEST_NO_C);
            };
        }
    }

}