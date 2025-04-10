//! Path finding algorithms for de Bruijn graphs
use crate::graph_old::Graph;
use crate::fasta_reader::DnaRecord;

use debruijn::{Dir, DnaBytes, Kmer};

use ahash::{AHashMap, AHashSet};
use std::collections::VecDeque;
use std::collections::hash_map::Entry::Vacant;
use std::error::Error;

/// Search the shortest path between two k-mers in a graph using the Breadth-First Search algorithm.
pub fn get_shortest_path<K: Kmer>(graph: &Graph<K>, start: K, end: K) -> Result<DnaBytes, Box<dyn Error>> {
    let mut parents= AHashMap::default();
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

/// Breaks the input haplotype into segments corresponding to shortest paths in the graph.
pub fn get_checkpoints<K: Kmer>(graph: &Graph<K>, haplo: &DnaRecord) -> Vec<(K, K)> {
    let mut chunks = Vec::new();
    let mut kmer_iter = haplo.iter_kmers::<K>(false);

    let mut start = kmer_iter.next();
    let mut nb_chunk = 0;

    // while not at the end of the haplotype: start a new BFS
    while start.is_some() {
        nb_chunk += 1;
        let mut size_chunk = 0;

        let start_kmer = start.unwrap();
        let mut FS = Vec::new();        // frontier set
        let mut NS = vec![start_kmer];  // next set
        let mut visited = AHashSet::default();

        let mut current_kmer = start_kmer;
        let mut next_kmer = kmer_iter.next();

        // while BFS coincides with the haplotype: explore next level of BFS
        while next_kmer.is_some() && !visited.contains(&next_kmer.unwrap()) {
            size_chunk += 1;

            FS = NS;
            NS = Vec::new();
            current_kmer = next_kmer.unwrap();
            next_kmer = kmer_iter.next();

            // while there are kmers in the frontier set: explore their neighbors
            while let Some(kmer) = FS.pop() {
                let exts = graph.get_exts(&kmer).expect("kmer was added to set so must be present");
                for &base in exts.get(Dir::Right).iter() {
                    let neigh = kmer.extend_right(base);
                    let new = visited.insert(neigh);
                    if new {
                        NS.push(neigh);
                    }
                }
            }
        }

        println!("Found chunk {} after {} steps", nb_chunk, size_chunk);
        chunks.push((start_kmer, current_kmer));
        start = next_kmer;
    }
    chunks
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use debruijn::kmer::Kmer3;
    use debruijn::{base_to_bits, bits_to_ascii, DnaBytes, Kmer, Vmer};

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

    #[test]
    fn test_get_checkpoints() {
        let record = DnaRecord::new(String::from("test"), SEQ.to_vec());
        for &canon in [false, true].iter() {
            let graph = make_graph(canon);
            let path = get_checkpoints(&graph, &record);
            println!("Path: {:?}", path);
        }
    }
}