pub use crate::dbg::Graph;

use debruijn::{Kmer, Dir, complement, DnaBytes};

use std::collections::VecDeque;
use std::error::Error;

/// Search the shortest path between two k-mers in a graph using the Breadth-First Search algorithm.
pub fn get_shortest_path<K: Kmer>(graph: &Graph<K>, start: K, end: K) -> Result<DnaBytes, Box<dyn Error>> {
    let mut parent_canon= vec![(usize::MAX, true); graph.len()];
    let mut parent_flip= if graph.canon {vec![(usize::MAX, true); graph.len()]} else {Vec::new()};
    let mut queue = VecDeque::new();

    let (start_c, start_flip) = if graph.canon {start.min_rc_flip()} else {(start, false)};
    let (end_c, end_flip) = if graph.canon {end.min_rc_flip()} else {(end, false)};

    // start kmer
    let start_id = graph.get_key_id_unsafe(&start_c);
    let start_flip = start_flip;

    // mark the start kmer as visited and add it to the queue
    (if start_flip {&mut parent_flip} else {&mut parent_canon})[start_id] = (start_id, start_flip);
    queue.push_back((start_id, start_flip));

    // perform BFS
    'kmer: while let Some((id, flip)) = queue.pop_front() {
        let kmer_c = graph.get_kmer(id);
        // println!("   Processing kmer: {:?} at id: {:?}, flip: {:?}", kmer_c, id, flip);
        let dir = if flip { Dir::Left } else { Dir::Right };
        let ext_bases =  graph.get_exts(id).get(dir);
        'neigh: for &base in ext_bases.iter() {
            let neigh = kmer_c.extend(base, dir);
            let (neigh_c, neigh_flip) = if graph.canon {neigh.min_rc_flip()} else {(neigh, false)};
            let neigh_id = graph.get_key_id_unsafe(&neigh_c);
            let neigh_flip = flip ^ neigh_flip;
            // println!("Neighbor: {:?} at id: {:?}, flip: {:?}", neigh, neigh_id, neigh_flip);

            let neigh_parent = &mut (if neigh_flip {&mut parent_flip} else {&mut parent_canon})[neigh_id];
            if *neigh_parent == (usize::MAX, true) {    // not visited yet
                *neigh_parent = (id, flip);
                queue.push_back((neigh_id, neigh_flip));
                if (neigh_c, neigh_flip) == (end_c, end_flip) {
                    // println!("Found end, exiting BFS");
                    break 'kmer;
                }
            }
        }
    }

    if queue.is_empty() {
        return Err("No path found".into());
    }
    // println!("Parent canon: {:?}", parent_canon);
    // println!("Parent flip: {:?}", parent_flip);

    // trace back to reconstruct the path
    let mut seq: Vec<u8> = end.iter().collect();
    seq.reverse();

    let mut id = graph.get_key_id_unsafe(&end_c);
    let mut flip = end_flip;
    while (id, flip) != (start_id, flip) {
        (id, flip) = (if flip { &parent_flip } else {& parent_canon})[id];
        let kmer_c = graph.get_kmer(id);
        if flip {
            seq.push(complement(kmer_c.get(K::k() - 1)));
        } else {
            seq.push(kmer_c.get(0));
        }
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

    static PATH: &str = "../data/input/chr1/AalbF5_splitN/AalbF5_splitN.part_006.fna";

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
    fn test_from_split_fasta() {
        let canon = true;
        let graph = Graph::<Kmer31>::from_fasta(PATH, canon);
        let mut reader = FastaReader::new(PATH).unwrap();
        let mut record = reader.next().unwrap();
        let seq = record.sequence();
        let seq = DnaBytes(seq.iter().map(|&b| base_to_bits(b)).collect());


        let (start, end) = seq.both_term_kmer::<Kmer31>();
        let path = get_shortest_path(&graph, start, end).unwrap();
        let path_ascii: Vec<u8> = path.0.iter().map(|&b| bits_to_ascii(b)).collect();
        let path_str = String::from_utf8_lossy(&path_ascii);
        println!("Path: {:?}", path_str);
    }

}