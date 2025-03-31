use crate::dbg::Graph;

use debruijn::{Mer, Kmer, Exts, Dir, kmer, complement};

use std::collections::{HashSet, HashMap, VecDeque};

/// Search the shortest path between two k-mers in a graph using the Breadth-First Search algorithm.
pub fn bfs<K: Kmer>(graph: &Graph<K>, start: K, end: K) -> Vec<u8> {
    let mut queue = VecDeque::new();
    let mut added_base: Vec<u8> = vec![u8::MAX; graph.len()];
    let mut parent: Vec<usize> = vec![usize::MAX; graph.len()];

    let (start_c, flip) = start.min_rc_flip();
    let dir = Dir::Right.cond_flip(flip);
    let start_id = graph.get_key_id_unsafe(&start_c);
    let end = end.min_rc();

    added_base[start_id] = 0;
    parent[start_id] = start_id;
    queue.push_back((start_id, dir));

    // perform BFS
    while let Some((id, dir)) = queue.pop_front() {
        let kmer = graph.get_kmer(id);
        let ext_bases =  graph.get_exts(id).get(dir);
        for &base in ext_bases.iter() {
            let neigh = kmer.extend(base, dir);
            let (neigh, flip) = neigh.min_rc_flip();
            let dir_neighbor = dir.cond_flip(flip);
            let neigh_id = graph.get_key_id_unsafe(&neigh);
            // add the neighbor to the queue if it has not been visited yet
            if added_base[neigh_id] == u8::MAX {
                added_base[neigh_id] = match dir {
                    Dir::Right => base,
                    Dir::Left => complement(base),
                };
                parent[neigh_id] = id;
                queue.push_back((neigh_id, dir_neighbor));
                if neigh == end {
                    break;
                }
            }
        }
    }

    // trace back to reconstruct the path
    let mut seq: Vec<u8> = Vec::new();
    let mut id = graph.get_key_id_unsafe(&end);
    while id != start_id {
        let base = added_base[id];
        seq.push(base);
        id = parent[id];
    }
    let mut start_bases: Vec<u8> = start.iter().collect();
    start_bases.reverse();
    seq.extend(start_bases);
    seq.reverse();
    // add the start k-mer to the path
    seq

}