use crate::dbg::Graph;

use debruijn::{Kmer, Exts, Dir, kmer};

use std::collections::{HashSet, HashMap, VecDeque};

/// Search the shortest path between two k-mers in a graph using the Breadth-First Search algorithm.
fn bfs<K: Kmer>(graph: &Graph<K>, start: K, end: K) {
    let mut queue = VecDeque::new();
    let mut parent: Vec<usize> = vec![usize::MAX; graph.len()];

    let (start, flip) = start.min_rc_flip();
    let dir = Dir::Right.cond_flip(flip);
    let start_id = graph.get_key_id_unsafe(&start);
    let end = end.min_rc();

    parent[start_id] = start_id;
    queue.push_back((start_id, dir));

    // perform BFS
    while let Some((id, dir)) = queue.pop_front() {
        let kmer = graph.get_kmer(id);
        let exts = graph.get_exts(id);

        for neigh in kmer.get_extensions(*exts, dir) {
            let (neigh, flip) = neigh.min_rc_flip();
            let neigh_id = graph.get_key_id_unsafe(&neigh);
            if neigh == end {
                parent[neigh_id] = id;
                break;
            }
            if parent[neigh_id] == usize::MAX {
                parent[neigh_id] = id;
                queue.push_back((neigh_id, dir.cond_flip(flip)));
            }
        }
    }

    // trace back to reconstruct the path

}