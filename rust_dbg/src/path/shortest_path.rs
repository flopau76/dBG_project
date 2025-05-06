//! Path finding algorithms for de Bruijn graphs

use crate::PathwayError;
use crate::graph::Graph;
use super::{node_iterator::NodeIterator, Extension};

use debruijn::{Dir, Vmer, Kmer};

use ahash::AHashMap;
use std::collections::hash_map::Entry;

/// A struct that represents a shortest path extension by adding the shortest path towards a target node.
pub struct ShortestPath {
    target_node: (usize, Dir),
}

impl<K: Kmer> Extension<K> for ShortestPath {
    fn extend_path(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>) {
        let start_node = *path.last().unwrap();
        let end_node = self.target_node;
        let extension = get_shortest_path(graph, start_node, end_node).unwrap();
        path.extend(extension);
    }

    fn next(graph: &Graph<K>, path: &Vec<(usize, Dir)>, position: usize) -> Self {
        todo!()
    }
}

// Advances the bfs by one depth and returns the new frontier set. Side indicates which end of the double-ended BFS is elongated.
fn advance_bfs<K: Kmer>(graph: &Graph<K>, queue: &mut Vec<(usize, Dir)>, parents: &mut AHashMap<(usize, Dir), (usize, Dir)>, side:Dir) -> Vec<(usize, Dir)> {
    let mut next_queue = Vec::new();
    while let Some(current_node) = queue.pop() {
        for (neigh_id, neigh_dir, _) in graph.get_node(current_node.0).edges(current_node.1.cond_flip(side==Dir::Left)) {
            let entry = parents.entry((neigh_id, neigh_dir.cond_flip(side!=Dir::Left)));
            if let Entry::Vacant(e) = entry {
                e.insert(current_node);
                next_queue.push((neigh_id, neigh_dir.cond_flip(side!=Dir::Left)));
            }
        }
    }
    next_queue
}

/// Get the shortest path between two nodes, using a double-ended BFS.
fn get_shortest_path<K: Kmer>(graph: &Graph<K>, start_node: (usize, Dir), end_node: (usize, Dir)) -> Result<Vec<(usize, Dir)>, PathwayError<K>> {

    // edge case: start and end are the same
    if start_node == end_node {
        return Ok(vec![start_node]);
    }

    // initialize the BFS
    let mut parents_left= AHashMap::default();
    let mut parents_right= AHashMap::default();
    let mut queue_left = Vec::new();
    let mut queue_right = Vec::new();
    parents_left.insert(start_node, start_node);
    parents_right.insert(end_node, end_node);
    queue_left.push(start_node);
    queue_right.push(end_node);

    // perform BFS
    let mut middle_node = None;
    'kmer: while !queue_left.is_empty() && !queue_right.is_empty() {
        if queue_left.len() < queue_right.len() {
            // elongate from the left
            queue_left = advance_bfs(graph, &mut queue_left, &mut parents_left, Dir::Left);
            // check if the two queues have met
            for node in &queue_left {
                if parents_right.contains_key(node) {
                    middle_node = Some(*node);
                    break 'kmer;
                }
            }
        } else {
            // elongate from the right
            queue_right = advance_bfs(graph, &mut queue_right, &mut parents_right, Dir::Right);
            // check if the two queues have met
            for node in &queue_right {
                if parents_left.contains_key(node) {
                    middle_node = Some(*node);
                    break 'kmer;
                }
            }
        }
    }

    if middle_node.is_none() {
        return Err(PathwayError::NoPathExists);
    }
    let middle_node = middle_node.unwrap();

    // traceback left side by baktracking from the middle node
    let mut node = middle_node;
    let mut path = Vec::new();
    path.push(node);

    while node != start_node {
        node = *parents_left.get(&node).expect("Node not found in parents map");
        path.push(node);
    }
    path.reverse();

    // traceback right side by baktracking from the middle node
    node = middle_node;
    while node != end_node {
        node = *parents_right.get(&node).expect("Node not found in parents map");
        path.push(node);
    }

    Ok(path)
}

// Perform a BFS to find the end of the longest path elongating allong the unitig iterator
// Before, the unitig iterator must peek at the node following start_node
// After, the unitig iterator peeks at the node following end_node
pub fn get_next_target_node<K: Kmer, D: Vmer>(graph: &Graph<K>, unitig_iter: &mut NodeIterator<K,D>, start_node: (usize, Dir)) -> Result<Option<ShortestPath>, PathwayError<K>> {
    let mut current_node = start_node;
    let mut next_node = match unitig_iter.peek() {
        Some(node) => node,
        None => {
            println!("1 unitigs:\t {:?}\t {:?}", start_node, start_node);
            return Ok(Some(ShortestPath{target_node: start_node}))
        },
    };

    // initialize the BFS
    let mut parents= AHashMap::default();
    let mut queue = Vec::new();
    parents.insert(start_node, start_node);
    queue.push(start_node);

    // perform the BFS
    let mut depth: usize = 1;
    while !parents.contains_key(&next_node) {
        // explore the next level of the BFS
        queue = advance_bfs(graph, &mut queue, &mut parents, Dir::Left);    // TODO: handle case with several paths of shortest length
        assert!(parents.contains_key(&next_node));

        // advance the kmer iterator to the next node
        current_node = next_node;
        unitig_iter.advance()?;
        next_node = match unitig_iter.peek() {
            Some(node) => node,
            None => break
        };
        depth += 1;
    }
    println!("{} unitigs:\t {:?}\t {:?}", depth, start_node, current_node);
    Ok(Some(ShortestPath{target_node: current_node}))
}