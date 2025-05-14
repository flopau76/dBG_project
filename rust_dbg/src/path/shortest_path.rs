//! Path finding algorithms for de Bruijn graphs

use crate::PathwayError;
use crate::graph::Graph;
use super::Extension;

use debruijn::{Dir, Kmer};

use ahash::{AHashMap, AHashSet};
use std::collections::hash_map::Entry;

pub const MAX_LENGTH: usize = 35;
pub const MIN_LENGTH: usize = 0;

/// A struct that represents a shortest path extension by adding the shortest path towards a target node.
#[derive(Debug, Copy, Clone)]
pub struct ShortestPath {
    pub target_node: (usize, Dir),
}

impl Extension for ShortestPath {
    fn extend_path<K: Kmer>(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>) {
        let start_node = *path.last().unwrap();
        let end_node = self.target_node;
        let extension = get_shortest_path(graph, start_node, end_node).unwrap();
        path.extend(extension);
    }
    fn get_next_ext<K: Kmer>(graph: &Graph<K>, path: &Vec<(usize, Dir)>, position: usize) -> Option<(Self, usize)> {
        // get_next_target_node(graph, path, position).unwrap()
        get_next_target_node_outer(graph, path, position).unwrap()
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
/// Note: results corresponds to (start_node ... end_node]
fn get_shortest_path<K: Kmer>(graph: &Graph<K>, start_node: (usize, Dir), end_node: (usize, Dir)) -> Result<Vec<(usize, Dir)>, PathwayError<K>> {

    // edge case: start and end are the same
    if start_node == end_node {
        return Ok(vec![end_node]);
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
        if queue_left.len() <= queue_right.len() {
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

    while node != start_node {
        path.push(node);
        node = *parents_left.get(&node).expect("Node not found in parents map");
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

// Perform a BFS to find the end of the longest path elongating from the start_position given path.  
// Return the next target node and the number of nodes in the path.
pub fn get_next_target_node<K: Kmer>(graph: &Graph<K>, path: &Vec<(usize, Dir)>, start_position: usize) -> Result<Option<(ShortestPath, usize)>, PathwayError<K>> {
    let mut node_iter = path.iter();
    let start_node = *node_iter.nth(start_position).unwrap();
    let mut current_node = start_node;
    let mut next_node = match node_iter.next() {
        Some(node) => node,
        None => return Ok(None)
    };

    // initialize the BFS
    let mut parents= AHashMap::default();
    let mut queue = Vec::new();
    queue.push(current_node);

    // perform the BFS
    let mut depth: usize = 0;
    while !parents.contains_key(next_node) && depth < MAX_LENGTH {
        depth += 1;
        // explore the next level of the BFS
        queue = advance_bfs(graph, &mut queue, &mut parents, Dir::Left);    // TODO: handle case with several paths of shortest length
        assert!(parents.contains_key(&next_node));

        // advance the kmer iterator to the next node
        current_node = *next_node;
        next_node = match node_iter.next() {
            Some(node) => node,
            None => break
        };
    }
    // println!("{} unitigs:\t {:?}\t {:?}", depth, start_node, current_node);
    Ok(Some((ShortestPath{target_node: current_node}, depth)))
}

// Recursive and double ended
pub fn get_next_target_node_outer<K: Kmer>(graph: &Graph<K>, path: &Vec<(usize, Dir)>, start_position: usize) -> Result<Option<(ShortestPath, usize)>, PathwayError<K>> {

    if start_position >= path.len()-1 {
        return Ok(None);
    }

    let mut end_position = std::cmp::min(start_position + 2*MAX_LENGTH, path.len()-1);
    // look for duplicates in the path: a shortest path cannot pass through a node twice
    let mut seen = AHashSet::new();
    for pos in start_position..end_position {
        if seen.contains(&path[pos]) {
            end_position = pos-1;
            break;
        }
        seen.insert(path[pos]);
    }

    let parents_left = AHashMap::from([(path[start_position], 0)]);
    let queue_left = vec!(path[start_position]);
    let depth_left = 0;
    get_next_target_node_inner(graph, path, start_position, end_position, parents_left, queue_left, depth_left)
}

// Advances the bfs by one depth and returns the new frontier set. Side indicates which end of the double-ended BFS is elongated.
fn advance_bfs_depth<K: Kmer>(graph: &Graph<K>, queue: &mut Vec<(usize, Dir)>, parents: &mut AHashMap<(usize, Dir), usize>, side:Dir, depth: usize) -> Vec<(usize, Dir)> {
    let mut next_queue = Vec::new();
    while let Some(current_node) = queue.pop() {
        for (neigh_id, neigh_dir, _) in graph.get_node(current_node.0).edges(current_node.1.cond_flip(side==Dir::Left)) {
            let entry = parents.entry((neigh_id, neigh_dir.cond_flip(side!=Dir::Left)));
            if let Entry::Vacant(e) = entry {
                e.insert(depth);
                next_queue.push((neigh_id, neigh_dir.cond_flip(side!=Dir::Left)));
            }
        }
    }
    next_queue
}

// Recursive and double ended
fn get_next_target_node_inner<K: Kmer> (
        graph: &Graph<K>,
        path: &Vec<(usize, Dir)>,
        start_position: usize,
        end_position: usize,
        mut parents_left: AHashMap<(usize, Dir), usize>,
        mut queue_left: Vec<(usize, Dir)>,
        mut depth_left: usize
    ) -> Result<Option<(ShortestPath, usize)>, PathwayError<K>> {

    // edge case: start and end are the same
    if start_position == end_position {
        return Ok(Some((ShortestPath{target_node: path[start_position]}, 1)));
    }

    // init BFS
    let mut parents_right = AHashMap::from([(path[end_position], 0)]);
    let mut queue_right = vec!(path[end_position]);
    let mut depth_right = 0;

    let mut shortcut = None;
    for pos in start_position + 1..=end_position {
        if let Some(depth) = parents_left.get(&path[pos]) {
            if start_position+*depth < pos {
                // panic!("Found a shortcut: {:?} -> {:?}", path[start_position+depth], path[pos]);
                shortcut = Some(pos);
                break;
            }
        }
    }

    while shortcut.is_none() && start_position+depth_left < end_position-depth_right {
        // advance the BFS
        if queue_left.len() <= queue_right.len() {
            // elongate from the left
            depth_left += 1;
            assert!(!queue_left.contains(&path[start_position+depth_left]));
            queue_left = advance_bfs_depth(graph, &mut queue_left, &mut parents_left, Dir::Left, depth_left);
            assert!(queue_left.contains(&path[start_position+depth_left]));
            // check if we found a shortcut
            shortcut = (start_position + depth_left + 1..=end_position)
                .find(|&pos| queue_left.contains(&path[pos]));
        } else {
            // elongate from the right
            depth_right += 1;
            assert!(!queue_right.contains(&path[end_position-depth_right]));
            queue_right = advance_bfs_depth(graph, &mut queue_right, &mut parents_right, Dir::Right, depth_right);
            assert!(queue_right.contains(&path[end_position-depth_right]));
            // check if we found a shortcut
            for pos in (start_position..=end_position-depth_right-1).rev() {
                if let Some(depth) = parents_right.get(&path[pos]) {
                    shortcut = Some(end_position - depth);
                    break;
                }
            }
        }
    }
    // if shortcut: recursive call to (start_pos, shortcut-1)
    if let Some(shortcut) = shortcut {
        let target = get_next_target_node_inner(graph, path, start_position, shortcut-1, parents_left, queue_left, depth_left )?;
        return Ok(target)
    }
    return Ok(Some((ShortestPath{target_node: path[end_position]}, end_position-start_position)))

}