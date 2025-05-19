//! Path finding algorithms for de Bruijn graphs

use super::MAX_PATH_LENGTH;
use crate::PathwayError;
use crate::graph::Graph;

use debruijn::{Dir, Kmer};

use ahash::{AHashMap, AHashSet};
use std::collections::hash_map::Entry;

// Advances the bfs by one depth and returns the new frontier set. Side indicates which end of the double-ended BFS is elongated.
fn advance_bfs_parents<K: Kmer>(
    graph: &Graph<K>,
    queue: &mut Vec<(usize, Dir)>,
    parents: &mut AHashMap<(usize, Dir), (usize, Dir)>,
    side: Dir,
) -> Vec<(usize, Dir)> {
    let mut next_queue = Vec::new();
    while let Some(current_node) = queue.pop() {
        for (neigh_id, neigh_dir, _) in graph
            .get_node(current_node.0)
            .edges(current_node.1.cond_flip(side == Dir::Left))
        {
            let entry = parents.entry((neigh_id, neigh_dir.cond_flip(side != Dir::Left)));
            if let Entry::Vacant(e) = entry {
                e.insert(current_node);
                next_queue.push((neigh_id, neigh_dir.cond_flip(side != Dir::Left)));
            }
        }
    }
    next_queue
}

/// Get the shortest path between two nodes, using a double-ended BFS.
/// Note: results corresponds to (start_node ... end_node]
pub fn get_shortest_path<K: Kmer>(
    graph: &Graph<K>,
    start_node: (usize, Dir),
    end_node: (usize, Dir),
) -> Result<Vec<(usize, Dir)>, PathwayError<K>> {
    // edge case: start and end are the same
    if start_node == end_node {
        return Ok(vec![end_node]);
    }

    // initialize the BFS
    let mut parents_left = AHashMap::default();
    let mut parents_right = AHashMap::default();
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
            queue_left = advance_bfs_parents(graph, &mut queue_left, &mut parents_left, Dir::Left);
            // check if the two queues have met
            for node in &queue_left {
                if parents_right.contains_key(node) {
                    middle_node = Some(*node);
                    break 'kmer;
                }
            }
        } else {
            // elongate from the right
            queue_right =
                advance_bfs_parents(graph, &mut queue_right, &mut parents_right, Dir::Right);
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
        node = *parents_left
            .get(&node)
            .expect("Node not found in parents map");
    }
    path.reverse();

    // traceback right side by baktracking from the middle node
    node = middle_node;
    while node != end_node {
        node = *parents_right
            .get(&node)
            .expect("Node not found in parents map");
        path.push(node);
    }

    Ok(path)
}

// Advances the bfs by one level and returns the new frontier set. Side indicates which end of the double-ended BFS is elongated.
fn advance_bfs_distances<K: Kmer>(
    graph: &Graph<K>,
    queue: &mut Vec<(usize, Dir)>,
    distances: &mut AHashMap<(usize, Dir), usize>,
    side: Dir,
    current_dist: usize,
) -> Vec<(usize, Dir)> {
    let mut next_queue = Vec::new();
    while let Some(current_node) = queue.pop() {
        for (neigh_id, neigh_dir, _) in graph
            .get_node(current_node.0)
            .edges(current_node.1.cond_flip(side == Dir::Left))
        {
            if let Entry::Vacant(e) =
                distances.entry((neigh_id, neigh_dir.cond_flip(side != Dir::Left)))
            {
                e.insert(current_dist);
                next_queue.push((neigh_id, neigh_dir.cond_flip(side != Dir::Left)));
            }
        }
    }
    next_queue
}

/// Perform a BFS to find the next target_node to elongate `path` from `start_position`.  
/// BFS is double-ended. If a shortcut is found, the right en is updated.
pub fn get_next_target_node<K: Kmer>(
    graph: &Graph<K>,
    path: &Vec<(usize, Dir)>,
    start_position: usize,
) -> Result<Option<((usize, Dir), usize)>, PathwayError<K>> {
    if start_position >= path.len() - 1 {
        return Ok(None);
    }

    // look for duplicates in the path: a shortest path cannot pass through the same node twice
    let mut end_position = std::cmp::min(start_position + MAX_PATH_LENGTH, path.len() - 1);
    let mut seen = AHashSet::new();
    for pos in start_position..=end_position {
        if seen.contains(&path[pos]) {
            end_position = pos - 1;
            break;
        }
        seen.insert(path[pos]);
    }

    // init BFS
    let mut distances_left = AHashMap::from([(path[start_position], 0)]);
    let mut queue_left = vec![path[start_position]];
    let mut dist_left = 0;
    let mut distances_right = AHashMap::from([(path[end_position], 0)]);
    let mut queue_right = vec![path[end_position]];
    let mut dist_right = 0;

    // perform BFS
    let mut shortcut = end_position + 1;
    while start_position + dist_left < end_position - dist_right
        && !queue_left.is_empty()
        && !queue_right.is_empty()
    {
        if queue_left.len() <= queue_right.len() {
            // elongate from the left
            dist_left += 1;
            queue_left = advance_bfs_distances(
                graph,
                &mut queue_left,
                &mut distances_left,
                Dir::Left,
                dist_left,
            );
            // check if we found a shortcut...
            // ...to a node visited from the right
            for node in queue_left.iter() {
                if let Some(&depth_right) = distances_right.get(node) {
                    if dist_left + depth_right < end_position - start_position {
                        shortcut = std::cmp::min(shortcut, end_position - depth_right);
                    }
                }
            }
            // ...to a node on the path
            for pos in start_position + dist_left + 1..=end_position {
                if queue_left.contains(&path[pos]) {
                    shortcut = std::cmp::min(shortcut, pos);
                    break;
                }
            }
        } else {
            // elongate from the right
            dist_right += 1;
            queue_right = advance_bfs_distances(
                graph,
                &mut queue_right,
                &mut distances_right,
                Dir::Right,
                dist_right,
            );
            // check if we found a shortcut...
            // ...to a node visited from the left
            for node in queue_right.iter() {
                if let Some(&depth_left) = distances_left.get(node) {
                    if depth_left + dist_right < end_position - start_position {
                        shortcut = std::cmp::min(shortcut, end_position - dist_right);
                    }
                }
            }
            // ...to a node on the path
            for pos in (start_position..=end_position - dist_right - 1).rev() {
                if queue_right.contains(&path[pos]) {
                    shortcut = std::cmp::min(shortcut, pos);
                    break;
                }
            }
        }
        // if shortcut: advance the end_position just before the shortcut
        if shortcut < end_position + 1 {
            end_position = shortcut - 1;
            distances_right = AHashMap::from([(path[end_position], 0)]);
            queue_right = vec![path[end_position]];
            dist_right = 0;
        }
    }
    // if the two queues have met
    if start_position + dist_left >= end_position - dist_right {
        // edge case: start and end are the same
        if start_position == end_position {
            return Ok(Some((path[start_position], 1)));
        }
        return Ok(Some((path[end_position], end_position - start_position)));
    }

    return Err(PathwayError::NoPathExists);
}
