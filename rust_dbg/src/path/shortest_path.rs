//! Path finding algorithms for de Bruijn graphs

use super::MAX_PATH_LENGTH;
use crate::graph::Graph;
use crate::PathwayError;

use debruijn::{Dir, Kmer};

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};


// Advances the bfs by one depth and returns the new frontier set. Side indicates which end of the double-ended BFS is elongated.
// The hasmap `visited` may contain additional information about the nodes, such as their parents or distances, provided the correct closure `f` is given.
fn advance_bfs<K: Kmer, D1: Copy, D2, F>(
    graph: &Graph<K>,
    side: Dir,
    queue: &mut Vec<((usize, Dir), D1)>,
    visited: &mut HashMap<(usize, Dir), D2>,
    f: F,
) -> Vec<((usize, Dir), D1)>  where F: Fn((usize, Dir), D1)-> (D1, D2) {
    let mut next_queue = Vec::new();
    while let Some((current_node, d1)) = queue.pop() {
        for neigh_node in graph
            .get_node(current_node.0)
            .edges(current_node.1.cond_flip(side == Dir::Left)).iter()
            .map(|(neigh_id, neigh_dir, _)| (*neigh_id, neigh_dir.cond_flip(side != Dir::Left)))
        {
            let entry = visited.entry(neigh_node);
            if let Entry::Vacant(e) = entry {
                let (neigh_d1, neigh_d2) = f(current_node, d1);
                e.insert(neigh_d2);
                next_queue.push((neigh_node, neigh_d1));
            }
        }
    }
    next_queue
}

// marks the nodes as visited
fn advance_bfs_simple<K: Kmer>(
    graph: &Graph<K>,
    side: Dir,
    queue: &mut Vec<((usize, Dir), ())>,
    visited: &mut HashMap<(usize, Dir), ()>,
) -> Vec<((usize, Dir), ())> {
    advance_bfs(graph, side, queue, visited, |_, _| {
        ((), ())
    })
}

// marks the nodes as visited and stores their parents
fn advance_bfs_parents<K: Kmer>(
    graph: &Graph<K>,
    side: Dir,
    queue: &mut Vec<((usize, Dir), ())>,
    parents: &mut HashMap<(usize, Dir), (usize, Dir)>,
) -> Vec<((usize, Dir), ())> {
    advance_bfs(graph, side, queue, parents, |current_node, _| {
        ((), current_node)
    })
}

// marks the nodes as visited and stores its ancestor (the same ancestor as its parent)
fn advance_bfs_ancestor<K: Kmer>(
    graph: &Graph<K>,
    side: Dir,
    queue: &mut Vec<((usize, Dir), usize)>,
    parents: &mut HashMap<(usize, Dir), usize>,
) -> Vec<((usize, Dir), usize)> {
    advance_bfs(graph, side, queue, parents, |_, id| {
        (id, id)
    })
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
    let mut parents_left = HashMap::default();
    let mut parents_right = HashMap::default();
    let mut queue_left = Vec::new();
    let mut queue_right = Vec::new();
    parents_left.insert(start_node, start_node);
    parents_right.insert(end_node, end_node);
    queue_left.push((start_node, ()));
    queue_right.push((end_node, ()));

    // perform BFS
    let mut middle_node = None;
    'kmer: while !queue_left.is_empty() && !queue_right.is_empty() {
        if queue_left.len() <= queue_right.len() {
            // elongate from the left
            queue_left = advance_bfs_parents(graph, Dir::Left, &mut queue_left, &mut parents_left);
            // check if the two queues have met
            for (node, _) in &queue_left {
                if parents_right.contains_key(node) {
                    middle_node = Some(*node);
                    break 'kmer;
                }
            }
        } else {
            // elongate from the right
            queue_right =
                advance_bfs_parents(graph, Dir::Right, &mut queue_right, &mut parents_right);
            // check if the two queues have met
            for (node, _) in &queue_right {
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
    let mut seen = HashSet::new();
    for pos in start_position..=end_position {
        if seen.contains(&path[pos]) {
            end_position = pos - 1;
            break;
        }
        seen.insert(path[pos]);
    }

    // init BFS
    let mut visited_left = HashMap::from([(path[start_position], ())]);
    let mut queue_left = vec![(path[start_position], ())];
    let mut pos_left = start_position;

    // let mut ancestors_right: HashMap<(usize, Dir), usize> = (start_position+1..end_position).map(|pos| (path[pos], pos)).collect();
    let mut ancestors_right= HashMap::from([(path[end_position], end_position)]);
    let mut queue_right = vec![(path[end_position], end_position)];
    let mut pos_right = end_position;

    // perform BFS
    let mut shortcut = end_position + 1;
    while pos_left < pos_right
        && !queue_left.is_empty()
        && !queue_right.is_empty()
    {
        if queue_left.len() <= queue_right.len() {
            // elongate from the left
            pos_left += 1;
            queue_left = advance_bfs_simple(
                graph,
                Dir::Left,
                &mut queue_left,
                &mut visited_left,
            );
            // check if we reach the right side of the BFS...
            // ... we found a shortcut
            for (node, _) in queue_left.iter() {
                if let Some(&pos_right) = ancestors_right.get(node) {
                    shortcut = std::cmp::min(shortcut, pos_right);
                }
            }
        } else {
            // elongate from the right
            pos_right -= 1;
            queue_right = advance_bfs_ancestor(
                graph,
                Dir::Right,
                &mut queue_right,
                &mut ancestors_right,
            );
            // correct the ancestor on the path
            assert!(ancestors_right[&path[pos_right]] == pos_right+1, "something is fishy: sveral possible paths ?");
            ancestors_right.insert(path[pos_right], pos_right);

            // check if we met the left search...
            // ... we found a shortcut
            for (node, ancestor) in queue_right.iter() {
                if visited_left.contains_key(node) {
                    shortcut = std::cmp::min(shortcut, *ancestor);
                }
            }
        }
        // if shortcut: advance the end_position just before the shortcut
        if shortcut < end_position + 1 {
            end_position = shortcut - 1;
            ancestors_right = HashMap::from([(path[end_position], end_position)]);
            queue_right = vec![(path[end_position], end_position)];
            pos_right = end_position;
        }
    }
    // if the two queues have met
    if start_position + pos_left >= end_position - pos_right {
        // edge case: start and end are the same
        if start_position == end_position {
            return Ok(Some((path[start_position], 1)));
        }
        return Ok(Some((path[end_position], end_position - start_position)));
    }

    return Err(PathwayError::NoPathExists);
}
