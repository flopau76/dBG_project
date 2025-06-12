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
) -> Vec<((usize, Dir), D1)>
where
    F: Fn((usize, Dir), D1) -> (D1, D2),
{
    let mut next_queue = Vec::new();
    while let Some((current_node, d1)) = queue.pop() {
        for neigh_node in graph
            .get_node(current_node.0)
            .edges(current_node.1.cond_flip(side == Dir::Left))
            .iter()
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
    advance_bfs(graph, side, queue, visited, |_, _| ((), ()))
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

fn get_ancestor_on_path(
    parents: &HashMap<(usize, Dir), (usize, Dir)>,
    path: &Vec<(usize, Dir)>,
    node: &(usize, Dir),
) -> usize {
    let mut parent = node;
    let mut pos = path.iter().position(|&x| x == *parent);
    while pos.is_none() {
        parent = parents
            .get(&parent)
            .expect("Parent not found in parents map");
        pos = path.iter().position(|&x| x == *parent);
    }
    pos.unwrap()
}

/// Perform a BFS to find the next target_node to elongate `path` from `source_pos`.  
/// BFS is double-ended. If a shortcut is found, the right en is updated.
pub fn get_next_target_node<K: Kmer>(
    graph: &Graph<K>,
    path: &Vec<(usize, Dir)>,
    source_pos: usize,
) -> Option<((usize, Dir), usize)> {
    if source_pos >= path.len() - 1 {
        return None;
    }

    // look for duplicates in the path: a shortest path cannot pass through the same node twice
    let mut shortcut = std::cmp::min(source_pos + MAX_PATH_LENGTH + 1, path.len());
    let mut seen = HashSet::new();
    for pos in source_pos..=shortcut - 1 {
        if seen.contains(&path[pos]) {
            shortcut = pos;
            break;
        }
        seen.insert(path[pos]);
    }

    // init left side of BFS
    let mut visited_left = HashMap::from([(path[source_pos], ())]);
    let mut queue_left = vec![(path[source_pos], ())];
    let mut pos_left = source_pos;

    loop {
        // update right extremity and (re)start BFS
        let target_pos = shortcut - 1;
        let mut parents_right = HashMap::from([(path[target_pos], path[target_pos])]);
        let mut queue_right = vec![(path[target_pos], ())];
        let mut pos_right = target_pos;

        while pos_left < pos_right - 1 && shortcut == target_pos + 1 {
            // advance BFS
            if queue_left.len() <= queue_right.len() {
                // elongate from the left
                pos_left += 1;
                queue_left =
                    advance_bfs_simple(graph, Dir::Left, &mut queue_left, &mut visited_left);
                // check if the two queues have met
                for (node, _) in &queue_left {
                    if parents_right.contains_key(node) {
                        // find the last ancestor on the path
                        let pos = get_ancestor_on_path(&parents_right, path, node);
                        shortcut = std::cmp::min(shortcut, pos);
                    }
                }
            } else {
                // elongate from the right
                pos_right -= 1;
                if pos_left == pos_right {
                    break;
                }
                queue_right =
                    advance_bfs_parents(graph, Dir::Right, &mut queue_right, &mut parents_right);
                // check if the two queues have met
                for (node, _) in queue_right.iter() {
                    if visited_left.contains_key(node) {
                        // find the last ancestor on the path
                        let pos = get_ancestor_on_path(&parents_right, path, node);
                        shortcut = std::cmp::min(shortcut, pos);
                    }
                }
            }
        }
        // if we found a shortcut, we have to restart the BFS with a new target position
        if shortcut <= target_pos {
            continue;
        }
        // otherwise both ends met and we can return the target node
        if source_pos == target_pos {
            // edge case: start and end are the same
            return Some((path[source_pos], 1));
        }
        return Some((path[target_pos], target_pos - source_pos));
    }
}
