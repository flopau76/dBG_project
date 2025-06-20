//! Path finding algorithms for de Bruijn graphs

use crate::graph::Graph;
use crate::path::MAX_PATH_LENGTH;
use crate::PathwayError;

use debruijn::{Dir, Kmer};

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::vec;

// Advances the bfs by one depth and returns the new frontier set. Side indicates which end of the double-ended BFS is elongated.
// The hasmap `visited` may contain additional information about the nodes, such as their parents or distances, provided the correct closure `f` is given.
fn advance_bfs<K: Kmer, D, F>(
    graph: &Graph<K>,
    side: Dir,
    queue: &mut Vec<(usize, Dir)>,
    visited: &mut HashMap<(usize, Dir), D>,
    f: F,
) -> Vec<(usize, Dir)>
where
    F: Fn((usize, Dir)) -> D,
{
    let mut next_queue = Vec::new();
    while let Some(current_node) = queue.pop() {
        for neigh_node in graph
            .get_node(current_node.0)
            .edges(current_node.1.cond_flip(side == Dir::Left))
            .iter()
            .map(|(neigh_id, neigh_dir, _)| (*neigh_id, neigh_dir.cond_flip(side != Dir::Left)))
        {
            let entry = visited.entry(neigh_node);
            if let Entry::Vacant(e) = entry {
                let neigh_d = f(current_node);
                e.insert(neigh_d);
                next_queue.push(neigh_node);
            }
        }
    }
    next_queue
}

// marks the nodes as visited
fn advance_bfs_simple<K: Kmer>(
    graph: &Graph<K>,
    side: Dir,
    queue: &mut Vec<(usize, Dir)>,
    visited: &mut HashMap<(usize, Dir), ()>,
) -> Vec<(usize, Dir)> {
    advance_bfs(graph, side, queue, visited, |_| ())
}

// marks the nodes as visited and stores their parents
fn advance_bfs_parents<K: Kmer>(
    graph: &Graph<K>,
    side: Dir,
    queue: &mut Vec<(usize, Dir)>,
    parents: &mut HashMap<(usize, Dir), (usize, Dir)>,
) -> Vec<(usize, Dir)> {
    advance_bfs(graph, side, queue, parents, |current_node| (current_node))
}

/// Get the shortest path between two nodes, using a double-ended BFS.
/// Note: results corresponds to ]start_node ... end_node]
pub fn get_shortest_path<K: Kmer>(
    graph: &Graph<K>,
    start_node: (usize, Dir),
    end_node: (usize, Dir),
) -> Result<Vec<(usize, Dir)>, PathwayError<K>> {
    // initialize the BFS
    let mut parents_left = HashMap::default();
    let mut parents_right = HashMap::default();
    let mut queue_left = vec![start_node];
    let mut queue_right = vec![end_node];

    // perform BFS
    let mut middle_node = None;
    'kmer: while !queue_left.is_empty() && !queue_right.is_empty() {
        if queue_left.len() <= queue_right.len() {
            // elongate from the left
            queue_left = advance_bfs_parents(graph, Dir::Left, &mut queue_left, &mut parents_left);
            // check if the two queues have met
            for node in &queue_left {
                if parents_right.contains_key(node) {
                    middle_node = Some(*node);
                    break 'kmer;
                }
                if node == &end_node {
                    middle_node = Some(*node);
                    break 'kmer;
                }
            }
        } else {
            // elongate from the right
            queue_right =
                advance_bfs_parents(graph, Dir::Right, &mut queue_right, &mut parents_right);
            // check if the two queues have met
            for node in &queue_right {
                if parents_left.contains_key(node) {
                    middle_node = Some(*node);
                    break 'kmer;
                }
                if node == &start_node {
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

    // edge case: start_node == middle_node == end_node
    if path.is_empty() {
        path.push(start_node);
    }

    Ok(path)
}

// Retuns the position-1 of the last ancestor of `node` on path[start_pos, end_pos].
fn get_ancestor_on_path(
    parents: &HashMap<(usize, Dir), (usize, Dir)>,
    path: &[(usize, Dir)],
    node: &(usize, Dir),
    start_pos: usize,
    end_pos: usize,
) -> usize {
    let mut parent = node;
    let mut pos = path[start_pos + 1..=end_pos]
        .iter()
        .position(|&x| x == *parent);
    while pos.is_none() {
        parent = parents
            .get(&parent)
            .expect("Parent not found in parents map");
        pos = path[start_pos + 1..=end_pos]
            .iter()
            .position(|&x| x == *parent);
    }
    pos.unwrap() + start_pos
}

/// Perform a BFS to find the next target_node to elongate `path` from `source_pos`.  
/// BFS is double-ended. If a shortcut is found, the right end is updated.
pub fn get_next_target_node<K: Kmer>(
    graph: &Graph<K>,
    path: &Vec<(usize, Dir)>,
    start_pos: usize,
) -> Option<((usize, Dir), usize)> {
    // no need to encode further nodes
    if start_pos >= path.len() - 1 {
        return None;
    }

    // look for duplicates in the path: a shortest path cannot pass through the same node twice (except the first node)
    let mut new_end_pos = std::cmp::min(start_pos + MAX_PATH_LENGTH, path.len() - 1);
    let mut seen = HashSet::new();
    for (pos, node) in path[start_pos..=new_end_pos].iter().enumerate() {
        if seen.contains(node) {
            new_end_pos = if node == &path[start_pos] {
                start_pos + pos
            } else {
                start_pos + pos - 1
            };
            break;
        }
        seen.insert(node);
    }

    // init left side of BFS
    let mut visited_left = HashMap::new();
    let mut queue_left = vec![path[start_pos]];
    let mut pos_left = start_pos;

    loop {
        // update right extremity and (re)start BFS
        let end_pos = new_end_pos;
        let mut parents_right = HashMap::new();
        let mut queue_right = vec![path[end_pos]];
        let mut pos_right = end_pos;

        while pos_left + 1 < pos_right && new_end_pos == end_pos {
            // advance BFS
            if queue_left.len() <= queue_right.len() {
                // elongate from the left
                pos_left += 1;
                queue_left =
                    advance_bfs_simple(graph, Dir::Left, &mut queue_left, &mut visited_left);
                // check if we found a shortcut
                // ... to the right side of the path
                for (i, node) in path[pos_left + 1..=end_pos].iter().enumerate() {
                    if visited_left.contains_key(node) {
                        new_end_pos = pos_left + i;
                        break;
                    }
                }
                // ... to the BFS from the right
                for node in &queue_left {
                    if parents_right.contains_key(node) {
                        // find the last ancestor on the path
                        let pos =
                            get_ancestor_on_path(&parents_right, path, node, start_pos, end_pos);
                        new_end_pos = std::cmp::min(new_end_pos, pos);
                    }
                }
            } else {
                // elongate from the right
                pos_right -= 1;
                queue_right =
                    advance_bfs_parents(graph, Dir::Right, &mut queue_right, &mut parents_right);
                // check if we found a shortcut
                // ... to the BFS from the left
                for node in queue_right.iter() {
                    if visited_left.contains_key(node) {
                        // find the last ancestor on the path
                        let pos = get_ancestor_on_path(
                            &parents_right,
                            path,
                            &parents_right[node],
                            start_pos,
                            end_pos,
                        );
                        new_end_pos = std::cmp::min(new_end_pos, pos);
                    }
                }
                // ... to the left side of the path
                for left_node in path[start_pos..pos_right].iter() {
                    if let Some(right_node) = parents_right.get(left_node) {
                        let pos = get_ancestor_on_path(
                            &parents_right,
                            path,
                            right_node,
                            start_pos,
                            end_pos,
                        );
                        new_end_pos = std::cmp::min(new_end_pos, pos);
                    }
                }
            }
        }
        // if we found a shortcut, we have to restart the BFS with a new target position
        if new_end_pos < end_pos {
            continue;
        }
        return Some((path[end_pos], end_pos - start_pos));
    }
}
