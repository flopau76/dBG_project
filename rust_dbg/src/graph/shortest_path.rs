//! Path finding algorithms for de Bruijn graphs

use crate::{Graph, KmerStorage, Node, PathwayError, Side};

use std::collections::{hash_map::Entry, HashMap, HashSet};

// Advances the bfs by one depth and updates `parents` with the first encountered parent.
// `dir` indicates in which direction to advance.
// Returns the new frontier set.
fn advance_bfs_parents(
    graph: &Graph<impl KmerStorage>,
    dir: Side,
    queue: &mut Vec<Node>,
    parents: &mut HashMap<Node, Node>,
) -> Vec<Node> {
    let mut next_queue = Vec::new();
    while let Some(current_node) = queue.pop() {
        for neigh_node in graph.node_neigh(current_node, dir) {
            let entry = parents.entry(neigh_node);
            if let Entry::Vacant(e) = entry {
                e.insert(current_node);
                next_queue.push(neigh_node);
            }
        }
    }
    next_queue
}

// Advances the bfs from the given `side` and updates `parents` with all encountered parents.
// Queue is modififed in place.
fn advance_bfs_parents_all(
    graph: &Graph<impl KmerStorage>,
    side: Side,
    queue: &mut Vec<Node>,
    parents: &mut HashMap<Node, Vec<Node>>,
) -> () {
    let mut next_parents = HashMap::new();
    while let Some(current_node) = queue.pop() {
        for neigh_node in graph.node_neigh(current_node, side.opposite()) {
            if !parents.contains_key(&neigh_node) {
                let entry = next_parents.entry(neigh_node);
                match entry {
                    Entry::Vacant(e) => {
                        e.insert(vec![current_node]);
                    }
                    Entry::Occupied(mut e) => {
                        e.get_mut().push(current_node);
                    }
                }
            }
        }
    }
    let mut next_queue = Vec::new();
    for (node, parent_nodes) in next_parents {
        parents.insert(node, parent_nodes);
        next_queue.push(node);
    }
    *queue = next_queue;
    ()
}

/// Get (one of) the shortest path between two nodes, using a double-ended BFS.
/// Returns ]start_node ... end_node]
pub fn get_shortest_path<K: KmerStorage>(
    graph: &Graph<K>,
    start_node: Node,
    end_node: Node,
) -> Result<Vec<Node>, PathwayError> {
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
            queue_left =
                advance_bfs_parents(graph, Side::Right, &mut queue_left, &mut parents_left);
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
                advance_bfs_parents(graph, Side::Left, &mut queue_right, &mut parents_right);
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

    let middle_node = middle_node.expect("No path exists between the two nodes");

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

// Find the id of the first ancestor of `node` belonging to `path``, given the map `parents_right`
// If `included` is true, the node itself is considered as an ancestor.
fn first_ancestor_on_path(
    node: Node,
    parents: &HashMap<Node, Vec<Node>>,
    path: &[Node],
    included: bool,
) -> usize {
    // if node is included and on the path, return its index
    if included {
        if let Some(idx) = path.iter().position(|&n| n == node) {
            return idx;
        }
    }
    // otherwise, recursive call to find the ancestor of its parents
    parents
        .get(&node)
        .expect(format!("Node {node:?} should have parents").as_str())
        .iter()
        .map(|&n| first_ancestor_on_path(n, parents, path, true))
        .min()
        .expect("Node should be visited from the right")
}

/// Perform a (double-ended) BFS to find the next target_node to elongate `path` from `start_pos`.  
/// Returns `(target_node, dist)`, such that:  
///  - `target_node`=`path[start_pos + dist]`  
///  - `path[start_pos..=start_pos+dist]` corresponds to the shortest path between `path[start_pos]`` and `target_node`.
pub fn get_next_target_node<K: KmerStorage>(
    graph: &Graph<K>,
    path: &Vec<Node>,
    start_pos: usize,
    max_depth: usize,
) -> Option<(Node, usize)> {
    // no need to encode further nodes
    if start_pos >= path.len() - 1 {
        return None;
    }

    // look for duplicates in the path: a shortest path cannot pass through the same node twice (except the first node)
    let mut new_end_pos = std::cmp::min(start_pos + max_depth, path.len() - 1);
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

    // init BFS
    let mut parents_left: HashMap<Node, Vec<Node>> = HashMap::new();
    let mut queue_left = vec![path[start_pos]];
    let mut pos_left = start_pos;

    let mut parents_right: HashMap<Node, Vec<Node>> = HashMap::new();
    let mut queue_right = Vec::new();

    loop {
        // update right extremity and (re)start BFS
        let end_pos = new_end_pos;
        parents_right.clear();
        queue_right.clear();
        queue_right.push(path[end_pos]);
        let mut pos_right = end_pos;

        // advance BFS
        while pos_left < pos_right && new_end_pos == end_pos {
            // chose which side to elongate
            // note: when a shortcut is found, right side is reset but not left side.
            // GAMMA indicates at wich point the left side is preferred.
            const GAMMA: usize = 3;
            let side = if queue_left.len() <= GAMMA * queue_right.len() {
                pos_left += 1;
                Side::Left
            } else {
                pos_right -= 1;
                Side::Right
            };
            // advance the BFS from the given side
            {
                let queue = side.choose(&mut queue_left, &mut queue_right);
                let parents = side.choose(&mut parents_left, &mut parents_right);
                advance_bfs_parents_all(graph, side, queue, parents);
            }
            // check if the path is unique: each node on the path should only have one parent on both sides
            match side {
                Side::Left => {
                    if parents_left[&path[pos_left]].len() > 1 {
                        new_end_pos = new_end_pos.min(pos_left - 1);
                    }
                }
                Side::Right => {
                    if parents_right[&path[pos_right]].len() > 1 {
                        for parent in parents_right[&path[pos_right]].iter() {
                            if parent == &path[pos_right + 1] {
                                continue;
                            }
                            new_end_pos = new_end_pos.min(
                                first_ancestor_on_path(
                                    *parent,
                                    &parents_right,
                                    &path[pos_left + 1..=end_pos],
                                    true,
                                ) + pos_left,
                            );
                        }
                    }
                }
            }
            // check if the two BFS have met
            let queue = side.choose(&queue_left, &queue_right);
            let parents_opposite = side.choose(&parents_right, &parents_left);
            if pos_left == pos_right {
                // the two BFS should have met at a single node, belonging to the path
                for node in queue.iter() {
                    // we found a node that is in the other BFS
                    if node != &path[pos_left] && parents_opposite.contains_key(node) {
                        new_end_pos = new_end_pos.min(
                            first_ancestor_on_path(
                                *node,
                                &parents_right,
                                &path[pos_left + 1..=end_pos],
                                true,
                            ) + pos_left,
                        );
                    }
                }
            } else {
                // the two BFS should not have met at all
                for node in queue.iter() {
                    // we found a node that is in the other BFS
                    if parents_opposite.contains_key(node) {
                        new_end_pos = new_end_pos.min(
                            first_ancestor_on_path(
                                *node,
                                &parents_right,
                                &path[pos_left + 1..=end_pos],
                                true,
                            ) + pos_left,
                        );
                    }
                    // we found a node between the two BFS, which should not be reached yet
                    if let Some(node) = path[pos_left + 1..pos_right].iter().find(|&&n| n == *node)
                    {
                        new_end_pos = new_end_pos.min(
                            first_ancestor_on_path(
                                *node,
                                &parents_right,
                                &path[pos_left + 1..=end_pos],
                                side == Side::Left,
                            ) + pos_left,
                        );
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

/// Perform a (single-ended) BFS to find the next target_node to elongate `path` from `start_pos`.  
/// Returns `(target_node, dist)`, such that:  
///  - `target_node`=`path[start_pos + dist]`  
///  - `path[start_pos..=start_pos+dist]` corresponds to the shortest path between `path[start_pos]`` and `target_node`.
pub fn get_next_target_node_naive<K: KmerStorage>(
    graph: &Graph<K>,
    path: &Vec<Node>,
    start_pos: usize,
    max_depth: usize,
) -> Option<(Node, usize)> {
    // no need to encode further nodes
    if start_pos >= path.len() - 1 {
        return None;
    }

    // look for duplicates in the path: a shortest path cannot pass through the same node twice (except the first node)
    let mut new_end_pos = std::cmp::min(start_pos + max_depth, path.len() - 1);
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

    // init BFS
    let mut parents: HashMap<Node, Vec<Node>> = HashMap::new();
    let mut queue = vec![path[start_pos]];
    let mut current_pos = start_pos;

    while current_pos < new_end_pos {
        // advance BFS
        current_pos += 1;
        advance_bfs_parents_all(graph, Side::Left, &mut queue, &mut parents);
        // check if multiple parents exist for the current node
        if parents[&path[current_pos]].len() > 1 {
            new_end_pos = current_pos - 1;
            break;
        }
        // check if we found a shortcut
        for pos in current_pos + 1..=new_end_pos {
            if parents.contains_key(&path[pos]) {
                new_end_pos = pos - 1;
                break;
            }
        }
    }
    return Some((path[new_end_pos], new_end_pos - start_pos));
}
