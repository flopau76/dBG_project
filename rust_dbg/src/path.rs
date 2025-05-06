//! Defines how to encode a path in a debruijn Graph. Pri
 
use debruijn::{dna_string::DnaString, Dir, Kmer};

use crate::graph::Graph;
use self::{node_iterator::NodeIterator, shortest_path::ShortestPath};

use std::fmt::Debug;

pub mod node_iterator;
pub mod shortest_path;

/// Some information about how to extend a path in a given graph
trait Extension<K: Kmer>: Debug {
    /// Get the next extension and the number of incoded nodes, given the full path and the current node in the path
    fn next(graph: &Graph<K>, path: &Vec<(usize, Dir)>, position: usize) -> Option<(Self, usize)> where Self: Sized;
    /// Extend the current path given the extension
    fn extend_path(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>);
}

/// The easiest way to extend the current path is to give the following node
impl<K: Kmer> Extension<K> for (usize, Dir) {
    fn next(_graph: &Graph<K>, path: &Vec<(usize, Dir)>, position: usize) -> Option<(Self, usize)> {
        path.get(position+1).map(|node| (*node, 1))
    }
    fn extend_path(&self, _graph: &Graph<K>, path: &mut Vec<(usize, Dir)>) {
        path.push(*self);
    }
}

enum MyExtension {
    ShortestPath(shortest_path::ShortestPath),
    NextNode((usize, Dir)),
    // Repetition ?
}
impl Debug for MyExtension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MyExtension::ShortestPath(sp) => write!(f, "SP{:?}", sp.target_node),
            MyExtension::NextNode(node) => write!(f, "NN{:?}", node),
        }
    }
}

pub struct MixedPath<'a, K: Kmer> {
    graph: &'a Graph<K>,
    start_node: (usize, Dir),
    extensions: Vec<MyExtension>,
}

impl<'a, K: Kmer> MixedPath<'a, K> {
    pub fn encode_seq(graph: &'a Graph<K>, seq: &DnaString) -> Self {
        let nodes = NodeIterator::new(graph, seq).unwrap().collect::<Vec<_>>();
        let start_node = nodes[0];
        println!("Start node: {:?}", start_node);
        let mut extensions = Vec::new();
        let mut position = 0;
        while let Some((path, path_length)) = ShortestPath::next(graph, &nodes, position) {
            if path_length >= 0 {
                position += path_length;
                let ext = MyExtension::ShortestPath(path);
                println!("{:?}: {:?}", path_length, ext);
                extensions.push(ext);
            } else {
                for _ in 0..path_length {
                    position += 1;
                    let ext = MyExtension::NextNode(nodes[position]);
                    println!("{:?}: {:?}",1, ext);
                    extensions.push(ext);
                }
            }
        }
        MixedPath {
            graph,
            start_node,
            extensions,
        }
    }

    pub fn decode_seq(&self) -> DnaString {
        let mut path = vec![self.start_node];
        for ext in self.extensions.iter() {
            match ext {
                MyExtension::ShortestPath(sp) => sp.extend_path(self.graph, &mut path),
                MyExtension::NextNode(node) => node.extend_path(self.graph, &mut path),
            }
        }
        self.graph.sequence_of_path(path.iter())
    }
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use crate::graph::Graph;

    use debruijn::kmer::Kmer3;
    use debruijn::{DnaSlice, dna_string::DnaString};

    const STRANDED : bool = true;
    const SEQ: DnaSlice = DnaSlice(&[2,2,2,1,1,1,1,2,2,2,0,0,0,0,0,1]);    // gggccccgggaaaaac
    const _SHORTEST: DnaSlice = DnaSlice(&[2,2,2,0,0,1]); // gggaac

    #[test]
    fn test_unitig_iterator() {
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);
        let mut unitig_iter = NodeIterator::new(&graph, &SEQ).unwrap();

        let mut path = Vec::new();
        while let Some(node) = unitig_iter.next().unwrap() {
            path.push(node);
        }
        let path_seq = graph.sequence_of_path(path.iter());
        let seq = DnaString::from_bytes(SEQ.0);
        assert!(path_seq == seq);

        // println!("Path: {:?}", path);

        // println!("{:?}", graph.base.sequences);
        // graph.print();
    }

    #[ignore]
    #[test]
    fn test_encode_seq() {
        let seq = DnaString::from_bytes(SEQ.0);
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);
        let path = MixedPath::encode_seq(&graph, &seq);
        for ext in path.extensions.iter() {
            println!("{:?}", ext);
        }
    }

    #[test]
    fn test_decode_seq() {
        let seq = DnaString::from_bytes(SEQ.0);
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);
        let path = MixedPath::encode_seq(&graph, &seq);
        let decoded_seq = path.decode_seq();

        assert_eq!(decoded_seq, seq);
    }
}