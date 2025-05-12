//! Defines how to encode a path in a debruijn Graph. Pri
 
use debruijn::{dna_string::DnaString, Dir, Kmer};

use crate::graph::Graph;
use self::{node_iterator::NodeIterator, shortest_path::ShortestPath};

use std::{fmt::Debug, path::Path};
use std::io::{BufReader, BufRead};

pub mod node_iterator;
pub mod shortest_path;

/// Some information about how to extend a path in a given graph
trait Extension : Debug {
    /// Get the next extension and the number of incoded nodes, given the full path and the current node in the path
    fn next<K: Kmer>(graph: &Graph<K>, path: &Vec<(usize, Dir)>, position: usize) -> Option<Self> where Self: Sized;
    /// Extend the current path given the extension
    fn extend_path<K: Kmer>(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>);
    /// Get the length of this extension
    fn length(&self) -> usize;
}

/// The easiest way to extend the current path is to give the following node
pub struct NextNode {
    next_node: (usize, Dir),
}

impl Debug for NextNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "NN{:?}", self.next_node)
    }
}

impl Extension for NextNode {
    fn next<K: Kmer>(_graph: &Graph<K>, path: &Vec<(usize, Dir)>, position: usize) -> Option<Self> {
        if position + 1 < path.len() {
            Some(NextNode { next_node: path[position + 1] })
        } else {
            None
        }
    }
    fn extend_path<K: Kmer>(&self, _graph: &Graph<K>, path: &mut Vec<(usize, Dir)>) {
        path.push(self.next_node);
    }
    fn length(&self) -> usize {
        1
    }
}

enum MyExtension {
    ShortestPath(ShortestPath),
    NextNode(NextNode),
}
impl Debug for MyExtension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MyExtension::ShortestPath(sp) => write!(f, "{:?}", sp),
            MyExtension::NextNode(nn) => write!(f, "{:?}", nn),
        }
    }
}
impl Extension for MyExtension {
    fn next<K: Kmer>(graph: &Graph<K>, path: &Vec<(usize, Dir)>, position: usize) -> Option<Self> {
        if let Some(sp) = ShortestPath::next(graph, path, position) {
            Some(MyExtension::ShortestPath(sp))
        } else if let Some(nn) = NextNode::next(graph, path, position) {
            Some(MyExtension::NextNode(nn))
        } else {
            None
        }
    }
    fn extend_path<K: Kmer>(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>) {
        match self {
            MyExtension::ShortestPath(sp) => sp.extend_path(graph, path),
            MyExtension::NextNode(nn) => nn.extend_path(graph, path),
        }
    }
    fn length(&self) -> usize {
        match self {
            MyExtension::ShortestPath(sp) => sp.length(),
            MyExtension::NextNode(nn) => nn.length(),
        }
    }
}

pub struct MixedPath<'a, K: Kmer> {
    name: Option<String>,
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
        while let Some(path) = ShortestPath::next(graph, &nodes, position) {
            let length = path.length();
            if length >= 0 {
                position += length;
                println!("{}: {:?}", length, path);
                extensions.push(MyExtension::ShortestPath(path));
            } else {
                for _ in 0..length {
                    position += 1;
                    let next_node = NextNode{next_node: nodes[position]};
                    println!("{}: {:?}", 1, next_node);
                    extensions.push(MyExtension::NextNode(next_node));
                }
            }
        }
        MixedPath {
            name: None,
            graph,
            start_node,
            extensions,
        }
    }

    pub fn decode_seq(&self) -> DnaString {
        let mut path = vec![self.start_node];
        for ext in self.extensions.iter() {
            ext.extend_path(self.graph, &mut path);
        }
        self.graph.sequence_of_path(path.iter())
    }

    pub fn from_text_file(graph: &'a Graph<K>, file: &Path) -> Vec<Self> {
        let file = BufReader::new(std::fs::File::open(file).unwrap());
        let mut lines_iter = file.lines().map(|l| l.unwrap());
        while let Some(line) = lines_iter.next() {
            if line.starts_with('>') {
                let name = line[1..].to_string();
                let first_line = lines_iter.next().unwrap();
                let start_node = first_line.split(":").nth(1).unwrap();
            }
        }
        todo!();
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