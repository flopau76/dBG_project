//! Defines how to encode a path in a debruijn Graph. Pri
 
use debruijn::{dna_string::DnaString, Dir, Kmer};

use crate::parse_node;
use crate::graph::Graph;
use self::{node_iterator::NodeIterator, shortest_path::ShortestPath};

use std::fmt::Debug;

pub mod node_iterator;
pub mod shortest_path;

/// Some information about how to extend a path in a given graph
trait Extension {
    /// Get the next extension and its length, given the full path and the current node in the path (note: path[position] is the last encoded node)
    fn get_next_ext<K: Kmer>(graph: &Graph<K>, path: &Vec<(usize, Dir)>, position: usize) -> Option<(Self, usize)> where Self: Sized;
    /// Extend the current path given the extension
    fn extend_path<K: Kmer>(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>);
}

#[derive(Copy, Clone)]
/// An enum containing all possible ways to encode a path extension.
enum MyExtension {
    ShortestPath(ShortestPath),
    NextNode((usize, Dir)),
}
impl Debug for MyExtension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MyExtension::ShortestPath(sp) => write!(f, "SP{:?}", sp.target_node),
            MyExtension::NextNode(nn) => write!(f, "NN{:?}", nn),
        }
    }
}
impl MyExtension {
    fn to_string(&self) -> String {
        match self {
            MyExtension::ShortestPath(sp) => format!("SP{:?}", sp.target_node),
            MyExtension::NextNode(nn) => format!("NN{:?}", nn),
        }
    }
    fn from_string(s: &str) -> Self {
        // let s = s.split_once(":").unwrap().1.trim();
        if s.starts_with("SP") {
            let target_node = parse_node(s[2..].trim()).unwrap();
            MyExtension::ShortestPath(ShortestPath {target_node})
        } else if s.starts_with("NN") {
            let target_node = parse_node(s[2..].trim()).unwrap();
            MyExtension::NextNode(target_node)
        } else {
            panic!("Unknown extension type")
        }
    }
    fn extend_path<K: Kmer>(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>) {
        match self {
            MyExtension::ShortestPath(sp) => sp.extend_path(graph, path),
            MyExtension::NextNode(nn) => path.push(*nn),
        }
    }
}

/// A struct that represents a path in a de Bruijn graph using a start node and a list of extensions.
pub struct MixedPath<'a, K: Kmer> {
    graph: &'a Graph<K>,
    start_node: (usize, Dir),
    extensions: Vec<MyExtension>,
}

impl<'a, K: Kmer> MixedPath<'a, K> {
    /// Encode the path into a MixedPath.
    pub fn encode_seq(graph: &'a Graph<K>, seq: &DnaString) -> Self {
        let nodes = NodeIterator::new(graph, seq).unwrap().collect::<Vec<_>>();
        let start_node = nodes[0];
        let mut extensions = Vec::new();
        let mut position = 0;
        while let Some((shortest_path, length)) = ShortestPath::get_next_ext(graph, &nodes, position) {
            if length >= shortest_path::MIN_LENGTH {
                extensions.push(MyExtension::ShortestPath(shortest_path));
            } else {
                extensions.extend(
                    nodes[position+1 .. position+1+length]
                        .iter()
                        .map(|&node| MyExtension::NextNode(node))
                );
            }
            position += length as usize;
        }
        MixedPath {
            graph,
            start_node,
            extensions,
        }
    }
    /// Decode the path into a DnaString.
    pub fn decode_seq(&self) -> DnaString {
        let mut path = vec![self.start_node];
        for ext in self.extensions.iter() {
            ext.extend_path(self.graph, &mut path);
        }
        self.graph.sequence_of_path(path.iter())
    }

    /// Transform the path into a string representation (to save to text file)
    pub fn to_string(&self) -> String {
        let mut result = String::new();
        result.push_str(&format!("Start node: {:?}", self.start_node));
        for ext in &self.extensions {
            result.push_str(&format!("\n{}", ext.to_string()));
        }
        result
    }

    /// Create a path from its string representation (to load from text file)
    pub fn from_string(s: &str, graph: &'a Graph<K>) -> Result<Self, Box<dyn std::error::Error>> {
        let mut lines = s.lines();
        let start_node_line = lines.next().ok_or("Missing start node line")?;
        let start_node = parse_node(&start_node_line[12..])?;

        let mut extensions = Vec::new();
        for line in lines {
            let ext = MyExtension::from_string(line);
            extensions.push(ext);
        }

        Ok(MixedPath {
            graph,
            start_node,
            extensions,
        })
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
        assert_eq!(seq, decoded_seq);
    }

    #[test]
    fn test_to_from_string() {
        let seq = DnaString::from_bytes(SEQ.0);
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);
        let path = MixedPath::encode_seq(&graph, &seq);

        let path_str = path.to_string();
        let loaded_path = MixedPath::from_string(&path_str, &graph).unwrap();

        assert_eq!(seq, loaded_path.decode_seq());
    }
}