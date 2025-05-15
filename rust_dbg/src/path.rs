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

#[derive(Debug, Copy, Clone)]
/// An enum containing all possible ways to encode a path extension.
enum MyExtension {
    ShortestPath(ShortestPath),
    NextNode((usize, Dir)),
    Repetition((u8, u16)),    // (size of pattern, number of nodes)
}

impl MyExtension {
    fn to_string(&self) -> String {
        match self {
            MyExtension::ShortestPath(sp) => format!("SP{:?}", sp.target_node),
            MyExtension::NextNode(nn) => format!("NN{:?}", nn),
            MyExtension::Repetition((size, count)) => format!("R{}x{}", size, count),
        }
    }
    fn from_string(s: &str) -> Self {
        if s.starts_with("SP") {
            let target_node = parse_node(s[2..].trim()).unwrap();
            MyExtension::ShortestPath(ShortestPath {target_node})
        } else if s.starts_with("NN") {
            let target_node = parse_node(s[2..].trim()).unwrap();
            MyExtension::NextNode(target_node)
        } else if s.starts_with("R") {
            let (size, count) = s[1..].split_once("x").unwrap();
            let size = size.parse::<u8>().unwrap();
            let count = count.parse::<u16>().unwrap();
            MyExtension::Repetition((size, count))
        } else {
            panic!("Unknown extension type: {}", s);
        }
    }
    fn extend_path<K: Kmer>(&self, graph: &Graph<K>, path: &mut Vec<(usize, Dir)>) {
        match self {
            MyExtension::ShortestPath(sp) => sp.extend_path(graph, path),
            MyExtension::NextNode(nn) => path.push(*nn),
            MyExtension::Repetition((size, count)) => {
                for _ in 0..*count {
                    let prev = path[path.len()-*size as usize];
                    path.push(prev);
                }
            }
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
    /// Encode the list of nodes from a sequence into a MixedPath.
    pub fn encode_seq(graph: &'a Graph<K>, seq: &DnaString) -> Self {
        let nodes = NodeIterator::new(graph, seq).unwrap().collect::<Vec<_>>();
        let start_node = nodes[0];
        let mut extensions = Vec::new();
        let mut position = 0;
        // TODO: handle repetitions
        while let Some((shortest_path, length)) = ShortestPath::get_next_ext(graph, &nodes, position) {
            // let mut test_path = vec!(nodes[position]);
            // shortest_path.extend_path(graph, &mut test_path);
            // if test_path != nodes[position ..= position+length as usize] {
            //     eprintln!("Beware: paths not matching");
            //     assert_eq!(test_path.len(), length as usize + 1, "Both paths do not even have the same length");
            // }

            // if the shortest path is long enough, we encode its target node (u32 for the whole path)
            if length >= shortest_path::MIN_LENGTH {
                extensions.push(MyExtension::ShortestPath(shortest_path));
            }
            // otherwise, we encode all its nodes directly (2 bits per node)
            else {
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
    fn test_node_iterator_offset() {
        let seq = DnaString::from_dna_string("aaccggtt");
        let graph = Graph::<Kmer3>::from_seq_serial(&seq, true);

        let seq_start = DnaString::from_dna_string("aacc");
        let unitig_iter_start = NodeIterator::new(&graph, &seq_start).unwrap();
        assert_eq!(unitig_iter_start.start_offset, 0);
        assert_eq!(unitig_iter_start.end_offset, Some(4));

        let seq_end = DnaString::from_dna_string("ggtt");
        let unitig_iter_end = NodeIterator::new(&graph, &seq_end).unwrap();
        assert_eq!(unitig_iter_end.start_offset, 4);
        assert_eq!(unitig_iter_end.end_offset, None);

        // let seq_middle = DnaString::from_dna_string("ccgg");
        // let unitig_iter_middle = NodeIterator::new(&graph, &seq_middle).unwrap();
    }

    #[test]
    fn test_node_iterator() {
        // seq -> node_list -> seq
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);
        let mut unitig_iter = NodeIterator::new(&graph, &SEQ).unwrap();

        let mut path = Vec::new();
        while let Some(node) = unitig_iter.next().unwrap() {
            path.push(node);
        }
        let path_seq = graph.sequence_of_path(path.iter());
        let seq = DnaString::from_bytes(SEQ.0);
        assert!(path_seq == seq);
    }

    #[ignore]
    #[test]
    fn test_encode_seq() {
        // seq -> (node_list) -> path
        let seq = DnaString::from_bytes(SEQ.0);
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);
        let path = MixedPath::encode_seq(&graph, &seq);
        for ext in path.extensions.iter() {
            println!("{:?}", ext);
        }
    }

    #[test]
    fn test_code_seq() {
        // seq -> (node_list) -> path -> (node_list) -> seq
        let seq = DnaString::from_bytes(SEQ.0);
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);

        let decoded_seq = MixedPath::encode_seq(&graph, &seq).decode_seq();
        assert_eq!(seq, decoded_seq);
    }

    #[test]
    fn test_path_to_string() {
        // seq -> (node_list) -> path -> string -> path -> (node_list) -> seq
        let seq = DnaString::from_bytes(SEQ.0);
        let graph = Graph::<Kmer3>::from_seq_serial(&SEQ, STRANDED);

        let path_str = MixedPath::encode_seq(&graph, &seq).to_string();
        let decoded_seq = MixedPath::from_string(&path_str, &graph).unwrap().decode_seq();
        assert_eq!(seq, decoded_seq);
    }
}



#[cfg(test)]
mod real_test {
    use std::path::PathBuf;

    use super::*;
    use crate::graph::Graph;
    use crate::fasta_reader::FastaReader;
    use debruijn::kmer;

    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    #[test]
    fn test_node_iterator_real() {
        // seq -> node_list -> seq
        let path_graph = PathBuf::from("/home/florence/Documents/dbg_project/data/output/chr1/AalbF5_k31.bin");
        let graph = Graph::<Kmer31>::load_from_binary(&path_graph).unwrap();

        let path_fasta = PathBuf::from("/home/florence/Documents/dbg_project/data/input/chr1/AalbF5_splitN.fna");
        let fasta_reader = FastaReader::new(path_fasta).unwrap();

        for record in fasta_reader {
            println!("Testing record {}", record.header());
            let seq = record.dna_string();
            let mut unitig_iter = NodeIterator::new(&graph, &seq).unwrap();

            let mut path = Vec::new();
            while let Some(node) = unitig_iter.next().unwrap() {
                path.push(node);
            }
            let path_seq = graph.sequence_of_path(path.iter());
            let path_seq = path_seq.slice(unitig_iter.start_offset, path_seq.len() - unitig_iter.end_offset.unwrap_or(0)).to_owned();
            assert_eq!(path_seq, seq);
        }
    }
}