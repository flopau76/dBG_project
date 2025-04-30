use rust_dbg::stats;
use rust_dbg::{fasta_reader::FastaReader, graph::Graph, path};

use debruijn::{kmer, Dir, Kmer};
use debruijn::dna_string::DnaString;

use regex::Regex;

use std::fs::File;
use std::io::{BufWriter, BufReader, BufRead, Write};

use serde::{Serialize, Deserialize};
use bincode::serde::{encode_into_std_write, decode_from_std_read};

use std::time::Instant;

fn main() {
    // input/output files
    let path_graph = "../data/output/chr1/graph_k31.fna";
    let path_bin = "../data/output/chr1/graph_k31.bin";
    let path_fasta = "../data/input/chr1/AalbF5_splitN.fna";
    let path_checkpoints = "../data/output/chr1/checkpoints_AalbF5_in_2.fna";

    // params used for kmer construction by ggcat
    let stranded = false;
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    // let graph = make_graph::<Kmer31>(path_graph, stranded);
    // save_graph(&graph, path_bin);
    let graph = load_graph::<Kmer31>(path_bin);

    stats::stats_depth_bfs(&graph, 30, 100);

    // make_checkpoints(&graph, path_fasta);
    // let checkpoints = load_checkpoints(path_checkpoints);
    // stats_checkpoints::<Kmer31>(&graph, checkpoints);
    // reconstruct_fasta(&graph, checkpoints);

    // let haplo = FastaReader::new(path_fasta).unwrap();
    // rust_dbg::stats::stats_haplo(&graph, haplo);
    // rust_dbg::stats::stats_graph(&graph);
}

/// Create a graph from a unitigs file.
pub fn make_graph<K: Kmer + Send + Sync>(path_graph: &str, stranded: bool) -> Graph<K> {
    eprint!("Creating graph (parallel)... ");
    std::io::stderr().flush().unwrap();
    let start = Instant::now();

    let graph = Graph::from_unitigs(path_graph, stranded);

    let duration = start.elapsed();
    eprintln!("done in {:?}", duration);

    graph
}

/// Save a graph to a binary file.
pub fn save_graph<K: Kmer + Serialize>(graph: &Graph<K>, path_bin: &str) {
    eprint!("Saving graph... ");
    std::io::stderr().flush().unwrap();
    let start = Instant::now();

    let mut f = BufWriter::new(File::create(path_bin).unwrap());
    let config = bincode::config::standard();
    let _ = encode_into_std_write(graph, &mut f, config).unwrap();

    let duration = start.elapsed();
    eprintln!("done in {:?}", duration);
}

/// Load a graph from a binary file.
pub fn load_graph<K: Kmer + for<'a> Deserialize<'a>>(path_bin: &str) -> Graph<K> {
    eprint!("Loading graph... ");
    std::io::stderr().flush().unwrap();
    let start = Instant::now();

    let mut f = BufReader::new(File::open(path_bin).unwrap());
    let config = bincode::config::standard();
    let graph = decode_from_std_read(&mut f, config).unwrap();

    let duration = start.elapsed();
    eprintln!("done in {:?}", duration);

    graph
}

/// Decomposes the fasta file into a suite of checkpoints.
/// Result is writen directly to stdout.
pub fn make_checkpoints<K: Kmer>(graph: &Graph<K>, path_fasta: &str) {
    let fasta_reader = FastaReader::new(path_fasta).unwrap();

    for record in fasta_reader {
        // dividing the record into chunks
        eprintln!("Processing record: {}", record.header());
        println!(">{}", record.header());
        let start = Instant::now();
        let path = path::get_checkpoints_bfs(&graph, &record).unwrap();
        let duration = start.elapsed();
        eprintln!("  - divided record into {} chunks\n  - time elapsed: {:?}", path.len(), duration);
    }
}

/// Reads the checkpoints computed by `get_checkpoints` and returns them as a vector of tuples.
/// Each tuple contains the header and a vector of checkpoints.
pub fn load_checkpoints(path_checkpoints: &str) -> Vec<(String, Vec<((usize, Dir), (usize, Dir))>)> {
    eprint!("Loading checkpoints... ");
    std::io::stderr().flush().unwrap();
    let start = Instant::now();

    let file = BufReader::new(File::open(path_checkpoints).unwrap());
    let mut current_header = String::new();
    let mut current_checkpoints = Vec::new();
    let mut result = Vec::new();

    let pattern = r"(\d+) unitigs, starting at (\d+):\s+\((\d+), (Left|Right)\)\s+\((\d+), (Left|Right)\)";
    let re = Regex::new(pattern).unwrap();

    for line in file.lines() {
        let line = line.unwrap();
        if line.starts_with('>') {
            if !current_checkpoints.is_empty() {
                result.push((current_header, current_checkpoints));
                current_checkpoints = Vec::new();
            }
            current_header = line[1..].to_string();
            continue;
        }
        // parse the line using regex
        let caps = re.captures(&line).unwrap_or_else(|| {
            panic!("Failed to parse line: {}", line);
        });
        let start_id = caps[3].parse::<usize>().unwrap();
        let start_dir = if &caps[4] == "Left" { Dir::Left } else { Dir::Right };
        let end_id = caps[5].parse::<usize>().unwrap();
        let end_dir = if &caps[6] == "Left" { Dir::Left } else { Dir::Right };
        current_checkpoints.push((
            (start_id, start_dir),
            (end_id, end_dir),
        ));
    }

    if !current_checkpoints.is_empty() {
        result.push((current_header, current_checkpoints));
    }

    let duration = start.elapsed();
    eprintln!("done in {:?}", duration);

    result
}

/// Reconstruct the records given their list of checkpoints.
/// Result is writen directly to stdout.
pub fn reconstruct_fasta<K: Kmer>(graph: &Graph<K>, checkpoints: Vec<(String, Vec<((usize, Dir), (usize, Dir))>)>) {
    for (header, current_checkpoints) in checkpoints {
        eprintln!("Processing {}", header);
        println!(">{}", header);
        let mut seq = DnaString::new();
        for (start_node, end_node) in current_checkpoints {
            let path = path::get_shortest_path_double_bfs(graph, start_node, end_node).unwrap();
            seq = graph.sequence_of_path(path.iter());
            println!("{}", seq.prefix(seq.len()-K::k()+1));
        }
        println!("{}", seq.suffix(K::k()-1));
    }
}

/// Counts the number of breakpoints per shortest path.
/// Result is writen directly to stdout.
pub fn stats_checkpoints<K: Kmer>(graph: &Graph<K>, checkpoints: Vec<(String, Vec<((usize, Dir), (usize, Dir))>)>) {

    for (header, current_checkpoints) in checkpoints {
        eprintln!("Processing {}", header);
        println!(">{}", header);
        for (start_node, end_node) in current_checkpoints {

            // get the path between the start and end nodes
            let path = path::get_shortest_path_double_bfs(graph, start_node, end_node).unwrap();

            // walk the path to get the number of breakpoints
            let mut breakpoints_l: usize = 0;
            let mut breakpoints_r: usize = 0;
            for node in path.iter() {
                let exts = graph.get_node(node.0).exts();
                let (l, r) = match node.1 {
                    Dir::Left => (exts.num_exts_l(), exts.num_exts_r()),
                    Dir::Right => (exts.num_exts_r(), exts.num_exts_l()),
                };
                if l > 1 {
                    breakpoints_l += 1;
                }
                if r > 1 {
                    breakpoints_r += 1;
                }
            }
            println!("{}\t{}\t{}", path.len(), breakpoints_l, breakpoints_r);
        }
    }
}