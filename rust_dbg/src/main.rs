use rust_dbg::{fasta_reader::FastaReader, graph, path};


use debruijn::{kmer, Kmer, Dir};
use debruijn::graph::DebruijnGraph;
use debruijn::dna_string::DnaString;

use std::fs::File;
use std::io::{BufWriter, BufReader, BufRead, Write};

use serde::{Serialize, Deserialize};
use bincode::serde::{encode_into_std_write, decode_from_std_read};

use std::time::Instant;

fn main() {
    // input/output files
    let path_graph = "../data/output/chr1/AalbF5_k31.fna";
    let path_bin = "../data/output/chr1/AalbF5_k31.bin";
    let path_fasta = "../data/input/chr1/AalbF5.fna";
    let path_chunks = "../data/output/chr1/chunks_nodes.AalbF5.fna";
    let path_reconstruct = "../data/output/chr1/reconstruct_nodes.AalbF5.fna";

    // params used for kmer construction by ggcat
    let stranded = false;
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    // let graph = make_graph::<Kmer31>(path_graph, stranded);
    // save_graph(&graph, path_bin);
    let graph = load_graph::<Kmer31>(path_bin);

    get_checkpoints(&graph, path_fasta, path_chunks);
    // reconstruct_fasta(&graph, path_chunks, path_reconstruct);
}

/// Create a graph from a unitigs file.
fn make_graph<K: Kmer + Send + Sync>(path_graph: &str, stranded: bool) -> DebruijnGraph<K,()> {
    print!("Creating graph (parallel)... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();

    let graph = graph::graph_from_unitigs(path_graph, stranded);

    let duration = start.elapsed();
    println!("done in {:?}", duration);

    graph
}

/// Save a graph to a binary file.
fn save_graph<K: Kmer + Serialize>(graph: &DebruijnGraph<K,()>, path_bin: &str) {
    print!("Saving graph... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();

    let mut f = BufWriter::new(File::create(path_bin).unwrap());
    let config = bincode::config::standard();
    let _ = encode_into_std_write(graph, &mut f, config).unwrap();

    let duration = start.elapsed();
    println!("done in {:?}", duration);
}

/// Load a graph from a binary file.
fn load_graph<K: Kmer + for<'a> Deserialize<'a>>(path_bin: &str) -> DebruijnGraph<K,()> {
    print!("Loading graph... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();

    let mut f = BufReader::new(File::open(path_bin).unwrap());
    let config = bincode::config::standard();
    let graph = decode_from_std_read(&mut f, config).unwrap();

    let duration = start.elapsed();
    println!("done in {:?}", duration);

    graph
}

/// Decompose the records in a fasta file into a suite a nodes in the graph and save them to a file.
fn get_checkpoints<K: Kmer>(graph: &DebruijnGraph<K,()>, path_fasta: &str, path_chunks: &str) {
    let fasta_reader = FastaReader::new(path_fasta).unwrap();
    let mut file_chunks = BufWriter::new(File::create(path_chunks).unwrap());

    let mut count: usize = 0;
    for record in fasta_reader.skip(count) {
        count += 1;
        // dividing the record into chunks
        println!("Processing record: {}", record.header());
        let start = Instant::now();
        let path = path::get_checkpoints_bfs(&graph, &record).unwrap();
        let duration = start.elapsed();
        println!("  - divided record into {} chunks\n  - time elapsed: {:?}", path.len(), duration);

        // saving them to the file
        writeln!(file_chunks, ">{}\tlen: {}\tdone in: {:?}", record.header(), path.len(), duration).unwrap();
        for chunk in path {
            let start_kmer = chunk.0;
            let end_kmer = chunk.1;
            writeln!(file_chunks, "{:?}\t\t{:?}", start_kmer, end_kmer).unwrap();
        }
        file_chunks.flush().unwrap();
    }
}

/// Reconstruct the records from a fasta file based on their decomposition into chunks.
fn reconstruct_fasta<K: Kmer>(graph: &DebruijnGraph<K,()>, file_chunks: &str, file_res: &str) {
    fn parse_node(kmer: &str) -> (usize,Dir) {
        let parts = kmer.trim_matches(|c| c == '(' || c == ')').split(',').map(str::trim).collect::<Vec<_>>();
        let id = parts[0].parse::<usize>().unwrap();
        let dir = match parts[1] {
            "Left" => Dir::Left,
            "Right" => Dir::Right,
            _ => panic!("Invalid direction"),
        };
        (id, dir)
    }

    let file_chunks = BufReader::new(File::open(file_chunks).unwrap());
    let mut file_res = BufWriter::new(File::create(file_res).unwrap());
    let mut header = String::new();
    let mut record = DnaString::new();
    let mut nb_chunks: usize = 0;
    for line in file_chunks.lines() {
        let line = line.unwrap();
        if line.starts_with('>') {
            // save the previous record and create a new one
            if nb_chunks > 0 {
                writeln!(file_res, ">{}\tnb chunks: {}\tlen: {}\t", header, nb_chunks, record.len()).unwrap();
                writeln!(file_res, "{}", record.to_string()).unwrap();
                file_res.flush().unwrap();
                record = DnaString::new();
                nb_chunks = 0;
            }
            header = line.trim_start_matches('>').split("\t").next().unwrap().to_owned();
            println!("Processing record: {}", header);
            continue;
        }
        // parse the line to get the start and end nodes
        let mut parts = line.split("\t\t");
        let start = parse_node(parts.next().unwrap());
        let end = parse_node(parts.next().unwrap());

        // get the path between the start and end nodes
        nb_chunks += 1;
        let chunk = path::get_shortest_path_bfs(graph, start, end).unwrap();
        if nb_chunks == 1 {
            record = chunk;
        } else {
            record.extend(chunk.to_bytes().into_iter().skip(K::k()-1));
        }
    }
    // save the last record
    if nb_chunks > 0 {
        writeln!(file_res, ">{}\tnb chunks: {}\tlen: {}\t", header, nb_chunks, record.len()).unwrap();
        writeln!(file_res, "{}", record.to_string()).unwrap();
        file_res.flush().unwrap();
    }
}