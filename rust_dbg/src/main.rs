mod dbg;
mod kmers;

use debruijn::{kmer, Kmer};
use dbg::Graph;

use std::fs::File;
use std::io::{BufWriter, BufReader, Write};

use serde::{Serialize, Deserialize};
use bincode::serde::{encode_into_std_write, decode_from_std_read};

use std::time::Instant;

fn main() {
    let path_graph = "data/unitigs_chr1_k31.fna";
    let path_haplo = "data/AalbF5_chr1.fna";
    let path_bin = "data/graph_k31.bin";

    let canon = true;
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    let graph = make_graph::<Kmer31>(path_graph, canon);
    save_graph(&graph, path_bin);

    let graph = load_graph::<Kmer31>(path_bin);

    let kmers_iter = kmers::KmerIterator::<Kmer31>::new(path_haplo, canon).unwrap().map(|b| b.0);
    walk_graph(&graph, kmers_iter);
}

fn make_graph<K: Kmer>(path: &str, canon: bool) -> Graph<K> {
    println!("Creating graph");
    let start = Instant::now();
    let graph = dbg::Graph::<K>::from_fasta(path, canon);
    let duration = start.elapsed();
    println!("Time to create graph: {:?}", duration);

    graph
}

fn save_graph<K: Kmer + Serialize>(graph: &Graph<K>, path_output: &str) {
    println!("Writing graph to file");
    let start = Instant::now();
    let mut f = BufWriter::new(File::create(path_output).unwrap());
    let config = bincode::config::standard();
    let nb_bytes = encode_into_std_write(graph, &mut f, config).unwrap();
    let duration = start.elapsed();
    println!("Time to save graph: {:?}", duration);

    println!("Wrote {} bytes to file {}", nb_bytes, path_output);
}

fn load_graph<K: Kmer + for<'a> Deserialize<'a>>(path: &str) -> Graph<K> {
    println!("Reading graph from file");
    let start = Instant::now();
    let mut f = BufReader::new(File::open(path).unwrap());
    let config = bincode::config::standard();
    let graph = decode_from_std_read(&mut f, config).unwrap();
    let duration = start.elapsed();
    println!("Time to read graph: {:?}", duration);

    graph
}

// Time to read graph: 78.713671409s
// Walking graph
// Time to walk graph: 216.341636567s
// Edges left: [444, 281_670_030, 26_059_907, 14_435_009, 15_533_432]
// Edges right: [497, 282_586_013, 25_754_030, 14_057_128, 15_301_154]
fn walk_graph<K: Kmer>(graph: &Graph<K>, kmers: impl Iterator<Item=K>) {
    println!("Walking graph");
    let start = Instant::now();
    let mut edges_histo: [usize; 8] = [0;8];

    let output_path = "data/walk_output.txt";
    let mut file = BufWriter::new(File::create(output_path).unwrap());

    for kmer in kmers {
        let id = graph.get_key_id(&kmer).expect("Kmer not found in graph");
        let exts = *graph.get_exts(id);
        let degree = exts.num_exts_l() + exts.num_exts_r();
        edges_histo[degree as usize] += 1;
        write!(file, "{} ", degree).unwrap();
    }
    let duration = start.elapsed();
    println!("Time to walk graph: {:?}", duration);
    println!("Edges distribution: {:?}", edges_histo);
    println!("Saved walk to {}", output_path);
}