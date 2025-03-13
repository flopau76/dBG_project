mod dbg;

use bio::io::fasta;

use debruijn::{kmer, Kmer, Exts, Vmer, dna_string::DnaString};
use boomphf::hashmap::BoomHashMap2;

use std::fs::File;
use std::io::{BufWriter, BufReader, Write};
use std::path;
use bincode::serde::{encode_into_std_write, decode_from_std_read};
use serde::{Serialize, Deserialize};

use std::time::Instant;

fn main() {
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;
    let path_input = "data/unitigs_chr1_k31.fna";
    let path_output = "data/graph_old_chr1_k31.serde";

    // let seq_iter = dbg::get_seq_iterator(path_input);
    // let _test: std::collections::HashSet<Kmer31> = dbg::get_kmers(seq_iter, true);

    // let graph = make_graph::<Kmer31>(path_input);    // ~1100s
    // save_graph(&graph, path_output); // ~40s

    // let graph = load_graph::<Kmer31>(path_output);  // ~80s

    // let haplo_path = "data/AalbF5_chr1.fna";
    // let path_output = "data/walk_output.txt";
    // let seq = load_seq(haplo_path);
    // walk_graph(&graph, seq, path_output);

}

// fn make_graph<K: Kmer>(path: &str) -> BoomHashMap2<K, Exts, ()> {
//     println!("Starting HashMap construction");
//     let start = Instant::now();
//     let graph = dbg::make_graph::<K>(path);
//     let duration = start.elapsed();
//     println!("Time to make map: {:?}", duration);

//     graph
// }

fn save_graph<K: Kmer + Serialize>(graph: &BoomHashMap2<K, Exts, ()>, path_output: &str) {
    println!("Writing map to file");
    let start = Instant::now();
    let mut f = BufWriter::new(File::create(path_output).unwrap());
    let config = bincode::config::standard();
    let nb_bytes = encode_into_std_write(graph, &mut f, config).unwrap();
    let duration = start.elapsed();
    println!("Time to write map: {:?}", duration);

    println!("Wrote {} bytes to file {}", nb_bytes, path_output);
}

fn load_graph<K: Kmer + for<'a> Deserialize<'a>>(path: &str) -> BoomHashMap2<K, Exts, ()> {
    println!("Reading map from file");
    let start = Instant::now();
    let mut f = BufReader::new(File::open(path).unwrap());
    let config = bincode::config::standard();
    let graph = decode_from_std_read(&mut f, config).unwrap();
    let duration = start.elapsed();
    println!("Time to read map: {:?}", duration);

    graph
}

fn walk_graph<K: Kmer>(graph: &BoomHashMap2<K, Exts, ()>, seq: impl Vmer, path_output: &str) {
    println!("Walking graph");
    let start = Instant::now();
    let mut file = BufWriter::with_capacity(65536, File::create(path_output).unwrap());
    for kmer in seq.iter_kmers() {
        if let Some(exts) = graph.get(&kmer) {
            writeln!(file, "{}\t{}", exts.0.num_exts_l(), exts.0.num_exts_r()).unwrap();
        }
    }
    file.flush().unwrap();
    let duration = start.elapsed();
    println!("Walk output written to {} in: {:?}", path_output, duration);
}

fn load_seq(path: &str) -> DnaString {
    let reader = fasta::Reader::from_file(path).unwrap();
    let mut records = reader.records();
    DnaString::from_acgt_bytes(records.next().unwrap().unwrap().seq())
}