pub mod fasta_reader;
pub mod dbg;
pub mod path;

use bincode::enc::write;
use fasta_reader::FastaReader;
use dbg::Graph;
use debruijn::{kmer, Kmer, bits_to_ascii, Mer};
use path::get_shortest_path;

use std::fs::File;
use std::io::{BufWriter, BufReader, Write};

use serde::{Serialize, Deserialize};
use bincode::serde::{encode_into_std_write, decode_from_std_read};

use std::time::{Duration, Instant};

#[allow(unused)]
fn main() {
    // input files
    let path_graph = "../data/output/chr1/AalbF5_k31.fna";
    let path_bin = "../data/output/chr1/AalbF5_k31.bin";
    let path_haplo = "../data/input/chr1/AalbF5_splitN.fna";

    let output_path = "../data/output/chr1/AalbF5_splitN.path.fna";

    // params used for kmer construction by ggcat
    let canon = true;
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    // let graph = make_graph::<Kmer31>(path_graph, canon);
    // save_graph(&graph, path_bin);
    let graph = load_graph::<Kmer31>(path_bin);
    println!("{}", graph);

    // let mut haplo = fasta_reader::FastaReader::new(path_haplo).unwrap()
    // make_report(&graph, haplo, output_walk_5);

    let mut haplo = fasta_reader::FastaReader::new(path_haplo).unwrap();
    get_path(&graph, &mut haplo, output_path);
}

fn make_graph<K: Kmer + Send + Sync>(path: &str, canon: bool) -> Graph<K> {
    print!("Creating graph (parallel)... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();
    let graph = dbg::Graph::<K>::from_unitigs(path, canon);
    let duration = start.elapsed();
    println!("done in {:?}", duration);

    graph
}

fn save_graph<K: Kmer + Serialize>(graph: &Graph<K>, path_output: &str) {
    print!("Saving graph... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();
    let mut f = BufWriter::new(File::create(path_output).unwrap());
    let config = bincode::config::standard();
    let _nb_bytes = encode_into_std_write(graph, &mut f, config).unwrap();
    let duration = start.elapsed();
    println!("done in {:?}", duration);
}

fn load_graph<K: Kmer + for<'a> Deserialize<'a>>(path: &str) -> Graph<K> {
    print!("Loading graph... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();
    let mut f = BufReader::new(File::open(path).unwrap());
    let config = bincode::config::standard();
    let graph = decode_from_std_read(&mut f, config).unwrap();
    let duration = start.elapsed();
    println!("done in {:?}", duration);

    graph
}

// Input graph: only AalbF3
//   - 255874974 nodes
//   - 259767078 edges
// Haplo contains:
//   - kmers: 373823221
//   - breakpoints: 52598113  (>1)      217  (<1)

// Input graph: only AalbF5
//   - 232971101 nodes
//   - 236446568 edges
// Haplo contains:
//   - kmers: 337698822
//   - breakpoints: 46619218  (>1)      862  (<1)

// Input graph: both haplo
//    - 376592961 nodes
//    - 383504450 edges
// Haplo AalbF3 contains:
//    - 373823221 kmers
//    - breakpoints: 60820242  (>1)     188  (<1)
// Haplo AalbF5 contains:
//    - 337698822 kmers
//    - breakpoints: 55574344  (>1)     470  (<1)
fn make_report<K: Kmer>(graph: &Graph<K>, mut haplo: FastaReader, output_walk: &str) {
    // overall graph stats
    print!("Iterating graph... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();
    let nb_nodes = graph.len();
    let mut nb_edges: usize = 0;
    let mut degree_histo = vec![0; 9];
    for id in 0..nb_nodes {
        let exts = *graph.get_exts(id);
        let degree = (exts.num_exts_l() + exts.num_exts_r()) as usize;
        nb_edges += degree;
        degree_histo[degree] += 1;
    }
    nb_edges /= 2;
    let duration = start.elapsed();
    println!("done in {:?}", duration);
    println!("Graph contains:\n  - {} nodes\n  - {} edges\n  degree histogramm: {:?}", nb_nodes, nb_edges, degree_histo);

    // haplo stats
    print!("Iterating haplo... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();
    let mut nb_kmers: usize = 0;
    let mut junctions: usize = 0;
    let mut dead_ends: usize = 0;

    // let file = BufWriter::new(File::create(output_walk).unwrap());
    while let Some(mut record) = haplo.next() {
        for (kmer, flip) in record.iter_kmers(graph.canon) {
            nb_kmers += 1;
            let id = graph.get_key_id(&kmer).expect("Kmer not found in graph");
            let exts = *graph.get_exts(id);
            let (mut l, mut r) = (exts.num_exts_l(), exts.num_exts_r());
            if flip { (l, r) = (r, l); }
            match r {
                0 => { dead_ends += 1; },
                1 => { },
                _ => { junctions += 1; },
            }
            // writeln!(file, "{},{}", l, r).unwrap();
        }
    }
    let duration = start.elapsed();
    println!("done in {:?}", duration);
    println!("Haplo contains:\n  - kmers: {}\n  - breakpoints: {}  (>1)\t\t{}  (<1)", nb_kmers, junctions, dead_ends);
}

fn get_path<K: Kmer>(graph: &Graph<K>, fasta_reader: &mut FastaReader, save_file: &str) {
    println!("Looking for path in graph... ");
    let mut file = BufWriter::new(File::create(save_file).unwrap());
    let mut count: usize = 1;
    while let Some(mut record) = fasta_reader.next() {
        println!("Processing record: {}", record.header());
        let start = Instant::now();
        let mut kmer_iter = record.iter_kmers::<K>(graph.canon);
        let (first_kmer, _) = kmer_iter.next().unwrap();
        let (last_kmer, _) = kmer_iter.last().unwrap();
        let path = get_shortest_path(graph, first_kmer, last_kmer);
        if path.is_err() {
            println!("  - no path found");
            continue;
        }
        let path = path.unwrap();
        let path_str = path.0.iter().map(|&b| bits_to_ascii(b) as char).collect::<String>();
        writeln!(file, ">path_{}", count).unwrap();
        writeln!(file, "{}", path_str).unwrap();
        count += 1;
        let duration = start.elapsed();
        println!("  - path length: {}\n  - time elapsed: {:?}", path.len(), duration);
        file.flush().unwrap();
    }
}