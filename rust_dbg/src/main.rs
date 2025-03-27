pub mod dbg;
pub mod fasta_reader;
pub mod path;

use debruijn::{kmer, Kmer};
use dbg::Graph;
use fasta_reader::FastaReader;

use std::fs::File;
use std::io::{BufWriter, BufReader, Write};

use serde::{Serialize, Deserialize};
use bincode::serde::{encode_into_std_write, decode_from_std_read};

use std::time::Instant;

#[allow(unused)]
fn main() {
    // input and output paths
    let path_graph = "../data/output/chr1/AalbF5_k31.fna";
    let path_bin = "../data/output/chr1/graph_k31.bin";

    let path_haplo_3 = "../data/input/chr1/AalbF3.fna";
    let path_haplo_5 = "../data/input/chr1/AalbF5.fna";
    let output_walk_3 = "../data/output/chr1_walk_AalbF3.txt";
    let output_walk_5 = "../data/output/chr1_walk_AalbF5.txt";

    // params used for kmer construction by ggcat
    let canon = true;
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    let graph = make_graph_parallel::<Kmer31>(path_graph, canon);
    // save_graph(&graph, path_bin);

    // let graph = load_graph::<Kmer31>(path_bin);

    // let mut haplo_3 = fasta_reader::FastaReader::new(path_haplo_3).unwrap();
    let mut haplo_5 = fasta_reader::FastaReader::new(path_haplo_5).unwrap();
    make_report(&graph, haplo_5, output_walk_5);
}

fn make_graph_parallel<K: Kmer + Send + Sync>(path: &str, canon: bool) -> Graph<K> {
    print!("Creating graph (parallel)... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();
    let graph = dbg::Graph::<K>::from_fasta_parallel(path, canon);
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
    for id in 0..nb_nodes {
        let exts = *graph.get_exts(id);
        nb_edges += (exts.num_exts_l() + exts.num_exts_r()) as usize;
    }
    nb_edges /= 2;
    let duration = start.elapsed();
    println!("done in {:?}", duration);
    println!("Graph contains:\n  - {} nodes\n  - {} edges\n", nb_nodes, nb_edges);

    // haplo stats
    print!("Iterating haplo... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();
    let mut nb_kmers: usize = 0;
    let mut junctions: usize = 0;
    let mut dead_ends: usize = 0;

    let file = BufWriter::new(File::create(output_walk).unwrap());

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