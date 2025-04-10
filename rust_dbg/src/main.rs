use rust_dbg::{fasta_reader::FastaReader, graph, path};


use debruijn::{kmer, Kmer, Vmer};
use debruijn::graph::DebruijnGraph;
use debruijn::dna_string::DnaString;

use std::fs::File;
use std::io::{BufWriter, BufReader, Write};

use serde::{Serialize, Deserialize};
use bincode::serde::{encode_into_std_write, decode_from_std_read};

use std::time::Instant;

fn main() {
    // input/output files
    let path_graph = "../data/output/chr1/AalbF5_k31.fna";
    let path_bin = "../data/output/chr1/AalbF5_k31.bin";
    let path_haplo = "../data/input/chr1/AalbF5_splitN.fna";
    let path_bfs = "../data/output/chr1/path.AalbF5_splitN.fna";
    let path_chunks = "../data/output/chr1/checkpoints.AalbF5_splitN.fna";

    // params used for kmer construction by ggcat
    let stranded = false;
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    let fasta_reader = FastaReader::new(path_haplo).unwrap();

    // let graph = make_graph::<Kmer31>(path_graph, stranded);
    // save_graph(&graph, path_bin);
    let graph = load_graph::<Kmer31>(path_bin);

    get_shortest_path(&graph, fasta_reader, path_bfs);
    // get_checkpoints(&graph, fasta_reader, path_chunks);
}

/// Create a graph from a unitigs file
fn make_graph<K: Kmer + Send + Sync>(path: &str, stranded: bool) -> DebruijnGraph<K,()> {
    print!("Creating graph (parallel)... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();
    let graph = graph::graph_from_unitigs(path, stranded);
    let duration = start.elapsed();
    println!("done in {:?}", duration);

    graph
}

/// Save a graph to a binary file
fn save_graph<K: Kmer + Serialize>(graph: &DebruijnGraph<K,()>, path_bin: &str) {
    print!("Saving graph... ");
    std::io::stdout().flush().unwrap();
    let start = Instant::now();
    let mut f = BufWriter::new(File::create(path_bin).unwrap());
    let config = bincode::config::standard();
    let _nb_bytes = encode_into_std_write(graph, &mut f, config).unwrap();
    let duration = start.elapsed();
    println!("done in {:?}", duration);
}

/// Load a graph from a binary file
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

// /// Print some stats about the graph
// fn stats_graph<K: Kmer>(graph: &DebruijnGraph<K,()>) {
//     print!("Iterating graph... ");
//     std::io::stdout().flush().unwrap();
//     let start = Instant::now();
//     let nb_nodes = graph.len();
//     let mut nb_edges: usize = 0;
//     let mut degree_histo = vec![0; 9];
//     for (_kmer, exts) in graph.into_iter() {
//         let degree = (exts.num_exts_l() + exts.num_exts_r()) as usize;
//         nb_edges += degree;
//         degree_histo[degree] += 1;
//     }
//     nb_edges /= 2;
//     let duration = start.elapsed();
//     println!("done in {:?}", duration);
//     println!("Graph contains:\n  - {} nodes\n  - {} edges\n  degree histogramm: {:?}", nb_nodes, nb_edges, degree_histo);
// }

// /// Count the number of breakpoints in a graph, for a given haplotype
// fn stats_haplo<K: Kmer>(graph: &Graph<K>, haplo: FastaReader) {
//     print!("Iterating haplo... ");
//     std::io::stdout().flush().unwrap();
//     let start = Instant::now();
//     let mut count: usize = 0;
//     let mut junctions: usize = 0;
//     let mut dead_ends: usize = 0;

//     let mut prev_base = 0;
//     for record in haplo {
//         for kmer in record.iter_kmers(false) {
//             count += 1;
//             let exts = graph.get_exts(&kmer);
//             if exts.is_none() {
//                 println!("Kmer {} not found in graph\n", count);
//                 continue;
//             }
//             // get the next possible kmers
//             let exts = exts.unwrap();
//             match exts.num_exts_r() {
//                 0 => { dead_ends += 1; },
//                 1 => { },
//                 _ => { junctions += 1; },
//             }

//             // check that there is an edge the previous kmer
//             if !exts.has_ext(Dir::Left, prev_base) && (count > 1) {
//                 println!("Kmer {} has no edge to previous the kmer", count);
//             }
//             prev_base = kmer.get(0);
//         }
//     }
//     let duration = start.elapsed();
//     println!("done in {:?}", duration);
//     println!("Haplo contains:\n  - kmers: {}\n  - breakpoints: {}  (>1)\t\t{}  (<1)", count, junctions, dead_ends);
// }

/// Get the shortest path between first and last kmer for all sequences
fn get_shortest_path<K: Kmer>(graph: &DebruijnGraph<K,()>, fasta_reader: FastaReader, save_file: &str) {
    println!("Looking for path in graph... ");
    let mut file = BufWriter::new(File::create(save_file).unwrap());
    // let mut file = BufWriter::new(File::options().append(true).create(true).open(save_file).unwrap());       // append instead of overwrite
    let mut count: usize = 0;
    for record in fasta_reader.skip(count) {
        count += 1;
        let seq = record.sequence();
        let seq = DnaString::from_acgt_bytes(seq);
        let (first_kmer, last_kmer) = seq.both_term_kmer();
        println!("  - first kmer: {:?}\tlast kmer: {:?}", first_kmer, last_kmer);
        println!("Processing record: {}", record.header());
        let start = Instant::now();
        let path = path::get_shortest_path(graph, first_kmer, last_kmer);
        if path.is_err() {
            println!("  - no path found");
            continue;
        }
        let path = path.unwrap();
        let duration = start.elapsed();
        writeln!(file, ">{}\tlen: {}\tdone in: {:?}", record.header(), path.len(), duration).unwrap();
        writeln!(file, "{:?}", path).unwrap();
        println!("  - path length: {}\n  - time elapsed: {:?}", path.len(), duration);
        file.flush().unwrap();
        // if count == 10 {
        //     break;
        // }
    }
}

// fn get_checkpoints<K: Kmer>(graph: &Graph<K>, fasta_reader: FastaReader, save_file: &str) {
//     println!("Cutting haplo into chunks... ");
//     let mut file = BufWriter::new(File::create(save_file).unwrap());

//     let mut count: usize = 0;
//     for record in fasta_reader.skip(count) {
//         count += 1;
//         println!("Processing record: {}", record.header());
//         let start = Instant::now();
//         let path = path::get_checkpoints(&graph, &record);
//         let duration = start.elapsed();
//         println!("  - divided record into {} chunks\n  - time elapsed: {:?}", path.len(), duration);

//         // saving to file
//         writeln!(file, ">{}\tlen: {}\tdone in: {:?}", record.header(), path.len(), duration).unwrap();
//         for chunk in path {
//             let start_kmer = chunk.0;
//             let end_kmer = chunk.1;
//             writeln!(file, "{:?}\t\t{:?}", start_kmer, end_kmer).unwrap();
//         }
//         file.flush().unwrap();
//     }
// }