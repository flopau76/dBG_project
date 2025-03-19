mod dbg;
mod kmers;

use debruijn::{kmer, Kmer};
use dbg::Graph;

use std::fs::File;
use std::hash::RandomState;
use std::io::{BufWriter, BufReader, Write};

use serde::{Serialize, Deserialize};
use bincode::serde::{encode_into_std_write, decode_from_std_read};

use std::time::Instant;

fn main() {
    let path_graph = "data/unitigs_chr1_k31.fna";
    let path_haplo = "data/AalbF5_chr1.fna";
    let path_bin = "data/graph_k31.bin";
    let path_test = "data/test.fna";

    let canon = true;
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    // let graph = make_graph::<Kmer31>(path_graph, canon);
    // save_graph(graph, path_bin);

    // let graph = load_graph::<Kmer31>(path_bin);
    // println!("Nb of kmers: {}", graph.len());
    // let mut nb_edges: usize = 0;
    // for i in 0..graph.len() {
    //     let exts = *graph.get_exts(i);
    //     nb_edges += (exts.num_exts_l() + exts.num_exts_r()) as usize;
    // }
    // println!("Nb of edges: {}", nb_edges);

    let kmers: Vec<Kmer31> = dbg::get_kmers(path_haplo, false);
    let kmers_canon: Vec<Kmer31> = dbg::get_kmers(path_haplo, true);
    println!("Total kmers: {}", kmers.len());
    let kmers_unique: std::collections::HashSet<&kmer::VarIntKmer<u64, kmer::K31>, RandomState> = std::collections::HashSet::from_iter(kmers.iter());
    let kmers_unique_canon: std::collections::HashSet<&kmer::VarIntKmer<u64, kmer::K31>, RandomState> = std::collections::HashSet::from_iter(kmers_canon.iter());
    println!("Nb non canon kmers: {}", kmers_unique.len());
    println!("Nb canon kmers: {}", kmers_unique_canon.len());

    // walk_graph(kmers, graph);

}

fn make_graph<K: Kmer>(path: &str, canon: bool) -> Graph<K> {
    println!("Creating graph");
    let start = Instant::now();
    let mut graph = dbg::Graph::<K>::from_fasta(path, canon);
    let duration = start.elapsed();
    println!("Time to create graph: {:?}", duration);

    graph
}

fn save_graph<K: Kmer + Serialize>(graph: Graph<K>, path_output: &str) {
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

// Time to read graph: 81.801018011s
// Time to walk graph: 92.534646519s
// Edges left: [251, 140836558, 13031541, 7215599, 7765715]
// Edges right: [240, 141292470, 12877727, 7031160, 7648067]
// Not found: 0
fn walk_graph<K: Kmer>(kmers: Vec<K> , graph: Graph<K>) {
    println!("Walking graph");
    let start = Instant::now();
    let mut edges_l = [0;5];
    let mut edges_r = [0;5];
    let mut not_found = 0;
    let mut tot = Vec::with_capacity(kmers.len());
    for kmer in kmers {
        let id = graph.get_key_id(&kmer);
        if id.is_none() {   // should not happen
            not_found += 1;
            tot.push(0);
            continue;
        }
        let id = id.unwrap();
        let exts = *graph.get_exts(id);
        edges_l[exts.num_exts_l() as usize] += 1;
        edges_r[exts.num_exts_r() as usize] += 1;
        tot.push(exts.num_exts_l() + exts.num_exts_r());
    }
    let duration = start.elapsed();
    println!("Time to walk graph: {:?}", duration);
    println!("Edges left: {:?}", edges_l);
    println!("Edges right: {:?}", edges_r);
    println!("Not found: {:?}", not_found);

    let output_path = "data/walk_output.txt";
    let mut file = BufWriter::new(File::create(output_path).unwrap());
    for value in tot {
        writeln!(file, "{}", value).unwrap();
    }
    println!("Saved tot to {}", output_path);
}