use debruijn::dna_string::DnaString;
use debruijn::kmer;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use rust_dbg::graph::Graph;
use rust_dbg::path::ContinuousPath;

fn generate_random_seq(n: usize) -> DnaString {
    let mut rng = rand::rng();
    let bytes: Vec<u8> = (0..n).map(|_| rng.random_range(0..=3)).collect();
    DnaString::from_bytes(&bytes)
}

#[allow(dead_code)]
fn generate_random_seq_seed(n: usize, seed: u64) -> DnaString {
    let mut rng = StdRng::seed_from_u64(seed);
    let bytes: Vec<u8> = (0..n).map(|_| rng.random_range(0..=3)).collect();
    DnaString::from_bytes(&bytes)
}

#[test]
fn main() {
    // let seed = 42;
    let n = 1000;
    let stranded = false;

    // type K = VarIntKmer<u64, K31>;
    type K = kmer::Kmer5;
    let seq = generate_random_seq(n);
    // println!("Generated sequence: {:?}", seq);
    let graph = Graph::<K>::from_seq_serial(&seq, stranded);
    graph.print_stats();

    println!("Encoding sequence into path...");
    let path = ContinuousPath::encode_seq(&graph, &seq);
    path.print_stats();

    println!("Decoding path back to sequence...");
    let seq_out = path.decode_seq(&graph);

    assert_eq!(seq, seq_out);
}
