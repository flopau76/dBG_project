// pub mod lookup;
// pub mod metadata;
// pub mod dbg;

pub mod debruijn;
mod debruijn_example;

use bio::io::fasta;

fn main() {
    // let path_graph = "data/graph_chr1_k31.fna";
    // let path_haplo = "data/AalbF5_chr1.fna";
    // let k = 31;

    // let dbg = dbg::DeBruijnGraph::from_fasta(path_graph, k).unwrap();
    // let records = read_fasta(path_haplo);
    // let haplo = records[0].as_slice();

    // let unitigs = dbg.get_kmers(haplo);
    // let edges = get_breakpoints(unitigs);
    // let mut count = 0;
    // for edge in edges {
    //     if edge == 1 {
    //         count += 1;
    //     } else {
    //         println!("{}:{}", edge, count);
    //         count = 0;
    //     }
    // }
}

// fn get_breakpoints(unitigs: Vec<&metadata::UnitigMetadata>) -> Vec<usize> {
//     unitigs.iter().map(|unitig| unitig.edges.len()/2).collect()
// }

// fn read_fasta(path: &str) -> Vec<Vec<u8>> {
//     let reader = fasta::Reader::from_file(path).unwrap();
//     let mut records = Vec::new();
//     for record in reader.records() {
//         let record = record.unwrap();
//         records.push(record.seq().to_ascii_uppercase());
//     }
//     records
// }
