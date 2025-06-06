use rust_dbg::fasta_reader::FastaReader;
use rust_dbg::graph::Graph;
use rust_dbg::path::{MixedPath, MAX_OFFSET, MAX_PATH_LENGTH, MIN_NB_REPEATS, MIN_PATH_LENGTH};
use rust_dbg::print_progress_bar;

use debruijn::{kmer, Kmer};

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::time::Instant;

use serde::{Deserialize, Serialize};

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// K-mer size
    #[arg(short)]
    k_size: usize,

    /// Path to the binary graph file
    #[arg(short, long)]
    graph: PathBuf,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Create a graph from a file of unitigs (as generated by ggcat)
    Build {
        /// Path to the unitigs file
        #[arg(short, long)]
        input: PathBuf,

        /// Treat reverse complementary kmers as different
        #[arg(short, long, default_value_t = false)]
        forward_only: bool,
    },
    /// Get some stats about a graph
    StatsG {},
    /// Encode a continuous sequence into a path in the graph
    Encode {
        /// Path to the sequence file (fasta format)
        #[arg(short, long)]
        input: PathBuf,

        /// Path to the output file (custom text file)
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Decode a path in the graph to retrieve the original sequence
    Decode {
        /// Path to the file containing the encoded path
        #[arg(short, long)]
        input: PathBuf,

        /// Path to the output file (reconstructed sequence in fasta format)
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Get some stats about the encoding of a path
    StatsP {
        /// Path to the file containing the encoded path
        #[arg(short, long)]
        input: PathBuf,
    },
}

impl Commands {
    fn run<K: Kmer + Send + Sync + Serialize + for<'a> Deserialize<'a>>(
        &self,
        path_graph: &PathBuf,
    ) {
        match self {
            Commands::Build {
                input,
                forward_only,
            } => {
                eprint!("Creating graph... ");
                std::io::stderr().flush().unwrap();
                let start = Instant::now();
                let graph = Graph::<K>::from_unitigs_serial(input, *forward_only);
                let duration = start.elapsed();
                eprintln!("done in {:?}", duration);

                eprint!("Saving graph... ");
                std::io::stderr().flush().unwrap();
                let start = Instant::now();
                graph.save_to_binary(path_graph).unwrap();
                let duration = start.elapsed();
                eprintln!("done in {:?}", duration);
            }
            Commands::StatsG {} => {
                let graph = Graph::<K>::load_from_binary(&path_graph).unwrap();
                graph.print_stats();
            }
            Commands::Encode { input, output } => {
                let graph = Graph::<K>::load_from_binary(path_graph).unwrap();
                let fasta_reader = FastaReader::new(input).unwrap();
                let mut output_writer = File::create(output).unwrap();

                // write the header
                let date = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
                writeln!(output_writer, "DATE:\n\t{}", date).unwrap();
                let git_hash = option_env!("GIT_COMMIT_HASH").unwrap_or("unknown");
                writeln!(output_writer, "GIT VERSION:\n\t{}", git_hash).unwrap();
                let command = std::env::args().collect::<Vec<_>>().join(" ");
                writeln!(output_writer, "COMMAND:\n\t{}", command).unwrap();
                writeln!(output_writer, "CONSTANTS:").unwrap();
                writeln!(output_writer, "\tMIN_PATH_LENGTH: {}", MIN_PATH_LENGTH).unwrap();
                writeln!(output_writer, "\tMAX_PATH_LENGTH: {}", MAX_PATH_LENGTH).unwrap();
                writeln!(output_writer, "\tMIN_NB_REPEATS: {}", MIN_NB_REPEATS).unwrap();
                writeln!(output_writer, "\tMAX_OFFSET: {}", MAX_OFFSET).unwrap();

                // encode the records
                for record in fasta_reader {
                    eprintln!("Encoding record {}", record.header());
                    let dna_strings = record.dna_strings();
                    let mut current = 0;
                    let total = dna_strings.len();
                    for dna_string in dna_strings {
                        print_progress_bar(current, total);
                        current += 1;
                        let path = MixedPath::encode_seq(&graph, &dna_string);
                        writeln!(output_writer, ">{}_{}", record.header(), current).unwrap();
                        writeln!(output_writer, "{}", path).unwrap();
                    }
                    print_progress_bar(current, total);
                    eprintln!();
                }
            }
            Commands::Decode { input, output } => {
                let graph = Graph::<K>::load_from_binary(path_graph).unwrap();
                let mut input_reader = BufReader::new(File::open(input).unwrap());
                let mut output_writer = File::create(output).unwrap();

                // go to the first '>'
                input_reader.skip_until(b'>').unwrap();
                let mut header = String::new();
                while input_reader.read_line(&mut header).unwrap() > 0 {
                    // read the path
                    let mut buffer = Vec::new();
                    input_reader.read_until(b'>', &mut buffer).unwrap();
                    if buffer.last() == Some(&b'>') {
                        buffer.pop();
                    }
                    let path_str = String::from_utf8_lossy(&buffer);
                    let path = MixedPath::from_string(&path_str, &graph).unwrap();

                    // proccess the path
                    eprintln!("Decoding record {}", header);
                    let seq = path.decode_seq();
                    writeln!(output_writer, ">{}", header).unwrap();
                    writeln!(output_writer, "{}", seq.to_string()).unwrap();
                }
            }
            Commands::StatsP { input } => {
                let graph = Graph::<K>::load_from_binary(path_graph).unwrap();
                let mut input_reader = BufReader::new(File::open(input).unwrap());
                let mut paths_list = Vec::new();

                // go to the first '>'
                input_reader.skip_until(b'>').unwrap();
                let mut header = String::new();
                while input_reader.read_line(&mut header).unwrap() > 0 {
                    // read the path
                    let mut buffer = Vec::new();
                    input_reader.read_until(b'>', &mut buffer).unwrap();
                    if buffer.last() == Some(&b'>') {
                        buffer.pop();
                    }
                    let path_str = String::from_utf8_lossy(&buffer);
                    let path = MixedPath::from_string(&path_str, &graph).unwrap();
                    paths_list.push(path);
                }
                MixedPath::print_stats(&paths_list);
            }
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Hash)]
struct MyKmerSize<const K: usize>;
impl<const K: usize> kmer::KmerSize for MyKmerSize<K> {
    fn K() -> usize {
        K
    }
}

pub fn main() {
    let cli = Cli::parse();
    let k_size = cli.k_size;
    let path_graph = cli.graph;
    let command = &cli.command;

    // macro to define the generic type Kmer, whih is known only at compile time
    macro_rules! kmer_size_match {
        ($($($n:expr)+ => $tp:ty$(,)?)+) => {
            match k_size {
                $($(
                    $n => {
                        type KmerSize = MyKmerSize<$n>;
                        type Kmer = kmer::VarIntKmer<$tp, KmerSize>;
                        command.run::<Kmer>(&path_graph);
                    },
                )+)+
                _ => unimplemented!()
            }
        };
    }
    // without generic constants, we have to match kmer_size one by one
    // this is uggly and dramatically increases compilation time
    // The only workaround would be to modify the debruijn crate
    kmer_size_match!(
        1 2 3 4 => u8,
        5 6 7 8 => u16,
        9 10 11 12 13 14 15 16 => u32,
        17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 => u64,
        // ...
    );
}
