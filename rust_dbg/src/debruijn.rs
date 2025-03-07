use bio::io::fasta;

use debruijn;
use debruijn::dna_string::DnaString;
use debruijn::Exts;

fn get_sequences(path: &str) -> Vec<(DnaString, Exts, ())> {
    let reader = fasta::Reader::from_file(path).unwrap();
    let mut sequences: Vec<(DnaString, Exts, ())> = Vec::new();
    for record in reader.records() {
        let record = record.unwrap();
        sequences.push((DnaString::from_dna_string(dna)), Exts::empty(), ()));
    }
    sequences

}