pub use bio::io::fasta;
pub use debruijn::{Exts, Dir};
pub use debruijn::{kmer, Kmer, Vmer, DnaBytes};

pub use std::collections::{HashSet, HashMap};
pub use boomphf::hashmap::BoomHashMap;


// const K: usize = 3;

// #[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
// pub struct SizeK;

// impl kmer::KmerSize for SizeK {
//     #[inline]
//     fn K() -> usize {
//         K
//     }
// }

// type KmerK = kmer::VarIntKmer<u64, SizeK>;
pub type KmerK = kmer::Kmer3;

fn get_seq_iterator(path: &str) -> impl Iterator< Item = DnaBytes > {
    let reader = fasta::Reader::from_file(path).unwrap();
    let records = reader.records();
    records.map(|record| DnaBytes(record.unwrap().seq().to_owned()))
}

// extract the kmers from a fasta file
fn get_kmers(seq_iter: impl Iterator<Item = DnaBytes>) -> HashSet<KmerK> {
    // Use a HashSet to automatically handle deduplication
    seq_iter
        .flat_map(|seq| seq.iter_kmers().collect::<Vec<KmerK>>())
        .collect()
}

// compute the edges of the de Bruijn graph
fn get_kmer_exts(kmer_set: HashSet<KmerK>) -> HashMap<KmerK, Exts> {
    let mut kmer_map: HashMap<KmerK, Exts> = HashMap::with_capacity(kmer_set.len());
    for kmer in &kmer_set {
        let mut exts = Exts::new(0);
        for i in 0..4 {
            let left_ext = kmer.extend_left(i);
            if kmer_set.contains(&left_ext) {
                exts = exts.set(Dir::Left, i);
            }
            let right_ext = kmer.extend_right(i);
            if kmer_set.contains(&right_ext) {
                exts = exts.set(Dir::Right, i);
            }
        }
        kmer_map.insert(*kmer, exts);
    }
    kmer_map
}


#[cfg(test)]
mod tests {
    use super::*;

    fn make_seq_iterator() -> impl Iterator<Item = DnaBytes> {
        let dna_vec = vec![DnaBytes(b"CGGTAAA".to_vec()), DnaBytes(b"ACGRAACCGGTT".to_vec())];
        dna_vec.into_iter()
    }

    #[test]
    fn test_kmer_set() {
        let seq_iter = make_seq_iterator();
        let kmer_set = get_kmers(seq_iter);
        let expected_kmer_set: HashSet<KmerK> = vec![
            KmerK::from_bytes(b"CGG"),
            KmerK::from_bytes(b"GGT"),
            KmerK::from_bytes(b"GTA"),
            KmerK::from_bytes(b"TAA"),
            KmerK::from_bytes(b"AAA"),
            KmerK::from_bytes(b"ACG"),
            KmerK::from_bytes(b"CGR"),
            KmerK::from_bytes(b"GRA"),
            KmerK::from_bytes(b"RAA"),
            KmerK::from_bytes(b"AAC"),
            KmerK::from_bytes(b"ACC"),
            KmerK::from_bytes(b"CCG"),
            KmerK::from_bytes(b"CGT"),
            KmerK::from_bytes(b"GTT"),
        ].into_iter().collect();
        assert_eq!(kmer_set, expected_kmer_set);
    }
}
