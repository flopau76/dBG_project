use bio::io::fasta;
use needletail::{parse_fastx_file, Sequence, FastxReader};

use debruijn::dna_string::DnaString;
use debruijn::{kmer, Kmer, Vmer, Mer};
use debruijn::{Exts, Dir};

use std::collections::HashSet;
use boomphf::hashmap::BoomHashMap2;

pub fn get_kmers_needletail<K: Kmer>(path: &str, canon: bool) -> HashSet<K> {
    let mut reader = parse_fastx_file(path).unwrap();
    let mut kmers = HashSet::new();

    while let Some(record) =  reader.next() {
        let seq = record.unwrap();
        for (k1, (kmer,k2), rc) in seq.bit_kmers(K::k() as u8, canon) {
            kmers.insert(K::from_u64(kmer));
        }
    }
    kmers
}

// compute the edges of a stranded de Bruijn graph
fn get_kmer_exts<K: Kmer>(kmer_set: HashSet<K>, canon: bool) -> BoomHashMap2<K, Exts, ()> {
    let nb_kmers = kmer_set.len();
    let mut exts_vec: Vec<Exts> = Vec::with_capacity(nb_kmers);

    let mut count: usize = 0;
    for kmer in &kmer_set {
        count += 1;
        if count % 100000 == 0 {
            println!("Processed {} kmers", count);
        }
        let mut new_exts = Exts::empty();
        for i in 0..4 {
            let left_ext = if canon {kmer.extend_left(i).min_rc()} else {kmer.extend_left(i)};
            if kmer_set.contains(&left_ext) {
                new_exts = new_exts.set(Dir::Left, i);
            }
            let right_ext = if canon {kmer.extend_right(i).min_rc()} else {kmer.extend_right(i)};
            if kmer_set.contains(&right_ext) {
                new_exts = new_exts.set(Dir::Right, i);
            }
        }
        exts_vec.push(new_exts);
    }
    let kmer_vec: Vec<K> = kmer_set.into_iter().collect();
    BoomHashMap2::new(kmer_vec, exts_vec, vec![(); nb_kmers])
}


#[cfg(test)]
mod unit_test {
    use super::*;

    use kmer::Kmer3;

    use std::io::Write;
    use std::fs::File;

    fn make_test_file(path: &str, seq: &str) {
        let mut writer = File::create(path).unwrap();
        writeln!(writer, ">test\n{}", seq).unwrap();
    }

    static PATH : &str = "data/test.fna";
    static SEQ : &str = "ACCGAT";
    fn expected_kmers(canon: bool) -> HashSet<Kmer3> {
        let mut kmers = HashSet::new();
        kmers.insert(Kmer3::from_ascii(b"ACC"));
        kmers.insert(Kmer3::from_ascii(b"CCG"));
        kmers.insert(Kmer3::from_ascii(b"CGA"));
        kmers.insert(if canon {Kmer3::from_ascii(b"GAT").rc()} else {Kmer3::from_ascii(b"GAT")});
        kmers
    }


    #[test]
    fn test_bio_non_canon() {
        let kmers = get_kmers::<Kmer3>(PATH, false);
        let expected = expected_kmers(false);
        assert_eq!(kmers, expected);
    }

    #[test]
    fn test_bio_canon() {
        let kmers = get_kmers::<Kmer3>(PATH, true);
        let expected = expected_kmers(true);
        assert_eq!(kmers, expected);
    }

    #[test]
    fn test_needletail_non_canon() {
        let kmers = get_kmers_needletail::<Kmer3>(PATH, false);
        let expected = expected_kmers(false);
        assert_eq!(kmers, expected);
    }

    #[test]
    fn test_needletail_canon() {
        let kmers = get_kmers_needletail::<Kmer3>(PATH, true);
        let expected = expected_kmers(true);
        assert_eq!(kmers, expected);
    }
}

#[cfg(test)]
mod tests_time {
    use super::*;

    use std::time::{Instant, Duration};

    static PATH: &str = "data/unitigs_chr1_k31.fna";
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    fn time_closure<F: Fn()->()> (f: F) -> Duration
    {
        let start = Instant::now();
        f();
        start.elapsed()
    }

    #[test]
    fn time_non_canon() {
        let duration = time_closure(||-> () {get_kmers::<Kmer31>(PATH, false);});
        println!("Time to get non canonical kmers: {:?}", duration);
    }

    #[test]
    fn time_canon() {
        let duration = time_closure(||-> () {get_kmers::<Kmer31>(PATH, true);});
        println!("Time to get canonical kmers: {:?}", duration);
    }

    #[test]
    fn time_needletail_non_canon() {
        let duration = time_closure(||-> () {get_kmers_needletail::<Kmer31>(PATH, false);});
        println!("Time to get canonical kmers with needletail: {:?}", duration);
    }

    #[test]
    fn time_needletail_canon() {
        let duration = time_closure(||-> () {get_kmers_needletail::<Kmer31>(PATH, true);});
        println!("Time to get canonical kmers with needletail: {:?}", duration);
    }
}
