use crate::kmers::KmerIterator;

use debruijn::{complement, Kmer};
use debruijn::{Exts, Dir};

use boomphf::Mphf;
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct Graph<K: Kmer> {
    mphf: Mphf<K>,
    kmers: Vec<K>,
    exts: Vec<Exts>,
    canon: bool,
}

/// Basic methods for graph manipulation
impl<K: Kmer> Graph<K> {
    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    /// returns the id of the kmer if it exists, None otherwise
    pub fn get_key_id(&self, kmer: &K) -> Option<usize> {
        let pos = self.mphf.try_hash(kmer)?;
        let hashed_kmer = &self.kmers[pos as usize];
        if kmer == hashed_kmer {
            Some(pos as usize)
        } else {
            None
        }
    }

    /// returns a ref to the data associated to the id
    pub fn get_exts(&self, id: usize) -> &Exts {
        &self.exts[id]
    }

    /// returns a mutable ref to the data associated to the id
    pub fn get_exts_mut(&mut self, id: usize) -> &mut Exts {
        &mut self.exts[id]
    }
}

/// Methods to construct the graph
impl<K: Kmer> Graph<K> {
    // reorder the keys according to the Mphf
    fn reorder_kmers(kmers: &mut Vec<K>, mphf: &Mphf<K>) {
        for i in 0..kmers.len() {
            loop {
                let kmer_slot = mphf.hash(&kmers[i]) as usize;
                if i == kmer_slot {
                    break;
                }
                kmers.swap(i, kmer_slot);
            }
        }
    }

    /// initialize the graph from a fasta file
    pub fn init_from_fasta(path: &str, canon: bool) -> Self {
        let mut kmers: Vec<K> = KmerIterator::<K>::new(path, canon).unwrap().map(|x|x.0).collect();
        let capacity = kmers.len();
        let mphf = Mphf::new(1.7, &kmers); // TODO: look at the paper for correct gamma factor
        Self::reorder_kmers(&mut kmers, &mphf);
        Self {mphf, kmers, exts: vec![Exts::empty(); capacity], canon}
    }

    /// compute all edges of a stranded de Bruijn graph
    pub fn update_exts(&mut self) {
        // iterate over the kmers
        for (kmer_id, kmer) in self.kmers.iter().enumerate() {
            let mut new_exts = *self.get_exts(kmer_id);
            // iterate over the 8 neighbors
            for &dir in [Dir::Left, Dir::Right].iter() {
                let removed_base = match dir {
                    Dir::Left => kmer.get(K::k()-1),
                    Dir::Right => kmer.get(0)
                };
                for added_base in 0..4 {
                    let (neighbor, flip) = if self.canon {kmer.extend(added_base, dir).min_rc_flip()} else {(kmer.extend(added_base, dir), false)};
                    // only lookup once each pair (kmer, neighbor)
                    if *kmer < neighbor { continue; }
                    if *kmer == neighbor { new_exts = new_exts.set(dir, added_base); continue; }
                    // check if the neighbor is in the set
                    if let Some(neighbor_id) = self.get_key_id(&neighbor) {
                        // update the kmer
                        new_exts = new_exts.set(dir, added_base);
                        // update the neighbor (beware if it is reversed)
                        let neighbor_exts = &mut self.exts[neighbor_id];
                        if !flip { *neighbor_exts = neighbor_exts.set(dir.flip(), removed_base); }
                        else { *neighbor_exts = neighbor_exts.set(dir, complement(removed_base)); }
                    }
                }
            }
            self.exts[kmer_id] = new_exts;
        }
    }

    /// initialize the graph from a fasta file and compute the edges
    pub fn from_fasta(path: &str, canon: bool) -> Self {
        let mut graph = Self::init_from_fasta(path, canon);
        graph.update_exts();
        graph
    }
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use boomphf::hashmap::BoomHashMap;

    use debruijn::Mer;
    use debruijn::kmer::Kmer3;

    use std::io::Write;
    use std::fs::File;

    // TODO: create the test file at beginning, then delete it at the end (beware of test parallelisation !)

    static PATH : &str = "data/test.fna";
    fn make_test_file(path: &str, seq: &str) {
        let mut writer = File::create(path).unwrap();
        writeln!(writer, ">test\n{}", seq).unwrap();
    }

    static SEQ : &str = "CGTAA";

    fn expected_result(canon: bool) -> BoomHashMap<Kmer3, Exts> {
        let kmers = vec![
            if !canon {Kmer3::from_ascii(b"CGT")} else {Kmer3::from_ascii(b"CGT").rc() },
            Kmer3::from_ascii(b"GTA"),
            Kmer3::from_ascii(b"TAA"),
        ];

        let exts = if !canon {
            vec![
                Exts::empty().set(Dir::Right, 0),
                Exts::empty().set(Dir::Left, 1).set(Dir::Right, 0),
                Exts::empty().set(Dir::Left, 2)
            ]
        } else {
            vec![
                Exts::empty().set(Dir::Right, 0).rc()                   .set(Dir::Right, 3),
                Exts::empty().set(Dir::Left, 1).set(Dir::Right, 0)      .set(Dir::Right, 1),
                Exts::empty().set(Dir::Left, 2)                         .set(Dir::Left, 3)
            ]
        };

        BoomHashMap::new(kmers, exts)
    }

    #[test]
    fn test_from_fasta() {
        for &canon in [true, false].iter() {
            let expected = expected_result(canon);
            let obtained = Graph::<Kmer3>::from_fasta(PATH, canon);
            assert_eq!(obtained.len(), expected.len());
            for (&kmer_e, &exts_e) in expected.iter() {
                // check that the kmer is present
                let id_o = obtained.get_key_id(&kmer_e).expect("kmer not found");
                // check that it is associated with the correct edges
                let exts_o = *obtained.get_exts(id_o);
                assert_eq!(exts_o, exts_e);
            }
        }
    }
}

#[cfg(test)]
mod tests_time {
    use super::*;
    use debruijn::kmer;

    use std::time::Instant;

    static PATH: &str = "data/unitigs_chr1_k31.fna";
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    fn time_from_fasta(canon: bool) {
        let start = Instant::now();
        let mut graph = Graph::<Kmer31>::init_from_fasta(PATH, canon);
        let duration = start.elapsed();
        println!("Time to initialize graph: {:?}", duration);

        let start = Instant::now();
        graph.update_exts();
        let duration = start.elapsed();
        println!("Time to update exts: {:?}", duration);
    }

    // Time to initialize graph: 576.703851699s
    // Time to update exts: 1192.2858803s
    #[test]
    fn time_no_canon() {
        time_from_fasta(false);
    }

    // Time to initialize graph: 591.076039135s
    // Time to update exts: 1209.704250992s
    #[test]
    fn time_canon() {
        time_from_fasta(true);
    }
}
