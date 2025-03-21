use needletail::{parse_fastx_file, Sequence};

use debruijn::{complement, Kmer};
use debruijn::{Exts, Dir};

use boomphf::hashmap::BoomHashMap;
use serde::{Deserialize, Serialize};

use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

/// get all kmers of a fasta file
pub fn get_kmers<K: Kmer>(path: &str, canon: bool) -> Vec<K> {
    let mut reader = parse_fastx_file(path).unwrap();
    let mut kmers = Vec::new();

    while let Some(record) =  reader.next() {
        let seq = record.unwrap();
        for (_k1, (kmer, _k2), _rc) in seq.bit_kmers(K::k() as u8, canon) {
            kmers.push(K::from_u64(kmer));
        }
    }
    kmers
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Graph<K: Kmer> {
    kmers: BoomHashMap<K, ()>,
    exts: Vec<Exts>,
    canon: bool,
}

impl<K: Kmer> Graph<K> {
    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    /// returns the id of the kmer if it exists, None otherwise
    pub fn get_key_id(&self, kmer: &K) -> Option<usize> {
        self.kmers.get_key_id(kmer)
    }

    /// returns a ref to the data associated to the id
    pub fn get_exts(&self, id: usize) -> &Exts {
        &self.exts[id]
    }

    /// returns a mutable ref to the data associated to the id
    pub fn get_exts_mut(&mut self, id: usize) -> &mut Exts {
        &mut self.exts[id]
    }

    // compute all edges of a stranded de Bruijn graph
    fn compute_exts(kmers: &BoomHashMap<K, ()>, canon: bool) -> Vec<Exts> {
        let mut exts = vec![Exts::empty(); kmers.len()];
        // iterate over the kmers
        for (kmer_id, (kmer, _)) in kmers.iter().enumerate() {
            let mut new_exts = exts[kmer_id];
            // iterate over the 8 neighbors
            for &dir in [Dir::Left, Dir::Right].iter() {
                let removed_base = match dir {
                    Dir::Left => kmer.get(K::k()-1),
                    Dir::Right => kmer.get(0)
                };
                for added_base in 0..4 {
                    let (neighbor, flip) = if canon {kmer.extend(added_base, dir).min_rc_flip()} else {(kmer.extend(added_base, dir), false)};
                    // only lookup once each pair (kmer, neighbor)
                    if *kmer < neighbor { continue; }
                    if *kmer == neighbor { new_exts = new_exts.set(dir, added_base); continue; }
                    // check if the neighbor is in the set
                    if let Some(neighbor_id) = kmers.get_key_id(&neighbor) {
                        // update the kmer
                        new_exts = new_exts.set(dir, added_base);
                        // update the neighbor (beware if it is reversed)
                        let mut neighbor_exts = &mut exts[neighbor_id];
                        if !flip { *neighbor_exts = neighbor_exts.set(dir.flip(), removed_base); }
                        else { *neighbor_exts = neighbor_exts.set(dir, complement(removed_base)); }
                    }
                }
            }
            exts[kmer_id] = new_exts;
        }
        exts
    }

    /// initialize the graph from a fasta file and compute the edges
    pub fn from_fasta(path: &str, canon: bool) -> Self {
        let kmers_vec = get_kmers(path, canon);
        let nb_kmers = kmers_vec.len();
        let kmers = BoomHashMap::new(kmers_vec, vec![(); nb_kmers]);
        let exts = Self::compute_exts(&kmers, canon);
        Self {kmers, exts, canon}
    }
}

// parallel implementation
impl <K:Kmer + Send + Sync> Graph<K> {
        /// Calculates the extensions of the kmers
        fn compute_exts_parallel(kmers: &BoomHashMap<K, ()>, canon: bool) -> Vec<Exts> {
            // iterate (in parallel) over the kmers and compute the extensions
            let new_exts: Vec<_> = kmers.par_iter()
                .map(|(kmer, _)| {
                    let mut new_exts = Exts::empty();
                    
                    // iterate over the 8 neighbors
                    for &dir in [Dir::Left, Dir::Right].iter() {
                        for added_base in 0..4 {
                            let neighbor = if canon { kmer.extend(added_base, dir).min_rc()}
                                else { kmer.extend(added_base, dir) };
                            // check if the neighbor is in the set
                            if kmers.get_key_id(&neighbor).is_some() {
                                new_exts = new_exts.set(dir, added_base);
                            }
                        }
                    }
                    new_exts
                })
                .collect();
            new_exts
        }


    /// initialize the graph from a fasta file and compute the edges
    pub fn from_fasta_parallel(path: &str, canon: bool) -> Self {
        let kmers_vec = get_kmers(path, canon);
        let nb_kmers = kmers_vec.len();
        let kmers = BoomHashMap::new_parallel(kmers_vec, vec![(); nb_kmers]);
        let exts = Self::compute_exts_parallel(&kmers, canon);
        Self {kmers, exts, canon}
    }
}

#[cfg(test)]
mod unit_test {
    use super::*;

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
            let mut obtained = Graph::<Kmer3>::from_fasta(PATH, canon);
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

    use std::time::{Instant, Duration};

    static PATH: &str = "data/unitigs_chr1_k31.fna";
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    fn time_closure<F: Fn()->()> (f: F) -> Duration
    {
        let start = Instant::now();
        f();
        start.elapsed()
    }

    // // 622s
    // #[test]
    // fn time_init_non_canon() {
    //     let duration = time_closure(||-> () {Graph::<Kmer31>::init_from_fasta(PATH, false);});
    //     println!("Time to get non canonical kmers: {:?}", duration);
    // }

    // // 627s
    // #[test]
    // fn time_init_canon() {
    //     let duration = time_closure(||-> () {Graph::<Kmer31>::init_from_fasta(PATH, true);});
    //     println!("Time to get canonical kmers: {:?}", duration);
    // }

    // // 1860s
    // #[test]
    // fn time_non_canon() {
    //     let duration = time_closure(||-> () {Graph::<Kmer31>::from_fasta(PATH, false);});
    //     println!("Time to get non canonical graph: {:?}", duration);
    // }

    // // 1890s
    // #[test]
    // fn time_canon() {
    //     let duration = time_closure(||-> () {Graph::<Kmer31>::from_fasta(PATH, true);});
    //     println!("Time to get canonical graph: {:?}", duration);
    // }
}
