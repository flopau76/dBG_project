use crate::fasta_reader::{get_kmers, get_kmer_exts};

use debruijn::{Kmer, Exts, Dir};

use boomphf::Mphf;
use serde::{Deserialize, Serialize};
use rayon::iter::{ParallelIterator, IntoParallelRefIterator, IndexedParallelIterator};

use std::collections::HashSet;
use std::fmt::{Display, Debug, Formatter};

const GAMMA: f64 = 3.5;    // gamma parameter for the mphf

/// A graph composed of a minimal perfect hash function, with a list of kmers and a list of their extensions
#[derive(Serialize, Deserialize)]
pub struct Graph<K: Kmer> {
    mphf: Mphf<K>,
    kmers: Vec<K>,
    exts: Vec<Exts>,
    pub canon: bool,
}

impl<K: Kmer> Display for Graph<K> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        if self.canon {
            writeln!(f, "Canonical graph with {} nodes", self.kmers.len())
        } else {
            writeln!(f, "Non canonical graph with {} nodes", self.kmers.len())
        }
    }
}
impl<K: Kmer> Debug for Graph<K> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.kmers.len() {
            let kmer = &self.kmers[i];
            let exts = self.exts[i];
            write!(f, "{:?}: {:?}\n", kmer, exts)?;
        }
        Ok(())
    }
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

    /// returns the id of a kmer
    /// The queried kmer must be in the graph. Querying an absent kmer will result in a panic or a wrong value being returned.
    pub fn get_key_id_unsafe(&self, kmer: &K) -> usize {
        self.mphf.try_hash(kmer).unwrap() as usize
    }

    /// returns a ref to the kmer associated with the id
    pub fn get_kmer(&self, id: usize) -> &K {
        &self.kmers[id]
    }

    /// returns a ref to the data associated with the id
    pub fn get_exts(&self, id: usize) -> &Exts {
        &self.exts[id]
    }

    /// returns a mutable ref to the data associated to the id
    pub fn get_exts_mut(&mut self, id: usize) -> &mut Exts {
        &mut self.exts[id]
    }
}

/// Methods for graph construction
impl<K: Kmer + Send + Sync> Graph<K> {
    // Reorder the kmers according to the Mphf
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
    // Reorder the kmers and their edges according to the Mphf
    fn reorder_data(kmers: &mut Vec<K>, exts: &mut Vec<Exts>, mphf: &Mphf<K>) {
        for i in 0..kmers.len() {
            loop {
                let kmer_slot = mphf.hash(&kmers[i]) as usize;
                if i == kmer_slot {
                    break;
                }
                kmers.swap(i, kmer_slot);
                exts.swap(i, kmer_slot);
            }
        }
    }

    /// Compute the edges of the graph
    /// If all=false, only recompute edges with in/out degree = 0
    fn update_exts(&mut self, all: bool) {
        // iterate over the kmers
        self.exts = self.kmers.par_iter().enumerate().map(|(kmer_id, kmer)| {
            let mut new_exts = *self.get_exts(kmer_id);
            // iterate over the 8 neighbors
            for &dir in [Dir::Left, Dir::Right].iter() {
                if !all && new_exts.num_ext_dir(dir) > 0 { continue; }
                for added_base in 0..4 {
                    let neighbor = if self.canon {kmer.extend(added_base, dir).min_rc()} else {kmer.extend(added_base, dir)};
                    // check if the neighbor is in the set
                    if self.get_key_id(&neighbor).is_some() {
                        // update the kmer
                        new_exts = new_exts.set(dir, added_base);
                    }
                }
            }
            new_exts
        }).collect();
    }

    /// Creates a graph from a list of unique kmers.
    /// Note: the kmers must be coherent with the argument 'canon'
    pub fn from_unique_kmers(mut kmers: Vec<K>, canon: bool) -> Self {
        let mphf = Mphf::new_parallel(GAMMA, &kmers, None);
        Self::reorder_kmers(&mut kmers, &mphf);
        let exts = vec![Exts::empty(); kmers.len()];
        let mut graph = Self{mphf, kmers, exts, canon};
        graph.update_exts(true);
        graph
    }

    /// Creates a graph from a list of (potentially duplicated) kmers.
    /// Note: the kmers must be coherent with the argument 'canon'
    pub fn from_kmers(kmers: Vec<K>, canon: bool) -> Self {
        let unique_kmers: Vec<K> = kmers.into_iter().collect::<HashSet<_>>().into_iter().collect();
        Self::from_unique_kmers(unique_kmers, canon)
    }

    /// Creates a graph from a fasta file.
    pub fn from_fasta(path: &str, canon: bool) -> Self {
        let kmers = get_kmers::<K>(path, canon);
        Self::from_kmers(kmers, canon)
    }

    /// Creates a graph from a unitig file, as output by ggcat for example.
    /// Compared to from_fasta, this version does not look for branching edges within sequences
    pub fn from_unitigs(path: &str, canon: bool) -> Self {
        let (mut kmers, mut exts) = get_kmer_exts::<K>(path, canon);
        let mphf = Mphf::new_parallel(GAMMA, &kmers, None);
        Self::reorder_data(&mut kmers, &mut exts, &mphf);
        let mut graph = Self {mphf, kmers, exts, canon};
        graph.update_exts(false);
        graph
    }
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use boomphf::hashmap::BoomHashMap;

    use debruijn::Mer;
    use debruijn::kmer::Kmer3;

    use core::hash;
    use std::io::Write;
    use std::fs::File;

    const PATH : &str = "test.fna";
    const SEQ : &str = "CGTAA";
    fn make_test_file() {
        let mut writer = File::create(PATH).unwrap();
        writeln!(writer, ">test\n{}", SEQ).unwrap();
    }

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
        make_test_file();
        for &canon in [false, true].iter() {
            // println!("\tcanon: {}", canon);
            let expected = expected_result(canon);
            let obtained = Graph::<Kmer3>::from_fasta(PATH, canon);
            // check that both contain the same umber of kmers
            assert_eq!(obtained.len(), expected.len());
            for (&kmer_e, &exts_e) in expected.iter() {
                // check that each kmer is present and associated with the correct edges
                let id_o = obtained.get_key_id(&kmer_e).expect("kmer not found");
                let exts_o = *obtained.get_exts(id_o);
                assert_eq!(exts_o, exts_e);
                // println!("{:?}: {:?} vs {:?}", kmer_e, exts_o, exts_e);
            }
        }
        std::fs::remove_file(PATH).unwrap();
    }

    #[test]
    fn compare_mphf() {
        let n: usize = 300_000_000;
        let nb_queries = 100_000_000;
        let keys: Vec<usize> = (0..n).collect();

        let start = std::time::Instant::now();
        let mphf = Mphf::new(GAMMA, &keys);
        let duration = start.elapsed();
        println!("Time to create mphf: {:?}", duration);

        let start = std::time::Instant::now();
        for i in 0..nb_queries {
            let key = i % n;
            let _ = mphf.try_hash(&key);
        }
        let duration = start.elapsed();
        println!("Time to query mphf (present): {:?}", duration);

        let start = std::time::Instant::now();
        for i in 0..nb_queries {
            let key = i + n;
            let _ = mphf.try_hash(&key);
        }
        let duration = start.elapsed();
        println!("Time to query mphf (absent): {:?}", duration);

        let start = std::time::Instant::now();
        let hashset = std::collections::HashSet::<usize>::from_iter(keys.iter().cloned());
        let end = start.elapsed();
        println!("Time to create hashset: {:?}", end);

        let start = std::time::Instant::now();
        for i in 0..nb_queries {
            let key = i % n;
            let _ = hashset.contains(&key);
        }
        let duration = start.elapsed();
        println!("Time to query hashset (present): {:?}", duration);

        let start = std::time::Instant::now();
        for i in 0..nb_queries {
            let key = i + n;
            let _ = hashset.contains(&key);
        }
        let duration = start.elapsed();
        println!("Time to query hashset (absent): {:?}", duration);

    }
}

// for integers keys between 0 and 300_000_000
// total time over 100_000_000 queries
// Time to create mphf: 133.744902417s
// Time to query mphf (present): 30.222365836s
// Time to query mphf (absent): 66.332115949s
// Time to create hashset: 93.400747101s
// Time to query hashset (present): 32.621489352s
// Time to query hashset (absent): 23.975642417s


    // GAMMA = 1.7 for mphf; graph on chr1 for both haplotypes:
// Without parallelization: iterating only kmers, then computing all edges
// Time to initialize graph: 591.076039135s
// Time to update exts: 1209.704250992s

// Without parallelization: iterating both kmers and exts, then recomputing only edges for kmers with in/out degree 0
// Time to initialize graph: 665.664836474s
// Time to update exts: 89.701427624s

// With parallelization: iterating both kmers and exts, then recomputing only edges for kmers with in/out degree 0
// Time to initialize graph: 461.868415931s
// Time to update exts: 9.453428725s
