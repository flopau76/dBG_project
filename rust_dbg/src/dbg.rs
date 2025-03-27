use crate::fasta_reader::get_kmer_exts;

use debruijn::{complement, Kmer};
use debruijn::{Exts, Dir};

use boomphf::Mphf;
use serde::{Deserialize, Serialize};
use rayon::iter::{ParallelIterator, IntoParallelRefIterator, IndexedParallelIterator};

const GAMMA: f64 = 3.5;    // gamma parameter for the mphf

/// A graph composed of a minimal perfect hash function, with a list of kmers and a list of their extensions
#[derive(Debug, Serialize, Deserialize)]
pub struct Graph<K: Kmer> {
    mphf: Mphf<K>,
    kmers: Vec<K>,
    exts: Vec<Exts>,
    pub canon: bool,
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

/// Methods to construct the graph
impl<K: Kmer> Graph<K> {
    // Reorder the data according to the Mphf
    fn reorder_kmers(kmers: &mut Vec<K>, exts: &mut Vec<Exts>, mphf: &Mphf<K>) {
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

    /// Initialize a graph from a fasta
    pub fn init_from_fasta(path: &str, canon: bool) -> Self {
        let (mut kmers, mut exts) = get_kmer_exts::<K>(path, canon);
        let mphf = Mphf::new(GAMMA, &kmers);
        Self::reorder_kmers(&mut kmers, &mut exts, &mphf);
        Self {mphf, kmers, exts, canon}
    }

    /// Recompute all edges of the graph
    pub fn update_all_exts(&mut self) {
        self.exts = vec![Exts::empty(); self.len()];
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

    /// Recompute only edges for kmers where in or out degree is 0
    pub fn update_exts(&mut self) {
        // iterate over the kmers
        for (kmer_id, kmer) in self.kmers.iter().enumerate() {
            let mut new_exts = *self.get_exts(kmer_id);
            // iterate over the 8 neighbors
            for &dir in [Dir::Left, Dir::Right].iter() {
                if new_exts.num_ext_dir(dir) > 0 { continue; }
                for added_base in 0..4 {
                    let neighbor = if self.canon {kmer.extend(added_base, dir).min_rc()} else {kmer.extend(added_base, dir)};
                    // check if the neighbor is in the set
                    if self.get_key_id(&neighbor).is_some() {
                        // update the kmer
                        new_exts = new_exts.set(dir, added_base);
                    }
                }
            }
            self.exts[kmer_id] = new_exts;
        }
    }

    /// Initialize the graph from a fasta file and compute the edges
    pub fn from_fasta(path: &str, canon: bool) -> Self {
        let mut graph = Self::init_from_fasta(path, canon);
        graph.update_exts();
        graph
    }
}

/// Methods to construct the graph in parallel
impl<K: Kmer + Send + Sync> Graph<K> {
    /// Initialize a graph
    pub fn init_from_fasta_parallel(path: &str, canon: bool) -> Self {
        let (mut kmers, mut exts) = get_kmer_exts::<K>(path, canon);
        let mphf = Mphf::new_parallel(GAMMA, &kmers, None);
        Self::reorder_kmers(&mut kmers, &mut exts, &mphf);
        Self {mphf, kmers, exts, canon}
    }

    /// Recompute edges of the graph
    pub fn update_exts_parallel(&mut self) {
        // iterate over the kmers
        self.exts = self.kmers.par_iter().enumerate().map(|(kmer_id, kmer)| {
            let mut new_exts = *self.get_exts(kmer_id);
            // iterate over the 8 neighbors
            for &dir in [Dir::Left, Dir::Right].iter() {
                if new_exts.num_ext_dir(dir) > 0 { continue; }
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

    /// Initialize the graph from a fasta file and compute the edges
    pub fn from_fasta_parallel(path: &str, canon: bool) -> Self {
        let mut graph = Self::init_from_fasta_parallel(path, canon);
        graph.update_exts_parallel();
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

    const PATH : &str = "data/test.fna";
    fn _make_test_file(path: &str, seq: &str) {
        let mut writer = File::create(path).unwrap();
        writeln!(writer, ">test\n{}", seq).unwrap();
    }

    const _SEQ : &str = "CGTAA";

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
        for &canon in [false, true].iter() {
            println!("\tcanon: {}", canon);
            let expected = expected_result(canon);
            let obtained = Graph::<Kmer3>::from_fasta(PATH, canon);
            // check that the number of kmers is correct
            assert_eq!(obtained.len(), expected.len());
            for (&kmer_e, &exts_e) in expected.iter() {
                // check that the kmer is present
                let id_o = obtained.get_key_id(&kmer_e).expect("kmer not found");
                // check that it is associated with the correct edges
                let exts_o = *obtained.get_exts(id_o);
                println!("{:?}: {:?} vs {:?}", kmer_e, exts_o, exts_e);
                // assert_eq!(exts_o, exts_e);
            }
        }
    }
}

#[cfg(test)]
mod tests_time {
    use super::*;
    use debruijn::kmer;

    use std::time::Instant;

    const PATH: &str = "../data/output/chr1/graph_k31.fna";
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    fn time_from_fasta(canon: bool, parallel: bool) {
        let start = Instant::now();
        let mut graph = if parallel {Graph::<Kmer31>::init_from_fasta_parallel(PATH, canon)} else {Graph::<Kmer31>::init_from_fasta(PATH, canon)};
        let duration = start.elapsed();
        println!("Time to initialize graph: {:?}", duration);

        let start = Instant::now();
        if parallel {graph.update_exts_parallel()} else {graph.update_exts()};
        let duration = start.elapsed();
        println!("Time to update exts: {:?}", duration);
    }
        // GAMMA = 1.7 for pmhf:
    // Without parallelization: iterating only kmers, then computing all edges
    // Time to initialize graph: 591.076039135s
    // Time to update exts: 1209.704250992s

    // Without parallelization: iterating both kmers and exts, then recomputing only edges for kmers with in/out degree 0
    // Time to initialize graph: 665.664836474s
    // Time to update exts: 89.701427624s

    // With parallelization: iterating both kmers and exts, then recomputing only edges for kmers with in/out degree 0
    // Time to initialize graph: 461.868415931s
    // Time to update exts: 9.453428725s
    #[test]
    fn time_canon_para() {
        time_from_fasta(true, false);
    }


    #[test]
    fn make_graph() {
        println!("Creating graph... ");
        let start = Instant::now();
        let graph = Graph::<Kmer31>::from_fasta(PATH, true);
        let duration = start.elapsed();
        println!("done in {:?}", duration);

    }
}
