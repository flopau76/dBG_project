//! To construct and manipulate De Bruijn graphs.

use crate::fasta_reader::{get_kmers, get_kmer_exts};

use debruijn::{Kmer, Exts, Dir};

use serde::{Deserialize, Serialize};
use rayon::iter::{ParallelIterator, IntoParallelRefIterator};

use std::collections::HashMap;
use std::fmt::{Display, Debug, Formatter};

/// A graph allowing to query the presence of a kmer, and its neighbors.
#[derive(Debug, Serialize, Deserialize)]
pub struct Graph<K: Kmer> {
    map: HashMap<K, Exts>,
    /// Wether two reverse complement kmers are considered equal or not.
    pub canon: bool,
}

impl<K: Kmer> Display for Graph<K> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        if self.canon {
            writeln!(f, "Canonical graph with {} nodes", self.map.len())
        } else {
            writeln!(f, "Non canonical graph with {} nodes", self.map.len())
        }
    }
}

/// Basic methods for graph manipulation
impl<K: Kmer> Graph<K> {
    /// Number of kmers in the graph.
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Check if the graph contains the given kmer.
    /// If self.canon is true, the kmer is first converted to its canonical form.
    pub fn contains_kmer(&self, kmer: &K) -> bool {
        let kmer = if self.canon {kmer.min_rc()} else {*kmer};
        self.map.contains_key(&kmer)
    }

    /// Returns the edges of the given kmer if present in the graph, None otherwise.
    /// If self.canon is true, the kmer is first converted to its canonical form, but the returned edges to correspond to its original orientation.
    pub fn get_exts(&self, kmer: &K) -> Option<Exts> {
        let (kmer_c, flip) = if self.canon {kmer.min_rc_flip()} else {(*kmer, false)};
        let exts = *self.map.get(&kmer_c)?;
        if flip {
            Some(exts.rc())
        } else {
            Some(exts)
        }
    }

    fn iter(&self) -> std::collections::hash_map::Iter<'_, K, Exts> {
        self.map.iter()
    }
}

impl<'a, K: Kmer> IntoIterator for &'a Graph<K> {
    type Item = (&'a K, &'a Exts);
    type IntoIter = std::collections::hash_map::Iter<'a, K, Exts>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// Methods for graph construction
impl<K: Kmer + Send + Sync> Graph<K> {
    
    /// Compute the edges of the graph
    /// If all=false, only recompute edges with in/out degree = 0
    fn update_exts(&mut self, all: bool) {
        // Create a vector of (K, Exts) pairs to update
        let updates: Vec<(K, Exts)> = self.map.par_iter()
            .filter_map(|(kmer, ext)| {
                let mut new_exts = ext.clone();
                let mut changed = false;
                // iterate over the 8 neighbors
                for &dir in [Dir::Left, Dir::Right].iter() {
                    if !all && new_exts.num_ext_dir(dir) > 0 { continue; }
                    for added_base in 0..4 {
                        let neighbor = if self.canon {
                            kmer.extend(added_base, dir).min_rc()
                        } else {
                            kmer.extend(added_base, dir)
                        };
                        // check if the neighbor is in the set
                        if self.map.contains_key(&neighbor) {
                            // update the kmer
                            new_exts = new_exts.set(dir, added_base);
                            changed = true;
                        }
                    }
                }
                
                // Only return pairs that actually changed
                if changed {
                    Some((kmer.clone(), new_exts))
                } else {
                    None
                }
            })
            .collect();
        
        // Apply all the updates to the original HashMap
        for (kmer, new_exts) in updates {
            if let Some(ext) = self.map.get_mut(&kmer) {
                *ext = new_exts;
            }
        }
    }

    /// Creates a graph from a list of (potentially duplicated) kmers.
    /// Note: the kmers must be coherent with the argument 'canon'.
    pub fn from_kmers(kmers: Vec<K>, canon: bool) -> Self {
        let map = kmers.into_iter().map(|kmer| (kmer, Exts::empty())).collect();
        let mut graph = Self{map, canon};
        graph.update_exts(true);
        graph
    }

    /// Creates a graph from a fasta file.
    pub fn from_fasta(path: &str, canon: bool) -> Self {
        let kmers = get_kmers::<K>(path, canon);
        Self::from_kmers(kmers, canon)
    }

    /// Creates a graph from a unitig file, as output by ggcat for example.
    /// Compared to `from_fasta`, this version does not look for branching edges within sequences.
    pub fn from_unitigs(path: &str, canon: bool) -> Self {
        let (kmers, exts) = get_kmer_exts::<K>(path, canon);
        let map = kmers.iter().zip(exts.iter()).map(|(kmer, ext)| (kmer.clone(), ext.clone())).collect();
        let mut graph = Self{map, canon};
        graph.update_exts(false);
        graph
    }
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use debruijn::Mer;
    use debruijn::kmer::Kmer3;

    use std::io::Write;
    use std::fs::File;

    const PATH : &str = "test.fna";
    const SEQ : &str = "CGTAA";
    fn make_test_file() {
        let mut writer = File::create(PATH).unwrap();
        writeln!(writer, ">test\n{}", SEQ).unwrap();
    }

    fn expected_result(canon: bool) -> HashMap<Kmer3, Exts> {
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
        kmers.into_iter().zip(exts.into_iter()).collect()
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
            for (kmer_e, exts_e) in expected.iter() {
                // check that each kmer is present and associated with the correct edges
                let exts_o = obtained.map.get(kmer_e).expect("kmer not found");
                assert_eq!(exts_o, exts_e);
                // println!("{:?}: {:?} vs {:?}", kmer_e, exts_o, exts_e);
            }
        }
        std::fs::remove_file(PATH).unwrap();
    }
}