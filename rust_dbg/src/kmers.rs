use std::io::{ BufReader, BufRead, Result};
use std::fs::{File};
use std::path::Path;

use debruijn::{complement, Kmer, dna_only_base_to_bits};
use debruijn::{Exts, Dir};

struct NucleoIterator {
    reader: BufReader<File>,
    buffer: Vec<u8>,
    cur_pos: usize,
}

impl NucleoIterator {
    fn new(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self {
            reader: BufReader::new(file),
            buffer: Vec::new(),
            cur_pos: 0,
        })
    }
}

impl Iterator for NucleoIterator {
    type Item = u8;
    
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // If we've reached the end of the current line, read a new one
            if self.cur_pos >= self.buffer.len() {
                self.buffer.clear();
                self.cur_pos = 0;
                
                // Read a new line into the buffer
                match self.reader.read_until(b'\n', &mut self.buffer) {
                    Ok(0) => return None, // EOF
                    Ok(_) => {}, 
                    Err(_) => return None, // Error
                }
                
                // Check if it is a header line
                if !self.buffer.is_empty() && self.buffer[0] == b'>' {
                    // Return a single char '>' and skip the rest of the line
                    self.cur_pos = self.buffer.len();
                    return Some(b'>');
                }
            }
            
            // Get the current byte
            if self.cur_pos < self.buffer.len() {
                let c = self.buffer[self.cur_pos];
                self.cur_pos += 1;
                
                // Skip whitespace and newlines
                if !c.is_ascii_whitespace() {
                    return Some(c);
                }
            } else {
                self.cur_pos += 1;
            }
        }
    }
}

pub struct KmerIterator<K: Kmer> {
    canonical: bool,
    nucleo_iter: NucleoIterator,
    cur_kmer: K,
    cur_kmer_rc: Option<K>,
    cur_count: usize,
}

impl<K: Kmer> KmerIterator<K> {
    fn new(path: impl AsRef<Path>, canonical: bool) -> Result<Self> {
        let nucleo_iter = NucleoIterator::new(path)?;
        Ok(Self {
            canonical,
            nucleo_iter,
            cur_kmer: K::empty(),
            cur_kmer_rc: if canonical { Some(K::empty()) } else { None },
            cur_count: 0,
        })
    }
}

impl<K: Kmer> Iterator for KmerIterator<K> {
    type Item = (K, bool);

    fn next(&mut self) -> Option<Self::Item> {
        // Elongate the current kmer until it reaches k
        while self.cur_count < K::k() {
            let base = self.nucleo_iter.next()?; {
            // convert the base to a 2-bit representation
            match dna_only_base_to_bits(base) {
                Some(nuc) => {
                    self.cur_kmer = self.cur_kmer.extend_right(nuc);
                    if let Some(rc) = self.cur_kmer_rc.as_mut() {
                        *rc = rc.extend_left(complement(nuc));
                    }
                    self.cur_count += 1;
                }
                // If we encounter a non-ACGT char (N or >), reset the current kmer
                None => {
                    self.cur_kmer = K::empty();
                    if let Some(rc) = self.cur_kmer_rc.as_mut() {
                        *rc = K::empty();
                    }
                    self.cur_count = 0;
                }
            }
            }
        }
        self.cur_count -= 1;
        if self.canonical {
            let rc = self.cur_kmer_rc.expect("Safe because canonical was computed");
            let (kmer, flip) = if self.cur_kmer < rc {
                (self.cur_kmer, false)
            } else {
                (rc, true)
            };
            Some((kmer, flip))
        } else {
            Some((self.cur_kmer, false))
        }
    }
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use debruijn::kmer::Kmer3;

    static PATH : &str = "data/test.fna";

    #[test]
    fn test_nucleo() {
        let iter = NucleoIterator::new(PATH).unwrap();
        let seq: Vec<char> = iter.map(|b  |{b as char}).collect();
        println!("{:?}", seq);
    }
    #[test]
    fn test_kmers() {
        let iter = KmerIterator::<Kmer3>::new(PATH, false).unwrap();
        let kmers: Vec<Kmer3> = iter.map(|b| b.0).collect();
        println!("{:?}", kmers);
    }
}

#[cfg(test)]
mod tests_time {
    use super::*;
    use debruijn::kmer;
    use std::time::Instant;

    static PATH: &str = "data/AalbF5_chr1.fna";
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    // On laptop:
    // Time to walk graph: 42.117445043s
    // Nb of kmers: 337_698_822
    #[test]
    fn walk_graph() {
        let start = Instant::now();
        let kmer_iter = KmerIterator::<Kmer31>::new(PATH, false).unwrap();
        let mut nb_kmers = 0;
        for kmer in kmer_iter {
            nb_kmers += 1;
        }
        let duration = start.elapsed();
        println!("Time to walk graph: {:?}", duration);
        println!("Nb of kmers: {}", nb_kmers);

    }
}