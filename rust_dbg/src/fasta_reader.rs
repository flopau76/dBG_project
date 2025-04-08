use std::io::{BufReader, BufRead, Result, Seek, SeekFrom};
use std::fs::File;
use std::path::Path;
use std::slice::Iter;

use debruijn::{Kmer, Exts, complement, dna_only_base_to_bits};

use rayon::prelude::*;
use std::sync::Arc;
use std::cell::RefCell;

/// Fasta reader
pub struct FastaReader {
    buffer: BufReader<File>,
    index: Vec<usize>,
    path: Arc<Path>,
}

impl FastaReader {
    pub fn new(path: impl AsRef<Path>) -> Result<Self> {
        let path = Arc::from(path.as_ref());

        let file = File::open(&path)?;
        let mut buffer = BufReader::new(file);
        
        // Iterate once over the file to make the index
        let mut index = Vec::new();
        let mut pos: usize = 0;
        let mut nb_bytes = buffer.skip_until(b'>').unwrap();
        while nb_bytes > 0 {
            pos += nb_bytes;
            index.push(pos);
            nb_bytes = buffer.read_until(b'>', &mut Vec::new()).unwrap();
        }

        // Rewind the buffer to the beginning
        buffer.seek(SeekFrom::Start(0)).unwrap();

        Ok(Self{buffer, index, path})
    }
    
    /// Reads the next fasta record from the file
    pub fn next(&mut self) -> Option<FastaRecord> {
        // Skip until the next header
        match self.buffer.skip_until(b'>') {
            Ok(0) => return None, // EOF
            Ok(_) => {}, 
            Err(_) => return None, // Error
        }
        // Read the header line
        let mut header = String::new();
        let header_size = self.buffer.read_line(&mut header).unwrap();
        if header_size == 0 {
            return None;
        }
        header = header.trim_end().to_string();

        // Read the sequence
        let mut sequence = Vec::new();
        self.buffer.read_until(b'>', &mut sequence).unwrap();
        if let Some(last) = sequence.last() {
            if *last == b'>' {sequence.pop();}
        }
        sequence.retain(|&c| c != b'\n');

        // Reposition the buffer at the start of the next seq header
        self.buffer.seek_relative(-1).unwrap();

        Some(FastaRecord {
            header,
            sequence,
        })
    }
}

impl Iterator for FastaReader {
    type Item = FastaRecord;
    
    fn next(&mut self) -> Option<Self::Item> {
        self.next()
    }
}

/// Parallel iterator for FastaReader
impl FastaReader {
       pub fn par_iter(self) -> impl ParallelIterator<Item = FastaRecord> {
        let path = self.path;
        let index = self.index[0..(self.index.len()-1)].to_vec();

        index.into_par_iter().map_init(
            // Initialize a thread-local file handle
            move || {
                let file = File::open(&*path).unwrap();
                RefCell::new(BufReader::new(file))
            },
            // Use the thread-local file handle for each chunk
            |reader, start| {
                // Get the file reader and position it at the start of the record
                let mut reader = reader.borrow_mut();
                reader.seek(SeekFrom::Start(start as u64)).unwrap();
                
                // Read the header line
                let mut header = String::new();
                reader.read_line(&mut header).unwrap();
                let header = header.trim_end().to_string();
                
                // Read the sequence
                let mut sequence = Vec::new();
                reader.read_until(b'>', &mut sequence).unwrap();
                if let Some(last) = sequence.last() {
                    if *last == b'>' {sequence.pop();}
                }
                sequence.retain(|&c| c != b'\n');
                                
                FastaRecord { header, sequence }
            }
        )
    }
}

/// Sequence record associated with a FastaReader
pub struct FastaRecord {
    header: String,
    sequence: Vec<u8>,
}

impl FastaRecord {
    pub fn iter_nucleotides(&self) -> Iter<'_, u8> {
        self.sequence.iter()
    }

    pub fn iter_kmers<K: Kmer>(&self, canonical: bool) -> impl Iterator<Item = (K, bool)> {
        KmerIterator::new(self, canonical)
    }

    pub fn iter_kmer_exts<K: Kmer>(&self, canonical: bool) -> impl Iterator<Item = (K, bool, Exts)> {
        KmerExtsIterator::new(self, canonical)
    }

    pub fn header(&self) -> &str {
        &self.header
    }

    pub fn sequence(&mut self) -> &Vec<u8> {
        &self.sequence
    }
}

/// Iterator over the kmers of a sequence
/// Skips all invalid kmers (containing non-ACGT characters)
pub struct KmerIterator<'a, K: Kmer> {
    nucleo_iter: Iter<'a, u8>,
    cur_kmer: K,
    cur_kmer_rc: Option<K>,
    cur_count: usize,   // nb of valid bases in the current kmer
}

impl<'a, K: Kmer> KmerIterator<'a, K> {
    fn new(record: &'a FastaRecord, canonical: bool) -> Self {
        Self {
            nucleo_iter: record.iter_nucleotides(),
            cur_kmer: K::empty(),
            cur_kmer_rc: if canonical {Some(K::empty())} else {None},
            cur_count: 0,
        }
    }
}

impl<K: Kmer> Iterator for KmerIterator<'_, K> {
    type Item = (K, bool);

    fn next(&mut self) -> Option<Self::Item> {
        // elongate the current kmer until it reaches k
        while self.cur_count < K::k() {
            let &base = self.nucleo_iter.next()?;
            // convert the base to a 2-bit representation
            match dna_only_base_to_bits(base) {
                Some(nuc) => {
                    self.cur_kmer = self.cur_kmer.extend_right(nuc);
                    if let Some(rc) = self.cur_kmer_rc.as_mut() {
                        *rc = rc.extend_left(complement(nuc));
                    }
                    self.cur_count += 1;
                }
                // If we encounter a non-ACGT base (N), skip and reset the current kmer to 0
                None => {
                    self.cur_kmer = K::empty();
                    if let Some(rc) = self.cur_kmer_rc.as_mut() {
                        *rc = K::empty();
                    }
                    self.cur_count = 0;
                }
            }
        }
        self.cur_count -= 1;
        if let Some(rc)  = self.cur_kmer_rc {
            if self.cur_kmer < rc {
                return Some((self.cur_kmer, false));
            } else {
                return Some((rc, true));
            };
        }
        return Some((self.cur_kmer, false));
    }
}

/// Iterator over the kmers and their extremities of a sequence
/// Note: Will not work properly if the sequence contains non-ACGT characters
pub struct KmerExtsIterator<'a, K: Kmer> {
    kmer_iter: KmerIterator<'a, K>,
    left_ext: Exts,
    current_kmer: Option<(K, bool)>,
}

impl<'a, K: Kmer> KmerExtsIterator<'a, K> {
    fn new(record: &'a FastaRecord, canonical: bool) -> Self {
        let mut kmer_iter = KmerIterator::new(record, canonical);
        let left_ext = Exts::empty();
        let current_kmer = kmer_iter.next();
        Self {
            kmer_iter,
            left_ext,
            current_kmer,
        }
    }
}

impl<K: Kmer> Iterator for KmerExtsIterator<'_, K> {
    type Item = (K, bool, Exts);

    fn next(&mut self) -> Option<Self::Item> {
        let (current_kmer, current_flip) = self.current_kmer?;

        // get the rightmost base of the next kmer
        let next_kmer = self.kmer_iter.next();
        let right_ext = match next_kmer {
            Some((next_kmer, next_flip)) => {
                let right_base = if !next_flip {next_kmer.get(K::k()-1)} else {complement(next_kmer.get(0))};
                Exts::mk_right(right_base)
            }
            None => Exts::empty(),
        };

        // merge the left and right extensions
        let mut exts = Exts::merge(self.left_ext, right_ext);
        if current_flip {exts = exts.rc();}

        // update the left_ext for the next kmer
        let left_base = if !current_flip {current_kmer.get(0)} else {complement(current_kmer.get(K::k()-1))};
        self.left_ext = Exts::mk_left(left_base);
        self.current_kmer = next_kmer;
        
        Some((current_kmer, current_flip, exts))
    }
}

pub fn get_kmers<K: Kmer>(path: &str, canon: bool) -> Vec<K> {
    let fasta_reader = FastaReader::new(path).unwrap();
    let mut kmers = Vec::new();

    for record in fasta_reader {
        let kmer_iter = record.iter_kmers::<K>(canon);
        kmers.extend(kmer_iter.map(|(kmer, _flip)| kmer));
    }
    kmers
}

pub fn get_kmer_exts<K: Kmer>(path: &str, canon: bool) -> (Vec<K>, Vec<Exts>) {
    let fasta_reader = FastaReader::new(path).unwrap();
    let mut kmers = Vec::new();

    for record in fasta_reader {
        let kmer_iter = record.iter_kmer_exts::<K>(canon);
        kmers.extend(kmer_iter.map(|(kmer, _flip, exts)| (kmer, exts)));
    }
    kmers.into_iter().unzip()
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use debruijn::kmer::Kmer3;

    static PATH : &str = "../data/input/test.fna";

    // TODO: proper unit tests

    #[test]
    fn test_fasta_reader() {
        let fasta_reader = FastaReader::new(PATH).unwrap();
        for record in fasta_reader {
            println!("Header: {}", record.header);
            println!("Sequence: {}", String::from_utf8_lossy(&record.sequence));
        }
    }

    #[test]
    fn test_par_iter_records() {
        let fasta_reader = FastaReader::new(PATH).unwrap();
        let par_iter = fasta_reader.par_iter();
        par_iter.for_each(|record| {
            println!("Header: {}", record.header);
            println!("Sequence: {}", String::from_utf8_lossy(&record.sequence));
        });

    }

    #[test]
    fn test_iter_kmers() {
        let fasta_reader = FastaReader::new(PATH).unwrap();
        for record in fasta_reader {
            println!(">{}", record.header);            
            let iter = record.iter_kmers::<Kmer3>(false);
            for kmer in iter {
                println!("{:?}", kmer);
            }
        }
    }

    #[test]
    fn test_iter_kmer_exts() {
        let fasta_reader = FastaReader::new(PATH).unwrap();
        for record in fasta_reader {
            println!(">{}", record.header);            
            let iter = record.iter_kmer_exts::<Kmer3>(true);
            for kmer in iter {
                println!("{:?}", kmer);
            }
        }
    }
}

#[cfg(test)]
mod tests_time {
    use super::*;
    use debruijn::kmer;
    use std::time::Instant;

    static PATH: &str = "../data/input/chr1/AalbF5.fna";
    type Kmer31 = kmer::VarIntKmer<u64, kmer::K31>;

    // Time to iter kmers: 16.733445004s
    // Nb of kmers: 337698822
    #[test]
    fn iter_kmers() {
        let start = Instant::now();
        let mut count: usize = 0;
        let fasta_reader = FastaReader::new(PATH).unwrap();
        for record in fasta_reader {
            let iter = record.iter_kmers::<Kmer31>(false);
            for _kmer in iter {
                count += 1;
            }
        }
        let duration = start.elapsed();
        println!("Time to iter kmers: {:?}", duration);
        println!("Nb of kmers: {}", count);
    }
    
    // Time to iter kmers with extensions: 29.078663587s
    // Nb of kmers: 337698822
    #[test]
    fn iter_kmer_exts() {
        let start = Instant::now();
        let mut reader = FastaReader::new(PATH).unwrap();
        let mut count: usize = 0;
        while let Some(record) = reader.next() {
            let iter = record.iter_kmer_exts::<Kmer31>(false);
            for _kmer in iter {
                count += 1;
            }
        }
        let duration = start.elapsed();
        println!("Time to iter kmers with extensions: {:?}", duration);
        println!("Nb of kmers: {}", count);
    }
}