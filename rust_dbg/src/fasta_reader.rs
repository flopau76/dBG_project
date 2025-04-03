use std::io::{BufReader, BufRead, Result, Seek};
use std::fs::File;
use std::path::Path;

use debruijn::{Kmer, Exts, complement, dna_only_base_to_bits};

/// Fasta reader
pub struct FastaReader {
    buffer: BufReader<File>,
}

impl FastaReader {
    pub fn new(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self{buffer: BufReader::new(file)})
    }

    /// Position the buffer at the start of the next sequence
    /// and return the corresponding record
    pub fn next(&mut self) -> Option<FastaRecord> {
        // Skip until the next header
        match self.buffer.skip_until(b'>') {
            Ok(0) => return None, // EOF
            Ok(_) => {}, 
            Err(_) => return None, // Error
        }
        let mut header = String::new();
        let header_size = self.buffer.read_line(&mut header).unwrap() as i64;
        if header_size == 0 {
            return None;
        }
        header = header.trim_end().to_string();
        let sequence_start = self.buffer.stream_position().unwrap();
        Some(FastaRecord {
            reader: self,
            header,
            sequence_start,
        })
    }
}

/// Sequence record associated with a FastaReader
pub struct FastaRecord<'a> {
    reader: &'a mut FastaReader,
    header: String,
    sequence_start: u64,
}

impl FastaRecord<'_> {
    pub fn iter_nucleotides(&mut self) -> NucleoIterator {
        if self.reader.buffer.stream_position().unwrap() != self.sequence_start {
            self.reader.buffer.seek(std::io::SeekFrom::Start(self.sequence_start)).unwrap();
        }
        NucleoIterator::new(self)
    }

    pub fn iter_kmers<K: Kmer>(&mut self, canonical: bool) -> KmerIterator<K> {
        if self.reader.buffer.stream_position().unwrap() != self.sequence_start {
            self.reader.buffer.seek(std::io::SeekFrom::Start(self.sequence_start)).unwrap();
        }
        KmerIterator::new(self, canonical)
    }

    pub fn iter_kmer_exts<K: Kmer>(&mut self, canonical: bool) -> KmerExtsIterator<K> {
        if self.reader.buffer.stream_position().unwrap() != self.sequence_start {
            self.reader.buffer.seek(std::io::SeekFrom::Start(self.sequence_start)).unwrap();
        }
        KmerExtsIterator::new(self, canonical)
    }

    pub fn header(&self) -> &str {
        &self.header
    }

    pub fn sequence(&mut self) -> Vec<u8> {
        self.iter_nucleotides().collect()
    }
}

/// Iterator over the nucleotides of a sequence
pub struct NucleoIterator<'a> {
    reader: &'a mut FastaReader,
    buffer: Vec<u8>,
    cur_pos: usize,
}

impl<'a> NucleoIterator<'a> {
    fn new(record: &'a mut FastaRecord) -> Self {
        Self {
            reader: record.reader,
            buffer: Vec::new(),
            cur_pos: 0,
        }
    }
}

impl<'a> Iterator for NucleoIterator<'a> {
    type Item = u8;
    
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // If we've reached the end of the current line, read a new one
            if self.cur_pos >= self.buffer.len() {
                self.buffer.clear();
                self.cur_pos = 0;
                
                // Read a new line into the buffer
                let len_header = match self.reader.buffer.read_until(b'\n', &mut self.buffer) {
                    Ok(0) => return None, // EOF
                    Ok(i) => i as i64, 
                    Err(_) => return None, // Error
                };
                
                // If we reach the end of this sequence
                if !self.buffer.is_empty() && self.buffer[0] == b'>' {
                    // position the buffer at the start of the next seq and return None
                    self.reader.buffer.seek_relative(-len_header).unwrap();
                    return None;
                }
            }
            
            // Get the current byte
            if self.cur_pos < self.buffer.len() {
                let c = self.buffer[self.cur_pos];
                self.cur_pos += 1;
                
                // Skip whitespaces and newlines
                if !c.is_ascii_whitespace() {
                    return Some(c);
                }
            } else {
                self.cur_pos += 1;
            }
        }
    }
}

/// Iterator over the kmers of a sequence
/// Skips all invalid kmers (containing non-ACGT characters)
pub struct KmerIterator<'a, K: Kmer> {
    canonical: bool,
    nucleo_iter: NucleoIterator<'a>,
    cur_kmer: K,
    cur_kmer_rc: Option<K>,
    cur_count: usize,   // nb of valid bases in the current kmer
}

impl<'a, K: Kmer> KmerIterator<'a, K> {
    fn new(record: &'a mut FastaRecord, canonical: bool) -> Self {
        Self {
            canonical,
            nucleo_iter: NucleoIterator::new(record),
            cur_kmer: K::empty(),
            cur_kmer_rc: if canonical {Some(K::empty())} else {None},
            cur_count: 0,
        }
    }
}

impl<K: Kmer> Iterator for KmerIterator<'_, K> {
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

/// Iterator over the kmers and their extremities of a sequence
/// Note: Will not work properly if the sequence contains non-ACGT characters
pub struct KmerExtsIterator<'a, K: Kmer> {
    kmer_iter: KmerIterator<'a, K>,
    current_kmer: Option<(K, bool)>,
    next_kmer: Option<(K, bool)>,
}

impl<'a, K: Kmer> KmerExtsIterator<'a, K> {
    fn new(record: &'a mut FastaRecord, canonical: bool) -> Self {
        let mut kmer_iter = KmerIterator::new(record, canonical);
        let current_kmer = None;
        let next_kmer = kmer_iter.next();
        Self {
            kmer_iter,
            current_kmer,
            next_kmer,
        }
    }
}

impl<K: Kmer> Iterator for KmerExtsIterator<'_, K> {
    type Item = (K, bool, Exts);

    fn next(&mut self) -> Option<Self::Item> {
        let prev_kmer: Option<(K, bool)>;  
        (prev_kmer, self.current_kmer, self.next_kmer) = (self.current_kmer, self.next_kmer, self.kmer_iter.next());
        let current_kmer = self.current_kmer?;


        let left_ext = match prev_kmer {
            Some((kmer, flip)) => {
                let base = if !flip {kmer.get(0)} else {complement(kmer.get(K::k()-1))};
                Exts::mk_left(base)
            }
            None => Exts::empty(),
        };
        let right_ext = match self.next_kmer {
            Some((kmer, flip)) => {
                let base = if !flip {kmer.get(K::k()-1)} else {complement(kmer.get(0))};
                Exts::mk_right(base)
            }
            None => Exts::empty(),
        };
        let exts = Exts::merge(left_ext, right_ext);
        if current_kmer.1 {
            return Some((current_kmer.0, true, exts.rc()))
        } else {
            return Some((current_kmer.0, false, exts))
        }
    }
}

pub fn get_kmers<K: Kmer>(path: &str, canon: bool) -> Vec<K> {
    let mut fasta_reader = FastaReader::new(path).unwrap();
    let mut kmers = Vec::new();

    while let Some(mut record) = fasta_reader.next() {
        let kmer_iter = record.iter_kmers::<K>(canon);
        for (kmer, _flip) in kmer_iter {
            kmers.push(kmer);
        }
    }
    kmers
}

pub fn get_kmer_exts<K: Kmer>(path: &str, canon: bool) -> (Vec<K>, Vec<Exts>) {
    let mut fasta_reader = FastaReader::new(path).unwrap();
    let mut kmers = Vec::new();
    let mut exts = Vec::new();

    while let Some(mut record) = fasta_reader.next() {
        let kmer_iter = record.iter_kmer_exts::<K>(canon);
        for (kmer, _flip, ext) in kmer_iter {
            exts.push(ext);
            kmers.push(kmer);
        }
    }
    (kmers, exts)
}

#[cfg(test)]
mod unit_test {
    use super::*;

    use debruijn::kmer::Kmer3;

    static PATH : &str = "test.fna";

    // TODO: proper unit tests

    #[test]
    fn test_fasta_reader() {
        let mut fasta_reader = FastaReader::new(PATH).unwrap();
        while let Some(record) = fasta_reader.next() {
            println!("Header: {}", record.header);
        }
    }

    #[test]
    fn test_iter_nucleo() {
        let mut fasta_reader = FastaReader::new(PATH).unwrap();
        while let Some(mut record) = fasta_reader.next() {
            let iter = record.iter_nucleotides();
            let seq = String::from_utf8( iter.collect()).unwrap();
            println!("{:?}", seq);
        }
    }

    #[test]
    fn test_iter_kmers() {
        let mut fasta_reader = FastaReader::new(PATH).unwrap();
        while let Some(mut record) = fasta_reader.next() {
            println!(">{}", record.header);            
            let iter = record.iter_kmers::<Kmer3>(false);
            for kmer in iter {
                println!("{:?}", kmer);
            }
        }
    }

    #[test]
    fn test_iter_kmer_exts() {
        let mut fasta_reader = FastaReader::new(PATH).unwrap();
        while let Some(mut record) = fasta_reader.next() {
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
        let mut reader = FastaReader::new(PATH).unwrap();
        let mut count: usize = 0;
        while let Some(mut record) = reader.next() {
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
        while let Some(mut record) = reader.next() {
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