//! To parse fasta files and get the DNA records.

use std::fs::File;
use std::io::{BufRead, BufReader, Result, Seek, SeekFrom};
use std::path::Path;

use debruijn::dna_string::DnaString;

/// Buffered reader for fasta files.
pub struct FastaReader {
    buffer: BufReader<File>,
}

impl FastaReader {
    /// Creates a new FastaReader from the path to a fasta file.  
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::open(&path)?;
        let mut buffer = BufReader::new(file);

        // Iterate once over the file to make the index
        let mut index = Vec::new();
        let mut pos: usize = 0;
        let mut nb_bytes = buffer.skip_until(b'>')?;
        while nb_bytes > 0 {
            pos += nb_bytes;
            index.push(pos);
            nb_bytes = buffer.read_until(b'>', &mut Vec::new())?;
        }

        // Rewind the buffer to the beginning
        buffer.seek(SeekFrom::Start(0))?;

        Ok(Self { buffer })
    }

    /// Reads the next record in the file.
    fn next(&mut self) -> Option<DnaRecord> {
        // Skip until the next header
        match self.buffer.skip_until(b'>') {
            Ok(0) => return None, // EOF
            Ok(_) => {}
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
            if *last == b'>' {
                sequence.pop();
            }
        }
        sequence.retain(|&c| c != b'\n');

        // Reposition the buffer at the start of the next seq header
        self.buffer.seek_relative(-1).unwrap();

        Some(DnaRecord { header, sequence })
    }
}

impl Iterator for FastaReader {
    type Item = DnaRecord;

    fn next(&mut self) -> Option<Self::Item> {
        self.next()
    }
}

/// DNA record composed of a header and the corresponding sequence.
pub struct DnaRecord {
    header: String,
    sequence: Vec<u8>,
}

impl DnaRecord {
    /// Creates a new DnaRecord from a header and a sequence.
    pub fn new(header: String, sequence: Vec<u8>) -> Self {
        Self { header, sequence }
    }
    /// Returns the header of the record.
    pub fn header(&self) -> &str {
        &self.header
    }
    /// Returns the sequence of the record (in ascii encoding).
    pub fn sequence(&self) -> &Vec<u8> {
        &self.sequence
    }
    /// Returns the sequence of the record as a DnaString. Non-acgt char will be converted to 'A'.
    pub fn dna_string(&self) -> DnaString {
        DnaString::from_acgt_bytes(&self.sequence)
    }
    /// Returns a list of DnaString, separated by 'N' in the sequence.
    pub fn dna_strings(&self) -> Vec<DnaString> {
        self.sequence
            .split(|&c| c == b'N')
            .filter(|chunk| !chunk.is_empty())
            .map(DnaString::from_acgt_bytes)
            .collect()
    }
}
