use crate::metadata::GraphMetadata;
use nthash::NtHashIterator;
use std::collections::HashMap;

pub struct KmerLookup {
    k: usize,
    lookup: HashMap<u64, usize>,
}

impl KmerLookup {
    /// Creates a new KmerLookup from a sequence (for testing mainly).
    /// The first occurence of a kmer in the sequence is stored in the table.
    pub fn from_seq(seq: &[u8], k: usize) -> Result<KmerLookup, Box<dyn std::error::Error>> {
        let mut lookup = HashMap::new();

        // Estimate capacity for better performance
        let estimated_capacity: usize = seq.len() - k + 1;
        lookup.reserve(estimated_capacity);

        // Populate lookup table
        let nthash_iter = NtHashIterator::new(seq, k)?;
        for (hash, i) in nthash_iter.zip(0..) {
            let old_value = lookup.insert(hash, i);
            if old_value.is_some() {
                panic!("Duplicate k-mer hash");  // TODO: handle hash collision
            }
        }

        Ok(KmerLookup {k, lookup })
    }

    /// Create a new KmerLookup from a GraphMetadata.
    /// The value of a kmer in the table correspond to its index in the GraphMetadata.
    pub fn from_graph(metadata: &GraphMetadata, k: usize) -> Result<KmerLookup, Box<dyn std::error::Error>> {
        let mut lookup = HashMap::new();

        // Estimate capacity for better performance
        let estimated_capacity: usize = metadata
            .iter()
            .map(|unitig| unitig.length as usize)
            .sum::<usize>();
        lookup.reserve(estimated_capacity);

        // Populate lookup table
        for (unitig_id, unitig) in metadata.iter().enumerate() {
            // Use NtHash for efficient canonical k-mer generation
            let nthash_iter = NtHashIterator::new(&unitig.seq, k)?;

            for hash in nthash_iter {
                let old_value = lookup.insert(hash, unitig_id);
                if old_value.is_some() {
                    panic!("Duplicate k-mer hash");  // TODO: handle hash collision
                }
            }
        }

        Ok(KmerLookup {k, lookup })
    }

    /// Lookup a kmer in the table.
    pub fn get_kmer(&self, kmer: &[u8]) -> Option<usize> {
        // note: no error checking to assert that kmer has the correct size
        let hash = nthash::ntc64(kmer, 0, self.k);
        self.lookup.get(&hash).copied()
    }

    /// Lookup all the kmers in the table.
    /// Returns a vector containing the unitig ids of all the present kmers.
    pub fn get_kmers(&self, seq: &[u8]) -> Vec<usize> {
        let mut unitig_ids = Vec::new();
        let nthash_iter = NtHashIterator::new(seq, self.k).unwrap();

        for hash in nthash_iter {
            if let Some(unitig_id) = self.lookup.get(&hash) {
                unitig_ids.push(*unitig_id);
            }
        }
        unitig_ids
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;

    #[test]
    fn test_kmer_lookup() {
        let seq = b"AAAGGGTAT";
        let k = 3;
        let lookup = KmerLookup::from_seq(seq, k).expect("k < seq.len()");

        let expected_id = 5;
        let obtained_id = lookup.get_kmer(seq[expected_id..expected_id+k].as_ref());
        assert_eq!(obtained_id, Some(expected_id));
    }

    #[test]
    fn test_kmers_lookup() {
        let seq = b"AAAGGGTAT";
        let k = 3;
        let lookup = KmerLookup::from_seq(seq, k).expect("k < seq.len()");

        let mut query = seq[..4+k].to_vec();
        query[1] = b'G';
        let expected: Vec<usize> = vec![2, 3, 4];

        let obtained = lookup.get_kmers(&query);
        assert_eq!(expected, obtained);
    }

    #[test]
    fn time_graph_construction() -> Result<(), Box<dyn std::error::Error>> {
        let path_graph = "data/graph_chr1_k31.fna";
        let k = 31;
        let kmer = b"GCCCCCCTTCATTTTCTAGCGATTCAAGCGC";
        let expected_id: usize = 4;

        // Parse metadata
        let start = Instant::now();
        let metadata = GraphMetadata::from_fasta(path_graph)?;
        let duration = start.elapsed();
        println!("  Time to parse metadata: {:?}", duration);

        // Create lookup table
        let start = Instant::now();
        let lookup = KmerLookup::from_graph(&metadata, k)?;
        let duration = start.elapsed();
        println!("  Time to create lookup table: {:?}", duration);

        // Perform lookup
        let start = Instant::now();
        let unitig_id = lookup.get_kmer(kmer).unwrap();
        let duration = start.elapsed();
        println!("  Time to lookup k-mer: {:?}", duration);
        assert_eq!(unitig_id, expected_id);
        Ok(())
    }
}
