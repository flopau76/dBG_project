use crate::metadata::GraphMetadata;
use nthash::NtHashIterator;
use std::collections::HashMap;

pub struct KmerLookup {
    k: usize,
    lookup: HashMap<u64, usize>,
}

impl KmerLookup {
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
                    return Err("Duplicate k-mer hash".into());  // TODO: handle hash collision
                }
            }
        }

        Ok(KmerLookup {k, lookup })
    }

    pub fn lookup(&self, kmer: &[u8]) -> Option<usize> {
        let hash = nthash::ntc64(kmer, 0, self.k);
        self.lookup.get(&hash).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;

    #[test]
    fn test_kmer_lookup() -> Result<(), Box<dyn std::error::Error>> {
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
        let unitig_id = lookup.lookup(kmer).unwrap();
        let duration = start.elapsed();
        println!("  Time to lookup k-mer: {:?}", duration);
        assert_eq!(unitig_id, expected_id);
        Ok(())
    }
}
