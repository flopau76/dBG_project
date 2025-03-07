use crate::metadata::{GraphMetadata, UnitigMetadata};
use crate::lookup::KmerLookup;
use nthash::NtHashIterator;

pub struct DeBruijnGraph {
    metadata: GraphMetadata,
    lookup: KmerLookup,
    k: usize,
}

impl DeBruijnGraph {
    pub fn from_fasta(path_graph: &str, k: usize) -> Result<DeBruijnGraph, Box<dyn std::error::Error>> {
        let metadata = GraphMetadata::from_fasta(path_graph)?;
        let lookup = KmerLookup::from_graph(&metadata, k)?;
        Ok(DeBruijnGraph { metadata, lookup, k})
    }

    pub fn get_kmer(&self, kmer: &[u8]) -> Option<&UnitigMetadata> {
        let unitig_id = self.lookup.lookup(kmer)?;
        self.metadata.get(unitig_id)
    }

    pub fn walk_path(&self, seq: &[u8]) -> Result<Vec<&UnitigMetadata>, Box<dyn std::error::Error>> {
        let mut unitigs = Vec::new();
        let nthash_iter = NtHashIterator::new(seq, self.k)?;

        for hash in nthash_iter {
            if let Some(unitig) = self.get_metadata(&hash.to_ne_bytes()) {
                unitigs.push(unitig);
            } else {
                return Err("Unitig not found".into());
            }
        }
    }
}
