use crate::metadata::{GraphMetadata, UnitigMetadata};
use crate::lookup::KmerLookup;

pub struct DeBruijnGraph {
    metadata: GraphMetadata,
    lookup: KmerLookup,
}

impl DeBruijnGraph {
    pub fn from_fasta(path_graph: &str, k: usize) -> Result<DeBruijnGraph, Box<dyn std::error::Error>> {
        let metadata = GraphMetadata::from_fasta(path_graph)?;
        let lookup = KmerLookup::from_graph(&metadata, k)?;
        Ok(DeBruijnGraph { metadata, lookup})
    }

    pub fn get_kmer(&self, kmer: &[u8]) -> Option<&UnitigMetadata> {
        let unitig_id = self.lookup.get_kmer(kmer)?;
        self.metadata.get(unitig_id)
    }

    pub fn get_kmers(&self, kmer: &[u8]) -> Vec<&UnitigMetadata> {
        let unitig_ids = self.lookup.get_kmers(kmer);
        unitig_ids.iter().filter_map(|&unitig_id| self.metadata.get(unitig_id)).collect()
    }
}