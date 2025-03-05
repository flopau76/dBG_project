use crate::metadata::{GraphMetadata, UnitigMetadata};
use crate::lookup::KmerLookup;

pub struct DeBruijnGraph {
    metadata: GraphMetadata,
    lookup: KmerLookup,
}

impl DeBruijnGraph {
    pub fn new(path_graph: &str, k: usize) -> Result<DeBruijnGraph, Box<dyn std::error::Error>> {
        let metadata = GraphMetadata::from_fasta(path_graph)?;
        let lookup = KmerLookup::from_graph(&metadata, k)?;
        Ok(DeBruijnGraph { metadata, lookup })
    }

    pub fn get_metadata(&self, kmer: &[u8]) -> Option<&UnitigMetadata> {
        let unitig_id = self.lookup.lookup(kmer)?;
        self.metadata.get(unitig_id)
    }
}
