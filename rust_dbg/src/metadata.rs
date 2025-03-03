use arrayvec::ArrayVec;
use bio::io::fasta;

type Sequence = Vec<u8>;

#[derive(Debug, PartialEq)]
pub struct Edge {
    pub neighbor: usize,
    pub direction: bool,
    pub neighbor_direction: bool,
}

#[derive(Debug, PartialEq)]
pub struct UnitigMetadata {
    pub length: u32,
    pub count: u32,
    pub edges: ArrayVec<Edge, 8>,
    pub seq: Sequence,
}
impl UnitigMetadata {
    fn from_unitig_header(line: &str) -> Result<UnitigMetadata, Box<dyn std::error::Error>> {
        // Initialize values
        let mut length: u32 = 0;
        let mut count: u32 = 0;
        let mut edges: ArrayVec<Edge, 8> = ArrayVec::new();

        let parts = line.split_whitespace();

        // Parse fields
        for part in parts {
            if let Some(ln_val) = part.strip_prefix("LN:i:") {
                length = ln_val.parse()?;
            } else if let Some(kc_val) = part.strip_prefix("KC:i:") {
                count = kc_val.parse()?;
            } else if let Some(l_val) = part.strip_prefix("L:") {
                let mut edge_parts = l_val.split(':');

                let direction = match edge_parts.next() {
                    Some("+") => true,
                    Some("-") => false,
                    _ => return Err("Invalid edge direction".into()),
                };

                let neighbor = match edge_parts.next() {
                    Some(n) => n.parse()?,
                    None => return Err("Missing neighbor ID".into()),
                };

                let neighbor_direction = match edge_parts.next() {
                    Some("+") => true,
                    Some("-") => false,
                    _ => return Err("Invalid neighbor direction".into()),
                };

                edges.push(Edge {
                    neighbor,
                    direction,
                    neighbor_direction,
                });
            }
        }
        Ok(UnitigMetadata {
            length,
            count,
            edges,
            seq: Sequence::new(),
        })
    }
}

pub type GraphMetadata = Vec<UnitigMetadata>;

pub fn metadata_from_fasta(file_path: &str) -> Result<GraphMetadata, Box<dyn std::error::Error>> {
    let mut metadata = GraphMetadata::new();
    let bio_reader = fasta::Reader::from_file(file_path)?;

    let mut current_id: usize = 0;

    // Parse BCALM file
    for result in bio_reader.records() {
        // Read record
        let record = result?;
        let unitig_id: usize = record.id().parse()?;
        let unitig_desc = record.desc().ok_or(format!(
            "Invalid Bcalm: Unitig {} does not have a description",
            unitig_id
        ))?;
        let unitig_seq = record.seq();

        if current_id != unitig_id {
            return Err(format!("Invalid bcalm: Unitigs must be indexed sequentially, starting at 0. Expected {}, got {}", 
                      current_id, unitig_id).into());
        }

        // Parse record into UnitigMetadata
        let mut unitig = UnitigMetadata::from_unitig_header(&unitig_desc)?;
        unitig.seq = Vec::with_capacity(unitig.length as usize);
        unitig.seq.extend_from_slice(unitig_seq);

        // Add metadata to graph
        metadata.push(unitig);
        current_id += 1;
    }
    Ok(metadata)
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrayvec::ArrayVec;
    use std::time::Instant;

    #[test]
    fn parse_header() -> Result<(), Box<dyn std::error::Error>> {
        let line = ">0 LN:i:100 KC:i:1 L:+:1:+";
        let unitig = UnitigMetadata::from_unitig_header(line)?;
        let expected = UnitigMetadata {
            length: 100,
            count: 1,
            edges: {
                let mut edges: ArrayVec<Edge, 8> = ArrayVec::new();
                edges.push(Edge {
                    neighbor: 1,
                    direction: true,
                    neighbor_direction: true,
                });
                edges
            },
            seq: Vec::new(),
        };
        assert_eq!(unitig, expected);
        Ok(())
    }

    #[test]
    fn time_for_whole_graph() -> Result<(), Box<dyn std::error::Error>> {
        let path_graph = "../data/ggcat_output/chr1.fna";
        let start = Instant::now();
        let _metadata = metadata_from_fasta(path_graph)?;
        let duration = start.elapsed();
        println!("Time taken to parse bcalm: {:?}", duration);
        Ok(())
    }
}
