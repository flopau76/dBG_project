use bincode::{Decode, Encode};
use derive_more::with_trait::{Add, AddAssign, Deref, DerefMut, Sum};
use std::fmt::Display;

use crate::{Graph, KmerStorage, Node};

use super::{Encoder, Encoding, EncodingStats};

#[derive(Default, Debug, Deref, DerefMut, Encode, Decode)]
pub struct VecNodes(pub Vec<Node>);

impl Encoding for VecNodes {
    type Stats = BasicStats;

    fn decode(&self, _graph: &Graph<impl KmerStorage>) -> Vec<Node> {
        self.0.clone()
    }
    fn get_encoding_stats(&self, _graph: &Graph<impl KmerStorage>) -> Self::Stats {
        BasicStats {
            num_nodes: self.0.len(),
        }
    }
}

#[derive(Default, Debug, Add, AddAssign, Sum)]
pub struct BasicStats {
    num_nodes: usize,
}
impl Display for BasicStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Number of nodes: {}", self.num_nodes)
    }
}
impl EncodingStats for BasicStats {}

#[derive(Default)]
pub struct BasicEncoder;
impl Encoder for BasicEncoder {
    type Encoding = VecNodes;
    fn encode_path(&self, nodes: Vec<Node>, _graph: &Graph<impl KmerStorage>) -> VecNodes {
        VecNodes(nodes)
    }
}
