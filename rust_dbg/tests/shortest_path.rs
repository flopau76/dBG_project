use rust_dbg::{BaseGraph, Graph, Node, NodeIterator};

use packed_seq::{PackedSeqVec, SeqVec};

fn generate_random_seq(n: usize) -> PackedSeqVec {
    PackedSeqVec::random(n)
}

#[test]
fn main() {
    // let seed = 42;
    let n = 50;
    let stranded = false;
    let k = 7;
    type KS = u16;

    let seq = generate_random_seq(n);
    let base = BaseGraph::from_seq(seq.clone(), k, stranded);
    let graph: Graph<KS> = base.finish();

    let path: Vec<Node> = NodeIterator::new(&graph, seq.clone()).unwrap().collect();

    let seq_out = graph.path_seq(&path);

    assert_eq!(seq.as_slice().unpack(), seq_out.as_slice().unpack());
}
