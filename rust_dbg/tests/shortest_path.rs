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
    let k = 3;
    type KS = u16;

    let seq = generate_random_seq(n);
    let seq_in = seq.as_slice().unpack();
    let seq_str = unsafe { String::from_utf8_unchecked(seq_in.clone()) };

    println!("Testing with n = {}, k = {}, stranded = {}", n, k, stranded);
    println!("seq: {:?}", seq_str);

    let base = BaseGraph::from_seq(seq.clone(), k, stranded);
    let graph: Graph<KS> = base.finish();

    let path: Vec<Node> = NodeIterator::new(&graph, seq.clone()).unwrap().collect();

    println!("path: {:?}", path);

    let seq_out = graph.path_seq(&path);

    assert_eq!(seq_in.len(), seq_out.len());
    assert_eq!(seq_in, seq_out);
}
