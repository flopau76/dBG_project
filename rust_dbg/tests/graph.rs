use rust_dbg::{BaseGraph, Graph, Node, NodeIterator};

use packed_seq::{PackedSeqVec, SeqVec};

#[ignore]
#[test]
fn node_iterator_offset() {
    let seq = b"aaattt";
    let base = BaseGraph::from_seq_ascii(seq, 3, true);
    let graph: Graph<u8> = base.finish();

    let seq_start = PackedSeqVec::from_ascii(b"aaa");
    let mut unitig_iter = NodeIterator::new(&graph, seq_start).unwrap();
    while let Some(_) = unitig_iter.next().unwrap() {
        // println!("Node: {:?}", node);
    }
    assert_eq!(unitig_iter.start_offset, Some(0));
    assert_eq!(unitig_iter.end_offset, Some(0));
}

#[test]
fn node_iterator() {
    let seq = PackedSeqVec::random(100);
    let seq_in = seq.as_slice().unpack();

    let base = BaseGraph::from_seq(&seq, 7, false);
    let graph: Graph<u32> = base.finish();

    let path: Vec<Node> = NodeIterator::new(&graph, seq).unwrap().collect();
    let seq_out = graph.path_seq(&path);

    assert_eq!(seq_in, seq_out);
}
