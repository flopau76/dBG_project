use rust_dbg::encoder::{Encoder, EncoderParams};
use rust_dbg::graph::SequenceSet;
use rust_dbg::{BaseGraph, Graph, Node, NodeIterator};

use packed_seq::{PackedSeqVec, SeqVec};

fn make_compacted_graph(stranded: bool) -> (PackedSeqVec, Graph<u8>) {
    let seq = PackedSeqVec::from_ascii(b"gggccccgggaaaaac");
    let mut sequences = SequenceSet::default();

    if stranded {
        sequences.push_seq(PackedSeqVec::from_ascii(b"GGAA").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"GGG").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"CCC").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"CCGG").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"AAA").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"AAC").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"GGCC").as_slice());
    } else {
        sequences.push_seq(PackedSeqVec::from_ascii(b"GGAA").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"CCC").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"AAA").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"CCG").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"AAC").as_slice());
        sequences.push_seq(PackedSeqVec::from_ascii(b"GCC").as_slice());
    }
    let base = BaseGraph {
        sequences,
        k: 3,
        stranded,
    };
    (seq, base.finish())
}

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
