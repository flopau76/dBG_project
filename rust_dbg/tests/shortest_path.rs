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

#[test]
fn easy_seq_encoding() {
    let (seq, graph) = make_compacted_graph(false);

    let encoder_params = EncoderParams {
        min_nb_repeats: 10 as u16,
        min_sp_length: 10,
        max_sp_length: 50,
        max_offset: 255,
    };
    let encoder = Encoder {
        params: encoder_params,
        graph: &graph,
    };

    let path_in: Vec<Node> = NodeIterator::new(&graph, seq.clone()).unwrap().collect();
    let encoding = encoder.encode_path(&path_in);
    let path_out = encoding.decode(&graph);

    println!("encoding: {:?}", encoding);

    assert_eq!(path_in.len(), path_out.len());
    assert_eq!(path_in, path_out);
}

#[test]
fn seq_encoding() {
    let seq = PackedSeqVec::random(100);

    let base = BaseGraph::from_seq(&seq, 5, false);
    let graph: Graph<u32> = base.finish();

    let encoder_params = EncoderParams {
        min_nb_repeats: 5 as u16,
        min_sp_length: 5,
        max_sp_length: 50,
        max_offset: 255,
    };
    let encoder = Encoder {
        params: encoder_params,
        graph: &graph,
    };

    let path_in: Vec<Node> = NodeIterator::new(&graph, seq.clone()).unwrap().collect();
    let encoding = encoder.encode_path(&path_in);
    let path_out = encoding.decode(&graph);

    assert_eq!(path_in.len(), path_out.len());
    // assert_eq!(path_in, path_out);
}
