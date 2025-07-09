use rust_dbg::embeddings::{Embedding, VecExtensions};
use rust_dbg::encoder::GreedyEncoder;
use rust_dbg::graph::SequenceSet;
use rust_dbg::{BaseGraph, Graph};

use packed_seq::{PackedSeqVec, SeqVec};

fn make_compacted_graph(stranded: bool) -> (Vec<u8>, Graph<u8>) {
    let seq = b"gggccccgggaaaaac".to_ascii_uppercase().to_vec();

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
fn easy_seq_encoding() {
    for stranded in [true, false] {
        let (seq_in, graph) = make_compacted_graph(stranded);

        let encoder_params = GreedyEncoder {
            min_nb_repeats: 10 as u16,
            min_sp_length: 1,
            max_sp_length: 50,
            max_offset: 255,
        };

        let encoding = VecExtensions::from_seq(&seq_in, &graph, &encoder_params).unwrap();
        let seq_out = encoding.get_seq(&graph);

        assert_eq!(seq_in, seq_out);
    }
}

#[test]
fn random_seq_encoding() {
    let seq_in = PackedSeqVec::random(10000);

    let base = BaseGraph::from_seq(&seq_in, 5, false);
    let graph: Graph<u32> = base.finish();

    let encoder = GreedyEncoder {
        min_nb_repeats: 500 as u16,
        min_sp_length: 1,
        max_sp_length: 50,
        max_offset: 255,
    };

    let seq_in = seq_in.as_slice().unpack();

    let encoding = VecExtensions::from_seq(&seq_in, &graph, &encoder).unwrap();
    let seq_out = encoding.get_seq(&graph);

    assert_eq!(seq_in, seq_out);
}
