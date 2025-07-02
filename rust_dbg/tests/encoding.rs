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
fn easy_seq_encoding() {
    for stranded in [true, false] {
        let (seq, graph) = make_compacted_graph(stranded);

        let encoder_params = EncoderParams {
            min_nb_repeats: 10 as u16,
            min_sp_length: 1,
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
        assert_eq!(path_in, path_out);
    }
}

#[test]
fn seq_encoding() {
    for stranded in [true, false] {
        let seq = PackedSeqVec::from_ascii(b"gggcggtggcggaaa");
        let graph = Graph::<u8>::from_seq(&seq, 3, stranded);

        let encoder_params = EncoderParams {
            min_nb_repeats: 10 as u16,
            min_sp_length: 1,
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

        println!("{:?}", path_out);

        assert_eq!(path_in.len(), path_out.len());
        assert_eq!(path_in, path_out);
    }
}

#[test]
fn random_seq_encoding() {
    let seq = PackedSeqVec::random(1000);

    let base = BaseGraph::from_seq(&seq, 5, false);
    let graph: Graph<u32> = base.finish();

    let encoder_params = EncoderParams {
        min_nb_repeats: 500 as u16,
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
