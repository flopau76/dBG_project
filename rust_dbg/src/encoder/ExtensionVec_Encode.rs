use super::{Extension, ExtensionVec, Node};

enum Patch {
    TargetNode(Vec<Node>),
    NextNucleotide(Vec<u8>),
    Repetition(Vec<(u16, u8)>),
}

impl Patch {
    fn new(ext: Extension) -> Self {
        match ext {
            Extension::TargetNode(node) => Patch::TargetNode(vec![node]),
            Extension::NextNucleotide(nuc) => Patch::NextNucleotide(vec![nuc]),
            Extension::Repetition { nb_repeats, offset } => {
                Patch::Repetition(vec![(nb_repeats, offset)])
            }
        }
    }
    fn push(&mut self, ext: Extension) {
        match (self, ext) {
            (Patch::TargetNode(nodes), Extension::TargetNode(node)) => nodes.push(node),
            (Patch::NextNucleotide(nucleotides), Extension::NextNucleotide(nuc)) => {
                nucleotides.push(nuc)
            }
            (Patch::Repetition(reps), Extension::Repetition { nb_repeats, offset }) => {
                reps.push((nb_repeats, offset))
            }
            _ => panic!("Incompatible extension type for patch"),
        }
    }
    fn len(&self) -> usize {
        match self {
            Patch::TargetNode(nodes) => nodes.len(),
            Patch::NextNucleotide(nucleotides) => nucleotides.len(),
            Patch::Repetition(reps) => reps.len(),
        }
    }

    fn from_extension_vec(extensions: ExtensionVec) -> Vec<Self> {
        let extensions = extensions.0;
        let mut patches = Vec::new();
        let mut current_patch = Patch::new(extensions[0]);
        let mut current_type = std::mem::discriminant(&extensions[0]);
        for ext in extensions.iter().skip(1) {
            let new_type = std::mem::discriminant(ext);
            if new_type == current_type && current_patch.len() < 63 {
                current_patch.push(*ext);
            } else {
                patches.push(current_patch);
                current_patch = Patch::new(*ext);
                current_type = new_type;
            }
        }
        patches.push(current_patch);
        patches
    }
    fn to_extension_vec(patches: Vec<Self>) -> ExtensionVec {
        let mut extensions = Vec::new();
        for patch in patches.into_iter() {
            match patch {
                Patch::TargetNode(nodes) => {
                    for node in nodes {
                        extensions.push(Extension::TargetNode(node));
                    }
                }
                Patch::NextNucleotide(nucleotides) => {
                    for nuc in nucleotides {
                        extensions.push(Extension::NextNucleotide(nuc));
                    }
                }
                Patch::Repetition(reps) => {
                    for (nb_repeats, offset) in reps {
                        extensions.push(Extension::Repetition { nb_repeats, offset });
                    }
                }
            }
        }
        ExtensionVec(extensions)
    }
}

// pack the two bit encoded bases into bytes
fn pack_bases(bases: &[u8]) -> Vec<u8> {
    let mut packed = Vec::with_capacity(bases.len().div_ceil(4));
    let mut current_byte = 0u8;
    for (i, &base) in bases.iter().enumerate() {
        current_byte |= base << (2 * (i % 4));
        if (i + 1) % 4 == 0 || i == bases.len() - 1 {
            packed.push(current_byte);
            current_byte = 0;
        }
    }
    packed
}

// unpack the two bit encoded bases from bytes
fn unpack_bases(bases: &[u8]) -> Vec<u8> {
    let mut unpacked = Vec::with_capacity(bases.len() * 4);
    for &byte in bases {
        for i in 0..4 {
            let base = (byte >> (2 * i)) & 0b11;
            unpacked.push(base);
        }
    }
    unpacked
}

use bincode::error::{AllowedEnumVariants, DecodeError, EncodeError};
use bincode::{
    de::{BorrowDecoder, Decoder},
    enc::{write::Writer, Encoder},
};
use bincode::{BorrowDecode, Decode, Encode};

impl Encode for Patch {
    fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        assert!(
            self.len() < (1 << 6),
            "Patch length exceeds maximum of 6 bits"
        );
        match self {
            Self::TargetNode(nodes) => {
                let id = (nodes.len() as u8) << 2 | 0;
                Encode::encode(&id, encoder)?;
                for node in nodes {
                    Encode::encode(node, encoder)?;
                }
                Ok(())
            }
            Self::NextNucleotide(nucleotides) => {
                let id = (nucleotides.len() as u8) << 2 | 1;
                Encode::encode(&id, encoder)?;
                let packed = pack_bases(nucleotides);
                encoder.writer().write(&packed)?;
                Ok(())
            }
            Self::Repetition(repetitions) => {
                let id = (repetitions.len() as u8) << 2 | 2;
                Encode::encode(&id, encoder)?;
                for (nb_repeats, offset) in repetitions {
                    Encode::encode(&(nb_repeats, offset), encoder)?;
                }
                Ok(())
            }
        }
    }
}

impl<Context> Decode<Context> for Patch {
    fn decode<D: Decoder>(decoder: &mut D) -> Result<Self, DecodeError> {
        let id: u8 = Decode::decode(decoder)?;
        let len = (id >> 2) as usize;
        match id & 0b11 {
            0 => {
                let mut nodes = Vec::with_capacity(len);
                for _ in 0..len {
                    nodes.push(Decode::decode(decoder)?);
                }
                Ok(Patch::TargetNode(nodes))
            }
            1 => {
                let mut packed_nucleotides: Vec<u8> = Vec::with_capacity(len.div_ceil(4));
                for _ in 0..len.div_ceil(4) {
                    packed_nucleotides.push(Decode::decode(decoder)?);
                }
                let mut nucleotides = unpack_bases(&packed_nucleotides);
                nucleotides.truncate(len);
                Ok(Patch::NextNucleotide(nucleotides))
            }
            2 => {
                let mut repetitions = Vec::with_capacity(len);
                for _ in 0..len {
                    repetitions.push(Decode::decode(decoder)?);
                }
                Ok(Patch::Repetition(repetitions))
            }
            other => Err(DecodeError::UnexpectedVariant {
                type_name: "Patch",
                found: other as u32,
                allowed: &AllowedEnumVariants::Range { min: 0, max: 2 },
            }),
        }
    }
}
impl<'d, Context> BorrowDecode<'d, Context> for Patch {
    fn borrow_decode<D: BorrowDecoder<'d>>(decoder: &mut D) -> Result<Self, DecodeError> {
        Ok(Decode::decode(decoder)?)
    }
}

impl Encode for ExtensionVec {
    fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        let patches = Patch::from_extension_vec(self.clone());
        Encode::encode(&patches, encoder)
    }
}
impl<Context> Decode<Context> for ExtensionVec {
    fn decode<D: Decoder>(decoder: &mut D) -> Result<Self, DecodeError> {
        let patches: Vec<Patch> = Decode::decode(decoder)?;
        Ok(Patch::to_extension_vec(patches))
    }
}
impl<'d, Context> BorrowDecode<'d, Context> for ExtensionVec {
    fn borrow_decode<D: BorrowDecoder<'d>>(decoder: &mut D) -> Result<Self, DecodeError> {
        Ok(Decode::decode(decoder)?)
    }
}
