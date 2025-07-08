use crate::embeddings::{Extension, VecExtensions};
use crate::Node;

use bincode::error::{AllowedEnumVariants, DecodeError, EncodeError};
use bincode::{
    de::{BorrowDecoder, Decoder},
    enc::{write::Writer, Encoder},
};
use bincode::{BorrowDecode, Decode, Encode};

//####################################################################################
//                             Serialize  Node                                      //
//####################################################################################

impl Encode for Node {
    fn encode<__E: Encoder>(&self, encoder: &mut __E) -> Result<(), EncodeError> {
        assert!(self.id < (1 << 63), "Node ID exceeds maximum of 63 bits");
        let value = ((self.id as u64) << 1) | (self.is_rc as u64);
        let value = vint64::encode(value);
        encoder.writer().write(value.as_ref())?;
        Ok(())
    }
}

impl<Context> Decode<Context> for Node {
    fn decode<D: Decoder>(decoder: &mut D) -> Result<Self, DecodeError> {
        let prefix: u8 = Decode::decode(decoder)?;
        let len = vint64::decoded_len(prefix);
        let mut buffer = vec![0; len];
        buffer[0] = prefix;
        for i in 1..len {
            buffer[i] = Decode::decode(decoder)?;
        }
        let value = vint64::decode(&mut buffer.as_slice())
            .map_err(|_| DecodeError::Other("Failed to parse var int in Node encoding"))?;
        let id = (value >> 1) as usize;
        let is_rc = (value & 1) != 0;
        Ok(Self { id, is_rc })
    }
}
impl<'d, Context> BorrowDecode<'d, Context> for Node {
    fn borrow_decode<D: BorrowDecoder<'d>>(decoder: &mut D) -> Result<Self, DecodeError> {
        Ok(Decode::decode(decoder)?)
    }
}

//####################################################################################
//                          Serialize  VecExtensions                                //
//####################################################################################

// pack the two bit encoded bases into bytes
fn pack_bases(bases: &[u8]) -> Vec<u8> {
    let mut packed = vec![0; bases.len().div_ceil(4)];
    for (i, &base) in bases.iter().enumerate() {
        packed[i / 4] |= (base & 0b11) << (2 * (i % 4));
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

// Group consecutive extensions of the same type
fn group_extensions(extensions: &Vec<Extension>) -> Vec<Vec<Extension>> {
    let mut result = Vec::new();

    let mut current_group = Vec::with_capacity(64);
    let mut current_discriminant = std::mem::discriminant(&current_group[0]);

    for ext in extensions {
        let next_discriminant = std::mem::discriminant(ext);
        if current_discriminant == next_discriminant && current_group.len() < 64 {
            current_group.push(*ext);
        } else {
            result.push(current_group);
            current_group = Vec::with_capacity(64);
            current_group.push(*ext);
            current_discriminant = next_discriminant;
        }
    }

    result.push(current_group);
    result
}

// Encode a vec of extensions composed only of TargetNode
fn encode_target_nodes(
    extensions: &Vec<Extension>,
    encoder: &mut impl Encoder,
) -> Result<(), EncodeError> {
    for ext in extensions {
        if let Extension::TargetNode(node) = ext {
            node.encode(encoder)?;
        } else {
            panic!("Expected only TargetNode in packed patch");
        }
    }
    Ok(())
}

// Encode a vec of extensions composed only of NextNucleotide
fn encode_next_nucleotides(
    extensions: &Vec<Extension>,
    encoder: &mut impl Encoder,
) -> Result<(), EncodeError> {
    let nucl = extensions
        .iter()
        .map(|ext| match ext {
            Extension::NextNucleotide(base) => *base,
            _ => panic!("Expected only NextNucleotide in packed patch"),
        })
        .collect::<Vec<u8>>();
    let packed_nucl = pack_bases(&nucl);
    encoder.writer().write(&packed_nucl)?;
    Ok(())
}

// Encode a vec of extensions composed only of Repetition
fn encode_repetitions(
    extensions: &Vec<Extension>,
    encoder: &mut impl Encoder,
) -> Result<(), EncodeError> {
    for ext in extensions {
        if let Extension::Repetition(nb_repeats, offset) = ext {
            (nb_repeats, offset).encode(encoder)?;
        } else {
            panic!("Expected only Repetition in packed patch");
        }
    }
    Ok(())
}

// Decode a vec of extensions composed only of TargetNode
fn decode_target_nodes<D: Decoder>(
    decoder: &mut D,
    len_patch: usize,
) -> Result<Vec<Extension>, DecodeError> {
    let mut nodes = Vec::with_capacity(len_patch);
    for _ in 0..len_patch {
        let node = Decode::decode(decoder)?;
        nodes.push(Extension::TargetNode(node));
    }
    Ok(nodes)
}

// Decode a vec of extensions composed only of NextNucleotide
fn decode_next_nucleotides<D: Decoder>(
    decoder: &mut D,
    len_patch: usize,
) -> Result<Vec<Extension>, DecodeError> {
    let mut packed_nucleotides: Vec<u8> = Vec::with_capacity(len_patch.div_ceil(4));
    for _ in 0..len_patch.div_ceil(4) {
        packed_nucleotides.push(Decode::decode(decoder)?);
    }
    let mut nucleotides = unpack_bases(&packed_nucleotides);
    nucleotides.truncate(len_patch);
    Ok(nucleotides
        .into_iter()
        .map(Extension::NextNucleotide)
        .collect())
}

// Decode a vec of extensions composed only of Repetition
fn decode_repetitions<D: Decoder>(
    decoder: &mut D,
    len_patch: usize,
) -> Result<Vec<Extension>, DecodeError> {
    let mut repetitions = Vec::with_capacity(len_patch);
    for _ in 0..len_patch {
        let (nb_repeats, offset) = Decode::decode(decoder)?;
        repetitions.push(Extension::Repetition(nb_repeats, offset));
    }
    Ok(repetitions)
}

impl Encode for VecExtensions {
    fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        let patches = group_extensions(&self.0);
        patches.len().encode(encoder)?;
        for patch in patches {
            let id = match patch[0] {
                Extension::TargetNode(_) => 0,
                Extension::NextNucleotide(_) => 1,
                Extension::Repetition(_, _) => 2,
            };
            let header = (patch.len() as u8 - 1) << 2 | id;
            header.encode(encoder)?;

            match id {
                0 => encode_target_nodes(&patch, encoder)?,
                1 => encode_next_nucleotides(&patch, encoder)?,
                2 => encode_repetitions(&patch, encoder)?,
                _ => {
                    panic!("Unreachable");
                }
            }
        }
        Ok(())
    }
}

impl<Context> Decode<Context> for VecExtensions {
    fn decode<D: Decoder>(decoder: &mut D) -> Result<Self, DecodeError> {
        let nb_patches: usize = Decode::decode(decoder)?;
        let mut result = Vec::new();
        for _ in 0..nb_patches {
            let id: u8 = Decode::decode(decoder)?;
            let len_patch = 1 + (id >> 2) as usize;
            let patch = match id & 0b11 {
                0 => decode_target_nodes(decoder, len_patch)?,
                1 => decode_next_nucleotides(decoder, len_patch)?,
                2 => decode_repetitions(decoder, len_patch)?,
                other => Err(DecodeError::UnexpectedVariant {
                    type_name: "ExtensionVec",
                    found: other as u32,
                    allowed: &AllowedEnumVariants::Range { min: 0, max: 2 },
                })?,
            };
            result.extend(patch);
        }
        Ok(VecExtensions(result))
    }
}
impl<'d, Context> BorrowDecode<'d, Context> for VecExtensions {
    fn borrow_decode<D: BorrowDecoder<'d>>(decoder: &mut D) -> Result<Self, DecodeError> {
        Ok(Decode::decode(decoder)?)
    }
}
