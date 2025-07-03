use super::Node;

use vint64;

use bincode::error::{DecodeError, EncodeError};
use bincode::{
    de::{BorrowDecoder, Decoder},
    enc::{write::Writer, Encoder},
};
use bincode::{BorrowDecode, Decode, Encode};

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
