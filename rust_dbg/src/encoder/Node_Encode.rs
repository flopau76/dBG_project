use super::Node;

use bincode::error::{DecodeError, EncodeError};
use bincode::{
    de::{BorrowDecoder, Decoder},
    enc::Encoder,
};
use bincode::{BorrowDecode, Decode, Encode};

impl Encode for Node {
    fn encode<__E: Encoder>(&self, encoder: &mut __E) -> Result<(), EncodeError> {
        let value: u128 = ((self.id as u128) << 1) | (self.is_rc as u128);
        Encode::encode(&value, encoder)?; // efficient if encoder uses variable int encoding
        Ok(())
    }
}

impl<Context> Decode<Context> for Node {
    fn decode<D: Decoder>(decoder: &mut D) -> Result<Self, DecodeError> {
        let value: u128 = Decode::decode(decoder)?;
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
