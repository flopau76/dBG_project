use std::fmt::Debug;
use std::hash::Hash;

use packed_seq::{unpack_base, PackedSeq, Seq};

use crate::Side;

pub trait IntHelp {
    /// Reverse the order of 2-bit units of the integer
    fn reverse_by_twos(&self) -> Self;
    fn higher_of_two() -> Self;
}

impl IntHelp for u128 {
    #[inline]
    fn reverse_by_twos(&self) -> u128 {
        // swap adjacent pairs
        let mut r = ((self & 0x33333333333333333333333333333333u128) << 2)
            | ((self >> 2) & 0x33333333333333333333333333333333u128);
        // swap nibbles
        r = ((r & 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0Fu128) << 4)
            | ((r >> 4) & 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0Fu128);
        // swap bytes
        r = ((r & 0x00FF00FF00FF00FF00FF00FF00FF00FFu128) << 8)
            | ((r >> 8) & 0x00FF00FF00FF00FF00FF00FF00FF00FFu128);
        // swap 2 bytes
        r = ((r & 0x0000FFFF0000FFFF0000FFFF0000FFFFu128) << 16)
            | ((r >> 16) & 0x0000FFFF0000FFFF0000FFFF0000FFFFu128);
        // swap 4 bytes
        r = ((r & 0x00000000FFFFFFFF00000000FFFFFFFFu128) << 32)
            | ((r >> 32) & 0x00000000FFFFFFFF00000000FFFFFFFFu128);
        // swap 8 bytes
        r = ((r & 0x0000000000000000FFFFFFFFFFFFFFFFu128) << 64)
            | ((r >> 64) & 0x0000000000000000FFFFFFFFFFFFFFFFu128);
        r
    }
    #[inline]
    fn higher_of_two() -> u128 {
        0xAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAu128
    }
}

impl IntHelp for u64 {
    #[inline]
    fn reverse_by_twos(&self) -> u64 {
        // swap adjacent pairs
        let mut r = ((self & 0x3333333333333333u64) << 2) | ((self >> 2) & 0x3333333333333333u64);
        // swap nibbles
        r = ((r & 0x0F0F0F0F0F0F0F0Fu64) << 4) | ((r >> 4) & 0x0F0F0F0F0F0F0F0Fu64);
        // swap bytes
        r = ((r & 0x00FF00FF00FF00FFu64) << 8) | ((r >> 8) & 0x00FF00FF00FF00FFu64);
        // swap 2 bytes
        r = ((r & 0x0000FFFF0000FFFFu64) << 16) | ((r >> 16) & 0x0000FFFF0000FFFFu64);
        // swap 4 bytes
        r = ((r & 0x00000000FFFFFFFFu64) << 32) | ((r >> 32) & 0x00000000FFFFFFFFu64);
        r
    }
    #[inline]
    fn higher_of_two() -> u64 {
        0xAAAAAAAAAAAAAAAAu64
    }
}

impl IntHelp for u32 {
    #[inline]
    fn reverse_by_twos(&self) -> u32 {
        // swap adjacent pairs
        let mut r = ((self & 0x33333333u32) << 2) | ((self >> 2) & 0x33333333u32);
        // swap nibbles
        r = ((r & 0x0F0F0F0Fu32) << 4) | ((r >> 4) & 0x0F0F0F0Fu32);
        // swap bytes
        r = ((r & 0x00FF00FFu32) << 8) | ((r >> 8) & 0x00FF00FFu32);
        // swap 2 bytes
        r = ((r & 0x0000FFFFu32) << 16) | ((r >> 16) & 0x0000FFFFu32);
        r
    }
    #[inline]
    fn higher_of_two() -> u32 {
        0xAAAAAAAAu32
    }
}

impl IntHelp for u16 {
    #[inline]
    fn reverse_by_twos(&self) -> u16 {
        // swap adjacent pairs
        let mut r = ((self & 0x3333u16) << 2) | ((self >> 2) & 0x3333u16);
        // swap nibbles
        r = ((r & 0x0F0Fu16) << 4) | ((r >> 4) & 0x0F0Fu16);
        // swap bytes
        r = ((r & 0x00FFu16) << 8) | ((r >> 8) & 0x00FFu16);
        r
    }
    #[inline]
    fn higher_of_two() -> u16 {
        0xAAAAu16
    }
}

impl IntHelp for u8 {
    #[inline]
    fn reverse_by_twos(&self) -> u8 {
        // swap adjacent pairs
        let mut r = ((self & 0x33u8) << 2) | ((self >> 2) & 0x33u8);
        // swap nibbles
        r = ((r & 0x0Fu8) << 4) | ((r >> 4) & 0x0Fu8);
        r
    }
    #[inline]
    fn higher_of_two() -> u8 {
        0xAAu8
    }
}

pub trait KmerStorage: Sized + IntHelp + Copy + Clone + Hash + Debug + Eq + Send + Sync {
    /// Maximal number of nucleotides that can be stored in this type.
    fn capacity() -> usize;
    /// Create a new empty k-mer storage.
    fn new() -> Self;
    /// Extend the k-mer to the left by adding a base, given its 2-bit encoding.
    fn extend_left(&mut self, k: usize, base: u8);
    /// Print the k-mer as a string of bases.
    fn print(&self, k: usize) -> String;
    /// Extend the k-mer to the right by adding a base, given its 2-bit encoding.
    fn extend_right(&mut self, k: usize, base: u8);
    /// Extend the k-mer to the given `side`.
    fn extend_side(&mut self, k: usize, base: u8, side: Side) {
        match side {
            Side::Left => self.extend_left(k, base),
            Side::Right => self.extend_right(k, base),
        }
    }
    /// Get the kmer starting at the given position in the sequence.
    fn get_kmer(k: usize, seq: PackedSeq, pos: usize) -> Self {
        let mut kmer = Self::new();
        for i in 0..k {
            let base = seq.get(pos + i);
            kmer.extend_right(k, base);
        }
        kmer
    }
    fn rc(self, k: usize) -> Self;
}

macro_rules! impl_kmer_storage {
    ($type:ty) => {
        impl KmerStorage for $type {
            fn capacity() -> usize {
                4 * std::mem::size_of::<Self>()
            }
            fn new() -> Self {
                0 as Self
            }
            fn extend_right(&mut self, k: usize, base: u8) {
                let bits = base as Self;
                *self = (*self << 2) | bits;
                if k < 4 * std::mem::size_of::<Self>() {
                    let mask: Self = (1 << (2 * k)) - 1;
                    *self &= mask;
                }
            }
            fn extend_left(&mut self, k: usize, base: u8) {
                let bits = base as Self;
                *self = (*self >> 2) | (bits << (2 * (k - 1)));
            }
            fn rc(self, k: usize) -> Self {
                let mut rc = self.reverse_by_twos() ^ Self::higher_of_two();
                let offset = 2 * (std::mem::size_of::<Self>() * 4 - k);
                if offset > 0 {
                    rc = rc >> offset;
                }
                rc
            }
            fn print(&self, k: usize) -> String {
                let mut bits = *self;
                let mut pairs = Vec::new();
                while bits > 0 {
                    pairs.push((bits & 0b11) as u8);
                    bits >>= 2;
                }
                while pairs.len() < k {
                    pairs.push(0);
                }
                pairs.reverse();
                let ascii = pairs.iter().map(|&x| unpack_base(x)).collect::<Vec<u8>>();
                unsafe { String::from_utf8_unchecked(ascii.clone()) }
            }
        }
    };
}

impl_kmer_storage!(u8);
impl_kmer_storage!(u16);
impl_kmer_storage!(u32);
impl_kmer_storage!(u64);
impl_kmer_storage!(u128);

#[cfg(test)]
mod tests {
    use super::*;
    use packed_seq::complement_base;
    #[test]
    fn extend_right() {
        let k = 3;
        let mut kmer: u8 = KmerStorage::new();
        kmer.extend_right(k, 3);
        kmer.extend_right(k, 1);
        kmer.extend_right(k, 2);
        kmer.extend_right(k, 3);
        assert_eq!(kmer, 0b00011011);
    }
    #[test]
    fn extend_left() {
        let k = 3; // k-mer size
        let mut kmer: u8 = KmerStorage::new();
        kmer.extend_left(k, 3);
        kmer.extend_left(k, 1);
        kmer.extend_left(k, 2);
        kmer.extend_left(k, 3);
        assert_eq!(kmer, 0b00111001);
    }
    #[test]
    fn reverse_complement() {
        let k = 3; // k-mer size
        let mut kmer: u8 = KmerStorage::new();
        kmer.extend_right(k, 1);
        kmer.extend_right(k, 2);
        kmer.extend_right(k, 3);
        let mut kmer_rc: u8 = KmerStorage::new();
        kmer_rc.extend_left(k, complement_base(1));
        kmer_rc.extend_left(k, complement_base(2));
        kmer_rc.extend_left(k, complement_base(3));
        assert_eq!(kmer.rc(k), kmer_rc);
    }
}
