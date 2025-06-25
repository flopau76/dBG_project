pub trait KmerStorage {
    /// Create a new empty k-mer storage.
    fn new() -> Self;
    /// Extend the k-mer to the left by adding a base, given its 2-bit encoding.
    fn extend_left(&mut self, k: usize, base: u8);
    /// Extend the k-mer to the right by adding a base, given its 2-bit encoding.
    fn extend_right(&mut self, k: usize, base: u8);
}

macro_rules! impl_kmer_storage {
    ($type:ty) => {
        impl KmerStorage for $type {
            fn new() -> Self {
                0 as $type
            }
            fn extend_left(&mut self, k: usize, base: u8) {
                let bits = base as $type;
                *self = (*self << 2) | bits;
                // If we've exceeded the k-mer size, mask off the extra bits
                if k < (std::mem::size_of::<$type>() * 8) {
                    let mask = (1 << (2 * k)) - 1;
                    *self &= mask as $type;
                }
            }

            fn extend_right(&mut self, k: usize, base: u8) {
                let bits = base as $type;
                *self = (*self >> 2) | (bits << (2 * (k - 1)));
            }
        }
    };
}

impl_kmer_storage!(u8);
impl_kmer_storage!(u16);
impl_kmer_storage!(u32);
impl_kmer_storage!(u64);
impl_kmer_storage!(u128);
