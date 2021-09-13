//! Iterator implementation for emitting all kmers of length k from a sequence
//!
//! Because the iterator is lazy, this implementation should be more
//! memory efficient, but at the cost of a less convenient programming
//! interface. Whereas the function version actually creates vectors
//! and holds them in memory, this iterator version can feed kmers
//! one-by-one.
//!
//! This Iterator is very similar to the built-in std::slice::Windows

struct KMers<'a> {
    sequence: &'a str,
    k: usize,
}

impl<'a> Iterator for KMers<'a> {
    type Item = &'a str;
    fn next(&mut self) -> Option<&'a str> {
        if self.k <= 0 {
            return None;
        }
        if self.sequence.len() < self.k {
            return None;
        }
        let kmer = &self.sequence[0..self.k];
        self.sequence = &self.sequence.get(1..).unwrap_or("");
        return Some(kmer);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmers_iterator_basic() {
        let kmers = KMers {
            sequence: "ABCD",
            k: 2,
        };
        assert_eq!(kmers.collect::<Vec<_>>(), ["AB", "BC", "CD"]);
    }
}
