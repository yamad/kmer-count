mod fasta_parser;

use log::warn;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::str;
use std::io;

#[derive(Eq, PartialEq, Debug)]
struct KMerRecord<'b> {
    seq: &'b str,
    count: u64,
}

impl<'a> Ord for KMerRecord<'a> {
    fn cmp(&self, other: &Self) -> Ordering {
        let ord = self.count.cmp(&other.count).reverse();
        if ord == Ordering::Equal {
            self.seq.cmp(&other.seq)
        } else {
            ord
        }
    }
}

impl<'a> PartialOrd for KMerRecord<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

type KMerCount<'a> = Vec<KMerRecord<'a>>;

fn main() -> io::Result<()> {
    let seq = "ABCDEFGHIJK";
    let kmer_count = count_kmers(seq, 3);
    save_kmer_count(kmer_count)?;
    Ok(())
}

/// Return all subsequences of length k from the given sequence
///
/// `sequence` must be an ASCII string, which is sufficient for sequencing data.
/// Multi-byte UTF-8 characters are not handled correctly.
///
/// For a team more familiar with imperative style, this implementation is probably clearest.
pub fn kmers(sequence: &str, k: usize) -> Vec<&str> {
    let mut res = vec![];
    let mut pos = 0;

    if k == 0 {
        return res;
    }

    if k > sequence.len() {
        warn!("requested kmer size {} is larger than sequence.", k);
    }

    while let Some(kmer) = sequence.get(pos..pos + k) {
        res.push(kmer);
        pos += 1;
    }
    return res;
}

///
/// ASCII strings are assumed. Multi-byte UTF-8 characters will not be handled correctly.
///
/// For a team familiar with functional style, I would probably implement the function this way
pub fn kmers2(sequence: &str, k: usize) -> impl Iterator<Item = &str> {
    sequence
        .as_bytes()
        .windows(k)
        .flat_map(|x| str::from_utf8(x)) // from string, so utf-8 conversion will always succeed
}

fn count_kmers(sequence: &str, k: usize) -> KMerCount {
    // count up kmer frequency
    let mut counter: HashMap<&str, u64> = HashMap::new();
    for kmer in kmers2(sequence, k) {
        *counter.entry(kmer).or_insert(0) += 1;
    }

    // order from most to least abundant
    let mut ordered: Vec<_> = counter
        .into_iter()
        .map(|(k, v)| KMerRecord { seq: k, count: v })
        .collect();

    // first, sort kmers alphabetically so order among equal counts is deterministic
    // second, sort by descending count
    //
    // n.b. this could be implemented by the Ord/PartialOrd traits on KMerRecord,
    // but for this simple program, this is more explicit and has less boilerplate.
    ordered.sort_by(|a, b| a.seq.cmp(b.seq));
    ordered.sort_by(|a, b| a.count.cmp(&b.count).reverse());
    ordered
}

use std::fs::File;
use std::io::Write;

fn save_kmer_count(kmer_count: KMerCount) -> std::io::Result<()> {
    let mut file = File::create("output.txt")?;
    for kmer in kmer_count {
        writeln!(file, "{}\t{}", kmer.seq, kmer.count)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmers_fn_basic() {
        assert_eq!(kmers("ABCD", 1), ["A", "B", "C", "D"]);
        assert_eq!(kmers("ABCD", 2), ["AB", "BC", "CD"]);
        assert_eq!(kmers("ABCD", 3), ["ABC", "BCD"]);
        assert_eq!(kmers("ABCD", 4), ["ABCD"]);
    }

    #[test]
    fn test_kmer_0() {
        let empty_vec: Vec<&str> = vec![];
        assert_eq!(kmers("ABCD", 0), empty_vec);
    }

    #[test]
    fn test_kmer_empty_string() {
        let empty_vec: Vec<&str> = vec![];
        assert_eq!(kmers("", 0), empty_vec);
        assert_eq!(kmers("", 10), empty_vec);
    }

    #[test]
    fn test_kmer_k_too_big() {
        let empty_vec: Vec<&str> = vec![];
        assert_eq!(kmers("ABCDE", 10), empty_vec);
    }

    #[test]
    fn test_kmer_windows() {
        println!("{:?}", kmers2("ABCD", 2).collect::<Vec<_>>());
    }

    impl<'a> From<(&'a str, u64)> for KMerRecord<'a> {
        fn from(item: (&'a str, u64)) -> Self {
            KMerRecord {
                seq: item.0,
                count: item.1,
            }
        }
    }

    fn kmer_count_from_tuples<'a>(item: Vec<(&'a str, u64)>) -> KMerCount<'a> {
        item.into_iter()
            .map(|x| KMerRecord {
                seq: x.0,
                count: x.1,
            })
            .collect()
    }

    #[test]
    fn test_count_kmers() {
        let sequence = "ATCGGATCG";
        let expected: KMerCount = kmer_count_from_tuples(vec![
            ("ATC", 2),
            ("TCG", 2),
            ("CGG", 1),
            ("GAT", 1),
            ("GGA", 1),
        ]);
        assert_eq!(count_kmers(sequence, 3), expected);
    }
}
