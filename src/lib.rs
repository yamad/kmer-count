use std::collections::HashMap;
use std::io::Write;
use std::fs;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::str;

use bio::io::fasta;

use anyhow::Result;
use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
enum KmerError {
    #[error("No valid kmers. kmer length is {k:?}, but must be 1 or greater")]
    KmerLengthTooSmall { k: usize },

    #[error(
        "No valid kmers. Sequence length {seq_len:?} is smaller than requested kmer length {k:?}"
    )]
    KmerLengthTooLong { k: usize, seq_len: usize },

    #[error("Suspect base(s) found: {bases:?}. Use only ATCG bases")]
    IncorrectBases { bases: String },
}

#[derive(Eq, PartialEq, Debug)]
struct KmerRecord<'b> {
    seq: &'b str,
    count: u64,
}

/// Aggregate count of all Kmers
type KmerCount<'a> = Vec<KmerRecord<'a>>;

/// Save counts for length `k` kmers from the fasta file at `fasta_path` at `output_path`
pub fn run_fasta_kmer_count(fasta_path: &PathBuf, k: usize, output_path: &PathBuf) -> Result<()> {
    let fasta_file = File::open(fasta_path)?;
    let reader = fasta::Reader::new(fasta_file);

    for record in reader.records() {
        let record = record?;

        if let Err(err) = check_bases(record.seq()) {
            println!("WARNING: {}", err);
        }

        match count_kmers(record.seq(), k) {
            Ok(kmer_count) => save_kmer_count(kmer_count, output_path)?,
            Err(err) => eprintln!("ERROR: {}", err),
        }
    }
    Ok(())
}

/// Return frequency of all kmers of length `k` in `sequence`, ordered from most to least abundant
fn count_kmers(sequence: &[u8], k: usize) -> Result<KmerCount, KmerError> {
    // calculate kmer frequencies
    let mut counter: HashMap<&str, u64> = HashMap::new();
    for kmer in kmers(sequence, k)? {
        *counter.entry(kmer).or_insert(0) += 1;
    }

    // order from most to least abundant
    let mut ordered: Vec<_> = counter
        .into_iter()
        .map(|(k, v)| KmerRecord { seq: k, count: v })
        .collect();

    // first, sort kmers alphabetically so order among equal counts is deterministic
    // second, sort by descending count
    //
    // n.b. this could be implemented by the Ord/PartialOrd traits on KmerRecord,
    // but for this simple program, putting the sorting logic here is clearer and
    // results in less boilerplate.
    ordered.sort_by(|a, b| a.seq.cmp(b.seq));
    ordered.sort_by(|a, b| a.count.cmp(&b.count).reverse());
    Ok(ordered)
}

/// Return all subsequences of length k from the given sequence
///
/// `sequence` must be an ASCII string, which is sufficient for sequencing data.
/// Multi-byte UTF-8 characters are not handled correctly.
fn kmers(sequence: &[u8], k: usize) -> Result<impl Iterator<Item = &str>, KmerError> {
    if k <= 0 {
        return Err(KmerError::KmerLengthTooSmall { k: k });
    }

    if sequence.len() < k {
        return Err(KmerError::KmerLengthTooLong {
            k: k,
            seq_len: sequence.len(),
        });
    }

    Ok(sequence.windows(k).flat_map(|x| str::from_utf8(x))) // from string, so utf-8 cast will always succeed
}

/// Check that all bases in `seq` are A, T, C, or G.
fn check_bases(seq: &[u8]) -> Result<(), KmerError> {
    let mut bad_bases = Vec::new();
    for base in seq {
        if !b"ATCG".contains(base) {
            bad_bases.push(*base);
        }
    }

    if bad_bases.is_empty() {
        Ok(())
    } else {
        let bases: String = String::from_utf8(bad_bases).unwrap();
        Err(KmerError::IncorrectBases { bases: bases })
    }
}

/// Derive an output file path from the suffix of the input path
pub fn output_path_from_input(
    input_path: &PathBuf,
    input_root: &PathBuf,
    output_root: &PathBuf,
) -> Result<PathBuf> {
    let path_stub = input_path.strip_prefix(&input_root)?;
    let mut output_path = output_root.join(path_stub);
    output_path.set_file_name(format!(
        "{}_kmer.txt",
        input_path.file_stem().unwrap().to_str().unwrap()
    ));
    Ok(output_path)
}

/// Find all files in directory `dir` with one of the given `extensions`
pub fn fs_find_files_with_extensions<T>(dir: &Path, extensions: &[T]) -> Result<Vec<PathBuf>>
where
    T: AsRef<str>,
{
    fn is_file_type<T: AsRef<str>>(p: &PathBuf, exts: &[T]) -> bool {
        p.is_file()
            && p.extension()
                .map(|s| exts.iter().any(|e| s == e.as_ref()))
                .unwrap_or(false)
    }

    let mut files = Vec::new();
    for entry in dir.read_dir()? {
        let entry = entry?;
        let path = fs::canonicalize(entry.path())?;
        if is_file_type(&path, &extensions) {
            files.push(path);
        }
    }
    Ok(files)
}

/// Save kmer count to `output_path`
fn save_kmer_count(kmer_count: KmerCount, output_path: &PathBuf) -> Result<()> {
    let mut file = File::create(output_path)?;

    writeln!(file, "kmer\tcount")?;
    for kmer in kmer_count {
        writeln!(file, "{}\t{}", kmer.seq, kmer.count)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    // test helper to run the kmers iterator
    fn kmers_vec(sequence: &[u8], k: usize) -> Vec<&str> {
        kmers(sequence, k).unwrap().collect()
    }

    #[test]
    /// kmers_vec simplifies comparisons
    fn test_kmers_vec_helper_demo() -> Result<(), KmerError> {
        assert_eq!(
            kmers(b"ABCD", 2)?.collect::<Vec<_>>(),
            kmers_vec(b"ABCD", 2)
        );
        Ok(())
    }

    #[test]
    fn test_kmers_fn_basic() {
        assert_eq!(kmers_vec(b"ABCD", 1), ["A", "B", "C", "D"]);
        assert_eq!(kmers_vec(b"ABCD", 2), ["AB", "BC", "CD"]);
        assert_eq!(kmers_vec(b"ABCD", 3), ["ABC", "BCD"]);
        assert_eq!(kmers_vec(b"ABCD", 4), ["ABCD"]);
    }

    #[test]
    fn test_kmer_0() -> Result<(), String>{
        match kmers(b"ABCD", 0) {
            Err(e) => Ok(assert_eq!(e, KmerError::KmerLengthTooSmall { k: 0 })),
            Ok(_) => Err(String::from("Should have generated error on k = 0")),
        }
    }

    #[test]
    fn test_kmer_empty_string() -> Result<(), String> {
        match kmers(b"", 10) {
            Err(err) => Ok(assert_eq!(err, KmerError::KmerLengthTooLong { k: 10, seq_len: 0 })),
            Ok(_) => Err(String::from("Should have generated error on empty string"))
        }
    }

    #[test]
    fn test_kmer_k_too_big() -> Result<(), String> {
        match kmers(b"ABC", 10) {
            Err(err) => Ok(assert_eq!(err, KmerError::KmerLengthTooLong { k: 10, seq_len: 3 })),
            Ok(_) => Err(String::from("Should have generated error when k > length of sequence"))
        }
    }

    /// test helper to convert tuple vector to KmerCount
    fn kmer_count_from_tuples<'a>(item: Vec<(&'a str, u64)>) -> KmerCount<'a> {
        item.into_iter()
            .map(|x| KmerRecord {
                seq: x.0,
                count: x.1,
            })
            .collect()
    }

    #[test]
    fn test_count_kmers() {
        let sequence = b"ATCGGATCG";
        let expected: KmerCount = kmer_count_from_tuples(vec![
            ("ATC", 2),
            ("TCG", 2),
            ("CGG", 1),
            ("GAT", 1),
            ("GGA", 1),
        ]);
        assert_eq!(count_kmers(sequence, 3).unwrap(), expected);
    }

    #[test]
    fn test_output_path_from_input() {
        let input_path = PathBuf::from("/a/input/dir/path.txt");
        let input_root = PathBuf::from("/a/input");
        let output_root = PathBuf::from("/output");
        assert_eq!(
            "/output/dir/path_kmer.txt",
            output_path_from_input(&input_path, &input_root, &output_root)
                .unwrap()
                .to_str()
                .unwrap()
        );
    }

    #[test]
    fn test_check_bases_success() {
        check_bases(b"ATCGATGCAAA").unwrap();
    }

    #[test]
    fn test_check_bases_bad_base() {
        assert_eq!(check_bases(b"ATCNTTZ").unwrap_err(),
        KmerError::IncorrectBases { bases: String::from("NZ") });
    }

    #[test]
    fn test_find_files() -> Result<()>{
        let dir = tempdir()?;

        let found_file_path = dir.path().join("foo.txt");
        File::create(&found_file_path)?;

        let missing_file_path = dir.path().join("bar.baz");
        File::create(&missing_file_path)?;

        let files = fs_find_files_with_extensions(dir.path(), &vec!["rs", "txt"])?;

        println!("{:?}", files);
        println!("{:?}", found_file_path);

        assert!(files.contains(&found_file_path.canonicalize()?));
        assert!(!files.contains(&missing_file_path));

        Ok(())
    }

    #[test]
    #[should_panic(expected = "Not a directory")]
    fn test_find_files_dir_is_file() {
        fs_find_files_with_extensions(&Path::new("./output.txt"), &vec!["rs", "txt"]).unwrap();
    }
}
