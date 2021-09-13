//! Calculate kmer frequency

#[macro_use]
extern crate log;

mod fasta_parser;

use std::collections::HashMap;
use std::fs;
use std::io;
use std::str;

use crate::KMerError::{KmerLengthTooLong, KmerLengthTooSmall};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use bio::io::fasta;
use structopt::StructOpt;

#[derive(Debug)]
enum KMerError {
    KmerLengthTooSmall,
    KmerLengthTooLong,
}

#[derive(Eq, PartialEq, Debug)]
struct KMerRecord<'b> {
    seq: &'b str,
    count: u64,
}

/// Aggregate count of all kmers
type KMerCount<'a> = Vec<KMerRecord<'a>>;

#[derive(Debug, StructOpt)]
#[structopt(
    name = "kmer count",
    about = "Count frequency of all kmers for all fasta files in directory"
)]
struct Opt {
    /// length of kmer
    #[structopt(short)]
    k: usize,

    /// input file extensions to find
    #[structopt(short, long, default_value = "fa, fasta")]
    extensions: Vec<String>,

    /// input directory
    #[structopt(parse(from_os_str), default_value = ".")]
    directory: PathBuf,

    /// output directory root
    #[structopt(parse(from_os_str), default_value = "./output")]
    output_root: PathBuf,
}

fn main() -> io::Result<()> {
    env_logger::init();

    let opt = Opt::from_args();

    let input_root = fs::canonicalize(opt.directory)?;

    let fasta_paths = fs_find_files_with_extensions(input_root.as_path(), &opt.extensions)?;

    println!("{:?}", input_root);
    println!("{:?}", fasta_paths);

    for fasta_path in fasta_paths {
        let output_path = output_path_from_input(&fasta_path, &input_root, &opt.output_root)?;
        info!("Counting in {:?}. Output to {:?}", fasta_path, output_path);
        run_fasta_kmer_count(&fasta_path, opt.k, &output_path)?
    }

    Ok(())
}

/// Derive output file path from input path
///
/// If the
fn output_path_from_input(
    input_path: &PathBuf,
    input_root: &PathBuf,
    output_root: &PathBuf,
) -> io::Result<PathBuf> {
    let path_stub = input_path.strip_prefix(&input_root).unwrap();
    let mut output_path = output_root.join(path_stub);
    output_path.set_file_name(format!(
        "{}_kmer.txt",
        input_path.file_stem().unwrap().to_str().unwrap()
    ));
    Ok(output_path)
}

fn run_fasta_kmer_count(fasta_path: &PathBuf, k: usize, output_path: &PathBuf) -> io::Result<()> {
    let fasta_file = File::open(fasta_path)?;
    let reader = fasta::Reader::new(fasta_file);
    for record in reader.records() {
        let record = record?;
        let kmer_count = count_kmers(record.seq(), k);
        save_kmer_count(kmer_count, output_path)?;
    }

    Ok(())
}

/// Find all files in directory `dir` with one of the given `extensions`
///
/// Examples:
///
///     fs_find_files_with_extension("./", ["rs"]);
fn fs_find_files_with_extensions<T>(dir: &Path, extensions: &[T]) -> io::Result<Vec<PathBuf>>
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
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = fs::canonicalize(entry.path())?;
        if is_file_type(&path, &extensions) {
            files.push(path);
        }
    }
    Ok(files)
}

/// Return all subsequences of length k from the given sequence
///
/// `sequence` must be an ASCII string, which is sufficient for sequencing data.
/// Multi-byte UTF-8 characters are not handled correctly.
fn kmers(sequence: &[u8], k: usize) -> Result<impl Iterator<Item = &str>, KMerError> {
    if k <= 0 {
        warn!("No valid kmers. kmer length must be 1 or greater");
        return Err(KmerLengthTooSmall);
    }

    if sequence.len() < k {
        warn!(
            "No valid kmers. Sequence length {} is smaller than requested kmer length {}",
            sequence.len(),
            k
        );
        return Err(KmerLengthTooLong);
    }

    Ok(sequence.windows(k).flat_map(|x| str::from_utf8(x))) // from string, so utf-8 cast will always succeed
}

/// Return frequency of all kmers of length `k` in `sequence`, ordered from most to least abundant
fn count_kmers(sequence: &[u8], k: usize) -> KMerCount {
    // calculate kmer frequencies
    let mut counter: HashMap<&str, u64> = HashMap::new();
    for kmer in kmers(sequence, k).unwrap() {
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
    // but for this simple program, putting the sorting logic here is clearer and
    // results in less boilerplate.
    ordered.sort_by(|a, b| a.seq.cmp(b.seq));
    ordered.sort_by(|a, b| a.count.cmp(&b.count).reverse());
    ordered
}

fn save_kmer_count(kmer_count: KMerCount, output_path: &PathBuf) -> io::Result<()> {
    //fs::create_dir_all(output_path.parent().unwrap_or(output_path))?;
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

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    // test helper to run the kmers iterator
    fn kmers_vec(sequence: &[u8], k: usize) -> Vec<&str> {
        kmers(sequence, k).unwrap().collect()
    }

    #[test]
    /// kmers_vec simplifies comparisons
    fn test_kmers_vec_helper_demo() -> Result<(), KMerError> {
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
    fn test_kmer_0() {
        log_init();

        let empty_vec: Vec<&str> = vec![];
        assert_eq!(kmers_vec(b"ABCD", 0), empty_vec);
    }

    #[test]
    fn test_kmer_empty_string() {
        let empty_vec: Vec<&str> = vec![];
        //  assert_eq!(kmers_vec("", 0), empty_vec);
        assert_eq!(kmers_vec(b"", 10), empty_vec);
    }

    #[test]
    fn test_kmer_k_too_big() {
        let empty_vec: Vec<&str> = vec![];
        assert_eq!(kmers_vec(b"ABCDE", 10), empty_vec);
    }

    /// test helper to convert tuple vector to KMerCount
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
        let sequence = b"ATCGGATCG";
        let expected: KMerCount = kmer_count_from_tuples(vec![
            ("ATC", 2),
            ("TCG", 2),
            ("CGG", 1),
            ("GAT", 1),
            ("GGA", 1),
        ]);
        assert_eq!(count_kmers(sequence, 3), expected);
    }

    #[test]
    fn test_find_files() {
        assert_eq!(
            fs_find_files_with_extensions(&Path::new("."), vec!["rs", "txt"]).unwrap(),
            vec![PathBuf::from("")]
        );
    }

    #[test]
    #[should_panic(expected = "Not a directory")]
    fn test_find_files_dir_is_file() {
        fs_find_files_with_extensions(&Path::new("./output.txt"), vec!["rs", "txt"]).unwrap();
    }
}
