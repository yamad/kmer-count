# kmer count

A kmer counter for directories of FASTA files.

Implemented in Rust as an exercise.

## Quickstart

  1. [Install Rust](https://www.rust-lang.org/tools/install)
  2. Clone repository and build
  
     ```
     git clone https://github.com/yamad/kmer-count
     cd kmer-count
     cargo build --release
     ```
  3. Run, e.g.
  
     ```
     target/release/kmer -vv -k 4 fasta-directory output-directory
     ```
      
  4. For help and usage:
  
     ```
     kmer --help
     ```
     
## Testing

Run:

```cargo test```

Current results:

```
running 11 tests
test tests::test_check_bases_success ... ok
test tests::test_check_bases_bad_base ... ok
test tests::test_count_kmers ... ok
test tests::test_find_files_dir_is_file - should panic ... ok
test tests::test_kmer_0 ... ok
test tests::test_kmer_empty_string ... ok
test tests::test_kmer_k_too_big ... ok
test tests::test_kmers_fn_basic ... ok
test tests::test_kmers_vec_helper_demo ... ok
test tests::test_output_path_from_input ... ok
test tests::test_find_files ... ok

test result: ok. 11 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
```


## Usage

```
kmer count 0.1.0
Count frequency of all kmers for all fasta files in directory

USAGE:
    kmer [FLAGS] [OPTIONS] -k <k> [--] [ARGS]

FLAGS:
    -h, --help
            Prints help information

    -q, --quiet
            Pass many times for less log output

    -V, --version
            Prints version information

    -v, --verbose
            Pass many times for more log output

            By default, it'll only report errors. Passing `-v` one time also prints warnings, `-vv` enables info
            logging, `-vvv` debug, and `-vvvv` trace.

OPTIONS:
    -e, --extensions <extensions>...
            input file extensions to find [default: fasta]

    -k <k>
            length of kmer


ARGS:
    <directory>
            input directory [default: .]

    <output-root>
            output directory root [default: ./output]

```