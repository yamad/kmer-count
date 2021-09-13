//! Read FASTA format files
//!
//! Implemented "by-hand" for the problem. In real life, prefer `bio::io::fasta` in the `bio` crate

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

use bio::io::fasta;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<impl BufRead>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;

    #[test]
    fn test_emit_fasta_record() -> Result<(), Box<dyn Error>> {
        let cursor = io::Cursor::new(
            b"\
        > First record\n\
        ATCG\n\
        ATCG\n\
        > Second record\n\
        GCTA\n\
        GCTA\n\
        ",
        );

        let expected_records = [
            ("First record", b"ATCGATCG"),
            ("Second record", b"GCTAGCTA"),
        ];

        let reader = fasta::Reader::new(cursor);
        let tests = expected_records.iter().zip(reader.records());
        for (expected, actual) in tests {
            let record = actual.unwrap();
            assert_eq!(record.desc().unwrap(), expected.0);
            assert_eq!(record.seq().to_vec(), expected.1);
        }

        Ok(())
    }
}
