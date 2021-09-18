//! Calculate kmer frequency

use log::info;

use anyhow::Result;

use std::fs;
use std::path::PathBuf;

use clap_verbosity_flag::Verbosity;
use structopt::StructOpt;

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
    #[structopt(short, long, default_value = "fasta")]
    extensions: Vec<String>,

    /// input directory
    #[structopt(parse(from_os_str), default_value = ".")]
    directory: PathBuf,

    /// output directory root
    #[structopt(parse(from_os_str), default_value = "./output")]
    output_root: PathBuf,

    /// verbosity
    #[structopt(flatten)]
    verbose: Verbosity,
}

fn main() -> Result<()> {
    let opt: Opt = Opt::from_args();
    opt.verbose.log_level().map(loggerv::init_with_level);

    let input_root = opt.directory.canonicalize()?;

    let fasta_paths = kmer::fs_find_files_with_extensions(input_root.as_path(), &opt.extensions)?;
    for fasta_path in fasta_paths {
        let output_path = kmer::output_path_from_input(&fasta_path, &input_root, &opt.output_root)?;
        fs::create_dir_all(output_path.parent().expect("Invalid paths"))
            .expect("Could not create directory");

        info!(
            "Counting kmers in {:?}. Output to {:?}",
            fasta_path, output_path
        );
        kmer::run_fasta_kmer_count(&fasta_path, opt.k, &output_path)?
    }

    Ok(())
}
