extern crate users;
extern crate clap;
extern crate time;
extern crate reqwest;
//use self::time::*;
//macro_use must be ahead of "extern crate" in order to import macros from accurity
#[macro_use]
extern crate accurity;
//use accurity; //useless
use clap::{Arg, App, SubCommand};
use std::collections::HashMap;
use std::io::prelude::*;
use std::process;
use std::process::Command;
use users::{get_user_by_uid, get_current_uid};
use reqwest::StatusCode;

fn main() {
    let matches = App::new("Accurity")
        .version("32acfd1e-debug")
        .author("www.yfish.org")
        .about("A program that infers tumor purity, ploidy from tumor-normal WGS data")
        .subcommand(SubCommand::with_name("gc_index")
            .about("GC indexing for a reference genome, counting the number of GCs in \
                    overlapping windows of 1/5/25/125bp")
            .version("32acfd1e-debug")
            .author("www.yfish.org")
            .arg(Arg::with_name("input_file")
                .short("i")
                .long("input_file")
                .value_name("INPUT FILE")
                .help("The input reference genome, fasta format")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("output_dir")
                .short("o")
                .long("outputdir")
                .value_name("OUTPUT DIR")
                .help("The output folder that will contain all gc_index files")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("debug")
                .short("d")
                .help("print debug information verbosely")
            )
        )
        .subcommand(SubCommand::with_name("normalize")
            .about("Smooth (window on either side also gets coverage) and normalize (divided by total fragment count) coverage of tumor and normal. output coverage ratio \
                    (tumor/normal).")
            .version("32acfd1e-debug")
            .author("www.yfish.org")
            .arg(Arg::with_name("read_len")
                .short("l")
                .long("read_len")
                .value_name("READ LENGTH")
                .help("Length of reads in sequencing")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("max_coverage")
                .short("x")
                .long("max_coverage")
                .value_name("MAXIMUM COVERAGE")
                .help("Coverage above this value is ignored.")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("smooth_window_half_size")
                .short("s")
                .long("smooth_window_half_size")
                .value_name("SMOOTH WINDOW HALF SIZE")
                .help("Number of windows on either side that will be used to smooth coverage and tumor/normal coverage ratio.")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("window_size")
                .short("w")
                .long("window_size")
                .value_name("WINDOW SIZE")
                .help("All fragments are grouped into windows to accrue coverage. Usually, it is equal to or smaller than the average insert size in paired end sequencing.")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("tumor_file_path")
                .short("t")
                .long("tumor_file_path")
                .value_name("TUMOR BAM FILE")
                .help("The tumor bam file.")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("normal_file_path")
                .short("n")
                .long("normal_file_path")
                .value_name("NORMAL BAM FILE")
                .help("The normal bam file")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("output_folder")
                .short("o")
                .long("output_folder")
                .value_name("OUTPUT FOLDER")
                .help("The output folder to contain regression data")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("debug")
                .short("d")
                .long("debug")
                .help("Debug mode. NOT used. no difference.")
                .takes_value(true)
            )
        )
        .subcommand(SubCommand::with_name("select_het_snp")
            .about("Select heterozygous SNPs")
            .version("32acfd1e-debug")
            .author("www.yfish.org")
            .arg(Arg::with_name("abp_max_tumor")
                .long("abp_max_tumor")
                .value_name("MAXIMUM ABP FOR TUMOR")
                .help("maximum ABP allelic balance phred allowed for the normal sample")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("abp_max_normal")
                .long("abp_max_normal")
                .value_name("MAXIMUM ABP FOR NORMAL")
                .help("maximum ABP allelic balance phred allowed for the normal sample")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("srp_max")
                .short("r")
                .long("srp_max")
                .value_name("MAXIMUM SRP")
                .help("maximum strand reference allele balance phred allowed")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("sap_max")
                .short("a")
                .long("sap_max")
                .value_name("MAXIMUM SAP")
                .help("maximum strand alternative allele balance phred allowed")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("max_coverage")
                .short("x")
                .long("max_coverage")
                .value_name("MAXIMUM COVERAGE")
                .help("Coverage above this value is ignored.")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("min_coverage")
                .short("m")
                .long("min_coverage")
                .value_name("MINIMUM COVERAGE")
                .help("Coverage below this value is ignored.")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("snp_file_path_tumor")
                .short("t")
                .long("snp_file_path_tumor")
                .value_name("TUMOR BAM FILE")
                .help("The tumor bam file.")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("snp_file_path_normal")
                .short("n")
                .long("snp_file_path_normal")
                .value_name("NORMAL BAM FILE")
                .help("The normal bam file")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("output_file_path")
                .short("o")
                .long("output_file_path")
                .value_name("OUTPUT FILE")
                .help("The output file to contain selected heterozygous SNPs")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("debug")
                .short("d")
                .long("debug")
                .help("Debug mode. NOT used. no difference.")
                .takes_value(true)
            )
        )
        .subcommand(SubCommand::with_name("infer")
            .about("infers tumor purity, ploidy from tumor-normal WGS data")
            .version("32acfd1e-debug")
            .author("www.yfish.org")
            .arg(Arg::with_name("config")
                .short("c")
                .long("config")
                .value_name("FILE")
                .help("Sets a custom config file")
                .required(true)
                .takes_value(true)
            )
        )
        .subcommand(SubCommand::with_name("recall_precision")
            .about("calculate recall and precision from truth result and predicted result")
            .version("32acfd1e-debug")
            .author("www.yfish.org")
            .arg(Arg::with_name("truth_result_file_path")
                .short("t")
                .long("truth_result_file_path")
                .value_name("FILE")
                .help("The file contains the truth result")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("predicted_result_file_path")
                .short("p")
                .long("predicted_result_file_path")
                .value_name("FILE")
                .help("The file contains the predicted result")
                .required(true)
                .takes_value(true)
            )
            .arg(Arg::with_name("output_file_path")
                .short("o")
                .long("output_file_path")
                .value_name("FILE")
                .help("The output file to contains recall and precision")
                .required(true)
                .takes_value(true)
            )
        )
        .get_matches();

    //set the default compile_time to future to guarantee that it will run.
    let compile_time = 1531365721;
    let run_time = time::now_utc().to_timespec().sec;
    //1296000 is 15-day in seconds.
    if run_time-compile_time>=1296000*12 {
        eprintln!("ERROR: your Accurity is out of date. Please go to https://www.yfish.org/display/PUB/Accurity to register and download one.");
        process::exit(1);
    }

    if let Some(matches) = matches.subcommand_matches("gc_index") {
        let input_filename = matches.value_of("input_file").unwrap();
        let output_dir = matches.value_of("output_dir").unwrap();
        let arguments = format!("-i {} -o {}", input_filename, output_dir);
        accurity::gc_index(input_filename, output_dir);
    } else if let Some(matches) = matches.subcommand_matches("normalize") {
        let tumor_file_path = matches.value_of("tumor_file_path").unwrap();
        let normal_file_path = matches.value_of("normal_file_path").unwrap();
        let output_folder = matches.value_of("output_folder").unwrap();
        let max_coverage: usize = matches.value_of("max_coverage").unwrap().parse().unwrap();
        let smooth_window_half_size: usize = matches.value_of("smooth_window_half_size").unwrap().parse().unwrap();
        let window_size: usize = matches.value_of("window_size").unwrap().parse().unwrap();
        let read_len: usize = matches.value_of("read_len").unwrap().parse().unwrap();
        let debug: i32 = matches.value_of("debug").unwrap().parse().unwrap();

        let arguments = format!("-t {} -n {} -w {} -l {} --smooth_window_half_size {} --max_coverage {} -d {} -o {}",
                                tumor_file_path, normal_file_path, window_size, read_len,
                                smooth_window_half_size, max_coverage, debug, output_folder);
        let ins = accurity::normalize::Normalize::new(tumor_file_path, normal_file_path, output_folder,
                                 window_size, read_len,
                                 max_coverage, smooth_window_half_size,
                                 debug);
        ins.run();
    } else if let Some(matches) = matches.subcommand_matches("select_het_snp") {
        let snp_file_path_tumor = matches.value_of("snp_file_path_tumor").unwrap();
        let snp_file_path_normal = matches.value_of("snp_file_path_normal").unwrap();
        let output_file_path = matches.value_of("output_file_path").unwrap();
        let min_coverage: usize = matches.value_of("min_coverage").unwrap().parse().unwrap();
        let max_coverage: usize = matches.value_of("max_coverage").unwrap().parse().unwrap();
        let abp_max_tumor: f32 = matches.value_of("abp_max_tumor").unwrap().parse().unwrap();
        let abp_max_normal: f32 = matches.value_of("abp_max_normal").unwrap().parse().unwrap();
        let srp_max: f32 = matches.value_of("srp_max").unwrap().parse().unwrap();
        let sap_max: f32 = matches.value_of("sap_max").unwrap().parse().unwrap();
        let debug: i32 = matches.value_of("debug").unwrap().parse().unwrap();

        let arguments = format!("-t {} -n {} --sap_max {} --srp_max {} --min_coverage {} --max_coverage {} --abp_max_normal {} --abp_max_tumor {} -d {} -o {}",
                                snp_file_path_tumor, snp_file_path_normal, sap_max, srp_max, min_coverage, max_coverage,
                                abp_max_normal, abp_max_tumor, debug, output_file_path);

        let ins = accurity::select_het_snp::SelectHetSNP::new(snp_file_path_tumor, snp_file_path_normal, output_file_path,
                                                              abp_max_tumor, abp_max_normal, srp_max, sap_max,
                                                              min_coverage, max_coverage, debug);
        ins.run();
    }else if let Some(matches) = matches.subcommand_matches("recall_precision") {
        let truth_result_file_path = matches.value_of("truth_result_file_path").unwrap();
        let predicted_result_file_path = matches.value_of("predicted_result_file_path").unwrap();
        let output_file_path = matches.value_of("output_file_path").unwrap();

        let arguments= format!("--truth_result_file_path {} --predicted_result_file_path {} -o {}",
                               truth_result_file_path, predicted_result_file_path, output_file_path);
        let ins = accurity::recall_precision::RecallPrecision::new(truth_result_file_path,
                                                                   predicted_result_file_path,output_file_path);
        ins.run();
    } else {
        eprintln!("Sub command not present. Try --help.");
    }

}
