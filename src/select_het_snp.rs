extern crate time;

use flate2;
use flate2::Compression;
use rust_htslib::bcf;
use std::cmp;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::path::{Path};
use std::str;


struct OneSNP {
    chr: String,
    start: u64,
    genotype: String,
    depth: usize,
    abp_value: f32,
    srp_value: f32,
    sap_value: f32,
    no_of_ref_obs: usize,
    no_of_alt_obs: usize,
}

impl OneSNP {
    fn new(chr: String,
           start: u64,
           genotype: String,
           depth: usize,
           abp_value: f32,
           srp_value: f32,
           sap_value: f32,
           no_of_ref_obs: usize,
           no_of_alt_obs: usize,) -> OneSNP {
        OneSNP {
            chr,
            start,
            genotype,
            depth,
            abp_value,
            srp_value,
            sap_value,
            no_of_ref_obs,
            no_of_alt_obs,
        }
    }
}

struct OneGenomeSNP{
    chr_start2snp: HashMap<(u64, u64), OneSNP>,
    no_of_total_records: usize,
    no_of_good_hets: usize,
}

impl OneGenomeSNP{
    fn new() -> OneGenomeSNP{
        OneGenomeSNP{
            chr_start2snp: HashMap::new(),
            no_of_total_records: 0,
            no_of_good_hets: 0,
        }
    }
}

pub struct SelectHetSNP<'a> {
    snp_file_path_tumor: &'a Path,
    snp_file_path_normal: &'a Path,
    output_file_path: &'a Path,
    abp_max_tumor: f32,
    abp_max_normal: f32,
    srp_max: f32,
    sap_max: f32,
    min_coverage: usize,
    max_coverage: usize,
    debug: i32,
}



impl<'a> SelectHetSNP<'a> {
    pub fn new(snp_file_path_tumor: &'a str,
           snp_file_path_normal: &'a str,
           output_file_path: &'a str,
           abp_max_tumor: f32,
           abp_max_normal: f32,
           srp_max: f32,
           sap_max: f32,
           min_coverage: usize,
           max_coverage: usize,
           debug: i32,
    ) -> SelectHetSNP<'a> {
        SelectHetSNP {
            snp_file_path_tumor: Path::new(snp_file_path_tumor),
            snp_file_path_normal: Path::new(snp_file_path_normal),
            output_file_path: Path::new(output_file_path),
            abp_max_tumor,
            abp_max_normal,
            srp_max,
            sap_max,
            min_coverage,
            max_coverage,
            debug,
        }
    }

    fn read_in_het_snp(&'a self, snp_file_path: &'a Path, abp_max: f32) -> OneGenomeSNP {
        println_stderr!("Reading from {:?} with abp_max={}, srp_max={}, sap_max={}, min_coverage={}, max_coverage={} ...",
            &snp_file_path, abp_max, self.srp_max, self.sap_max, self.min_coverage, self.max_coverage);
        let mut vcf = bcf::Reader::from_path(&snp_file_path).ok().expect("Error opening file.");
        let mut one_genome_snp = OneGenomeSNP::new();
        let vcf_header = vcf.header().clone();
        for (i, rec) in vcf.records().enumerate() {
            one_genome_snp.no_of_total_records += 1;
            let mut record = rec.ok().expect("Error reading record.");
            let abp_value = record.info(b"ABP").float().ok().expect("Error reading ABP float.").expect("Missing tag ABP")[0];
            let srp_value = record.info(b"SRP").float().ok().expect("Error reading SRP float.").expect("Missing tag SRP")[0];
            let sap_value = record.info(b"SAP").float().ok().expect("Error reading SAP float.").expect("Missing tag SAP")[0];
            let sample_1_genotype: String;
            {
                //a separate scope due to conflict between mutable borrow, record.genotypes(), and immutable borrows, record.rid(), etc..
                let genotypes = record.genotypes().expect("Error reading genotypes");
                sample_1_genotype = format!("{}", genotypes.get(0));
            }
            if sample_1_genotype == "0/1" && abp_value<=abp_max {
                //sample_1_genotype could be . (uncalled), then RO and AO are undefined (a big minus number)
                //  and integer().ok() still works.
                let no_of_ref_obs = record.format(b"RO").integer().ok().expect("Error reading RO integer.")[0][0] as usize;
                let no_of_alt_obs = record.format(b"AO").integer().ok().expect("Error reading AO integer.")[0][0] as usize;
                let depth = no_of_ref_obs + no_of_alt_obs;
                if depth >= self.min_coverage && depth <= self.max_coverage {
                    one_genome_snp.no_of_good_hets += 1;
                    let ref_id = record.rid().expect("Error reading rid.") as u64;
                    let snp_key = (ref_id, record.pos() as u64);
                    let chr = String::from_utf8_lossy(vcf_header.rid2name(ref_id as u32)).to_string();
                    let one_snp = OneSNP::new(chr, snp_key.1, sample_1_genotype,
                                              depth, abp_value, srp_value, sap_value, no_of_ref_obs, no_of_alt_obs);
                    one_genome_snp.chr_start2snp.insert(snp_key, one_snp);
                }
            }
        }
        println_stderr!("{} good hets out of {} SNPs in total.", one_genome_snp.no_of_good_hets, one_genome_snp.no_of_total_records);
        return one_genome_snp;
    }

    fn intersect_snp(&self, one_genome_snp_tumor: &OneGenomeSNP, one_genome_snp_normal: &OneGenomeSNP){
        let output_f = File::create(&self.output_file_path)
            .expect(&format!("Error in creating output file {:?}", &self.output_file_path));
        let mut gz_writer = flate2::GzBuilder::new()
            .filename(self.output_file_path.file_stem().unwrap().to_str().unwrap())
            .comment("Comment")
            .write(output_f, Compression::default());
        gz_writer.write_fmt(format_args!("#abp_max_tumor={}, abp_max_normal={}, srp_max={}, \
                sap_max={}, min_coverage={}, max_coverage={}\n",
                                         self.abp_max_tumor, self.abp_max_normal, self.srp_max,
                                         self.sap_max, self.min_coverage, self.max_coverage)).unwrap();
        gz_writer.write_fmt(format_args!("#tumor snp:{:?}\n", &self.snp_file_path_tumor)).unwrap();
        gz_writer.write_fmt(format_args!("#tumor no_of_total_records: {}\n", one_genome_snp_tumor.no_of_total_records)).unwrap();
        gz_writer.write_fmt(format_args!("#tumor no_of_good_hets: {}\n", one_genome_snp_tumor.no_of_good_hets)).unwrap();
        gz_writer.write_fmt(format_args!("#normal snp:{:?}\n", &self.snp_file_path_tumor)).unwrap();
        gz_writer.write_fmt(format_args!("#normal no_of_total_records: {}\n", one_genome_snp_normal.no_of_total_records)).unwrap();
        gz_writer.write_fmt(format_args!("#normal no_of_good_hets: {}\n", one_genome_snp_normal.no_of_good_hets)).unwrap();
        gz_writer.write_fmt(format_args!("chr\tpos\ttumor_maf_normalized\ttumor_depth\t\
                tumor_maf\ttumor_ro\ttumor_ao\t\
                normal_maf\tnormal_ro\tnormal_ao\n")).unwrap();

        let mut no_of_intersect = 0usize;

        let mut snp_key_tumor_vec: Vec<(u64, u64)> = one_genome_snp_tumor.chr_start2snp.keys().map(|snp_key| *snp_key).collect();
        //let mut snp_key_tumor_vec = one_genome_snp_tumor.chr_start2snp.iter().map(|(snp_key, _)| *snp_key).collect::<Vec<(u64, u64)>>();
        //snp_key_tumor_vec.sort_by(|a, b| ((a.0<<30)+a.1).cmp(&((b.0<<30)+b.1)));
        //plus operator has precedence over bit operator. Without () for (k.0<<31), shift will cause overflow.
        snp_key_tumor_vec.sort_by_key(|k| (k.0<<31)+k.1);
        for snp_key in snp_key_tumor_vec.iter() {
            let snp_tumor = &one_genome_snp_tumor.chr_start2snp[snp_key];
            if one_genome_snp_normal.chr_start2snp.contains_key(snp_key) {
                no_of_intersect += 1;
                //let snp_tumor = &one_genome_snp_tumor.chr_start2snp[snp_key];
                let snp_normal = &one_genome_snp_normal.chr_start2snp[snp_key];
                let tumor_maf = cmp::max(snp_tumor.no_of_alt_obs, snp_tumor.no_of_ref_obs) as f32/snp_tumor.depth as f32;
                let normal_maf = cmp::max(snp_normal.no_of_alt_obs, snp_normal.no_of_ref_obs) as f32/snp_normal.depth as f32;
                let ratio_ref = snp_tumor.no_of_ref_obs as f32/snp_normal.no_of_ref_obs as f32;
                let ratio_alt  = snp_tumor.no_of_alt_obs as f32/snp_normal.no_of_alt_obs as f32;
                let tumor_maf_normalized: f32;
                if ratio_alt>ratio_ref{
                    tumor_maf_normalized = ratio_alt/(ratio_alt + ratio_ref);
                } else {
                    tumor_maf_normalized = ratio_ref/(ratio_alt + ratio_ref);
                }
                gz_writer.write_fmt(format_args!("{}\t{}\t{}\t{}\t\
                        {}\t{}\t{}\t\
                        {}\t{}\t{}\n",
                                                 snp_tumor.chr, snp_tumor.start, tumor_maf_normalized, snp_tumor.depth,
                                                 tumor_maf, snp_tumor.no_of_ref_obs, snp_tumor.no_of_alt_obs,
                                                 normal_maf, snp_normal.no_of_ref_obs, snp_normal.no_of_alt_obs)
                ).unwrap();
            }

        }

        gz_writer.write_fmt(format_args!("#no_of_intersect: {}\n", no_of_intersect)).unwrap();
        gz_writer.finish()
            .expect(&format!("ERROR finish() failure for gz_writer of {:?}.", &self.output_file_path));
        println_stderr!("{} intersect SNPs.", no_of_intersect);
    }

    pub fn run(&self){
        let one_genome_snp_tumor =
              self.read_in_het_snp(self.snp_file_path_tumor, self.abp_max_tumor);
        let one_genome_snp_normal =
              self.read_in_het_snp(self.snp_file_path_normal, self.abp_max_normal);
        self.intersect_snp(&one_genome_snp_tumor, &one_genome_snp_normal);
    }
}
