use std::fs::File;
use std::io::prelude::*;

const LEN: usize = 22;

static CHR_LIST: [&str;22] = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
    "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
    "chr21","chr22"];

struct Segment {
    chr_id: u32,
    start: u64,
    end: u64,
    copynumber: f64,
}

struct Chromosomes {
    segment_list: Vec<Segment>,
}

struct GenomicSegments {
    chromosomes_list: Vec<Chromosomes>,
}

impl Segment {
    fn new(chr_id: u32, start: u64, end: u64, copynumber: f64) -> Segment {
        Segment { chr_id, start, end, copynumber }
    }
}

impl Chromosomes {
    fn add(&mut self, seg: Segment) {
        let mut index = 0;
        for each_seg in self.segment_list.iter() {
            if each_seg.start > seg.start && each_seg.end > seg.end {
                break;
            }
            index += 1;
        }
        self.segment_list.insert(index,seg);
    }

    fn len(&self) -> u32 {
        self.segment_list.len() as u32
    }
}

impl GenomicSegments {
    fn new() -> GenomicSegments {
        let chromosomes_list: Vec<Chromosomes> = vec![];
        let mut genoms = GenomicSegments{chromosomes_list};
        for _ in 0..22 {
            let segment_list: Vec<Segment> = vec![];
            genoms.chromosomes_list.push(Chromosomes { segment_list });
        }
        genoms
    }

    fn add(&mut self, seg: Segment) {
        self.chromosomes_list[seg.chr_id as usize - 1].add(seg);
    }
}

pub struct RecallPrecision {
    resultfile: String,
    accurityfile: String,
    outputfile: String,
}

impl RecallPrecision {
    pub fn new(resultfile: &str, accurityfile: &str, outputfile: &str) -> RecallPrecision {
        RecallPrecision {
            resultfile: resultfile.to_string(),
            accurityfile: accurityfile.to_string(),
            outputfile: outputfile.to_string(),
        }
    }

    pub fn run(&self) {
        let mut actual = GenomicSegments::new();
        let mut accurity = GenomicSegments::new();

        read_from_actual_result(&self.resultfile, &mut actual);
        read_from_accurity_result(&self.accurityfile, &mut accurity);
        let (recall, precision) = call_recall_and_precision(&actual, &accurity);
        println!("recall: {:.*}, precision: {:.*}", 6, recall, 6, precision);

        let mut outputfile = File::create(&self.outputfile).expect("Cannot create file");
        outputfile.write(b"recall\tprecision\n");
        outputfile.write_fmt(format_args!("{:.*}\t{:.*}\n", 6, recall, 6, precision));
    }
}

fn print_genomic_segment(gs: &GenomicSegments) {
    for chromosome in gs.chromosomes_list.iter() {
        for segment in chromosome.segment_list.iter() {
            println!("chr{}: {}--{} {}", segment.chr_id, segment.start, segment.end, segment.copynumber);
        }
    }
}

fn read_from_actual_result(resultfile: &String, actual: &mut GenomicSegments)
{
    let mut f = File::open(resultfile).expect("file not found");
    let mut contents = String::new();
    f.read_to_string(&mut contents);

    let mut isheader = true;
    for line in contents.lines() {
        if isheader {
            isheader = false;
            continue;
        }
        let mut data_list = line.split('\t');
        let chr = data_list.next().unwrap();
        if CHR_LIST.contains(&&chr) == false {
            continue;
        }
        let chr_id: u32 = chr[3..].parse().unwrap();
        let start: u64 = data_list.next().unwrap().parse().unwrap();
        let end: u64 = data_list.next().unwrap().parse().unwrap();
        let copynumber: f64 = data_list.next().unwrap().parse().unwrap();
        data_list.next();
        let seg = Segment { chr_id, start, end, copynumber };
        actual.add(seg);
    }
}

fn read_from_accurity_result(accurityfile: &String, accurity: &mut GenomicSegments)
{
    let mut f = File::open(accurityfile).expect("file not found");
    let mut contents = String::new();
    f.read_to_string(&mut contents);

    let mut isheader = true;
    for line in contents.lines() {
        if isheader {
            isheader = false;
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let mut data_list = line.split('\t');
        let chr_id: u32 = data_list.next().unwrap().parse().unwrap();
        data_list.next();
        data_list.next();
        let copynumber: f64 = data_list.next().unwrap().parse().unwrap();
        data_list.next();
        data_list.next();
        data_list.next();
        data_list.next();
        data_list.next();
        data_list.next();
        let start: u64 = data_list.next().unwrap().parse().unwrap();
        let end: u64 = data_list.next().unwrap().parse().unwrap();
        if copynumber == 2.0 {
            continue;
        }
        let seg = Segment { chr_id, start, end, copynumber };
        accurity.add(seg);
    }
}

fn call_recall_and_precision(actual: &GenomicSegments, accurity: &GenomicSegments) -> (f64,f64) {
    let mut xony = GenomicSegments::new();
    y_segments_map_on_x(&actual, &accurity, &mut xony, LEN);
    let similarity = similarity_of_x_to_y(&xony, LEN);
    //print_genomic_segment(&xony);
    //println!("similarity: {}",similarity);
    let recall = similarity / sum_of_area(&actual, LEN);
    let precision = similarity / sum_of_area(&accurity, LEN);
    (recall, precision)
}

fn max(x: u64, y: u64) -> u64 {
    if x > y {
        x
    } else {
        y
    }
}

fn min(x: u64, y: u64) -> u64 {
    if x > y {
        y
    } else {
        x
    }
}

fn y_segments_map_on_x(x: &GenomicSegments, y: &GenomicSegments, xony: &mut GenomicSegments, length: usize) {
    for chr_index in 0..length {
        if x.chromosomes_list[chr_index].len() == 0 || y.chromosomes_list[chr_index].len() == 0 {
            continue;
        }
        for y_segment in y.chromosomes_list[chr_index].segment_list.iter() {
            for x_segment in x.chromosomes_list[chr_index].segment_list.iter() {
                let start = max(y_segment.start, x_segment.start);
                let end = min(y_segment.end, x_segment.end);
                if start < end + 1 {
                    let seg = Segment {
                        chr_id: chr_index as u32 + 1,
                        start,
                        end,
                        copynumber: (x_segment.copynumber - y_segment.copynumber).abs()
                    };
                    xony.add(seg);
                }
            }
        }
    }
}

fn similarity_of_x_to_y(xony: &GenomicSegments, length: usize) -> f64 {
    let mut similarity: f64 = 0.0;
    for chr_idx in 0..length {
        if xony.chromosomes_list[chr_idx].len() == 0 {
            continue;
        }
        for segment in xony.chromosomes_list[chr_idx].segment_list.iter() {
            similarity += (segment.end - segment.start + 1) as f64 * (0.0 - segment.copynumber).exp();
        }
    }
    return similarity;
}

fn sum_of_area(gs: &GenomicSegments, length: usize) -> f64 {
    let mut sum: u64 = 0;
    for chr_idx in 0..length {
        if gs.chromosomes_list[chr_idx].len() == 0 {
            continue;
        }
        for segments in gs.chromosomes_list[chr_idx].segment_list.iter() {
            sum += segments.end - segments.start + 1;
        }
    }
    return sum as f64;
}