extern crate clap;
extern crate rust_htslib;
extern crate bio;

use clap::{Arg, App};
use rust_htslib::bam;
use rust_htslib::prelude::*;
use bio::io::fasta;

#[derive(Clone)]
pub struct GenomicInterval {
    pub tid: u32,
    pub chrom: String,
    // chromosome name
    pub start_pos: u32,
    // start of interval
    pub end_pos: u32,
    // end of interval (inclusive)
}

pub fn parse_target_names(bam_file: &String) -> Vec<String> {
    let bam = bam::Reader::from_path(bam_file).unwrap();
    let header_view = bam.header();
    let target_names_dec: Vec<&[u8]> = header_view.target_names();
    let mut target_names: Vec<String> = vec![];

    for t_name_dec in target_names_dec {
        let mut name_vec: Vec<char> = vec![];
        for decr in t_name_dec {
            let dec: u8 = *decr;
            name_vec.push(dec as char);
        }
        let name_string: String = name_vec.into_iter().collect();
        target_names.push(name_string);
    }

    target_names
}

pub fn u8_to_string(u: &[u8]) -> String {
    String::from_utf8(u.to_vec()).unwrap()
}

pub fn dna_vec(u: &[u8]) -> (Vec<char>) {
    let mut v: Vec<char> = Vec::with_capacity(u.len());
    for cu in u.to_ascii_uppercase() {
        let c = cu as char;
        //assert!(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
        if c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' {
            v.push(c);
        } else {
            eprintln!("Warning: Unexpected base \"{}\" encountered. Replaced with \"N\".",
                      c);
            v.push('N');
        }
    }
    v
}

pub fn get_whole_genome_intervals(bam_file: &String) -> Vec<GenomicInterval> {
    let bam = bam::Reader::from_path(bam_file).unwrap();
    let header_view = bam.header();
    let target_names_dec: Vec<&[u8]> = header_view.target_names();
    let mut intervals: Vec<GenomicInterval> = vec![];

    for (tid, t_name_dec) in target_names_dec.iter().enumerate() {
        let mut name_vec: Vec<char> = vec![];
        for decr in t_name_dec.iter() {
            let dec: u8 = *decr;
            name_vec.push(dec as char);
        }
        let name_string: String = name_vec.into_iter().collect();
        intervals.push(GenomicInterval{
            tid: tid as u32,
            chrom: name_string,
            start_pos: 0,
            end_pos: header_view.target_len(tid as u32).unwrap()-1
        });
    }

    intervals
}

// given a bam file name and a possible genomic interval,
// if the interval exists then just return a vector holding that lone interval
// otherwise, if the interval is None,
// return a vector holding GenomicIntervals representing the whole genome.
pub fn get_interval_lst(bam_file: &String, interval: &Option<GenomicInterval>) -> Vec<GenomicInterval> {
    match interval {
        &Some(ref iv) => {
            vec![iv.clone()]
        }
        &None => {
            get_whole_genome_intervals(bam_file)
        }
    }
}

// this is really ugly. TODO a less verbose implementation
pub fn parse_chrom_string(chrom_string: Option<&str>,
                           bamfile_name: &String)
                           -> Option<GenomicInterval> {
    let bam = bam::Reader::from_path(bamfile_name).unwrap();

    match chrom_string {
        Some(r) if r.contains(":") && r.contains("-") => {
            panic!("Invalid chromosome string");
        }
        Some(r) => {
            let r_str = r.to_string();
            let mut tid: u32 = 0;
            for name in bam.header().target_names() {
                if u8_to_string(name) == r_str {
                    break;
                }
                tid += 1;
            }
            if tid as usize == bam.header().target_names().len() {
                panic!("Chromosome name for region is not in BAM file.");
            }

            let tlen = bam.header().target_len(tid).unwrap();
            Some(GenomicInterval {
                tid: tid,
                chrom: r_str,
                start_pos: 0,
                end_pos: tlen - 1,
            })
        }
        None => None,
    }
}

pub fn count_mapped_reads(bam_file: &String,
                          fasta_file: &String,
                           interval: &Option<GenomicInterval>,
                           min_coverage: u32,
                           min_mapq: u8,
                           min_map_frac: f64) {

    let target_names = parse_target_names(&bam_file);

    let mut fasta = fasta::IndexedReader::from_file(&fasta_file).unwrap();

    // pileup over all covered sites
    let mut ref_seq: Vec<char> = vec![];
    let mut prev_tid = 4294967295;

    let a_str = "A".to_string();
    let c_str = "C".to_string();
    let g_str = "G".to_string();
    let t_str = "T".to_string();

    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bam_file, interval);
    let mut bam_ix = bam::IndexedReader::from_path(bam_file).unwrap();

    let mut count = 0;

    for iv in interval_lst {
        bam_ix.fetch(iv.tid as u32, iv.start_pos as u32, iv.end_pos as u32 + 1).ok().expect("Error seeking BAM file while extracting fragments.");
        let bam_pileup = bam_ix.pileup();

        for p in bam_pileup {
            let pileup = p.unwrap();

            let tid: usize = pileup.tid() as usize;
            let chrom: String = target_names[tid].clone();

            if tid != prev_tid {
                let mut ref_seq_u8: Vec<u8> = vec![];
                fasta.read_all(&chrom, &mut ref_seq_u8).expect("Failed to read fasta sequence record.");
                ref_seq = dna_vec(&ref_seq_u8);
            }

            let ref_base_str = (ref_seq[pileup.pos() as usize]).to_string();

            if ref_base_str.contains("N") {
                continue;
            }

            assert!(ref_base_str == a_str || ref_base_str == c_str || ref_base_str == g_str || ref_base_str == t_str);

            let mut depth: usize = 0;
            let mut well_mapped: usize = 0;

            // pileup the bases for a single position and count number of each base
            for alignment in pileup.alignments() {
                let record = alignment.record();

                // may be faster to implement this as bitwise operation on raw flag in the future?
                if record.is_secondary() || record.is_quality_check_failed() ||
                    record.is_duplicate() || record.is_supplementary() {
                    continue;
                }

                depth += 1;

                if record.is_unmapped() || record.mapq() < min_mapq {
                    continue;
                }

                well_mapped += 1;
            }

            let well_mapped_frac = well_mapped as f64 / depth as f64;

            if depth >= min_coverage as usize && well_mapped_frac >= min_map_frac {
                count += 1;
            }

            prev_tid = tid;
        }
    }

    println!("{}",count);
}


fn main() {
    let input_args = App::new("Map Counter")
        .version("0.1")
        .author("Peter Edge <edge.peterj@gmail.com>")
        .about("Given a bam, count the number of positions exceeding a given min coverage and \"well-mapped\" fraction.")
    .arg(Arg::with_name("Input BAM")
        .short("b")
        .long("bam")
        .value_name("BAM")
        .help("sorted, indexed BAM file.")
        .display_order(10)
        .required(true)
        .takes_value(true))
    .arg(Arg::with_name("Input FASTA")
        .short("r")
        .long("ref")
        .value_name("FASTA")
        .help("indexed fasta reference that BAM file is aligned to")
        .display_order(20)
        .required(true)
        .takes_value(true))
    .arg(Arg::with_name("Chrom")
        .short("C")
        .long("chrom")
        .value_name("string")
        .help("Chromosome to limit analysis to.")
        .display_order(30)
        .takes_value(true))
    .arg(Arg::with_name("Min coverage")
        .short("c")
        .long("min_cov")
        .value_name("int")
        .help("Minimum coverage (of reads passing filters) to consider position as a potential SNV.")
        .display_order(40)
        .required(true)
        .default_value("0"))
    .arg(Arg::with_name("Well-mapped fraction")
        .short("f")
        .long("map_frac")
        .value_name("float")
        .help("Minimum fraction of mapped reads with mapq >= MAPQ_CUTOFF.")
        .display_order(50)
        .required(true)
        .default_value("0"))
    .arg(Arg::with_name("Min mapq")
         .short("q")
         .long("min_mapq")
         .value_name("int")
         .help("Map quality cutoff (for calculating well-mapped fraction).")
         .display_order(60)
         .default_value("60"))
    .get_matches();

    let bamfile_name = input_args.value_of("Input BAM").unwrap().to_string();
    let fastafile_name = input_args.value_of("Input FASTA").unwrap().to_string();

    let interval: Option<GenomicInterval> = parse_chrom_string(input_args.value_of("Chrom"),
                                                                &bamfile_name);

    let min_mapq = input_args.value_of("Min mapq")
        .unwrap()
        .parse::<u8>()
        .expect("Argument min_mapq must be an int!");

    let min_cov = input_args.value_of("Min coverage")
        .unwrap()
        .parse::<u32>()
        .expect("Argument min_cov must be an int!");

    let min_map_frac = input_args.value_of("Well-mapped fraction")
        .unwrap()
        .parse::<f64>()
        .expect("Argument map_frac must be a positive float!");

    count_mapped_reads(&bamfile_name,
                       &fastafile_name,
                       &interval,
                       min_cov,
                       min_mapq,
                       min_map_frac);

}
