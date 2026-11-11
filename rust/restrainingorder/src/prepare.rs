//! `reStraining`: the genome preparation.
//!
//! Reads the Mouse Genomes Project VCF, keeps the positions where at least one strain has a
//! high-confidence homozygous ALT call, and writes:
//!
//! - `MGPv8_SNP_matrix_chr1.txt.gz` (or `MGPv7_...` with `--v7_VCF`): one row per kept
//!   position on chromosome 1, a 0/1 call per strain. This is what `reStrainingOrder` reads.
//! - `SNPs_directory/chr<name>.txt`: the same rows, per chromosome.
//! - `MGP_strains_N-masked/chr<name>.N-masked.fa`: every chromosome the VCF header names, with
//!   the kept positions replaced by `N`.
//! - `reStrainingOrder_SNP_filtering_report.txt` and `..._genome_preparation_report.txt`.
//!
//! As in the Perl, filtering stops at the first record that is not on chromosome `1`, so the
//! other chromosomes are written out unchanged (upper-cased). The matrix only ever needs
//! chromosome 1, and the masked genome only needs to hide the positions the matrix scores.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;

use crate::Tool;
use crate::cli::{self, Kind, Spec};
use crate::textio::{self, next_line, perl_split};
use crate::tool::{banner, die, print_help};

pub const OPTIONS: &[Spec] = &[
    ("help", Kind::Flag),
    ("versions", Kind::Flag),
    ("vcf_file", Kind::Str),
    ("nmasking", Kind::Flag),
    ("genome|reference_genome", Kind::Str),
    ("v7_VCF", Kind::Flag),
];

const GENOME_BUILD: &str = "GRCm39";
const SNP_DIR: &str = "SNPs_directory";
const MASKED_DIR: &str = "MGP_strains_N-masked";

pub fn run(args: &[String]) -> ExitCode {
    let opts = match cli::parse(args, OPTIONS) {
        Ok(o) => o,
        Err(e) => return die(format!("{e}\nPlease respecify command line options\n")),
    };
    if opts.flag("help") {
        return print_help(Tool::Prepare);
    }
    if opts.flag("versions") {
        print!("{}", banner(Tool::Prepare));
        return ExitCode::from(1); // the Perl exits 1 after the banner
    }
    match prepare(
        opts.value("vcf_file"),
        opts.value("genome"),
        opts.flag("v7_VCF"),
    ) {
        Ok(code) => code,
        Err(e) => die(format!("{e:#}\n")),
    }
}

fn prepare(vcf: Option<&str>, genome: Option<&str>, v7: bool) -> Result<ExitCode> {
    let Some(vcf) = vcf else {
        bail!(
            "\nYou need to provide a VCF file detailing SNPs positions with '--vcf_file your.file' \
             (e.g.: '--vcf mgp_REL2021_snps.vcf.gz'). Please respecify!\n"
        );
    };
    if !Path::new(vcf).exists() {
        bail!(
            "Input VCF file '{vcf}' doesn't exist in the folder. Please check filenames and try again!\n"
        );
    }
    if v7 {
        if vcf.rsplit('/').next() == Some("mgp_REL2005_snps_indels.vcf.gz") {
            eprintln!(
                "Using v7 MGP file 'mgp_REL2005_snps_indels.vcf.gz'. This file is not marked as Current_SNPs, so consider this approach experimental.\n"
            );
        } else {
            eprintln!(
                "Version 7 input file selected. Input VCF file '{vcf}' doesn't appear to be file 'mgp_REL2005_snps_indels.vcf.gz'. If something goes wrong, I won't accept responsibility....!\n"
            );
        }
    }
    let gzipped = vcf.ends_with("gz");
    let header = read_vcf_header(vcf, gzipped)?;

    let Some(genome) = genome else {
        // The Perl warns and exits successfully here.
        eprintln!(
            "Reference genome folder was not specified! Please use --reference_genome </genome/folder/>"
        );
        return Ok(ExitCode::SUCCESS);
    };
    let genome_dir = fs::canonicalize(genome).with_context(|| {
        format!("Failed to move to genome folder > {genome}/ <\n\n'reStraining --help' for more details\n")
    })?;
    if !genome_dir.is_dir() {
        bail!("Failed to move to genome folder > {genome}/ <: Not a directory\n");
    }

    eprintln!("\nSummarising reStraining parameters - Genome Preparation");
    eprintln!("{}", "=".repeat(50));
    eprintln!("Processing SNPs from VCF file:\t\t{vcf}");
    eprintln!("Reference genome:\t\t\t{}/\n", genome_dir.display());

    println!("Generating N-masking output for all of the following strains:");
    println!("{}", "=".repeat(37));
    for strain in &header.strains {
        print!("{strain}  //  ");
    }
    println!("\n{}\n", "=".repeat(37));
    println!(
        "Using the following chromosome(s):\t{}\n",
        header.chroms.join("\t")
    );

    let chromosomes = read_genome(&genome_dir)?;

    filter_vcf(vcf, gzipped, &header, v7)?;
    eprintln!("Finished filtering and writing out SNPs\n");

    let report_name = "reStrainingOrder_genome_preparation_report.txt";
    let mut report = BufWriter::new(
        File::create(report_name)
            .with_context(|| format!("Failed to write to file {report_name}"))?,
    );
    let mut new_n_total = 0usize;
    for chr in &header.chroms {
        new_n_total += mask_chromosome(chr, &chromosomes, &mut report)?;
    }
    let summary = format!(
        "\nSummary\n{new_n_total} Ns were newly introduced into the N-masked genome for all current MGP strains in total\n"
    );
    eprint!("\n{summary}\n");
    report.write_all(summary.as_bytes())?;
    report.flush()?;

    eprintln!(
        "All done. Genome(s) are now ready to be indexed with your favourite aligner!\nFYI, aligners known to work with reStrainingOrder are Bowtie2, Tophat, STAR, HISAT2, HiCUP and Bismark (STAR and HISAT2 require disabling soft-clipping, please check the manual for details)\n"
    );
    Ok(ExitCode::SUCCESS)
}

struct VcfHeader {
    /// Chromosomes named by `##contig` lines, in header order (the Perl uses hash order).
    chroms: Vec<String>,
    /// Sample columns of the `#CHROM` line.
    strains: Vec<String>,
}

fn read_vcf_header(vcf: &str, gzipped: bool) -> Result<VcfHeader> {
    let mut reader =
        textio::open(vcf, gzipped).with_context(|| format!("Failed to open file '{vcf}'"))?;
    let mut chroms: Vec<String> = Vec::new();
    let mut strains = Vec::new();
    let mut buf = Vec::new();
    while next_line(&mut *reader, &mut buf)? {
        let line = strip_eol(&buf);
        if !line.starts_with(b"#") {
            break;
        }
        if line.starts_with(b"##contig") {
            // /ID=(.+),/: greedy, so the name runs to the last comma on the line.
            let text = String::from_utf8_lossy(line);
            if let Some(start) = text.find("ID=")
                && let Some(end) = text.rfind(',')
                && end > start + 3
            {
                let chr = text[start + 3..end].to_string();
                if !chroms.contains(&chr) {
                    chroms.push(chr);
                }
            }
        } else if line.starts_with(b"#CHROM") {
            strains = perl_split(line, b'\t')
                .iter()
                .skip(9)
                .map(|s| String::from_utf8_lossy(s).into_owned())
                .collect();
        }
    }
    Ok(VcfHeader { chroms, strains })
}

/// `s/(\r|\n)//g`
fn strip_eol(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && (line[end - 1] == b'\r' || line[end - 1] == b'\n') {
        end -= 1;
    }
    &line[..end]
}

/// Read every `*.fa` file in the folder (or every `*.fasta` if there is none), upper-cased.
fn read_genome(dir: &Path) -> Result<HashMap<String, Vec<u8>>> {
    eprintln!(
        "Now reading in and storing sequence information of the genome specified in: {}/\n",
        dir.display()
    );
    let mut files = list_with_extension(dir, ".fa")?;
    if files.is_empty() {
        files = list_with_extension(dir, ".fasta")?;
    }
    if files.is_empty() {
        bail!(
            "The specified reference genome folder {}/ does not contain any sequence files in FastA format (with .fa or .fasta file extensions)",
            dir.display()
        );
    }

    let mut chromosomes: HashMap<String, Vec<u8>> = HashMap::new();
    for name in files {
        let path = dir.join(&name);
        let path_str = path.to_string_lossy().into_owned();
        let mut reader = textio::open(&path_str, false)
            .with_context(|| format!("Failed to read from sequence file {name}"))?;
        let mut buf = Vec::new();
        let mut current: Option<(String, Vec<u8>)> = None;
        let mut first = true;
        while next_line(&mut *reader, &mut buf)? {
            let line = remove_first_cr(&buf);
            let line: &[u8] = &line;
            if first || line.starts_with(b">") {
                if !line.starts_with(b">") {
                    bail!(
                        "The specified chromosome ({}) file doesn't seem to be in FASTA format as required!",
                        String::from_utf8_lossy(line)
                    );
                }
                if let Some((chr, seq)) = current.take() {
                    store(&mut chromosomes, chr, seq, &name)?;
                }
                current = Some((chromosome_name(line), Vec::new()));
                first = false;
            } else if let Some((_, seq)) = current.as_mut() {
                seq.extend(line.iter().map(u8::to_ascii_uppercase));
            }
        }
        if let Some((chr, seq)) = current.take() {
            store(&mut chromosomes, chr, seq, &name)?;
        }
    }
    println!();
    Ok(chromosomes)
}

fn store(
    chromosomes: &mut HashMap<String, Vec<u8>>,
    chr: String,
    seq: Vec<u8>,
    file: &str,
) -> Result<()> {
    println!("chr {chr} ({} bp)", seq.len());
    if chromosomes.contains_key(&chr) {
        bail!(
            "Exiting because chromosome name already exists. Please make sure all chromosomes have a unique name!"
        );
    }
    if seq.is_empty() {
        eprintln!("Chromosome {chr} in the file {file} did not contain any sequence information!");
    }
    chromosomes.insert(chr, seq);
    Ok(())
}

/// Files in `dir` ending in `ext`, sorted, without dotfiles: what Perl's `<*.fa>` matches.
fn list_with_extension(dir: &Path, ext: &str) -> Result<Vec<String>> {
    let mut names: Vec<String> = fs::read_dir(dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(ext) && !n.starts_with('.'))
        .collect();
    names.sort();
    Ok(names)
}

/// `s/\r//`: the first carriage return, wherever it is.
fn remove_first_cr(line: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    match line.iter().position(|&b| b == b'\r') {
        Some(i) => {
            let mut v = line.to_vec();
            v.remove(i);
            v.into()
        }
        None => line.into(),
    }
}

/// The name Bowtie takes from a FASTA header: the first word after `>`.
fn chromosome_name(header: &[u8]) -> String {
    let text = String::from_utf8_lossy(&header[1..]);
    // split /\s+/ keeps a leading empty field, so a header like "> 1" names chromosome "".
    text.split(|c: char| c.is_ascii_whitespace())
        .next()
        .unwrap_or("")
        .to_string()
}

#[derive(Default)]
struct FilterCounts {
    read: u64,
    skipped: u64,
    indels: u64,
    discarded: u64,
    matrix: u64,
}

fn filter_vcf(vcf: &str, gzipped: bool, header: &VcfHeader, v7: bool) -> Result<()> {
    eprintln!("Now reading from VCF file {vcf}");
    let mut reader =
        textio::open(vcf, gzipped).with_context(|| format!("Failed to open file '{vcf}'"))?;

    let matrix_name = if v7 {
        "MGPv7_SNP_matrix_chr1.txt.gz"
    } else {
        "MGPv8_SNP_matrix_chr1.txt.gz"
    };
    eprintln!("Writing a SNP matrix for chromosome 1 for all strains to {matrix_name}\n");
    let mut matrix = GzEncoder::new(
        BufWriter::new(
            File::create(matrix_name)
                .with_context(|| format!("Failed to write to file '{matrix_name}'"))?,
        ),
        Compression::default(),
    );
    write!(matrix, "Chromosome\tPosition\tREF\tALT")?;
    for strain in &header.strains {
        write!(matrix, "\t{strain}")?;
    }
    writeln!(matrix)?;

    if !Path::new(SNP_DIR).is_dir() {
        eprintln!("Folder '{SNP_DIR}' doesn't exist. Creating it for you...\n");
        fs::create_dir(SNP_DIR)
            .with_context(|| format!("Failed to created directory {SNP_DIR}"))?;
    }
    let mut per_chr: HashMap<&str, BufWriter<File>> = HashMap::new();
    for chr in &header.chroms {
        let mut fh = BufWriter::new(
            File::create(format!("{SNP_DIR}/chr{chr}.txt")).context("Couldn't open filehandle")?,
        );
        writeln!(fh, ">{chr}")?;
        per_chr.insert(chr, fh);
    }

    let fi_index = if v7 { 7 } else { 5 };
    let mut counts = FilterCounts::default();
    let mut last_chr: Option<Vec<u8>> = None;
    let mut calls: Vec<u8> = Vec::with_capacity(header.strains.len());
    let mut buf = Vec::new();

    while next_line(&mut *reader, &mut buf)? {
        let line = strip_eol(&buf);
        if line.starts_with(b"##") {
            continue;
        }
        if line.starts_with(b"#CHROM") {
            let fields = perl_split(line, b'\t');
            if !fields.contains(&&b"FORMAT"[..]) {
                bail!("Failed to extract index of field 'FORMAT'. Hmmm...");
            }
            if !fields.contains(&&b"INFO"[..]) {
                bail!("Failed to extract index of field 'INFO'. Hmmm...");
            }
            continue;
        }
        counts.read += 1;
        if counts.read % 1_000_000 == 0 {
            eprintln!("processed {} lines", counts.read);
        }

        let fields = perl_split(line, b'\t');
        let field = |i: usize| fields.get(i).copied().unwrap_or(b"");
        let (chr, pos, reference, alt, info) = (field(0), field(1), field(3), field(4), field(7));

        if last_chr.as_deref() != Some(chr) {
            match last_chr {
                None => eprintln!(
                    "Processing first chromosome (chr {})",
                    String::from_utf8_lossy(chr)
                ),
                Some(_) => eprintln!(
                    "Reached a new chromosome (chr {})",
                    String::from_utf8_lossy(chr)
                ),
            }
            last_chr = Some(chr.to_vec());
        }
        if chr != b"1" {
            eprintln!(
                "Reached a new chromosome ({}), but we are currently limiting filtering to chromosome 1",
                String::from_utf8_lossy(chr)
            );
            break;
        }
        if !all_bases(reference) {
            counts.skipped += 1;
            continue;
        }
        if v7 && contains(info, b"INDEL") {
            counts.indels += 1;
            continue;
        }
        if !all_bases(alt) {
            counts.skipped += 1;
            continue;
        }

        calls.clear();
        for sample in fields.iter().skip(9) {
            let parts = perl_split(sample, b':');
            let gt = parts.first().copied().unwrap_or(b"");
            let confident = parts.get(fi_index).copied() == Some(b"1");
            match gt {
                b"1/1" => calls.push(if confident { b'1' } else { b'0' }),
                b"2/2" | b"3/3" => {
                    eprintln!("{}", String::from_utf8_lossy(gt));
                    bail!("Should never happen");
                }
                _ => calls.push(b'0'),
            }
        }

        if !calls.contains(&b'1') {
            counts.discarded += 1;
            continue;
        }
        let mut row = Vec::with_capacity(line.len() / 4);
        for part in [chr, pos, reference, alt] {
            row.extend_from_slice(part);
            row.push(b'\t');
        }
        for (i, c) in calls.iter().enumerate() {
            if i > 0 {
                row.push(b'\t');
            }
            row.push(*c);
        }
        row.push(b'\n');

        matrix.write_all(&row)?;
        counts.matrix += 1;
        let chr_name = String::from_utf8_lossy(chr);
        match per_chr.get_mut(chr_name.as_ref()) {
            Some(fh) => fh.write_all(&row)?,
            None => bail!(
                "Can't use an undefined value as a symbol reference (no ##contig line for chromosome {chr_name})"
            ),
        }
    }

    matrix.finish()?.flush()?;
    for (_, mut fh) in per_chr {
        fh.flush()?;
    }

    let mut report = format!(
        "SNP position summary for all MGP strains (based on mouse genome build {GENOME_BUILD})\n{}\n\nPositions read in total:\t{}\n",
        "=".repeat(75),
        counts.read
    );
    if v7 {
        report += &format!(
            "{}\tPositions were INDELs (and hence skipped)\n",
            counts.indels
        );
    }
    report += &format!(
        "Positions skipped because the REF/ALT bases were not well defined:\t{}\n\
         Positions discarded as no strain had a high confidence call:\t{}\n\n\
         Positions printed to THE MATRIX in total:\t{}\n",
        counts.skipped, counts.discarded, counts.matrix
    );
    eprint!("\n{report}");
    let name = "reStrainingOrder_SNP_filtering_report.txt";
    fs::write(name, report).with_context(|| format!("Failed to write to file {name}"))?;
    Ok(())
}

/// `$s !~ /[^ATCG]/`: true for the empty string, as in the Perl.
fn all_bases(s: &[u8]) -> bool {
    s.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// N-mask one chromosome from its SNP file and write it out. Returns the number of Ns added.
fn mask_chromosome(
    chr: &str,
    chromosomes: &HashMap<String, Vec<u8>>,
    report: &mut impl Write,
) -> Result<usize> {
    eprintln!("Processing chromosome {chr}");
    // `unless ($chromosomes{$chr})`: an empty sequence counts as missing too.
    let Some(sequence) = chromosomes
        .get(chr)
        .filter(|s| !s.is_empty() && s.as_slice() != b"0")
    else {
        eprintln!(
            "\nThe chromosome name given in the VCF file was '{chr}' and was not found in the reference genome.\n\
             A rather common mistake might be that the VCF file was downloaded from Ensembl (who use chromosome names such as 1, 2, X, MT)\n\
             but the genome from UCSC (who use chromosome names such as chr1, chr2, chrX, chrM)"
        );
        eprintln!("The chromosome names in the reference genome folder were:");
        let mut names: Vec<&String> = chromosomes.keys().collect();
        names.sort();
        for c in names {
            eprintln!("{c}");
        }
        bail!(
            "[FATAL ERROR] Please ensure that the same version of the genome is used for both VCF annotations and reference genome (FastA files). Exiting...\n"
        );
    };
    let mut sequence = sequence.clone();

    let snps = read_snps(chr)?;
    let mut last_pos: i64 = 0;
    let mut new_n = 0usize;
    for (pos, reference, alt) in &snps {
        if *pos == last_pos {
            continue; // duplicate SNP
        }
        last_pos = *pos;
        let base = perl_substr(&sequence, *pos - 1);
        if base != reference.as_slice() {
            eprintln!(
                "Skipping {pos} {}/{} since the reference base didn't match",
                String::from_utf8_lossy(reference),
                String::from_utf8_lossy(alt)
            );
            continue;
        }
        let i = index(&sequence, *pos - 1).expect("perl_substr found a base");
        sequence[i] = b'N';
        new_n += 1;
    }

    eprintln!("Writing modified chromosome (N-masking)");
    if !Path::new(MASKED_DIR).is_dir() {
        fs::create_dir(MASKED_DIR)?;
    }
    let out_name = format!("{MASKED_DIR}/chr{chr}.N-masked.fa");
    let mut out = BufWriter::new(
        File::create(&out_name).with_context(|| format!("Failed to write to file {out_name}"))?,
    );
    writeln!(out, ">{chr}")?;
    // `while ($pos < length - 100)`: whole 100-base lines, then the rest (1 to 100 bases,
    // or an empty line for an empty sequence) on a line of its own.
    let mut start = 0usize;
    while start + 100 < sequence.len() {
        out.write_all(&sequence[start..start + 100])?;
        out.write_all(b"\n")?;
        start += 100;
    }
    out.write_all(&sequence[start..])?;
    out.write_all(b"\n")?;
    out.flush()?;

    eprintln!("{} SNPs total for chromosome {chr}", snps.len());
    let line = format!("{new_n} positions on chromosome {chr} were changed to 'N'\n");
    eprintln!("{line}");
    report.write_all(line.as_bytes())?;
    Ok(new_n)
}

/// Where `substr($s, $offset, 1)` reads, if anywhere. A negative offset counts from the end.
fn index(sequence: &[u8], offset: i64) -> Option<usize> {
    let len = sequence.len() as i64;
    let i = if offset < 0 { len + offset } else { offset };
    (0..len).contains(&i).then_some(i as usize)
}

fn perl_substr(sequence: &[u8], offset: i64) -> &[u8] {
    match index(sequence, offset) {
        Some(i) => &sequence[i..=i],
        None => b"",
    }
}

/// The SNPs for one chromosome, sorted by position (stable, like Perl's mergesort).
/// A SNP file row: position, REF, ALT.
type Snp = (i64, Vec<u8>, Vec<u8>);

fn read_snps(chr: &str) -> Result<Vec<Snp>> {
    if !Path::new(SNP_DIR).is_dir() {
        bail!(
            "Folder >>{SNP_DIR}<< did not exist, I'm afraid something seems to have gone wrong...\n"
        );
    }
    let file = format!("{SNP_DIR}/chr{chr}.txt");
    let mut snps = Vec::new();
    if !Path::new(&file).exists() {
        eprintln!(
            "Couldn't find SNP file for chromosome '{chr}' '{file}' didn't exist. Skipping..."
        );
        return Ok(snps);
    }
    eprintln!("Reading SNPs from file {file}");
    let mut reader = textio::open(&file, false)?;
    let mut buf = Vec::new();
    next_line(&mut *reader, &mut buf)?; // the ">chr" header
    while next_line(&mut *reader, &mut buf)? {
        let line = remove_first_cr(&buf);
        if line.is_empty() {
            continue;
        }
        let fields = perl_split(&line, b'\t');
        let get = |i: usize| fields.get(i).map(|f| f.to_vec()).unwrap_or_default();
        snps.push((perl_int(&get(1)), get(2), get(3)));
    }
    snps.sort_by_key(|s| s.0);
    Ok(snps)
}

/// Perl's numeric value of a string of digits: its leading integer, or 0.
fn perl_int(s: &[u8]) -> i64 {
    let digits: Vec<u8> = s.iter().copied().take_while(u8::is_ascii_digit).collect();
    std::str::from_utf8(&digits)
        .ok()
        .and_then(|d| d.parse().ok())
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_base_check_like_the_perl_regex() {
        assert!(all_bases(b"A"));
        assert!(all_bases(b"GT")); // two bases pass /[^ATCG]/, as in the Perl
        assert!(all_bases(b""));
        assert!(!all_bases(b"N"));
        assert!(!all_bases(b"A,C"));
    }

    #[test]
    fn chromosome_names_are_the_first_word() {
        assert_eq!(chromosome_name(b">1 dna:chromosome"), "1");
        assert_eq!(chromosome_name(b">MT"), "MT");
        assert_eq!(chromosome_name(b"> 1"), "");
    }

    #[test]
    fn substr_reads_one_base_or_nothing() {
        assert_eq!(perl_substr(b"ACGT", 0), b"A");
        assert_eq!(perl_substr(b"ACGT", 4), b"");
        assert_eq!(perl_substr(b"ACGT", -1), b"T"); // position 0 in a SNP file
    }
}
