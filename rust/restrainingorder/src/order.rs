//! `reStrainingOrder`: scoring alignments against every strain and strain pair.
//!
//! Reads that cover an N-masked position on chromosome 1 are looked up in the SNP matrix
//! written by `reStraining`. Each covered position collects how many reads showed the REF
//! base, the ALT base, or something else, and those counts are then compared with every
//! strain's genotype (strain scores), every pair of strains (hybrid scores), and the share of
//! informative reads going to each strain of a pair (allele ratios).
//!
//! The arithmetic, the rounding and the order of the output are the Perl's, including where
//! the Perl is surprising. Two of those are worth knowing up front:
//!
//! - State carries over between input files, as the Perl's globals do: the per-position
//!   counts accumulate, and the strain tallies are added to again for each file.
//! - Only reads whose CIGAR is a single `M` are scored. The Perl has a CIGAR walk for the
//!   rest, but returns before it, so soft-clipped, spliced or gapped reads are unassignable.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result, bail};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;

use crate::Tool;
use crate::alignments::{Alignments, header_lines};
use crate::cli::{self, Kind, Spec};
use crate::textio::{self, next_line, perl_split};
use crate::tool::{banner, die, print_help};

pub const OPTIONS: &[Spec] = &[
    ("help|man", Kind::Flag),
    ("versions", Kind::Flag),
    ("output_dir", Kind::Str),
    ("SNP_file", Kind::Str),
    ("no_sorting", Kind::Flag),
    ("verbose", Kind::Flag),
    ("samtools_path", Kind::Str),
    ("sam", Kind::Flag),
    ("paired", Kind::Flag),
    ("hic", Kind::Flag),
    ("conflicting|weird", Kind::Flag),
    ("singletons", Kind::Flag),
    ("bisulfite", Kind::Flag),
    ("skip_hybrid_detection", Kind::Flag),
];

pub fn run(args: &[String]) -> ExitCode {
    let opts = match cli::parse(args, OPTIONS) {
        Ok(o) => o,
        Err(e) => return die(format!("{e}\nPlease respecify command line options\n")),
    };
    if opts.flag("help") {
        return print_help(Tool::Order);
    }
    if opts.flag("versions") {
        print!("{}", banner(Tool::Order));
        return ExitCode::SUCCESS;
    }
    if opts.args.is_empty() {
        eprintln!(
            "You need to provide one or more SAM/BAM files to start the allele-specific pipeline. Please respecify!"
        );
        return print_help(Tool::Order);
    }
    match score(&opts) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => die(format!("{e:#}\n")),
    }
}

struct Config {
    snp_file: String,
    output_dir: String,
    bisulfite: bool,
    skip_hybrid: bool,
}

fn configure(opts: &cli::Parsed) -> Result<Config> {
    let Some(snp_file) = opts.value("SNP_file") else {
        bail!(
            "You need to provide a matrix file detailing SNPs on chromosome 1 with e.g. '--SNP_file MGPv8_SNP_matrix_chr1.txt.gz'. Please respecify!"
        );
    };
    if opts.flag("no_sorting") {
        for file in &opts.args {
            if !file.ends_with(".bam") {
                bail!(
                    "The option --no_sorting requires you to supply a BAM file (ending in .bam), but at least one of the file names did not look like a BAM file: '{file}'. Please respecify!\n"
                );
            }
        }
    }
    for file in &opts.args {
        if !Path::new(file).exists() {
            bail!(
                "Input file '{file}' doesn't exist in the input folder. Please check filenames and try again!\n"
            );
        }
    }
    if opts.flag("singletons") && !opts.flag("paired") {
        eprintln!(
            "The option --singletons only makes sense for paired-end files. Simply ignoring it in single-end mode..."
        );
    }

    let output_dir = match opts.value("output_dir") {
        None | Some("") => String::new(),
        Some(dir) => {
            if !Path::new(dir).is_dir() {
                fs::create_dir(dir).with_context(|| format!("Unable to create directory {dir}"))?;
                eprintln!("Created output directory {dir}!\n");
            }
            let abs = fs::canonicalize(dir).with_context(|| format!("Failed to move to {dir}"))?;
            eprintln!(
                "Output will be written into the directory: {}/",
                abs.display()
            );
            format!("{}/", abs.display())
        }
    };

    let first = &opts.args[0];
    eprintln!("Testing if input file '{first}' looks like a Bisulfite-Seq file");
    let mut bisulfite = opts.flag("bisulfite");
    if is_bismark(first)? {
        bisulfite = true;
        eprintln!("File looks like a Bismark file. Setting '--bisulfite'...");
    } else {
        eprintln!("File looks like standard genomic file (aligned with e.g. Bowtie2 or HISAT2)");
    }

    Ok(Config {
        snp_file: snp_file.to_string(),
        output_dir,
        bisulfite,
        skip_hybrid: opts.flag("skip_hybrid_detection"),
    })
}

/// The Perl's `check_for_bs`: the last `@PG` line decides, so a Bismark file that another
/// program has since written (and added its own `@PG` to) is not recognised, there or here.
///
/// Whether the file is paired-end is not needed: in the Perl it only decides whether to sort
/// by name first, and no output depends on the order reads are scored in.
fn is_bismark(file: &str) -> Result<bool> {
    if !(file.ends_with(".sam") || file.ends_with(".bam")) {
        return Ok(false); // the main loop rejects it with the Perl's message
    }
    let reader = Alignments::open(Path::new(file))?;
    let mut bismark = false;
    for line in header_lines(reader.header())? {
        if line.starts_with("@PG") {
            bismark = line.contains("Bismark");
        }
    }
    Ok(bismark)
}

// ------------------------------------------------------------------ the SNP table

/// A genotype call from the matrix, as the Perl compares it: numerically (`== 0`, `== 1`)
/// in most places, as a string (`ne`) in the allele ratios.
const CALL_REF: u8 = 0;
const CALL_ALT: u8 = 1;
const CALL_OTHER: u8 = 2; // any other value: the Perl dies if it is ever scored
const CALL_MISSING: u8 = 3; // a row shorter than the header: undef, numerically 0

#[derive(Default)]
struct Snp {
    reference: Box<[u8]>,
    alt: Box<[u8]>,
    calls: Box<[u8]>,
    /// Reads showing REF, ALT, or another base. Zero means the Perl's key is undefined.
    counts: [u64; 3],
}

impl Snp {
    fn call(&self, index: usize) -> u8 {
        self.calls.get(index).copied().unwrap_or(CALL_MISSING)
    }
    fn covered(&self) -> bool {
        self.counts.iter().any(|&c| c > 0)
    }
}

/// `$x == 0` / `$x == 1` on a call.
fn numeric(call: u8) -> Option<u8> {
    match call {
        CALL_REF | CALL_MISSING => Some(0),
        CALL_ALT => Some(1),
        _ => None,
    }
}

const REF: usize = 0;
const ALT: usize = 1;
const OTHER: usize = 2;

type SnpTable = HashMap<Vec<u8>, HashMap<i64, Snp>>;

fn read_snps(snp_file: &str) -> Result<(Vec<String>, SnpTable, u64)> {
    eprintln!("Storing SNP positions provided in '{snp_file}'");
    let mut reader = textio::open(snp_file, snp_file.ends_with(".gz"))
        .with_context(|| format!("Failed to read from {snp_file}"))?;
    let mut buf = Vec::new();
    next_line(&mut *reader, &mut buf)?; // chomp removes the \n only
    let strains: Vec<String> = perl_split(&buf, b'\t')
        .iter()
        .skip(4)
        .map(|s| String::from_utf8_lossy(s).into_owned())
        .collect();

    let mut table = SnpTable::new();
    let mut count = 0u64;
    while next_line(&mut *reader, &mut buf)? {
        count += 1;
        buf.retain(|&b| b != b'\r' && b != b'\n');
        let fields = perl_split(&buf, b'\t');
        let field = |i: usize| fields.get(i).copied().unwrap_or(b"");
        let calls: Box<[u8]> = fields
            .iter()
            .skip(4)
            .map(|c| match *c {
                b"0" => CALL_REF,
                b"1" => CALL_ALT,
                b"" => CALL_MISSING, // "" == 0 in Perl
                _ => CALL_OTHER,
            })
            .collect();
        let snp = table
            .entry(field(0).to_vec())
            .or_default()
            .entry(leading_int(field(1)))
            .or_default();
        snp.reference = field(2).into();
        snp.alt = field(3).into();
        snp.calls = calls;
    }
    eprintln!("Stored {count} positions in total\n");
    Ok((strains, table, count))
}

/// Perl's numeric value of a string: its leading integer, or 0.
fn leading_int(s: &[u8]) -> i64 {
    let (sign, digits) = match s.first() {
        Some(b'-') => (-1, &s[1..]),
        Some(b'+') => (1, &s[1..]),
        _ => (1, s),
    };
    let end = digits
        .iter()
        .position(|b| !b.is_ascii_digit())
        .unwrap_or(digits.len());
    std::str::from_utf8(&digits[..end])
        .ok()
        .and_then(|d| d.parse::<i64>().ok())
        .unwrap_or(0)
        * sign
}

// ------------------------------------------------------------------ the run

/// What the Perl keeps in globals across input files.
struct State {
    strains: Vec<String>,
    snps: SnpTable,
    snps_stored: u64,
    /// Per strain name, in the order names first appear in the header.
    tallies: Vec<(String, Tally)>,
    snp_found: u64,
    no_snp_found: u64,
}

#[derive(Default, Clone)]
struct Tally {
    positions_covered: u64,
    agree: u64,
    disagree: u64,
}

fn score(opts: &cli::Parsed) -> Result<()> {
    let config = configure(opts)?;
    let mut state: Option<State> = None;

    for file in &opts.args {
        eprintln!("Now starting to process file <<< '{file}' >>> \n");
        let names = OutputNames::for_input(file);

        // The Perl opens its outputs before it checks the extension, so an input named
        // `x.txt` is overwritten by its own report. Checking first keeps the input.
        if !(file.ends_with(".sam") || file.ends_with(".bam")) {
            bail!(
                "Please supply a file in either SAM or BAM format (ending with .sam or .bam). File name given was '{file}'\n"
            );
        }

        let state = match &mut state {
            Some(s) => {
                eprintln!("Using already stored SNP information");
                s
            }
            None => {
                let (strains, snps, snps_stored) = read_snps(&config.snp_file)?;
                if snps.is_empty() {
                    bail!(
                        "SNP file doesn't appear to contain readable SNP positions. Please check the filename/file format and try again\n"
                    );
                }
                let mut tallies: Vec<(String, Tally)> = Vec::new();
                for s in &strains {
                    if !tallies.iter().any(|(n, _)| n == s) {
                        tallies.push((s.clone(), Tally::default()));
                    }
                }
                state.insert(State {
                    strains,
                    snps,
                    snps_stored,
                    tallies,
                    snp_found: 0,
                    no_snp_found: 0,
                })
            }
        };

        let mut report = names.report_header(config.skip_hybrid);
        report += &process_alignments(Path::new(file), &config, state)?;
        let strain_scores = strain_scores(state)?;

        let mut summary = format!("STRAIN SCORES\n{strain_scores}");
        let dir = &config.output_dir;
        write(dir, &names.strain_scores, &strain_scores)?;
        if !config.skip_hybrid {
            let (hybrids, ratios) = hybrid_scores(state)?;
            write(dir, &names.hybrid_scores, &hybrids)?;
            write(dir, &names.hybrid_ratios, &ratios)?;
            summary += &format!("\nHYBRID SCORES\n{hybrids}\nALLELE RATIOS\n{ratios}");
        }
        write(dir, &names.report, &report)?;
        write(dir, &names.summary, &summary)?;
        eprintln!("reStrainingOrder processing complete.\n");

        // The Perl runs reStrainingReport with no arguments, which reports on every summary
        // in the current directory. With --output_dir that finds nothing (the summaries are
        // elsewhere), so there the report is made for this summary, next to it.
        let report_args: Vec<String> = if dir.is_empty() {
            Vec::new()
        } else {
            vec![
                "--summary_report".into(),
                format!("{dir}{}", names.summary),
                "--dir".into(),
                dir.clone(),
            ]
        };
        let _ = crate::report::run(&report_args); // the Perl ignores its exit status too
        eprintln!("reStrainingOrder run complete...\n");
    }
    Ok(())
}

fn write(dir: &str, name: &str, content: &str) -> Result<()> {
    let path = PathBuf::from(format!("{dir}{name}"));
    fs::write(&path, content)
        .with_context(|| format!("Unable to write to file '{}'", path.display()))
}

struct OutputNames {
    input: String,
    report: String,
    strain_scores: String,
    hybrid_scores: String,
    hybrid_ratios: String,
    summary: String,
}

impl OutputNames {
    /// The Perl's `open_output_filehandles`, substitution for substitution.
    fn for_input(file: &str) -> Self {
        let base = file.rsplit('/').next().unwrap_or(file).to_string();
        let derive = |suffix: &str| sub_end(&sub_end(&base, "sam", "bam"), "bam", suffix);
        // s/sorted\.//: the first "sorted." anywhere in the name, not just a suffix.
        let unsorted = base.replacen("sorted.", "", 1);
        Self {
            report: sub_end(
                &sub_end(&unsorted, "sam", "bam"),
                "bam",
                "reStrainingOrder_report.txt",
            ),
            strain_scores: derive("reStrainingOrder.strain_scores.txt"),
            hybrid_scores: derive("reStrainingOrder.hybrid_scores.txt"),
            hybrid_ratios: derive("reStrainingOrder.hybrid_ratios.txt"),
            summary: derive("reStrainingOrder.summary_stats.txt"),
            input: base,
        }
    }

    fn report_header(&self, skip_hybrid: bool) -> String {
        let mut s = format!("Input file:\t\t\t\t\t\t'{}'\n", self.input);
        let _ = writeln!(
            s,
            "Writing pure strain compatibility scores to:\t\t'{}'",
            self.strain_scores
        );
        if !skip_hybrid {
            let _ = writeln!(
                s,
                "Writing possible hybrid compatibility scores to:\t'{}'",
                self.hybrid_scores
            );
            let _ = writeln!(
                s,
                "Writing hybrid allele ratio scores to:\t\t\t'{}'",
                self.hybrid_ratios
            );
        }
        let _ = writeln!(
            s,
            "Writing summary statistics file to:\t\t\t'{}'\n",
            self.summary
        );
        eprint!("{s}");
        s
    }
}

/// `s/from$/to/`
fn sub_end(s: &str, from: &str, to: &str) -> String {
    match s.strip_suffix(from) {
        Some(stem) => format!("{stem}{to}"),
        None => s.to_string(),
    }
}

// ------------------------------------------------------------------ per read

#[derive(Default)]
struct ReadCounts {
    count: u64,
    unmapped: u64,
    unassigned: u64,
    genome1: u64,
    genome2: u64,
    no_snp: u64,
    bizarre: u64,
    unassigned_but_ct: u64,
    ct_snp: u64,
    n_containing: u64,
    non_n: u64,
    n_deletion: u64,
    multi_n_deletion: u64,
}

fn process_alignments(file: &Path, config: &Config, state: &mut State) -> Result<String> {
    let mut reader = Alignments::open(file)?;
    eprintln!("Reading from sorted mapping file '{}'", file.display());
    let names: Vec<Vec<u8>> = reader
        .header()
        .reference_sequences()
        .keys()
        .map(|k| k.to_vec())
        .collect();

    let md_tag = Tag::new(b'M', b'D');
    let mut c = ReadCounts::default();
    let mut record = RecordBuf::default();
    while reader.read(&mut record)? {
        c.count += 1;
        if c.count % 1_000_000 == 0 {
            eprintln!("Processed {} lines so far", c.count);
        }
        if record.flags().is_unmapped() {
            c.unmapped += 1;
            continue;
        }
        let md = match record.data().get(&md_tag) {
            Some(Value::String(s)) if !s.is_empty() => s.to_vec(),
            _ => bail!(
                "Failed to extract MD:Z tag from read {}\n",
                String::from_utf8_lossy(record.name().map_or(&b"*"[..], |n| n.as_ref()))
            ),
        };

        if !md.contains(&b'N') {
            c.non_n += 1;
            c.unassigned += 1; // covers no SNP position
            continue;
        }
        let n_dels = deleted_n_runs(&md);
        if n_dels > 0 {
            c.n_deletion += 1;
            if n_dels > 1 {
                c.multi_n_deletion += 1;
            }
            c.bizarre += 1; // a deleted N is called conflicting
            continue;
        }
        c.n_containing += 1;

        let chr = record
            .reference_sequence_id()
            .and_then(|i| names.get(i))
            .map_or(&b""[..], |n| n.as_slice());
        let start = record.alignment_start().map_or(0, usize::from) as i64;
        let ops: Vec<_> = record.cigar().as_ref().to_vec();
        let simple = ops.len() == 1 && ops[0].kind() == CigarKind::Match;
        let (g1, g2, ct_unassigned) = if simple {
            let read = Read {
                chr,
                start,
                record: &record,
                bisulfite: config.bisulfite,
            };
            score_read(&read, &md, state, &mut c)?
        } else {
            (0, 0, 0) // the Perl returns before its CIGAR walk
        };

        if g1 > 0 && g2 > 0 {
            c.bizarre += 1;
        } else if g1 > g2 {
            c.genome1 += 1;
        } else if g2 > g1 {
            c.genome2 += 1;
        } else {
            if config.bisulfite && ct_unassigned > 0 {
                c.unassigned_but_ct += 1;
            } else {
                c.no_snp += 1;
            }
            c.unassigned += 1;
        }
    }
    Ok(tagging_report(&c, config, state))
}

/// `while ($md =~ /\^\D*N\D*\d+/g)`: deletions in the MD tag that include an N.
fn deleted_n_runs(md: &[u8]) -> usize {
    let mut n = 0;
    let mut i = 0;
    while i < md.len() {
        if md[i] == b'^' {
            let mut j = i + 1;
            let mut has_n = false;
            while j < md.len() && !md[j].is_ascii_digit() {
                has_n |= md[j] == b'N';
                j += 1;
            }
            // \D* is greedy but backtracks, so a run that contains an N and ends at a digit
            // matches; one that runs to the end of the string has no digit and does not.
            if has_n && j < md.len() {
                n += 1;
                while j < md.len() && md[j].is_ascii_digit() {
                    j += 1;
                }
                i = j;
                continue;
            }
        }
        i += 1;
    }
    n
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Strand {
    Ot,
    Ctot,
    Ctob,
    Ob,
}

impl Strand {
    fn top(self) -> bool {
        matches!(self, Strand::Ot | Strand::Ctot)
    }
}

/// The Bismark strand from the XR (read conversion) and XG (genome conversion) tags.
fn strand(record: &RecordBuf) -> Result<Strand> {
    let tag = |t: &[u8; 2]| match record.data().get(&Tag::new(t[0], t[1])) {
        Some(Value::String(s)) => s
            .iter()
            .copied()
            .filter(|&b| b != b'\r')
            .collect::<Vec<u8>>(),
        _ => Vec::new(),
    };
    let (xr, xg) = (tag(b"XR"), tag(b"XG"));
    Ok(match (xr.as_slice(), xg.as_slice()) {
        (b"CT", b"CT") => Strand::Ot,
        (b"GA", b"CT") => Strand::Ctot,
        (b"GA", b"GA") => Strand::Ctob,
        (b"CT", b"GA") => Strand::Ob,
        _ => bail!(
            "Unexpected combination of read and genome conversion: {} / {}\n",
            String::from_utf8_lossy(&xr),
            String::from_utf8_lossy(&xg)
        ),
    })
}

struct Read<'a> {
    chr: &'a [u8],
    start: i64,
    record: &'a RecordBuf,
    bisulfite: bool,
}

/// The Perl's `determine_N_position` for a read whose CIGAR is a single M: walk the MD tag,
/// find each N's offset in the read, and score the read's base there against the SNP.
///
/// Returns the genome 1 (REF) and genome 2 (ALT) evidence, and how many N positions were
/// set aside as bisulfite C>T (or G>A) SNPs.
fn score_read(
    read: &Read,
    md: &[u8],
    state: &mut State,
    c: &mut ReadCounts,
) -> Result<(u64, u64, u64)> {
    // split /N/, which drops trailing empty fields; the last element is what follows the
    // last N and is not visited.
    let ns = perl_split(md, b'N');
    if ns.len() < 2 {
        bail!(
            "There were fewer than 2 elements in the MD tag {}! number of elements: {}\n",
            String::from_utf8_lossy(md),
            ns.len()
        );
    }
    let (mut g1, mut g2, mut ct_unassigned) = (0u64, 0u64, 0u64);
    let mut pos: i64 = 0;
    let mut strand_of_read: Option<Strand> = None; // read from XR/XG once a SNP needs it
    let sequence = read.record.sequence().as_ref();

    for element in &ns[..ns.len() - 1] {
        if element.iter().any(|b| !b.is_ascii_digit()) {
            // split /[ATCG]/, again dropping trailing empty fields
            let parts = perl_split_any(element, b"ATCG");
            for part in &parts {
                pos += leading_int(part);
            }
            pos += parts.len() as i64 - 1;
        } else {
            pos += leading_int(element);
        }

        let genomic = pos + read.start;
        match state
            .snps
            .get_mut(read.chr)
            .and_then(|t| t.get_mut(&genomic))
        {
            Some(snp) => {
                state.snp_found += 1;
                let base = perl_substr(sequence, pos);
                match read.bisulfite {
                    true => {
                        let strand = match strand_of_read {
                            Some(s) => s,
                            None => *strand_of_read.insert(strand(read.record)?),
                        };
                        if score_bisulfite(snp, base, strand, &mut g1, &mut g2) {
                            c.ct_snp += 1;
                            ct_unassigned += 1;
                        }
                    }
                    false => {
                        if base == &*snp.reference {
                            g1 += 1;
                            snp.counts[REF] += 1;
                        } else if base == &*snp.alt {
                            g2 += 1;
                            snp.counts[ALT] += 1;
                        } else {
                            snp.counts[OTHER] += 1;
                        }
                    }
                }
            }
            None => state.no_snp_found += 1,
        }
        pos += 1; // the N itself
    }
    Ok((g1, g2, ct_unassigned))
}

/// `split /[set]/`: Perl's split with a character class, dropping trailing empty fields.
fn perl_split_any<'a>(s: &'a [u8], set: &[u8]) -> Vec<&'a [u8]> {
    let mut fields: Vec<&[u8]> = s.split(|b| set.contains(b)).collect();
    while fields.last().is_some_and(|f| f.is_empty()) {
        fields.pop();
    }
    fields
}

/// `substr($seq, $pos, 1)`: one base, or nothing past the end.
fn perl_substr(s: &[u8], pos: i64) -> &[u8] {
    let len = s.len() as i64;
    let i = if pos < 0 { len + pos } else { pos };
    if (0..len).contains(&i) {
        &s[i as usize..=i as usize]
    } else {
        b""
    }
}

/// The Perl's `score_bisulfite_SNPs`, branch for branch. Returns true when the position is a
/// C>T (or G>A) SNP on the strand where conversion makes REF and ALT indistinguishable.
///
/// Some branches leave the OTHER count alone where their neighbours raise it. That is the
/// Perl's doing, kept because the strain scores depend on it.
fn score_bisulfite(snp: &mut Snp, base: &[u8], strand: Strand, g1: &mut u64, g2: &mut u64) -> bool {
    let (r, a) = (&*snp.reference, &*snp.alt);
    let is = |s: &[u8], b: u8| s == [b];
    let mut genome1 = |snp_counts: &mut [u64; 3]| {
        *g1 += 1;
        snp_counts[REF] += 1;
    };
    let mut genome2 = |snp_counts: &mut [u64; 3]| {
        *g2 += 1;
        snp_counts[ALT] += 1;
    };
    let counts = &mut snp.counts;

    if is(r, b'C') {
        if is(a, b'T') && strand.top() {
            return true;
        }
        if strand.top() {
            if is(base, b'C') || is(base, b'T') {
                genome1(counts);
            } else if base == a {
                genome2(counts);
            } // else: irrelevant, and OTHER is not counted here
        } else if base == r {
            genome1(counts);
        } else if is(a, b'G') {
            if is(base, b'G') || is(base, b'A') {
                genome2(counts);
            } else {
                counts[OTHER] += 1;
            }
        } else if base == a {
            genome2(counts);
        } else {
            counts[OTHER] += 1;
        }
    } else if is(a, b'C') {
        if is(r, b'T') && strand.top() {
            return true;
        }
        if strand.top() {
            if is(base, b'C') || is(base, b'T') {
                genome2(counts);
            } else if base == r {
                genome1(counts);
            } else {
                counts[OTHER] += 1;
            }
        } else if base == a {
            genome2(counts);
        } else if is(r, b'G') {
            if is(base, b'G') || is(base, b'A') {
                genome1(counts);
            } else {
                counts[OTHER] += 1;
            }
        } else if base == r {
            genome1(counts);
        } else {
            counts[OTHER] += 1;
        }
    } else if is(r, b'G') {
        if is(a, b'A') && !strand.top() {
            return true;
        }
        if strand.top() {
            if base == r {
                genome1(counts);
            } else if is(a, b'C') {
                if is(base, b'C') || is(base, b'T') {
                    genome2(counts);
                } else {
                    counts[OTHER] += 1;
                }
            } else if base == a {
                genome2(counts);
            } else {
                counts[OTHER] += 1;
            }
        } else if is(base, b'G') || is(base, b'A') {
            genome1(counts);
        } else if base == a {
            genome2(counts);
        } else {
            counts[OTHER] += 1;
        }
    } else if is(a, b'G') {
        if is(r, b'A') && !strand.top() {
            return true;
        }
        if strand.top() {
            if base == a {
                genome2(counts);
            } else if is(r, b'C') {
                if is(base, b'C') || is(base, b'T') {
                    genome1(counts);
                } else {
                    counts[OTHER] += 1;
                }
            } else if base == r {
                genome1(counts);
            } else {
                counts[OTHER] += 1;
            }
        } else if is(base, b'G') || is(base, b'A') {
            genome2(counts);
        } else if base == r {
            genome1(counts);
        } else {
            counts[OTHER] += 1;
        }
    } else if base == r {
        genome1(counts);
    } else if base == a {
        genome2(counts);
    } else {
        counts[OTHER] += 1;
    }
    false
}

fn pct(n: u64, total: u64) -> String {
    if total == 0 {
        "N/A".into()
    } else {
        format!("{:.2}", n as f64 * 100.0 / total as f64)
    }
}

fn tagging_report(c: &ReadCounts, config: &Config, state: &State) -> String {
    let n = c.count;
    let mut s = String::new();
    let _ = write!(
        s,
        "\nAllele-tagging report\n{}\nProcessed {n} read alignments in total\n\
         Reads were unaligned and hence skipped: {} ({}%)\n\
         {} reads were unassignable ({}%)\n\
         {} reads were specific for genome 1 ({}%)\n\
         {} reads were specific for genome 2 ({}%)\n",
        "=".repeat(21),
        c.unmapped,
        pct(c.unmapped, n),
        c.unassigned,
        pct(c.unassigned, n),
        c.genome1,
        pct(c.genome1, n),
        c.genome2,
        pct(c.genome2, n),
    );
    if config.bisulfite {
        let _ = writeln!(
            s,
            "{} reads that were unassignable contained C>T SNPs preventing the assignment",
            c.unassigned_but_ct
        );
    }
    let _ = write!(
        s,
        "{} reads did not contain one of the expected bases at known SNP positions ({}%)\n\
         {} contained conflicting allele-specific SNPs ({}%)\n\n\n",
        c.no_snp,
        pct(c.no_snp, n),
        c.bizarre,
        pct(c.bizarre, n)
    );

    let found = state.snp_found + state.no_snp_found;
    let (pc_found, pc_not_found) = if n == 0 || found == 0 {
        ("N/A".to_string(), "N/A".to_string())
    } else {
        (pct(state.snp_found, found), pct(state.no_snp_found, found))
    };
    let _ = write!(
        s,
        "SNP coverage report\n===================\nSNP annotation file:\t{}\nSNPs stored in total:\t{}\n\
         N-containing reads:\t{}\nnon-N:\t\t\t{}\ntotal:\t\t\t{n}\n\
         Reads had a deletion of the N-masked position (and were thus dropped):\t{} ({}%)\n\
         Of which had multiple deletions of N-masked positions within the same read:\t{}\n\n\
         Of valid N containing reads,\nN was present in the list of known SNPs:\t{} ({pc_found}%)\n",
        config.snp_file,
        state.snps_stored,
        c.n_containing,
        c.non_n,
        c.n_deletion,
        pct(c.n_deletion, n),
        c.multi_n_deletion,
        state.snp_found,
    );
    if config.bisulfite {
        let _ = writeln!(
            s,
            "Positions were skipped since they involved C>T SNPs:\t{}",
            c.ct_snp
        );
    }
    let _ = write!(
        s,
        "N was not present in the list of SNPs:\t\t{} ({pc_not_found}%)\n\n",
        state.no_snp_found
    );
    eprint!("{s}");
    s
}

// ------------------------------------------------------------------ scores

fn covered(state: &State) -> impl Iterator<Item = &Snp> {
    state
        .snps
        .values()
        .flat_map(|m| m.values())
        .filter(|s| s.covered())
}

/// Add this file's evidence to the strain tallies and render the strain score table.
fn strain_scores(state: &mut State) -> Result<String> {
    eprintln!("STAGE I: Calculating strain tallies to determine agreement percentage");
    let mut delta: Vec<Tally> = vec![Tally::default(); state.strains.len()];
    let mut extra: HashMap<usize, Tally> = HashMap::new(); // rows longer than the header
    for snp in covered(state) {
        let [r, a, o] = snp.counts;
        for (i, &call) in snp.calls.iter().enumerate() {
            let t = match delta.get_mut(i) {
                Some(t) => t,
                None => extra.entry(i).or_default(),
            };
            t.positions_covered += 1;
            match numeric(call) {
                Some(0) => {
                    t.agree += r;
                    t.disagree += a + o;
                }
                Some(_) => {
                    t.disagree += r + o;
                    t.agree += a;
                }
                None => bail!("This should never happen!"),
            }
        }
    }
    for (i, d) in delta.iter().enumerate() {
        let slot = state
            .tallies
            .iter()
            .position(|(n, _)| *n == state.strains[i])
            .expect("every strain has a tally");
        let t = &mut state.tallies[slot].1;
        t.positions_covered += d.positions_covered;
        t.agree += d.agree;
        t.disagree += d.disagree;
    }
    if !extra.is_empty() {
        // `$strains[$index]` is undef past the header: those calls tally under the name "".
        let slot = match state.tallies.iter().position(|(n, _)| n.is_empty()) {
            Some(s) => s,
            None => {
                state.tallies.push((String::new(), Tally::default()));
                state.tallies.len() - 1
            }
        };
        for d in extra.values() {
            let t = &mut state.tallies[slot].1;
            t.positions_covered += d.positions_covered;
            t.agree += d.agree;
            t.disagree += d.disagree;
        }
    }

    // Only strains that ever covered a position exist in the Perl's %strain_tallies.
    let mut rows: Vec<(&str, &Tally, String)> = state
        .tallies
        .iter()
        .filter(|(_, t)| t.positions_covered > 0)
        .map(|(n, t)| {
            let total = t.agree + t.disagree;
            let agreeable = if total == 0 {
                "NA".to_string()
            } else {
                format!("{:.2}", t.agree as f64 / total as f64 * 100.0)
            };
            (n.as_str(), t, agreeable)
        })
        .collect();
    // sort { $b->{agreeable} cmp $a->{agreeable} }: a string sort, so "100.00" comes after
    // "99.50" and "NA" before everything. Ties keep header order (the Perl's are random).
    rows.sort_by(|x, y| y.2.as_bytes().cmp(x.2.as_bytes()));

    let mut out =
        String::from("Strain\tPositions covered\tAgreeing GT\tDisagreeing GT\tagreement\n");
    for (name, t, agreeable) in rows {
        let _ = writeln!(
            out,
            "{name}\t{}\t{}\t{}\t{agreeable}",
            t.positions_covered, t.agree, t.disagree
        );
    }
    Ok(out)
}

struct Hybrid {
    name: String,
    yes: u64,
    no: u64,
    percentage: String,
    s1: usize,
    s2: usize,
    count1: u64,
    count2: u64,
    pct1: String,
    pct2: String,
}

/// Every pair of strains: how well the pair explains the reads, and how the informative
/// reads split between the two. Rendered as the hybrid score and allele ratio tables.
fn hybrid_scores(state: &State) -> Result<(String, String)> {
    let positions: Vec<&Snp> = covered(state).collect();
    let n = state.strains.len();
    let mut hybrids = Vec::with_capacity(n * n.saturating_sub(1) / 2);

    for i in 0..n {
        for j in i + 1..n {
            let (mut yes, mut no, mut c1, mut c2) = (0u64, 0u64, 0u64, 0u64);
            for snp in &positions {
                let [r, a, o] = snp.counts;
                let (x, y) = (snp.call(i), snp.call(j));
                match (numeric(x), numeric(y)) {
                    (Some(0), Some(0)) => {
                        yes += r;
                        no += a + o;
                    }
                    (Some(1), Some(1)) => {
                        no += r + o;
                        yes += a;
                    }
                    (Some(_), Some(_)) => {
                        yes += r + a;
                        no += o;
                    }
                    _ => bail!("We require to have one of the above combinations!"),
                }
                // Allele ratios compare the calls as strings (`ne`).
                if x != y {
                    if numeric(x) == Some(0) {
                        c1 += r;
                        c2 += a;
                    } else if numeric(y) == Some(0) {
                        c2 += r;
                        c1 += a;
                    } else {
                        bail!("Either strain 1 or strain 2 should look like the reference!");
                    }
                }
            }
            let percentage = if yes + no > 0 {
                format!("{:.2}", yes as f64 / (yes + no) as f64 * 100.0)
            } else {
                "NA ".to_string()
            };
            let (pct1, pct2) = if c1 + c2 > 0 {
                (
                    format!("{:.1}", c1 as f64 / (c1 + c2) as f64 * 100.0),
                    format!("{:.1}", c2 as f64 / (c1 + c2) as f64 * 100.0),
                )
            } else {
                ("NA ".to_string(), "NA ".to_string())
            };
            hybrids.push(Hybrid {
                name: format!("{}/{}", state.strains[i], state.strains[j]),
                yes,
                no,
                percentage,
                s1: i,
                s2: j,
                count1: c1,
                count2: c2,
                pct1,
                pct2,
            });
        }
    }

    // sort { $b->{percentage} <=> $a->{percentage} }: numeric, so "NA " counts as 0.
    let key = |h: &Hybrid| h.percentage.trim().parse::<f64>().unwrap_or(0.0);
    hybrids.sort_by(|a, b| key(b).total_cmp(&key(a)));

    let mut table = String::from(
        "Potential Hybrid\tAgreeing\tDisagreeing\tPercentage agreement\tStrain1 Index\tStrain2 Index\n",
    );
    let mut ratios = String::from(
        "Strain1\tStrain1 Count\tStrain1 Percentage\tStrain2\tStrain2 Count\tStrain2 Percentage\n",
    );
    for h in &hybrids {
        let _ = writeln!(
            table,
            "{}\t{}\t{}\t{}\t{}\t{}",
            h.name, h.yes, h.no, h.percentage, h.s1, h.s2
        );
        let _ = writeln!(
            ratios,
            "{}\t{}\t{}\t{}\t{}\t{}",
            state.strains[h.s1], h.count1, h.pct1, state.strains[h.s2], h.count2, h.pct2
        );
    }
    Ok((table, ratios))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deleted_ns_are_counted_like_the_perl_regex() {
        assert_eq!(deleted_n_runs(b"20^N5N0N6N16"), 1);
        assert_eq!(deleted_n_runs(b"5^N20^N25"), 2);
        assert_eq!(deleted_n_runs(b"5^AC20N25"), 0);
        assert_eq!(deleted_n_runs(b"5^ANC20"), 1);
        assert_eq!(deleted_n_runs(b"10N39"), 0);
    }

    #[test]
    fn output_names_follow_the_perl_substitutions() {
        let n = OutputNames::for_input("dir/sample.sam");
        assert_eq!(n.input, "sample.sam");
        assert_eq!(n.report, "sample.reStrainingOrder_report.txt");
        assert_eq!(n.summary, "sample.reStrainingOrder.summary_stats.txt");
        // "sorted." goes wherever it is, as `s/sorted\.//` does.
        let n = OutputNames::for_input("x_possorted.sam");
        assert_eq!(n.report, "x_posreStrainingOrder_report.txt");
        assert_eq!(
            n.strain_scores,
            "x_possorted.reStrainingOrder.strain_scores.txt"
        );
    }

    #[test]
    fn md_splits_drop_trailing_empty_fields() {
        assert_eq!(perl_split_any(b"5A3", b"ATCG"), vec![&b"5"[..], b"3"]);
        assert_eq!(perl_split(b"0N2N46", b'N'), vec![&b"0"[..], b"2", b"46"]);
        assert_eq!(perl_split(b"N10", b'N'), vec![&b""[..], b"10"]);
    }

    fn snp(r: &[u8], a: &[u8]) -> Snp {
        Snp {
            reference: r.into(),
            alt: a.into(),
            ..Default::default()
        }
    }

    #[test]
    fn bisulfite_ct_snps_are_set_aside_on_the_converted_strand_only() {
        let (mut g1, mut g2) = (0, 0);
        let mut s = snp(b"C", b"T");
        assert!(score_bisulfite(&mut s, b"T", Strand::Ot, &mut g1, &mut g2));
        assert!(!score_bisulfite(&mut s, b"C", Strand::Ob, &mut g1, &mut g2));
        assert_eq!((g1, g2, s.counts), (1, 0, [1, 0, 0]));

        // A converted C reads as T on the top strand and still counts as REF.
        let mut s = snp(b"C", b"A");
        score_bisulfite(&mut s, b"T", Strand::Ctot, &mut g1, &mut g2);
        assert_eq!(s.counts, [1, 0, 0]);
        // ...and an unexpected base there is not counted as OTHER: the Perl's quirk.
        score_bisulfite(&mut s, b"G", Strand::Ot, &mut g1, &mut g2);
        assert_eq!(s.counts, [1, 0, 0]);
    }
}
