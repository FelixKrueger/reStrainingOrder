//! `reStrainingReport`: the HTML report for a `reStrainingOrder` summary.
//!
//! The template, plotly.js and the logos are compiled in from `plotly/` at the repository
//! root, so the binary needs nothing beside it. The Perl also builds an allele-ratio table
//! and a confidence table, but the template has no place for either, so they never reach
//! the page and are not built here.

use std::fs;
use std::process::ExitCode;

use anyhow::{Context, Result, bail};

use crate::Tool;
use crate::cli::{self, Kind, Spec};
use crate::tool::{banner, die, print_help};

pub const OPTIONS: &[Spec] = &[
    ("help|man", Kind::Flag),
    ("dir", Kind::Str),
    ("o|output", Kind::Str),
    ("summary_report", Kind::Str),
    ("version", Kind::Flag),
    ("candidates", Kind::Int),
    ("verbose", Kind::Flag),
];

const REPORT_VERSION: &str = "v0.4.0";
const TEMPLATE: &str = include_str!("../../../plotly/plotly_template.tpl");
const PLOTLY: &str = include_str!("../../../plotly/plot.ly");
const LOGO: &str = include_str!("../../../plotly/reStrainingOrder.logo");
const BIOINF_LOGO: &str = include_str!("../../../plotly/bioinf.logo");

pub fn run(args: &[String]) -> ExitCode {
    let opts = match cli::parse(args, OPTIONS) {
        Ok(o) => o,
        Err(e) => return die(format!("{e}\nPlease re-specify command line options\n")),
    };
    if opts.flag("help") {
        return print_help(Tool::Report);
    }
    if opts.flag("version") {
        print!("{}", banner(Tool::Report));
        return ExitCode::SUCCESS;
    }

    let output_dir = match opts.value("dir") {
        None | Some("") => String::new(),
        Some(d) if d.ends_with('/') => d.to_string(),
        Some(d) => format!("{d}/"),
    };
    let summaries = match opts.value("summary_report") {
        Some(s) => vec![s.to_string()],
        None => match summaries_in_cwd() {
            Ok(found) if found.is_empty() => {
                eprintln!(
                    "Found no potential reStrainingOrder summary reports in the current directory. Please specify a single reStrainingOrder summary report file using the option '--summary_report FILE'\n"
                );
                return print_help(Tool::Report);
            }
            Ok(found) => {
                eprintln!(
                    "Found {} summary reports in current directory. Now processing all of them.",
                    found.len()
                );
                found
            }
            Err(e) => return die(format!("{e:#}\n")),
        },
    };
    let manual = opts.value("o");
    if summaries.len() > 1 && manual.is_some() {
        return die(
            "You cannot run reStrainingReport on more than 1 file while specifying a single output file. Either lose the option -o to derive the output filenames automatically, or specify a single reStrainingOrder summary report file using the option '--summary_report FILE'\n\n",
        );
    }
    let candidates = match opts.value("candidates") {
        Some(c) => {
            let n: i64 = c.parse().unwrap_or(0);
            if !(1..=30).contains(&n) {
                return die(
                    "Please select a number of candidates between 1 and 30 (the default is 10) and try again\n",
                );
            }
            n as usize
        }
        None => 30,
    };

    for summary in &summaries {
        if let Err(e) = write_report(summary, &output_dir, manual, candidates) {
            return die(format!("{e:#}\n"));
        }
    }
    ExitCode::SUCCESS
}

/// `<*reStrainingOrder.summary_stats.txt>`: sorted, dotfiles excluded.
fn summaries_in_cwd() -> Result<Vec<String>> {
    let mut found: Vec<String> = fs::read_dir(".")?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with("reStrainingOrder.summary_stats.txt") && !n.starts_with('.'))
        .collect();
    found.sort();
    Ok(found)
}

fn write_report(
    summary: &str,
    output_dir: &str,
    manual: Option<&str>,
    candidates: usize,
) -> Result<()> {
    let name = match manual {
        Some(m) => m.to_string(),
        None => {
            let base = summary.rsplit('/').next().unwrap_or(summary);
            format!("{}.html", base.strip_suffix(".txt").unwrap_or(base))
        }
    };
    let output = format!("{output_dir}{name}");
    eprintln!("\nWriting reStrainingOrder HTML report to >> {output} <<\n");

    let mut doc = as_read_by_perl(TEMPLATE).replace("{{report_version}}", REPORT_VERSION);
    doc = inject(&doc, "{{logo_goes_here}}", &as_read_by_perl(LOGO)).unwrap_or(doc);
    // s/reStrainingOrder.summary_stats.txt/bam/: the first occurrence, path included.
    let bam_name = summary.replacen("reStrainingOrder.summary_stats.txt", "bam", 1);
    doc = doc.replace("{{filename}}", &bam_name);
    let Some(with_plotly) = inject(&doc, "{{plotly_goes_here}}", &as_read_by_perl(PLOTLY)) else {
        bail!(
            "Plot.ly incjection not working, won't be able to construct any meaningful HTML reports in this case....\n"
        );
    };
    doc = with_plotly;
    doc = inject(
        &doc,
        "{{bioinf_logo_goes_here}}",
        &as_read_by_perl(BIOINF_LOGO),
    )
    .unwrap_or(doc);
    let (date, time) = local_date_time();
    doc = doc.replace("{{date}}", &date).replace("{{time}}", &time);

    let text = fs::read(summary).with_context(|| format!("Couldn't read from file {summary}"))?;
    let tables = Tables::parse(&String::from_utf8_lossy(&text))?;
    doc = fill(doc, &tables, candidates);

    fs::write(&output, doc).with_context(|| format!("Failed to write to output file {output}"))
}

/// What the Perl's `read_report_template` makes of a file: each line with its `\r`s removed
/// and a newline after it, so a file without a final newline gains one.
fn as_read_by_perl(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 1);
    for line in text.split_inclusive('\n') {
        out.extend(
            line.strip_suffix('\n')
                .unwrap_or(line)
                .chars()
                .filter(|&c| c != '\r'),
        );
        out.push('\n');
    }
    out
}

/// `s/\{\{x\}\}.*\{\{x\}\}/$content/s`: replace from the first marker to the last.
fn inject(doc: &str, marker: &str, content: &str) -> Option<String> {
    let start = doc.find(marker)?;
    let end = doc.rfind(marker).filter(|&e| e > start)? + marker.len();
    Some(format!("{}{}{}", &doc[..start], content, &doc[end..]))
}

/// The local date and time, as the Perl's `localtime` formats them.
fn local_date_time() -> (String, String) {
    // SAFETY: localtime_r only writes into the struct we own.
    let tm = unsafe {
        let now = libc::time(std::ptr::null_mut());
        let mut tm: libc::tm = std::mem::zeroed();
        libc::localtime_r(&now, &mut tm);
        tm
    };
    (
        format!(
            "{:04}-{:02}-{:02}",
            tm.tm_year + 1900,
            tm.tm_mon + 1,
            tm.tm_mday
        ),
        format!("{:02}:{:02}:{:02}", tm.tm_hour, tm.tm_min, tm.tm_sec),
    )
}

#[derive(Default)]
struct Tables {
    /// strain, covered, agree, disagree, percent
    strains: Vec<[String; 5]>,
    /// hybrid, agree, disagree, percent (the indices are not used)
    hybrids: Vec<[String; 4]>,
    /// strain1, strain2, strain1 count, strain2 count, strain1 %, strain2 %
    ratios: Vec<[String; 6]>,
}

impl Tables {
    /// The Perl's `read_rSO_summary_report`.
    fn parse(text: &str) -> Result<Self> {
        let mut t = Tables::default();
        let mut section = 0;
        for (i, line) in text.split('\n').enumerate() {
            if i == 0 {
                if line != "STRAIN SCORES" {
                    bail!(
                        "This file does not look like a reStrainingOrder summary report. Aborting... (First line was: >{line}<\n"
                    );
                }
                section = 1;
                continue;
            }
            if line.is_empty() {
                continue;
            }
            // split /\t/ drops trailing empty fields; a missing field reads as "".
            let f: Vec<&str> = {
                let mut v: Vec<&str> = line.split('\t').collect();
                while v.last() == Some(&"") {
                    v.pop();
                }
                v
            };
            let get = |k: usize| f.get(k).copied().unwrap_or("").to_string();
            match section {
                1 if line == "HYBRID SCORES" => section = 2,
                1 if line.starts_with("Strain") => {}
                1 => t.strains.push([get(0), get(1), get(2), get(3), get(4)]),
                2 if line == "ALLELE RATIOS" => section = 3,
                2 if line.starts_with("Potential") => {}
                2 => t.hybrids.push([get(0), get(1), get(2), get(3)]),
                _ if line.starts_with("Strain1") => {}
                _ => t
                    .ratios
                    .push([get(0), get(3), get(1), get(4), get(2), get(5)]),
            }
        }
        if text.is_empty() {
            bail!(
                "This file does not look like a reStrainingOrder summary report. Aborting... (First line was: ><\n"
            );
        }
        t.ratios = t.ratios_in_hybrid_order();
        Ok(t)
    }

    /// The Perl's `adjust_allele_ratio_order`: one ratio row per hybrid row, matched by name.
    fn ratios_in_hybrid_order(&self) -> Vec<[String; 6]> {
        let mut out = Vec::new();
        for h in &self.hybrids {
            // A later ratio row with the same name overwrites an earlier one in the Perl hash.
            match self
                .ratios
                .iter()
                .rev()
                .find(|r| format!("{}/{}", r[0], r[1]) == h[0])
            {
                Some(r) => out.push(r.clone()),
                None => eprintln!("We've got a problem with the hybrid name!"),
            }
        }
        out
    }
}

/// Perl's numeric value of a string: its leading number, or 0.
fn perl_num(s: &str) -> f64 {
    let s = s.trim_start();
    let mut end = 0;
    let bytes = s.as_bytes();
    if matches!(bytes.first(), Some(b'+' | b'-')) {
        end = 1;
    }
    let digits_start = end;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    if end < bytes.len() && bytes[end] == b'.' {
        end += 1;
        while end < bytes.len() && bytes[end].is_ascii_digit() {
            end += 1;
        }
    }
    if end == digits_start || &s[digits_start..end] == "." {
        return 0.0;
    }
    s[..end].parse().unwrap_or(0.0)
}

fn fill(doc: String, t: &Tables, candidates: usize) -> String {
    let mut covered = String::from("<tr>\n<td>Different positions on chr1 covered:</td>\n");
    if let Some(first) = t.strains.first() {
        covered += &format!("\t<td>{}</td>\n\t</tr>\n", first[1]);
    }

    let mut strains = String::from(
        "<tbody>\n\t<tr>\n\t\t<th>Strain</th>\n\t\t<th>Agreeing Calls</th>\n\t\t<th>Disagreeing Calls</th>\n\t\t<th>Compatibility Percentage</th>\n\t</tr>\n</tbody>\n<tbody>\n",
    );
    let (mut strain_names, mut strain_pcts) = (String::new(), String::new());
    for [strain, _, agree, disagree, percent] in &t.strains {
        strains += &format!(
            "\t<tr>\n\t\t<td>{strain}</td>\n\t\t<td>{agree}</td>\n\t\t<td>{disagree}</td>\n\t\t<td>{percent}</td>\n\t</tr>\n"
        );
        strain_names += &format!(" '{strain}' ,");
        strain_pcts += &format!(" {percent} ,");
    }
    strains += "</tbody>\n";

    let mut hybrids = String::from(
        "<tbody>\n\t<tr>\n\t\t<th>Potential Hybrid</th>\n\t\t<th>Agreeing Calls</th>\n\t\t<th>Disagreeing Calls</th>\n\t\t<th>Compatibility %</th>\n\t\t<th>Strain1 %</th>\n\t\t<th>Strain2 %</th>\n\t</tr>\n</tbody>\n<tbody>\n",
    );
    let (mut hybrid_names, mut s1_pcts, mut s2_pcts) =
        (String::new(), String::new(), String::new());
    for (i, [hybrid, agree, disagree, percent]) in t.hybrids.iter().take(candidates).enumerate() {
        let (p1, p2) = t
            .ratios
            .get(i)
            .map_or(("", ""), |r| (r[4].as_str(), r[5].as_str()));
        hybrids += &format!(
            "\t<tr>\n\t\t<td>{hybrid}</td>\n\t\t<td>{agree}</td>\n\t\t<td>{disagree}</td>\n\t\t<td>{percent}</td>\n\t\t<td>{p1}</td>\n\t\t<td>{p2}</td>\n\t</tr>\n"
        );
        let p = perl_num(percent);
        hybrid_names += &format!(" '{hybrid}' ,");
        s1_pcts += &format!(" {:.2} ,", p * perl_num(p1) / 100.0);
        s2_pcts += &format!(" {:.2} ,", p * perl_num(p2) / 100.0);
    }
    hybrids += "</tbody>\n";

    // s/,$//: the trailing comma of each list goes.
    let trim = |s: String| s.strip_suffix(',').map(str::to_string).unwrap_or(s);
    doc.replace("{{number_of_positions_covered}}", &covered)
        .replace("{{strain_compatibility_content}}", &strains)
        .replace("{{strain_scores_strains_plotly}}", &trim(strain_names))
        .replace("{{strain_scores_percentages_plotly}}", &trim(strain_pcts))
        .replace("{{hybrid_compatibility_content}}", &hybrids)
        .replace(
            "{{hybrid_compatibility_strains_plotly}}",
            &trim(hybrid_names),
        )
        .replace("{{hybrid_percentages_strain1_plotly}}", &trim(s1_pcts))
        .replace("{{hybrid_percentages_strain2_plotly}}", &trim(s2_pcts))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn templates_gain_a_final_newline_and_lose_carriage_returns() {
        assert_eq!(as_read_by_perl("a\r\nb"), "a\nb\n");
        assert_eq!(as_read_by_perl("a\n"), "a\n");
    }

    #[test]
    fn injection_spans_the_first_to_the_last_marker() {
        assert_eq!(inject("x{{m}}old{{m}}y", "{{m}}", "new").unwrap(), "xnewy");
        assert!(inject("x{{m}}y", "{{m}}", "new").is_none());
    }

    #[test]
    fn perl_numbers() {
        assert_eq!(perl_num("75.03"), 75.03);
        assert_eq!(perl_num("NA "), 0.0);
        assert_eq!(perl_num(""), 0.0);
        assert_eq!(perl_num("12abc"), 12.0);
    }

    #[test]
    fn the_embedded_assets_are_the_ones_in_plotly() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../plotly");
        assert_eq!(fs::read_to_string(root.join("plot.ly")).unwrap(), PLOTLY);
    }
}
