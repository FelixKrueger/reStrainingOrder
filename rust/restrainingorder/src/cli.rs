//! Command-line parsing with Getopt::Long's defaults, which the Perl scripts rely on.
//!
//! Three of those defaults are visible to users and none is what a modern parser does:
//!
//! - `ignore_case`: `--snp_file` and `--SNP_file` are the same option. The help text of
//!   `reStrainingOrder` documents the first, the code registers the second.
//! - `auto_abbrev`: any unambiguous prefix works, so `--vcf` and `--genome` are what the
//!   user guide shows for `--vcf_file` and `--reference_genome`, and `--version` reaches the
//!   option the Perl registers as `versions`.
//! - `permute`: options and file names can be mixed in any order.
//!
//! The rules follow SNPsplit's port (`rust/snpsplit/src/optmatch.rs` there), extended with
//! aliases: `genome|reference_genome` is one option under two names, so a prefix matching
//! both is not ambiguous.

use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// A switch: `--paired`.
    Flag,
    /// `=s`: takes a value, inline (`--x=v`) or as the next argument.
    Str,
    /// `=i`: like `Str`, but the value must be an integer.
    Int,
}

/// One option: its names as Getopt::Long spells them (`"genome|reference_genome"`) and kind.
pub type Spec = (&'static str, Kind);

#[derive(Debug, Default)]
pub struct Parsed {
    values: HashMap<&'static str, String>,
    /// The non-option arguments, in order.
    pub args: Vec<String>,
}

impl Parsed {
    /// Whether the option was given. `name` is the first of its names.
    pub fn flag(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }

    /// The option's value, if it was given. The last occurrence wins, as in Getopt::Long.
    pub fn value(&self, name: &str) -> Option<&str> {
        self.values.get(name).map(String::as_str)
    }
}

/// Parse `args` (without the program name) against `specs`.
///
/// The error is the message Getopt::Long would warn with, without the trailing newline.
pub fn parse(args: &[String], specs: &[Spec]) -> Result<Parsed, String> {
    let mut parsed = Parsed::default();
    let mut rest = args.iter();

    while let Some(arg) = rest.next() {
        if arg == "--" {
            parsed.args.extend(rest.by_ref().cloned());
            break;
        }
        let Some(body) = arg.strip_prefix("--").or_else(|| arg.strip_prefix('-')) else {
            parsed.args.push(arg.clone());
            continue;
        };
        if body.is_empty() {
            parsed.args.push(arg.clone()); // a lone "-" is a file name
            continue;
        }
        let (name, inline) = match body.split_once('=') {
            Some((n, v)) => (n, Some(v.to_string())),
            None => (body, None),
        };
        let (canonical, kind) = resolve(name, specs)?;
        let value = match kind {
            Kind::Flag => {
                if inline.is_some() {
                    return Err(format!("Option {name} does not take an argument"));
                }
                String::new()
            }
            Kind::Str | Kind::Int => match inline.or_else(|| rest.next().cloned()) {
                Some(v) => v,
                None => return Err(format!("Option {name} requires an argument")),
            },
        };
        if kind == Kind::Int && !is_integer(&value) {
            return Err(format!(
                "Value \"{value}\" invalid for option {name} (number expected)"
            ));
        }
        parsed.values.insert(canonical, value);
    }
    Ok(parsed)
}

fn is_integer(s: &str) -> bool {
    let digits = s.strip_prefix(['-', '+']).unwrap_or(s);
    !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit())
}

/// Resolve a name as typed to its option: the first of its names, and its kind.
fn resolve(input: &str, specs: &[Spec]) -> Result<(&'static str, Kind), String> {
    let lowered = input.to_ascii_lowercase();
    let canonical = |spec: &Spec| spec.0.split('|').next().unwrap_or(spec.0);

    // An exact name wins over being a prefix of another one.
    for spec in specs {
        if spec.0.split('|').any(|n| n.to_ascii_lowercase() == lowered) {
            return Ok((canonical(spec), spec.1));
        }
    }

    let mut matched_specs: Vec<&Spec> = Vec::new();
    let mut matched_names: Vec<String> = Vec::new();
    for spec in specs {
        let mut hit = false;
        for n in spec.0.split('|') {
            if n.to_ascii_lowercase().starts_with(&lowered) {
                matched_names.push(n.to_ascii_lowercase());
                hit = true;
            }
        }
        if hit {
            matched_specs.push(spec);
        }
    }

    match matched_specs.as_slice() {
        [one] => Ok((canonical(one), one.1)),
        [] => Err(format!("Unknown option: {lowered}")),
        _ => {
            matched_names.sort();
            Err(format!(
                "Option {input} is ambiguous ({})",
                matched_names.join(", ")
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ORDER: &[Spec] = &[
        ("help|man", Kind::Flag),
        ("versions", Kind::Flag),
        ("output_dir", Kind::Str),
        ("SNP_file", Kind::Str),
        ("verbose", Kind::Flag),
        ("samtools_path", Kind::Str),
        ("sam", Kind::Flag),
    ];

    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn the_documented_lowercase_spelling_reaches_the_registered_one() {
        let p = parse(&args(&["--snp_file", "m.txt"]), ORDER).unwrap();
        assert_eq!(p.value("SNP_file"), Some("m.txt"));
    }

    #[test]
    fn a_prefix_works_and_an_exact_name_beats_a_longer_one() {
        let p = parse(&args(&["--version", "--sam"]), ORDER).unwrap();
        assert!(p.flag("versions"));
        assert!(p.flag("sam"));
    }

    #[test]
    fn options_and_files_mix_and_inline_values_work() {
        let p = parse(&args(&["a.bam", "--snp=m.txt", "-verbose", "b.bam"]), ORDER).unwrap();
        assert_eq!(p.args, args(&["a.bam", "b.bam"]));
        assert_eq!(p.value("SNP_file"), Some("m.txt"));
        assert!(p.flag("verbose"));
    }

    #[test]
    fn aliases_are_one_option() {
        let specs: &[Spec] = &[
            ("genome|reference_genome", Kind::Str),
            ("v7_VCF", Kind::Flag),
        ];
        let p = parse(&args(&["--reference", "g/"]), specs).unwrap();
        assert_eq!(p.value("genome"), Some("g/"));
        let p = parse(&args(&["--man"]), ORDER).unwrap();
        assert!(p.flag("help"));
    }

    #[test]
    fn errors_read_like_getopt_long() {
        assert_eq!(
            parse(&args(&["--ver"]), ORDER).unwrap_err(),
            "Option ver is ambiguous (verbose, versions)"
        );
        assert_eq!(
            parse(&args(&["--Frob"]), ORDER).unwrap_err(),
            "Unknown option: frob"
        );
        assert_eq!(
            parse(&args(&["--snp_file"]), ORDER).unwrap_err(),
            "Option snp_file requires an argument"
        );
        let int: &[Spec] = &[("candidates", Kind::Int)];
        assert_eq!(
            parse(&args(&["--candidates", "ten"]), int).unwrap_err(),
            "Value \"ten\" invalid for option candidates (number expected)"
        );
    }

    #[test]
    fn a_double_dash_ends_the_options() {
        let p = parse(&args(&["--sam", "--", "--verbose"]), ORDER).unwrap();
        assert!(p.flag("sam"));
        assert_eq!(p.args, args(&["--verbose"]));
    }
}
