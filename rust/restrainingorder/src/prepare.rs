//! `reStraining`: the genome preparation.

use std::process::ExitCode;

use crate::Tool;
use crate::cli::{self, Kind, Spec};
use crate::tool::{banner, die, print_help};

pub const OPTIONS: &[Spec] = &[
    ("help", Kind::Flag),
    ("versions", Kind::Flag),
    ("vcf_file", Kind::Str),
    ("nmasking", Kind::Flag),
    ("genome|reference_genome", Kind::Str),
    ("v7_VCF", Kind::Flag),
];

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
    die("reStraining is not ported yet\n")
}
