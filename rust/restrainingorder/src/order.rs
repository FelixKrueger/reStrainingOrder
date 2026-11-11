//! `reStrainingOrder`: scoring alignments against every strain and strain pair.

use std::process::ExitCode;

use crate::Tool;
use crate::cli::{self, Kind, Spec};
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
    die("reStrainingOrder is not ported yet\n")
}
