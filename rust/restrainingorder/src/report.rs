//! `reStrainingReport`: the HTML report.

use std::process::ExitCode;

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
    die("reStrainingReport is not ported yet\n")
}
