//! Which of the three tools an invocation means, and the text each prints verbatim.

use std::process::ExitCode;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tool {
    /// `reStraining`: filter the MGP VCF, write the SNP matrix and the N-masked genome.
    Prepare,
    /// `reStrainingOrder`: score alignments against every strain and strain pair.
    Order,
    /// `reStrainingReport`: turn a summary into an HTML report.
    Report,
}

/// The tool a classic-name invocation selects. The path is stripped first, since pipelines
/// call the tools by absolute path and the fixture runner through symlinks.
pub fn from_argv0(argv0: &str) -> Option<Tool> {
    match argv0.rsplit('/').next().unwrap_or(argv0) {
        "reStraining" => Some(Tool::Prepare),
        "reStrainingOrder" => Some(Tool::Order),
        "reStrainingReport" => Some(Tool::Report),
        _ => None,
    }
}

/// The tool a `restrainingorder <subcommand>` invocation selects.
pub fn from_subcommand(name: &str) -> Option<Tool> {
    match name {
        "prepare" | "reStraining" => Some(Tool::Prepare),
        "score" | "reStrainingOrder" => Some(Tool::Order),
        "report" | "reStrainingReport" => Some(Tool::Report),
        _ => None,
    }
}

/// `--version` output, captured byte for byte from the Perl. Regenerate, do not edit.
pub fn banner(tool: Tool) -> &'static str {
    match tool {
        Tool::Prepare => include_str!("../banners/reStraining.txt"),
        Tool::Order => include_str!("../banners/reStrainingOrder.txt"),
        Tool::Report => include_str!("../banners/reStrainingReport.txt"),
    }
}

/// `--help` output, captured byte for byte from the Perl. Regenerate, do not edit.
pub fn help(tool: Tool) -> &'static str {
    match tool {
        Tool::Prepare => include_str!("../help/reStraining.txt"),
        Tool::Order => include_str!("../help/reStrainingOrder.txt"),
        Tool::Report => include_str!("../help/reStrainingReport.txt"),
    }
}

/// Print the help text. Every Perl `print_helpfile` ends in `exit 1`.
pub fn print_help(tool: Tool) -> ExitCode {
    print!("{}", help(tool));
    ExitCode::from(1)
}

/// A fatal error, the way the Perl's `die` ends a run. The message goes to stderr as is.
pub fn die(message: impl std::fmt::Display) -> ExitCode {
    eprint!("{message}");
    ExitCode::from(255)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classic_names_resolve_through_a_path() {
        assert_eq!(from_argv0("/opt/bin/reStraining"), Some(Tool::Prepare));
        assert_eq!(from_argv0("reStrainingOrder"), Some(Tool::Order));
        assert_eq!(from_argv0("reStrainingReport"), Some(Tool::Report));
        assert_eq!(from_argv0("restrainingorder"), None);
    }

    #[test]
    fn every_banner_states_the_crate_version() {
        for tool in [Tool::Prepare, Tool::Order, Tool::Report] {
            assert!(
                banner(tool).contains(env!("CARGO_PKG_VERSION")),
                "{tool:?} banner does not mention {}; regenerate banners/",
                env!("CARGO_PKG_VERSION")
            );
        }
    }
}
