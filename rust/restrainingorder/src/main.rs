use std::process::ExitCode;

use restrainingorder::{Tool, order, prepare, report, tool};

fn main() -> ExitCode {
    let mut args: Vec<String> = std::env::args().collect();
    let argv0 = args.first().cloned().unwrap_or_default();

    // A classic name selects the tool; otherwise the first argument is a subcommand.
    let selected = tool::from_argv0(&argv0).or_else(|| {
        let t = args.get(1).and_then(|s| tool::from_subcommand(s));
        if t.is_some() {
            args.remove(1);
        }
        t
    });

    match selected {
        Some(Tool::Prepare) => prepare::run(&args[1..]),
        Some(Tool::Order) => order::run(&args[1..]),
        Some(Tool::Report) => report::run(&args[1..]),
        None => {
            eprintln!("Usage: restrainingorder <prepare|score|report> [options]");
            eprintln!("   or: invoke as reStraining, reStrainingOrder or reStrainingReport");
            ExitCode::from(2)
        }
    }
}
