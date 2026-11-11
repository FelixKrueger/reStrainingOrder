//! Where the port deliberately differs from the Perl, and what the fixtures cannot see.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test/data")
}

fn score(cwd: &Path, extra: &[&str], input: &Path) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_restrainingorder"))
        .arg("score")
        .arg("--SNP_file")
        .arg(data().join("MGPv8_SNP_matrix_chr1.txt.gz"))
        .args(extra)
        .arg(input)
        .current_dir(cwd)
        .output()
        .unwrap()
}

const OUTPUTS: &[&str] = &[
    "pure_CAST.reStrainingOrder_report.txt",
    "pure_CAST.reStrainingOrder.strain_scores.txt",
    "pure_CAST.reStrainingOrder.hybrid_scores.txt",
    "pure_CAST.reStrainingOrder.hybrid_ratios.txt",
    "pure_CAST.reStrainingOrder.summary_stats.txt",
];

/// The Perl dies with --output_dir: it collates the summary from the current directory.
#[test]
fn output_dir_gets_the_same_files_as_the_current_directory() {
    let here = tempfile::tempdir().unwrap();
    let there = tempfile::tempdir().unwrap();
    let input = data().join("pure_CAST.sam");
    assert!(score(here.path(), &[], &input).status.success());
    let out = there.path().join("results");
    let run = score(
        there.path(),
        &["--output_dir", out.to_str().unwrap()],
        &input,
    );
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    for name in OUTPUTS {
        assert_eq!(
            fs::read(here.path().join(name)).unwrap(),
            fs::read(out.join(name)).unwrap(),
            "{name}"
        );
    }
}

/// The Perl opens its outputs before checking the extension, and for `reads.txt` the derived
/// report name is `reads.txt` itself: the input is truncated to a few lines of report.
#[test]
fn an_input_that_is_neither_sam_nor_bam_is_refused_and_left_alone() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.txt");
    fs::copy(data().join("pure_CAST.sam"), &input).unwrap();
    let before = fs::read(&input).unwrap();
    assert!(!score(dir.path(), &[], &input).status.success());
    assert_eq!(fs::read(&input).unwrap(), before);
}

/// Ties are broken in header order, so two runs write identical files.
#[test]
fn runs_are_reproducible() {
    let (a, b) = (tempfile::tempdir().unwrap(), tempfile::tempdir().unwrap());
    let input = data().join("hybrid_DBA_PWK.bam");
    assert!(score(a.path(), &[], &input).status.success());
    assert!(score(b.path(), &[], &input).status.success());
    for name in OUTPUTS {
        let name = name.replace("pure_CAST", "hybrid_DBA_PWK");
        assert_eq!(
            fs::read(a.path().join(&name)).unwrap(),
            fs::read(b.path().join(&name)).unwrap(),
            "{name}"
        );
    }
}
