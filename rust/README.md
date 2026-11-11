# reStrainingOrder in Rust

A port of the three Perl tools to one Rust binary. This file is the status journal: it says
what is ported and at what fidelity, and it is updated in the same PR as the work it
describes.

## Building

```sh
./rust/build-impl.sh
```

Builds `rust/target/release/restrainingorder` and exposes it under the three classic names
in `rust/target/impl/`. Called by one of those names it behaves as that tool; called as
`restrainingorder` it takes a subcommand (`prepare`, `score`, `report`).

## Fidelity

`test/run_tests.py` runs the fixture suite (see `test/README.md`) against any directory
holding the three tools. Which fixtures the Rust build passes is not a claim made here but a
checked file, `test/rust_fixtures.txt`, enforced in both directions:

```sh
python3 test/run_tests.py --impl rust/target/impl --gate test/rust_fixtures.txt
```

## Status

| Tool | State |
|---|---|
| multicall dispatch, option parsing, `--help`, `--version` | done |
| `reStraining` | done: all `prepare_*` fixtures |
| `reStrainingOrder` | done: all `order_*` fixtures but `order_skip_hybrid`, which needs the HTML report |
| `reStrainingReport` | not started |

## Verified equivalences

- `reStraining`: the matrix, the per-chromosome SNP files, the N-masked FASTA files and both
  reports are identical to the Perl for the v8 and v7 VCFs. Every filter rule is exercised:
  a REF that is not a base, a two-base REF (which the Perl's `/[^ATCG]/` lets through, and so
  does this), two ALTs, a `*` ALT, REFs that disagree with the genome, low-confidence and
  heterozygous calls, and v7 INDELs.

- The three `--version` banners and `--help` texts are byte-identical to the Perl v0.4.0
  output (captured into `banners/` and `help/`, checked by the `*_version` and `*_help`
  fixtures), and a unit test fails if the crate version moves without the banners.
- Options are resolved as Getopt::Long does by default: case-insensitive, any unambiguous
  prefix, options and file names in any order. Every option each Perl script declares is
  declared here under the same names and aliases.

- `reStrainingOrder`: the report, strain scores, hybrid scores, allele ratios and summary are
  identical to the Perl (up to the order of tied rows, which is random in the Perl) for pure,
  hybrid, BAM, plain-text matrix, Bismark single-end and paired-end input, and for two files
  in one run. Every branch of the per-read classification has a read in the fixtures.

## Differences from the Perl

- `reStraining` processes chromosomes in VCF header order. The Perl uses hash order, which
  changes from run to run; only the order of lines in the genome preparation report and on
  stderr depends on it.
- No `gzip`/`gunzip` subprocesses: gzip is read and written in-process.
- `reStrainingOrder` reads SAM and BAM itself, so samtools is not needed and `--samtools_path`
  is accepted and ignored. It writes no BAM files: the Perl converts SAM input to BAM next to
  the input, and sorts paired-end BAM by name into `*.sortedByName.bam`. Neither changes a
  score (reads are counted, so their order does not matter), so neither step is done.
- Ties in the strain and hybrid tables keep the header order, so runs are reproducible. The
  Perl's order of tied rows comes from hash order and changes between runs.
- `--output_dir` works. The Perl writes the per-file tables there but then reads them back from
  the current directory to build the summary, and dies. The HTML report is written next to the
  summary.
- An input whose name ends in neither `.sam` nor `.bam` is refused before anything is written.
  The Perl opens its outputs first, and for such a name the report file is the input itself,
  so the input is truncated before the error.
- `--verbose` is accepted but prints nothing extra.
