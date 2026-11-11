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
| `reStraining` | not started |
| `reStrainingOrder` | not started |
| `reStrainingReport` | not started |

## Verified equivalences

- The three `--version` banners and `--help` texts are byte-identical to the Perl v0.4.0
  output (captured into `banners/` and `help/`, checked by the `*_version` and `*_help`
  fixtures), and a unit test fails if the crate version moves without the banners.
- Options are resolved as Getopt::Long does by default: case-insensitive, any unambiguous
  prefix, options and file names in any order. Every option each Perl script declares is
  declared here under the same names and aliases.
