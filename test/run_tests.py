#!/usr/bin/env python3
"""Run the fixture suite against an implementation of the three tools.

    test/run_tests.py --impl .                    # the Perl scripts in the repository root
    test/run_tests.py --impl rust/target/impl     # any directory holding the three tools
    test/run_tests.py --impl . --record           # rewrite expected/ from this implementation

Each fixture is a directory under test/fixtures/ with:

    cmd         one line: the tool and its arguments (paths are relative to the work dir)
    input/      optional, copied into the work dir before the run
    ignore      optional, glob patterns of outputs not to compare (one per line)
    stdout      optional marker: also compare what the tool prints on stdout
    expected/   what the run must produce

The run happens in a fresh directory holding a copy of test/data/ as data/. Every file the
run creates outside data/ is an output, and outputs are compared in both directions: a file
the implementation writes that expected/ lacks is a failure too. EXIT holds the exit status.

Some of what the Perl writes depends on Perl's hash order, which changes from run to run.
Those outputs are put in a canonical form before comparing, never the expected ones alone:

  * strain, hybrid and allele-ratio tables: within each run of rows that tie on the sort key,
    rows are ordered by name (the runs themselves are not moved)
  * the genome preparation report: its lines are sorted (chromosomes come out in hash order)
  * HTML reports: the date and time are blanked, and each embedded asset line (plotly.js, the
    logos) is replaced by its length and SHA-256 so expected/ stays small
  * gzip files are compared decompressed; BAM files, which only the Perl writes as a side
    effect of calling samtools, are ignored
  * an exit status other than 0 or 1 is recorded as "error": Perl's `die` exits with errno
"""

import argparse
import difflib
import fnmatch
import gzip
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(HERE, "fixtures")
DATA = os.path.join(HERE, "data")


# ------------------------------------------------------------------ canonical forms

def perl_num(s):
    """Perl's numeric value of a string: its leading number, or 0."""
    m = re.match(r"\s*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?", s)
    return float(m.group(0)) if m else 0.0


def within_ties(indices, key, name):
    """Order rows by name inside each run of consecutive rows with an equal sort key.

    Only runs are reordered, never the runs themselves, so an output that is not sorted by
    its key still differs from the expected one.
    """
    out, run = [], []
    for i in indices:
        if run and key(i) != key(run[-1]):
            out += sorted(run, key=name)
            run = []
        run.append(i)
    return out + sorted(run, key=name)


def canon_strains(lines):
    if not lines:
        return lines
    header, rows = lines[0], [line.split("\t") for line in lines[1:] if line]
    order = within_ties(range(len(rows)), lambda i: rows[i][4] if len(rows[i]) > 4 else "", lambda i: rows[i][0])
    return [header] + ["\t".join(rows[i]) for i in order]


def canon_hybrids(hybrid_lines, ratio_lines):
    """Order tied hybrid rows by name, moving each allele-ratio row with its hybrid row."""
    h_head, h_rows = hybrid_lines[:1], [line for line in hybrid_lines[1:] if line]
    r_head, r_rows = ratio_lines[:1], [line for line in ratio_lines[1:] if line]
    if r_rows and len(r_rows) != len(h_rows):
        return hybrid_lines, ratio_lines  # not pairable: compare as written
    cols = [row.split("\t") for row in h_rows]
    order = within_ties(range(len(h_rows)), lambda i: perl_num(cols[i][3]) if len(cols[i]) > 3 else 0.0,
                        lambda i: cols[i][0])
    return h_head + [h_rows[i] for i in order], r_head + ([r_rows[i] for i in order] if r_rows else [])


def canon_summary(text):
    parts = re.split(r"^(STRAIN SCORES|HYBRID SCORES|ALLELE RATIOS)\n", text, flags=re.M)
    sections = dict(zip(parts[1::2], parts[2::2]))
    strain = canon_strains(sections.get("STRAIN SCORES", "").rstrip("\n").split("\n"))
    out = ["STRAIN SCORES"] + strain
    if "HYBRID SCORES" in sections:
        h, r = canon_hybrids(sections["HYBRID SCORES"].rstrip("\n").split("\n"),
                             sections.get("ALLELE RATIOS", "").rstrip("\n").split("\n"))
        out += ["", "HYBRID SCORES"] + h + ["", "ALLELE RATIOS"] + r
    return "\n".join(out) + "\n"


def canon_html(text):
    text = re.sub(r"\d{4}-\d{2}-\d{2}", "YYYY-MM-DD", text)
    text = re.sub(r"\d{2}:\d{2}:\d{2}", "HH:MM:SS", text)
    lines = []
    for line in text.split("\n"):
        if len(line) > 1000:
            digest = hashlib.sha256(line.encode()).hexdigest()
            line = f"[embedded asset: {len(line)} bytes, sha256 {digest}]"
        lines.append(line)
    return "\n".join(lines)


def canonical(outputs):
    """outputs: {relative path: bytes}. Returns {relative path: text} in canonical form."""
    out = {}
    for path, raw in outputs.items():
        if path.endswith(".gz"):
            raw = gzip.decompress(raw)
        out[path] = raw.decode("utf-8", errors="replace")

    for path, text in list(out.items()):
        name = os.path.basename(path)
        if name.endswith("genome_preparation_report.txt"):
            out[path] = "\n".join(sorted(text.split("\n")))
        elif name.endswith(".html"):
            out[path] = canon_html(text)
        elif name.endswith(".strain_scores.txt"):
            out[path] = "\n".join(canon_strains(text.split("\n")))
        elif name.endswith(".summary_stats.txt"):
            out[path] = canon_summary(text)
        elif name.endswith(".hybrid_scores.txt"):
            ratio_path = path[: -len("hybrid_scores.txt")] + "hybrid_ratios.txt"
            h, r = canon_hybrids(text.split("\n"), out.get(ratio_path, "").split("\n"))
            out[path] = "\n".join(h)
            if ratio_path in out:
                out[ratio_path] = "\n".join(r)
    return out


# ------------------------------------------------------------------ running

def collect(root, ignore):
    outputs = {}
    for dirpath, dirnames, filenames in os.walk(root):
        rel_dir = os.path.relpath(dirpath, root)
        if rel_dir == "data" or rel_dir.startswith("data" + os.sep):
            continue
        for f in filenames:
            rel = os.path.normpath(os.path.join(rel_dir, f))
            if rel.endswith(".bam") or any(fnmatch.fnmatch(rel, pat) for pat in ignore):
                continue
            with open(os.path.join(dirpath, f), "rb") as fh:
                outputs[rel] = fh.read()
    return outputs


def run_fixture(name, impl, record, keep):
    fdir = os.path.join(FIXTURES, name)
    with open(os.path.join(fdir, "cmd")) as fh:
        argv = shlex.split(fh.read().strip())
    ignore = []
    if os.path.exists(os.path.join(fdir, "ignore")):
        with open(os.path.join(fdir, "ignore")) as fh:
            ignore = [line.strip() for line in fh if line.strip()]
    want_stdout = os.path.exists(os.path.join(fdir, "stdout"))

    work = tempfile.mkdtemp(prefix=f"rso_{name}_")
    shutil.copytree(DATA, os.path.join(work, "data"))
    if os.path.isdir(os.path.join(fdir, "input")):
        shutil.copytree(os.path.join(fdir, "input"), work, dirs_exist_ok=True)
    before = set(collect(work, []))

    tool = os.path.join(os.path.abspath(impl), argv[0])
    proc = subprocess.run([tool] + argv[1:], cwd=work, capture_output=True, timeout=600)

    outputs = {k: v for k, v in collect(work, ignore).items() if k not in before}
    code = proc.returncode
    outputs["EXIT"] = (str(code) if code in (0, 1) else "error").encode() + b"\n"
    if want_stdout:
        outputs["STDOUT"] = proc.stdout

    got = canonical(outputs)
    expected_dir = os.path.join(fdir, "expected")

    if record:
        if os.path.isdir(expected_dir):
            shutil.rmtree(expected_dir)
        for rel, text in got.items():
            dest = os.path.join(expected_dir, rel[:-3] if rel.endswith(".gz") else rel)
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with open(dest, "w") as fh:
                fh.write(text)
        shutil.rmtree(work)
        return True, []

    expected = {}
    for dirpath, _, filenames in os.walk(expected_dir):
        for f in filenames:
            rel = os.path.relpath(os.path.join(dirpath, f), expected_dir)
            with open(os.path.join(dirpath, f)) as fh:
                expected[rel] = fh.read()
    got = {(k[:-3] if k.endswith(".gz") else k): v for k, v in got.items()}

    problems = []
    for rel in sorted(set(expected) | set(got)):
        if rel not in got:
            problems.append(f"missing output: {rel}")
        elif rel not in expected:
            problems.append(f"unexpected output: {rel}")
        elif got[rel] != expected[rel]:
            diff = difflib.unified_diff(expected[rel].split("\n"), got[rel].split("\n"),
                                        f"expected/{rel}", f"got/{rel}", lineterm="", n=1)
            problems.append("\n".join(list(diff)[:40]))
    ok = not problems
    if not ok:
        problems.append("stderr:\n" + proc.stderr.decode(errors="replace")[-2000:])
    if keep or not ok:
        problems.append(f"(work dir kept: {work})")
    else:
        shutil.rmtree(work)
    return ok, problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--impl", required=True, help="directory holding reStraining, reStrainingOrder and reStrainingReport")
    ap.add_argument("--record", action="store_true", help="rewrite expected/ from this implementation")
    ap.add_argument("--list", help="file naming the fixtures to run, one per line (default: all)")
    ap.add_argument("--keep", action="store_true", help="keep work directories")
    ap.add_argument("fixtures", nargs="*")
    args = ap.parse_args()

    names = args.fixtures or sorted(os.listdir(FIXTURES))
    if args.list:
        with open(args.list) as fh:
            names = [line.split("#")[0].strip() for line in fh if line.split("#")[0].strip()]

    failed = []
    for name in names:
        ok, problems = run_fixture(name, args.impl, args.record, args.keep)
        print(f"{'ok  ' if ok else 'FAIL'} {name}")
        for p in problems:
            print("     " + p.replace("\n", "\n     "))
        if not ok:
            failed.append(name)
    print(f"\n{len(names) - len(failed)} of {len(names)} fixtures passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
