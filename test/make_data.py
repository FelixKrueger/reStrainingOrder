#!/usr/bin/env python3
"""Generate the small synthetic inputs in test/data/.

Everything is derived from one seed, so running this again reproduces the committed files
byte for byte. The data is tiny on purpose: what it has to cover is the branches, not the
scale. Each edge case below exists because a line of the Perl handles it.

    python3 test/make_data.py
"""

import gzip
import os
import random
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

STRAINS = ["129S1_SvImJ", "A_J", "CAST_EiJ", "DBA_2J", "PWK_PhJ", "WSB_EiJ"]
READ_LEN = 50
BASES = "ACGT"

rng = random.Random(20181107)


def random_seq(n):
    return "".join(rng.choice(BASES) for _ in range(n))


# ---------------------------------------------------------------- genome

chr1 = list(random_seq(4000))
chr1[1000:1005] = "NNNNN"  # genuine Ns in the reference: an N in MD that is no known SNP
chr2 = random_seq(1500)
chrmt = random_seq(250)


def fasta(name, seq, header_extra="", width=60, lower=None):
    lines = [f">{name}{header_extra}"]
    for i in range(0, len(seq), width):
        chunk = seq[i : i + width]
        if lower and lower[0] <= i < lower[1]:
            chunk = chunk.lower()  # soft-masked lines: the Perl upper-cases them
        lines.append(chunk)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- VCF

def genotype(v7):
    r = rng.random()
    if r < 0.50:
        gt, fi = "0/0", "1"
    elif r < 0.75:
        gt, fi = "1/1", "1"
    elif r < 0.82:
        gt, fi = "1/1", "0"  # low confidence: must not count as a call
    elif r < 0.88:
        gt, fi = "0/1", "0"
    elif r < 0.94:
        gt, fi = "0/0", "0"
    else:
        gt, fi = "./.", "."
    if v7:
        return f"{gt}:0,0,0:9:0,0:0,0:0,9:30:{fi}"
    return f"{gt}:255,30,0:9:0,9:30:{fi}"


def vcf_records(v7):
    """Records for chr1, then a few on chr2 (where the Perl stops filtering)."""
    positions = sorted(rng.sample([p for p in range(10, 3990) if not 995 <= p <= 1010], 160))
    records = []
    for i, pos in enumerate(positions):
        ref = chr1[pos - 1]
        alt = rng.choice([b for b in BASES if b != ref])
        info = "DP=99"
        if i == 7:
            ref = "N"  # REF not a base: skipped
        elif i == 13:
            ref, alt = ref + "T", alt  # multi-base REF of bases: kept, then never matches the genome
        elif i == 21:
            alt = alt + "," + rng.choice([b for b in BASES if b not in (ref, alt)])  # two ALTs: skipped
        elif i == 34:
            alt = "*"  # ALT not a base: skipped
        elif i in (40, 77):
            ref = rng.choice([b for b in BASES if b != chr1[pos - 1]])  # REF disagrees with the genome
            alt = rng.choice([b for b in BASES if b != ref])
        elif v7 and i in (50, 51):
            info = "INDEL;DP=12"  # v7 only: INDELs are skipped
        calls = [genotype(v7) for _ in STRAINS]
        if i == 13:
            calls[0] = calls[0].split(":", 1)[0].replace(calls[0].split(":", 1)[0], "1/1") + ":" + calls[0].split(":", 1)[1][:-1] + "1"
        if i == 60:
            calls = ["1/1:255,30,0:9:0,9:30:1" if not v7 else "1/1:0,0,0:9:0,0:0,0:0,9:30:1" for _ in STRAINS]
        records.append(("1", pos, ref, alt, info, calls))
    # One position inside the genuine-N stretch: REF never matches, so it is not masked.
    calls = [("1/1:255,30,0:9:0,9:30:1" if not v7 else "1/1:0,0,0:9:0,0:0,0:0,9:30:1") for _ in STRAINS]
    records.append(("1", 1003, "A", "G", "DP=99", calls))
    records.sort(key=lambda r: r[1])
    for pos in (100, 200, 300):
        ref = chr2[pos - 1]
        alt = rng.choice([b for b in BASES if b != ref])
        records.append(("2", pos, ref, alt, "DP=99", [genotype(v7) for _ in STRAINS]))
    return records


def write_vcf(path, v7):
    fmt = "GT:PL:DP:ADF:ADR:AD:GQ:FI" if v7 else "GT:PL:DP:AD:GQ:FI"
    out = ["##fileformat=VCFv4.2", "##source=reStrainingOrder test data"]
    out += [f"##contig=<ID=1,length={len(chr1)}>", f"##contig=<ID=2,length={len(chr2)}>",
            f"##contig=<ID=MT,length={len(chrmt)}>"]
    out.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + STRAINS))
    records = vcf_records(v7)
    for chrom, pos, ref, alt, info, calls in records:
        out.append("\t".join([chrom, str(pos), ".", ref, alt, "99", "PASS", info, fmt] + calls))
    with gzip.GzipFile(path, "wb", mtime=0) as fh:
        fh.write(("\n".join(out) + "\n").encode())
    return records


# ---------------------------------------------------------------- matrix, as reStraining builds it

def matrix_rows(records, v7):
    rows = []
    for chrom, pos, ref, alt, info, calls in records:
        if chrom != "1":
            break
        if any(c not in BASES for c in ref):  # the Perl's test: "AT" passes, as it does there
            continue
        if v7 and "INDEL" in info:
            continue
        if any(c not in BASES for c in alt):
            continue
        bits = []
        for c in calls:
            f = c.split(":")
            gt, fi = f[0], f[-1]
            bits.append("1" if gt == "1/1" and fi == "1" else "0")
        if "1" in bits:
            rows.append((pos, ref, alt, bits))
    return rows


# ---------------------------------------------------------------- reads

def masked_chr1(rows):
    seq = chr1[:]
    for pos, ref, _, _ in rows:
        if seq[pos - 1] == ref:
            seq[pos - 1] = "N"
    return seq


def haplotype(rows, strain):
    idx = STRAINS.index(strain)
    seq = chr1[:]
    for pos, ref, alt, bits in rows:
        seq[pos - 1] = alt if bits[idx] == "1" else chr1[pos - 1]
    return seq


def md_for(read, ref):
    """MD:Z for an ungapped alignment of `read` against `ref` (same length)."""
    md, run = [], 0
    for r, g in zip(read, ref):
        if r == g and g != "N":
            run += 1
        else:
            md.append(str(run))
            md.append(g)
            run = 0
    md.append(str(run))
    return "".join(md)


def sam_line(name, flag, chrom, pos, cigar, seq, md, extra=()):
    qual = "I" * len(seq)
    fields = [name, str(flag), chrom, str(pos), "42", cigar, "*", "0", "0", seq, qual]
    if md is not None:
        fields.append(f"MD:Z:{md}")
    fields += list(extra)
    return "\t".join(fields)


def mutate(seq, rate):
    out = []
    for b in seq:
        if rng.random() < rate:
            out.append(rng.choice([x for x in BASES if x != b]))
        else:
            out.append(b)
    return out


def simulate(rows, strains, n, prefix, error=0.01):
    masked = masked_chr1(rows)
    haps = {s: haplotype(rows, s) for s in strains}
    lines = []
    for i in range(n):
        s = strains[i % len(strains)]
        start = rng.randint(1, len(chr1) - READ_LEN + 1)
        read = mutate(haps[s][start - 1 : start - 1 + READ_LEN], error)
        ref = masked[start - 1 : start - 1 + READ_LEN]
        flag = rng.choice([0, 16])
        lines.append(sam_line(f"{prefix}_{i:04d}", flag, "1", start, f"{READ_LEN}M", "".join(read), md_for(read, ref)))
    return lines


def edge_cases(rows):
    """One read per branch of the Perl's per-read classification."""
    masked = masked_chr1(rows)
    npos = [i + 1 for i, b in enumerate(masked) if b == "N" and not 1001 <= i + 1 <= 1005]
    lines = []

    lines.append(sam_line("edge_unmapped", 4, "*", 0, "*", "ACGT" * 5, None))

    # Soft clip: a CIGAR that is not a single M, so the read is unassignable.
    p = npos[5] - 10
    ref = masked[p - 1 : p - 1 + 47]
    read = ["T", "T", "T"] + ["A" if b == "N" else b for b in ref]
    lines.append(sam_line("edge_softclip", 0, "1", p, "3S47M", "".join(read), md_for(read[3:], ref)))

    # Deletion of an N-masked position: conflicting.
    p = npos[10] - 20
    ref = masked[p - 1 : p - 1 + 51]
    read = ref[:20] + ref[21:]
    read = ["C" if b == "N" else b for b in read]
    md = md_for(read[:20], ref[:20]) + "^N" + md_for(read[20:], ref[21:])
    lines.append(sam_line("edge_del_one_N", 0, "1", p, "20M1D30M", "".join(read), md.replace("^N0", "^N0")))

    # Two N deletions in one read.
    a, b = npos[15], npos[16]
    if b - a > 3 and b - a < 40:
        p = a - 5
        ref = masked[p - 1 : p - 1 + 60]
        da, db = a - p, b - p
        read = ref[:da] + ref[da + 1 : db] + ref[db + 1 :]
        read = read[:50]
        read = ["G" if x == "N" else x for x in read]
        m1 = md_for(read[:da], ref[:da])
        m2 = md_for(read[da : db - 1], ref[da + 1 : db])
        m3 = md_for(read[db - 1 :], ref[db + 1 : db + 1 + 50 - (db - 1)])
        lines.append(sam_line("edge_del_two_N", 0, "1", p, f"{da}M1D{db - 1 - da}M1D{50 - (db - 1)}M", "".join(read), f"{m1}^N{m2}^N{m3}"))

    # Insertion: not a simple match, unassignable.
    p = npos[20] - 10
    ref = masked[p - 1 : p - 1 + 48]
    read = ref[:20] + ["A", "C"] + ref[20:]
    read = ["T" if x == "N" else x for x in read]
    lines.append(sam_line("edge_insertion", 0, "1", p, "20M2I28M", "".join(read), md_for(read[:20] + read[22:], ref)))

    # A read across the genuine N stretch: Ns in MD that are not known SNPs.
    p = 990
    ref = masked[p - 1 : p - 1 + 50]
    read = ["A" if x == "N" else x for x in ref]
    lines.append(sam_line("edge_genuine_N", 0, "1", p, "50M", "".join(read), md_for(read, ref)))

    # Mismatch letters before an N: exercises the [ATCG] split in the MD walk.
    p = npos[25] - 30
    ref = masked[p - 1 : p - 1 + 50]
    read = list(ref)
    for k in (3, 11):
        if read[k] != "N":
            read[k] = rng.choice([x for x in BASES if x != read[k]])
    for row in rows:
        if p <= row[0] < p + 50 and masked[row[0] - 1] == "N":
            read[row[0] - p] = row[2]
    read = [rng.choice(BASES) if x == "N" else x for x in read]
    lines.append(sam_line("edge_mismatches", 16, "1", p, "50M", "".join(read), md_for(read, ref)))

    # Both alleles in one read: conflicting by allele.
    pair = None
    for k in range(len(rows) - 1):
        p1, p2 = rows[k][0], rows[k + 1][0]
        if masked[p1 - 1] == "N" and masked[p2 - 1] == "N" and 5 < p2 - p1 < 40:
            pair = (rows[k], rows[k + 1])
            break
    if pair:
        p = pair[0][0] - 5
        ref = masked[p - 1 : p - 1 + 50]
        read = list(chr1[p - 1 : p - 1 + 50])
        read[pair[0][0] - p] = pair[0][1]  # REF allele
        read[pair[1][0] - p] = pair[1][2]  # ALT allele
        read = [rng.choice(BASES) if x == "N" else x for x in read]
        lines.append(sam_line("edge_both_alleles", 0, "1", p, "50M", "".join(read), md_for(read, ref)))

    # Reads on chr2: no Ns there (the Perl only masks chromosome 1's SNPs).
    for k, p in enumerate((50, 700)):
        seq = chr2[p - 1 : p - 1 + 50]
        lines.append(sam_line(f"edge_chr2_{k}", 0, "2", p, "50M", seq, "50"))

    # MD starting with N: the empty leading element of split /N/.
    first = next(r for r in rows if masked[r[0] - 1] == "N" and r[0] > 1)
    p = first[0]
    ref = masked[p - 1 : p - 1 + 50]
    read = [first[2]] + [("A" if x == "N" else x) for x in ref[1:]]
    lines.append(sam_line("edge_md_leading_N", 0, "1", p, "50M", "".join(read), md_for(read, ref)))
    return lines


# ---------------------------------------------------------------- bisulfite

def convert(seq, frm, to, rate=0.8):
    return [to if b == frm and rng.random() < rate else b for b in seq]


def bisulfite_reads(rows, strain, n, prefix, paired):
    masked = masked_chr1(rows)
    hap = haplotype(rows, strain)
    # (strand, flag R1, flag R2, XR R1, XR R2, XG, genome conversion)
    strands = [("OT", 99, 147, "CT", "GA", "CT"), ("OB", 83, 163, "CT", "GA", "GA"),
               ("CTOT", 147, 99, "GA", "CT", "CT"), ("CTOB", 163, 83, "GA", "CT", "GA")]
    lines = []
    for i in range(n):
        strand, f1, f2, xr1, xr2, xg = strands[i % 4] if paired else strands[i % 4]
        conv = ("C", "T") if xg == "CT" else ("G", "A")
        mates = [(f1, xr1)] + ([(f2, xr2)] if paired else [])
        start = rng.randint(1, len(chr1) - 2 * READ_LEN)
        for m, (flag, xr) in enumerate(mates):
            s = start + m * 30
            read = convert(hap[s - 1 : s - 1 + READ_LEN], *conv)
            ref = convert(masked[s - 1 : s - 1 + READ_LEN], *conv, rate=0)  # MD against the unconverted reference
            ref_conv = [conv[1] if (a == conv[0] and b == conv[1]) else a for a, b in zip(ref, read)]
            if not paired:
                flag = {"OT": 0, "OB": 16, "CTOT": 16, "CTOB": 0}[strand]
            lines.append(sam_line(f"{prefix}_{i:04d}", flag, "1", s, f"{READ_LEN}M", "".join(read),
                                  md_for(read, ref_conv), ("NM:i:0", f"XM:Z:{'.' * READ_LEN}", f"XR:Z:{xr}", f"XG:Z:{xg}")))
    return lines


# ---------------------------------------------------------------- write

def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def sam(header_pg, lines):
    head = ["@HD\tVN:1.0\tSO:unsorted", f"@SQ\tSN:1\tLN:{len(chr1)}", f"@SQ\tSN:2\tLN:{len(chr2)}",
            f"@SQ\tSN:MT\tLN:{len(chrmt)}"]
    if header_pg:
        head.append(header_pg)
    return "\n".join(head + lines) + "\n"


def main():
    os.makedirs(os.path.join(DATA, "genome"), exist_ok=True)
    write(os.path.join(DATA, "genome", "1.fa"), fasta("1", "".join(chr1), " dna:chromosome chromosome:GRCm39:1", lower=(1980, 2160)))
    write(os.path.join(DATA, "genome", "other.fa"), fasta("2", chr2, " dna:chromosome") + fasta("MT", chrmt, width=70))
    # A genome missing chromosomes the VCF names: the Perl stops with a FATAL ERROR.
    os.makedirs(os.path.join(DATA, "genome_chr1"), exist_ok=True)
    shutil.copy(os.path.join(DATA, "genome", "1.fa"), os.path.join(DATA, "genome_chr1", "1.fa"))

    records = write_vcf(os.path.join(DATA, "mgp_REL2021_snps.vcf.gz"), v7=False)
    rng.seed(7)
    records7 = write_vcf(os.path.join(DATA, "mgp_REL2005_snps_indels.vcf.gz"), v7=True)

    rows = matrix_rows(records, v7=False)
    header = "\t".join(["Chromosome", "Position", "REF", "ALT"] + STRAINS)
    body = "".join("\t".join(["1", str(p), r, a] + bits) + "\n" for p, r, a, bits in rows)
    with gzip.GzipFile(os.path.join(DATA, "MGPv8_SNP_matrix_chr1.txt.gz"), "wb", mtime=0) as fh:
        fh.write((header + "\n" + body).encode())
    write(os.path.join(DATA, "MGPv8_SNP_matrix_chr1.txt"), header + "\n" + body)
    _ = records7

    bowtie_pg = "@PG\tID:bowtie2\tPN:bowtie2\tVN:2.5.4\tCL:\"bowtie2 -x MGP.N-masked -U reads.fq.gz\""
    rng.seed(1)
    write(os.path.join(DATA, "pure_CAST.sam"), sam(bowtie_pg, simulate(rows, ["CAST_EiJ"], 400, "cast")))
    rng.seed(2)
    hybrid = simulate(rows, ["DBA_2J", "PWK_PhJ"], 600, "hyb") + edge_cases(rows)
    write(os.path.join(DATA, "hybrid_DBA_PWK.sam"), sam(bowtie_pg, hybrid))
    # BAM input, written by samtools so it is the BAM real pipelines hand over.
    subprocess.run(["samtools", "view", "-b", "--no-PG", "-o", os.path.join(DATA, "hybrid_DBA_PWK.bam"),
                    os.path.join(DATA, "hybrid_DBA_PWK.sam")], check=True)

    se_pg = "@PG\tID:Bismark\tVN:v0.24.2\tCL:\"bismark --genome MGP_strains_N-masked reads.fq.gz\""
    pe_pg = "@PG\tID:Bismark\tVN:v0.24.2\tCL:\"bismark --genome MGP_strains_N-masked -1 r1.fq.gz -2 r2.fq.gz\""
    rng.seed(3)
    write(os.path.join(DATA, "bs_WSB_SE.sam"), sam(se_pg, bisulfite_reads(rows, "WSB_EiJ", 300, "bsse", paired=False)))
    rng.seed(4)
    pe = bisulfite_reads(rows, "A_J", 200, "bspe", paired=True)
    write(os.path.join(DATA, "bs_AJ_PE.sam"), sam(pe_pg, pe))
    # The same pairs with the mates of each pair apart, as after a positional sort.
    shuffled = pe[0::2] + pe[1::2]
    write(os.path.join(DATA, "bs_AJ_PE_possorted.sam"), sam(pe_pg, shuffled))


if __name__ == "__main__":
    main()
