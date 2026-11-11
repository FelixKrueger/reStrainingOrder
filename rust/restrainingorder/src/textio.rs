//! Line-oriented text I/O shared by the tools.

use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;

/// Open `path` for reading lines, decompressing it if `gzipped`.
///
/// The Perl pipes gzipped input through `gunzip -c`, which reads every member of a
/// multi-member file; bgzip output (the MGP VCFs) is one, so this must be the multi-member
/// decoder, not `GzDecoder`, which stops after the first block.
pub fn open(path: &str, gzipped: bool) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("Failed to open file '{path}'"))?;
    Ok(if gzipped {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    })
}

/// Read the next line into `buf` without its line ending. Returns false at end of input.
pub fn next_line(reader: &mut dyn BufRead, buf: &mut Vec<u8>) -> Result<bool> {
    buf.clear();
    if reader.read_until(b'\n', buf)? == 0 {
        return Ok(false);
    }
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    Ok(true)
}

/// `split /\t/` as Perl does it: trailing empty fields are dropped.
pub fn perl_split(line: &[u8], sep: u8) -> Vec<&[u8]> {
    let mut fields: Vec<&[u8]> = line.split(|&b| b == sep).collect();
    while fields.last().is_some_and(|f| f.is_empty()) {
        fields.pop();
    }
    fields
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perl_split_drops_trailing_empty_fields_only() {
        assert_eq!(perl_split(b"a\t\tb\t\t", b'\t'), vec![&b"a"[..], b"", b"b"]);
        assert!(perl_split(b"", b'\t').is_empty());
    }
}
