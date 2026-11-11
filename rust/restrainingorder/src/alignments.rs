//! Reading alignments from SAM or BAM in-process, where the Perl pipes `samtools view -h`.
//!
//! Adapted from SNPsplit's port (`rust/snpsplit/src/io/reader.rs` there).

use std::fs::File;
use std::io::{BufReader, Read};
use std::num::NonZero;
use std::path::Path;

use anyhow::{Context, Result};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;

enum Inner {
    Sam(noodles_sam::io::Reader<BufReader<File>>),
    Bam(noodles_bam::io::Reader<noodles_bgzf::io::MultithreadedReader<File>>),
}

/// Streams owned alignment records out of a SAM or BAM file.
pub struct Alignments {
    inner: Inner,
    header: Header,
}

impl Alignments {
    /// Open `path`, deciding SAM or BAM from its first bytes rather than its name.
    pub fn open(path: &Path) -> Result<Self> {
        let opening = || format!("Failed to read from BAM file {}", path.display());
        let mut magic = [0u8; 2];
        let is_bgzf = File::open(path)
            .with_context(opening)?
            .read_exact(&mut magic)
            .is_ok()
            && magic == [0x1f, 0x8b];

        let file = File::open(path).with_context(opening)?;
        let (inner, header) = if is_bgzf {
            // samtools view decompresses on one thread; a few workers are free speed.
            let workers = std::thread::available_parallelism()
                .map_or(1, |n| n.get().min(4))
                .try_into()
                .unwrap_or(NonZero::<usize>::MIN);
            let decoder = noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, file);
            let mut reader = noodles_bam::io::Reader::from(decoder);
            let header = reader.read_header().with_context(opening)?;
            (Inner::Bam(reader), header)
        } else {
            let mut reader = noodles_sam::io::Reader::new(BufReader::new(file));
            let header = reader.read_header().with_context(opening)?;
            (Inner::Sam(reader), header)
        };
        Ok(Self { inner, header })
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Read the next record into `record`. Returns false at the end of the file.
    pub fn read(&mut self, record: &mut RecordBuf) -> Result<bool> {
        let n = match &mut self.inner {
            Inner::Sam(r) => r.read_record_buf(&self.header, record)?,
            Inner::Bam(r) => r.read_record_buf(&self.header, record)?,
        };
        Ok(n > 0)
    }
}

/// The header as SAM text, one line per record, as `samtools view -h` prints it.
pub fn header_lines(header: &Header) -> Result<Vec<String>> {
    let mut writer = noodles_sam::io::Writer::new(Vec::new());
    writer.write_header(header)?;
    let text = String::from_utf8_lossy(writer.get_ref()).into_owned();
    Ok(text.lines().map(str::to_string).collect())
}
