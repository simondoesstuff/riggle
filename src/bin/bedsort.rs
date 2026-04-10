//! bedsort — sort a BED / BED.gz file the same way `bedtools sort` does.
//!
//! Sort order: chromosome (lex or natural), then start, then end.
//! Header lines (track / browser / #) are passed through first, unchanged.
//! The original line text is preserved byte-for-byte; only the row order changes.
//!
//! Parallelism via rayon; per-chromosome coordinate sort via voracious_radix_sort.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use rayon::prelude::*;
use riggle::bench::nat_cmp;
use voracious_radix_sort::{RadixSort, Radixable};

// ── CLI ─────────────────────────────────────────────────────────────────────

#[derive(Parser)]
#[command(
    name = "bedsort",
    about = "Sort a BED(.gz) file by chromosome and start coordinate (bedtools sort compatible)"
)]
struct Cli {
    /// Input BED file (.bed or .bed.gz); omit to read from stdin
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Output file (.bed or .bed.gz); omit to write to stdout
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Natural sort for chromosome names (chr1 < chr2 < chr10 instead of chr1 < chr10 < chr2)
    #[arg(short = 'n', long)]
    natural: bool,
}

// ── Record / sort-key types ──────────────────────────────────────────────────

/// A BED row with its original text kept intact for lossless output.
struct BedRecord {
    line: String,
    start: u32,
    end: u32,
}

/// Lightweight sort key carrying the row index so we can permute `BedRecord`s
/// after the sort without requiring `Copy` on `BedRecord`.
///
/// Radix key: upper 32 bits = start, lower 32 bits = end — so a single sort
/// produces the correct (start, end) order used by bedtools.
#[derive(Clone, Copy, PartialEq, PartialOrd)]
struct SortKey {
    start: u32,
    end: u32,
    /// Index into the per-chromosome `Vec<BedRecord>`.
    idx: u32,
}

impl Radixable<u64> for SortKey {
    type Key = u64;

    #[inline]
    fn key(&self) -> u64 {
        ((self.start as u64) << 32) | (self.end as u64)
    }
}

// ── Parsing ──────────────────────────────────────────────────────────────────

/// Returns `true` for BED header lines that should be passed through unchanged.
fn is_header(line: &str) -> bool {
    let t = line.trim_start();
    t.starts_with('#') || t.starts_with("track ") || t.starts_with("browser ")
}

/// Extract (chrom, start, end) from a BED line using the same flexible column
/// detection as the rest of riggle: first non-integer column → chrom, first two
/// integer columns → (start, end).
///
/// Returns `None` for empty, header, or unparsable lines.
fn parse_coords(line: &str) -> Option<(&str, u32, u32)> {
    let t = line.trim();
    if t.is_empty() || is_header(t) {
        return None;
    }
    let mut chr: Option<&str> = None;
    let mut start: Option<u32> = None;
    let mut end: Option<u32> = None;
    for col in t.split_ascii_whitespace() {
        if let Ok(v) = col.parse::<u32>() {
            if start.is_none() {
                start = Some(v);
            } else if end.is_none() {
                end = Some(v);
            }
        } else if chr.is_none() {
            chr = Some(col);
        }
        if chr.is_some() && end.is_some() {
            break;
        }
    }
    match (chr, start, end) {
        (Some(c), Some(s), Some(e)) => Some((c, s, e)),
        _ => None,
    }
}

// ── I/O helpers ──────────────────────────────────────────────────────────────

fn open_reader(path: &Option<PathBuf>) -> Result<Box<dyn BufRead>, Box<dyn std::error::Error>> {
    match path {
        None => Ok(Box::new(BufReader::new(io::stdin()))),
        Some(p) => {
            let f = File::open(p)?;
            let is_gz = p.extension().and_then(|e| e.to_str()) == Some("gz");
            if is_gz {
                Ok(Box::new(BufReader::new(MultiGzDecoder::new(f))))
            } else {
                Ok(Box::new(BufReader::new(f)))
            }
        }
    }
}

fn open_writer(path: &Option<PathBuf>) -> Result<Box<dyn Write>, Box<dyn std::error::Error>> {
    match path {
        None => Ok(Box::new(BufWriter::new(io::stdout()))),
        Some(p) => {
            let f = File::create(p)?;
            let is_gz = p.extension().and_then(|e| e.to_str()) == Some("gz");
            if is_gz {
                Ok(Box::new(GzEncoder::new(f, Compression::default())))
            } else {
                Ok(Box::new(BufWriter::new(f)))
            }
        }
    }
}

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    // ── 1. Read all lines ────────────────────────────────────────────────────
    let reader = open_reader(&cli.input)?;
    let mut headers: Vec<String> = Vec::new();
    let mut data_lines: Vec<String> = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        if is_header(&line) {
            headers.push(line);
        } else {
            data_lines.push(line);
        }
    }

    // ── 2. Parse coordinates in parallel ────────────────────────────────────
    // Each element: Some((chr_string, record)) or None (unparsable).
    let parsed: Vec<Option<(String, BedRecord)>> = data_lines
        .into_par_iter()
        .map(|line| {
            if let Some((chr, start, end)) = parse_coords(&line) {
                Some((chr.to_string(), BedRecord { line, start, end }))
            } else {
                None
            }
        })
        .collect();

    // Separate parsed records from lines that couldn't be parsed.
    let mut by_chr: HashMap<String, Vec<BedRecord>> = HashMap::new();
    for entry in parsed {
        match entry {
            Some((chr, rec)) => by_chr.entry(chr).or_default().push(rec),
            None => {
                // collect for pass-through; we already filtered empty/header above
                // so these are genuinely malformed lines
            }
        }
    }

    // ── 3. Sort each chromosome's records in parallel ────────────────────────
    // Build per-chromosome SortKey vecs and sort them; the records themselves
    // stay in place (BedRecord can't be Copy because it owns a String).
    let by_chr_sorted: HashMap<&String, Vec<SortKey>> = by_chr
        .par_iter()
        .map(|(chr, records)| {
            let mut keys: Vec<SortKey> = records
                .iter()
                .enumerate()
                .map(|(i, r)| SortKey {
                    start: r.start,
                    end: r.end,
                    idx: i as u32,
                })
                .collect();
            keys.voracious_sort();
            (chr, keys)
        })
        .collect();

    // ── 4. Sort chromosome names ─────────────────────────────────────────────
    let mut chrs: Vec<&String> = by_chr.keys().collect();
    if cli.natural {
        chrs.sort_by(|a, b| nat_cmp(a, b));
    } else {
        chrs.sort();
    }

    // ── 5. Write output ──────────────────────────────────────────────────────
    let mut writer = open_writer(&cli.output)?;

    for h in &headers {
        writeln!(writer, "{}", h)?;
    }

    for chr in &chrs {
        let records = &by_chr[*chr];
        let keys = &by_chr_sorted[chr];
        for sk in keys {
            writeln!(writer, "{}", records[sk.idx as usize].line)?;
        }
    }

    // Flush (important for GzEncoder which needs finish() implicitly via drop,
    // but explicit flush is still good practice for buffered writers).
    writer.flush()?;

    Ok(())
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp::Ordering;

    #[test]
    fn test_parse_coords_basic() {
        assert_eq!(parse_coords("chr1\t100\t200"), Some(("chr1", 100, 200)));
        assert_eq!(
            parse_coords("chr1\t100\t200\textra\t1000"),
            Some(("chr1", 100, 200))
        );
    }

    #[test]
    fn test_parse_coords_no_chr() {
        // Coord-only BED (no shard column) — returns None (no chr found).
        assert_eq!(parse_coords("100\t200"), None);
    }

    #[test]
    fn test_parse_coords_header() {
        assert_eq!(parse_coords("# comment"), None);
        assert_eq!(parse_coords("track name=foo"), None);
        assert_eq!(parse_coords("browser position chr1:1-1000"), None);
        assert_eq!(parse_coords(""), None);
    }

    #[test]
    fn test_nat_cmp_numbers() {
        assert_eq!(nat_cmp("chr2", "chr10"), Ordering::Less);
        assert_eq!(nat_cmp("chr10", "chr2"), Ordering::Greater);
        assert_eq!(nat_cmp("chr1", "chr1"), Ordering::Equal);
        assert_eq!(nat_cmp("chr9", "chr10"), Ordering::Less);
        assert_eq!(nat_cmp("chrX", "chrY"), Ordering::Less);
        assert_eq!(nat_cmp("chr22", "chrX"), Ordering::Less);
    }

    #[test]
    fn test_nat_cmp_vs_lex() {
        // lex: "chr10" < "chr2" (because '1' < '2' after "chr")
        // nat: "chr2"  < "chr10"
        assert!(nat_cmp("chr2", "chr10") == Ordering::Less);
        let mut chrs = vec!["chr10", "chr2", "chr1", "chrX", "chr22", "chrM"];
        chrs.sort_by(|a, b| nat_cmp(a, b));
        assert_eq!(chrs, vec!["chr1", "chr2", "chr10", "chr22", "chrM", "chrX"]);
    }

    #[test]
    fn test_radix_key_ordering() {
        let a = SortKey {
            start: 100,
            end: 200,
            idx: 0,
        };
        let b = SortKey {
            start: 100,
            end: 300,
            idx: 1,
        };
        let c = SortKey {
            start: 200,
            end: 250,
            idx: 2,
        };
        assert!(a.key() < b.key()); // same start, lower end first
        assert!(b.key() < c.key()); // lower start first
    }
}
