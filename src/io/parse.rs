use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use thiserror::Error;

use crate::core::Interval;

/// Default shard name when no non-integer column is found
pub const DEFAULT_SHARD: &str = "default";

/// Check if a path is a BED file (.bed or .bed.gz)
pub fn is_bed_file(path: &Path) -> bool {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    name.ends_with(".bed") || name.ends_with(".bed.gz")
}

/// Errors that can occur during BED parsing
#[derive(Debug, Error)]
pub enum BedParseError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Line {line}: invalid format - expected at least 2 integer columns")]
    InvalidFormat { line: usize },

    #[error("Line {line}: invalid start coordinate '{value}'")]
    InvalidStart { line: usize, value: String },

    #[error("Line {line}: invalid end coordinate '{value}'")]
    InvalidEnd { line: usize, value: String },

    #[error("Line {line}: start ({start}) must be less than end ({end})")]
    InvalidRange { line: usize, start: u32, end: u32 },
}

/// Result of parsing a single BED line
struct ParsedLine<'a> {
    shard: &'a str,
    interval: Interval,
}

/// Parse a single BED line with flexible column detection.
///
/// - First non-integer column → shard name
/// - First two integer columns → (start, end) coordinates
/// - If no non-integer column found → use DEFAULT_SHARD
///
/// Zero allocations per line: borrows shard from the input string.
fn parse_bed_line<'a>(
    line_num: usize,
    line: &'a str,
    sid: u32,
) -> Option<Result<ParsedLine<'a>, BedParseError>> {
    let trimmed = line.trim();

    if trimmed.is_empty() || trimmed.starts_with('#') {
        return None;
    }

    let mut shard: Option<&'a str> = None;
    let mut start: Option<u32> = None;
    let mut end: Option<u32> = None;
    let mut col_count = 0;

    for part in trimmed.split_whitespace() {
        col_count += 1;
        if let Ok(val) = part.parse::<u32>() {
            if start.is_none() {
                start = Some(val);
            } else if end.is_none() {
                end = Some(val);
            }
        } else if shard.is_none() {
            shard = Some(part);
        }
        if shard.is_some() && end.is_some() {
            break;
        }
    }

    if col_count < 2 {
        return Some(Err(BedParseError::InvalidFormat { line: line_num }));
    }

    let (start, end) = match (start, end) {
        (Some(s), Some(e)) => (s, e),
        _ => return Some(Err(BedParseError::InvalidFormat { line: line_num })),
    };

    if start >= end {
        return Some(Err(BedParseError::InvalidRange {
            line: line_num,
            start,
            end,
        }));
    }

    Some(Ok(ParsedLine {
        shard: shard.unwrap_or(DEFAULT_SHARD),
        interval: Interval::new(start, end, sid),
    }))
}

/// Parse a BED file into intervals grouped by shard.
///
/// Supports plain `.bed` and gzip-compressed `.bed.gz` files.
/// Skips comment lines (starting with `#`) and empty lines.
/// The `sid` is embedded into every returned interval.
///
/// Performance: only allocates one `String` per unique shard name.
pub fn parse_bed_file(
    path: &Path,
    sid: u32,
) -> Result<HashMap<String, Vec<Interval>>, BedParseError> {
    let file = File::open(path)?;
    let is_gzipped = path.extension().and_then(OsStr::to_str) == Some("gz");

    let reader: Box<dyn BufRead> = if is_gzipped {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut shards: HashMap<String, Vec<Interval>> = HashMap::new();

    for (line_idx, line_result) in reader.lines().enumerate() {
        let line_num = line_idx + 1;
        let line = line_result?;

        if let Some(result) = parse_bed_line(line_num, &line, sid) {
            let parsed = result?;
            if let Some(intervals) = shards.get_mut(parsed.shard) {
                intervals.push(parsed.interval);
            } else {
                shards.insert(parsed.shard.to_string(), vec![parsed.interval]);
            }
        }
    }

    Ok(shards)
}

/// Parse BED data from a string (useful for testing).
pub fn parse_bed_string(
    content: &str,
    sid: u32,
) -> Result<HashMap<String, Vec<Interval>>, BedParseError> {
    let mut shards: HashMap<String, Vec<Interval>> = HashMap::new();

    for (line_idx, line) in content.lines().enumerate() {
        let line_num = line_idx + 1;

        if let Some(result) = parse_bed_line(line_num, line, sid) {
            let parsed = result?;
            if let Some(intervals) = shards.get_mut(parsed.shard) {
                intervals.push(parsed.interval);
            } else {
                shards.insert(parsed.shard.to_string(), vec![parsed.interval]);
            }
        }
    }

    Ok(shards)
}

/// Flatten parsed shards into a single vec (for tests that don't need shard grouping).
#[cfg(test)]
pub fn parse_bed_string_flat(
    content: &str,
    sid: u32,
) -> Result<Vec<Interval>, BedParseError> {
    let shards = parse_bed_string(content, sid)?;
    Ok(shards.into_values().flatten().collect())
}

/// Flatten a parsed BED file into a single vec (for tests).
#[cfg(test)]
pub fn parse_bed_file_flat(path: &Path, sid: u32) -> Result<Vec<Interval>, BedParseError> {
    let shards = parse_bed_file(path, sid)?;
    Ok(shards.into_values().flatten().collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parse_basic_bed() {
        let content = "chr1\t100\t200\nchr1\t300\t400\n";
        let shards = parse_bed_string(content, 1).unwrap();

        assert_eq!(shards.len(), 1);
        assert!(shards.contains_key("chr1"));
        let intervals = &shards["chr1"];
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].start, 100);
        assert_eq!(intervals[0].end, 200);
        assert_eq!(intervals[0].sid, 1);
        assert_eq!(intervals[1].start, 300);
        assert_eq!(intervals[1].end, 400);
    }

    #[test]
    fn test_parse_multi_shard() {
        let content = "chr1\t100\t200\nchr2\t300\t400\nchr1\t500\t600\n";
        let shards = parse_bed_string(content, 1).unwrap();

        assert_eq!(shards.len(), 2);
        assert_eq!(shards["chr1"].len(), 2);
        assert_eq!(shards["chr2"].len(), 1);
    }

    #[test]
    fn test_parse_extended_bed() {
        let content = "chr1\t100\t200\tpeak1\t500\t+\n";
        let shards = parse_bed_string(content, 2).unwrap();

        assert_eq!(shards["chr1"].len(), 1);
        assert_eq!(shards["chr1"][0].start, 100);
        assert_eq!(shards["chr1"][0].end, 200);
    }

    #[test]
    fn test_parse_no_shard_column() {
        let content = "100\t200\n300\t400\n";
        let shards = parse_bed_string(content, 1).unwrap();

        assert_eq!(shards.len(), 1);
        assert!(shards.contains_key(DEFAULT_SHARD));
        assert_eq!(shards[DEFAULT_SHARD].len(), 2);
    }

    #[test]
    fn test_parse_with_comments() {
        let content = "# This is a comment\nchr1\t100\t200\n# Another comment\nchr1\t300\t400\n";
        let intervals = parse_bed_string_flat(content, 1).unwrap();
        assert_eq!(intervals.len(), 2);
    }

    #[test]
    fn test_parse_with_empty_lines() {
        let content = "chr1\t100\t200\n\n\nchr1\t300\t400\n   \nchr1\t500\t600\n";
        let intervals = parse_bed_string_flat(content, 1).unwrap();
        assert_eq!(intervals.len(), 3);
    }

    #[test]
    fn test_parse_space_delimited() {
        let content = "chr1 100 200\nchr1   300   400\n";
        let shards = parse_bed_string(content, 1).unwrap();
        let intervals = &shards["chr1"];
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].start, 100);
        assert_eq!(intervals[1].start, 300);
    }

    #[test]
    fn test_parse_invalid_format() {
        let content = "100\n";
        let result = parse_bed_string(content, 1);
        assert!(matches!(result, Err(BedParseError::InvalidFormat { line: 1 })));
    }

    #[test]
    fn test_parse_invalid_range() {
        let content = "chr1\t200\t100\n";
        let result = parse_bed_string(content, 1);
        assert!(matches!(result, Err(BedParseError::InvalidRange { .. })));
    }

    #[test]
    fn test_parse_equal_start_end() {
        let content = "chr1\t100\t100\n";
        let result = parse_bed_string(content, 1);
        assert!(matches!(result, Err(BedParseError::InvalidRange { .. })));
    }

    #[test]
    fn test_parse_bed_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "chr1\t100\t200").unwrap();
        writeln!(file, "chr2\t300\t400").unwrap();
        file.flush().unwrap();

        let shards = parse_bed_file(file.path(), 42).unwrap();
        assert_eq!(shards.len(), 2);
        assert_eq!(shards["chr1"][0].sid, 42);
        assert_eq!(shards["chr2"][0].sid, 42);
    }

    #[test]
    fn test_parse_bed_file_flat() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "chr1\t100\t200").unwrap();
        writeln!(file, "chr2\t300\t400").unwrap();
        file.flush().unwrap();

        let intervals = parse_bed_file_flat(file.path(), 42).unwrap();
        assert_eq!(intervals.len(), 2);
    }
}
