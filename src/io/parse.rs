use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use thiserror::Error;

use crate::core::TaggedInterval;

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
/// Uses borrowed string slice for shard to avoid per-line allocations
struct ParsedLine<'a> {
    /// Shard name (first non-integer column, or DEFAULT_SHARD if none)
    shard: &'a str,
    /// The interval
    interval: TaggedInterval,
}

/// Parse a single BED line with flexible column detection
///
/// Logic:
/// - Scan columns left-to-right
/// - First non-integer column -> shard name
/// - First two integer columns -> (start, end) coordinates
/// - If no non-integer column found -> use DEFAULT_SHARD
///
/// Performance: Zero allocations per line - borrows shard from input,
/// uses scalar variables for coordinates instead of Vec.
///
/// Returns:
/// - None if the line should be skipped (empty or comment)
/// - Some(Ok(parsed)) for valid intervals
/// - Some(Err(error)) for invalid lines
fn parse_bed_line<'a>(
    line_num: usize,
    line: &'a str,
    sid: u32,
) -> Option<Result<ParsedLine<'a>, BedParseError>> {
    let trimmed = line.trim();

    // Skip empty lines and comments
    if trimmed.is_empty() || trimmed.starts_with('#') {
        return None;
    }

    // Scan columns to find shard and coordinates
    // Use scalar variables instead of Vec to avoid allocation
    let mut shard: Option<&'a str> = None;
    let mut start: Option<u32> = None;
    let mut end: Option<u32> = None;
    let mut col_count = 0;

    // Iterate directly instead of collecting to Vec
    for part in trimmed.split_whitespace() {
        col_count += 1;
        if let Ok(val) = part.parse::<u32>() {
            // Integer column - use as coordinate (take first two)
            if start.is_none() {
                start = Some(val);
            } else if end.is_none() {
                end = Some(val);
            }
        } else if shard.is_none() {
            // First non-integer column - borrow as shard name (no allocation!)
            shard = Some(part);
        }
        // Stop once we have both shard and 2 coordinates
        if shard.is_some() && end.is_some() {
            break;
        }
    }

    // Need at least 2 columns to be valid
    if col_count < 2 {
        return Some(Err(BedParseError::InvalidFormat { line: line_num }));
    }

    // Need at least 2 coordinates
    let (start, end) = match (start, end) {
        (Some(s), Some(e)) => (s, e),
        _ => return Some(Err(BedParseError::InvalidFormat { line: line_num })),
    };

    // Validate range
    if start >= end {
        return Some(Err(BedParseError::InvalidRange {
            line: line_num,
            start,
            end,
        }));
    }

    Some(Ok(ParsedLine {
        shard: shard.unwrap_or(DEFAULT_SHARD),
        interval: TaggedInterval::new(start, end, sid),
    }))
}

/// Parse a BED file into intervals grouped by shard
///
/// Features:
/// - Automatically infers and decompresses .gz files based on extension
/// - Skips comment lines (starting with #)
/// - Skips empty or whitespace-only lines
/// - Flexible column detection: first non-integer = shard, first two integers = coords
/// - Supports tab or space delimiters
/// - Returns intervals grouped by shard name, tagged with the given source ID
///
/// Performance: Only allocates a String once per unique shard (~24 for human genome)
/// instead of once per line (millions of allocations avoided).
pub fn parse_bed_file(
    path: &Path,
    sid: u32,
) -> Result<HashMap<String, Vec<TaggedInterval>>, BedParseError> {
    let file = File::open(path)?;

    // Check if the file extension is .gz
    let is_gzipped = path.extension().and_then(OsStr::to_str) == Some("gz");

    // Use Box<dyn BufRead> to dynamically dispatch the reader type
    // This allows us to use the same iteration logic below
    let reader: Box<dyn BufRead> = if is_gzipped {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut shards: HashMap<String, Vec<TaggedInterval>> = HashMap::new();

    for (line_idx, line_result) in reader.lines().enumerate() {
        let line_num = line_idx + 1;
        let line = line_result?;

        if let Some(result) = parse_bed_line(line_num, &line, sid) {
            let parsed = result?;
            // Check using &str first - only allocate String if shard is new
            if let Some(intervals) = shards.get_mut(parsed.shard) {
                intervals.push(parsed.interval);
            } else {
                shards.insert(parsed.shard.to_string(), vec![parsed.interval]);
            }
        }
    }

    Ok(shards)
}

/// Parse BED data from a string (useful for testing)
/// Returns intervals grouped by shard
///
/// Performance: Only allocates a String once per unique shard.
pub fn parse_bed_string(
    content: &str,
    sid: u32,
) -> Result<HashMap<String, Vec<TaggedInterval>>, BedParseError> {
    let mut shards: HashMap<String, Vec<TaggedInterval>> = HashMap::new();

    for (line_idx, line) in content.lines().enumerate() {
        let line_num = line_idx + 1;

        if let Some(result) = parse_bed_line(line_num, line, sid) {
            let parsed = result?;
            // Check using &str first - only allocate String if shard is new
            if let Some(intervals) = shards.get_mut(parsed.shard) {
                intervals.push(parsed.interval);
            } else {
                shards.insert(parsed.shard.to_string(), vec![parsed.interval]);
            }
        }
    }

    Ok(shards)
}

/// Parse BED data from a string, returning a flat vector of intervals (ignores shard)
/// This is a convenience function for tests that don't need shard information
#[cfg(test)]
pub fn parse_bed_string_flat(
    content: &str,
    sid: u32,
) -> Result<Vec<TaggedInterval>, BedParseError> {
    let shards = parse_bed_string(content, sid)?;
    Ok(shards.into_values().flatten().collect())
}

/// Parse a BED file, returning a flat vector of intervals (ignores shard)
/// This is a convenience function for cases that don't need shard information
#[cfg(test)]
pub fn parse_bed_file_flat(path: &Path, sid: u32) -> Result<Vec<TaggedInterval>, BedParseError> {
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

        // All intervals should be in chr1 shard
        assert_eq!(shards.len(), 1);
        assert!(shards.contains_key("chr1"));
        let intervals = &shards["chr1"];
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].iv.start, 100);
        assert_eq!(intervals[0].iv.end, 200);
        assert_eq!(intervals[0].sid, 1);
        assert_eq!(intervals[1].iv.start, 300);
        assert_eq!(intervals[1].iv.end, 400);
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
        // BED with extra columns (name, score, strand, etc.)
        let content = "chr1\t100\t200\tpeak1\t500\t+\n";
        let shards = parse_bed_string(content, 2).unwrap();

        assert_eq!(shards["chr1"].len(), 1);
        assert_eq!(shards["chr1"][0].iv.start, 100);
        assert_eq!(shards["chr1"][0].iv.end, 200);
    }

    #[test]
    fn test_parse_no_shard_column() {
        // When all columns are integers, use DEFAULT_SHARD
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
        assert_eq!(intervals[0].iv.start, 100);
        assert_eq!(intervals[1].iv.start, 300);
    }

    #[test]
    fn test_parse_invalid_format() {
        let content = "100\n"; // Only 1 column
        let result = parse_bed_string(content, 1);

        assert!(matches!(
            result,
            Err(BedParseError::InvalidFormat { line: 1 })
        ));
    }

    #[test]
    fn test_parse_invalid_range() {
        let content = "chr1\t200\t100\n"; // start > end
        let result = parse_bed_string(content, 1);

        assert!(matches!(result, Err(BedParseError::InvalidRange { .. })));
    }

    #[test]
    fn test_parse_equal_start_end() {
        let content = "chr1\t100\t100\n"; // start == end (empty interval)
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
