use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use thiserror::Error;

use crate::core::TaggedInterval;

/// Default shard name when no non-integer column is found
pub const DEFAULT_SHARD: &str = "default";

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
struct ParsedLine {
    /// Shard name (first non-integer column, or DEFAULT_SHARD if none)
    shard: String,
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
/// Returns:
/// - None if the line should be skipped (empty or comment)
/// - Some(Ok(parsed)) for valid intervals
/// - Some(Err(error)) for invalid lines
fn parse_bed_line(line_num: usize, line: &str, sid: u32) -> Option<Result<ParsedLine, BedParseError>> {
    let trimmed = line.trim();

    // Skip empty lines and comments
    if trimmed.is_empty() || trimmed.starts_with('#') {
        return None;
    }

    // Split by whitespace (handles both tabs and spaces)
    let parts: Vec<&str> = trimmed.split_whitespace().collect();

    if parts.len() < 2 {
        return Some(Err(BedParseError::InvalidFormat { line: line_num }));
    }

    // Scan columns to find shard and coordinates
    let mut shard: Option<String> = None;
    let mut coords: Vec<u32> = Vec::with_capacity(2);

    for part in &parts {
        if let Ok(val) = part.parse::<u32>() {
            // Integer column - use as coordinate (take first two)
            if coords.len() < 2 {
                coords.push(val);
            }
        } else if shard.is_none() {
            // First non-integer column - use as shard name
            shard = Some(part.to_string());
        }
        // Stop once we have both shard and 2 coordinates
        if shard.is_some() && coords.len() == 2 {
            break;
        }
    }

    // Need at least 2 coordinates
    if coords.len() < 2 {
        return Some(Err(BedParseError::InvalidFormat { line: line_num }));
    }

    let start = coords[0];
    let end = coords[1];

    // Validate range
    if start >= end {
        return Some(Err(BedParseError::InvalidRange {
            line: line_num,
            start,
            end,
        }));
    }

    let shard_name = shard.unwrap_or_else(|| DEFAULT_SHARD.to_string());

    Some(Ok(ParsedLine {
        shard: shard_name,
        interval: TaggedInterval::new(start, end, sid),
    }))
}

/// Parse a BED file into intervals grouped by shard
///
/// Features:
/// - Skips comment lines (starting with #)
/// - Skips empty or whitespace-only lines
/// - Flexible column detection: first non-integer = shard, first two integers = coords
/// - Supports tab or space delimiters
/// - Returns intervals grouped by shard name, tagged with the given source ID
pub fn parse_bed_file(path: &Path, sid: u32) -> Result<HashMap<String, Vec<TaggedInterval>>, BedParseError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut shards: HashMap<String, Vec<TaggedInterval>> = HashMap::new();

    for (line_idx, line_result) in reader.lines().enumerate() {
        let line_num = line_idx + 1;
        let line = line_result?;

        if let Some(result) = parse_bed_line(line_num, &line, sid) {
            let parsed = result?;
            shards.entry(parsed.shard).or_default().push(parsed.interval);
        }
    }

    Ok(shards)
}

/// Parse BED data from a string (useful for testing)
/// Returns intervals grouped by shard
pub fn parse_bed_string(content: &str, sid: u32) -> Result<HashMap<String, Vec<TaggedInterval>>, BedParseError> {
    let mut shards: HashMap<String, Vec<TaggedInterval>> = HashMap::new();

    for (line_idx, line) in content.lines().enumerate() {
        let line_num = line_idx + 1;

        if let Some(result) = parse_bed_line(line_num, line, sid) {
            let parsed = result?;
            shards.entry(parsed.shard).or_default().push(parsed.interval);
        }
    }

    Ok(shards)
}

/// Parse BED data from a string, returning a flat vector of intervals (ignores shard)
/// This is a convenience function for tests that don't need shard information
#[cfg(test)]
pub fn parse_bed_string_flat(content: &str, sid: u32) -> Result<Vec<TaggedInterval>, BedParseError> {
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
