use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use thiserror::Error;

use crate::core::TaggedInterval;

/// Errors that can occur during BED parsing
#[derive(Debug, Error)]
pub enum BedParseError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Line {line}: invalid format - expected at least 3 columns")]
    InvalidFormat { line: usize },

    #[error("Line {line}: invalid start coordinate '{value}'")]
    InvalidStart { line: usize, value: String },

    #[error("Line {line}: invalid end coordinate '{value}'")]
    InvalidEnd { line: usize, value: String },

    #[error("Line {line}: start ({start}) must be less than end ({end})")]
    InvalidRange { line: usize, start: u32, end: u32 },
}

/// Parse a BED file into a vector of TaggedIntervals
///
/// Features:
/// - Skips comment lines (starting with #)
/// - Skips empty or whitespace-only lines
/// - Parses only first 3 columns (chrom is ignored, we use sid for chromosome identity)
/// - Supports tab or space delimiters
/// - Returns intervals tagged with the given source ID
pub fn parse_bed_file(path: &Path, sid: u32) -> Result<Vec<TaggedInterval>, BedParseError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();

    for (line_idx, line_result) in reader.lines().enumerate() {
        let line_num = line_idx + 1;
        let line = line_result?;
        let trimmed = line.trim();

        // Skip empty lines and comments
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        // Split by whitespace (handles both tabs and spaces)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() < 3 {
            return Err(BedParseError::InvalidFormat { line: line_num });
        }

        // Parse start coordinate
        let start: u32 = parts[1].parse().map_err(|_| BedParseError::InvalidStart {
            line: line_num,
            value: parts[1].to_string(),
        })?;

        // Parse end coordinate
        let end: u32 = parts[2].parse().map_err(|_| BedParseError::InvalidEnd {
            line: line_num,
            value: parts[2].to_string(),
        })?;

        // Validate range
        if start >= end {
            return Err(BedParseError::InvalidRange {
                line: line_num,
                start,
                end,
            });
        }

        intervals.push(TaggedInterval::new(start, end, sid));
    }

    Ok(intervals)
}

/// Parse BED data from a string (useful for testing)
pub fn parse_bed_string(content: &str, sid: u32) -> Result<Vec<TaggedInterval>, BedParseError> {
    let mut intervals = Vec::new();

    for (line_idx, line) in content.lines().enumerate() {
        let line_num = line_idx + 1;
        let trimmed = line.trim();

        // Skip empty lines and comments
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        // Split by whitespace
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() < 3 {
            return Err(BedParseError::InvalidFormat { line: line_num });
        }

        let start: u32 = parts[1].parse().map_err(|_| BedParseError::InvalidStart {
            line: line_num,
            value: parts[1].to_string(),
        })?;

        let end: u32 = parts[2].parse().map_err(|_| BedParseError::InvalidEnd {
            line: line_num,
            value: parts[2].to_string(),
        })?;

        if start >= end {
            return Err(BedParseError::InvalidRange {
                line: line_num,
                start,
                end,
            });
        }

        intervals.push(TaggedInterval::new(start, end, sid));
    }

    Ok(intervals)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parse_basic_bed() {
        let content = "chr1\t100\t200\nchr1\t300\t400\n";
        let intervals = parse_bed_string(content, 1).unwrap();

        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].iv.start, 100);
        assert_eq!(intervals[0].iv.end, 200);
        assert_eq!(intervals[0].sid, 1);
        assert_eq!(intervals[1].iv.start, 300);
        assert_eq!(intervals[1].iv.end, 400);
    }

    #[test]
    fn test_parse_extended_bed() {
        // BED with extra columns (name, score, strand, etc.)
        let content = "chr1\t100\t200\tpeak1\t500\t+\n";
        let intervals = parse_bed_string(content, 2).unwrap();

        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].iv.start, 100);
        assert_eq!(intervals[0].iv.end, 200);
    }

    #[test]
    fn test_parse_with_comments() {
        let content = "# This is a comment\nchr1\t100\t200\n# Another comment\nchr1\t300\t400\n";
        let intervals = parse_bed_string(content, 1).unwrap();

        assert_eq!(intervals.len(), 2);
    }

    #[test]
    fn test_parse_with_empty_lines() {
        let content = "chr1\t100\t200\n\n\nchr1\t300\t400\n   \nchr1\t500\t600\n";
        let intervals = parse_bed_string(content, 1).unwrap();

        assert_eq!(intervals.len(), 3);
    }

    #[test]
    fn test_parse_space_delimited() {
        let content = "chr1 100 200\nchr1   300   400\n";
        let intervals = parse_bed_string(content, 1).unwrap();

        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].iv.start, 100);
        assert_eq!(intervals[1].iv.start, 300);
    }

    #[test]
    fn test_parse_invalid_format() {
        let content = "chr1\t100\n"; // Only 2 columns
        let result = parse_bed_string(content, 1);

        assert!(matches!(
            result,
            Err(BedParseError::InvalidFormat { line: 1 })
        ));
    }

    #[test]
    fn test_parse_invalid_start() {
        let content = "chr1\tabc\t200\n";
        let result = parse_bed_string(content, 1);

        assert!(matches!(result, Err(BedParseError::InvalidStart { .. })));
    }

    #[test]
    fn test_parse_invalid_end() {
        let content = "chr1\t100\txyz\n";
        let result = parse_bed_string(content, 1);

        assert!(matches!(result, Err(BedParseError::InvalidEnd { .. })));
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
    fn test_parse_negative_values() {
        // u32 parse will fail on negative numbers
        let content = "chr1\t-100\t200\n";
        let result = parse_bed_string(content, 1);

        assert!(matches!(result, Err(BedParseError::InvalidStart { .. })));
    }

    #[test]
    fn test_parse_bed_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "chr1\t100\t200").unwrap();
        writeln!(file, "chr1\t300\t400").unwrap();
        file.flush().unwrap();

        let intervals = parse_bed_file(file.path(), 42).unwrap();
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].sid, 42);
    }
}
