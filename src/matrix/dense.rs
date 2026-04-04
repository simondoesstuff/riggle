use bitvec::prelude::*;

/// Dense row-major matrix for accumulating intersection counts
/// Layout: row i = query i, column j = database source j
#[derive(Debug)]
pub struct DenseMatrix {
    data: Vec<u32>,
    num_rows: usize,
    num_cols: usize,
}

impl DenseMatrix {
    /// Allocate a new dense matrix initialized to zero
    pub fn new(num_rows: usize, num_cols: usize) -> Self {
        Self {
            data: vec![0; num_rows * num_cols],
            num_rows,
            num_cols,
        }
    }

    /// Get dimensions
    #[inline]
    pub fn num_rows(&self) -> usize {
        self.num_rows
    }

    #[inline]
    pub fn num_cols(&self) -> usize {
        self.num_cols
    }

    /// Get value at (row, col)
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> u32 {
        debug_assert!(row < self.num_rows && col < self.num_cols);
        self.data[row * self.num_cols + col]
    }

    /// Set value at (row, col)
    #[inline]
    pub fn set(&mut self, row: usize, col: usize, value: u32) {
        debug_assert!(row < self.num_rows && col < self.num_cols);
        self.data[row * self.num_cols + col] = value;
    }

    /// Increment value at (row, col) by delta
    #[inline]
    pub fn add(&mut self, row: usize, col: usize, delta: u32) {
        debug_assert!(row < self.num_rows && col < self.num_cols);
        self.data[row * self.num_cols + col] += delta;
    }

    /// Get a row as a slice
    #[inline]
    pub fn row(&self, row: usize) -> &[u32] {
        let start = row * self.num_cols;
        &self.data[start..start + self.num_cols]
    }

    /// Get a mutable row as a slice
    #[inline]
    pub fn row_mut(&mut self, row: usize) -> &mut [u32] {
        let start = row * self.num_cols;
        &mut self.data[start..start + self.num_cols]
    }

    /// Get raw data slice
    pub fn as_slice(&self) -> &[u32] {
        &self.data
    }

    /// Resize and zero the matrix (only reallocates if larger)
    pub fn resize_and_zero(&mut self, num_rows: usize, num_cols: usize) {
        let new_size = num_rows * num_cols;
        if new_size > self.data.len() {
            self.data.resize(new_size, 0);
        }
        // Zero the used portion
        self.data[..new_size].fill(0);
        self.num_rows = num_rows;
        self.num_cols = num_cols;
    }
}

/// Bitwise mask tracking non-zero columns per row
/// Used to efficiently skip zero regions during sparse condensation
#[derive(Debug)]
pub struct BitwiseMask {
    data: BitVec<u64, Lsb0>,
    num_rows: usize,
    num_cols: usize,
}

impl BitwiseMask {
    /// Allocate a new mask initialized to all false (no flags)
    pub fn new(num_rows: usize, num_cols: usize) -> Self {
        Self {
            data: bitvec![u64, Lsb0; 0; num_rows * num_cols],
            num_rows,
            num_cols,
        }
    }

    /// Get dimensions
    #[inline]
    pub fn num_rows(&self) -> usize {
        self.num_rows
    }

    #[inline]
    pub fn num_cols(&self) -> usize {
        self.num_cols
    }

    /// Flag a position as containing a non-zero value
    #[inline]
    pub fn flag(&mut self, row: usize, col: usize) {
        debug_assert!(row < self.num_rows && col < self.num_cols);
        self.data.set(row * self.num_cols + col, true);
    }

    /// Check if a position is flagged
    #[inline]
    pub fn is_flagged(&self, row: usize, col: usize) -> bool {
        debug_assert!(row < self.num_rows && col < self.num_cols);
        self.data[row * self.num_cols + col]
    }

    /// Clear all flags (reset to zero)
    pub fn clear_all(&mut self) {
        self.data.fill(false);
    }

    /// Iterate over flagged positions in a row
    pub fn flagged_cols(&self, row: usize) -> impl Iterator<Item = usize> + '_ {
        let start = row * self.num_cols;
        let end = start + self.num_cols;
        self.data[start..end].iter_ones()
    }

    /// Count flagged positions in a row
    pub fn count_flagged(&self, row: usize) -> usize {
        let start = row * self.num_cols;
        let end = start + self.num_cols;
        self.data[start..end].count_ones()
    }

    /// Resize and clear the mask (only reallocates if larger)
    pub fn resize_and_clear(&mut self, num_rows: usize, num_cols: usize) {
        let new_size = num_rows * num_cols;
        if new_size > self.data.len() {
            self.data.resize(new_size, false);
        }
        // Clear the used portion
        self.data[..new_size].fill(false);
        self.num_rows = num_rows;
        self.num_cols = num_cols;
    }
}

/// Allocate a paired dense matrix and bitmask for accumulation
pub fn allocate_dense_accumulator(
    num_queries: usize,
    num_db_sids: usize,
) -> (DenseMatrix, BitwiseMask) {
    (
        DenseMatrix::new(num_queries, num_db_sids),
        BitwiseMask::new(num_queries, num_db_sids),
    )
}

/// Zero only the flagged regions in the dense matrix (for reuse)
pub fn zero_flagged_regions(dense: &mut DenseMatrix, mask: &BitwiseMask) {
    for row in 0..dense.num_rows() {
        for col in mask.flagged_cols(row) {
            dense.set(row, col, 0);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dense_matrix_basic() {
        let mut mat = DenseMatrix::new(3, 4);
        assert_eq!(mat.num_rows(), 3);
        assert_eq!(mat.num_cols(), 4);
        assert_eq!(mat.get(0, 0), 0);

        mat.set(1, 2, 42);
        assert_eq!(mat.get(1, 2), 42);

        mat.add(1, 2, 8);
        assert_eq!(mat.get(1, 2), 50);
    }

    #[test]
    fn test_dense_matrix_row_access() {
        let mut mat = DenseMatrix::new(2, 3);
        mat.set(0, 0, 1);
        mat.set(0, 1, 2);
        mat.set(0, 2, 3);

        assert_eq!(mat.row(0), &[1, 2, 3]);
        assert_eq!(mat.row(1), &[0, 0, 0]);

        mat.row_mut(1)[0] = 10;
        assert_eq!(mat.get(1, 0), 10);
    }

    #[test]
    fn test_bitmask_basic() {
        let mut mask = BitwiseMask::new(3, 4);
        assert!(!mask.is_flagged(0, 0));

        mask.flag(1, 2);
        assert!(mask.is_flagged(1, 2));
        assert!(!mask.is_flagged(1, 3));
    }

    #[test]
    fn test_bitmask_flagged_cols() {
        let mut mask = BitwiseMask::new(2, 5);
        mask.flag(0, 1);
        mask.flag(0, 3);
        mask.flag(0, 4);

        let flagged: Vec<_> = mask.flagged_cols(0).collect();
        assert_eq!(flagged, vec![1, 3, 4]);
        assert_eq!(mask.count_flagged(0), 3);

        let row1_flagged: Vec<_> = mask.flagged_cols(1).collect();
        assert!(row1_flagged.is_empty());
    }

    #[test]
    fn test_bitmask_clear() {
        let mut mask = BitwiseMask::new(2, 3);
        mask.flag(0, 0);
        mask.flag(1, 2);

        mask.clear_all();
        assert!(!mask.is_flagged(0, 0));
        assert!(!mask.is_flagged(1, 2));
    }

    #[test]
    fn test_zero_flagged_regions() {
        let mut mat = DenseMatrix::new(2, 4);
        let mut mask = BitwiseMask::new(2, 4);

        // Set some values and flag them
        mat.set(0, 1, 10);
        mat.set(0, 3, 20);
        mat.set(1, 0, 30);
        mask.flag(0, 1);
        mask.flag(0, 3);
        mask.flag(1, 0);

        // Also set an unflagged value
        mat.set(1, 2, 100);

        // Zero flagged regions
        zero_flagged_regions(&mut mat, &mask);

        // Flagged regions should be zero
        assert_eq!(mat.get(0, 1), 0);
        assert_eq!(mat.get(0, 3), 0);
        assert_eq!(mat.get(1, 0), 0);

        // Unflagged region should remain
        assert_eq!(mat.get(1, 2), 100);
    }

    #[test]
    fn test_allocate_dense_accumulator() {
        let (mat, mask) = allocate_dense_accumulator(10, 20);
        assert_eq!(mat.num_rows(), 10);
        assert_eq!(mat.num_cols(), 20);
        assert_eq!(mask.num_rows(), 10);
        assert_eq!(mask.num_cols(), 20);
    }
}
