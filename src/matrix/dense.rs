/// Dense row-major matrix for accumulating intersection counts.
/// Layout: row i = query i, column j = database source j.
#[derive(Debug, Clone)]
pub struct DenseMatrix {
    data: Vec<u32>,
    num_rows: usize,
    num_cols: usize,
}

impl DenseMatrix {
    pub fn new(num_rows: usize, num_cols: usize) -> Self {
        Self {
            data: vec![0; num_rows * num_cols],
            num_rows,
            num_cols,
        }
    }

    #[inline]
    pub fn num_rows(&self) -> usize {
        self.num_rows
    }

    #[inline]
    pub fn num_cols(&self) -> usize {
        self.num_cols
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> u32 {
        debug_assert!(row < self.num_rows && col < self.num_cols);
        self.data[row * self.num_cols + col]
    }

    #[inline]
    pub fn set(&mut self, row: usize, col: usize, value: u32) {
        debug_assert!(row < self.num_rows && col < self.num_cols);
        self.data[row * self.num_cols + col] = value;
    }

    #[inline]
    pub fn add(&mut self, row: usize, col: usize, delta: u32) {
        debug_assert!(row < self.num_rows && col < self.num_cols);
        self.data[row * self.num_cols + col] += delta;
    }

    #[inline]
    pub fn row(&self, row: usize) -> &[u32] {
        let start = row * self.num_cols;
        &self.data[start..start + self.num_cols]
    }

    /// SIMD-friendly element-wise addition; the compiler auto-vectorizes this loop.
    pub fn add_dense(&mut self, other: &DenseMatrix) {
        debug_assert_eq!(self.num_rows, other.num_rows);
        debug_assert_eq!(self.num_cols, other.num_cols);
        for (a, b) in self.data.iter_mut().zip(other.data.iter()) {
            *a += b;
        }
    }

    /// Resize the matrix to new dimensions, zeroing all elements.
    /// Only reallocates if the new size exceeds the current allocation.
    pub fn resize_and_zero(&mut self, num_rows: usize, num_cols: usize) {
        let new_size = num_rows * num_cols;
        if new_size > self.data.len() {
            self.data.resize(new_size, 0);
        }
        self.data[..new_size].fill(0);
        self.num_rows = num_rows;
        self.num_cols = num_cols;
    }
}

/// Row-major f32 matrix for accumulating exact bin-overlap totals.
/// Mirrors `DenseMatrix` layout but stores base-pair overlaps in 100 bp bins.
#[derive(Debug, Clone)]
pub struct OverlapMatrix {
    data: Vec<f32>,
    num_rows: usize,
    num_cols: usize,
}

impl OverlapMatrix {
    pub fn new(num_rows: usize, num_cols: usize) -> Self {
        Self {
            data: vec![0.0; num_rows * num_cols],
            num_rows,
            num_cols,
        }
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> f32 {
        debug_assert!(row < self.num_rows && col < self.num_cols);
        self.data[row * self.num_cols + col]
    }

    #[inline]
    pub fn add(&mut self, row: usize, col: usize, delta: f32) {
        debug_assert!(row < self.num_rows && col < self.num_cols);
        self.data[row * self.num_cols + col] += delta;
    }

    /// Element-wise f32 addition; auto-vectorizes to SIMD on supported targets.
    pub fn add_dense(&mut self, other: &OverlapMatrix) {
        debug_assert_eq!(self.num_rows, other.num_rows);
        debug_assert_eq!(self.num_cols, other.num_cols);
        let n = self.num_rows * self.num_cols;
        for (a, b) in self.data[..n].iter_mut().zip(other.data[..n].iter()) {
            *a += b;
        }
    }

    pub fn resize_and_zero(&mut self, num_rows: usize, num_cols: usize) {
        let new_size = num_rows * num_cols;
        if new_size > self.data.len() {
            self.data.resize(new_size, 0.0);
        }
        self.data[..new_size].fill(0.0);
        self.num_rows = num_rows;
        self.num_cols = num_cols;
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
    fn test_dense_matrix_row() {
        let mut mat = DenseMatrix::new(2, 3);
        mat.set(0, 0, 1);
        mat.set(0, 1, 2);
        mat.set(0, 2, 3);

        assert_eq!(mat.row(0), &[1, 2, 3]);
        assert_eq!(mat.row(1), &[0, 0, 0]);
    }

    #[test]
    fn test_dense_matrix_add_dense() {
        let matrices: Vec<DenseMatrix> = (1u32..=5)
            .map(|i| {
                let mut d = DenseMatrix::new(2, 2);
                d.set(0, 0, i);
                d
            })
            .collect();
        let folded = matrices
            .into_iter()
            .fold(DenseMatrix::new(2, 2), |mut acc, m| {
                acc.add_dense(&m);
                acc
            });
        assert_eq!(folded.get(0, 0), 15);
        assert_eq!(folded.get(1, 0), 0);
    }

    #[test]
    fn test_dense_matrix_resize_and_zero() {
        let mut mat = DenseMatrix::new(2, 2);
        mat.set(0, 0, 99);
        mat.resize_and_zero(3, 4);
        assert_eq!(mat.num_rows(), 3);
        assert_eq!(mat.num_cols(), 4);
        for row in 0..3 {
            for col in 0..4 {
                assert_eq!(mat.get(row, col), 0);
            }
        }
    }
}
