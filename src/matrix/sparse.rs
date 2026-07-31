use sprs::CsMat;

use super::dense::DenseMatrix;

/// Sparse matrix type alias using CSR format.
pub type SparseMatrix = CsMat<u32>;

/// Convert a dense matrix to CSR sparse format.
///
/// Used once at the end of the query pipeline after all chunk matrices have
/// been folded together via `add_dense`.
pub fn condense_to_sparse_no_mask(dense: &DenseMatrix, num_rows: usize, num_cols: usize) -> SparseMatrix {
    let mut indptr = Vec::with_capacity(num_rows + 1);
    let mut indices = Vec::new();
    let mut data = Vec::new();

    indptr.push(0);
    for row in 0..num_rows {
        let row_slice = dense.row(row);
        for (col, &val) in row_slice.iter().enumerate() {
            if val > 0 {
                indices.push(col);
                data.push(val);
            }
        }
        indptr.push(indices.len());
    }

    CsMat::new((num_rows, num_cols), indptr, indices, data)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_condense_to_sparse_no_mask() {
        let mut dense = DenseMatrix::new(3, 4);
        dense.set(0, 1, 5);
        dense.set(1, 0, 10);
        dense.set(1, 3, 15);
        dense.set(2, 2, 20);

        let sparse = condense_to_sparse_no_mask(&dense, 3, 4);

        assert_eq!(sparse.rows(), 3);
        assert_eq!(sparse.cols(), 4);
        assert_eq!(sparse.nnz(), 4);

        assert_eq!(sparse.get(0, 1), Some(&5));
        assert_eq!(sparse.get(1, 0), Some(&10));
        assert_eq!(sparse.get(1, 3), Some(&15));
        assert_eq!(sparse.get(2, 2), Some(&20));
        assert_eq!(sparse.get(0, 0), None);
    }

    #[test]
    fn test_condense_skips_zero_values() {
        let mut dense = DenseMatrix::new(2, 2);
        dense.set(0, 1, 10);

        let sparse = condense_to_sparse_no_mask(&dense, 2, 2);
        assert_eq!(sparse.nnz(), 1);
        assert_eq!(sparse.get(0, 0), None);
        assert_eq!(sparse.get(0, 1), Some(&10));
    }
}
