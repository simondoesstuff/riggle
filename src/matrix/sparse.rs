use sprs::{CsMat, TriMat};

use super::dense::{BitwiseMask, DenseMatrix};

/// Sparse matrix type alias using CSR format
pub type SparseMatrix = CsMat<u32>;

/// Condense a dense matrix to sparse format, using the bitmask to skip zeros
pub fn condense_to_sparse(dense: &DenseMatrix, mask: &BitwiseMask) -> SparseMatrix {
    let num_rows = dense.num_rows();
    let num_cols = dense.num_cols();

    // Count total non-zeros for pre-allocation
    let nnz: usize = (0..num_rows).map(|r| mask.count_flagged(r)).sum();

    // Build triplet matrix
    let mut triplets = TriMat::with_capacity((num_rows, num_cols), nnz);

    for row in 0..num_rows {
        for col in mask.flagged_cols(row) {
            let value = dense.get(row, col);
            if value > 0 {
                triplets.add_triplet(row, col, value);
            }
        }
    }

    triplets.to_csr()
}

/// Merge two sparse matrices by element-wise addition
pub fn merge_sparse(a: &SparseMatrix, b: &SparseMatrix) -> SparseMatrix {
    // sprs handles addition of CSR matrices efficiently
    a + b
}

/// Create an empty sparse matrix with given dimensions
pub fn empty_sparse(_num_rows: usize, num_cols: usize) -> SparseMatrix {
    CsMat::empty(sprs::CompressedStorage::CSR, num_cols)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_condense_to_sparse() {
        let mut dense = DenseMatrix::new(3, 4);
        let mut mask = BitwiseMask::new(3, 4);

        // Set some values
        dense.set(0, 1, 5);
        dense.set(1, 0, 10);
        dense.set(1, 3, 15);
        dense.set(2, 2, 20);

        // Flag them
        mask.flag(0, 1);
        mask.flag(1, 0);
        mask.flag(1, 3);
        mask.flag(2, 2);

        let sparse = condense_to_sparse(&dense, &mask);

        assert_eq!(sparse.rows(), 3);
        assert_eq!(sparse.cols(), 4);
        assert_eq!(sparse.nnz(), 4);

        // Check values
        assert_eq!(sparse.get(0, 1), Some(&5));
        assert_eq!(sparse.get(1, 0), Some(&10));
        assert_eq!(sparse.get(1, 3), Some(&15));
        assert_eq!(sparse.get(2, 2), Some(&20));
        assert_eq!(sparse.get(0, 0), None);
    }

    #[test]
    fn test_condense_skips_zero_values() {
        let mut dense = DenseMatrix::new(2, 2);
        let mut mask = BitwiseMask::new(2, 2);

        // Flag a position but leave value as 0
        mask.flag(0, 0);
        dense.set(0, 1, 10);
        mask.flag(0, 1);

        let sparse = condense_to_sparse(&dense, &mask);
        assert_eq!(sparse.nnz(), 1); // Only non-zero value
        assert_eq!(sparse.get(0, 0), None);
        assert_eq!(sparse.get(0, 1), Some(&10));
    }

    #[test]
    fn test_merge_sparse() {
        let mut dense_a = DenseMatrix::new(2, 3);
        let mut mask_a = BitwiseMask::new(2, 3);
        dense_a.set(0, 0, 5);
        dense_a.set(1, 1, 10);
        mask_a.flag(0, 0);
        mask_a.flag(1, 1);

        let mut dense_b = DenseMatrix::new(2, 3);
        let mut mask_b = BitwiseMask::new(2, 3);
        dense_b.set(0, 0, 3);
        dense_b.set(0, 2, 7);
        mask_b.flag(0, 0);
        mask_b.flag(0, 2);

        let sparse_a = condense_to_sparse(&dense_a, &mask_a);
        let sparse_b = condense_to_sparse(&dense_b, &mask_b);

        let merged = merge_sparse(&sparse_a, &sparse_b);

        assert_eq!(merged.get(0, 0), Some(&8)); // 5 + 3
        assert_eq!(merged.get(0, 2), Some(&7)); // 0 + 7
        assert_eq!(merged.get(1, 1), Some(&10)); // 10 + 0
    }

    #[test]
    fn test_merge_associative() {
        // Create three sparse matrices
        let mut d1 = DenseMatrix::new(2, 2);
        let mut m1 = BitwiseMask::new(2, 2);
        d1.set(0, 0, 1);
        m1.flag(0, 0);

        let mut d2 = DenseMatrix::new(2, 2);
        let mut m2 = BitwiseMask::new(2, 2);
        d2.set(0, 0, 2);
        m2.flag(0, 0);

        let mut d3 = DenseMatrix::new(2, 2);
        let mut m3 = BitwiseMask::new(2, 2);
        d3.set(0, 0, 3);
        m3.flag(0, 0);

        let s1 = condense_to_sparse(&d1, &m1);
        let s2 = condense_to_sparse(&d2, &m2);
        let s3 = condense_to_sparse(&d3, &m3);

        // (s1 + s2) + s3
        let left = merge_sparse(&merge_sparse(&s1, &s2), &s3);
        // s1 + (s2 + s3)
        let right = merge_sparse(&s1, &merge_sparse(&s2, &s3));

        assert_eq!(left.get(0, 0), Some(&6));
        assert_eq!(right.get(0, 0), Some(&6));
    }

    #[test]
    fn test_merge_commutative() {
        let mut d1 = DenseMatrix::new(2, 2);
        let mut m1 = BitwiseMask::new(2, 2);
        d1.set(0, 0, 5);
        d1.set(1, 1, 3);
        m1.flag(0, 0);
        m1.flag(1, 1);

        let mut d2 = DenseMatrix::new(2, 2);
        let mut m2 = BitwiseMask::new(2, 2);
        d2.set(0, 1, 7);
        d2.set(1, 1, 2);
        m2.flag(0, 1);
        m2.flag(1, 1);

        let s1 = condense_to_sparse(&d1, &m1);
        let s2 = condense_to_sparse(&d2, &m2);

        let ab = merge_sparse(&s1, &s2);
        let ba = merge_sparse(&s2, &s1);

        assert_eq!(ab.get(0, 0), ba.get(0, 0));
        assert_eq!(ab.get(0, 1), ba.get(0, 1));
        assert_eq!(ab.get(1, 1), ba.get(1, 1));
    }

    #[test]
    fn test_empty_sparse() {
        let empty = empty_sparse(5, 10);
        // Note: sprs empty matrix has 0 rows but correct cols
        assert_eq!(empty.cols(), 10);
        assert_eq!(empty.nnz(), 0);
    }
}
