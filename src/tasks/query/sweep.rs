use crate::matrix::SparseMatrix;

/// Build a CSR sparse matrix from entries that are already sorted by (row, col).
pub(super) fn build_csr_from_sorted_entries(
    entries: &[(usize, usize, u32)],
    num_rows: usize,
    num_cols: usize,
) -> SparseMatrix {
    let mut indptr = vec![0usize; num_rows + 1];
    let mut indices = Vec::with_capacity(entries.len());
    let mut data = Vec::with_capacity(entries.len());

    for &(row, col, val) in entries {
        indices.push(col);
        data.push(val);
        indptr[row + 1] += 1;
    }
    for i in 0..num_rows {
        indptr[i + 1] += indptr[i];
    }

    sprs::CsMat::new((num_rows, num_cols), indptr, indices, data)
}
