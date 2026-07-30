use std::cell::RefCell;

use crate::core::Interval;
use crate::io::MappedJumpTable;
use crate::matrix::{DenseMatrix, OverlapMatrix, SparseMatrix};
use crate::sweep::query_sweep;

thread_local! {
    static THREAD_LOCAL_BUFFER: RefCell<Option<(DenseMatrix, OverlapMatrix)>> =
        const { RefCell::new(None) };
}

/// Run a single sweep block using thread-local scratch buffers.
///
/// Returns `(count_matrix, overlap_matrix)` where `overlap_matrix` accumulates
/// the exact base-pair overlap in 100 bp bins for each (query, DB-source) pair.
pub(super) fn run_sweep_block(
    db_intervals: &[Interval],
    layer_max_size: u32,
    jump_table: &MappedJumpTable,
    query_block: &[Interval],
    num_queries: usize,
    num_sources: usize,
) -> (DenseMatrix, OverlapMatrix) {
    THREAD_LOCAL_BUFFER.with(|cell| {
        let mut borrow = cell.borrow_mut();
        let (counts, overlap) = borrow.get_or_insert_with(|| {
            (
                DenseMatrix::new(num_queries, num_sources),
                OverlapMatrix::new(num_queries, num_sources),
            )
        });

        if counts.num_rows() != num_queries || counts.num_cols() != num_sources {
            counts.resize_and_zero(num_queries, num_sources);
            overlap.resize_and_zero(num_queries, num_sources);
        }

        query_sweep(db_intervals, layer_max_size, jump_table, query_block, counts, overlap);

        let result_counts = counts.clone();
        let result_overlap = overlap.clone();
        counts.resize_and_zero(num_queries, num_sources);
        overlap.resize_and_zero(num_queries, num_sources);
        (result_counts, result_overlap)
    })
}

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
