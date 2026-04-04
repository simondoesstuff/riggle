mod dense;
mod sparse;

pub use dense::{BitwiseMask, DenseMatrix, allocate_dense_accumulator, zero_flagged_regions};
pub use sparse::{SparseMatrix, condense_to_sparse, empty_sparse, merge_sparse};
