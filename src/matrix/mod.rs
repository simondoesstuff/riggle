mod dense;
mod sparse;

pub use dense::{DenseMatrix, OverlapMatrix};
pub use sparse::{SparseMatrix, condense_to_sparse_no_mask};
