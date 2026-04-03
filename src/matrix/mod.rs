mod dense;
mod sparse;

pub use dense::{allocate_dense_accumulator, zero_flagged_regions, BitwiseMask, DenseMatrix};
pub use sparse::{condense_to_sparse, empty_sparse, merge_sparse, SparseMatrix};
