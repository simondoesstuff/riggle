//! eigenbed: Genomic interval similarity via truncated SVD embedding.
//!
//! Each BED file is encoded as a sparse binary vector over a 100 bp tile
//! space, the corpus is stacked into a sparse matrix X (files × tiles),
//! and truncated SVD is computed via eigendecomposition of the gram matrix
//! X @ X^T (files × files, typically small).  The resulting k-dimensional
//! embeddings support cosine-similarity search at query time with no
//! statistical model required.
//!
//! # Commands
//! - `index`: build the index from a corpus of BED files
//! - `query`: project one or more BED files and rank corpus similarity

use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use clap::{Parser, Subcommand};
use faer::Side;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sprs::TriMat;

use riggle::bench::nat_cmp;
use riggle::io::{BedParseError, parse_bed_file};

// ─── CLI ──────────────────────────────────────────────────────────────────────

#[derive(Parser)]
#[command(
    name = "eigenbed",
    about = "Genomic interval similarity via SVD embedding"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build an SVD embedding index from a corpus of BED files.
    ///
    /// Parses each BED file in parallel, encodes them as L2-normalized sparse
    /// binary vectors in 100 bp tile space, then decomposes the gram matrix to
    /// find the top-k eigenbeds (principal genomic variation directions).
    /// Saves the projection matrix V and pre-computed corpus embeddings.
    Index {
        /// BED files to index (plain or gzip-compressed).
        #[arg(required = true)]
        beds: Vec<PathBuf>,

        /// Output directory for the index.
        #[arg(short, long)]
        output: PathBuf,

        /// Tile size in base pairs.
        #[arg(short, long, default_value_t = 100)]
        tile_size: u32,

        /// Maximum number of SVD components to retain (k).
        #[arg(short = 'k', long, default_value_t = 128)]
        components: usize,

        /// Retain the fewest components that explain at least this fraction of
        /// variance (0 < v ≤ 1).  When set, --components acts as an upper cap.
        #[arg(short = 'v', long, value_parser = parse_variance, default_value = "0.95")]
        variance: Option<f32>,
    },

    /// Project BED file(s) through the index and emit their k-dimensional embeddings.
    ///
    /// Output is TSV with one line per file: path followed by k float values.
    /// Embeddings are L2-normalised unit vectors; cosine similarity = dot product.
    Embed {
        /// BED file(s) to embed (plain or gzip-compressed).
        #[arg(required = true)]
        beds: Vec<PathBuf>,

        /// Index directory produced by `index`.
        #[arg(short, long, default_value = "eigenbed_index")]
        index: PathBuf,
    },

    /// Project BED file(s) through the index and rank corpus by cosine similarity.
    Query {
        /// Query BED file(s) (plain or gzip-compressed).
        #[arg(required = true)]
        queries: Vec<PathBuf>,

        /// Index directory produced by `index`.
        #[arg(short, long, default_value = "eigenbed_index")]
        index: PathBuf,

        /// Maximum number of results to show per query.
        #[arg(short = 'n', long, default_value_t = 10)]
        top: usize,
    },
}

// ─── Index metadata ───────────────────────────────────────────────────────────

#[derive(Serialize, Deserialize)]
struct IndexMeta {
    tile_size: u32,
    components: usize,
    /// chromosome name → first global tile offset for that chromosome
    chrom_offsets: HashMap<String, usize>,
    total_tiles: usize,
    num_active_tiles: usize,
    /// paths of indexed BED files, in row order (row i = bed_files[i])
    bed_files: Vec<String>,
}

// ─── Error ────────────────────────────────────────────────────────────────────

#[derive(Debug)]
enum EigenError {
    Io(std::io::Error),
    Parse(BedParseError),
    Json(serde_json::Error),
    Other(String),
}

impl std::fmt::Display for EigenError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "IO error: {e}"),
            Self::Parse(e) => write!(f, "parse error: {e}"),
            Self::Json(e) => write!(f, "JSON error: {e}"),
            Self::Other(s) => write!(f, "{s}"),
        }
    }
}

impl From<std::io::Error> for EigenError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}
impl From<BedParseError> for EigenError {
    fn from(e: BedParseError) -> Self {
        Self::Parse(e)
    }
}
impl From<serde_json::Error> for EigenError {
    fn from(e: serde_json::Error) -> Self {
        Self::Json(e)
    }
}

fn parse_variance(s: &str) -> Result<f32, String> {
    let v: f32 = s
        .parse()
        .map_err(|_| format!("'{s}' is not a valid number"))?;
    if v > 0.0 && v <= 1.0 {
        Ok(v)
    } else {
        Err(format!("variance must be in (0, 1], got {v}"))
    }
}

// ─── Entry point ─────────────────────────────────────────────────────────────

fn main() {
    let cli = Cli::parse();
    let result = match cli.command {
        Commands::Index {
            beds,
            output,
            tile_size,
            components,
            variance,
        } => cmd_index(beds, output, tile_size, components, variance),
        Commands::Embed { beds, index } => cmd_embed(beds, index),
        Commands::Query {
            queries,
            index,
            top,
        } => cmd_query(queries, index, top),
    };
    if let Err(e) = result {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

// ─── Index ───────────────────────────────────────────────────────────────────

fn cmd_index(
    beds: Vec<PathBuf>,
    output: PathBuf,
    tile_size: u32,
    max_components: usize,
    variance: Option<f32>,
) -> Result<(), EigenError> {
    let m = beds.len();
    if m < 2 {
        return Err(EigenError::Other(
            "need at least 2 BED files to build an index".into(),
        ));
    }
    let max_k = max_components.min(m - 1);

    let mp = MultiProgress::new();
    let bar_style = ProgressStyle::with_template(
        "{spinner:.cyan} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}",
    )
    .unwrap()
    .progress_chars("█▉▊▋▌▍▎▏  ");

    let pb_parse = mp.add(ProgressBar::new(m as u64));
    pb_parse.set_style(bar_style.clone());
    pb_parse.set_message("parsing BED files");

    // Parse every BED file in parallel; retain the file-index→intervals mapping.
    let parsed: Vec<(usize, HashMap<String, Vec<riggle::core::Interval>>)> = beds
        .par_iter()
        .enumerate()
        .map(|(sid, path)| {
            let ivs = parse_bed_file(path, sid as u32)?;
            pb_parse.inc(1);
            Ok::<_, BedParseError>((sid, ivs))
        })
        .collect::<Result<_, _>>()?;
    pb_parse.finish_with_message(format!("parsed {m} files"));

    let pb_chrom = mp.add(ProgressBar::new_spinner());
    pb_chrom.set_message("computing chromosome tile offsets");

    // Max observed endpoint per chromosome → size each chromosome's tile span.
    let mut chrom_max: HashMap<String, u32> = HashMap::new();
    for (_, ivs_by_chrom) in &parsed {
        for (chrom, ivs) in ivs_by_chrom {
            let max_end = ivs.iter().map(|iv| iv.end).max().unwrap_or(0);
            chrom_max
                .entry(chrom.clone())
                .and_modify(|e| *e = (*e).max(max_end))
                .or_insert(max_end);
        }
    }

    // Assign global tile offsets in natural chromosome order.
    let mut chroms: Vec<String> = chrom_max.keys().cloned().collect();
    chroms.sort_by(|a, b| nat_cmp(a, b));

    let mut chrom_offsets: HashMap<String, usize> = HashMap::new();
    let mut total_tiles = 0usize;
    for chrom in &chroms {
        chrom_offsets.insert(chrom.clone(), total_tiles);
        total_tiles += (chrom_max[chrom] as usize).div_ceil(tile_size as usize);
    }
    pb_chrom.finish_with_message(format!(
        "{} chromosomes → {} total tiles",
        chroms.len(),
        total_tiles
    ));

    let pb_tiles = mp.add(ProgressBar::new(2 * m as u64));
    pb_tiles.set_style(bar_style.clone());
    pb_tiles.set_message("building sparse tile matrix");

    // For each file, collect the sorted unique set of global tile indices it covers.
    let file_tiles: Vec<Vec<usize>> = parsed
        .par_iter()
        .map(|(_, ivs_by_chrom)| {
            let mut set: HashSet<usize> = HashSet::new();
            for (chrom, ivs) in ivs_by_chrom {
                let Some(&chrom_off) = chrom_offsets.get(chrom) else {
                    continue;
                };
                for iv in ivs {
                    let t0 = iv.start as usize / tile_size as usize;
                    let t1 = (iv.end as usize - 1) / tile_size as usize;
                    for t in t0..=t1 {
                        set.insert(chrom_off + t);
                    }
                }
            }
            let mut v: Vec<usize> = set.into_iter().collect();
            v.sort_unstable();
            pb_tiles.inc(1);
            v
        })
        .collect();

    // Union of all covered tiles → compact index space for the stored matrices.
    let mut active_set: HashSet<usize> = HashSet::new();
    for tiles in &file_tiles {
        active_set.extend(tiles.iter().copied());
    }
    let mut active_tiles: Vec<usize> = active_set.into_iter().collect();
    active_tiles.sort_unstable();
    let n_active = active_tiles.len();

    // Map global tile → compact row index (for O(log n) lookup at query time).
    let tile_to_compact: HashMap<usize, usize> = active_tiles
        .iter()
        .enumerate()
        .map(|(i, &t)| (t, i))
        .collect();

    // Build the (m × n_active) sparse matrix X in triplet form.
    // Each row is L2-normalised: val = 1/sqrt(|tiles for this file|).
    let mut tri: TriMat<f32> = TriMat::new((m, n_active));
    for (i, tiles) in file_tiles.iter().enumerate() {
        if tiles.is_empty() {
            continue;
        }
        let val = 1.0_f32 / (tiles.len() as f32).sqrt();
        for &t in tiles {
            if let Some(&ci) = tile_to_compact.get(&t) {
                tri.add_triplet(i, ci, val);
            }
        }
        pb_tiles.inc(1);
    }
    let x = tri.to_csr::<usize>();
    pb_tiles.finish_with_message(format!(
        "{n_active} active tiles ({:.2}% of genome covered)",
        100.0 * n_active as f64 / total_tiles as f64
    ));

    let pb_gram = mp.add(ProgressBar::new_spinner());
    pb_gram.set_message(format!("gram matrix ({m}×{m}) + eigendecomposition"));
    pb_gram.enable_steady_tick(std::time::Duration::from_millis(80));

    // gram = X @ X^T  (m × m).
    // Since rows of X are L2-normalised, gram[i,j] = cosine_similarity(file_i, file_j).
    // Computed as: X (CSR, m×n) × X^T (CSR, n×m).
    let gram_sparse = {
        let xt_csr = x.transpose_view().to_csr(); // n_active × m, CSR
        &x * &xt_csr // m × m
    };

    // Fill only the lower triangle into a column-major Vec, then wrap as a faer
    // MatRef.  self_adjoint_eigen reads only one triangle, so no symmetry
    // enforcement loop is needed.
    let mut gram_data = vec![0.0f32; m * m]; // column-major: data[row + col*m]
    for (&val, (r, c)) in gram_sparse.iter() {
        if r >= c {
            gram_data[r + c * m] = val;
        }
    }
    let gram = faer::mat::MatRef::<f32>::from_column_major_slice(&gram_data, m, m);

    let eig = gram
        .self_adjoint_eigen(Side::Lower)
        .map_err(|_| EigenError::Other("eigendecomposition failed".into()))?;

    // faer returns eigenvalues in ascending order; sort descending for variance selection.
    let mut order: Vec<(usize, f32)> = eig
        .S()
        .column_vector()
        .iter()
        .copied()
        .enumerate()
        .collect();
    order.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let total_var: f32 = order.iter().map(|(_, v)| v.max(0.0)).sum();

    // Determine k: either the smallest number of components that meets the variance
    // target, or max_k if no target is given—whichever is smaller.
    let k = match variance {
        Some(target) => {
            let mut cum = 0.0f32;
            let mut selected = 1usize;
            for (j, &(_, λ)) in order.iter().enumerate() {
                cum += λ.max(0.0);
                selected = j + 1;
                if cum / total_var.max(1e-10) >= target {
                    break;
                }
            }
            selected.min(max_k)
        }
        None => max_k,
    };
    order.truncate(k);

    let top_var: f32 = order.iter().map(|(_, v)| v.max(0.0)).sum();
    pb_gram.finish_with_message(format!(
        "k={k} components explain {:.1}% of variance{}",
        100.0 * top_var / total_var.max(1e-10),
        variance
            .map(|v| format!("  (target {:.0}%)", v * 100.0))
            .unwrap_or_default(),
    ));

    // U_m: (m × k) row-major flat vec; u_m[i*k + comp] = eigenvector comp, file i.
    // sigma: corresponding singular values (√eigenvalues).
    let evecs = eig.U();
    let mut u_m = vec![0.0f32; m * k];
    let sigma: Vec<f32> = order
        .iter()
        .enumerate()
        .map(|(j, &(idx, λ))| {
            for i in 0..m {
                u_m[i * k + j] = evecs[(i, idx)];
            }
            λ.max(0.0_f32).sqrt()
        })
        .collect();

    let pb_proj = mp.add(ProgressBar::new(m as u64));
    pb_proj.set_style(bar_style.clone());
    pb_proj.set_message(format!("projection matrix ({n_active} tiles × {k})"));

    // Projection matrix V = X^T @ U_m @ diag(1/σ)  (n_active × k).
    //
    // Computed via a single pass over the non-zeros of X:
    //   V[j, :] += X[i, j] · U_m[i, :]   for every non-zero (i, j) in X
    // then each column is scaled by 1/σ.
    let mut proj = vec![0.0f32; n_active * k];
    for (i, row) in x.outer_iterator().enumerate() {
        let u_row = &u_m[i * k..(i + 1) * k];
        for (j, &x_val) in row.iter() {
            let proj_row = &mut proj[j * k..(j + 1) * k];
            for comp in 0..k {
                proj_row[comp] += x_val * u_row[comp];
            }
        }
        pb_proj.inc(1);
    }
    pb_proj.finish_with_message("projection matrix done");
    for comp in 0..k {
        let scale = if sigma[comp] > 1e-10 {
            1.0 / sigma[comp]
        } else {
            0.0
        };
        for tile in 0..n_active {
            proj[tile * k + comp] *= scale;
        }
    }

    // Corpus embeddings = U_m @ diag(σ)  (m × k), then L2-normalised per row.
    // This is equivalent to X @ V, the projection of each corpus file.
    let mut embeddings = vec![0.0f32; m * k];
    for i in 0..m {
        let mut norm_sq = 0.0f32;
        for comp in 0..k {
            let v = u_m[i * k + comp] * sigma[comp];
            embeddings[i * k + comp] = v;
            norm_sq += v * v;
        }
        let norm = norm_sq.sqrt();
        if norm > 1e-10 {
            for comp in 0..k {
                embeddings[i * k + comp] /= norm;
            }
        }
    }

    let pb_write = mp.add(ProgressBar::new_spinner());
    pb_write.set_message(format!("writing index to {}", output.display()));
    fs::create_dir_all(&output)?;

    // meta.json
    let meta = IndexMeta {
        tile_size,
        components: k,
        chrom_offsets,
        total_tiles,
        num_active_tiles: n_active,
        bed_files: beds
            .iter()
            .map(|p| p.to_string_lossy().into_owned())
            .collect(),
    };
    fs::write(
        output.join("meta.json"),
        serde_json::to_string_pretty(&meta)?,
    )?;

    // active_tiles.bin — sorted list of global tile indices (u64 LE)
    {
        let mut f = BufWriter::new(File::create(output.join("active_tiles.bin"))?);
        for &t in &active_tiles {
            f.write_all(&(t as u64).to_le_bytes())?;
        }
    }

    // projection.bin — V matrix (n_active × k, f32 LE row-major)
    write_f32_binary(&output.join("projection.bin"), &proj)?;

    // embeddings.bin — corpus embeddings (m × k, f32 LE row-major, L2-normalised rows)
    write_f32_binary(&output.join("embeddings.bin"), &embeddings)?;

    pb_write.finish_with_message(format!(
        "index written  (projection: {:.1} MB, embeddings: {:.1} KB)",
        (n_active * k * 4) as f64 / 1e6,
        (m * k * 4) as f64 / 1e3,
    ));
    Ok(())
}

// ─── Shared index loading + projection ───────────────────────────────────────

struct LoadedIndex {
    meta: IndexMeta,
    /// Sorted global tile indices present in the corpus (binary-searchable).
    active_tiles: Vec<u64>,
    /// Projection matrix V (n_active × k, row-major f32).
    proj: Vec<f32>,
    /// L2-normalised corpus embeddings (m × k, row-major f32); empty if not loaded.
    embeddings: Vec<f32>,
}

impl LoadedIndex {
    fn load(index_dir: &PathBuf, with_embeddings: bool) -> Result<Self, EigenError> {
        let meta: IndexMeta =
            serde_json::from_str(&fs::read_to_string(index_dir.join("meta.json"))?)?;

        let at_bytes = fs::read(index_dir.join("active_tiles.bin"))?;
        let active_tiles: Vec<u64> = at_bytes
            .chunks_exact(8)
            .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
            .collect();

        let proj = read_f32_binary(&index_dir.join("projection.bin"))?;

        let embeddings = if with_embeddings {
            read_f32_binary(&index_dir.join("embeddings.bin"))?
        } else {
            Vec::new()
        };

        Ok(Self {
            meta,
            active_tiles,
            proj,
            embeddings,
        })
    }

    /// Project a BED file into the k-dimensional embedding space.
    ///
    /// Returns `(n_tiles, hits, embedding)` where `n_tiles` is the total number
    /// of tiles the file covers, `hits` is how many of those appear in the index,
    /// and `embedding` is the L2-normalised projection vector (length k).
    fn project_bed(&self, path: &PathBuf) -> Result<(usize, usize, Vec<f32>), BedParseError> {
        let k = self.meta.components;
        let tile_size = self.meta.tile_size;
        let ivs_by_chrom = parse_bed_file(path, 0)?;

        let mut tiles: HashSet<u64> = HashSet::new();
        for (chrom, ivs) in &ivs_by_chrom {
            let Some(&chrom_off) = self.meta.chrom_offsets.get(chrom) else {
                continue;
            };
            for iv in ivs {
                let t0 = iv.start as usize / tile_size as usize;
                let t1 = (iv.end as usize - 1) / tile_size as usize;
                for t in t0..=t1 {
                    tiles.insert((chrom_off + t) as u64);
                }
            }
        }

        let mut embedding = vec![0.0f32; k];
        let mut hits = 0usize;
        for &global_tile in &tiles {
            if let Ok(pos) = self.active_tiles.binary_search(&global_tile) {
                let base = pos * k;
                for comp in 0..k {
                    embedding[comp] += self.proj[base + comp];
                }
                hits += 1;
            }
        }

        let norm: f32 = embedding.iter().map(|&v| v * v).sum::<f32>().sqrt();
        if norm > 1e-10 {
            for v in &mut embedding {
                *v /= norm;
            }
        }

        Ok((tiles.len(), hits, embedding))
    }
}

// ─── Embed ────────────────────────────────────────────────────────────────────

fn cmd_embed(beds: Vec<PathBuf>, index_dir: PathBuf) -> Result<(), EigenError> {
    let idx = LoadedIndex::load(&index_dir, false)?;
    let k = idx.meta.components;

    // Header
    let header_comps: Vec<String> = (1..=k).map(|i| i.to_string()).collect();
    println!("# file\t{}", header_comps.join("\t"));

    let pb = ProgressBar::new(beds.len() as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.cyan} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}",
        )
        .unwrap()
        .progress_chars("█▉▊▋▌▍▎▏  "),
    );
    pb.set_message("embedding");

    let results: Vec<_> = beds
        .par_iter()
        .map(|path| {
            let name = path.to_string_lossy().into_owned();
            let (_, _, embedding) = idx.project_bed(path)?;
            pb.inc(1);
            Ok::<_, BedParseError>((name, embedding))
        })
        .collect::<Result<_, _>>()?;
    pb.finish_and_clear();

    for (name, embedding) in &results {
        let vals: Vec<String> = embedding.iter().map(|v| format!("{v:.6}")).collect();
        println!("{}\t{}", name, vals.join("\t"));
    }

    Ok(())
}

// ─── Query ────────────────────────────────────────────────────────────────────

fn cmd_query(queries: Vec<PathBuf>, index_dir: PathBuf, top: usize) -> Result<(), EigenError> {
    let idx = LoadedIndex::load(&index_dir, true)?;
    let k = idx.meta.components;
    let m = idx.meta.bed_files.len();

    let pb = ProgressBar::new(queries.len() as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.cyan} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}",
        )
        .unwrap()
        .progress_chars("█▉▊▋▌▍▎▏  "),
    );
    pb.set_message("querying");

    let results: Vec<_> = queries
        .par_iter()
        .map(|qpath| {
            let name = qpath.to_string_lossy().into_owned();
            let (n_tiles, hits, embedding) = idx.project_bed(qpath)?;

            let mut sims: Vec<(usize, f32)> = (0..m)
                .map(|i| {
                    let base = i * k;
                    let sim = (0..k)
                        .map(|c| embedding[c] * idx.embeddings[base + c])
                        .sum();
                    (i, sim)
                })
                .collect();
            sims.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            sims.truncate(top);
            pb.inc(1);

            Ok::<_, BedParseError>((name, n_tiles, hits, sims))
        })
        .collect::<Result<_, _>>()?;
    pb.finish_and_clear();

    for (name, n_tiles, hits, sims) in &results {
        println!("# query: {}  ({} tiles, {} in index)", name, n_tiles, hits);
        for &(idx_i, sim) in sims {
            println!("{:.4}\t{}", sim, idx.meta.bed_files[idx_i]);
        }
    }

    Ok(())
}

// ─── I/O helpers ─────────────────────────────────────────────────────────────

fn write_f32_binary(path: &std::path::Path, data: &[f32]) -> Result<(), EigenError> {
    let mut f = BufWriter::new(File::create(path)?);
    for &v in data {
        f.write_all(&v.to_le_bytes())?;
    }
    Ok(())
}

fn read_f32_binary(path: &std::path::Path) -> Result<Vec<f32>, EigenError> {
    let bytes = fs::read(path)?;
    Ok(bytes
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
        .collect())
}
