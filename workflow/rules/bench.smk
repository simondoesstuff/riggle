# At-scale benchmark: index and query timing for Chuckle, IGD, and Giggle.
#
# Each trial is uniquely identified by (src_files × src_records, qry_records) and
# produces an independent JSON artifact.  Add new configurations to BENCH_CONFIGS;
# Snakemake will only run the trials that do not yet have outputs.
#
# Run all configured trials:
#   snakemake bench_all
#
# Run a specific trial:
#   snakemake data/bench/10f_50000n_20000m.json
#
# Produce plots from completed trials:
#   snakemake data/plots/bench_index.png

# Quick-run configs (few minutes total; scale up when stable).
# Total source intervals = src_files × src_records.
BENCH_CONFIGS = [
    {"src_files": 2,   "src_records": 1_000,   "qry_records": 500},
    {"src_files": 5,   "src_records": 5_000,   "qry_records": 2_000},
    {"src_files": 10,  "src_records": 20_000,  "qry_records": 10_000},
    {"src_files": 20,  "src_records": 50_000,  "qry_records": 25_000},
    {"src_files": 50,  "src_records": 100_000, "qry_records": 50_000},
    {"src_files": 100, "src_records": 200_000, "qry_records": 100_000},
]

# Full-scale configs (hours; uncomment when ready to scale up):
# BENCH_CONFIGS = [
#     {"src_files": 2,    "src_records": 5_000,      "qry_records": 1_000},
#     {"src_files": 2,    "src_records": 50_000,     "qry_records": 10_000},
#     {"src_files": 5,    "src_records": 200_000,    "qry_records": 100_000},
#     {"src_files": 10,   "src_records": 200_000,    "qry_records": 200_000},
#     {"src_files": 20,   "src_records": 200_000,    "qry_records": 400_000},
#     {"src_files": 50,   "src_records": 200_000,    "qry_records": 1_000_000},
#     {"src_files": 100,  "src_records": 200_000,    "qry_records": 2_000_000},
#     {"src_files": 200,  "src_records": 200_000,    "qry_records": 4_000_000},
#     {"src_files": 500,  "src_records": 200_000,    "qry_records": 10_000_000},
#     {"src_files": 750,  "src_records": 200_000,    "qry_records": 150_000_000},
#     {"src_files": 1000, "src_records": 1_000_000,  "qry_records": 100_000_000},
#     {"src_files": 2000, "src_records": 1_000_000,  "qry_records": 200_000_000},
# ]

BENCH_TIMEOUT = 60   # seconds per phase (increase for full-scale)
BENCH_BATCH_SIZE = 400

BENCH_ALL = [
    f"data/bench/{c['src_files']}f_{c['src_records']}n_{c['qry_records']}m.json"
    for c in BENCH_CONFIGS
]


rule bench_trial:
    """One at-scale benchmark trial → structured JSON with per-tool timing and index size.

    Indexes are built in a temp directory and deleted after each trial — only the JSON
    result is kept.  This avoids accumulating large synthetic indexes on disk.
    """
    output:
        "data/bench/{src_files}f_{src_records}n_{qry_records}m.json",
    wildcard_constraints:
        src_files=r"\d+",
        src_records=r"\d+",
        qry_records=r"\d+",
    params:
        timeout=BENCH_TIMEOUT,
        batch_size=BENCH_BATCH_SIZE,
    script:
        "../scripts/bench_trial.py"


rule bench_all:
    """Trigger all configured benchmark trials."""
    input:
        BENCH_ALL,


rule bench_plot:
    """Asymptotic complexity plots (index time, query time, index size) from all trials."""
    input:
        BENCH_ALL,
    output:
        index="data/plots/bench_index.png",
        query="data/plots/bench_query.png",
        size="data/plots/bench_size.png",
    shell:
        "uv run workflow/scripts/bench_plot.py {input}"
        " --index {output.index}"
        " --query {output.query}"
        " --size {output.size}"
