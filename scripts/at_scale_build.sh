#!/bin/bash
set -e

# Default parameters
NUM_FILES=2
NUM_RECORDS=100000
TIMEOUT=300 # Default timeout in seconds

# Parse command-line arguments (-t added for timeout)
while getopts f:n:t: flag; do
	case "${flag}" in
	f) NUM_FILES=${OPTARG} ;;
	n) NUM_RECORDS=${OPTARG} ;;
	t) TIMEOUT=${OPTARG} ;;
	esac
done

DATA_DIR="data"
PREFIX="${DATA_DIR}/at_scale_bench_"
BED_DIR="${PREFIX}beds"
IDX_NATIVE="${PREFIX}riggle_idx"
IDX_IGD="${PREFIX}igd_idx"
IDX_GIGGLE="${PREFIX}giggle_idx"

cleanup() {
	rm -rf "$BED_DIR" "$IDX_NATIVE" "$IDX_IGD" "$IDX_GIGGLE"
}
trap cleanup EXIT

# Initial cleanup & setup
cleanup
mkdir -p "$BED_DIR"

echo "RUN: files=${NUM_FILES}, records=${NUM_RECORDS}, timeout=${TIMEOUT}s"
just gen -f "$NUM_FILES" -n "$NUM_RECORDS" --sort -c -o "$BED_DIR"

echo
echo -e "\t ================= Native (Riggle) ================="
time timeout "$TIMEOUT" just run build -i "$BED_DIR" -o "$IDX_NATIVE" || echo "[!] Native timed out or failed"
[ -d "$IDX_NATIVE" ] && du -sh "$IDX_NATIVE"
rm -rf "$IDX_NATIVE"

echo
echo -e "\t ================= IGD ================="
time timeout "$TIMEOUT" igd create "$BED_DIR" "$IDX_IGD" -s 0 || echo "[!] IGD timed out or failed"
[ -d "$IDX_IGD" ] && du -sh "$IDX_IGD"
rm -rf "$IDX_IGD"

echo
echo -e "\t ================= Giggle ================="
time timeout "$TIMEOUT" giggle index -s -i "${BED_DIR}/*" -o "$IDX_GIGGLE" || echo "[!] Giggle timed out or failed"
[ -d "$IDX_GIGGLE" ] && du -sh "$IDX_GIGGLE"
rm -rf "$IDX_GIGGLE"
