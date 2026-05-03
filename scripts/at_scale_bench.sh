#!/bin/bash
set -e

# Default parameters
SRC_FILES=10       # Number of files indexed into the database
SRC_RECORDS=100000 # Records per source file
QRY_RECORDS=10000  # Records in the single query file
TIMEOUT=300        # Default timeout in seconds
BATCH_SIZE=''      # Query batch size, None = all

# Parse command-line arguments
# -f / -n : source (database) shape
# -m      : query file record count
# -t      : timeout
# -b      : riggle query batch size
while getopts f:n:m:t:b: flag; do
	case "${flag}" in
	f) SRC_FILES=${OPTARG} ;;
	n) SRC_RECORDS=${OPTARG} ;;
	m) QRY_RECORDS=${OPTARG} ;;
	t) TIMEOUT=${OPTARG} ;;
	b) BATCH_SIZE=${OPTARG} ;;
	esac
done

DATA_DIR="data"
PREFIX="${DATA_DIR}/at_scale_bench_"
SRC_BED_DIR="${PREFIX}src_beds"
QRY_FILE="${PREFIX}query.bed.gz"
IDX_NATIVE="${PREFIX}riggle_idx"
IDX_IGD="${PREFIX}igd_idx"
IDX_GIGGLE="${PREFIX}giggle_idx"

echo "RUN: src_files=${SRC_FILES}, src_records=${SRC_RECORDS}, qry_records=${QRY_RECORDS}, timeout=${TIMEOUT}s"

# Generate and sort source BED files (seed 42)
if [ ! -d "$SRC_BED_DIR" ]; then
	mkdir -p "$SRC_BED_DIR"
	just gen -f "$SRC_FILES" -n "$SRC_RECORDS" --min-len 1 --max-len 50000 -c --sort -o "$SRC_BED_DIR"
fi

# Generate and sort single query BED file (seed 99 to differ from source)
if [ ! -f "$QRY_FILE" ]; then
	just gen -n "$QRY_RECORDS" -s 99 --min-len 1 --max-len 50000 -c --sort -o "$QRY_FILE"
fi

cleanup() {
	# not cleaning up indices for fast repeated runs
	rm -f "$QRY_FILE"
}
trap cleanup EXIT

echo
echo "==================================  Native  ================================== "
batchArg=''
[ -n "$BATCH_SIZE" ] && batchArg="--batch-size $BATCH_SIZE"
if [ ! -d "$IDX_NATIVE" ]; then
	echo ">>>   Building index   <<<"
	time timeout "$TIMEOUT" just run add -i "$SRC_BED_DIR" -d "$IDX_NATIVE" ${batchArg} || {
		rm -rf "$IDX_NATIVE"
		echo "[!] Native build timed out or failed"
	}
fi
[ -d "$IDX_NATIVE" ] && {
	du -sh "$IDX_NATIVE"
	echo ">>>   Querying   <<<"
	time timeout "$TIMEOUT" just run query -d "$IDX_NATIVE" -q "$QRY_FILE" -o /dev/null ${batchArg} || echo "[!] Native query timed out or failed"
}

echo
echo "==================================  IGD  ===================================== "
if [ ! -d "$IDX_IGD" ]; then
	echo ">>>   Building index   <<<"
	time timeout "$TIMEOUT" igd create "$SRC_BED_DIR" "$IDX_IGD" -s 0 || {
		rm -rf "$IDX_IGD"
		echo "[!] IGD build timed out or failed"
	}
fi
[ -d "$IDX_IGD" ] && {
	du -sh "$IDX_IGD"
	IGD_FILE=$(ls "$IDX_IGD"/*.igd 2>/dev/null | head -1)
	if [ -n "$IGD_FILE" ]; then
		echo ">>>   Querying   <<<"
		time timeout "$TIMEOUT" igd search "$IGD_FILE" -q "$QRY_FILE" >/dev/null || echo "[!] IGD query timed out or failed"
	else
		echo "[!] No .igd file found in $IDX_IGD"
	fi
}

# echo
# echo "==================================  Giggle  ================================== "
# if [ ! -f "$IDX_GIGGLE/root_ids.dat" ]; then
# 	echo ">>>   Building index   <<<"
# 	time timeout "$TIMEOUT" giggle index -s -f -i "${SRC_BED_DIR}/*" -o "$IDX_GIGGLE" || {
# 		rm -rf "$IDX_GIGGLE"
# 		echo "[!] Giggle build timed out or failed"
# 	}
# fi
# [ -f "$IDX_GIGGLE/root_ids.dat" ] && {
# 	du -sh "$IDX_GIGGLE"
# 	echo ">>>   Querying   <<<"
# 	time timeout "$TIMEOUT" giggle search -i "$IDX_GIGGLE" -q "$QRY_FILE" >/dev/null || echo "[!] Giggle query timed out or failed"
# }
