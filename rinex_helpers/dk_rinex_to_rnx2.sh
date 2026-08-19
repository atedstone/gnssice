#!/usr/bin/env bash
# Convert a folder of 1-second RINEX3 observation files to 10-second RINEX2
# files using gfzrnx, with automatic RINEX2 filename convention output.
#
# Supports all combinations of input file types:
#   .crx.gz  (Hatanaka + gzip)  → gunzip to tmp, crx2rnx to tmp, gfzrnx
#   .crx     (Hatanaka only)    → crx2rnx to tmp, gfzrnx
#   .rnx.gz  (gzip only)        → gunzip to tmp, gfzrnx
#   .rnx     (plain)            → gfzrnx directly
#
# Usage:
#   chmod +x convert_rinex3_to_rinex2_10s.sh
#   ./convert_rinex3_to_rinex2_10s.sh [input_dir] [output_dir]
#
# Defaults:
#   input_dir  = .
#   output_dir = .
#
# Author: Andrew Tedstone, March 2026 using Claude AI

set -euo pipefail

INPUT_DIR="${1:-.}"
OUTPUT_DIR="${2:-.}"
SAMPLING=10   # output sampling rate in seconds
TMP_DIR=$(mktemp -d)

cleanup() { rm -rf "$TMP_DIR"; }
trap cleanup EXIT

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: input directory '$INPUT_DIR' not found." >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Collect all RINEX3 observation files (compressed or plain)
mapfile -t files < <(find "$INPUT_DIR" -maxdepth 1 -type f \
  \( -name "*_MO.crx.gz" -o -name "*_MO.rnx.gz" \
     -o -name "*_MO.crx"  -o -name "*_MO.rnx"  \) | sort)

if [[ ${#files[@]} -eq 0 ]]; then
  echo "No RINEX3 observation files found in '$INPUT_DIR'." >&2
  exit 1
fi

total=${#files[@]}
success=0
failed=0
failed_files=()

echo "Found ${total} file(s) in '$INPUT_DIR'."
echo "Converting to ${SAMPLING}s RINEX2 → '$OUTPUT_DIR'"
echo ""

for filepath in "${files[@]}"; do
  filename=$(basename "$filepath")
  echo "[$(( success + failed + 1 ))/${total}] Processing: $filename"

  # Decompress to a temp file so gfzrnx always receives a real file path.
  # After each iteration the tmp dir is cleared, keeping disk use low.
  gfzrnx_input="$filepath"   # default: use file directly (plain .rnx)

  case "$filename" in
    *.crx.gz)
      tmp_rnx="$TMP_DIR/${filename%.crx.gz}.rnx"
      if ! gunzip -c "$filepath" | crx2rnx - > "$tmp_rnx"; then
        echo "    ✗ Failed to decompress/decode: $filename" >&2
        failed=$(( failed + 1 )); failed_files+=("$filename"); continue
      fi
      gfzrnx_input="$tmp_rnx"
      ;;
    *.crx)
      tmp_rnx="$TMP_DIR/${filename%.crx}.rnx"
      if ! crx2rnx - < "$filepath" > "$tmp_rnx"; then
        echo "    ✗ Failed to decode Hatanaka: $filename" >&2
        failed=$(( failed + 1 )); failed_files+=("$filename"); continue
      fi
      gfzrnx_input="$tmp_rnx"
      ;;
    *.rnx.gz)
      tmp_rnx="$TMP_DIR/${filename%.gz}"
      if ! gunzip -c "$filepath" > "$tmp_rnx"; then
        echo "    ✗ Failed to decompress: $filename" >&2
        failed=$(( failed + 1 )); failed_files+=("$filename"); continue
      fi
      gfzrnx_input="$tmp_rnx"
      ;;
    *.rnx)
      gfzrnx_input="$filepath"   # no decompression needed
      ;;
    *)
      echo "    ✗ Unrecognised file type, skipping: $filename" >&2
      failed=$(( failed + 1 )); failed_files+=("$filename"); continue
      ;;
  esac

  if gfzrnx \
      -finp  "$gfzrnx_input" \
      -fout  "${OUTPUT_DIR}/::RX2::" \
      -vo    2 \
      -smp   "$SAMPLING" \
      -f; then
    echo "    ✓ Done"
    success=$(( success + 1 ))
  else
    echo "    ✗ gfzrnx failed: $filename" >&2
    failed=$(( failed + 1 ))
    failed_files+=("$filename")
  fi

  # Clear tmp dir between iterations to avoid accumulating large files
  rm -f "$TMP_DIR"/*

done

echo ""
echo "Done. ${success}/${total} file(s) converted successfully to '$OUTPUT_DIR'."

if [[ ${#failed_files[@]} -gt 0 ]]; then
  echo ""
  echo "Failed files (${failed}):"
  printf '  %s\n' "${failed_files[@]}"
  exit 1
fi
