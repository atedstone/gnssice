#!/usr/bin/env bash
# Download files listed in urls.txt (one URL per line) using curl.
# Skips blank lines and comments (lines starting with #).
#
# Usage:
#   chmod +x download_files.sh
#   ./download_files.sh [urls_file] [output_dir]
#
# Defaults:
#   urls_file  = urls.txt        (in the same directory as this script)
#   output_dir = ./downloads
# Insert your username-password combination below or save in a .netrc file
# Written by Andrew Tedstone with Claude AI, March 2026

set -euo pipefail

URLS_FILE="${1:-urls.txt}"
OUTPUT_DIR="${2:-./downloads}"

if [[ ! -f "$URLS_FILE" ]]; then
  echo "Error: URL file '$URLS_FILE' not found." >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

total=0
success=0
failed=0
failed_urls=()

while IFS= read -r url || [[ -n "$url" ]]; do
  # Skip blank lines and comments
  [[ -z "$url" || "$url" == \#* ]] && continue

  total=$((total + 1))
  filename=$(basename "$url")
  dest="$OUTPUT_DIR/$filename"

  echo "[${total}] Downloading: $filename"

  if curl \
      --silent \
      --show-error \
      --fail \
      --ftp-ssl \
      --user "user:password" \
      --output "$dest" \
      "$url"; then
    echo "    ✓ Saved to $dest"
    success=$((success + 1))
  else
    echo "    ✗ Failed: $url" >&2
    failed=$((failed + 1))
    failed_urls+=("$url")
    rm -f "$dest"   # remove incomplete file
  fi

done < "$URLS_FILE"

echo ""
echo "Done. ${success}/${total} files downloaded successfully to '$OUTPUT_DIR'."

if [[ ${#failed_urls[@]} -gt 0 ]]; then
  echo ""
  echo "Failed URLs (${failed}):"
  printf '  %s\n' "${failed_urls[@]}"
  exit 1
fi
