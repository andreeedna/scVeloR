#!/usr/bin/env bash
# Copy the Rust core into the R package so that source tarballs
# (R CMD build, devtools::install, r-universe) are self-contained.
# Run after every change in velocore/ (the R tests check the copy is in sync).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
dest="$root/r-package/src/rust/velocore"
rm -rf "$dest"
mkdir -p "$dest"
cp -r "$root/velocore/src" "$dest/src"
sed -e '/^\[dev-dependencies\]/,$d' "$root/velocore/Cargo.toml" > "$dest/Cargo.toml"
echo "velocore copied to $dest"
