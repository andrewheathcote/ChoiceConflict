#!/bin/sh
# Install the EMC2 / DMC2 builds the archive uses, each into its own library: <CC_LIB, default ../lib>/<sha8>/EMC2.
#   sh builds/install_builds.sh [sha8 ...]      (default: every tarball in builds/)
# Needs R (>= 4.4), a C++17 compiler and EMC2's CRAN dependencies (see the archive's README.md); OpenMP is used when the
# compiler supports it (chains run in parallel threads). Each library gets EMC2_SHA.txt, which config.R's use_build() checks.
cd "$(dirname "$0")"; LIB=${CC_LIB:-$(cd .. && pwd)/lib}
[ $# -eq 0 ] && set -- $(ls EMC2_*.tar.gz | sed 's/EMC2_\(.*\)\.tar\.gz/\1/')
for s in "$@"; do
  full=$(awk -F, -v s="$s" '$1 == s {print $2}' BUILDS.csv); [ -n "$full" ] || { echo "unknown build $s"; continue; }
  want=$(awk -v f="EMC2_$s.tar.gz" '$2 == f {print $1}' MD5SUMS)
  got=$(Rscript -e "cat(tools::md5sum('EMC2_$s.tar.gz'))")
  [ -n "$want" ] && [ "$want" = "$got" ] || { echo "EMC2_$s.tar.gz: md5 does not match MD5SUMS, skipped"; continue; }
  dest=$LIB/$s; mkdir -p "$dest"; tmp=$(mktemp -d)
  tar xzf "EMC2_$s.tar.gz" -C "$tmp" && R CMD INSTALL --no-test-load -l "$dest" "$tmp/EMC2" > "$dest/install.log" 2>&1 \
    && echo "$full" > "$dest/EMC2_SHA.txt" && echo "installed $s" || echo "FAILED $s (see $dest/install.log)"
  rm -rf "$tmp"
done
