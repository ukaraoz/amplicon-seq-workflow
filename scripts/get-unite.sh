#!/bin/bash
# Download the UNITE general FASTA release (ITS reference for DADA2 assignTaxonomy) into <refdir>/UNITE_ITS/<release date>/
#
# usage: get-unite.sh <refdir> [doi]
#   refdir : reference folder; the release goes to <refdir>/UNITE_ITS/<date>/ (config: taxonomy.unite_dir)
#   doi    : default 10.15156/BIO/3301231 = UNITE 10.0 (19.02.2025), all eukaryotes, singletons set as RefS
#            (fungi only: 10.15156/BIO/3301229; list of releases: https://unite.ut.ee/repository.php)
#
# The DOI landing page (doi.plutof.ut.ee) is a JavaScript app; the file link comes from the PlutoF API behind it.
# The tarball holds two files: sh_general_release_dynamic_all_<date>.fasta (ITS region only) and ..._dev.fasta
# (untrimmed, keeps flanking 18S/5.8S); dada2-asv.R uses the _dev file. Re-running skips a release already present
# (README.<tarball name>.txt written last, as the completion marker).

set -euo pipefail

refdir=${1:?usage: get-unite.sh <refdir> [doi]}
doi=${2:-10.15156/BIO/3301231}
api="https://api.plutof.ut.ee/v1/public/dois/?identifier=$doi"

# release file name + download URL from the DOI record
meta=$(curl -fsSL "$api")
name=$(printf '%s' "$meta" | grep -o '"name": *"[^"]*\.tgz"' | head -1 | sed 's/.*"\([^"]*\)"$/\1/')
url=$(printf '%s' "$meta" | grep -o '"url": *"https://[^"]*\.tgz"' | head -1 | sed 's/.*"\([^"]*\)"$/\1/')
title=$(printf '%s' "$meta" | grep -o '"title": *"[^"]*"' | head -1 | sed 's/.*"\([^"]*\)"$/\1/')
[ -n "$name" ] && [ -n "$url" ] || { echo "no .tgz found in the DOI record: $api" >&2; exit 1; }

# sh_general_release_all_19.02.2025.tgz -> 2025-02-19 (folder name, as for the existing UNITE_ITS releases)
reldate=$(printf '%s' "$name" | grep -o '[0-9][0-9]\.[0-9][0-9]\.[0-9][0-9][0-9][0-9]' | head -1)
[ -n "$reldate" ] || { echo "no release date in file name: $name" >&2; exit 1; }
outdir=$refdir/UNITE_ITS/$(echo "$reldate" | awk -F. '{print $3"-"$2"-"$1}')
echo "$title, release $reldate (doi:$doi)"
echo "  $url -> $outdir"

# one README per downloaded tarball: the same date folder can hold the fungi-only and all-eukaryotes releases
readme=$outdir/README.${name%.tgz}.txt
if [ -s "$readme" ]; then
  echo "already present, nothing to do:"; ls -l "$outdir"
  exit 0
fi

mkdir -p "$outdir"
tmp=$(mktemp -d "$outdir/.download.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
curl -fsSL --retry 3 -o "$tmp/$name" "$url"
tar -xzf "$tmp/$name" -C "$tmp"
md5=$(md5sum "$tmp/$name" 2>/dev/null | cut -d' ' -f1 || md5 -q "$tmp/$name")

# sanity check and move into place
{
  echo "$title"
  echo "doi: https://doi.org/$doi"
  echo "url: $url"
  echo "file: $name  md5: $md5"
  echo "downloaded: $(date '+%Y-%m-%d') by get-unite.sh"
  echo
  echo "sequences per kingdom:"
} > "$tmp/README"
n=0
for f in "$tmp"/*.fasta; do
  nseq=$(grep -c '^>' "$f")
  [ "$nseq" -gt 0 ] || { echo "empty fasta: $(basename "$f")" >&2; exit 1; }
  echo "  $(basename "$f"): $nseq" >> "$tmp/README"
  grep '^>' "$f" | cut -d'|' -f5 | cut -d';' -f1 | sort | uniq -c | sort -rn | sed 's/^/    /' >> "$tmp/README"
  mv "$f" "$outdir/"
  n=$((n + 1))
done
[ "$n" -gt 0 ] || { echo "no fasta in $name" >&2; exit 1; }
mv "$tmp/README" "$readme"
chmod g+r "$outdir"/*
cat "$readme"
ls -l "$outdir"
