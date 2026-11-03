#!/bin/bash
# Place DADA2 16S V4 ASVs into the SILVA 138.2 reference alignment (mothur MiSeq SOP approach) and build a tree.
#
# usage: REFDIR=<mothur SILVA folder> [V4START=13862 V4END=23444] align-asv-v4.sh <ASV.fa> [outdir] [processors]
#   (workflow.py tree step sets these from the config)
#   ASV.fa : DADA2 ASV fasta, primer-free V4 (515F/806R), any sequence names
#
# steps (https://mothur.org/wiki/miseq_sop/):
#   1. one-time: trim the SILVA nr alignment to primer-free V4 (pcr.seqs). Coordinates were found by aligning
#      E. coli J01859 cut between 515F and 806R (253 bp) to silva.nr_v138_2.align -> columns 13862-23444
#      (with primers: 11895-25318), the same columns the SOP gives.
#   2. align.seqs: kmer search (8-mers) for the closest template, Needleman-Wunsch, NAST re-gapping
#   3. screen.seqs: keep ASVs that span the whole V4 (columns 1-9583 of the trimmed reference)
#   4. filter.seqs(vertical=T, trump=.): drop all-gap columns
#   no unique.seqs: ASVs that become identical after trimming stay separate rows in the ASV table
#   5. FastTree (GTR) on the filtered alignment

set -euo pipefail

asvfa=$(realpath "$1")
outdir=${2:-$(dirname "$asvfa")}
procs=${3:-8}

mothur=${MOTHUR:-mothur}
fasttree=${FASTTREE:-FastTreeMP}
# reference folder (set by workflow.py from the config's tree.reference): holds silva.nr_v138_2.tgz
# (https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_2.tgz); the V4 alignment is built there once
refdir=${REFDIR:?set REFDIR to the folder with silva.nr_v138_2.tgz}
silva=$refdir/silva.nr_v138_2.align
v4ref=$refdir/silva.nr_v138_2.v4.align
# V4 columns of the SILVA alignment (config tree.v4_start / v4_end; defaults: primer-free 515F/806R)
v4start=${V4START:-13862}
v4end=${V4END:-23444}
v4ncol=$((v4end - v4start + 1))           # 9583

key=$(basename "$asvfa" .fa)
key=${key%.fasta}
mkdir -p "$outdir"
cd "$outdir"

# 1. V4 reference (one-time)
if [ ! -s "$v4ref" ]; then
  # full alignment (8 GB) is only needed here; unpack it from the tarball if it was deleted
  [ -s "$silva" ] || tar -xzf "$refdir/silva.nr_v138_2.tgz" -C "$refdir" silva.nr_v138_2.align
  (cd "$refdir" && "$mothur" "#pcr.seqs(fasta=$silva, start=$v4start, end=$v4end, keepdots=F, processors=$procs); rename.file(input=silva.nr_v138_2.pcr.align, new=$(basename "$v4ref"))")
fi

# short IDs for mothur; map back to the ASV sequence (= DADA2 ASV table row name)
# (sequences may be wrapped over several lines: writeXStringSet wraps at 80 columns)
awk -v map="$key.ids.txt" 'function out() { if(seq != "") { n++; id = sprintf("ASV%05d", n); print ">" id; print seq; print id "\t" seq > map } }
  BEGIN{print "asv_id\tasv_seq" > map}
  /^>/{out(); seq = ""; next} {seq = seq $0} END{out()}' "$asvfa" > "$key.ids.fa"

# 2-4. align, screen, filter
"$mothur" "#set.current(processors=$procs);
  align.seqs(fasta=$key.ids.fa, reference=$v4ref);
  summary.seqs(fasta=$key.ids.align);
  screen.seqs(fasta=$key.ids.align, start=1, end=$v4ncol);
  filter.seqs(fasta=$key.ids.good.align, vertical=T, trump=.);
  summary.seqs(fasta=$key.ids.good.filter.fasta)" > "$key.mothur.log"

nin=$(grep -c '^>' "$key.ids.fa")
nbad=$( [ -f "$key.ids.bad.accnos" ] && wc -l < "$key.ids.bad.accnos" || echo 0)
echo "ASVs in: $nin; removed by screen.seqs (not spanning V4): $nbad; see $key.ids.bad.accnos"
echo "low-similarity alignments (SimBtwnQuery&Template < 80, likely off-target):"
awk -F'\t' 'NR>1 && $16 < 80 {n++} END{print "  " n+0}' "$key.ids.align_report"

# 5. tree
"$fasttree" -nt -gtr "$key.ids.good.filter.fasta" > "$key.ids.good.filter.tree" 2> "$key.fasttree.log"
echo "tree: $outdir/$key.ids.good.filter.tree (tip labels = asv_id; map: $key.ids.txt)"

# same tree with tip labels = ASV sequence, matching the ASV table / phyloseq taxa names
Rscript -e "map = read.table('$key.ids.txt', sep = '\t', header = T, stringsAsFactors = F);
  tr = ape::read.tree('$key.ids.good.filter.tree');
  tr\$tip.label = map\$asv_seq[match(tr\$tip.label, map\$asv_id)];
  ape::write.tree(tr, '$key.tree')"
echo "tree with sequence labels: $outdir/$key.tree"
