# DADA2 workflow: paired-end amplicon reads -> ASV table + taxonomy
#
# version 2.18
#
# usage: Rscript dada2-asv.R --config=<key>.config.json [--threads=N] [--nsamples=N] [--no-taxonomy]
#                            [--resume[=outputs]] [--pdf]
#   normally run by workflow.py (dada2 and phyloseq steps), which writes the resolved config
#   --config=FILE      resolved per-locus config (workflow.py): paths, sample naming, locus parameters, references
#   --threads=N        threads (default: config params.threads, else $SLURM_CPUS_PER_TASK, else 4)
#   --nsamples=N       test run on the first N samples (of each run); outputs go to <results>/test-nN and <work>/test-nN
#   --no-taxonomy      stop after the seqtab checkpoint (before assignTaxonomy)
#   --resume           load the seqtab checkpoint (<key>.dada2.seqtab.RData) and run only taxonomy + outputs
#   --resume=outputs   load the final <key>.dada2.RData (with taxonomy) and only rewrite the outputs
#                      (tables, fasta, phyloseq, failed samples), e.g. after building the tree
#   --pdf              also write the raw read quality profiles and error-model plots as PDF (default: PNG only)
#   R1 only or merged pairs: config params.merge (workflow.py --merge / --no-merge); R1-only keys end in ".R1"
#
# failed samples (<key>.failed_samples.txt, with the reason): < mincount reads after quality filtering (dropped
#   before denoising), < minfinal reads AFTER the config's keep_taxa / remove_taxa / require_ranks filter (with
#   keep_taxa Fungi: fewer than minfinal fungal reads), or less than min_taxa_kept_fraction of the sample's reads
#   kept by that filter (default 0 = no check; e.g. 0.5 with keep_taxa Fungi: < 50% fungal reads). Failed samples
#   are kept in <key>.phyloseq.rds but removed from <key>.phyloseq.filtered.rds
#
# inputs  : paths.metadata (tab-separated; samples.id_column = sample id); fastq from samples.manifest (R1, R2
#           columns) or paths.raw/<samples.r1|r2>
# outputs : paths.results (QC plots, read tracking, ASV table/fasta, taxonomy, RData, phyloseq rds; the filtered phyloseq
#           also as text files per slot, <key>.filtered.*, with <key>.filtered.README.txt; the package step bundles them
#           into <key>.filtered.tar.gz)
#           paths.work    (cutadapt + filtered reads)

packages = c("dada2", "ShortRead", "Biostrings", "ggplot2", "parallel", "phyloseq")
invisible(lapply(packages, library, character.only = TRUE))

allargs = commandArgs(trailingOnly = TRUE)
opts = allargs[grepl("^--", allargs)]
optval = function(name, default) {
  hit = sub(paste0("^--", name, "="), "", opts[grepl(paste0("^--", name, "="), opts)])
  if(length(hit) > 0) hit[1] else default
}
config_file = optval("config", NA)
if(is.na(config_file)) stop("usage: Rscript dada2-asv.R --config=<key>.config.json [options] (see workflow.py)")
cfg = jsonlite::read_json(config_file, simplifyVector = TRUE)
locus = cfg$locus
run = cfg$run
key = cfg$key

cutadapt = Sys.which("cutadapt")
# threads: --threads, else the config (dada2.threads / locus_settings.<locus>.threads), else the SLURM allocation, else 4
nthreads = as.integer(optval("threads", if(!is.null(cfg$params$threads)) cfg$params$threads else Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
nsamples = as.integer(optval("nsamples", NA))
run_taxonomy = !("--no-taxonomy" %in% opts)
write_pdf = "--pdf" %in% opts
resume_mode = if("--resume" %in% opts) "taxonomy" else optval("resume", NA)
stopifnot(resume_mode %in% c(NA, "taxonomy", "outputs"))
resume = !is.na(resume_mode)
mincount = cfg$params$mincount   # samples with fewer filtered reads are dropped before denoising
minfinal = cfg$params$minfinal   # samples with fewer reads after the taxa filter (e.g. fungal reads) fail
# samples keeping a smaller fraction of their reads after the taxa filter fail; default 0 (no check)
min_kept = if(is.null(cfg$params$min_taxa_kept_fraction)) 0 else cfg$params$min_taxa_kept_fraction
set.seed(100)
# run per-sample jobs in parallel: as many forked workers as fit in nthreads, splitting the threads between them
# (dada/mergePairs on one sample use only a fraction of 128 threads, so parallelizing across samples is what scales)
mcrun = function(X, FUN) {
  workers = max(1, min(length(X), nthreads))
  res = mclapply(X, FUN, threads = max(1, nthreads %/% workers), mc.cores = workers, mc.preschedule = FALSE)
  failed = sapply(res, inherits, "try-error")
  if(any(failed)) stop("parallel job failed for: ", paste(names(X)[failed], collapse = ", "), "\n", res[failed][[1]])
  res
}
script_version = "2.18"   # bump on every modification (2.0, 2.1, ...)
cat("dada2-asv.R version", script_version, "\n")

# ---- locus settings (config params) ----
# parameters follow the DADA2 tutorial (https://benjjneb.github.io/dada2/tutorial.html) unless noted in the config.
# run_cutadapt: remove primers / read-through with cutadapt (variable-length amplicons such as ITS1, where reads run
#   into the opposite primer); otherwise primers must already be gone and truncLen handles read-through.
#   A primer at the read start is removed only when complete (anchored); read-through at the 3' end is removed.
# max_read_length: with cutadapt, shorten every read to at most this many bases after the read-through is removed
#   (cutadapt --length; default: no limit). Reads without read-through (amplicon longer than the read) keep the full
#   cycle count, which can differ between runs (2 x 300 vs 2 x 301), so the same long-ITS1 sequence would give
#   ASVs of different lengths (286 vs 287 bp after trimLeft 14) that never match across runs. A common maximum makes
#   them identical; reads already shorter (trimmed read-through) are unchanged.
# merge: denoise R1 and R2 and merge the pairs; FALSE = R1 only (maxEE[2] is used only when merging).
# asv_lengths: c(min, max) ASV length window ("cut a band"); null keeps all lengths.
# trimLeft: bases removed from the start of R1 (and R2) in filterAndTrim, after cutadapt (default 0), e.g. to drop
#   early cycles with N calls that maxN = 0 would otherwise discard whole reads for.
# runs: per-sequencing-run overrides of truncLen / maxEE / trimLeft, keyed by the run names in the manifest's
#   samples.run_column (see "samples"). With a run column, each run is filtered, error-modelled, denoised and merged
#   on its own; the per-run tables are combined before the length window, chimera removal and taxonomy.
#   Merged pairs span the whole amplicon whatever the truncation, so per-run truncLen keeps ASVs comparable; in R1-only
#   mode the ASV ends at truncLen, so all runs must then share it.
p = cfg$params
params = list(fwd_primer = p$fwd_primer, rev_primer = p$rev_primer, run_cutadapt = isTRUE(p$run_cutadapt),
              merge = isTRUE(p$merge), truncLen = p$truncLen, maxEE = p$maxEE, minLen = p$minLen,
              max_read_length = p$max_read_length,
              minOverlap = p$minOverlap, trimLeft = rep_len(if(is.null(p$trimLeft)) 0 else p$trimLeft, 2),
              asv_lengths = if(length(p$asv_lengths) == 2) p$asv_lengths[1]:p$asv_lengths[2] else NULL,
              taxdb = p$taxonomy$db, taxLevels = p$taxonomy$levels, speciesdb = p$taxonomy$species_db,
              tryRC = isTRUE(p$taxonomy$try_rc), runs = p$runs)
# filtering settings of one run: the locus settings with that run's overrides
run_params = function(r) {
  o = params$runs[[r]]
  list(truncLen = if(!is.null(o$truncLen)) o$truncLen else params$truncLen,
       maxEE = if(!is.null(o$maxEE)) o$maxEE else params$maxEE,
       trimLeft = rep_len(if(!is.null(o$trimLeft)) o$trimLeft else params$trimLeft, 2))
}
cat("threads:", nthreads, " nsamples:", ifelse(is.na(nsamples), "all", nsamples), " taxonomy:", run_taxonomy,
    " resume:", ifelse(resume, resume_mode, FALSE), " reads:", ifelse(params$merge, "R1+R2 merged", "R1 only"), "\n")

rawdir = cfg$paths$raw
workdir = cfg$paths$work
outdir = cfg$paths$results
if(!is.na(nsamples)) {
  workdir = file.path(workdir, paste0("test-n", nsamples))
  outdir = file.path(outdir, paste0("test-n", nsamples))
}
dir.create(workdir, showWarnings = F, recursive = T)
dir.create(outdir, showWarnings = F, recursive = T)

# fail early on missing tools / reference files rather than hours into the run
if(!resume && params$run_cutadapt && cutadapt == "") stop("cutadapt not found on PATH")
needs_taxonomy = run_taxonomy && !identical(resume_mode, "outputs")
if(needs_taxonomy && !file.exists(params$taxdb)) stop("taxonomy reference not found: ", params$taxdb)
if(needs_taxonomy && !is.null(params$speciesdb) && !file.exists(params$speciesdb)) {
  cat("species reference not found, skipping addSpecies:", params$speciesdb, "\n")
  params$speciesdb = NULL
}

# ---- samples ----
# sample names: the metadata id column, optionally renamed (samples.rename: e.g. 16S library names -> ITS names)
# read as text, then convert columns to numbers / logicals, except codes with leading zeros (e.g. SamplingCode 058),
# which stay text; no quote or comment characters, so values may contain ' and #
metadata = read.table(cfg$paths$metadata, sep = "\t", header = T, check.names = F, colClasses = "character",
                      quote = "", comment.char = "", na.strings = "")
metadata[] = lapply(metadata, function(x) if(any(grepl("^0[0-9]", x))) x else type.convert(x, as.is = TRUE))

# ---- site tables (paths.site_tables) joined to the sample metadata ----
# samples.site_key: {metadata: <metadata column>, tables: <table column>, pattern: <regex>}: the first group of the
# pattern applied to both columns is the join key (e.g. "^([0-9]{3})": 058 for 058-ER18 / 058-ER22 / 058-ER18_A), so
# every sample of a site gets the site's row. Tables are .csv (comma) or tab-separated, read as text and converted
# like the metadata. Per table column, in order:
#   - same name as a metadata column (also when equal after make.names, e.g. "Sampling.Area" = "Sampling Area"):
#     fills only that column's empty values; stops if a non-empty value differs
#   - otherwise added under its own name; when its values equal a column already there (e.g. Site_Name and
#     "Site Name"), the report says so
# <key>.metadata_columns.txt lists every sample-data column with its source file and what was done.
metadata_sources = data.frame(column = names(metadata), source = basename(cfg$paths$metadata), action = "metadata",
                              stringsAsFactors = FALSE)
site_tables = unlist(cfg$paths$site_tables)
if(length(site_tables)) {
  sk = cfg$samples$site_key
  if(is.null(sk$metadata) || is.null(sk$tables)) stop("paths.site_tables needs samples.site_key: {metadata, tables, pattern}")
  pat = if(is.null(sk$pattern)) "^(.*)$" else sk$pattern
  keyof = function(x) ifelse(grepl(pat, x), sub(paste0(pat, ".*$"), "\\1", x), NA)
  same = function(a, b) { a = as.character(a); b = as.character(b); identical(is.na(a), is.na(b)) && all(a == b, na.rm = TRUE) }
  mkey = keyof(as.character(metadata[[sk$metadata]]))
  for(f in site_tables) {
    if(!file.exists(f)) stop("site table not found: ", f)
    tb = read.table(f, sep = if(grepl("\\.csv$", f, ignore.case = TRUE)) "," else "\t", header = TRUE, check.names = FALSE,
                    colClasses = "character", quote = "\"", comment.char = "", na.strings = c("", "NA"))
    tb[] = lapply(tb, function(x) if(any(grepl("^0[0-9]", x))) x else type.convert(x, as.is = TRUE))
    tkey = keyof(as.character(tb[[sk$tables]]))
    if(anyNA(tkey) || anyDuplicated(tkey)) stop("site table ", basename(f), ": ", sk$tables, " values must give unique keys with ", pat)
    i = match(mkey, tkey)
    cat(sprintf("site table %s: %d rows, %d columns; %d of %d samples matched\n", basename(f), nrow(tb), ncol(tb), sum(!is.na(i)), nrow(metadata)))
    for(col in setdiff(names(tb), sk$tables)) {
      v = tb[[col]][i]
      hit = names(metadata)[names(metadata) == col | make.names(names(metadata)) == make.names(col)]
      if(length(hit)) {
        old = metadata[[hit[1]]]
        clash = !is.na(old) & !is.na(v) & as.character(old) != as.character(v)
        if(is.numeric(old) && is.numeric(v)) clash = !is.na(old) & !is.na(v) & abs(old - v) > 1e-9 * pmax(1, abs(old))
        if(any(clash)) stop("site table ", basename(f), " column ", col, " differs from metadata column ", hit[1], " for ",
                            sum(clash), " samples, e.g. ", metadata[[1]][which(clash)[1]])
        fill = is.na(old) & !is.na(v)
        if(any(fill)) metadata[[hit[1]]][fill] = v[fill]
        metadata_sources = rbind(metadata_sources, data.frame(column = hit[1], source = basename(f),
          action = sprintf("same values as %s; filled %d empty", hit[1], sum(fill))))
        next
      }
      dup = names(metadata)[vapply(metadata, same, TRUE, b = v)]
      metadata[[col]] = v
      metadata_sources = rbind(metadata_sources, data.frame(column = col, source = basename(f),
        action = if(length(dup)) paste("added; same values as", dup[1]) else "added"))
    }
  }
  acts = table(sub(":.*|;.*", "", metadata_sources$action[metadata_sources$action != "metadata"]))
  cat("site tables:", paste(names(acts), acts, sep = " ", collapse = ", "), "->", ncol(metadata), "metadata columns\n")
}
sample_ids = function(md) {
  ids = as.character(md[[cfg$samples$id_column]])
  if(!is.null(cfg$samples$rename)) ids = sub(cfg$samples$rename$pattern, cfg$samples$rename$replacement, ids)
  ids
}
# fastq files, two ways:
#   samples.manifest: tab-separated file with a header, one line per sample, columns R1 and R2 (paths; relative ones
#     are relative to the manifest's folder) and optionally the sample name (column samples.id_column or "sample";
#     without one, the R1 file name minus _L001_R1_001.fastq.gz). The manifest sets which samples are run; each must
#     be in the metadata. samples.run_column (optional) names a manifest column with the sequencing run of each
#     sample: runs are then processed separately (see "runs" above).
#   otherwise: every metadata sample, paths.raw/<samples.r1|r2> with {sample} filled in (one run)
fastq_suffix = "(_L[0-9]{3})?_R1(_[0-9]{3})?[.](fastq|fq)([.]gz)?$"
read_manifest = function(path) {
  m = read.table(path, sep = "\t", header = T, quote = "", comment.char = "", colClasses = "character",
                 check.names = F, strip.white = T)
  col = function(names) { hit = which(tolower(colnames(m)) %in% tolower(names)); if(length(hit)) colnames(m)[hit[1]] else NA }
  r1 = col("R1"); r2 = col("R2"); id = col(c(cfg$samples$id_column, "sample"))
  if(is.na(r1) || is.na(r2)) stop("manifest needs columns R1 and R2: ", path, " (columns: ", paste(colnames(m), collapse = ", "), ")")
  m = m[m[[r1]] != "" | m[[r2]] != "", , drop = F]
  if(any(m[[r1]] == "" | m[[r2]] == "")) stop("manifest lines with only one of R1 / R2: ", path)
  abs_path = function(p) ifelse(grepl("^(/|~)", p), path.expand(p), file.path(dirname(path), p))
  ids = if(!is.na(id)) m[[id]] else sub(fastq_suffix, "", basename(m[[r1]]))
  if(anyDuplicated(ids)) stop("duplicated sample names in manifest ", path, ": ", paste(unique(ids[duplicated(ids)]), collapse = ", "))
  runcol = cfg$samples$run_column
  if(!is.null(runcol) && !(runcol %in% colnames(m))) stop("manifest has no run column ", runcol, ": ", path)
  runs = if(!is.null(runcol)) m[[runcol]] else rep("all", nrow(m))
  if(any(runs == "")) stop("manifest lines with an empty ", runcol, ": ", path)
  list(samples = ids, r1 = abs_path(m[[r1]]), r2 = abs_path(m[[r2]]), runs = runs)
}
if(!is.null(cfg$samples$manifest)) {
  man = read_manifest(cfg$samples$manifest)
  nometa = setdiff(man$samples, sample_ids(metadata))
  if(length(nometa) > 0) stop(length(nometa), " manifest samples not in the metadata (", cfg$samples$id_column, "): ",
                              paste(head(nometa, 10), collapse = ", "))
  samplenames = man$samples; Ffiles = man$r1; Rfiles = man$r2; sampleruns = man$runs
  cat("fastq manifest:", cfg$samples$manifest, "\n")
} else {
  samplenames = sample_ids(metadata)
  fastq = function(pattern) file.path(rawdir, vapply(samplenames, function(s) gsub("{sample}", s, pattern, fixed = TRUE), ""))
  Ffiles = fastq(cfg$samples$r1)
  Rfiles = fastq(cfg$samples$r2)
  sampleruns = rep("all", length(samplenames))
}
if(!is.na(nsamples)) {
  # test runs: the first N samples of each run, so every run (and its settings) is exercised
  pick = unlist(lapply(split(seq_along(samplenames), factor(sampleruns, unique(sampleruns))), head, nsamples))
  samplenames = samplenames[pick]; Ffiles = Ffiles[pick]; Rfiles = Rfiles[pick]; sampleruns = sampleruns[pick]
}
missing = c(Ffiles, Rfiles)[!file.exists(c(Ffiles, Rfiles))]
if(!resume && length(missing) > 0) stop("missing fastq files:\n", paste(missing, collapse = "\n"))
names(Ffiles) = names(Rfiles) = names(sampleruns) = samplenames
runnames = unique(sampleruns)
unknown = setdiff(names(params$runs), runnames)
if(length(unknown) > 0) cat("note: run settings for runs not in this sample set:", unknown, "\n")
if(!params$merge && length(unique(lapply(runnames, function(r) run_params(r)$truncLen[1]))) > 1)
  stop("R1-only mode (merge: false) needs the same truncLen in all runs (the ASVs end at truncLen)")
cat(locus, run, ":", length(samplenames), "samples in", length(runnames), "run(s)\n")
for(r in runnames) {
  rp = run_params(r)
  cat(sprintf("  %-32s %4d samples  truncLen %s  maxEE %s  trimLeft %s\n", r, sum(sampleruns == r),
              paste(rp$truncLen, collapse = "/"), paste(rp$maxEE, collapse = "/"), paste(rp$trimLeft, collapse = "/")))
}

# ---- cutadapt tracking ----
# per-sample numbers from the cutadapt reports (<work>/cutadapt/<sample>.cutadapt.log): read pairs in, R1 / R2 reads
# with a primer or read-through removed, pairs dropped as too short, pairs written, bases in / out per read, and
# how often each primer was found. Written to <key>.cutadapt_summary.txt; readtracking.txt gets raw (pairs given
# to cutadapt) and cutadapt_R1_trimmed / cutadapt_R2_trimmed before input (= pairs written by cutadapt).
parse_cutadapt_logs = function(dir, samples) {
  num = function(pat, txt) { m = regmatches(txt, regexec(pat, txt))[[1]]; if(length(m) < 2) NA_real_ else as.numeric(gsub(",", "", m[2])) }
  rows = lapply(samples, function(s) {
    f = file.path(dir, paste0(s, ".cutadapt.log"))
    if(!file.exists(f)) return(NULL)
    txt = paste(readLines(f, warn = FALSE), collapse = "\n")
    ad = regmatches(txt, gregexpr("Sequence: [ACGTNRYSWKMBDHVacgtn]+; Type: [^;]+;[^\n]*Trimmed: [0-9]+ times", txt))[[1]]
    adn = as.numeric(sub(".*Trimmed: ([0-9]+) times.*", "\\1", ad))
    adt = sub(".*Type: ([^;]+);.*", "\\1", ad)
    data.frame(sample = s,
      pairs_in = num("Total read pairs processed:\\s+([0-9,]+)", txt),
      R1_trimmed = num("Read 1 with adapter:\\s+([0-9,]+)", txt), R2_trimmed = num("Read 2 with adapter:\\s+([0-9,]+)", txt),
      pairs_too_short = num("Pairs that were too short:\\s+([0-9,]+)", txt),
      pairs_written = num("Pairs written \\(passing filters\\):\\s+([0-9,]+)", txt),
      R1_bp_in = num("Total basepairs processed:[^\n]*\n\\s+Read 1:\\s+([0-9,]+)", txt),
      R1_bp_out = num("Total written \\(filtered\\):[^\n]*\n\\s+Read 1:\\s+([0-9,]+)", txt),
      R2_bp_in = num("Total basepairs processed:[^\n]*\n[^\n]*\n\\s+Read 2:\\s+([0-9,]+)", txt),
      R2_bp_out = num("Total written \\(filtered\\):[^\n]*\n[^\n]*\n\\s+Read 2:\\s+([0-9,]+)", txt),
      # adapters in command-line order: -g (R1 5'), -a (R1 3'), -G (R2 5'), -A (R2 3')
      R1_5prime_primer = if(length(adn) >= 1) adn[1] else NA, R1_3prime_readthrough = if(length(adn) >= 2) adn[2] else NA,
      R2_5prime_primer = if(length(adn) >= 3) adn[3] else NA, R2_3prime_readthrough = if(length(adn) >= 4) adn[4] else NA,
      stringsAsFactors = FALSE)
  })
  rows = do.call(rbind, rows)
  if(is.null(rows)) return(NULL)
  rows$pct_R1_trimmed = round(100 * rows$R1_trimmed / rows$pairs_in, 2)
  rows$pct_R1_bp_kept = round(100 * rows$R1_bp_out / rows$R1_bp_in, 2)
  rows
}
add_cutadapt_tracking = function(track, cs) {
  if(is.null(cs) || "raw" %in% colnames(track)) return(track)
  i = match(rownames(track), cs$sample)
  cbind(raw = cs$pairs_in[i], cutadapt_R1_trimmed = cs$R1_trimmed[i], cutadapt_R2_trimmed = cs$R2_trimmed[i], track)
}
write_cutadapt_summary = function(cs) {
  if(is.null(cs)) return(invisible())
  cs$run = sampleruns[cs$sample]
  write.table(cs[, c("sample", "run", setdiff(names(cs), c("sample", "run")))],
              file.path(outdir, paste0(key, ".cutadapt_summary.txt")), sep = "\t", quote = F, row.names = F)
  for(r in unique(cs$run)) { z = cs[cs$run == r, ]
    cat(sprintf("cutadapt %-28s %4d samples: %s pairs in, R1 trimmed %.1f%%, R2 trimmed %.1f%%, too short %d, written %s, R1 bases kept %.1f%%\n",
                r, nrow(z), format(sum(z$pairs_in), big.mark = ","), 100 * sum(z$R1_trimmed) / sum(z$pairs_in),
                100 * sum(z$R2_trimmed) / sum(z$pairs_in), as.integer(sum(z$pairs_too_short)),
                format(sum(z$pairs_written), big.mark = ","), 100 * sum(z$R1_bp_out) / sum(z$R1_bp_in))) }
}
write_track = function(track) {
  tw = tibble::rownames_to_column(as.data.frame(track), "sample")
  if(length(runnames) > 1) tw = cbind(tw[1], run = sampleruns[tw$sample], tw[-1])
  write.table(tw, file = file.path(outdir, paste0(key, ".readtracking.txt")), sep = "\t", row.names = F, col.names = T, quote = F)
}
cutstats = NULL

checkpoint = file.path(outdir, paste0(key, ".dada2.seqtab.RData"))
if(resume) {
  # ---- resume from the seqtab checkpoint (or, for --resume=outputs, from the final RData with taxonomy) ----
  if(resume_mode == "outputs") checkpoint = file.path(outdir, paste0(key, ".dada2.RData"))
  if(!file.exists(checkpoint)) stop("no checkpoint to resume from: ", checkpoint)
  current_params = params   # keep this run's reference paths, not the saved ones
  load(checkpoint)
  params = current_params
  cat("resumed from", checkpoint, ":", ncol(seqtab.nochim), "ASVs in", nrow(seqtab.nochim), "samples\n")
} else {
  # ---- raw quality profiles (a few samples) ----
  qc.samples = sample(samplenames, min(6, length(samplenames)))
  # PNG (300 dpi); PDF too with --pdf
  for(rd in c("R1", "R2")) {
    qp = plotQualityProfile(if(rd == "R1") Ffiles[qc.samples] else Rfiles[qc.samples])
    f = file.path(outdir, paste0(key, ".quality.raw.", rd))
    ggsave(paste0(f, ".png"), qp, width = 12, height = 8, dpi = 300)
    if(write_pdf) ggsave(paste0(f, ".pdf"), qp, width = 12, height = 8)
  }

  # ---- primer check ----
  allOrients = function(primer) {
    dna = DNAString(primer)
    orients = c(Forward = dna, Complement = Biostrings::complement(dna),
                Reverse = Biostrings::reverse(dna), RevComp = reverseComplement(dna))
    sapply(orients, toString)
  }
  primerHits = function(primer, fn) {
    nhits = vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    sum(nhits > 0)
  }
  primerCheck = function(Ffile, Rfile) {
    rbind(FWD.R1 = sapply(allOrients(params$fwd_primer), primerHits, fn = Ffile),
          FWD.R2 = sapply(allOrients(params$fwd_primer), primerHits, fn = Rfile),
          REV.R1 = sapply(allOrients(params$rev_primer), primerHits, fn = Ffile),
          REV.R2 = sapply(allOrients(params$rev_primer), primerHits, fn = Rfile))
  }
  # filter Ns first, otherwise vcountPattern with IUPAC codes is unreliable
  prefiltF = file.path(workdir, "prefilt", basename(Ffiles))
  prefiltR = file.path(workdir, "prefilt", basename(Rfiles))
  filterAndTrim(Ffiles[qc.samples[1]], prefiltF[1], Rfiles[qc.samples[1]], prefiltR[1], maxN = 0, multithread = nthreads)
  cat("primer hits before trimming (", qc.samples[1], "):\n"); print(primerCheck(prefiltF[1], prefiltR[1]))

  # ---- primer / read-through removal (ITS) ----
  if(params$run_cutadapt) {
    cutdir = file.path(workdir, "cutadapt")
    dir.create(cutdir, showWarnings = F)
    cutFfiles = file.path(cutdir, basename(Ffiles))
    cutRfiles = file.path(cutdir, basename(Rfiles))
    fwd.rc = dada2:::rc(params$fwd_primer)
    rev.rc = dada2:::rc(params$rev_primer)
    invisible(mcrun(setNames(seq_along(samplenames), samplenames), function(i, threads) {
      # R1: fwd primer anchored at the 5' end, rc(rev) read-through at 3'; R2: rev primer anchored at 5', rc(fwd) at 3'.
      # Anchored (^primer): only a complete primer at the read start is removed. Unanchored, cutadapt also took 3-6
      # bases off the start of 0.1-0.3% of NEON ITS reads by chance (their primers are never sequenced). The 3'
      # read-through keeps cutadapt's default minimum overlap (3): short matches there are mostly real read-through.
      system2(cutadapt, args = c("-g", paste0("^", params$fwd_primer), "-a", rev.rc,
                                 "-G", paste0("^", params$rev_primer), "-A", fwd.rc,
                                 "-n", 2, "--minimum-length", params$minLen, "--cores", threads,
                                 # same maximum read length in every run (see max_read_length above); applied after
                                 # adapter trimming, to R1 and R2
                                 if(!is.null(params$max_read_length)) c("--length", params$max_read_length),
                                 "-o", cutFfiles[i], "-p", cutRfiles[i], Ffiles[i], Rfiles[i]),
              stdout = file.path(cutdir, paste0(samplenames[i], ".cutadapt.log")))
    }))
    names(cutFfiles) = names(cutRfiles) = samplenames
    cutstats = parse_cutadapt_logs(cutdir, samplenames)
    write_cutadapt_summary(cutstats)
    Ffiles = cutFfiles
    Rfiles = cutRfiles
    cat("primer hits after cutadapt (", qc.samples[1], "):\n"); print(primerCheck(Ffiles[qc.samples[1]], Rfiles[qc.samples[1]]))
  }

  # ---- quality filtering ----
  filtFfiles = file.path(workdir, "filtered", paste0(samplenames, "_R1.filt.fastq.gz"))
  filtRfiles = file.path(workdir, "filtered", paste0(samplenames, "_R2.filt.fastq.gz"))
  names(filtFfiles) = names(filtRfiles) = samplenames
  # per run: its own truncLen / maxEE / trimLeft (run_params)
  out = do.call(rbind, lapply(runnames, function(r) {
    s = samplenames[sampleruns == r]; rp = run_params(r)
    o = if(params$merge) {
      filterAndTrim(Ffiles[s], filtFfiles[s], Rfiles[s], filtRfiles[s],
                    truncLen = rp$truncLen, trimLeft = rp$trimLeft, maxN = 0, maxEE = rp$maxEE, truncQ = 2,
                    minLen = params$minLen, rm.phix = TRUE, compress = TRUE, multithread = nthreads, verbose = TRUE)
    } else {
      filterAndTrim(Ffiles[s], filtFfiles[s],
                    truncLen = rp$truncLen[1], trimLeft = rp$trimLeft[1], maxN = 0, maxEE = rp$maxEE[1], truncQ = 2,
                    minLen = params$minLen, rm.phix = TRUE, compress = TRUE, multithread = nthreads, verbose = TRUE)
    }
    rownames(o) = s
    o
  }))
  out = out[samplenames, , drop = F]

  lowcount = rownames(out)[out[, "reads.out"] < mincount]
  if(length(lowcount) > 0) cat("dropping", length(lowcount), "low-count samples:", lowcount, "\n")
  keep = samplenames[!(samplenames %in% lowcount)]
  filtFfiles = filtFfiles[keep]
  filtRfiles = filtRfiles[keep]
  keepruns = runnames[runnames %in% sampleruns[keep]]

  # ---- error models: one per run (and read) ----
  # all run x read models are independent: learn them concurrently, splitting the threads
  jobs = list()
  for(r in keepruns) {
    s = keep[sampleruns[keep] == r]
    jobs[[paste(r, "R1", sep = "\t")]] = filtFfiles[s]
    if(params$merge) jobs[[paste(r, "R2", sep = "\t")]] = filtRfiles[s]
  }
  errs = mcrun(jobs, function(fns, threads) { set.seed(100); learnErrors(fns, multithread = threads, randomize = TRUE) })
  errF = setNames(lapply(keepruns, function(r) errs[[paste(r, "R1", sep = "\t")]]), keepruns)
  errR = if(params$merge) setNames(lapply(keepruns, function(r) errs[[paste(r, "R2", sep = "\t")]]), keepruns) else NULL
  for(r in keepruns) {
    tag = if(length(runnames) > 1) paste0(".", r) else ""
    # PNG (300 dpi); PDF too with --pdf
    for(rd in if(params$merge) c("R1", "R2") else "R1") {
      ep = plotErrors(if(rd == "R1") errF[[r]] else errR[[r]], nominalQ = TRUE)
      f = file.path(outdir, paste0(key, tag, ".errors.", rd))
      ggsave(paste0(f, ".png"), ep, width = 10, height = 8, dpi = 300)
      if(write_pdf) ggsave(paste0(f, ".pdf"), ep, width = 10, height = 8)
    }
  }
  save(out, keep, sampleruns, errF, errR, params, file = file.path(outdir, paste0(key, ".dada2.errors.RData")))

  # ---- denoise, merge, chimeras ----
  # pool = "pseudo" done by hand so that samples (R1 and R2) run in parallel; same algorithm as dada():
  #   pass 1 per sample, ASVs found in >= PSEUDO_PREVALENCE samples (or >= PSEUDO_ABUNDANCE reads) become priors,
  #   pass 2 per sample with those priors
  dadaPseudo = function(fns, err) {
    pass = function(priors) mcrun(as.list(fns), function(fn, threads)
      dada(fn, err = err, priors = priors, multithread = threads, verbose = FALSE))
    first = pass(character(0))
    st = makeSequenceTable(first)
    priors = colnames(st)[colSums(st > 0) >= getDadaOpt("PSEUDO_PREVALENCE") | colSums(st) >= getDadaOpt("PSEUDO_ABUNDANCE")]
    cat(length(priors), "pseudo-pooling priors\n")
    pass(priors)
  }
  # each run with its own error model, pseudo-pooled within the run
  dadaRuns = function(fns, errs) {
    d = do.call(c, lapply(keepruns, function(r) { cat(r, ": "); dadaPseudo(fns[sampleruns[names(fns)] == r], errs[[r]]) }))
    d[names(fns)]
  }
  dadaFs = dadaRuns(filtFfiles, errF)
  if(params$merge) {
    dadaRs = dadaRuns(filtRfiles, errR)
    # trimOverhang: drop any bases where one read extends past the start of the other (primer read-through
    # in amplicons shorter than truncLen). mergePairs has no multithread option -> one sample per worker
    mergers = mcrun(setNames(as.list(keep), keep), function(s, threads)
      mergePairs(dadaFs[[s]], filtFfiles[[s]], dadaRs[[s]], filtRfiles[[s]], minOverlap = params$minOverlap,
                 trimOverhang = TRUE, verbose = TRUE))
    seqtab = makeSequenceTable(mergers)
  } else {
    dadaRs = mergers = NULL
    seqtab = makeSequenceTable(dadaFs)
  }
  cat("ASV length distribution:\n"); print(table(nchar(getSequences(seqtab))))
  if(!is.null(params$asv_lengths)) {
    inrange = nchar(getSequences(seqtab)) %in% params$asv_lengths
    cat("keeping", sum(inrange), "of", ncol(seqtab), "ASVs (", round(100 * sum(seqtab[, inrange]) / sum(seqtab), 1),
        "% of reads ) with length in", range(params$asv_lengths), "\n")
    seqtab = seqtab[, inrange, drop = F]
  }

  seqtab.nochim = removeBimeraDenovo(seqtab, method = "consensus", multithread = nthreads, verbose = TRUE)
  cat("fraction of reads kept after chimera removal:", sum(seqtab.nochim) / sum(seqtab), "\n")

  # ---- read tracking ----
  getN = function(x) sum(getUniques(x))
  track = cbind(out[keep, , drop = F], denoisedF = sapply(dadaFs, getN))
  if(params$merge) track = cbind(track, denoisedR = sapply(dadaRs, getN), merged = sapply(mergers, getN))
  track = cbind(track, nonchim = rowSums(seqtab.nochim))
  colnames(track)[1:2] = c("input", "filtered")
  track = add_cutadapt_tracking(track, cutstats)
  write_track(track)
  # checkpoint before taxonomy: assignTaxonomy against SILVA is memory hungry and can get the job killed
  save(out, keep, sampleruns, errF, errR, dadaFs, dadaRs, mergers, seqtab, seqtab.nochim, track, params, cutstats,
       file = checkpoint)
  if(!run_taxonomy) {
    cat("stopping before taxonomy (--no-taxonomy):", ncol(seqtab.nochim), "ASVs in", nrow(seqtab.nochim), "samples ->", outdir, "\n")
    quit(save = "no", status = 0)
  }
}

# ---- taxonomy ----
if(!(resume && resume_mode == "outputs")) {
  taxa = assignTaxonomy(seqtab.nochim, params$taxdb, taxLevels = params$taxLevels,
                        multithread = nthreads, tryRC = params$tryRC)
  if(!is.null(params$speciesdb)) taxa = addSpecies(taxa, params$speciesdb)
}

# ---- write outputs (ASV sequence as identifier) ----
asvtable.towrite = tibble::rownames_to_column(as.data.frame(t(seqtab.nochim)), "asv_seq")
write.table(asvtable.towrite, file = file.path(outdir, paste0(key, ".ASVtable.txt")),
            sep = "\t", row.names = F, col.names = T, quote = F)

taxa.towrite = tibble::rownames_to_column(as.data.frame(taxa), "asv_seq")
write.table(taxa.towrite, file = file.path(outdir, paste0(key, ".taxonomy.txt")),
            sep = "\t", row.names = F, col.names = T, quote = F)

asvs = DNAStringSet(colnames(seqtab.nochim))
names(asvs) = colnames(seqtab.nochim)
writeXStringSet(asvs, file.path(outdir, paste0(key, ".ASV.fa")))

# cutadapt tracking for runs resumed from a checkpoint made before it existed: read the reports if still in <work>
if(params$run_cutadapt && is.null(cutstats)) {
  cutstats = parse_cutadapt_logs(file.path(workdir, "cutadapt"), rownames(out))
  if(is.null(cutstats)) cat("note: no cutadapt reports in", file.path(workdir, "cutadapt"), "- readtracking without cutadapt columns\n")
  write_cutadapt_summary(cutstats)
}
if(!"raw" %in% colnames(track) && !is.null(cutstats)) { track = add_cutadapt_tracking(track, cutstats); write_track(track) }

save(out, keep, errF, errR, dadaFs, dadaRs, mergers, seqtab, seqtab.nochim, taxa, track, params, cutstats,
     file = file.path(outdir, paste0(key, ".dada2.RData")))

# sample-data columns and where they come from (metadata file, site tables)
write.table(metadata_sources, file.path(outdir, paste0(key, ".metadata_columns.txt")), sep = "\t", quote = F, row.names = F)

# ---- phyloseq object ----
# sample_data = metadata (+ read tracking as reads_* columns); refseq holds the ASV sequences, taxa renamed ASV1..N
# (seqtab columns are in decreasing abundance). The tree is added if <key>.ASV.tree exists in outdir or in
# outdir/v4align (where the tree step, align-asv-v4.sh, writes it); its tip labels are the ASV sequences, as in ASV.fa
sampledata = metadata
rownames(sampledata) = sample_ids(metadata)
sampledata = sampledata[rownames(seqtab.nochim), , drop = F]
sampledata = cbind(sampledata, setNames(as.data.frame(track[rownames(sampledata), , drop = F]), paste0("reads_", colnames(track))))
ps = phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(taxa), sample_data(sampledata), refseq(asvs))
treefile = file.path(c(outdir, file.path(outdir, "v4align")), paste0(key, ".ASV.tree"))
treefile = treefile[file.exists(treefile)][1]
if(!is.na(treefile)) {
  # merge_phyloseq keeps only ASVs that are tips of the tree (align-asv-v4.sh drops ASVs not spanning V4)
  ntaxa_before = ntaxa(ps); nreads_before = sum(taxa_sums(ps))
  ps = merge_phyloseq(ps, phy_tree(ape::read.tree(treefile)))
  cat("added tree:", treefile, "-", ntaxa_before - ntaxa(ps), "ASVs not in the tree dropped (",
      round(100 * (1 - sum(taxa_sums(ps)) / nreads_before), 2), "% of reads )\n")
}
taxa_names(ps) = paste0("ASV", seq_len(ntaxa(ps)))
# metadata column names are kept as they are (e.g. "Soil pH"): phyloseq rebuilds sample_data with data.frame() in
# many operations ($<-, prune_samples, merge_phyloseq), which turns them into "Soil.pH". sample_df() reads the table
# back with the original names, and with_metadata_names() puts them back before an object is saved.
md_names = c(colnames(sampledata))
sample_df = function(x) { s = sample_data(x); d = as.data.frame(s@.Data, col.names = seq_along(s@names))
  names(d) = c(md_names, s@names[-seq_along(md_names)]); rownames(d) = sample_names(x); d }
with_metadata_names = function(x) { sample_data(x) = sample_data(sample_df(x)); x }
ps = with_metadata_names(ps)
saveRDS(ps, file.path(outdir, paste0(key, ".phyloseq.rds")))
print(ps)
# filtered copy: keeps only the config's keep_taxa (rank: names, e.g. Kingdom: Fungi; unassigned at that rank is
# dropped too), drops remove_taxa (rank: names, e.g. Order: Chloroplast, Family: Mitochondria) and ASVs unassigned at
# any of require_ranks (e.g. Phylum). Names are matched without UNITE's rank prefixes (k__Fungi = Fungi).
# Every removed ASV gets a reason (the first that applies, in that order); <key>.taxa_filter.txt summarizes the
# removed ASVs by reason and Kingdom / Phylum / Class (ASVs, reads, % of all reads, samples), plus what is kept.
psf = ps
tt = as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
tt[] = lapply(tt, function(x) sub("^[a-z]__", "", x))
reason = rep(NA_character_, ntaxa(ps))
for(rank in names(p$keep_taxa)) {
  out_k = !(tt[[rank]] %in% p$keep_taxa[[rank]])
  reason[is.na(reason) & out_k] = paste0("not ", rank, " ", paste(p$keep_taxa[[rank]], collapse = "/"))
}
for(rank in names(p$remove_taxa)) {
  hit = tt[[rank]] %in% p$remove_taxa[[rank]]
  reason[is.na(reason) & hit] = paste0(rank, " ", tt[[rank]][is.na(reason) & hit])
}
for(rank in p$require_ranks) reason[is.na(reason) & is.na(tt[[rank]])] = paste0("unassigned at ", rank)
drop = !is.na(reason)
reads = taxa_sums(ps); total = sum(reads)
otu = as(otu_table(ps), "matrix")                                  # samples x taxa (taxa_are_rows = FALSE)
nsamp = colSums(otu > 0)
lvl = function(x) ifelse(is.na(x), "unassigned", x)
grp = data.frame(status = ifelse(drop, "removed", "kept"), reason = ifelse(drop, reason, ""),
                 Kingdom = lvl(tt[[1]]), Phylum = if(ncol(tt) > 1) lvl(tt[[2]]) else "",
                 Class = if(ncol(tt) > 2) lvl(tt[[3]]) else "", stringsAsFactors = FALSE)
rep_tab = do.call(rbind, lapply(split(seq_len(ntaxa(ps)), grp, drop = TRUE), function(i) cbind(grp[i[1], ],
  ASVs = length(i), reads = sum(reads[i]), pct_reads = round(100 * sum(reads[i]) / total, 3),
  samples = sum(rowSums(otu[, i, drop = F]) > 0), max_samples_one_ASV = max(nsamp[i]))))
rep_tab = rep_tab[order(rep_tab$status != "removed", -rep_tab$reads), ]
tot_row = function(s, sel) data.frame(status = s, reason = "TOTAL", Kingdom = "", Phylum = "", Class = "",
  ASVs = sum(sel), reads = sum(reads[sel]), pct_reads = round(100 * sum(reads[sel]) / total, 3),
  samples = sum(rowSums(otu[, sel, drop = F]) > 0), max_samples_one_ASV = if(any(sel)) max(nsamp[sel]) else 0)
rep_tab = rbind(tot_row("removed", drop), tot_row("kept", !drop), rep_tab)
write.table(rep_tab, file = file.path(outdir, paste0(key, ".taxa_filter.txt")), sep = "\t", row.names = F, quote = F)
if(any(drop)) {
  cat("taxa filter: removing", sum(drop), "of", ntaxa(ps), "ASVs (", round(100 * sum(reads[drop]) / total, 2),
      "% of reads ) ->", paste0(key, ".taxa_filter.txt"), "\n")
  byreason = tapply(reads[drop], reason[drop], sum)
  for(r in names(sort(byreason, decreasing = TRUE)))
    cat(sprintf("  %-40s %6d ASVs  %6.2f%% of reads\n", r, sum(reason[drop] == r), 100 * byreason[[r]] / total))
  psf = prune_taxa(!drop, ps)
}
sd_f = sample_df(psf)
sd_f$reads_filtered_taxa = sample_sums(psf)
sample_data(psf) = sample_data(sd_f)
reads_before_taxa_filter = sample_sums(ps)
# read tracking: reads kept by the taxa filter (what minfinal applies to) and their % of the sample's nonchim reads
# (with a tree, reads of ASVs not in the tree count as not kept). Column taxa_filtered_<label>: params
# taxa_filter_label if set, else from the filter: keep_taxa names (Fungi -> fungal, Bacteria -> bacterial,
# Archaea -> archaeal, others lowercased), else no_<remove_taxa names> and <require_ranks>_assigned, e.g.
# taxa_filtered_fungal, taxa_filtered_no_chloroplast_mitochondria_phylum_assigned; taxa_filtered_<label>_pct_nonchim
filter_label = function(p) {
  if(length(p$taxa_filter_label)) return(p$taxa_filter_label)
  adj = c(Fungi = "fungal", Bacteria = "bacterial", Archaea = "archaeal")
  nm = function(x) { x = unlist(x); ifelse(x %in% names(adj), adj[x], tolower(x)) }
  if(length(p$keep_taxa)) return(paste(nm(p$keep_taxa), collapse = "_"))
  lab = c(if(length(p$remove_taxa)) paste0("no_", paste(nm(p$remove_taxa), collapse = "_")),
          if(length(p$require_ranks)) paste0(paste(tolower(p$require_ranks), collapse = "_"), "_assigned"))
  if(length(lab)) paste(lab, collapse = "_") else "none"
}
tf_col = paste0("taxa_filtered_", gsub("[^A-Za-z0-9]+", "_", filter_label(p)))
track = track[, !grepl("^taxa_filtered", colnames(track)), drop = F]
tf_reads = as.numeric(sample_sums(psf)[rownames(track)])
track = cbind(track, setNames(data.frame(tf_reads, round(100 * tf_reads / track[, "nonchim"], 2)),
                              c(tf_col, paste0(tf_col, "_pct_nonchim"))))
write_track(track)

# ---- failed samples ----
# context for each reason: the sample's reads vs the median (and the lowest passing sample) at the same stage
final = sample_sums(psf)   # reads after the taxa filter (e.g. fungal reads): what minfinal applies to
lowcount = setdiff(rownames(out), keep)
lowfinal = names(final)[final < minfinal]
kept_frac = final / reads_before_taxa_filter[names(final)]
lowfrac = setdiff(names(final)[kept_frac < min_kept], lowfinal)
pass = setdiff(names(final), c(lowfinal, lowfrac))
nextlow = pass[which.min(final[pass])]
med_in = median(out[, "reads.in"])
failed = rbind(
  data.frame(sample = lowcount, stage = rep("quality filtering", length(lowcount)),
             reads_input = out[lowcount, "reads.in"], reads_at_failure = out[lowcount, "reads.out"],
             reason = sprintf("%d of %d input reads passed filterAndTrim (< %d); median input across samples %d (sample has %.1f%% of it)",
                              out[lowcount, "reads.out"], out[lowcount, "reads.in"], mincount, as.integer(med_in),
                              100 * out[lowcount, "reads.in"] / med_in)),
  data.frame(sample = lowfinal, stage = rep("final read count", length(lowfinal)),
             reads_input = out[lowfinal, "reads.in"], reads_at_failure = as.integer(final[lowfinal]),
             reason = sprintf("low number of reads: %d reads after the taxa filter (< %d), %.1f%% of the median (%d); lowest passing sample %s has %d; input %d reads vs median %d",
                              as.integer(final[lowfinal]), minfinal, 100 * final[lowfinal] / median(final), as.integer(median(final)),
                              nextlow, as.integer(final[nextlow]), out[lowfinal, "reads.in"], as.integer(med_in))),
  data.frame(sample = lowfrac, stage = rep("taxa filter fraction", length(lowfrac)),
             reads_input = out[lowfrac, "reads.in"], reads_at_failure = as.integer(final[lowfrac]),
             reason = sprintf("%.1f%% of the sample's reads kept by the taxa filter (%d of %d; < %.0f%%); median across samples %.1f%%",
                              100 * kept_frac[lowfrac], as.integer(final[lowfrac]), as.integer(reads_before_taxa_filter[lowfrac]),
                              100 * min_kept, 100 * median(kept_frac))))
# reads before the taxa filter: shows which samples fail only because of it
failed$reads_before_taxa_filter = as.integer(reads_before_taxa_filter[failed$sample])
only_filter = failed$stage == "final read count" & failed$reads_before_taxa_filter >= minfinal
if(any(only_filter)) cat(sum(only_filter), "sample(s) fail only after the taxa filter:",
                         paste(failed$sample[only_filter], collapse = ", "), "\n")
write.table(failed, file = file.path(outdir, paste0(key, ".failed_samples.txt")), sep = "\t", row.names = F, col.names = T, quote = F)
cat(nrow(failed), "failed samples ->", file.path(outdir, paste0(key, ".failed_samples.txt")), "\n")
if(nrow(failed) > 0) print(failed[, c("sample", "stage", "reads_input", "reads_at_failure")], row.names = F)
psf = with_metadata_names(prune_samples(pass, psf))
saveRDS(psf, file.path(outdir, paste0(key, ".phyloseq.filtered.rds")))
print(psf)

# ---- text exports of the filtered phyloseq object, and the README of their tarball ----
# <key>.filtered.otu_table.txt (ASVs x samples, asv_id = the object's taxa names), .tax_table.txt (asv_id, ranks,
# asv_seq), .sample_data.txt (one row per sample, original column names), .refseq.fa (names = asv_id), .tree (when the
# object has one), and <key>.filtered.README.txt describing them. The package step (package-filtered.py) bundles
# them with <key>.readtracking.txt, metadata_columns.txt, taxa_filter.txt, failed_samples.txt and, for ITS, the guild tables into
# <key>.filtered.tar.gz (README.txt = this README plus a guild section)
ex = file.path(outdir, paste0(key, ".filtered"))
otu_x = as(otu_table(psf), "matrix"); if(!taxa_are_rows(psf)) otu_x = t(otu_x)          # ASVs x samples
write.table(data.frame(asv_id = rownames(otu_x), otu_x, check.names = FALSE), paste0(ex, ".otu_table.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
tt_x = as(tax_table(psf), "matrix")
write.table(data.frame(asv_id = rownames(tt_x), tt_x, asv_seq = as.character(refseq(psf))[rownames(tt_x)], check.names = FALSE),
            paste0(ex, ".tax_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
sd_x = sample_df(psf)
if(!identical(as.character(sd_x[[1]]), rownames(sd_x))) sd_x = data.frame(sample_id = rownames(sd_x), sd_x, check.names = FALSE)
write.table(sd_x, paste0(ex, ".sample_data.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
writeXStringSet(refseq(psf), paste0(ex, ".refseq.fa"))
has_tree = !is.null(phy_tree(psf, errorIfNULL = FALSE))
if(has_tree) ape::write.tree(phy_tree(psf), paste0(ex, ".tree"))
exports = paste0(basename(ex), c(".otu_table.txt", ".tax_table.txt", ".sample_data.txt", ".refseq.fa", if(has_tree) ".tree"))
extra = paste0(key, c(".metadata_columns.txt", ".taxa_filter.txt", ".failed_samples.txt"))
extra = extra[file.exists(file.path(outdir, extra))]
filt = c(if(length(p$keep_taxa)) paste0("keep_taxa ", paste(names(p$keep_taxa), sapply(p$keep_taxa, paste, collapse = "/"), collapse = "; ")),
         if(length(p$remove_taxa)) paste0("remove_taxa ", paste(names(p$remove_taxa), sapply(p$remove_taxa, paste, collapse = "/"), collapse = "; ")),
         if(length(p$require_ranks)) paste0("require_ranks ", paste(p$require_ranks, collapse = ", ")))
readme = c(
  paste0(key, ".filtered: the filtered phyloseq object as text files"),
  paste0("written ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), " by dada2-asv.R ", script_version, " (amplicon-seq-workflow)"),
  "",
  paste0("Content of ", key, ".phyloseq.filtered.rds: ", nsamples(psf), " samples, ", ntaxa(psf), " ASVs, ",
         format(sum(otu_x), big.mark = ","), " reads", if(has_tree) ", with a tree" else ""),
  paste0("Taxa filter: ", if(length(filt)) paste(filt, collapse = "; ") else "none", " (", key, ".taxa_filter.txt)"),
  paste0("Failed samples removed: fewer than ", minfinal, " reads after the taxa filter",
         if(min_kept > 0) paste0(", or less than ", 100 * min_kept, "% of their reads kept by it") else "",
         " (", key, ".failed_samples.txt; ", nrow(failed), " samples)"),
  paste0("Taxonomy: ", basename(params$taxdb)),
  "",
  "Files",
  paste0("  ", basename(ex), ".otu_table.txt     read counts: one row per ASV (asv_id), one column per sample"),
  paste0("  ", basename(ex), ".tax_table.txt     taxonomy: asv_id, ", paste(colnames(tt_x), collapse = ", "), ", asv_seq (the ASV sequence)"),
  paste0("  ", basename(ex), ".sample_data.txt   one row per sample: sample metadata (", ncol(sd_x), " columns; reads_* = read tracking)"),
  paste0("  ", basename(ex), ".refseq.fa         ASV sequences, named by asv_id"),
  if(has_tree) paste0("  ", basename(ex), ".tree              phylogenetic tree (Newick), tips named by asv_id"),
  if(length(extra)) c("  provenance: ", paste0("    ", extra)),
  "",
  "IDs: asv_id (ASV1, ASV2, ...) are the taxa names of the phyloseq object, numbered by decreasing total abundance",
  "before filtering; asv_seq in the tax table maps them to the sequences used as IDs in the unfiltered ASVtable.txt,",
  "taxonomy.txt and ASV.fa. Missing values are empty cells.",
  "",
  "In R:  library(phyloseq); ps = readRDS(\"<key>.phyloseq.filtered.rds\")   (the same data)")
writeLines(readme, paste0(ex, ".README.txt"))
cat("filtered exports:", paste(c(exports, paste0(basename(ex), ".README.txt")), collapse = ", "),
    "(bundled into", paste0(basename(ex), ".tar.gz"), "by the package step)\n")
cat("done:", ncol(seqtab.nochim), "ASVs in", nrow(seqtab.nochim), "samples ->", outdir, "\n")
