# DADA2 workflow: paired-end amplicon reads -> ASV table + taxonomy
#
# version 2.6
#
# usage: Rscript dada2-asv.R --config=<key>.config.json [--threads=N] [--nsamples=N] [--no-taxonomy]
#                            [--resume[=outputs]]
#   normally run by workflow.py (dada2 and phyloseq steps), which writes the resolved config
#   --config=FILE      resolved per-locus config (workflow.py): paths, sample naming, locus parameters, references
#   --threads=N        threads (default: config params.threads, else $SLURM_CPUS_PER_TASK, else 4)
#   --nsamples=N       test run on the first N samples (of each run); outputs go to <results>/test-nN and <work>/test-nN
#   --no-taxonomy      stop after the seqtab checkpoint (before assignTaxonomy)
#   --resume           load the seqtab checkpoint (<key>.dada2.seqtab.RData) and run only taxonomy + outputs
#   --resume=outputs   load the final <key>.dada2.RData (with taxonomy) and only rewrite the outputs
#                      (tables, fasta, phyloseq, failed samples), e.g. after building the tree
#   R1 only or merged pairs: config params.merge (workflow.py --merge / --no-merge); R1-only keys end in ".R1"
#
# failed samples (<key>.failed_samples.txt, with the reason): < mincount reads after quality filtering (dropped
#   before denoising), or < minfinal reads in the final table (after removing the config's remove_taxa /
#   require_ranks ASVs). Failed samples are kept in <key>.phyloseq.rds but removed from <key>.phyloseq.filtered.rds
#
# inputs  : paths.metadata (tab-separated; samples.id_column = sample id); fastq from samples.manifest (R1, R2
#           columns) or paths.raw/<samples.r1|r2>
# outputs : paths.results (QC plots, read tracking, ASV table/fasta, taxonomy, RData, phyloseq rds)
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
resume_mode = if("--resume" %in% opts) "taxonomy" else optval("resume", NA)
stopifnot(resume_mode %in% c(NA, "taxonomy", "outputs"))
resume = !is.na(resume_mode)
mincount = cfg$params$mincount   # samples with fewer filtered reads are dropped before denoising
minfinal = cfg$params$minfinal   # samples with fewer reads in the final table are reported as failed
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
script_version = "2.6"   # bump on every modification (2.0, 2.1, ...)
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
metadata = read.table(cfg$paths$metadata, sep = "\t", header = T, check.names = F)
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
  ggsave(file.path(outdir, paste0(key, ".quality.raw.R1.pdf")), plotQualityProfile(Ffiles[qc.samples]), width = 12, height = 8)
  ggsave(file.path(outdir, paste0(key, ".quality.raw.R2.pdf")), plotQualityProfile(Rfiles[qc.samples]), width = 12, height = 8)

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
    ggsave(file.path(outdir, paste0(key, tag, ".errors.R1.pdf")), plotErrors(errF[[r]], nominalQ = TRUE), width = 10, height = 8)
    if(params$merge) ggsave(file.path(outdir, paste0(key, tag, ".errors.R2.pdf")), plotErrors(errR[[r]], nominalQ = TRUE),
                            width = 10, height = 8)
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
  track.towrite = tibble::rownames_to_column(as.data.frame(track), "sample")
  if(length(runnames) > 1) track.towrite = cbind(track.towrite[1], run = sampleruns[track.towrite$sample], track.towrite[-1])
  write.table(track.towrite, file = file.path(outdir, paste0(key, ".readtracking.txt")),
              sep = "\t", row.names = F, col.names = T, quote = F)
  # checkpoint before taxonomy: assignTaxonomy against SILVA is memory hungry and can get the job killed
  save(out, keep, sampleruns, errF, errR, dadaFs, dadaRs, mergers, seqtab, seqtab.nochim, track, params,
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

save(out, keep, errF, errR, dadaFs, dadaRs, mergers, seqtab, seqtab.nochim, taxa, track, params,
     file = file.path(outdir, paste0(key, ".dada2.RData")))

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
saveRDS(ps, file.path(outdir, paste0(key, ".phyloseq.rds")))
print(ps)
# filtered copy: drops the config's remove_taxa (rank: names, e.g. Order: Chloroplast, Family: Mitochondria) and
# ASVs unassigned at any of require_ranks (e.g. Phylum)
psf = ps
tt = as.data.frame(as(tax_table(ps), "matrix"))
drop = rep(FALSE, ntaxa(ps))
for(rank in names(p$remove_taxa)) drop = drop | tt[[rank]] %in% p$remove_taxa[[rank]]
for(rank in p$require_ranks) drop = drop | is.na(tt[[rank]])
if(any(drop)) {
  what = c(unlist(lapply(names(p$remove_taxa), function(r) paste(p$remove_taxa[[r]], collapse = "/"))),
           if(length(p$require_ranks)) paste0(paste(p$require_ranks, collapse = "/"), "-NA"))
  cat("removing", sum(drop), paste(what, collapse = "/"), "ASVs (", round(100 * sum(taxa_sums(ps)[drop]) / sum(taxa_sums(ps)), 2), "% of reads )\n")
  psf = prune_taxa(!drop, ps)
}
sample_data(psf)$reads_filtered_taxa = sample_sums(psf)

# ---- failed samples ----
# context for each reason: the sample's reads vs the median (and the lowest passing sample) at the same stage
final = sample_sums(psf)
lowcount = setdiff(rownames(out), keep)
lowfinal = names(final)[final < minfinal]
pass = setdiff(names(final), lowfinal)
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
             reason = sprintf("low number of reads: %d reads in the final table (< %d), %.1f%% of the median (%d); lowest passing sample %s has %d; input %d reads vs median %d",
                              as.integer(final[lowfinal]), minfinal, 100 * final[lowfinal] / median(final), as.integer(median(final)),
                              nextlow, as.integer(final[nextlow]), out[lowfinal, "reads.in"], as.integer(med_in))))
write.table(failed, file = file.path(outdir, paste0(key, ".failed_samples.txt")), sep = "\t", row.names = F, col.names = T, quote = F)
cat(nrow(failed), "failed samples ->", file.path(outdir, paste0(key, ".failed_samples.txt")), "\n")
if(nrow(failed) > 0) print(failed[, c("sample", "stage", "reads_input", "reads_at_failure")], row.names = F)
psf = prune_samples(pass, psf)
saveRDS(psf, file.path(outdir, paste0(key, ".phyloseq.filtered.rds")))
print(psf)
cat("done:", ncol(seqtab.nochim), "ASVs in", nrow(seqtab.nochim), "samples ->", outdir, "\n")
