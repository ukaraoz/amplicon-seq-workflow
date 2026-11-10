# Read tracking figure: distribution of each readtracking.txt metric per sample, by sequencing run
#
# version 1.0
#
# usage: Rscript plot-readtracking.R --config=<key>.config.json        (workflow.py `figures` step)
#   reads   <results>/<key>.readtracking.txt and the resolved config (processing parameters for the subtitle)
#   writes  <results>/figures/<key>.readtracking_by_run.png and .pdf
# One panel per metric, in pipeline order (cutadapt_R2_trimmed left out; R1-only runs have no R2 metrics anyway):
# box (median, quartiles) and points (samples) per run on a log scale, the median per run as a number, and the
# failed-sample threshold (minfinal) on the last panel. The subtitle lists the cutadapt, filterAndTrim and DADA2
# settings from the config, including per-run overrides.
suppressPackageStartupMessages({library(ggplot2); library(jsonlite)})
args = commandArgs(TRUE)
cfgfile = sub("^--config=", "", args[grepl("^--config=", args)])
if(length(cfgfile) != 1) stop("usage: Rscript plot-readtracking.R --config=<key>.config.json")
cfg = fromJSON(cfgfile, simplifyVector = TRUE); p = cfg$params; key = cfg$key
tr = read.delim(file.path(cfg$paths$results, paste0(key, ".readtracking.txt")), check.names = FALSE)
if(!"run" %in% names(tr)) tr$run = "all"
figdir = file.path(cfg$paths$results, "figures"); dir.create(figdir, showWarnings = FALSE)

# ---- metrics, in pipeline order ----
desc = c(raw = "pairs into cutadapt", cutadapt_R1_trimmed = "R1 with primer / read-through removed",
         input = if("raw" %in% names(tr)) "pairs after cutadapt" else "input pairs", filtered = "after filterAndTrim",
         denoisedF = "R1 after DADA2", denoisedR = "R2 after DADA2", merged = "merged pairs",
         nonchim = "after chimera removal")
metrics = intersect(names(desc), setdiff(names(tr), c("sample", "run", "cutadapt_R2_trimmed")))
mlab = setNames(paste0(metrics, "\n(", desc[metrics], ")"), metrics)

# ---- runs: short labels (year + flowcell when the name has them), sample counts, colours ----
runs = unique(tr$run)
short = vapply(runs, function(r) { y = regmatches(r, regexpr("(19|20)[0-9]{2}", r)); f = sub(".*_", "", r)
  if(length(y) && f != r) paste(y, f) else paste(strwrap(gsub("_", " ", r), 14), collapse = "\n") }, "")
n = table(factor(tr$run, levels = runs))
rlab = setNames(paste0(short, "\n(n = ", as.integer(n), ")"), runs)
pal = c("#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
cols = setNames(pal[(seq_along(runs) - 1) %% length(pal) + 1], rlab)

long = do.call(rbind, lapply(metrics, function(m) data.frame(sample = tr$sample, run = rlab[tr$run], metric = mlab[[m]], reads = tr[[m]])))
long$metric = factor(long$metric, levels = mlab); long$run = factor(long$run, levels = rlab)
long = long[!is.na(long$reads) & long$reads > 0, ]
med = aggregate(reads ~ run + metric, long, median)
med$lab = formatC(med$reads, format = "d", big.mark = ",")

# ---- subtitle: processing parameters from the config ----
v = function(x, i = 1) if(is.null(x)) NA else x[min(i, length(x))]
r1only = !isTRUE(p$merge)
fmt_tl = function(tl) if(r1only) (if(v(tl) == 0) "0 (none)" else v(tl)) else paste(tl, collapse = "/")
lines = character(0)
if(isTRUE(p$run_cutadapt)) {
  # version: from a cutadapt report in <work>/cutadapt ("This is cutadapt X.Y ..."), else from cutadapt on PATH
  logs = Sys.glob(file.path(cfg$paths$work, "cutadapt", "*.cutadapt.log"))
  cv = if(length(logs)) sub("^This is cutadapt ([^ ]+).*", "\\1", readLines(logs[1], n = 1)) else ""
  if(!nzchar(cv) && nzchar(Sys.which("cutadapt"))) cv = tryCatch(system2(Sys.which("cutadapt"), "--version", stdout = TRUE, stderr = FALSE)[1], error = function(e) "")
  lines = c(lines, paste0("cutadapt", if(nzchar(cv)) paste0(" ", cv) else "", " on raw pairs: ", p$fwd_primer, " / ", p$rev_primer,
    " anchored at the read starts, read-through to the opposite primer removed at the 3' ends",
    if(!is.null(p$max_read_length)) paste0(", reads cut to ", p$max_read_length, " bp") else "",
    ", pairs < ", p$minLen, " bp dropped"))
}
lines = c(lines, paste0("filterAndTrim, ", if(r1only) "R1 only" else "R1 + R2", ": truncLen ", fmt_tl(p$truncLen),
  ", trimLeft ", if(is.null(p$trimLeft)) 0 else paste(if(r1only) v(p$trimLeft) else p$trimLeft, collapse = "/"), ", maxEE ",
  paste(if(r1only) v(p$maxEE) else p$maxEE, collapse = "/"), ", maxN 0, truncQ 2, minLen ", p$minLen,
  "; samples < ", format(p$mincount, big.mark = ","), " filtered reads dropped before denoising"))
if(length(p$runs)) lines = c(lines, paste0("per-run settings: ", paste(vapply(names(p$runs), function(r) {
  o = p$runs[[r]]; paste0(gsub("\n", " ", short[r]), " ", paste(vapply(names(o), function(k) paste(k, paste(o[[k]], collapse = "/")), ""), collapse = ", ")) }, ""),
  collapse = "; ")))
lines = c(lines, paste0("DADA2: error model learned per run, pseudo-pooled denoising per run",
  if(!r1only) paste0(", mergePairs (minOverlap ", if(is.null(p$minOverlap)) 12 else p$minOverlap, ")") else "",
  "; removeBimeraDenovo (consensus) on the combined table; ",
  if(length(p$asv_lengths) == 2) paste0("ASV length window ", p$asv_lengths[1], "-", p$asv_lengths[2], " bp") else "no ASV length window",
  "; samples < ", format(p$minfinal, big.mark = ","), " final reads reported as failed (dashed line)"))

last = mlab[[tail(metrics, 1)]]
g = ggplot(long, aes(run, reads, colour = run)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.18, height = 0, size = 0.9, alpha = 0.45) +
  geom_hline(data = data.frame(metric = factor(last, levels = mlab), y = p$minfinal), aes(yintercept = y),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_text(data = med, aes(y = reads, label = lab), colour = "black", size = 2.6, nudge_x = 0.42, hjust = 0) +
  facet_wrap(~ metric, nrow = 1) +
  scale_y_log10(labels = function(x) formatC(x, format = "d", big.mark = ","),
                breaks = as.vector(outer(c(1, 3), 10^(1:7)))) +
  scale_x_discrete(expand = expansion(add = c(0.6, 1.1))) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = NULL, y = "reads per sample (log scale)",
       title = paste0(cfg$project, " ", cfg$locus, " read tracking per sample, by sequencing run (", key, ", ", nrow(tr), " samples)"),
       subtitle = paste(lines, collapse = "\n"), caption = "box: median and quartiles; points: samples; numbers: median per run") +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_text(size = 7.6, lineheight = 1.1), strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 7.5), panel.grid.minor = element_blank())
w = max(8, length(metrics) * (0.85 * length(runs) + 0.9)); h = 5.4 + 0.18 * length(lines)
out = file.path(figdir, paste0(key, ".readtracking_by_run"))
ggsave(paste0(out, ".png"), g, width = w, height = h, dpi = 140)
ggsave(paste0(out, ".pdf"), g, width = w, height = h)
cat("plot-readtracking.R 1.0:", paste0(out, ".png / .pdf"), "\n")
