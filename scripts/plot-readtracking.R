# Read tracking figure: distribution of each readtracking.txt metric per sample, by sequencing run
#
# version 1.5
#
# usage: Rscript plot-readtracking.R --config=<key>.config.json [--pdf]        (workflow.py `figures` step)
#   reads   <results>/<key>.readtracking.txt and the resolved config
#   writes  <results>/<key>.readtracking_by_run.png          reads per sample
#           <results>/<key>.readtracking_pct_by_run.png      % of the sample's input reads
#           with --pdf also the same as .pdf
#           (in the results folder itself, next to the other plots; no figures/ subfolder)
# One panel per metric, in pipeline order (cutadapt_R2_trimmed left out; R1-only runs have no R2 metrics anyway):
# box (median, quartiles) and points (samples) per run on a log scale, the median per run as a number, and the
# failed-sample threshold (minfinal) on the last panel. No subtitle.
suppressPackageStartupMessages({library(ggplot2); library(jsonlite)})
args = commandArgs(TRUE)
cfgfile = sub("^--config=", "", args[grepl("^--config=", args)])
if(length(cfgfile) != 1) stop("usage: Rscript plot-readtracking.R --config=<key>.config.json [--pdf]")
want_pdf = "--pdf" %in% args
cfg = fromJSON(cfgfile, simplifyVector = TRUE); p = cfg$params; key = cfg$key
tr = read.delim(file.path(cfg$paths$results, paste0(key, ".readtracking.txt")), check.names = FALSE)
if(!"run" %in% names(tr)) tr$run = "all"
figdir = cfg$paths$results

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
       caption = "box: median and quartiles; points: samples; numbers: median per run") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 7.5), panel.grid.minor = element_blank())
w = max(8, length(metrics) * (0.85 * length(runs) + 0.9)); h = 5.4
out = file.path(figdir, paste0(key, ".readtracking_by_run"))
ggsave(paste0(out, ".png"), g, width = w, height = h, dpi = 300)
if(want_pdf) ggsave(paste0(out, ".pdf"), g, width = w, height = h)
cat("plot-readtracking.R 1.5:", paste0(out, if(want_pdf) ".png / .pdf" else ".png"), "\n")

# ---- the same metrics as % of each sample's input reads (pairs after cutadapt, or input pairs) ----
pm = setdiff(metrics, c("raw", "input"))   # raw / input is ~100% (cutadapt drops almost no pairs; counts in the other figure)
plong = do.call(rbind, lapply(pm, function(m) data.frame(sample = tr$sample, run = rlab[tr$run], metric = mlab[[m]],
                                                         pct = 100 * tr[[m]] / tr$input)))
plong$metric = factor(plong$metric, levels = mlab[pm]); plong$run = factor(plong$run, levels = rlab)
plong = plong[is.finite(plong$pct), ]
pmed = aggregate(pct ~ run + metric, plong, median)
pmed$lab = sprintf("%.1f%%", pmed$pct)
ref = if("raw" %in% names(tr)) "input = pairs after cutadapt" else "input = input pairs"
g2 = ggplot(plong, aes(run, pct, colour = run)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.18, height = 0, size = 0.9, alpha = 0.45) +
  geom_text(data = pmed, aes(y = pct, label = lab), colour = "black", size = 2.6, nudge_x = 0.42, hjust = 0) +
  facet_wrap(~ metric, nrow = 1) +
  scale_y_continuous(limits = c(0, max(100, plong$pct) * 1.02), breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%")) +
  scale_x_discrete(expand = expansion(add = c(0.6, 1.1))) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = NULL, y = paste0("% of the sample's input reads (", ref, ")"),
       title = paste0(cfg$project, " ", cfg$locus, " read tracking per sample as % of input reads, by sequencing run (", key, ", ", nrow(tr), " samples)"),
       caption = paste0("each metric divided by the same sample's input reads (", ref, "); box: median and quartiles; points: samples; numbers: median per run")) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 7.5), panel.grid.minor = element_blank())
w2 = max(8, length(pm) * (0.85 * length(runs) + 0.9))
out2 = file.path(figdir, paste0(key, ".readtracking_pct_by_run"))
ggsave(paste0(out2, ".png"), g2, width = w2, height = h, dpi = 300)
if(want_pdf) ggsave(paste0(out2, ".pdf"), g2, width = w2, height = h)
cat("plot-readtracking.R 1.5:", paste0(out2, if(want_pdf) ".png / .pdf" else ".png"), "\n")
