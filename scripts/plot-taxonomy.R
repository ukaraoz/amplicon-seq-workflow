# Taxonomic composition per sample, one figure per rank
#
# version 1.0
#
# usage: Rscript plot-taxonomy.R --config=<key>.config.json [--top=15] [--pdf]     (workflow.py `figures` step)
#   reads   <results>/<key>.phyloseq.filtered.rds (after the taxa filter, failed samples removed)
#   writes  <results>/<key>.taxonomy_<Rank>.png for every rank below the first (e.g. Phylum ... Species for UNITE);
#           with --pdf also .pdf
# Stacked bars: relative abundance per sample (share of its reads), the --top taxa by mean relative abundance across
# samples, the rest as "Other", reads not assigned at that rank (or only to a placeholder such as
# Ascomycota_gen_Incertae_sedis) as "unassigned". One panel per sampling year (metadata column `year`), else per
# sequencing run (manifest run column); samples ordered within a panel by the abundance of the most abundant taxa. Names without UNITE's rank prefixes; species as
# "Genus species".
suppressPackageStartupMessages({library(phyloseq); library(ggplot2); library(jsonlite)})
args = commandArgs(TRUE)
cfgfile = sub("^--config=", "", args[grepl("^--config=", args)])
if(length(cfgfile) != 1) stop("usage: Rscript plot-taxonomy.R --config=<key>.config.json [--top=15] [--pdf]")
top = as.integer(sub("^--top=", "", c(args[grepl("^--top=", args)], "--top=15")[1]))
want_pdf = "--pdf" %in% args
cfg = fromJSON(cfgfile, simplifyVector = TRUE); key = cfg$key; outdir = cfg$paths$results
ps = readRDS(file.path(outdir, paste0(key, ".phyloseq.filtered.rds")))

otu = as(otu_table(ps), "matrix"); if(taxa_are_rows(ps)) otu = t(otu)          # samples x taxa
tt = as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
tt[] = lapply(tt, function(x) { x = sub("^[a-z]__", "", x); x[grepl("Incertae_sedis|_gen_|_fam_|_ord_|_cls_|_phy_|^unidentified", x)] = NA; x })
ranks = names(tt)[-1]                                                           # below Kingdom
if("Species" %in% ranks && "Genus" %in% ranks)
  tt$Species = ifelse(is.na(tt$Species) | is.na(tt$Genus), NA, paste(tt$Genus, gsub("_", " ", tt$Species)))
rel = otu / rowSums(otu)

# panels: sampling year (metadata `year`), else sequencing run (manifest run column), else one group
s = sample_data(ps); sdn = s@names
run = rep("all samples", nsamples(ps))
if("year" %in% sdn) {
  run = as.character(s@.Data[[match("year", sdn)]])
} else if(!is.null(cfg$samples$manifest) && !is.null(cfg$samples$run_column) && file.exists(cfg$samples$manifest)) {
  man = read.delim(cfg$samples$manifest, check.names = FALSE, colClasses = "character")
  idc = if(cfg$samples$id_column %in% names(man)) cfg$samples$id_column else "sample"
  r = man[[cfg$samples$run_column]][match(sample_names(ps), man[[idc]])]
  if(!anyNA(r)) run = r
}
short = vapply(unique(run), function(r) { y = regmatches(r, regexpr("(19|20)[0-9]{2}", r)); f = sub(".*_", "", r)
  if(length(y) && f != r) paste(y, f) else r }, "")
nrun = table(run)
runs = sort(unique(run))
runlab = setNames(paste0(short[runs], " (n = ", as.integer(nrun[runs]), ")"), runs)

pal = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf",
        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5",
        "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#e6550d")
desc = paste0(cfg$project, " ", cfg$locus, " (", key, "): ", nsamples(ps), " samples, ", ntaxa(ps), " ASVs, ",
              format(sum(otu), big.mark = ","), " reads (phyloseq.filtered.rds",
              if(length(cfg$params$keep_taxa)) paste0(", ", paste(names(cfg$params$keep_taxa), unlist(cfg$params$keep_taxa), collapse = ", "), " only") else "", ")")

for(rk in ranks) {
  lab = ifelse(is.na(tt[[rk]]), "unassigned", tt[[rk]])
  m = t(rowsum(t(rel), lab))                                                    # samples x taxa at this rank
  mean_ab = sort(colMeans(m[, colnames(m) != "unassigned", drop = FALSE]), decreasing = TRUE)
  keep = names(mean_ab)[seq_len(min(top, length(mean_ab)))]
  other = setdiff(colnames(m), c(keep, "unassigned"))
  mm = cbind(m[, keep, drop = FALSE],
             Other = if(length(other)) rowSums(m[, other, drop = FALSE]) else 0,
             unassigned = if("unassigned" %in% colnames(m)) m[, "unassigned"] else 0)
  # sample order within run: by the most abundant kept taxon, then the next
  o = do.call(order, c(list(run), lapply(keep[seq_len(min(3, length(keep)))], function(k) -mm[, k])))
  smp = rownames(mm)[o]
  long = data.frame(sample = rep(rownames(mm), ncol(mm)), taxon = rep(colnames(mm), each = nrow(mm)), ab = as.vector(mm))
  long$run = factor(runlab[run[match(long$sample, rownames(mm))]], levels = runlab)
  long$sample = factor(long$sample, levels = smp)
  lv = c(keep, "Other", "unassigned")
  long$taxon = factor(long$taxon, levels = rev(lv))
  cols = c(setNames(pal[seq_along(keep)], keep), Other = "#b3b3b3", unassigned = "#e8e8e8")
  pct_un = 100 * mean(mm[, "unassigned"]); pct_oth = 100 * mean(mm[, "Other"])
  g = ggplot(long, aes(sample, ab, fill = taxon)) +
    geom_col(width = 1) +
    facet_grid(. ~ run, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cols, breaks = lv, name = rk) +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%"), expand = c(0, 0)) +
    labs(x = paste0("samples, by ", if("year" %in% sdn) "sampling year" else "sequencing run", " (ordered within each panel by the most abundant taxa)"), y = "relative abundance (share of sample reads)",
         title = paste0(cfg$locus, " taxonomic composition per sample: ", rk),
         subtitle = paste0(desc, "\n", length(keep), " most abundant ", tolower(rk), " taxa by mean relative abundance (",
                           length(mean_ab), " in total); Other = the remaining ", length(other), " (mean ", sprintf("%.1f", pct_oth),
                           "%); unassigned = not assigned at ", tolower(rk), " level (mean ", sprintf("%.1f", pct_un), "%)")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
          panel.spacing = unit(3, "pt"), strip.text = element_text(size = 9), plot.subtitle = element_text(size = 8, lineheight = 1.1),
          legend.key.size = unit(10, "pt"), legend.text = element_text(size = 8,
            face = if(rk %in% c("Genus", "Species")) "italic" else "plain"))
  out = file.path(outdir, paste0(key, ".taxonomy_", rk))
  ggsave(paste0(out, ".png"), g, width = 14, height = 6.5, dpi = 140)
  if(want_pdf) ggsave(paste0(out, ".pdf"), g, width = 14, height = 6.5)
  cat(sprintf("%-8s %5d taxa; top %d cover %.1f%% on average, Other %.1f%%, unassigned %.1f%% -> %s.png%s\n", rk, length(mean_ab),
              length(keep), 100 * mean(rowSums(mm[, keep, drop = FALSE])), pct_oth, pct_un, basename(out), if(want_pdf) " + .pdf" else ""))
}
cat("plot-taxonomy.R 1.0 done\n")
