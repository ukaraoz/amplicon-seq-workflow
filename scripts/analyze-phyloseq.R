# Community analysis of the DADA2 phyloseq objects (output of dada2-asv.R), for survey designs such as NEON ER18+ER22
#
# version 1.1
#
# usage: Rscript analyze-phyloseq.R --config=<key>.config.json [--pdf]
#   normally run by workflow.py (analyze step, analysis.script), which writes the resolved config
#   --pdf  also write PDFs of the figures (default: PNG only)
# input : paths.results/<key>.phyloseq.filtered.rds (taxa filter applied, failed samples removed; tree from the
#         tree step, if any)
# output: paths.analysis/<key>.<analysis>....png (and .pdf with --pdf)
#
# design (config analysis):
#   sample_order: {by: [columns], levels: {column: [values]}}: the order of the samples in every figure and table:
#     sorted by the `by` sample-data columns in turn (numbers numerically, text alphabetically or in the given
#     `levels`; missing values last), then by sample name; written to <key>.sample_order.txt
#   groups: the sample-data column that splits the figures into panels (in sample order); without it (or when the
#     column is missing) the panels are the sequencing runs (samples.run_column), else one group
# Every analysis has its own config section and runs only when the section is present.
#
# analyses:
#   1. taxonomic composition (taxonomic_composition): for every rank below the first, relative abundance per sample
#      (share of its reads) as stacked bars: the ntop taxa by mean relative abundance across samples, the rest as
#      "Other", reads not assigned at that rank (or only to a placeholder such as Ascomycota_gen_Incertae_sedis) as
#      "unassigned"; one panel per group, samples in sample order (order: abundance: by the most abundant taxa
#      within each panel instead); band: a sample-data column shown as a coloured strip under the bars (e.g. year);
#      names without UNITE's rank prefixes, species as "Genus species"
#
# All analyses of the phyloseq objects go in this script, as further numbered sections.

packages = c("phyloseq", "ggplot2")
invisible(lapply(packages, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))
script_version = "1.1"

allargs = commandArgs(trailingOnly = TRUE)
write_pdf = "--pdf" %in% allargs
config_file = sub("^--config=", "", allargs[grepl("^--config=", allargs)])
if(length(config_file) != 1) stop("usage: Rscript analyze-phyloseq.R --config=<key>.config.json [--pdf] (see workflow.py)")
cfg = jsonlite::read_json(config_file, simplifyVector = TRUE)
acfg = if(is.null(cfg$analysis)) list() else cfg$analysis
key = cfg$key
resdir = cfg$paths$results
outdir = if(is.null(cfg$paths$analysis)) resdir else cfg$paths$analysis
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat("analyze-phyloseq.R version", script_version, ":", key, "->", outdir, "\n")

# figures: always PNG; PDF only with --pdf
save_plot = function(p, name, width, height, dpi = 140) {
  f = file.path(outdir, paste0(key, ".", name))
  ggsave(paste0(f, ".png"), p, width = width, height = height, dpi = dpi)
  if(write_pdf) ggsave(paste0(f, ".pdf"), p, width = width, height = height)
  basename(f)
}

# ---- data ----
ps = readRDS(file.path(resdir, paste0(key, ".phyloseq.filtered.rds")))
otu = as(otu_table(ps), "matrix"); if(taxa_are_rows(ps)) otu = t(otu)       # samples x taxa
# sample data with its original column names (phyloseq's data.frame conversions would turn "Soil pH" into "Soil.pH")
sample_df = function(x) { s = sample_data(x); d = as.data.frame(s@.Data, col.names = seq_along(s@names))
  names(d) = s@names; rownames(d) = sample_names(x); d }
sd = sample_df(ps)
# taxonomy without UNITE's rank prefixes; placeholder names count as unassigned; species as "Genus species"
tax = as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
tax[] = lapply(tax, function(x) { x = sub("^[a-z]__", "", x)
  x[grepl("Incertae_sedis|_gen_|_fam_|_ord_|_cls_|_phy_|^unidentified", x)] = NA; x })
if(all(c("Genus", "Species") %in% names(tax)))
  tax$Species = ifelse(is.na(tax$Species) | is.na(tax$Genus), NA, paste(tax$Genus, gsub("_", " ", tax$Species)))
desc = paste0(cfg$project, " ", cfg$locus, " (", key, "): ", nsamples(ps), " samples, ", ntaxa(ps), " ASVs, ",
              format(sum(otu), big.mark = ","), " reads (phyloseq.filtered.rds",
              if(length(cfg$params$keep_taxa)) paste0("; ", paste(names(cfg$params$keep_taxa), unlist(cfg$params$keep_taxa), collapse = ", "), " only") else "", ")")
cat(desc, "\n")

# ---- sample order (analysis.sample_order) ----
# all figures and tables below follow this order; ps and otu are reordered, sd keeps the original column names
so = acfg$sample_order
by = if(is.null(so$by)) character(0) else unlist(so$by)
miss = setdiff(by, names(sd))
if(length(miss)) stop("analysis.sample_order.by: no sample-data column ", paste(miss, collapse = ", "))
sort_key = function(cn) {
  x = sd[[cn]]; lv = so$levels[[cn]]
  if(!is.null(lv)) { k = match(as.character(x), as.character(lv)); if(anyNA(k[!is.na(x)])) stop("sample_order.levels.", cn, " misses ",
                       paste(setdiff(unique(as.character(x[!is.na(x)])), lv), collapse = ", ")); k }
  else if(is.numeric(x)) x else as.character(x)
}
sample_order = rownames(sd)[do.call(order, c(lapply(by, sort_key), list(rownames(sd)), list(na.last = TRUE)))]
sd = sd[sample_order, , drop = FALSE]
otu = otu[sample_order, , drop = FALSE]
ps = prune_samples(sample_order, ps)
otu_table(ps) = otu_table(if(taxa_are_rows(ps)) t(otu) else otu, taxa_are_rows = taxa_are_rows(ps))
stopifnot(identical(sample_names(ps), sample_order), identical(rownames(otu), sample_order))
write.table(data.frame(order = seq_along(sample_order), sample = sample_order, sd[, by, drop = FALSE], check.names = FALSE),
            file.path(outdir, paste0(key, ".sample_order.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
cat("sample order:", if(length(by)) paste(by, collapse = " > ") else "sample name", "->", paste0(key, ".sample_order.txt"), "\n")

# ---- groups (figure panels) ----
# analysis.groups: a sample-data column (e.g. year); else the sequencing run (manifest run column); else one group
run_short = function(r) { y = regmatches(r, regexpr("(19|20)[0-9]{2}", r)); f = sub(".*_", "", r)
  if(length(y) && f != r) paste(y, f) else r }
group = rep("all samples", nsamples(ps)); group_name = "all samples"
if(!is.null(acfg$groups) && acfg$groups %in% names(sd)) {
  group = as.character(sd[[acfg$groups]]); group_name = acfg$groups
} else if(!is.null(cfg$samples$manifest) && !is.null(cfg$samples$run_column) && file.exists(cfg$samples$manifest)) {
  man = read.delim(cfg$samples$manifest, check.names = FALSE, colClasses = "character")
  idc = if(cfg$samples$id_column %in% names(man)) cfg$samples$id_column else "sample"
  r = man[[cfg$samples$run_column]][match(rownames(sd), man[[idc]])]
  if(!anyNA(r)) { group = vapply(r, run_short, ""); group_name = "sequencing run" }
}
group[is.na(group)] = "NA"
glev = unique(group)                       # panels in sample order (first appearance)
glab = setNames(paste0(glev, " (n = ", as.integer(table(group)[glev]), ")"), glev)
cat("groups (", group_name, "):", paste(glab, collapse = ", "), "\n")

pal = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf",
        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5",
        "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#e6550d")

# ---- 1. taxonomic composition (analysis.taxonomic_composition) ----
if(!is.null(acfg$taxonomic_composition)) {
  tc = acfg$taxonomic_composition
  ntop = if(is.null(tc$ntop)) 15 else as.integer(tc$ntop)
  ranks = if(!is.null(tc$ranks)) intersect(tc$ranks, names(tax)) else names(tax)[-1]     # below the first rank
  by_abundance = identical(tc$order, "abundance")
  band = if(!is.null(tc$band) && tc$band %in% names(sd)) tc$band else NULL
  band_cols = c("#4d4d4d", "#bdbdbd", "#d95f02", "#1b9e77", "#7570b3", "#e7298a")
  rel = otu / rowSums(otu)
  for(rk in ranks) {
    lab = ifelse(is.na(tax[[rk]]), "unassigned", tax[[rk]])
    m = t(rowsum(t(rel), lab))                                                # samples x taxa at this rank
    mean_ab = sort(colMeans(m[, colnames(m) != "unassigned", drop = FALSE]), decreasing = TRUE)
    keep = names(mean_ab)[seq_len(min(ntop, length(mean_ab)))]
    other = setdiff(colnames(m), c(keep, "unassigned"))
    mm = cbind(m[, keep, drop = FALSE],
               Other = if(length(other)) rowSums(m[, other, drop = FALSE]) else 0,
               unassigned = if("unassigned" %in% colnames(m)) m[, "unassigned"] else 0)
    o = if(by_abundance) do.call(order, c(list(match(group, glev)), lapply(keep[seq_len(min(3, length(keep)))], function(k) -mm[, k])))
        else seq_len(nrow(mm))                                                  # sample order
    long = data.frame(sample = factor(rep(rownames(mm), ncol(mm)), levels = rownames(mm)[o]),
                      taxon = factor(rep(colnames(mm), each = nrow(mm)), levels = rev(c(keep, "Other", "unassigned"))),
                      group = factor(rep(glab[group], ncol(mm)), levels = glab), ab = as.vector(mm))
    lv = c(keep, "Other", "unassigned")
    bd = if(!is.null(band)) data.frame(sample = factor(rownames(mm), levels = rownames(mm)[o]),
                                       group = factor(glab[group], levels = glab),
                                       band = factor(as.character(sd[rownames(mm), band]))) else NULL
    cols = c(setNames(pal[seq_along(keep)], keep), Other = "#b3b3b3", unassigned = "#e8e8e8")
    pct_un = 100 * mean(mm[, "unassigned"]); pct_oth = 100 * mean(mm[, "Other"])
    g = ggplot(long, aes(sample, ab, fill = taxon)) +
      geom_col(width = 1) +
      (if(!is.null(bd)) geom_linerange(data = bd, aes(x = sample, ymin = -0.045, ymax = -0.01, colour = band), inherit.aes = FALSE,
                                       linewidth = 1.2) else NULL) +
      (if(!is.null(bd)) scale_colour_manual(values = setNames(band_cols[seq_along(levels(bd$band))], levels(bd$band)), name = band) else NULL) +
      facet_grid(. ~ group, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = cols, breaks = lv, name = rk) +
      scale_y_continuous(labels = function(x) paste0(100 * x, "%"), breaks = seq(0, 1, 0.25), expand = c(0, 0),
                         limits = c(if(!is.null(bd)) -0.05 else 0, 1)) +
      labs(x = paste0("samples, by ", group_name, if(by_abundance) " (ordered within each panel by the most abundant taxa)"
                      else paste0(" (in sample order: ", paste(by, collapse = " > "), ")")),
           y = "relative abundance (share of sample reads)",
           title = paste0(cfg$locus, " taxonomic composition per sample: ", rk),
           subtitle = paste0(desc, "\n", length(keep), " most abundant ", tolower(rk), " taxa by mean relative abundance (",
                             length(mean_ab), " in total); Other = the remaining ", length(other), " (mean ", sprintf("%.1f", pct_oth),
                             "%); unassigned = not assigned at ", tolower(rk), " level (mean ", sprintf("%.1f", pct_un), "%)")) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
            panel.spacing = unit(3, "pt"), strip.text = element_text(size = 8, angle = if(length(glev) > 6) 90 else 0), plot.subtitle = element_text(size = 8, lineheight = 1.1),
            legend.key.size = unit(10, "pt"),
            legend.text = element_text(size = 8, face = if(rk %in% c("Genus", "Species")) "italic" else "plain"))
    f = save_plot(g, paste0("taxonomy_", rk), width = 14, height = 6.5)
    cat(sprintf("1. taxonomic composition %-8s %5d taxa; top %d cover %.1f%% on average, Other %.1f%%, unassigned %.1f%% -> %s.png%s\n",
                rk, length(mean_ab), length(keep), 100 * mean(rowSums(mm[, keep, drop = FALSE])), pct_oth, pct_un, f,
                if(write_pdf) " + .pdf" else ""))
  }
}

cat("analyze-phyloseq.R done\n")
