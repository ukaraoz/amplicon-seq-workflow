# Community analysis of the DADA2 phyloseq objects (output of dada2-asv.R)
#
# version 2.0
#
# usage: Rscript analyze-phyloseq.R --config=<key>.config.json [--pdf]
#   normally run by workflow.py (analyze step), which writes the resolved config
#   --pdf  also write PDFs of the figures (default: PNG only)
# input : paths.results/<key>.phyloseq.filtered.rds (remove_taxa / require_ranks ASVs and failed samples removed;
#         tree from the tree step, if any)
# output: paths.analysis/
#
# experimental design (config analysis): metadata columns are mapped to the internal names habitat, treatment,
# pot (replicate), is_T0 (baseline flag) and is_pellet (source flag); habitat and treatment levels, labels and colors,
# the baseline and source sample labels, and the settings of each analysis come from the config. Optional analyses
# (venn, source_sharing, reference_distances, betapart) run only when their config section is present.
#
# analyses:
#   1. taxonomic composition: relative-abundance bar plots per rank, in sample order
#   2. alpha diversity (rarefied counts): observed, Shannon, inverse Simpson, Pielou, Faith's PD (with a tree)
#   3. ASV overlap among a set of treatments (venn), per habitat: n-way Venn diagrams and pairwise reference vs each other
#   4. NMDS (with a tree: weighted and unweighted UniFrac; without: Bray-Curtis and Jaccard presence/absence), per
#      ordination subset of habitats
#   5. variance partitioning (PERMANOVA, adonis2), formula per ordination subset
#   6. pairwise PERMANOVA between treatments within each habitat (BH-adjusted)
#   7. source establishment (source_sharing): ASVs shared between the source sample (e.g. an inoculum) and each
#      treatment sample, per habitat (rarefied per habitat, n_rare repeats)
#   8. distances between a reference treatment and every other treatment, per habitat (reference_distances)
#   9. source colonization vs pruning/boosting (betapart): betapart partition of with- vs without-source pot pairs
#      into turnover and nestedness (Sorensen, presence/absence), per habitat (rarefied per habitat, n_rare repeats)

packages = c("phyloseq", "vegan", "phangorn", "ggplot2", "ggrepel", "dplyr", "tibble")
invisible(lapply(packages, library, character.only = TRUE))
pdf(NULL)   # null graphics device: functions that open one (vegan::ordiellipse) do not leave an Rplots.pdf behind

allargs = commandArgs(trailingOnly = TRUE)
write_pdf = "--pdf" %in% allargs
config_file = sub("^--config=", "", allargs[grepl("^--config=", allargs)])
if(length(config_file) != 1) stop("usage: Rscript analyze-phyloseq.R --config=<key>.config.json [--pdf] (see workflow.py)")
cfg = jsonlite::read_json(config_file, simplifyVector = TRUE)
acfg = cfg$analysis
locus = cfg$locus
run = cfg$run
key = cfg$key
resdir = cfg$paths$results
outdir = cfg$paths$analysis
dir.create(outdir, showWarnings = F, recursive = T)
# figures: always PNG; PDF only with --pdf
save_plot = function(p, name, width, height, dpi = 130) {
  ggsave(file.path(outdir, paste0(key, ".", name, ".png")), p, width = width, height = height, dpi = dpi)
  if(write_pdf) ggsave(file.path(outdir, paste0(key, ".", name, ".pdf")), p, width = width, height = height)
}
set.seed(100)

# ASV counts -> relative abundances (each sample sums to 1)
to_relative = function(ps) {
  stopifnot(all(sample_sums(ps) > 0))
  rel = transform_sample_counts(ps, function(x) x / sum(x))
  stopifnot(isTRUE(all.equal(unname(sample_sums(rel)), rep(1, nsamples(rel)))))
  rel
}

ps = readRDS(file.path(resdir, paste0(key, ".phyloseq.filtered.rds")))
print(ps)
# all analyses below use relative abundances; the read counts stay in ps_counts (for analyses that need counts,
# e.g. richness or rarefaction)
ps_counts = ps
ps = to_relative(ps)

# ---- design (config analysis) ----
# metadata columns -> internal names; the baseline and source flags default to 0 when the config has no column
sd = data.frame(sample_data(ps), check.names = F)
cols = acfg$columns
for(nm in c("habitat", "treatment", "pot")) sd[[nm]] = sd[[cols[[c(habitat = "habitat", treatment = "treatment", pot = "replicate")[[nm]]]]]]
sd$is_T0 = if(!is.null(cols$baseline)) sd[[cols$baseline]] else 0
sd$is_pellet = if(!is.null(cols$source)) sd[[cols$source]] else 0
baseline_label = if(!is.null(acfg$baseline$label)) acfg$baseline$label else "T0"
baseline_habitat = if(!is.null(acfg$baseline$habitat)) acfg$baseline$habitat else acfg$habitats[1]
# the source sample (e.g. an inoculum) has no pot/treatment/habitat in the metadata; label it so it shows in plots
has_source = !is.null(acfg$source) && any(sd$is_pellet %in% 1)
source_label = if(has_source) acfg$source$label else character(0)      # e.g. "SMT pellet"
source_short = if(has_source) acfg$source$short else character(0)      # panel label, e.g. "SMT"
if(has_source) {
  sd$habitat[is.na(sd$habitat) & sd$is_pellet %in% 1] = source_label
  sd$treatment[is.na(sd$treatment) & sd$is_pellet %in% 1] = source_label
}
# baseline samples carry the baseline label as their habitat
sd$habitat[sd$is_T0 %in% 1] = baseline_label

# ---- sample order ----
# habitat (config order), then treatment, then within each treatment its baseline (T0) first and the pots by
# replicate; the source sample last. The baseline goes at the head of its treatment in its habitat block.
# habitat and treatment become factors with these levels, so plots and tables follow the same order.
habitat_levels = c(acfg$habitats, baseline_label, source_label)
treatment_levels = c(acfg$treatments$levels, source_label)
stopifnot(all(sd$habitat %in% habitat_levels), all(sd$treatment %in% treatment_levels))
# treatment names shown in the figures (data, tables and file names keep the codes). treatment_name(x, width):
# width = NA breaks long names as set in the config (treatments.wrap: pattern -> replacement with a newline; axes,
# legends); a number wraps at that many characters (narrow label blocks)
treatment_names = c(unlist(acfg$treatments$labels), setNames(source_label, source_label))
treatment_name = function(x, width = Inf) {
  x = as.character(x)
  y = ifelse(x %in% names(treatment_names), treatment_names[x], x)
  if(is.na(width)) { if(!is.null(acfg$treatments$wrap)) y = sub(acfg$treatments$wrap$pattern, acfg$treatments$wrap$replacement, y) }
  else if(is.finite(width)) y = vapply(y, function(t) paste(strwrap(t, width), collapse = "\n"), "")
  unname(y)
}
treatment_axis_labels = function(x) treatment_name(x, NA)
# x labels of "<block>|<treatment>" or "<block>|<treatment> <baseline label>" groups
group_axis_labels = function(x) {
  g = sub("^[^|]*[|]", "", x); suffix = paste0(" ", baseline_label)
  is_bl = endsWith(g, suffix); g[is_bl] = substr(g[is_bl], 1, nchar(g[is_bl]) - nchar(suffix))
  paste0(treatment_axis_labels(g), ifelse(is_bl, suffix, ""))
}
habitat_block = ifelse(sd$habitat == baseline_label, baseline_habitat, sd$habitat)
sample_order = rownames(sd)[order(sd$is_pellet %in% 1,
                                  match(habitat_block, habitat_levels),
                                  match(sd$treatment, treatment_levels),
                                  !(sd$is_T0 %in% 1),
                                  sd$pot)]
sd = sd[sample_order, ]
sd$habitat = factor(sd$habitat, levels = habitat_levels)
sd$treatment = factor(sd$treatment, levels = treatment_levels)
otu = as(otu_table(ps), "matrix")
otu_table(ps) = otu_table(if(taxa_are_rows(ps)) otu[, sample_order] else otu[sample_order, ], taxa_are_rows = taxa_are_rows(ps))
sample_data(ps) = sample_data(sd)
stopifnot(identical(sample_names(ps), sample_order))
ps_counts = prune_samples(sample_order, ps_counts)
otu_counts = as(otu_table(ps_counts), "matrix")
otu_table(ps_counts) = otu_table(if(taxa_are_rows(ps_counts)) otu_counts[, sample_order] else otu_counts[sample_order, ],
                                 taxa_are_rows = taxa_are_rows(ps_counts))
sample_data(ps_counts) = sample_data(sd)
write.table(data.frame(order = seq_along(sample_order), sample = sample_order, sd[, c("habitat", "treatment", "pot")]),
            file.path(outdir, paste0(key, ".sample_order.txt")), sep = "\t", row.names = F, quote = F)
treatment_colors = if(!is.null(acfg$treatments$colors)) unlist(acfg$treatments$colors) else
  setNames(RColorBrewer::brewer.pal(max(3, length(acfg$treatments$levels)), "Dark2")[seq_along(acfg$treatments$levels)],
           acfg$treatments$levels)
treatment_colors = c(treatment_colors[acfg$treatments$levels], setNames(rep("grey30", length(source_label)), source_label))

# ---- tree rooting ----
# FastTree trees are unrooted and phyloseq::UniFrac would then root them at a random tip (different distances
# on every run). Root on an outgroup instead (config tree_root.outgroup, a name at the first taxonomic rank; 16S:
# Archaea): on the branch between the outgroup ASVs (when they form a clade) and the rest. Rooting on a single
# outgroup tip would split the outgroup. Fallbacks: the most abundant outgroup ASV, else the midpoint.
has_tree = !is.null(phy_tree(ps, errorIfNULL = FALSE))
outgroup = if(!is.null(acfg$tree_root$outgroup)) acfg$tree_root$outgroup else "Archaea"
if(has_tree) {
  tt = as.data.frame(as(tax_table(ps), "matrix"))
  arch = taxa_names(ps)[tt[, 1] %in% outgroup]   # first rank (SILVA 144: Domain)
  tr = phy_tree(ps)
  if(length(arch) > 0 && length(arch) < ntaxa(ps) && is.monophyletic(tr, arch, reroot = TRUE)) {
    rooted = root(tr, outgroup = arch, resolve.root = TRUE)
    cat("tree rooted on the", outgroup, "clade:", length(arch), "ASVs\n")
  } else if(length(arch) > 0) {
    a1 = arch[which.max(taxa_sums(ps)[arch])]
    rooted = root(tr, outgroup = a1, resolve.root = TRUE)
    cat(outgroup, "not monophyletic: tree rooted on the most abundant", outgroup, "ASV", a1, "\n")
  } else {
    rooted = midpoint(tr)
    cat("no", outgroup, "ASVs: tree rooted at the midpoint\n")
  }
  phy_tree(ps) = rooted
  phy_tree(ps_counts) = rooted
}

# ---- 1. taxonomic composition ----
# stacked bars of relative abundance, one figure per taxonomic rank, samples in sample_order and split into
# habitat blocks (baseline with its habitat). Per rank: the ntop taxa with the highest mean relative abundance across samples,
# the rest summed as "Other"; ASVs without an assignment at that rank are summed as "Unassigned"
ntop = if(!is.null(acfg$taxa_barplots$ntop)) acfg$taxa_barplots$ntop else 15
rank_table = function(ps, rank) {
  tt = as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = F)
  taxon = sub("^[a-z]__", "", tt[[rank]])   # UNITE names carry k__/p__/... prefixes
  # SILVA species are epithets only: prefix the genus
  if(rank == "Species" && "Genus" %in% names(tt) && !any(grepl("_", taxon[!is.na(taxon)])))
    taxon = ifelse(is.na(taxon), NA, paste(sub("^[a-z]__", "", tt$Genus), taxon))
  taxon[is.na(taxon)] = "Unassigned"
  m = as(otu_table(ps), "matrix"); if(!taxa_are_rows(ps)) m = t(m)
  rowsum(m, taxon)[, sample_names(ps), drop = F]   # taxa x samples
}
# panel labels: the baseline goes with its habitat; the one-sample source panel gets the short label so it fits
block_map = c(setNames(baseline_habitat, baseline_label), setNames(source_short, source_label))
sample_block = as.character(sd$habitat)
sample_block = factor(ifelse(sample_block %in% names(block_map), block_map[sample_block], sample_block),
                      levels = c(acfg$habitats, source_short))
names(sample_block) = rownames(sd)
sample_label = if(!is.null(cfg$samples$label_strip)) sub(cfg$samples$label_strip, "", sample_names(ps)) else sample_names(ps)
names(sample_label) = sample_names(ps)
# x positions: one unit per sample, plus a space (treatment_gap) between treatments within a habitat panel; a
# vertical line after each T0 separates it from the pots of its treatment. Positions run on across habitat panels
# (not reset per panel), so each free-x panel only shows the labels of its own samples
treatment_gap = 0.15
n = nsamples(ps)
same_block = c(FALSE, sample_block[-1] == sample_block[-n])
new_treatment = same_block & c(FALSE, as.character(sd$treatment)[-1] != as.character(sd$treatment)[-n])
after_t0 = same_block & !new_treatment & c(FALSE, (sd$is_T0 %in% 1)[-n] & !(sd$is_T0 %in% 1)[-1])
# (+5 at each new habitat panel: panels never share x values, so widening the one-bar SMT panel below cannot pull
# a neighbouring panel's label into it)
sample_x = cumsum(1 + treatment_gap * new_treatment + 5 * c(FALSE, !same_block[-1]))
names(sample_x) = sample_names(ps)
group_lines = data.frame(x = sample_x[after_t0] - 0.5, block = sample_block[after_t0])
# treatment labels in a band under the bars (y < 0): a bracket line spanning the treatment's samples (T0 and pots)
# and the name below it, wrapped to fit its block (12 characters per line; 9 in blocks of <= 3 samples, e.g. 2
# rhizosphere pots of Conv_BC)
label_band = 0.2
treatment_labels = data.frame(block = sample_block, treatment = as.character(sd$treatment), x = sample_x) %>%
  group_by(block, treatment) %>%
  summarise(xmin = min(x) - 0.45, xmax = max(x) + 0.45, x = mean(x), n = n(), .groups = "drop") %>%
  filter(!(treatment %in% source_label)) %>%   # one-bar panel, already named by its strip
  mutate(label = ifelse(n <= 3, treatment_name(treatment, 9), treatment_name(treatment, 12)))
# the source panel holds one bar: widen it to 3 bar widths so its strip title fits
pellet_x = sample_x[sd$habitat %in% source_label]
panel_padding = data.frame(x = c(pellet_x - 1.5, pellet_x + 1.5), block = factor(rep(source_short, 2 * length(pellet_x)), levels = levels(sample_block)))
ranks = rank_names(ps)
if(write_pdf) pdf(file.path(outdir, paste0(key, ".taxa_barplots.pdf")), width = 16, height = 7)
for(rank in ranks) {
  rt = rank_table(ps, rank)
  write.table(rownames_to_column(as.data.frame(rt), rank), file.path(outdir, paste0(key, ".relabund.", rank, ".txt")),
              sep = "\t", row.names = F, quote = F)
  assigned = setdiff(rownames(rt), "Unassigned")
  top = head(assigned[order(-rowMeans(rt[assigned, , drop = F]))], ntop)
  other = setdiff(assigned, top)
  rt2 = rbind(rt[top, , drop = F],
              Other = colSums(rt[other, , drop = F]),
              Unassigned = if("Unassigned" %in% rownames(rt)) rt["Unassigned", ] else 0)
  long = as.data.frame(rt2) %>% rownames_to_column("taxon") %>%
    tidyr::pivot_longer(-taxon, names_to = "sample", values_to = "relab") %>%
    mutate(taxon = factor(taxon, levels = rev(c(top, "Other", "Unassigned"))),
           x = sample_x[sample],
           block = sample_block[as.character(sample)])
  pal = c(setNames(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(top)), top),
          Other = "grey75", Unassigned = "grey40")
  p = ggplot(long, aes(x, relab, fill = taxon)) +
    geom_col(width = 1) +
    geom_blank(data = panel_padding, aes(x = x), inherit.aes = F) +
    geom_segment(data = group_lines, aes(x = x, xend = x, y = 0, yend = 1), inherit.aes = F,
                 color = "black", linewidth = 0.15) +
    geom_segment(data = treatment_labels, aes(x = xmin, xend = xmax, y = -0.012, yend = -0.012), inherit.aes = F,
                 linewidth = 0.3) +
    geom_text(data = treatment_labels, aes(x = x, y = -0.02, label = label), inherit.aes = F,
              vjust = 1, size = 3, lineheight = 0.85) +   # ~8.5 pt; habitat strips are 12 pt
    facet_grid(. ~ block, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal, breaks = c(top, "Other", "Unassigned")) +
    scale_x_continuous(breaks = sample_x, labels = sample_label, expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    # zoom (not scale limits, which would drop stacked bars ending a rounding error above 1) to make room for labels
    coord_cartesian(ylim = c(-label_band, 1)) +
    labs(title = paste(locus, "relative abundance:", rank),
         subtitle = sprintf("top %d of %d %s by mean relative abundance", length(top), length(assigned), tolower(rank)),
         x = NULL, y = "relative abundance", fill = rank) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5), panel.grid = element_blank(),
          panel.spacing.x = unit(0.6, "cm"),
          strip.text.x = element_text(size = 12, face = "bold"),
          legend.key.size = unit(0.4, "cm"))
  if(write_pdf) print(p)
  ggsave(file.path(outdir, paste0(key, ".taxa_barplot.", rank, ".png")), p, width = 16, height = 7, dpi = 130)
}
if(write_pdf) invisible(dev.off())
cat("taxa barplots:", paste(ranks, collapse = ", "), "->", file.path(outdir, paste0(key, ".taxa_barplot.<rank>.png")), "\n")

# ---- 2. alpha diversity ----
# from read counts (ps_counts), rarefied once to the smallest library (seed 100) so richness is compared at equal
# depth. Metrics: observed ASVs, Shannon, inverse Simpson, Pielou evenness (Shannon / log(observed)) and, with a
# tree, Faith's phylogenetic diversity (sum of branch lengths connecting a sample's ASVs to the root).
# DADA2 removes singletons, so singleton-based estimators (Chao1, ACE) are not used.
depth = min(sample_sums(ps_counts))
ps_rare = rarefy_even_depth(ps_counts, sample.size = depth, rngseed = 100, replace = FALSE, trimOTUs = FALSE, verbose = FALSE)
cat("alpha diversity: rarefied to", depth, "reads per sample\n")
alpha_div = estimate_richness(ps_rare, measures = c("Observed", "Shannon", "InvSimpson"))
rownames(alpha_div) = sample_names(ps_rare)
alpha_div$Pielou = alpha_div$Shannon / log(alpha_div$Observed)
faith_pd = function(ps) {
  tr = phy_tree(ps)
  m = as(otu_table(ps), "matrix"); if(taxa_are_rows(ps)) m = t(m)
  m = m[, tr$tip.label, drop = F] > 0
  tips_below = phangorn::Descendants(tr, tr$edge[, 2], type = "tips")   # tips below each edge
  apply(m, 1, function(present) sum(tr$edge.length[vapply(tips_below, function(t) any(present[t]), logical(1))]))
}
if(has_tree) alpha_div$FaithPD = faith_pd(ps_rare)[rownames(alpha_div)]
metric_names = c(Observed = "Observed ASVs", Shannon = "Shannon", InvSimpson = "Inverse Simpson",
                 Pielou = "Pielou evenness", FaithPD = "Faith's PD")
metric_names = metric_names[names(metric_names) %in% names(alpha_div)]
write.table(rownames_to_column(cbind(sd[rownames(alpha_div), c("habitat", "treatment", "pot")], alpha_div), "sample"),
            file.path(outdir, paste0(key, ".alpha_diversity.txt")), sep = "\t", row.names = F, quote = F)

# one group per treatment (pots) and per treatment T0, within each habitat panel, in sample order
alpha_long = rownames_to_column(alpha_div, "sample") %>%
  tidyr::pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric_names[metric], levels = metric_names),
         block = sample_block[sample],
         treatment = sd[sample, "treatment"],
         is_T0 = sd[sample, "is_T0"] %in% 1,
         group = ifelse(sd[sample, "habitat"] %in% source_label, as.character(sd[sample, "habitat"]),
                        ifelse(is_T0, paste(treatment, baseline_label), as.character(treatment))))
group_levels = unique(paste(sample_block, alpha_long$group[match(sample_names(ps), alpha_long$sample)], sep = "|"))
alpha_long$x = factor(paste(alpha_long$block, alpha_long$group, sep = "|"), levels = group_levels)
p = ggplot(alpha_long, aes(x, value, color = treatment)) +
  # all groups first: the x scale takes its order from the first layer, and the boxes below leave out T0
  geom_blank() +
  geom_boxplot(data = filter(alpha_long, !is_T0, !(block %in% source_short)), aes(fill = treatment), alpha = 0.15,
               outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_point(aes(shape = is_T0), size = 1.6) +
  facet_grid(metric ~ block, scales = "free", space = "free_x", switch = "y") +
  scale_x_discrete(labels = group_axis_labels) +
  scale_color_manual(values = treatment_colors, labels = treatment_axis_labels, drop = FALSE) +
  scale_fill_manual(values = treatment_colors, guide = "none") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8), labels = c("FALSE" = "pot", "TRUE" = baseline_label), name = NULL) +
  labs(title = paste(locus, "alpha diversity"),
       subtitle = sprintf("rarefied to %d reads per sample; boxes: pots A-D per treatment", depth),
       x = NULL, y = NULL, color = "treatment") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), panel.grid.major.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"), strip.text.y = element_text(size = 9),
        strip.placement = "outside", strip.background.y = element_blank(), panel.spacing.x = unit(0.6, "cm"))
save_plot(p, "alpha_diversity", width = 15, height = 2.2 * length(metric_names) + 1.5)

# ---- 3. ASV overlap among a set of treatments (Venn), per habitat (config venn) ----
# sets: venn.sets (pots only, no baseline). On the rarefied counts (ps_rare, equal depth), an ASV belongs to a treatment
# if it is detected in at least half of that treatment's pots (2 of 4, 2 of 3, 1 of 2), so single-pot detections of
# rare ASVs do not count as presence. Pairwise diagrams: venn.reference vs each other set.
rare_m = as(otu_table(ps_rare), "matrix"); if(taxa_are_rows(ps_rare)) rare_m = t(rare_m)
if(!is.null(acfg$venn)) {
venn_treatments = acfg$venn$sets
venn_habitats = acfg$venn$habitats
venn_ref = acfg$venn$reference
venn_name = paste0("venn_", acfg$venn$name)
venn_plots = list(); venn_membership = list(); venn_summary = list(); venn_sets = list(); venn_npots = list()
for(h in venn_habitats) {
  sets = list(); npots = c()
  for(trt in venn_treatments) {
    pots = rownames(sd)[sd$habitat %in% h & sd$treatment %in% trt & !(sd$is_T0 %in% 1)]
    npots[trt] = length(pots)
    detected = colSums(rare_m[pots, , drop = F] > 0)
    sets[[trt]] = colnames(rare_m)[detected >= ceiling(length(pots) / 2)]
  }
  venn_sets[[h]] = sets; venn_npots[[h]] = npots
  set_labels = sprintf("%s\n(%d pots, %d ASVs)", treatment_name(names(sets), 22), npots, lengths(sets))
  venn_plots[[h]] = ggvenn::ggvenn(setNames(sets, set_labels), fill_color = unname(treatment_colors[venn_treatments]),
                                   fill_alpha = 0.35, stroke_size = 0.4, set_name_size = 3.2, text_size = 3.2,
                                   show_percentage = TRUE, digits = 0) +
    ggtitle(h) + theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  # per-ASV membership with its region, mean relative abundance in the habitat's pots of the sets and taxonomy
  all_asv = unique(unlist(sets))
  memb = sapply(sets, function(x) all_asv %in% x); rownames(memb) = all_asv
  region = apply(memb, 1, function(r) paste(venn_treatments[r], collapse = " & "))
  set_pots = rownames(sd)[sd$habitat %in% h & sd$treatment %in% venn_treatments & !(sd$is_T0 %in% 1)]
  rel_m = as(otu_table(ps), "matrix"); if(taxa_are_rows(ps)) rel_m = t(rel_m)
  mean_relab = colMeans(rel_m[set_pots, all_asv, drop = F])
  tt_venn = as.data.frame(as(tax_table(ps), "matrix"))[all_asv, , drop = F]
  venn_membership[[h]] = data.frame(habitat = h, asv = all_asv, region = region, n_sets = rowSums(memb), memb,
                                    mean_relab = mean_relab, tt_venn, check.names = F)
  total_relab = sum(mean_relab)   # share of the set pots' reads carried by ASVs in any set
  venn_summary[[h]] = venn_membership[[h]] %>% group_by(habitat, region, n_sets) %>%
    summarise(n_asv = n(), pct_asv = 100 * n() / length(all_asv), pct_relab = 100 * sum(mean_relab) / total_relab,
              .groups = "drop") %>%
    arrange(desc(n_sets), desc(n_asv))
  cat(sprintf("Venn %s: ASVs in any set carry %.1f%% of the mean relative abundance of the set pots\n", h, 100 * total_relab))
}
write.table(do.call(rbind, venn_membership), file.path(outdir, paste0(key, ".", venn_name, ".asv_membership.txt")),
            sep = "\t", row.names = F, quote = F)
venn_summary = do.call(rbind, venn_summary)
write.table(venn_summary, file.path(outdir, paste0(key, ".", venn_name, ".regions.txt")), sep = "\t", row.names = F, quote = F)
for(h in venn_habitats) {
  core = venn_summary %>% filter(habitat == h, n_sets == length(venn_treatments))
  cat(sprintf("Venn %s: %d ASVs shared by all %d sets (%.0f%% of ASVs, %.0f%% of mean relative abundance)\n",
              h, sum(core$n_asv), length(venn_treatments), sum(core$pct_asv), sum(core$pct_relab)))
}
p = patchwork::wrap_plots(venn_plots, nrow = 1) +
  patchwork::plot_annotation(title = paste(locus, "ASV overlap:", acfg$venn$title),
                             subtitle = "ASV counted in a treatment if detected in >= half of its pots (rarefied counts); % of ASVs in the union")
save_plot(p, venn_name, width = 14, height = 7.5)

# pairwise: the reference set vs each other set (same sets as the n-way diagrams), rows = habitat.
# Drawn directly (ggvenn ignores fill_color for 2 sets): two circles, counts and % of the union in each region
venn2 = function(a, b, labels, colors, title) {
  circle = function(x0, id) { t = seq(0, 2 * pi, length.out = 200); data.frame(x = x0 + cos(t), y = sin(t), id = id) }
  circles = rbind(circle(-0.55, "a"), circle(0.55, "a2"))
  n = c(length(setdiff(a, b)), length(intersect(a, b)), length(setdiff(b, a)))
  counts = data.frame(x = c(-1.05, 0, 1.05), y = 0, label = sprintf("%d\n(%.0f%%)", n, 100 * n / sum(n)))
  ggplot() +
    geom_polygon(data = circles, aes(x, y, group = id, fill = id), alpha = 0.35, color = "black", linewidth = 0.4) +
    scale_fill_manual(values = setNames(colors, c("a", "a2")), guide = "none") +
    geom_text(data = counts, aes(x, y, label = label), size = 3.4, lineheight = 0.9) +
    annotate("text", x = c(-0.9, 0.9), y = 1.3, label = labels, size = 3.2, lineheight = 0.9, vjust = 0) +
    coord_equal(xlim = c(-1.8, 1.8), ylim = c(-1.1, 1.75), clip = "off") +
    ggtitle(title) + theme_void() + theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
}
rel_all = as(otu_table(ps), "matrix"); if(taxa_are_rows(ps)) rel_all = t(rel_all)
pair_plots = list(); pair_summary = list()
for(h in venn_habitats) for(trt in setdiff(venn_treatments, venn_ref)) {
  pair = venn_sets[[h]][c(venn_ref, trt)]
  np = venn_npots[[h]][c(venn_ref, trt)]
  pair_plots[[paste(h, trt)]] = venn2(pair[[1]], pair[[2]],
                                      labels = sprintf("%s\n(%d pots, %d ASVs)", treatment_name(names(pair), 22), np, lengths(pair)),
                                      colors = unname(treatment_colors[c(venn_ref, trt)]),
                                      title = h)
  # share of each treatment's reads (mean relative abundance over its pots) carried by the shared ASVs
  shared = intersect(pair[[1]], pair[[2]])
  read_share = sapply(names(pair), function(t) {
    pots = rownames(sd)[sd$habitat %in% h & sd$treatment %in% t & !(sd$is_T0 %in% 1)]
    m = colMeans(rel_all[pots, , drop = F]); 100 * sum(m[shared]) / sum(m)
  })
  pair_summary[[paste(h, trt)]] = data.frame(habitat = h, reference = venn_ref, compared = trt,
    n_reference = length(pair[[1]]), n_compared = length(pair[[2]]), n_shared = length(shared),
    only_reference = length(setdiff(pair[[1]], pair[[2]])), only_compared = length(setdiff(pair[[2]], pair[[1]])),
    jaccard = length(shared) / length(union(pair[[1]], pair[[2]])),
    pct_reads_shared_reference = read_share[[1]], pct_reads_shared_compared = read_share[[2]])
}
pair_summary = do.call(rbind, pair_summary)
write.table(pair_summary, file.path(outdir, paste0(key, ".", venn_name, "_pairwise.txt")), sep = "\t", row.names = F, quote = F)
print(pair_summary[, c("habitat", "compared", "n_shared", "only_reference", "only_compared", "jaccard",
                       "pct_reads_shared_reference", "pct_reads_shared_compared")], row.names = F, digits = 3)
p = patchwork::wrap_plots(pair_plots, nrow = length(venn_habitats)) +
  patchwork::plot_annotation(title = paste(locus, "ASV overlap:", acfg$venn$pairwise_title),
                             subtitle = "ASV counted in a treatment if detected in >= half of its pots (rarefied counts); % of ASVs in the union")
save_plot(p, paste0(venn_name, "_pairwise"), width = 15, height = 10)
}

# ---- 4. NMDS ----
# with a tree: weighted UniFrac (rooted tree); without (e.g. ITS): Bray-Curtis. Both on relative abundances (ps).
# One ordination per config ordination subset (a set of habitats; e.g. bulk + rhizosphere together and root alone,
# when root and soil communities are too far apart to share one NMDS). The baseline (T0) is left out.
# The source sample is ordinated with the subset named in source.ordinate_with but is not in the PERMANOVA or in the
# distances passed to the downstream sections (pots only).
dist_method = if(has_tree) "wunifrac" else "bray"
dist_name = c(wunifrac = "weighted UniFrac", bray = "Bray-Curtis")[[dist_method]]
# also an NMDS (with PERMANOVA) on a presence/absence distance: 16S unweighted UniFrac (presence/absence of lineages),
# ITS binary Jaccard (presence/absence of ASVs; no tree); figures and tables only, the downstream sections use dist_method
pa_method = if(has_tree) "unifrac" else "jaccard"
nmds_methods = c(dist_method, pa_method)
nmds_names = c(wunifrac = "weighted UniFrac", unifrac = "unweighted UniFrac", bray = "Bray-Curtis",
               jaccard = "Jaccard (presence/absence)")
calc_distance = function(physeq, method)
  if(method == "jaccard") phyloseq::distance(physeq, method = "jaccard", binary = TRUE) else phyloseq::distance(physeq, method = method)
nmds_subsets = lapply(acfg$ordination$subsets, as.character)
habitat_shapes = setNames(c(16, 17, 15, 18, 3, 4)[seq_along(acfg$habitats)], acfg$habitats)
habitat_shapes = c(habitat_shapes, setNames(rep(8, length(source_label)), source_label))
nmds_pellet = rownames(sd)[sd$is_pellet %in% 1]
source_subset = if(has_source && !is.null(acfg$source$ordinate_with)) acfg$source$ordinate_with else NA
# PERMANOVA (section 5) runs here too, on the same distance, so each NMDS figure can carry its table.
# formula per subset from the config (e.g. bulk + rhizosphere: habitat * treatment; root: treatment). Sequential
# (type I) sums of squares (by = "terms"), 999 permutations, seed 100; R2 = share of the total variance explained by
# each term. Note: habitats sampled from the same pots are not independent, but permutations are unrestricted
# (pots are not used as strata), so p-values treat the samples as independent.
adonis_formulas = lapply(acfg$ordination$formulas, function(rhs) as.formula(paste("d ~", rhs)))
run_adonis = function(d, formula) {
  meta = sd[labels(d), , drop = F] %>% mutate(habitat = droplevels(habitat), treatment = droplevels(treatment))
  environment(formula) = environment()
  set.seed(100)
  as.data.frame(adonis2(formula, data = meta, permutations = 999, by = "terms")) %>% rownames_to_column("term")
}
# table grob for the NMDS inset: term, df, R2 (%), F, p (Total row dropped)
adonis_grob = function(res) {
  tab = res %>% filter(term != "Total") %>%
    transmute(term = dplyr::recode(term, "habitat:treatment" = "habitat x treatment", "Residual" = "residual"),
              df = Df, `R2 (%)` = sprintf("%.1f", 100 * R2), `F` = ifelse(is.na(`F`), "", sprintf("%.2f", `F`)),
              p = ifelse(is.na(`Pr(>F)`), "", sprintf("%.3f", `Pr(>F)`)))
  gridExtra::tableGrob(tab, rows = NULL, theme = gridExtra::ttheme_minimal(
    base_size = 8, padding = unit(c(3, 2), "mm"),
    core = list(bg_params = list(fill = "white", alpha = 0.8)),
    colhead = list(bg_params = list(fill = "grey90", alpha = 0.9), fg_params = list(fontface = "bold"))))
}
# corner of the panel with the fewest points (scores scaled to 0-1); returns inset_element coordinates
inset_corner = function(x, y, w = 0.46, h = 0.3) {
  sx = (x - min(x)) / diff(range(x)); sy = (y - min(y)) / diff(range(y))
  corners = list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))   # bottom-left, bottom-right, top-left, top-right
  n = sapply(corners, function(cn) sum(abs(sx - cn[1]) < w & abs(sy - cn[2]) < h))
  cn = corners[[which.min(n)]]
  left = if(cn[1] == 0) 0.01 else 0.99 - w; bottom = if(cn[2] == 0) 0.01 else 0.99 - h
  c(left = left, bottom = bottom, right = left + w, top = bottom + h)
}
nmds_dist = list(); adonis_results = list(); adonis_extra = list(); nmds_plots = list()
for(nm in nmds_methods) for(subset in names(nmds_subsets)) {
  nm_name = nmds_names[[nm]]
  sn = rownames(sd)[sd$habitat %in% nmds_subsets[[subset]] & !(sd$is_T0 %in% 1)]
  sn_ord = if(identical(subset, source_subset)) c(sn, nmds_pellet) else sn   # samples in the ordination
  ps_sub = prune_samples(sn_ord, ps)
  ps_sub = prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  d_ord = calc_distance(ps_sub, nm)
  d = as.dist(as.matrix(d_ord)[sn, sn])   # pots only (pairwise distances, so unchanged by the pellet)
  saveRDS(d, file.path(outdir, paste0(key, ".", nm, ".", subset, ".dist.rds")))
  res = run_adonis(d, adonis_formulas[[subset]]) %>%
    mutate(subset = subset, distance = nm_name, n_samples = attr(d, "Size"), .before = 1)
  if(nm == dist_method) { nmds_dist[[subset]] = d; adonis_results[[subset]] = res } else adonis_extra[[paste(nm, subset)]] = res
  cat(sprintf("PERMANOVA (%s), %s:\n", nm_name, subset))
  print(res[, c("term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")], row.names = F, digits = 3)
  nmds = metaMDS(d_ord, k = 2, try = 100, trymax = 500, trace = 0)
  cat(sprintf("NMDS (%s), %s: %d samples, stress %.3f, converged: %s\n", nm_name, subset, length(sn_ord), nmds$stress,
              nmds$converged))
  scores = as.data.frame(vegan::scores(nmds, display = "sites")) %>%
    rownames_to_column("sample") %>%
    left_join(rownames_to_column(sd, "sample"), by = "sample") %>%
    mutate(habitat = droplevels(habitat))
  write.table(scores, file.path(outdir, paste0(key, ".nmds.", nm, ".", subset, ".txt")), sep = "\t", row.names = F, quote = F)

  # ellipses (vegan::ordiellipse, kind = "sd"): one standard deviation of the group's site scores around its
  # centroid, per treatment x habitat; groups with fewer than 3 pots get none
  scores$group = paste(scores$treatment, scores$habitat, sep = "|")
  big = names(which(table(scores$group) >= 3))
  in_big = scores$group %in% big
  ell = vegan::ordiellipse(as.matrix(scores[in_big, c("NMDS1", "NMDS2")]), scores$group[in_big], kind = "sd", draw = "none")
  ellipses = do.call(rbind, lapply(names(ell), function(g) {
    xy = vegan:::veganCovEllipse(ell[[g]]$cov, ell[[g]]$center, ell[[g]]$scale)
    data.frame(NMDS1 = xy[, 1], NMDS2 = xy[, 2], group = g,
               treatment = factor(sub("[|].*", "", g), levels = levels(scores$treatment)),
               habitat = factor(sub(".*[|]", "", g), levels = levels(scores$habitat)))
  }))
  p = ggplot(scores, aes(NMDS1, NMDS2, color = treatment, shape = habitat)) +
    geom_polygon(data = ellipses, mapping = aes(NMDS1, NMDS2, group = group, color = treatment, fill = treatment,
                                                linetype = habitat), inherit.aes = F, alpha = 0.08, linewidth = 0.5) +
    geom_point(size = 3) +
    scale_fill_manual(values = treatment_colors, guide = "none") +
    # ellipse line type per habitat, in subset order (solid, dashed, ...)
    scale_linetype_manual(values = setNames(c("solid", "dashed", "dotted", "dotdash")[seq_along(nmds_subsets[[subset]])],
                                            nmds_subsets[[subset]]),
                          guide = if(length(unique(scores$habitat)) > 1) "legend" else "none") +
    scale_color_manual(values = treatment_colors, labels = treatment_axis_labels, drop = TRUE) +
    scale_shape_manual(values = habitat_shapes, drop = TRUE) +
    coord_equal() +
    labs(title = sprintf("%s NMDS, %s: %s", locus, nm_name, gsub("_", " + ", subset)),
         subtitle = sprintf("%d pots%s, %d ASVs, stress %.3f; ellipses: 1 SD per treatment x habitat",
                            length(sn), if(length(sn_ord) > length(sn)) paste(" +", source_label) else "", ntaxa(ps_sub), nmds$stress)) +
    theme_bw()
  # PERMANOVA table in the emptiest corner of the panel
  h = 0.07 + 0.055 * (nrow(res) - 1)
  pos = inset_corner(scores$NMDS1, scores$NMDS2, h = h)
  p = p + patchwork::inset_element(adonis_grob(res), left = pos[["left"]], bottom = pos[["bottom"]],
                                   right = pos[["right"]], top = pos[["top"]], align_to = "panel", clip = FALSE)
  save_plot(p, paste0("nmds.", nm, ".", subset), width = 8, height = 6, dpi = 150)

  # an ordination of several habitats drawn one habitat per panel, on the axes of the whole ordination
  if(length(nmds_subsets[[subset]]) > 1) {
    xr = range(scores$NMDS1, ellipses$NMDS1); yr = range(scores$NMDS2, ellipses$NMDS2)
    for(hb in nmds_subsets[[subset]]) {
      sc = scores[scores$habitat == hb, ]
      pel = scores[scores$sample %in% nmds_pellet, ]   # source sample: shown in every habitat panel
      q = ggplot(sc, aes(NMDS1, NMDS2, color = treatment)) +
        geom_polygon(data = ellipses[ellipses$habitat == hb, ], mapping = aes(NMDS1, NMDS2, group = group, color = treatment,
                     fill = treatment), inherit.aes = F, alpha = 0.08, linewidth = 0.5, show.legend = FALSE) +
        geom_point(size = 5, shape = 16) +
        geom_point(data = pel, aes(NMDS1, NMDS2), inherit.aes = F, shape = 16, size = 6, color = "grey25") +
        geom_text(data = pel, aes(NMDS1, NMDS2), inherit.aes = F, label = source_short[seq_len(nrow(pel) > 0)], size = 6, fontface = "bold",
                  vjust = -1.1) +
        scale_fill_manual(values = treatment_colors, guide = "none") +
        scale_color_manual(values = treatment_colors, labels = treatment_axis_labels, drop = FALSE,
                           limits = setdiff(levels(droplevels(scores$treatment)), source_label)) +
        coord_equal(xlim = xr, ylim = yr) +
        guides(color = guide_legend(override.aes = list(size = 7))) +
        theme_bw() + theme(legend.text = element_text(size = 18), legend.title = element_text(size = 20),
                           legend.key.size = unit(1.3, "cm"), axis.title = element_text(size = 18),
                           axis.text = element_text(size = 14))
      nmds_plots[[subset]][[nm]][[hb]] = q
    }
  }
}
# multi-habitat ordinations by habitat: rows = distance (with a tree: weighted, unweighted UniFrac; without:
# Bray-Curtis, Jaccard), columns = habitat, with large column (habitat) and row (distance) headers
header = function(label, rot = 0) patchwork::wrap_elements(full = grid::textGrob(label, rot = rot,
  gp = grid::gpar(fontsize = 26, fontface = "bold")))
for(subset in names(nmds_plots)) {
  hb_cols = nmds_subsets[[subset]]
  grid_cells = c(list(patchwork::plot_spacer()), lapply(tools::toTitleCase(hb_cols), header))
  for(nm in nmds_methods) grid_cells = c(grid_cells, list(header(tools::toTitleCase(nmds_names[[nm]]), rot = 90)),
                                         unname(nmds_plots[[subset]][[nm]][hb_cols]))
  p = patchwork::wrap_plots(grid_cells, ncol = length(hb_cols) + 1, byrow = TRUE,
                            widths = c(0.1, rep(1, length(hb_cols))), heights = c(0.1, rep(1, length(nmds_methods)))) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = sprintf("%s NMDS of %s pots: one ordination per distance, shown by habitat", locus,
                                               paste(hb_cols, collapse = " + ")),
                               subtitle = "ellipses: 1 SD per treatment x habitat")
  save_plot(p, paste0("nmds.", paste(nmds_methods, collapse = "_"), ".", subset, ".by_habitat"),
            width = 8 * length(hb_cols), height = 5.5 * length(nmds_methods) + 1.5, dpi = 150)
}
# PERMANOVA tables of the presence/absence NMDS distances (unweighted UniFrac or Jaccard); sections 5-6 and 8 use dist_method
if(length(adonis_extra) > 0) {
  adonis_extra = do.call(rbind, adonis_extra)
  for(nm in unique(names(nmds_names)[match(adonis_extra$distance, nmds_names)]))
    write.table(adonis_extra[adonis_extra$distance == nmds_names[[nm]], ], file.path(outdir, paste0(key, ".permanova.", nm, ".txt")),
                sep = "\t", row.names = F, quote = F)
}

# ---- 5. variance partitioning (PERMANOVA, vegan::adonis2) ----
# computed with the NMDS above (same distances and samples; settings there); here: table and bar plot
adonis_results = do.call(rbind, adonis_results)
write.table(adonis_results, file.path(outdir, paste0(key, ".permanova.", dist_method, ".txt")), sep = "\t", row.names = F, quote = F)

# variance partition plot: R2 per term (residual in grey), stars for p <= 0.05 / 0.01 / 0.001 (with 999
# permutations the smallest possible p is 0.001)
vp = adonis_results %>% filter(term != "Total") %>%
  mutate(term = dplyr::recode(term, "habitat:treatment" = "habitat x treatment", "Residual" = "residual"),
         term = factor(term, levels = c("habitat", "treatment", "habitat x treatment", "residual")),
         subset = factor(gsub("_", " + ", subset), levels = gsub("_", " + ", names(adonis_formulas))),
         stars = ifelse(is.na(`Pr(>F)`), "", ifelse(`Pr(>F)` <= 0.001, "***", ifelse(`Pr(>F)` <= 0.01, "**",
                        ifelse(`Pr(>F)` <= 0.05, "*", "ns")))),
         label = ifelse(term == "residual", sprintf("%.0f%%", 100 * R2), sprintf("%.1f%% %s", 100 * R2, stars)))
p = ggplot(vp, aes(subset, R2, fill = term)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.3) +
  scale_fill_manual(values = c("habitat" = "#8c6d31", "treatment" = "#3182bd", "habitat x treatment" = "#9ecae1",
                               "residual" = "grey85")) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(title = paste(locus, "variance partitioning (PERMANOVA),", dist_name),
       subtitle = "R2 per term; 999 permutations; *** p<=0.001, ** p<=0.01, * p<=0.05",
       x = NULL, y = "variance explained (R2)", fill = NULL) +
  theme_bw() + theme(panel.grid.major.x = element_blank())
save_plot(p, paste0("permanova.", dist_method), width = 7, height = 6)

# ---- 6. pairwise PERMANOVA between treatments, within each habitat ----
# for each habitat and each pair of treatments: adonis2(distance ~ treatment) on those pots only (same distances
# as the NMDS; baseline excluded), 999 permutations, seed 100; p adjusted (Benjamini-Hochberg) over the pairs of a
# habitat. With 4 vs 4 pots there are only 35 distinct ways to split the samples, so the smallest attainable
# p is ~0.03 (2 vs 4 pots: ~0.07), and few pairs can stay significant after adjustment.
# each habitat's distances come from the ordination subset that contains it
pair_habitats = as.list(setNames(rep(names(nmds_subsets), lengths(nmds_subsets)), unlist(nmds_subsets)))
pairwise = list()
for(h in names(pair_habitats)) {
  dm = as.matrix(nmds_dist[[pair_habitats[[h]]]])
  pots = rownames(sd)[sd$habitat %in% h & !(sd$is_T0 %in% 1) & rownames(sd) %in% rownames(dm)]
  trts = intersect(treatment_levels, as.character(sd[pots, "treatment"]))
  for(pr in combn(trts, 2, simplify = FALSE)) {
    sn = pots[sd[pots, "treatment"] %in% pr]
    d = as.dist(dm[sn, sn])
    res = run_adonis(d, d ~ treatment)
    n = table(as.character(sd[sn, "treatment"]))
    pairwise[[length(pairwise) + 1]] = data.frame(habitat = h, treatment1 = pr[1], treatment2 = pr[2],
      n1 = as.integer(n[pr[1]]), n2 = as.integer(n[pr[2]]), R2 = res$R2[1], F = res$F[1], p = res$`Pr(>F)`[1])
  }
}
pairwise = do.call(rbind, pairwise) %>% group_by(habitat) %>% mutate(p_adj = p.adjust(p, method = "BH")) %>% ungroup() %>%
  mutate(distance = dist_name, .before = 1)
write.table(pairwise, file.path(outdir, paste0(key, ".permanova_pairwise.", dist_method, ".txt")), sep = "\t", row.names = F, quote = F)
cat(sprintf("pairwise PERMANOVA (%s): %d tests; p <= 0.05: %d, p_adj <= 0.05: %d\n", dist_name, nrow(pairwise),
            sum(pairwise$p <= 0.05), sum(pairwise$p_adj <= 0.05)))

# heatmap per habitat: R2 of each pair (lower triangle), label R2 and p (p_adj in brackets), bold outline p <= 0.05
hm = pairwise %>%
  mutate(habitat = factor(habitat, levels = names(pair_habitats)),
         t1 = factor(treatment1, levels = treatment_levels), t2 = factor(treatment2, levels = rev(treatment_levels)),
         label = sprintf("%.2f\np=%.3f\n(%.3f)", R2, p, p_adj), sig = p <= 0.05)
p = ggplot(hm, aes(t1, t2, fill = R2)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_tile(data = filter(hm, sig), color = "black", linewidth = 0.7, fill = NA) +
  geom_text(aes(label = label), size = 2.4, lineheight = 0.9) +
  facet_wrap(~ habitat, nrow = 1, drop = TRUE) +
  scale_x_discrete(labels = treatment_axis_labels) + scale_y_discrete(labels = treatment_axis_labels) +
  scale_fill_gradient(low = "white", high = "#3182bd", limits = c(0, max(hm$R2)), name = "R2") +
  labs(title = paste(locus, "pairwise PERMANOVA between treatments,", dist_name),
       subtitle = "cell: R2, p, (BH-adjusted p within habitat); black outline: p <= 0.05; 999 permutations",
       x = NULL, y = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
                     strip.text = element_text(size = 12, face = "bold"))
save_plot(p, paste0("permanova_pairwise.", dist_method), width = 15, height = 5.5)


# ---- rarefaction depth per habitat block (sections 7 and 9) ----
# each habitat block (with its baseline samples) and the source sample, if any, are rarefied together to the smallest
# library in that block (or the source's, if smaller), so soil samples are not cut to the depth of the shallowest root
# library. Repeated n_rare times (seeds 100, 101, ...) and averaged.
n_rare = c(acfg$source_sharing$n_rare, acfg$betapart$n_rare, 100)[1]
pellet = rownames(sd)[sd$is_pellet %in% 1]
counts_m = as(otu_table(ps_counts), "matrix"); if(taxa_are_rows(ps_counts)) counts_m = t(counts_m)
rel_m = as(otu_table(ps), "matrix"); if(taxa_are_rows(ps)) rel_m = t(rel_m)
smt_samples = rownames(sd)[!(sd$is_pellet %in% 1)]
smt_blocks = intersect(levels(sample_block), as.character(sample_block[smt_samples]))
smt_depth = sapply(smt_blocks, function(b) min(rowSums(counts_m[c(smt_samples[sample_block[smt_samples] == b], pellet), , drop = F])))

# ---- 7. source establishment: ASVs shared with the source sample (config source_sharing) ----
# for every treatment sample (pots and baseline) the number of ASVs it shares with the source sample (e.g. an
# inoculum). Treatments without the source give the background: source ASVs already in the soil.
# Rarefied per habitat block as above, n_rare times, averaged per sample, so the result does not depend on a single
# random subsample.
# Per sample: ASVs (sample_asvs), ASVs shared with the source (shared_asvs), shared as % of the source's ASVs and of
# the sample's ASVs, and shared_relab: % of the sample's reads (unrarefied, ps) in ASVs detected in the rarefied source.
if(has_source && !is.null(acfg$source_sharing)) {
stopifnot(length(pellet) == 1)
source_name = acfg$source$name
source_noun = if(!is.null(acfg$source$noun)) acfg$source$noun else "source"   # e.g. "pellet"
smt_pellet_asvs = c(); smt_shared = list()
for(b in smt_blocks) {
  sn = smt_samples[sample_block[smt_samples] == b]
  res = lapply(seq_len(n_rare), function(i) {
    set.seed(99 + i)
    r = vegan::rrarefy(counts_m[c(sn, pellet), ], smt_depth[b])
    in_pellet = r[pellet, ] > 0
    cbind(pellet_asvs = sum(in_pellet),
          sample_asvs = rowSums(r[sn, , drop = F] > 0),
          shared_asvs = rowSums(r[sn, in_pellet, drop = F] > 0),
          shared_relab = 100 * rowSums(rel_m[sn, in_pellet, drop = F]))
  })
  avg = Reduce(`+`, res) / n_rare   # samples x metrics, mean over the rarefactions
  smt_pellet_asvs[b] = avg[1, "pellet_asvs"]
  smt_shared[[b]] = data.frame(sample = sn, sd[sn, c("habitat", "treatment", "pot")], is_T0 = sd[sn, "is_T0"] %in% 1,
    block = b, depth = smt_depth[b], pellet_asvs = avg[, "pellet_asvs"], sample_asvs = avg[, "sample_asvs"],
    shared_asvs = avg[, "shared_asvs"],
    shared_asvs_sd = apply(sapply(res, function(x) x[, "shared_asvs"]), 1, stats::sd),   # spread over rarefactions
    shared_relab = avg[, "shared_relab"], row.names = NULL)
  cat(sprintf("%s %s, %s: rarefied %d x to %d reads; %s %.0f ASVs (mean)\n", source_label, pellet, b, n_rare,
              smt_depth[b], source_noun, smt_pellet_asvs[b]))
}
smt_shared = do.call(rbind, smt_shared)
smt_shared$pct_pellet_asvs = 100 * smt_shared$shared_asvs / smt_shared$pellet_asvs
smt_shared$pct_sample_asvs = 100 * smt_shared$shared_asvs / smt_shared$sample_asvs
write.table(smt_shared, file.path(outdir, paste0(key, ".", source_name, "_shared_asvs.txt")), sep = "\t", row.names = F, quote = F)

# per habitat and treatment (pots only): mean and sd (stats::sd: sd is the sample data frame here)
smt_summary = smt_shared %>% filter(!is_T0) %>% group_by(habitat, treatment) %>%
  summarise(n = n(), across(c(shared_asvs, pct_pellet_asvs, shared_relab), list(mean = mean, sd = stats::sd)), .groups = "drop")
write.table(smt_summary, file.path(outdir, paste0(key, ".", source_name, "_shared_asvs.summary.txt")), sep = "\t", row.names = F, quote = F)
print(as.data.frame(smt_summary[, c("habitat", "treatment", "n", "shared_asvs_mean", "shared_asvs_sd")]))

# plot: same layout as the alpha diversity figure (a box per treatment, baseline as stars), but the baseline samples
# in their own panel (they were rarefied with their habitat above)
smt_panel_levels = c(baseline_label, acfg$habitats)
smt_panel = setNames(ifelse(smt_shared$is_T0, baseline_label, as.character(smt_shared$block)), smt_shared$sample)
asv_metric = sprintf("%% of %s ASVs shared\n(right axis: number of ASVs)", source_noun)
reads_metric = "shared ASVs, % of reads"
smt_long = smt_shared %>%
  tidyr::pivot_longer(c(shared_asvs, shared_relab), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(c(shared_asvs = asv_metric, shared_relab = reads_metric)[metric], levels = c(asv_metric, reads_metric)),
         block = factor(smt_panel[sample], levels = smt_panel_levels),
         group = as.character(treatment))
smt_levels = unique(paste(smt_panel[smt_samples], smt_shared$treatment[match(smt_samples, smt_shared$sample)], sep = "|"))
smt_long$x = factor(paste(smt_long$block, smt_long$group, sep = "|"), levels = smt_levels)
# per-panel y scales (ggh4x): the shared-ASV panels show % of the source's ASVs on the left and the number of shared
# ASVs on the right (values are counts; the count axis is drawn on the right, the % axis is its secondary). The source
# was rarefied to a different depth per habitat (baseline with its habitat), so the conversion differs per panel.
# Shared-ASV panels span their own data range (baseline and its habitat share one: same rarefaction); the % of reads
# row keeps one range (over the panels shown). Versions: all panels, and the config's subset_figure habitats.
pellet_n = c(setNames(smt_pellet_asvs[[baseline_habitat]], baseline_label), smt_pellet_asvs)
# t0_in_bulk: each baseline sample drawn in its habitat panel on top of its treatment's box (same x position, star
# symbol); otherwise in their own panel. metrics: the rows to draw; txt: text size multiplier
smt_plot = function(panels, t0_in_bulk = FALSE, metrics = levels(smt_long$metric), txt = 1) {
  dat = filter(smt_long, block %in% panels, metric %in% metrics) %>% mutate(metric = droplevels(metric))
  if(t0_in_bulk) {
    dat = mutate(dat, block = factor(ifelse(is_T0, baseline_habitat, as.character(block)), levels = smt_panel_levels))
    dat$x = factor(paste(dat$block, dat$group, sep = "|"), levels = levels(smt_long$x))
    dat = droplevels(dat)
  }
  shown = filter(dat, metric == asv_metric)
  asv_lim = tapply(shown$value, as.character(shown$block), range, simplify = FALSE)
  bl = c(baseline_label, baseline_habitat)
  if(all(bl %in% names(asv_lim))) asv_lim[bl] = list(range(unlist(asv_lim[bl])))
  y_scales = c(lapply(names(asv_lim), function(b) rlang::new_formula(
      bquote(metric == .(asv_metric) & block == .(b)),
      bquote(scale_y_continuous(limits = .(asv_lim[[b]]), position = "right", sec.axis = sec_axis(~ . * .(100 / pellet_n[[b]])))))),
    list(rlang::new_formula(bquote(metric != .(asv_metric)), bquote(scale_y_continuous(limits = .(range(dat$value[dat$metric != asv_metric])))))))
  ggplot(dat, aes(x, value)) +
    geom_blank() +
    # all treatments in one color (treatments are named on the x axis): unfilled boxes, points in the same color
    geom_boxplot(data = filter(dat, !is_T0), color = "grey25", fill = NA, outlier.shape = NA, width = 0.65, linewidth = 0.5) +
    geom_point(aes(shape = is_T0), color = "grey25", size = 4.2 * txt, position = position_jitter(width = 0.15, height = 0, seed = 100)) +   # horizontal only
    ggh4x::facet_grid2(metric ~ block, scales = "free", space = "free_x", switch = "y", independent = "y") +
    ggh4x::facetted_pos_scales(y = y_scales) +
    scale_x_discrete(labels = group_axis_labels) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8), labels = c("FALSE" = "pot", "TRUE" = baseline_label), name = NULL) +
    guides(shape = guide_legend(override.aes = list(size = 6 * txt))) +
    labs(title = sprintf("%s ASVs shared with the %s\n(each habitat rarefied with the %s to its smallest library; mean of %d rarefactions)",
                         locus, source_label, source_noun, n_rare),
         x = NULL, y = NULL) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 24 * txt), axis.text.y = element_text(size = 24 * txt),
          panel.grid.major.x = element_blank(), plot.title = element_text(size = 26 * txt),
          strip.text.x = element_text(size = 32 * txt, face = "bold"), strip.text.y = element_text(size = 28 * txt),
          legend.text = element_text(size = 22 * txt), legend.key.size = unit(1.2 * txt, "cm"),
          strip.placement = "outside", strip.background.y = element_blank(), panel.spacing.x = unit(1.1, "cm"))
}
smt_file = paste0(source_name, "_shared_asvs")
save_plot(smt_plot(smt_panel_levels), smt_file, width = 30, height = 15)
# subset of habitats (e.g. soil only), baseline drawn in its habitat panel; and its % of reads row only, larger text
sub_fig = acfg$source_sharing$subset_figure
if(!is.null(sub_fig)) {
  sub_panels = c(baseline_label, sub_fig$habitats)
  save_plot(smt_plot(sub_panels, t0_in_bulk = TRUE), paste0(smt_file, ".", sub_fig$name), width = 19, height = 16)
  save_plot(smt_plot(sub_panels, t0_in_bulk = TRUE, metrics = reads_metric, txt = 1.25),
            paste0(smt_file, ".", sub_fig$name, ".reads"), width = 20, height = 12)
}
}

# ---- 8. distances: a reference treatment vs every other treatment, per habitat (config reference_distances) ----
# every distance between a reference pot and a pot of each other treatment in the same habitat (pots only, no
# baseline; each habitat from the matrix of its ordination subset). Distances: the NMDS ones (with a tree: weighted
# UniFrac, without: Bray-Curtis) and the presence/absence one on the same samples (unweighted UniFrac or binary
# Jaccard; it depends on sequencing depth more than the abundance-weighted distances).
if(!is.null(acfg$reference_distances)) {
rd = acfg$reference_distances
dist_ref = rd$reference
dist_file = paste0("dist_", rd$name, "_vs_treatments")
dist_others = setdiff(treatment_levels, c(dist_ref, source_label))
dist_sets = setNames(list(nmds_dist), dist_method)
dist_sets[[pa_method]] = lapply(nmds_dist, function(d) {
  ps_sub = prune_samples(labels(d), ps)
  calc_distance(prune_taxa(taxa_sums(ps_sub) > 0, ps_sub), pa_method)
})
dist_names = nmds_names
dist_plot_habitats = rd$plot_habitats
all_ref_dist = list()
for(dm_method in names(dist_sets)) {
  dm_name = dist_names[[dm_method]]
  ref_dist = list()
  for(h in names(pair_habitats)) {
    dm = as.matrix(dist_sets[[dm_method]][[pair_habitats[[h]]]])
    pots = rownames(sd)[sd$habitat %in% h & !(sd$is_T0 %in% 1) & rownames(sd) %in% rownames(dm)]
    ref = pots[sd[pots, "treatment"] %in% dist_ref]
    for(trt in dist_others) {
      other = pots[sd[pots, "treatment"] %in% trt]
      if(length(ref) == 0 || length(other) == 0) next
      pr = expand.grid(sample1 = ref, sample2 = other, stringsAsFactors = FALSE)
      ref_dist[[length(ref_dist) + 1]] = data.frame(habitat = h, treatment = trt, comparison = paste(dist_ref, "vs", trt),
        sample1 = pr$sample1, sample2 = pr$sample2, distance = dm[cbind(pr$sample1, pr$sample2)])
    }
  }
  ref_dist = do.call(rbind, ref_dist) %>% mutate(distance_type = dm_name, .before = 1)
  write.table(ref_dist, file.path(outdir, paste0(key, ".", dist_file, ".", dm_method, ".txt")), sep = "\t", row.names = F, quote = F)
  ref_dist_summary = ref_dist %>% group_by(distance_type, habitat, treatment, comparison) %>%
    summarise(n_pairs = n(), mean = mean(distance), sd = stats::sd(distance), median = median(distance),
              min = min(distance), max = max(distance), .groups = "drop") %>%
    mutate(habitat = factor(habitat, levels = names(pair_habitats)),
           treatment = factor(treatment, levels = dist_others)) %>% arrange(habitat, treatment)
  write.table(ref_dist_summary, file.path(outdir, paste0(key, ".", dist_file, ".", dm_method, ".summary.txt")), sep = "\t", row.names = F, quote = F)
  print(as.data.frame(ref_dist_summary), digits = 3)

  all_ref_dist[[dm_method]] = ref_dist

  # plot: one panel per habitat (plot_habitats), a box per comparison (pairs as points), same style as the
  # source-sharing figure
  ref_plot = ref_dist %>% filter(habitat %in% dist_plot_habitats) %>%
    mutate(habitat = factor(habitat, levels = dist_plot_habitats),
           group = factor(paste("vs", treatment_axis_labels(treatment)), levels = paste("vs", treatment_axis_labels(dist_others))))
  p = ggplot(ref_plot, aes(group, distance)) +
    geom_boxplot(color = "grey25", fill = NA, outlier.shape = NA, width = 0.65, linewidth = 0.5) +
    geom_point(color = "grey25", size = 3.2, position = position_jitter(width = 0.15, height = 0, seed = 100)) +
    facet_wrap(~ habitat, nrow = 1, scales = "free_y") +
    labs(title = sprintf("%s %s distances from %s pots, per habitat", locus, dm_name, treatment_name(dist_ref)),
         x = NULL, y = paste(dm_name, "distance")) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), axis.text.y = element_text(size = 18),
          axis.title.y = element_text(size = 20), panel.grid.major.x = element_blank(), plot.title = element_text(size = 18),
          strip.text = element_text(size = 24, face = "bold"), panel.spacing.x = unit(1, "cm"))
  save_plot(p, paste0(dist_file, ".", dm_method), width = 14, height = 10)
}

# combined figure: abundance (top: weighted UniFrac or Bray-Curtis) and presence/absence (bottom: unweighted UniFrac
# or Jaccard) distances, plot_habitats; each panel its own y range
{
  comb_plot = bind_rows(all_ref_dist[nmds_methods]) %>% filter(habitat %in% dist_plot_habitats) %>%
    mutate(habitat = factor(habitat, levels = dist_plot_habitats),
           distance_type = factor(distance_type, levels = dist_names[nmds_methods]),
           group = factor(paste("vs", treatment_axis_labels(treatment)), levels = paste("vs", treatment_axis_labels(dist_others))))
  p = ggplot(comb_plot, aes(group, distance)) +
    geom_boxplot(color = "grey25", fill = NA, outlier.shape = NA, width = 0.65, linewidth = 0.5) +
    geom_point(color = "grey25", size = 3.2, position = position_jitter(width = 0.15, height = 0, seed = 100)) +
    ggh4x::facet_grid2(distance_type ~ habitat, scales = "free_y", independent = "y", switch = "y") +
    labs(title = sprintf("%s %s distances from %s pots, per habitat", locus, if(has_tree) "UniFrac" else "Bray-Curtis and Jaccard",
                         treatment_name(dist_ref)),
         x = NULL, y = NULL) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), axis.text.y = element_text(size = 18),
          panel.grid.major.x = element_blank(), plot.title = element_text(size = 18),
          strip.text.x = element_text(size = 24, face = "bold"), strip.text.y = element_text(size = 20),
          strip.placement = "outside", strip.background.y = element_blank(), panel.spacing = unit(1, "cm"))
  save_plot(p, paste0(dist_file, ".", if(has_tree) "unifrac_weighted_unweighted" else paste(nmds_methods, collapse = "_")),
            width = 14, height = 15)
}

# selected contrasts (reference vs each of rd$selected); table with every distance (all habitats; both distances).
# Figure laid out like the by-habitat NMDS grid: rows = distance, columns = plot_habitats, a box plot per contrast in
# each panel (own y range)
if(length(rd$selected) > 0) {
sel_trts = rd$selected
sel_file = paste0("dist_", rd$name, "_vs_", rd$selected_name)
sel_dist = bind_rows(all_ref_dist) %>% filter(treatment %in% sel_trts)
write.table(sel_dist, file.path(outdir, paste0(key, ".", sel_file, ".txt")), sep = "\t", row.names = F, quote = F)
print(as.data.frame(sel_dist %>% group_by(distance_type, habitat, treatment) %>%
  summarise(n_pairs = n(), mean = mean(distance), median = median(distance), .groups = "drop")), digits = 3)
# y axis label: a double-headed arrow, "more distant" at the top and "less distant" at the bottom, with the axis
# title in a gap in the middle of the arrow (drawn to the left of each panel, aligned to it)
# (lengths in absolute units: the end labels sit at the panel's top and bottom edges, the arrow fills the rest)
arrow_gp = grid::gpar(fontsize = 18)
arrow_text = function(label, ...) grid::textGrob(label, x = 0.5, rot = 90, gp = arrow_gp, ...)
arrow_title = arrow_text(paste("distance from", rd$label), y = 0.5)
arrow_more = arrow_text("more distant", y = unit(1, "npc"), hjust = 1)   # rotated: hjust = 1 ends the text at y
arrow_less = arrow_text("less distant", y = unit(0, "npc"), hjust = 0)
pad = unit(2.5, "mm")
text_len = function(g) grid::grobWidth(grid::editGrob(g, rot = 0))   # length along the text (rotated grobs report their thickness)
arrow_seg = function(y0, y1, ends) grid::segmentsGrob(0.5, y0, 0.5, y1, gp = grid::gpar(lwd = 2.5),
  arrow = grid::arrow(ends = ends, length = unit(0.3, "cm"), type = "open"))
distance_arrow = grid::gTree(children = grid::gList(
  arrow_seg(text_len(arrow_less) + pad, unit(0.5, "npc") - 0.5 * text_len(arrow_title) - pad, "first"),
  arrow_seg(unit(0.5, "npc") + 0.5 * text_len(arrow_title) + pad, unit(1, "npc") - text_len(arrow_more) - pad, "last"),
  arrow_title, arrow_more, arrow_less))
grid_cells = c(list(patchwork::plot_spacer()), lapply(tools::toTitleCase(dist_plot_habitats), header))
for(dm_method in names(all_ref_dist)) {
  grid_cells = c(grid_cells, list(header(tools::toTitleCase(dist_names[[dm_method]]), rot = 90)))
  for(hb in dist_plot_habitats) {
    sp = sel_dist %>% filter(distance_type == dist_names[[dm_method]], habitat == hb) %>%
      mutate(group = factor(treatment_name(treatment, 14), levels = treatment_name(sel_trts, 14)))
    grid_cells[[length(grid_cells) + 1]] = ggplot(sp, aes(group, distance)) +
      geom_boxplot(color = "grey25", fill = NA, outlier.shape = NA, width = 0.6, linewidth = 0.5) +
      geom_point(color = "grey25", size = 4, position = position_jitter(width = 0.12, height = 0, seed = 100)) +
      labs(x = NULL, y = NULL) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18),
            panel.grid.major.x = element_blank(), plot.margin = margin(5.5, 5.5, 5.5, 80)) +
      patchwork::inset_element(distance_arrow, left = -0.3, right = -0.2, bottom = 0, top = 1,
                               align_to = "panel", clip = FALSE)
  }
}
p = patchwork::wrap_plots(grid_cells, ncol = length(dist_plot_habitats) + 1, byrow = TRUE,
                          widths = c(0.1, rep(1, length(dist_plot_habitats))), heights = c(0.1, rep(1, length(all_ref_dist)))) +
  patchwork::plot_annotation(title = sprintf("%s distances from %s pots", locus, treatment_name(dist_ref)))
save_plot(p, sel_file, width = 8 * length(dist_plot_habitats), height = 9 * length(all_ref_dist) + 1.5)
}
}

# ---- 9. source colonization (turnover) or pruning/boosting of native ASVs (nestedness), betapart (config betapart) ----
# Sorensen dissimilarity between two pots split into turnover (beta.sim: ASVs replaced by other ASVs) and nestedness
# (beta.sne: the poorer pot holds a subset of the richer pot's ASVs); betapart::beta.pair on presence/absence. The
# abundance partition (beta.pair.abund) is not used: its gradient component is the difference in total counts, which
# is 0 for rarefied samples (and only sequencing depth without rarefaction).
# contrasts: pairs (treatment with the source, the same treatment without it; e.g. Conv_BC-SMT vs Conv_BC). If the
# source colonizes the recipient soil, with- vs without-source pots should show more turnover than replicate pots
# within a treatment; if it only prunes or boosts native ASVs, the extra dissimilarity is nestedness (with-source pots
# poorer: pruning; richer: gain of ASVs nested within the richer pot's set). Pots only (no baseline), per habitat;
# counts rarefied to the habitat block's depth (above), n_rare times (seeds 100, 101, ...), indices averaged.
if(!is.null(acfg$betapart)) {
bp_cfg = acfg$betapart
smt_contrasts = bp_cfg$contrasts
if(is.matrix(smt_contrasts)) smt_contrasts = lapply(seq_len(nrow(smt_contrasts)), function(i) smt_contrasts[i, ])
bp_label = bp_cfg$label              # e.g. "SMT vs no SMT"
bp_pots_label = bp_cfg$pots_label    # e.g. "SMT vs no-SMT pots"
bp_file = paste0("betapart_", bp_cfg$name)
bp_trts = unique(unlist(smt_contrasts))
bp_metrics = c("beta.sor", "beta.sim", "beta.sne")
bp_pairs = list()
for(h in intersect(names(pair_habitats), smt_blocks)) {
  sn = rownames(sd)[sd$habitat %in% h & !(sd$is_T0 %in% 1) & sd$treatment %in% bp_trts]
  res = lapply(seq_len(n_rare), function(i) {
    set.seed(99 + i)
    r = vegan::rrarefy(counts_m[sn, ], smt_depth[h])
    r = r[, colSums(r) > 0, drop = F]
    c(lapply(betapart::beta.pair(1 * (r > 0), index.family = "sorensen"), as.matrix),
      list(richness = rowSums(r > 0)))
  })
  avg = setNames(lapply(c(bp_metrics, "richness"), function(m) Reduce(`+`, lapply(res, `[[`, m)) / n_rare),
                 c(bp_metrics, "richness"))
  pr = t(combn(sn, 2))
  t1 = as.character(sd[pr[, 1], "treatment"]); t2 = as.character(sd[pr[, 2], "treatment"])
  for(k in seq_len(nrow(pr))) {
    s1 = pr[k, 1]; s2 = pr[k, 2]
    if(t1[k] == t2[k]) {
      type = "replicates"; trt = t1[k]; label = paste("within", t1[k])
    } else {
      hit = Filter(function(x) setequal(x, c(t1[k], t2[k])), smt_contrasts)
      if(length(hit) == 0) next   # e.g. Conv_BC vs Conv_BC-MPSS: not a contrast
      type = bp_label; trt = hit[[1]][1]; label = paste(hit[[1]], collapse = " vs ")
      if(t1[k] != hit[[1]][1]) { s1 = pr[k, 2]; s2 = pr[k, 1] }   # sample1 = the with-source pot
    }
    bp_pairs[[length(bp_pairs) + 1]] = data.frame(habitat = h, depth = smt_depth[h], type = type, treatment = trt,
      comparison = label, sample1 = s1, sample2 = s2, as.list(sapply(bp_metrics, function(m) avg[[m]][s1, s2])),
      richness1 = avg$richness[s1], richness2 = avg$richness[s2], row.names = NULL)
  }
  cat(sprintf("betapart, %s: %d pots rarefied %d x to %d reads\n", h, length(sn), n_rare, smt_depth[h]))
}
bp_pairs = do.call(rbind, bp_pairs) %>% mutate(nestedness_share = beta.sne / beta.sor)
write.table(bp_pairs, file.path(outdir, paste0(key, ".", bp_file, ".txt")), sep = "\t", row.names = F, quote = F)
bp_summary = bp_pairs %>% group_by(habitat, type, treatment, comparison) %>%
  summarise(n_pairs = n(), across(c(all_of(bp_metrics), nestedness_share), mean),
            richness_diff = mean(richness1 - richness2),                 # with- minus without-source pot (replicates: arbitrary sign)
            pct_sample1_richer = 100 * mean(richness1 > richness2), .groups = "drop") %>%
  mutate(habitat = factor(habitat, levels = names(pair_habitats)),
         treatment = factor(treatment, levels = treatment_levels)) %>% arrange(habitat, treatment, type)
bp_summary$richness_diff[bp_summary$type == "replicates"] = NA
bp_summary$pct_sample1_richer[bp_summary$type == "replicates"] = NA
write.table(bp_summary, file.path(outdir, paste0(key, ".", bp_file, ".summary.txt")), sep = "\t", row.names = F, quote = F)
print(as.data.frame(bp_summary), digits = 3)

# plot: rows = Sorensen components (turnover, nestedness), columns
# = habitats; per contrast, the two treatments' replicate pairs next to the with- vs without-source pairs
bp_components = c(beta.sim = "turnover (Simpson)", beta.sne = "nestedness")
bp_groups = unlist(lapply(smt_contrasts, function(x)
  c(paste("within", treatment_axis_labels(x[2])), paste("within", treatment_axis_labels(x[1])),
    paste0(bp_label, "\n(", treatment_name(x[2]), ")"))))
without_source = setNames(vapply(smt_contrasts, `[`, "", 2), vapply(smt_contrasts, `[`, "", 1))
bp_plot = bp_pairs %>%
  mutate(group = ifelse(type == "replicates", paste("within", treatment_axis_labels(treatment)),
                        paste0(bp_label, "\n(", treatment_name(without_source[treatment]), ")"))) %>%
  tidyr::pivot_longer(all_of(names(bp_components)), names_to = "component", values_to = "dissimilarity") %>%
  mutate(habitat = factor(habitat, levels = names(pair_habitats)), group = factor(group, levels = bp_groups),
         component = factor(bp_components[component], levels = bp_components),
         contrast = ifelse(type == "replicates", "replicate pots", bp_pots_label))
p = ggplot(bp_plot, aes(group, dissimilarity)) +
  geom_boxplot(aes(fill = contrast), color = "grey25", outlier.shape = NA, width = 0.65, linewidth = 0.5) +
  geom_point(color = "grey25", size = 2.4, position = position_jitter(width = 0.15, height = 0, seed = 100)) +
  scale_fill_manual(values = setNames(c(NA, "grey80"), c("replicate pots", bp_pots_label)), name = NULL) +
  facet_grid(component ~ habitat, scales = "free_y", switch = "y") +
  labs(title = sprintf("%s beta diversity partition (betapart): %s and replicate pots (rarefied per habitat; mean of %d rarefactions)",
                       locus, bp_pots_label, n_rare), x = NULL, y = NULL) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15),
        panel.grid.major.x = element_blank(), plot.title = element_text(size = 16),
        strip.text.x = element_text(size = 22, face = "bold"), strip.text.y = element_text(size = 17),
        strip.placement = "outside", strip.background.y = element_blank(), legend.position = "top",
        legend.text = element_text(size = 18), legend.key.size = unit(1, "cm"), panel.spacing = unit(0.6, "cm"),
        plot.margin = margin(10, 10, 10, 70))
save_plot(p, bp_file, width = 22, height = 12)
}
