# Fungal guilds / lifestyles for ITS ASVs from their UNITE taxonomy: FungalTraits and FUNGuild
#
# version 1.2
#
# usage: Rscript fungal-guilds.R --config=<key>.config.json     (workflow.py `guilds` step, after `phyloseq`)
#   reads   <results>/<key>.phyloseq.rds (all ASVs) and <key>.phyloseq.filtered.rds
#           params.guilds.fungaltraits (FungalTraits_1.2.tsv) and params.guilds.funguild (FUNGuild.tsv),
#           made by get-guild-dbs.py
#   writes  <key>.guilds.txt              one row per ASV (all ASVs): taxonomy, reads, samples, FungalTraits columns
#                                         (FT_*), FUNGuild columns (FG_*), how each database was matched
#           <key>.guild_summary.txt       per database / field / value: ASVs, reads, % of fungal reads, samples, over
#                                         the filtered phyloseq (Fungi kept by the taxa filter, passing samples),
#                                         with rows for fungal reads that could not be assigned and why
#           <key>.guild_abundance_FungalTraits.txt   wide: one row per sample (filtered phyloseq), fungal_reads,
#                                         then the relative abundance (share of the sample's fungal reads) of each
#                                         FungalTraits primary lifestyle; columns sum to 1 per sample
#           <key>.guild_abundance.txt                FUNGuild, wide: one row per sample, fungal_reads, then the relative
#                                         abundance of each FUNGuild trophic mode ("trophicMode: ...") and each guild
#                                         ("guild: ..."); each of the two sets sums to 1 per sample
#           ("unassigned": fungal reads without a match in that database)
#           <key>.phyloseq.filtered.guilds.rds   the filtered phyloseq with guild columns added to the tax table
#
# matching (names without UNITE's rank prefixes; placeholder names such as Ascomycota_gen_Incertae_sedis count as
# unassigned):
#   FungalTraits: genus-level database -> exact genus name
#   FUNGuild: the most specific rank that has a record: species ("Genus epithet"), genus, family, order, class,
#   phylum, as FUNGuild's own script does. Records keep their confidence ranking (Highly Probable / Probable /
#   Possible); FUNGuild's authors suggest relying on Probable or better.
# A guild from a genus-level match is a hypothesis for the genus; species in a genus can differ.
suppressPackageStartupMessages({library(phyloseq); library(jsonlite)})

args = commandArgs(TRUE)
cfgfile = sub("^--config=", "", args[grepl("^--config=", args)])
if(length(cfgfile) != 1) stop("usage: Rscript fungal-guilds.R --config=<key>.config.json")
cfg = fromJSON(cfgfile, simplifyVector = TRUE)
key = cfg$key; outdir = cfg$paths$results; g = cfg$params$guilds
if(is.null(g$fungaltraits) || is.null(g$funguild)) stop("params.guilds needs fungaltraits and funguild (get-guild-dbs.py tables)")
cat("fungal-guilds.R version 1.2 :", key, "\n")

ps = readRDS(file.path(outdir, paste0(key, ".phyloseq.rds")))
psf = readRDS(file.path(outdir, paste0(key, ".phyloseq.filtered.rds")))

# ---- taxonomy without prefixes; placeholders -> NA ----
tt0 = as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
tt0[] = lapply(tt0, function(x) sub("^[a-z]__", "", x))             # names as assigned (for the output table)
tt = tt0                                                              # for matching: placeholders -> NA
tt[] = lapply(tt, function(x) { x[grepl("Incertae_sedis|_gen_|_fam_|_ord_|_cls_|_phy_|^unidentified", x)] = NA; x })
rank = function(r) if(r %in% names(tt)) tt[[r]] else rep(NA_character_, nrow(tt))
genus = rank("Genus")
species = ifelse(is.na(genus) | is.na(rank("Species")), NA, paste(genus, gsub("_", " ", rank("Species"))))
is_fungus = !is.na(rank("Kingdom")) & rank("Kingdom") == "Fungi"

# ---- FungalTraits (genus) ----
ft = read.delim(g$fungaltraits, quote = "", stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
ft = ft[!duplicated(ft$GENUS), ]
ftcols = setdiff(names(ft), c("jrk", "Phylum", "Class", "Order", "Family", "GENUS", "COMMENT on genus"))
fti = match(genus, ft$GENUS)
FT = setNames(as.data.frame(lapply(ftcols, function(cn) ft[[cn]][fti]), stringsAsFactors = FALSE), paste0("FT_", ftcols))
FT_match = ifelse(!is_fungus, "not Fungi", ifelse(is.na(genus), "no genus assigned",
                  ifelse(is.na(fti), "genus not in FungalTraits", "genus")))

# ---- FUNGuild (most specific rank with a record) ----
fg = read.delim(g$funguild, quote = "", stringsAsFactors = FALSE, na.strings = c("", "NULL"))
fg$guild = gsub("\\|", "", fg$guild)          # FUNGuild marks uncertain guild components with |...|
fg = fg[!duplicated(fg$taxon), ]
fg_ranks = list(species = species, genus = genus, family = rank("Family"), order = rank("Order"),
                class = rank("Class"), phylum = rank("Phylum"))
fgi = rep(NA_integer_, nrow(tt)); fglevel = rep(NA_character_, nrow(tt))
for(r in names(fg_ranks)) {
  m = match(fg_ranks[[r]], fg$taxon)
  new = is.na(fgi) & !is.na(m) & is_fungus
  fgi[new] = m[new]; fglevel[new] = r
}
FG = data.frame(FG_taxon = fg$taxon[fgi], FG_match = ifelse(!is_fungus, "not Fungi", ifelse(is.na(fgi), "no FUNGuild record", fglevel)),
                FG_trophicMode = fg$trophicMode[fgi], FG_guild = fg$guild[fgi], FG_confidence = fg$confidenceRanking[fgi],
                FG_growthForm = fg$growthForm[fgi], FG_trait = fg$trait[fgi], FG_citation = fg$citationSource[fgi],
                stringsAsFactors = FALSE)

# ---- per-ASV table (all ASVs) ----
otu_all = as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)) otu_all = t(otu_all)
asv = data.frame(asv_id = taxa_names(ps), asv_seq = as.character(refseq(ps)), tt0, reads = colSums(otu_all),
                 samples = colSums(otu_all > 0), in_filtered = taxa_names(ps) %in% taxa_names(psf),
                 FT_match = FT_match, FT, FG, stringsAsFactors = FALSE, check.names = FALSE)
write.table(asv, file.path(outdir, paste0(key, ".guilds.txt")), sep = "\t", quote = F, row.names = F, na = "")

# ---- summaries over the filtered phyloseq (fungal reads in passing samples) ----
otu = as(otu_table(psf), "matrix")
if(taxa_are_rows(psf)) otu = t(otu)
a = asv[match(colnames(otu), asv$asv_id), ]
total = sum(otu)
summ = function(db, field, value) {
  value[is.na(value)] = "unassigned"
  do.call(rbind, lapply(split(seq_along(value), value), function(i) data.frame(database = db, field = field,
    value = value[i[1]], ASVs = length(i), reads = sum(otu[, i]), pct_fungal_reads = round(100 * sum(otu[, i]) / total, 3),
    samples = sum(rowSums(otu[, i, drop = F]) > 0), stringsAsFactors = FALSE)))
}
lst = list(summ("FungalTraits", "match", a$FT_match),
           summ("FungalTraits", "primary_lifestyle", a$FT_primary_lifestyle),
           summ("FungalTraits", "Secondary_lifestyle", a$FT_Secondary_lifestyle),
           summ("FungalTraits", "Ectomycorrhiza_exploration_type", ifelse(a$FT_primary_lifestyle %in% "ectomycorrhizal", a$FT_Ectomycorrhiza_exploration_type, "(not ectomycorrhizal)")),
           summ("FungalTraits", "Decay_type", a$FT_Decay_type),
           summ("FungalTraits", "Plant_pathogenic_capacity", a$FT_Plant_pathogenic_capacity),
           summ("FUNGuild", "match", a$FG_match),
           summ("FUNGuild", "confidence", a$FG_confidence),
           summ("FUNGuild", "trophicMode", a$FG_trophicMode),
           summ("FUNGuild", "trophicMode (Probable or better)", ifelse(a$FG_confidence %in% c("Probable", "Highly Probable"), a$FG_trophicMode, NA)),
           summ("FUNGuild", "guild", a$FG_guild))
S = do.call(rbind, lst)
S = S[order(S$database, match(S$field, unique(S$field)), -S$reads), ]
write.table(S, file.path(outdir, paste0(key, ".guild_summary.txt")), sep = "\t", quote = F, row.names = F)

# ---- per-sample relative abundance, wide: one file per database, one row per sample ----
# columns ordered by total abundance (largest first), "unassigned" last; prefix (e.g. "guild: ") when a file has
# more than one field
wide = function(value, prefix = "") {
  value[is.na(value)] = "unassigned"
  grp = split(seq_along(value), value)
  m = matrix(sapply(grp, function(i) rowSums(otu[, i, drop = F])), nrow = nrow(otu), dimnames = list(rownames(otu), names(grp)))
  o = order(colnames(m) == "unassigned", -colSums(m))
  m = round(m[, o, drop = F] / rowSums(otu), 6)
  colnames(m) = paste0(prefix, colnames(m))
  m
}
write_wide = function(m, db, suffix) {
  w = data.frame(sample = rownames(otu), fungal_reads = as.integer(rowSums(otu)), m, check.names = FALSE)
  f = file.path(outdir, paste0(key, ".guild_abundance", suffix, ".txt"))
  write.table(w, f, sep = "\t", quote = F, row.names = F)
  cat(sprintf("%s: %d samples x %d categories -> %s\n", db, nrow(w), ncol(m), basename(f)))
}
write_wide(wide(a$FT_primary_lifestyle), "FungalTraits", "_FungalTraits")
write_wide(cbind(wide(a$FG_trophicMode, "trophicMode: "), wide(a$FG_guild, "guild: ")), "FUNGuild", "")
unlink(file.path(outdir, paste0(key, ".guild_abundance_FUNGuild.txt")))   # name in version 1.1

# ---- filtered phyloseq with guild columns ----
add = c("FT_match", "FT_primary_lifestyle", "FT_Secondary_lifestyle", "FT_Ectomycorrhiza_exploration_type",
        "FT_Decay_type", "FT_Plant_pathogenic_capacity", "FG_match", "FG_trophicMode", "FG_guild", "FG_confidence")
tx = cbind(as(tax_table(psf), "matrix"), as.matrix(a[, add]))
tax_table(psf) = tax_table(tx)
saveRDS(psf, file.path(outdir, paste0(key, ".phyloseq.filtered.guilds.rds")))

# ---- log ----
pct = function(sel) round(100 * sum(otu[, sel]) / total, 1)
cat(sprintf("fungal reads in the filtered phyloseq: %d (%d ASVs, %d samples)\n", total, ncol(otu), nrow(otu)))
cat(sprintf("FungalTraits: genus match %.1f%% of fungal reads; no genus %.1f%%; genus not in FungalTraits %.1f%%\n",
            pct(a$FT_match == "genus"), pct(a$FT_match == "no genus assigned"), pct(a$FT_match == "genus not in FungalTraits")))
cat(sprintf("FUNGuild: record found %.1f%% of fungal reads (species %.1f%%, genus %.1f%%, higher %.1f%%); Probable or better %.1f%%\n",
            pct(!is.na(a$FG_trophicMode) | a$FG_match %in% names(fg_ranks)), pct(a$FG_match == "species"), pct(a$FG_match == "genus"),
            pct(a$FG_match %in% c("family", "order", "class", "phylum")), pct(a$FG_confidence %in% c("Probable", "Highly Probable"))))
top = S[S$database == "FungalTraits" & S$field == "primary_lifestyle", ][1:min(8, sum(S$field == "primary_lifestyle")), ]
cat("FungalTraits primary lifestyle (% of fungal reads):\n"); print(top[, c("value", "ASVs", "pct_fungal_reads", "samples")], row.names = F)
cat("done ->", outdir, "\n")
