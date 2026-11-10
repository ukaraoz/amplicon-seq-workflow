# Join the 2018 soil physical/chemical measurements to the NEON ER amplicon metadata by SamplingCode, the 3-digit
# site number (058 for 058-ER18, 058-ER22, 058-ER18_A ...): all samples of a site, both years, get the same row.
# usage: Rscript merge_soil.R <metadata.tsv> <soil.csv> <out.tsv>
a = commandArgs(TRUE)
md = read.delim(a[1], check.names = FALSE, colClasses = "character")
soil = read.csv(a[2], check.names = FALSE, colClasses = c(SampleSiteCode = "character", IGSN = "character"))
code = function(x) sub("^([0-9]{3}).*", "\\1", x)
md$SamplingCode = code(md$sample)
soil$SamplingCode = code(soil$SampleSiteCode)
stopifnot(all(grepl("^[0-9]{3}$", md$SamplingCode)), all(grepl("^[0-9]{3}$", soil$SamplingCode)), !anyDuplicated(soil$SamplingCode))
# soil SampleSiteCode is the join source only (the metadata has its own SampleSiteCode); other columns unchanged
soil$SampleSiteCode = NULL
out = merge(md, soil, by = "SamplingCode", all.x = TRUE, sort = FALSE)
out = out[match(md$sample, out$sample), c(names(md), setdiff(names(soil), "SamplingCode"))]
# keep the metadata's column order, SamplingCode after SampleSiteCode
cols = names(out); cols = c(setdiff(cols[seq_len(match("SampleSiteCode", cols))], "SamplingCode"), "SamplingCode",
                            setdiff(cols[-seq_len(match("SampleSiteCode", cols))], "SamplingCode"))
out = out[, cols]
write.table(out, a[3], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(nrow(out), "samples x", ncol(out), "columns ->", a[3], "\n")
cat("samples with a soil row:", sum(!is.na(out$IGSN) & out$IGSN != ""), "of", nrow(out), "\n")
nomatch = out$sample[is.na(out$IGSN) | out$IGSN == ""]
cat("without:", if(length(nomatch)) paste(nomatch, collapse = " ") else "none", "\n")
cat("soil rows used by no sample:", sum(!soil$SamplingCode %in% md$SamplingCode), "\n")
cat("example, site 058:\n"); print(out[out$SamplingCode == "058", c("sample", "year", "SampleSiteCode", "SamplingCode", "Sampling Area", "Soil pH", "Soil TOC")], row.names = FALSE)
