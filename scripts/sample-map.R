# Map of all sampling locations of a dataset (all loci together) on a base map
#
# version 1.0
#
# usage: Rscript sample-map.R --out=<dir> [--name=<prefix>] [--basemap=Esri.WorldImagery] [--zoom=14]
#                             [--cache=<dir>] [--pdf] <results dir> [<results dir> ...]
#   e.g. Rscript scripts/sample-map.R --out=config/neon_er_ER18_ER22/metadata --name=neon_er_ER18_ER22 \
#          results/neon_er_ER18_ER22.16S results/neon_er_ER18_ER22.ITS.R1
# input : <results dir>/<key>.filtered.sample_data.txt of each locus (filtered phyloseq objects; failed samples
#         are not in them); a sample in several loci counts once
# output: <out>/<name>.sample_map.png (.pdf with --pdf): one point per sampling location (Latitude / Longitude;
#           the 2018 and 2022 samples of a location share it), fill = Vegetation.Type, shape = years sampled
#           there, white outline; one label per site (Site Name) with its number of samples
#         <out>/<name>.sample_locations.txt: one row per sample: site, vegetation type, year, coordinates, and the
#           loci it is in
# Base map: map tiles from --basemap (a maptiles provider; default Esri.WorldImagery) at --zoom (default 14,
# ~10 m per pixel), downloaded into --cache (default: a temporary folder, deleted after the run); none = no base
# map. Web Mercator, axes in degrees. Samples without coordinates are listed in the subtitle.
suppressPackageStartupMessages({library(ggplot2)})
script_version = "1.0"

args = commandArgs(TRUE)
opt = function(name, default = NULL) { a = args[grepl(paste0("^--", name, "="), args)]
  if(length(a)) sub(paste0("^--", name, "="), "", a[1]) else default }
dirs = args[!grepl("^--", args)]
outdir = opt("out"); if(is.null(outdir) || !length(dirs)) stop("usage: Rscript sample-map.R --out=<dir> [--name=] <results dir> ...")
name = opt("name", basename(outdir))
provider = opt("basemap", "Esri.WorldImagery")
zoom = as.integer(opt("zoom", "14"))
cachedir = opt("cache", file.path(tempdir(), "tiles"))
write_pdf = "--pdf" %in% args
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
col_site = "Site Name"; col_veg = "Vegetation.Type"; col_year = "year"

# ---- samples of all loci ----
sd = do.call(rbind, lapply(dirs, function(d) {
  f = Sys.glob(file.path(d, "*.filtered.sample_data.txt"))
  if(length(f) != 1) stop("no (or several) *.filtered.sample_data.txt in ", d)
  x = read.delim(f, check.names = FALSE, colClasses = "character")
  locus = sub("^[^.]+\\.([^.]+).*$", "\\1", basename(dirname(f[1])))   # neon_er_ER18_ER22.ITS.R1 -> ITS
  data.frame(sample = x[[1]], site = x[[col_site]], veg = x[[col_veg]], year = x[[col_year]],
             lat = as.numeric(x$Latitude), lon = as.numeric(x$Longitude), locus = locus, stringsAsFactors = FALSE)
}))
loci = unique(sd$locus)
smp = do.call(rbind, lapply(split(sd, sd$sample), function(z) cbind(z[1, c("sample", "site", "veg", "year", "lat", "lon")],
  loci = paste(intersect(loci, z$locus), collapse = "+"))))
smp = smp[order(smp$site, smp$sample), ]
write.table(smp, file.path(outdir, paste0(name, ".sample_locations.txt")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
nocoord = smp$sample[is.na(smp$lat) | is.na(smp$lon)]
m = smp[!(smp$sample %in% nocoord), ]
cat("sample-map.R", script_version, ":", nrow(smp), "samples (", paste(names(table(smp$loci)), table(smp$loci), collapse = ", "),
    ");", length(nocoord), "without coordinates\n")

# ---- locations and sites, in Web Mercator ----
merc = function(lon, lat) data.frame(x = 6378137 * lon * pi / 180, y = 6378137 * log(tan(pi / 4 + lat * pi / 360)))
m = cbind(m, merc(m$lon, m$lat))
years = sort(unique(m$year))
loc = do.call(rbind, lapply(split(m, paste(m$lat, m$lon)), function(z) data.frame(x = z$x[1], y = z$y[1], veg = z$veg[1],
  sampled = if(length(unique(z$year)) > 1) paste(sort(unique(z$year)), collapse = " + ") else paste(z$year[1], "only"))))
loc$sampled = factor(loc$sampled, levels = intersect(c(paste(years, "only"), paste(years, collapse = " + ")), unique(loc$sampled)))
sites = do.call(rbind, lapply(split(m, m$site), function(z) data.frame(x = mean(range(z$x)), y = mean(range(z$y)),
  text = paste0(z$site[1], " (", nrow(z), ")"))))
pad = function(r, f) r + c(-1, 1) * diff(r) * f
xr = pad(range(m$x), 0.08); yr = pad(range(m$y), 0.08)

# ---- base map ----
bm = NULL
if(!identical(provider, "none")) {
  if(!requireNamespace("maptiles", quietly = TRUE) || !requireNamespace("terra", quietly = TRUE)) {
    cat("no base map: R packages maptiles and terra are needed\n")
  } else bm = tryCatch({
    box = terra::vect(cbind(rep(xr, each = 2), c(yr, rev(yr))), type = "polygons", crs = "EPSG:3857")
    tl = maptiles::get_tiles(box, provider = provider, zoom = zoom, crop = TRUE, cachedir = cachedir)
    df = terra::as.data.frame(tl, xy = TRUE)
    df$fill = grDevices::rgb(df[[3]], df[[4]], df[[5]], maxColorValue = 255)
    list(df = df[, c("x", "y", "fill")], credit = maptiles::get_credit(provider))
  }, error = function(e) { cat("no base map (", conditionMessage(e), ")\n"); NULL })
  if(is.null(opt("cache"))) unlink(cachedir, recursive = TRUE)   # temporary tile folder
}

# ---- plot ----
veg_lv = sort(unique(loc$veg))
veg_cols = setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")[seq_along(veg_lv)], veg_lv)  # as in the composition plots
shapes = setNames(c(21, 24, 22, 23)[seq_along(levels(loc$sampled))], levels(loc$sampled))
lonb = pretty(range(m$lon), 5); latb = pretty(range(m$lat), 5)
g = ggplot() +
  (if(!is.null(bm)) geom_raster(data = bm$df, aes(x, y, fill = I(fill)), show.legend = FALSE) else NULL) +
  geom_point(data = loc[order(loc$sampled), ], aes(x, y, fill = veg, shape = sampled), colour = "white", size = 2.6, stroke = 0.5) +
  scale_fill_manual(values = veg_cols, name = "Vegetation type", guide = guide_legend(override.aes = list(shape = 21, size = 3.5))) +
  scale_shape_manual(values = shapes, name = "Sampled", guide = guide_legend(override.aes = list(fill = "grey50", size = 3.5))) +
  ggrepel::geom_label_repel(data = sites, aes(x, y, label = text), size = 3.2, label.size = 0, fill = "white", alpha = 0.9,
                            label.padding = unit(2.5, "pt"), min.segment.length = 0, segment.colour = "white",
                            segment.size = 0.5, box.padding = 1, seed = 1, max.overlaps = Inf) +
  scale_x_continuous(breaks = merc(lonb, 0)$x, labels = paste0(abs(lonb), "°W"), expand = c(0, 0)) +
  scale_y_continuous(breaks = merc(0, latb)$y, labels = paste0(latb, "°N"), expand = c(0, 0)) +
  coord_fixed(xlim = xr, ylim = yr) +
  labs(x = NULL, y = NULL,
       title = paste0(name, ": ", nrow(m), " samples at ", nrow(loc), " locations, ", nrow(sites), " sites"),
       subtitle = paste(strwrap(paste0("samples in ", paste(loci, collapse = " or "), " (filtered objects; ",
                                       paste(names(table(smp$loci)), table(smp$loci), collapse = ", "), "); labels: site (samples)",
                                       if(length(nocoord)) paste0("; no coordinates (not shown): ", paste(sort(nocoord), collapse = ", ")) else ""),
                                130), collapse = "\n"),
       caption = if(!is.null(bm)) paste(strwrap(paste0("Base map: ", provider, " - ", bm$credit), 160), collapse = "\n") else NULL) +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_text(size = 8), plot.caption = element_text(size = 6.5, colour = "grey30"),
        panel.grid = element_blank(), legend.key = element_rect(fill = NA))
asp = diff(yr) / diff(xr)
f = file.path(outdir, paste0(name, ".sample_map"))
h = min(14, 7.2 * asp + 1.3)
ggsave(paste0(f, ".png"), g, width = 9, height = h, dpi = 300)
if(write_pdf) ggsave(paste0(f, ".pdf"), g, width = 9, height = h)
cat("sample map:", nrow(m), "samples,", nrow(loc), "locations,", nrow(sites), "sites ->", paste0(f, if(write_pdf) ".png / .pdf" else ".png"), "\n")
