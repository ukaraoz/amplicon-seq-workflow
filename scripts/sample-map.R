# Map of all sampling locations of a dataset (all loci together) on a base map
#
# version 1.9
#
# usage: Rscript sample-map.R --out=<dir> [--name=<prefix>] [--site=<Site Name>] [--basemap=Esri.WorldImagery]
#                             [--zoom=<level>] [--cache=<dir>] [--pdf] <results dir> [<results dir> ...]
#   e.g. Rscript scripts/sample-map.R --out=config/neon_er_ER18_ER22/metadata --name=neon_er_ER18_ER22 \
#          results/neon_er_ER18_ER22.16S results/neon_er_ER18_ER22.ITS.R1
# input : <results dir>/<key>.filtered.sample_data.txt of each locus (filtered phyloseq objects; failed samples
#         are not in them); a sample in several loci counts once
# output: <out>/<name>.sample_map.png (.pdf with --pdf): one point per sampling location (Latitude / Longitude;
#           the 2018 and 2022 samples of a location share it), fill = Vegetation.Type, shape = years sampled
#           there, white outline; one label per site (Site Name) with its number of samples
#         <out>/<name>.sample_locations.txt: one row per sample: site, vegetation type, year, coordinates, and the
#           loci it is in
#   --site=<Site Name>: only that site's samples -> <out>/<name>.<site>.sample_map.png, each location labelled
#           with its sampling code (the NNN of NNN-ER18), no locations table
#   every map has a scale bar (bottom left on site maps, bottom right on the map of all sites)
#   --slope: terrain slope over the base map, from the USGS 3DEP elevation service (1-10 m DEM, slope computed in
#           UTM), in 10-degree classes, semi-transparent, with contour lines (--contour=<m>, default 20; 0 = none);
#           writes <out>/<name>[.<site>].slope_map.png instead of .sample_map.png
#   --colour=family: points coloured by the plant family of the location (column Family) within its vegetation type:
#           each vegetation type keeps its hue (Meadow greens, Shrubs oranges, Trees purples, in alphabetical order of
#           the types) and its most common families (counted over all samples given, so the colours are the same on
#           every map) get shades of it: as many as the hue has shades (4, 3, 2), the last shade for "other <type>"
#           when there are more families; no family (litter) grey. Shape = vegetation type (circle, triangle, square),
#           small star inside the points of locations sampled in both years.
#           --trees=species: Trees by species instead of family (pine, spruce, fir: conifer, darker shades; aspen:
#           broadleaf, lightest), e.g. "Trees: pine (conifer)". Default (--colour=vegetation): fill =
#           vegetation type, shape = years sampled
# Base map: map tiles from --basemap (a maptiles provider; default Esri.WorldImagery) at --zoom (default: from the
# extent, ~2700 pixels across, at most 18; all sites: 14, ~10 m per pixel; one site: 18), downloaded into --cache (default: a temporary folder, deleted after the run); none = no base
# map. Web Mercator, axes in degrees. Title only (no subtitle or caption): samples without coordinates are the rows
# without Latitude / Longitude in sample_locations.txt; the base-map credit is printed (cite it with the figure).
suppressPackageStartupMessages({library(ggplot2)})
script_version = "1.9"

args = commandArgs(TRUE)
opt = function(name, default = NULL) { a = args[grepl(paste0("^--", name, "="), args)]
  if(length(a)) sub(paste0("^--", name, "="), "", a[1]) else default }
dirs = args[!grepl("^--", args)]
outdir = opt("out"); if(is.null(outdir) || !length(dirs)) stop("usage: Rscript sample-map.R --out=<dir> [--name=] <results dir> ...")
name = opt("name", basename(outdir))
provider = opt("basemap", "Esri.WorldImagery")
zoom = opt("zoom"); site = opt("site")
cachedir = opt("cache", file.path(tempdir(), "tiles"))
write_pdf = "--pdf" %in% args
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
col_site = "Site Name"; col_veg = "Vegetation.Type"; col_year = "year"; col_family = "Family"
colour_by = opt("colour", "vegetation")

# ---- samples of all loci ----
sd = do.call(rbind, lapply(dirs, function(d) {
  f = Sys.glob(file.path(d, "*.filtered.sample_data.txt"))
  if(length(f) != 1) stop("no (or several) *.filtered.sample_data.txt in ", d)
  x = read.delim(f, check.names = FALSE, colClasses = "character")
  locus = sub("^[^.]+\\.([^.]+).*$", "\\1", basename(dirname(f[1])))   # neon_er_ER18_ER22.ITS.R1 -> ITS
  data.frame(sample = x[[1]], site = x[[col_site]], veg = x[[col_veg]], year = x[[col_year]],
             family = if(col_family %in% names(x)) x[[col_family]] else NA,
             species = if("PlantGenusSpecies" %in% names(x)) x[["PlantGenusSpecies"]] else NA,
             lat = as.numeric(x$Latitude), lon = as.numeric(x$Longitude), locus = locus, stringsAsFactors = FALSE)
}))
loci = unique(sd$locus)
smp = do.call(rbind, lapply(split(sd, sd$sample), function(z) cbind(z[1, c("sample", "site", "veg", "year", "family", "species", "lat", "lon")],
  loci = paste(intersect(loci, z$locus), collapse = "+"))))
smp = smp[order(smp$site, smp$sample), ]
# plant-family classes (--colour=family), from all samples so that every map uses the same colours
ramps = list(c("#00441b", "#2b9348", "#80c95b", "#d9f0a3"), c("#7f2704", "#d94801", "#fdae6b"), c("#3f007d", "#b4a7e0"))
fam_none = "#8c8c8c"
smp$family[is.na(smp$family) | smp$family == ""] = NA
fam_class = rep(NA_character_, nrow(smp)); fam_cols = c()
for(i in seq_along(vt <- sort(unique(smp$veg)))) {
  sel = smp$veg == vt[i]; r = ramps[[(i - 1) %% length(ramps) + 1]]
  cnt = sort(table(smp$family[sel & !is.na(smp$family)]), decreasing = TRUE)
  named = if(length(cnt) <= length(r)) names(cnt) else names(cnt)[seq_len(length(r) - 1)]
  lab = ifelse(smp$family[sel] %in% named, paste0(vt[i], ": ", smp$family[sel]), paste0(vt[i], ": other"))
  lab[is.na(smp$family[sel])] = paste0(vt[i], ": no plant (litter)")
  fam_class[sel] = lab
  lv = c(paste0(vt[i], ": ", named), if(length(cnt) > length(named)) paste0(vt[i], ": other"))
  fam_cols = c(fam_cols, setNames(r[seq_along(lv)], lv), if(any(is.na(smp$family[sel]))) setNames(fam_none, paste0(vt[i], ": no plant (litter)")))
}
# --trees=species: Trees classed by species (genus -> common name, leaf type), fixed shades of the Trees hue
tree_sp = list(Pinus = c("pine", "conifer", "#3f007d"), Abies = c("fir", "conifer", "#5b3a9e"),
               Picea = c("spruce", "conifer", "#7b6bc4"), Populus = c("aspen", "broadleaf", "#cbc0f0"))
if(identical(opt("trees"), "species") && "Trees" %in% smp$veg) {
  sel = smp$veg == "Trees"; genus = sub(" .*", "", smp$species[sel])
  known = genus %in% names(tree_sp)
  fam_class[sel] = ifelse(known, paste0("Trees: ", sapply(genus, function(g) if(g %in% names(tree_sp)) tree_sp[[g]][1] else ""),
                                        " (", sapply(genus, function(g) if(g %in% names(tree_sp)) tree_sp[[g]][2] else ""), ")"),
                          "Trees: other")
  fam_cols = fam_cols[!startsWith(names(fam_cols), "Trees: ")]
  for(g in names(tree_sp)) fam_cols[paste0("Trees: ", tree_sp[[g]][1], " (", tree_sp[[g]][2], ")")] = tree_sp[[g]][3]
  if(any(!known)) fam_cols["Trees: other"] = "#e0dcf2"
}
smp$fam_class = fam_class
if(is.null(site)) {
  write.table(smp, file.path(outdir, paste0(name, ".sample_locations.txt")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
} else {
  if(!(site %in% smp$site)) stop("no samples of site '", site, "'; sites: ", paste(sort(unique(smp$site)), collapse = ", "))
  smp = smp[smp$site == site, ]
}
nocoord = smp$sample[is.na(smp$lat) | is.na(smp$lon)]
m = smp[!(smp$sample %in% nocoord), ]
cat("sample-map.R", script_version, ":", nrow(smp), "samples (", paste(names(table(smp$loci)), table(smp$loci), collapse = ", "),
    ");", length(nocoord), "without coordinates\n")

# ---- locations and sites, in Web Mercator ----
merc = function(lon, lat) data.frame(x = 6378137 * lon * pi / 180, y = 6378137 * log(tan(pi / 4 + lat * pi / 360)))
m = cbind(m, merc(m$lon, m$lat))
years = sort(unique(m$year))
loc = do.call(rbind, lapply(split(m, paste(m$lat, m$lon)), function(z) data.frame(x = z$x[1], y = z$y[1], veg = z$veg[1],
  fam_class = z$fam_class[1],
  code = paste(unique(sub("-.*", "", z$sample)), collapse = "/"),
  sampled = if(length(unique(z$year)) > 1) paste(sort(unique(z$year)), collapse = " + ") else paste(z$year[1], "only"))))
loc$sampled = factor(loc$sampled, levels = intersect(c(paste(years, "only"), paste(years, collapse = " + ")), unique(loc$sampled)))
sites = do.call(rbind, lapply(split(m, m$site), function(z) data.frame(x = mean(range(z$x)), y = mean(range(z$y)),
  text = paste0(z$site[1], " (", nrow(z), ")"))))
# labels: one per site, or (one site) one per location with its sampling code
labs_df = if(is.null(site)) sites else data.frame(x = loc$x, y = loc$y, text = loc$code)
pad = function(r, f) { w = max(diff(r), 200); mean(r) + c(-1, 1) * w * (0.5 + f) }   # at least 200 m wide
xr = pad(range(m$x), 0.08); yr = pad(range(m$y), 0.08)
# frame at least 1:2 (height:width) and at most 2:1, widening the shorter side around the centre
if(diff(yr) < 0.5 * diff(xr)) yr = mean(yr) + c(-1, 1) * 0.25 * diff(xr)
if(diff(xr) < 0.5 * diff(yr)) xr = mean(xr) + c(-1, 1) * 0.25 * diff(yr)
# tile zoom: ~2700 pixels across the map (256-pixel tiles; 156543 m per pixel at zoom 0), at most 18
zoom = if(is.null(zoom)) min(18L, as.integer(floor(log2(2700 * 156543.03 / diff(xr))))) else as.integer(zoom)

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
    cat("base map credit:", maptiles::get_credit(provider), "\n")
    list(df = df[, c("x", "y", "fill")], credit = maptiles::get_credit(provider))
  }, error = function(e) { cat("no base map (", conditionMessage(e), ")\n"); NULL })
  if(is.null(opt("cache"))) unlink(cachedir, recursive = TRUE)   # temporary tile folder
}

# ---- scale bar: bottom left (site maps) or bottom right (all sites), a round length near 1/5 of the map width ----
# Web Mercator stretches distances by 1/cos(latitude): ground metres = map units * cos(lat)
scalebar = NULL
if(TRUE) {
  k = cos(mean(m$lat) * pi / 180)
  target = diff(xr) * k / 5
  nice = as.vector(outer(c(1, 2, 5), 10^(0:5))); len = max(nice[nice <= target])
  seg = len / 2 / k                                                    # map units of one half
  # bottom left on site maps, bottom right on the map of all sites
  x0 = if(is.null(site)) xr[2] - 0.06 * diff(xr) - 2 * seg else xr[1] + 0.04 * diff(xr); y0 = yr[1] + 0.05 * diff(yr); bh = 0.012 * diff(yr)
  unit_lab = function(v) if(len >= 1000) paste(v / 1000) else paste(v)
  scalebar = list(
    annotate("rect", xmin = x0 - 0.015 * diff(xr), xmax = x0 + 2 * seg + 0.045 * diff(xr), ymin = y0 - 0.025 * diff(yr),
             ymax = y0 + bh + 0.07 * diff(yr), fill = "white", alpha = 0.85),
    annotate("rect", xmin = x0 + c(0, seg), xmax = x0 + c(seg, 2 * seg), ymin = y0, ymax = y0 + bh,
             fill = c("black", "white"), colour = "black", linewidth = 0.3),
    annotate("text", x = x0 + c(0, seg, 2 * seg), y = y0 + bh + 0.02 * diff(yr), vjust = 0, size = 3,
             label = c("0", unit_lab(len / 2), paste(unit_lab(len), if(len >= 1000) "km" else "m"))))
}

# ---- slope and contours from USGS 3DEP (--slope) ----
want_slope = "--slope" %in% args
contour_m = as.numeric(opt("contour", "20"))
slope_df = NULL; contour_df = NULL
if(want_slope) slope_df = tryCatch({
  zone = floor((mean(m$lon) + 180) / 6) + 1; utm = paste0("EPSG:", 32600 + zone)
  ext_utm = terra::ext(terra::project(terra::vect(cbind(rep(xr, each = 2), c(yr, rev(yr))), type = "polygons", crs = "EPSG:3857"), utm))
  ext_utm = ext_utm + 50                                                # margin: slope needs neighbours at the edges
  px = 2                                                                # metres per pixel requested
  nx = min(4000, ceiling((ext_utm$xmax - ext_utm$xmin) / px)); ny = min(4000, ceiling((ext_utm$ymax - ext_utm$ymin) / px))
  url = paste0("https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage?bbox=",
               paste(c(ext_utm$xmin, ext_utm$ymin, ext_utm$xmax, ext_utm$ymax), collapse = ","), "&bboxSR=", 32600 + zone,
               "&imageSR=", 32600 + zone, "&size=", nx, ",", ny, "&format=tiff&pixelType=F32",
               "&interpolation=RSP_BilinearInterpolation&f=image")
  tf = tempfile(fileext = ".tif"); utils::download.file(url, tf, mode = "wb", quiet = TRUE)
  dem = terra::rast(tf); terra::crs(dem) = utm
  slp = terra::terrain(dem, v = "slope", unit = "degrees")
  # back to Web Mercator, only the map frame
  slp_m = terra::project(slp, "EPSG:3857", method = "bilinear"); dem_m = terra::project(dem, "EPSG:3857", method = "bilinear")
  fr = terra::ext(xr[1], xr[2], yr[1], yr[2]); slp_m = terra::crop(slp_m, fr); dem_m = terra::crop(dem_m, fr)
  sdf = terra::as.data.frame(slp_m, xy = TRUE, na.rm = TRUE); names(sdf)[3] = "slope"
  if(contour_m > 0) contour_df <<- { d = terra::as.data.frame(dem_m, xy = TRUE, na.rm = TRUE); names(d)[3] = "z"; d }
  # check: slope at the sample locations vs the lidar slope in the sample data, if present
  at = terra::extract(slp, terra::project(terra::vect(cbind(loc$x, loc$y), crs = "EPSG:3857"), utm))[, 2]
  cat(sprintf("slope (3DEP): %.1f-%.1f degrees in the frame; at the %d locations median %.1f\n",
              min(sdf$slope), max(sdf$slope), nrow(loc), median(at, na.rm = TRUE)))
  sdf
}, error = function(e) { cat("no slope (", conditionMessage(e), ")\n"); NULL })

# ---- plot ----
veg_lv = sort(unique(loc$veg))
veg_cols = setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")[seq_along(veg_lv)], veg_lv)  # as in the composition plots
shapes = setNames(c(21, 24, 22, 23)[seq_along(levels(loc$sampled))], levels(loc$sampled))
if(identical(colour_by, "family")) {
  # fill = family class (shade of the vegetation hue), shape = vegetation type, white outline; a small star inside
  # the points of locations sampled in both years (black)
  fc = fam_cols[names(fam_cols) %in% loc$fam_class]
  veg_shape = setNames(c(21, 24, 22, 23, 25)[seq_along(veg_lv)], veg_lv)
  both = levels(loc$sampled)[grepl("\\+", levels(loc$sampled))]
  stars = loc[loc$sampled %in% both, ]
  lum = function(h) { v = grDevices::col2rgb(h) / 255; v = ifelse(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055)^2.4)
    0.2126 * v[1, ] + 0.7152 * v[2, ] + 0.0722 * v[3, ] }
  stars$dy = ifelse(stars$veg == names(veg_shape)[veg_shape == 24][1] & !is.na(names(veg_shape)[veg_shape == 24][1]),
                    -0.12, 0) * diff(yr) / 60            # triangles: centre of mass a little below the point centre
  star_lab = paste0("sampled in ", paste(years, collapse = " and "), " (", nrow(stars), ")")
  points_layer = list(
    geom_point(data = loc, aes(x, y, fill = fam_class, shape = veg), colour = "white", size = 3.1, stroke = 0.5),
    scale_fill_manual(values = fc, breaks = names(fc), name = "Vegetation type: plant family (locations)",
                      labels = paste0(names(fc), " (", as.integer(table(factor(loc$fam_class, levels = names(fc)))), ")"),
                      guide = guide_legend(order = 1, ncol = 2, override.aes = list(shape = unname(veg_shape[sub(":.*", "", names(fc))]),
                                                                           colour = "white", size = 3.9))),
    scale_shape_manual(values = veg_shape, guide = "none"),
    if(nrow(stars)) geom_point(data = stars, aes(x, y + dy), colour = "black", shape = 8, size = 0.65, stroke = 0.35),
    # legend entry for the star: an invisible point mapped to a one-level size scale (size is not used otherwise)
    if(nrow(stars)) geom_point(data = stars[1, ], aes(x, y, size = star_lab), shape = 8, alpha = 0),
    if(nrow(stars)) scale_size_manual(values = setNames(1.6, star_lab), name = "Years (locations)",
                                      guide = guide_legend(order = 2, override.aes = list(alpha = 1, stroke = 0.5))))
} else {
  points_layer = list(
    geom_point(data = loc[order(loc$sampled), ], aes(x, y, fill = veg, shape = sampled), colour = "white", size = 2.6, stroke = 0.5),
    scale_fill_manual(values = veg_cols, name = "Vegetation type", guide = guide_legend(override.aes = list(shape = 21, size = 3.5))),
    scale_shape_manual(values = shapes, name = "Sampled", guide = guide_legend(override.aes = list(fill = "grey50", size = 3.5))))
}
lonb = pretty(range(m$lon), 5); latb = pretty(range(m$lat), 5)
deg = function(v) format(v, nsmall = max(0, -floor(log10(diff(range(v)))) + 1), trim = TRUE)   # enough decimals
g = ggplot() +
  (if(!is.null(bm)) geom_raster(data = bm$df, aes(x, y, fill = I(fill)), show.legend = FALSE) else NULL) +
  (if(!is.null(slope_df)) list(
     ggnewscale::new_scale_fill(),
     geom_raster(data = slope_df, aes(x, y, fill = cut(slope, c(0, 10, 20, 30, 90), right = FALSE,
                                                            labels = c("0-10", "10-20", "20-30", "30+"))), alpha = 0.5),
     scale_fill_manual(values = c("#f7fbff", "#9ecae1", "#3182bd", "#08306b"), name = "Slope (degrees)", drop = FALSE,
                       guide = guide_legend(order = 3, ncol = 1, override.aes = list(alpha = 0.8))),
     ggnewscale::new_scale_fill()) else NULL) +
  (if(!is.null(contour_df)) geom_contour(data = contour_df, aes(x, y, z = z), breaks = seq(0, 9000, contour_m),
                                          colour = "white", linewidth = 0.25, alpha = 0.8) else NULL) +
  points_layer +
  ggrepel::geom_label_repel(data = labs_df, aes(x, y, label = text), size = if(is.null(site)) 3.2 else 2.4, label.size = 0, fill = "white", alpha = 0.9,
                            label.padding = unit(2.5, "pt"), min.segment.length = 0, segment.colour = "white",
                            segment.size = if(is.null(site)) 0.5 else 0.3, box.padding = if(is.null(site)) 1 else 0.35, seed = 1, max.overlaps = Inf) +
  scale_x_continuous(breaks = merc(lonb, 0)$x, labels = paste0(deg(abs(lonb)), "°W"), expand = c(0, 0)) +
  scale_y_continuous(breaks = merc(0, latb)$y, labels = paste0(deg(latb), "°N"), expand = c(0, 0)) +
  scalebar +
  coord_fixed(xlim = xr, ylim = yr) +
  labs(x = NULL, y = NULL,
       title = paste0(if(is.null(site)) name else site, " sampling locations (", nrow(loc), ")")) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), legend.key = element_rect(fill = NA)) +
  (if(identical(colour_by, "family")) theme(legend.position = "bottom", legend.box = "horizontal",
                                            legend.box.just = "top", legend.justification = "left",
                                            legend.title.position = "top", legend.title = element_text(size = 9),
                                            legend.text = element_text(size = 8), legend.key.size = unit(11, "pt"),
                                            legend.key.spacing.y = unit(1, "pt"), legend.spacing.x = unit(18, "pt"),
                                            legend.margin = margin(0, 0, 0, 0)) else NULL)
asp = diff(yr) / diff(xr)
f = file.path(outdir, paste0(name, if(!is.null(site)) paste0(".", gsub("[^A-Za-z0-9]+", "_", site)) else "",
                              if(!is.null(slope_df)) ".slope_map" else ".sample_map"))
h = min(14, 7.2 * asp + 0.8)
if(identical(colour_by, "family")) h = min(14, 9 * asp + 0.8 + 1.3)   # map full width, legends below
ggsave(paste0(f, ".png"), g, width = 9, height = h, dpi = 300)
if(write_pdf) ggsave(paste0(f, ".pdf"), g, width = 9, height = h)
cat("sample map:", nrow(m), "samples,", nrow(loc), "locations,", nrow(sites), "sites ->", paste0(f, if(write_pdf) ".png / .pdf" else ".png"), "\n")
