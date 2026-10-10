# =============================================================================
# STEP 10: METHODS SCHEMATICS FOR SLIDES
#
# Two slide-width figures drawn on the study provinces exactly as analysed
# (Basin and Range, Coastal Plain, Great Plains with the Central Lowland west
# of the Mississippi; Fenneman & Johnson 1946 boundaries from Step 3):
#   methods_workflow.png/.svg/.pdf   five panels: fossil occurrences (FAUNMAP
#                                    and PBDB logos) -> spatial bins -> time
#                                    bins -> province x species matrix ->
#                                    multisite beta_SIM equation
#   beta_sim_schematic.png/.svg/.pdf the three provinces on a map, linked
#                                    pairwise, beside the beta_SIM equation
#
# Logos: save the FAUNMAP and PBDB logos as PNG files and set 'faunmap_logo'
# and 'pbdb_logo' below (paths relative to work_dir or absolute). A missing
# logo is replaced by the database name in plain type.
#
# Font: Helvetica ('font_family'). Windows has no Helvetica, so the PNG and PDF
# use Arial there (the same letter widths); the SVG asks for Helvetica, which
# PowerPoint / Illustrator on a Mac will use.
#
# Inputs: Outputs/3_spatial/spatial_layers.rds (Step 3) and the Mississippi
#         line cached in Outputs/maps/basemap (Step 5b; downloaded if missing).
# Outputs: Outputs/10_methods_schematics/
# Needs: ggplot2, sf, maps, png; svglite for the SVG files (optional).
# =============================================================================

for (pkg in c("ggplot2", "sf", "maps", "png")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}
library(grid)
library(ggplot2)
library(sf)
suppressMessages(sf_use_s2(FALSE))

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

layers_file <- file.path("Outputs", "3_spatial", "spatial_layers.rds")
basemap_dir <- file.path(work_dir, "Outputs", "maps", "basemap")
output_dir  <- file.path(work_dir, "Outputs", "10_methods_schematics")
faunmap_logo <- "faunmap_logo.png"     # PNG files (relative to work_dir or absolute)
pbdb_logo    <- "pbdb_logo.png"

regions <- c("Basin and Range", "Great Plains", "Coastal Plain")       # map order 1, 2, 3
western_central_lowland <- TRUE
rivers_url <- paste0("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/",
                     "geojson/ne_10m_rivers_lake_centerlines.geojson")

font_family <- "Helvetica"
region_fill <- c("Basin and Range" = "#E69F00", "Coastal Plain" = "#56B4E9", "Great Plains" = "#009E73")
region_alpha <- 0.6
ink <- "#1E2A3E"; soft <- "#5F6B80"; faint <- "#C9CFD8"
land_fill <- "#F4EFE4"; sea_fill <- "#E4EEF7"; border_col <- "#C9BFA9"; river_col <- "#2C6FB7"
map_crs <- "+proj=laea +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
margin_deg <- 1.5
fig_dpi <- 300

# =============================================================================
# 1. STUDY PROVINCES
# =============================================================================

setwd(work_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (.Platform$OS.type == "windows" && font_family == "Helvetica") {
  grDevices::windowsFonts(Helvetica = grDevices::windowsFont("Arial"))
}

polygons_only <- function(x) {
  if (any(st_geometry_type(x) == "GEOMETRYCOLLECTION")) x <- st_collection_extract(x, "POLYGON", warn = FALSE)
  x <- x[st_geometry_type(x) %in% c("POLYGON", "MULTIPOLYGON"), ]
  st_cast(x, "MULTIPOLYGON", warn = FALSE)
}
repair_geom <- function(x) {
  x <- suppressWarnings(suppressMessages(polygons_only(st_make_valid(x))))
  bad <- !st_is_valid(x)
  if (any(bad, na.rm = TRUE)) {
    st_geometry(x)[which(bad)] <- suppressWarnings(st_buffer(st_geometry(x)[which(bad)], 0))
    x <- suppressWarnings(suppressMessages(polygons_only(st_make_valid(x))))
  }
  x[!st_is_empty(x), ]
}
safe_op <- function(f, a, b) tryCatch(suppressWarnings(suppressMessages(f(a, b))),
                                      error = function(e) suppressWarnings(suppressMessages(
                                        f(st_buffer(st_make_valid(a), 0), st_buffer(st_make_valid(b), 0)))))
union_units <- function(prov, u) {
  x <- prov[prov$Unit %in% u, ]
  if (!nrow(x)) return(NULL)
  suppressWarnings(suppressMessages(st_union(st_make_valid(x))))
}

if (!file.exists(layers_file)) stop("File not found: ", file.path(work_dir, layers_file), "\nRun Step 3 first.")
prov <- repair_geom(readRDS(layers_file)$province) |> st_transform(4326)
missing_units <- setdiff(regions, prov$Unit)
if (length(missing_units)) stop("Province(s) not in the Step 3 layers: ", paste(missing_units, collapse = ", "))
geoms <- lapply(regions, function(r) union_units(prov, r)); names(geoms) <- regions

dir.create(basemap_dir, recursive = TRUE, showWarnings = FALSE)
rivers_file <- file.path(basemap_dir, basename(rivers_url))
if (!file.exists(rivers_file)) {
  cat("  Downloading Natural Earth rivers (one time, ~7 MB)...\n")
  download.file(rivers_url, rivers_file, mode = "wb", quiet = TRUE)
}
river <- st_read(rivers_file, quiet = TRUE)
river <- river[river$name %in% "Mississippi", ]

if (western_central_lowland) {
  aea <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +datum=WGS84 +units=m +no_defs"
  rxy <- st_coordinates(suppressWarnings(st_transform(st_segmentize(st_transform(river, aea), 1000), 4326)))[, c("X", "Y")]
  bands <- seq(floor(min(rxy[, 2]) * 20) / 20, ceiling(max(rxy[, 2]) * 20) / 20, by = 0.05)
  band_lon <- vapply(bands, function(b) { inb <- abs(rxy[, 2] - b) <= 0.05
                                          if (any(inb)) min(rxy[inb, 1]) else NA_real_ }, numeric(1))
  river_lon_at <- approxfun(bands[!is.na(band_lon)], band_lon[!is.na(band_lon)], rule = 2)
  lat_seq <- seq(20, 60, by = 0.05)
  west_poly <- st_sfc(st_polygon(list(rbind(cbind(river_lon_at(lat_seq), lat_seq),
                                            c(-150, max(lat_seq)), c(-150, min(lat_seq)),
                                            c(river_lon_at(min(lat_seq)), min(lat_seq))))), crs = 4326)
  cl <- union_units(prov, "Central Lowland")
  if (!is.null(cl)) geoms[["Great Plains"]] <- st_union(c(geoms[["Great Plains"]], safe_op(st_intersection, cl, west_poly)))
}
study <- st_sf(Region = factor(regions, levels = regions), geometry = do.call(c, unname(geoms)), crs = 4326)

bb <- st_bbox(study)
frame_ll <- st_as_sfc(st_bbox(c(xmin = bb[["xmin"]] - margin_deg, xmax = bb[["xmax"]] + margin_deg,
                                ymin = bb[["ymin"]] - margin_deg, ymax = bb[["ymax"]] + margin_deg), crs = 4326))
clip_frame <- function(x) suppressWarnings(suppressMessages(
  st_intersection(st_set_agr(st_make_valid(x), "constant"), st_buffer(frame_ll, 3))))
countries <- st_as_sf(maps::map("world", regions = c("Canada", "USA", "Mexico"), fill = TRUE, plot = FALSE))
st_crs(countries) <- NA; st_crs(countries) <- 4326
countries_m <- clip_frame(countries) |> st_transform(map_crs)
prov_m  <- clip_frame(prov) |> st_transform(map_crs)
study_m <- st_transform(study, map_crs)
river_m <- suppressWarnings(st_transform(clip_frame(river), map_crs))
lims <- st_bbox(st_transform(st_segmentize(frame_ll, 0.25), map_crs))

# One label point inside each province: the centre of the largest circle that fits in its
# largest part (falls back to a point on the surface).
label_point <- function(g) {
  parts <- st_cast(st_sfc(g, crs = map_crs), "POLYGON")
  big <- parts[which.max(st_area(parts))]
  pt <- tryCatch(st_centroid(st_inscribed_circle(big, nQuadSegs = 8)), error = function(e) st_point_on_surface(big))
  st_coordinates(pt)[1, 1:2]
}
nodes <- do.call(rbind, lapply(seq_along(regions), function(i) label_point(st_geometry(study_m)[[i]])))
nodes <- data.frame(Region = regions, x = nodes[, 1], y = nodes[, 2], lab = as.character(seq_along(regions)))

base_map <- function(alpha = region_alpha, outline = 0.4) {
  ggplot() +
    geom_sf(data = countries_m, fill = land_fill, colour = NA) +
    geom_sf(data = prov_m, fill = NA, colour = border_col, linewidth = 0.2) +
    geom_sf(data = study_m, fill = land_fill, colour = NA) +
    geom_sf(data = study_m, aes(fill = Region), colour = ink, linewidth = outline, alpha = alpha) +
    geom_sf(data = countries_m, fill = NA, colour = soft, linewidth = 0.3) +
    geom_sf(data = river_m, colour = river_col, linewidth = 0.5) +
    scale_fill_manual(values = region_fill, guide = "none")
}
map_coord <- function() coord_sf(crs = map_crs, xlim = lims[c("xmin", "xmax")], ylim = lims[c("ymin", "ymax")],
                                 expand = FALSE, datum = NA)

beta_eq  <- quote(frac(sum(min*group("(", list(b[ij], b[ji]), ")"), "i<j", ""),
                       (sum(S[i], i, "") - S[T]) + sum(min*group("(", list(b[ij], b[ji]), ")"), "i<j", "")))
txt <- function(lab, x, y, size, col = ink, face = "plain", just = "left", lineheight = 1.2)
  grid.text(lab, x, y, just = just, gp = gpar(fontfamily = font_family, fontsize = size, col = col,
                                              fontface = face, lineheight = lineheight))

read_logo <- function(path) {
  if (is.na(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) { message("  Logo not found: ", path, " (using the name instead)"); return(NULL) }
  png::readPNG(path)
}
logos <- list(FAUNMAP = read_logo(faunmap_logo), PBDB = read_logo(pbdb_logo))

save_all <- function(draw_fun, name, w, h) {
  png(file.path(output_dir, paste0(name, ".png")), width = w, height = h, units = "in", res = fig_dpi,
      type = "cairo")
  draw_fun(); invisible(dev.off())
  grDevices::cairo_pdf(file.path(output_dir, paste0(name, ".pdf")), width = w, height = h)
  draw_fun(); invisible(dev.off())
  if (requireNamespace("svglite", quietly = TRUE)) {
    svg_path <- file.path(output_dir, paste0(name, ".svg"))
    svglite::svglite(svg_path, width = w, height = h)
    draw_fun(); invisible(dev.off())
    svg <- readLines(svg_path, warn = FALSE)     # name the font as Helvetica in the SVG
    svg <- gsub("font-family: *[^;\"']+", paste0("font-family: ", font_family, ", Arial, sans-serif"), svg)
    writeLines(svg, svg_path)
  }
  cat("  Saved", name, "\n")
}

# =============================================================================
# 2. METHODS WORKFLOW
# =============================================================================

glyph_logos <- function() {
  for (k in 1:2) {
    x <- c(0.27, 0.73)[k]; nm <- names(logos)[k]; img <- logos[[k]]
    if (is.null(img)) {
      grid.roundrect(x, 0.55, 0.42, 0.42, r = unit(0.08, "snpc"), gp = gpar(fill = "#EEF2F6", col = faint))
      txt(nm, x, 0.55, 15, ink, "bold", "centre")
    } else {
      asp <- dim(img)[2] / dim(img)[1]
      pushViewport(viewport(x, 0.55, 0.44, 0.6))
      vw <- convertWidth(unit(1, "npc"), "in", TRUE); vh <- convertHeight(unit(1, "npc"), "in", TRUE)
      hh <- min(vh, vw / asp)
      grid.raster(img, width = unit(hh * asp, "in"), height = unit(hh, "in"), interpolate = TRUE)
      popViewport()
    }
    txt(nm, x, 0.12, 10.5, soft, just = "centre")
  }
}
glyph_map <- function() {
  bbs <- st_bbox(study_m); pad <- 0.04 * (bbs[["xmax"]] - bbs[["xmin"]])
  p <- base_map(alpha = 0.75, outline = 0.3) +
    coord_sf(crs = map_crs, xlim = c(bbs[["xmin"]] - pad, bbs[["xmax"]] + pad),
             ylim = c(bbs[["ymin"]] - pad, bbs[["ymax"]] + pad), expand = FALSE, datum = NA) +
    theme_void()
  print(p, newpage = FALSE)
}
glyph_time <- function() {
  b <- c(3.25, 2.5, 1.75, 1.0, 0.25, 0)
  sx <- function(a) 0.08 + (3.25 - a) / 3.25 * 0.84
  for (i in 1:5) grid.rect(sx(b[i]), 0.55, sx(b[i + 1]) - sx(b[i]), 0.26, just = c("left", "centre"),
                           gp = gpar(fill = if (i %% 2) "#E9DEC6" else "#F6F1E6", col = "#FFFFFF", lwd = 1.5))
  grid.lines(c(0.08, 0.92), c(0.36, 0.36), gp = gpar(col = soft, lwd = 1))
  for (a in c(3.25, 2, 1, 0)) {
    grid.lines(c(sx(a), sx(a)), c(0.36, 0.33), gp = gpar(col = soft, lwd = 1))
    txt(format(a), sx(a), 0.24, 10, soft, just = "centre")
  }
  txt("Ma", 0.92, 0.12, 9.5, soft, just = "right")
}
glyph_matrix <- function() {
  set.seed(7); nc <- 9; m <- matrix(runif(3 * nc) < 0.5, 3)
  cols <- region_fill[regions]
  for (i in 1:3) for (j in 1:nc) {
    grid.circle(0.20 + (j - 1) * 0.085, 0.72 - (i - 1) * 0.2, unit(4.2, "pt"),
                gp = gpar(fill = if (m[i, j]) cols[i] else "#FFFFFF", col = if (m[i, j]) NA else faint, lwd = 1))
  }
  txt("species \u2192", 0.20, 0.18, 9.5, soft)
}
glyph_beta <- function() {
  grid.text(expression(beta[SIM] == ""), 0.02, 0.86, just = "left", gp = gpar(fontfamily = font_family, fontsize = 14, col = ink))
  grid.text(as.expression(beta_eq), 0.02, 0.46, just = "left", gp = gpar(fontfamily = font_family, fontsize = 12, col = ink))
}

panels <- list(
  list("Fossil occurrences", "FAUNMAP and PBDB records,\nlinked and taxonomically cleaned", glyph_logos),
  list("Spatial bins", "Three physiographic provinces\n(Fenneman & Johnson 1946)", glyph_map),
  list("Time bins", "Five 0.75-Myr bins,\n3.25 Ma to 11.7 ka", glyph_time),
  list("Province \u00d7 species", "Large mammals (>1 kg),\npresence\u2013absence per bin", glyph_matrix),
  list("Beta diversity", "Multisite Simpson dissimilarity\n(Baselga 2010), equal species", glyph_beta))

draw_workflow <- function() {
  grid.newpage(); grid.rect(gp = gpar(fill = "#FFFFFF", col = NA))
  n <- length(panels); left <- 0.025; gap <- 0.03; pw <- (1 - 2 * left - gap * (n - 1)) / n
  for (i in seq_len(n)) {
    p <- panels[[i]]; x0 <- left + (i - 1) * (pw + gap)
    pushViewport(viewport(x = x0, y = 0.42, width = pw, height = 0.54, just = c("left", "bottom")))
    p[[3]](); popViewport()
    txt(p[[1]], x0, 0.31, 15, ink, "bold")
    txt(p[[2]], x0, 0.18, 12, soft)
    if (i < n) grid.lines(c(x0 + pw + 0.005, x0 + pw + gap - 0.005), c(0.70, 0.70),
                          arrow = arrow(length = unit(0.09, "in"), type = "open"), gp = gpar(col = soft, lwd = 1.3))
  }
}
save_all(draw_workflow, "methods_workflow", 13.33, 3.5)

# =============================================================================
# 3. beta_SIM SCHEMATIC
# =============================================================================

shrink <- function(a, b, d = 2.1e5) { v <- b - a; a + v * d / sqrt(sum(v^2)) }
pairs <- data.frame(from = c(1, 2, 1), to = c(2, 3, 3), curv = c(-0.18, -0.18, 0.22))
seg <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(k) {
  a <- unlist(nodes[pairs$from[k], c("x", "y")]); b <- unlist(nodes[pairs$to[k], c("x", "y")])
  s <- shrink(a, b); e <- shrink(b, a)
  data.frame(x = s[1], y = s[2], xend = e[1], yend = e[2])
}))
map_plot <- base_map()
for (k in seq_len(nrow(seg)))
  map_plot <- map_plot + geom_curve(data = seg[k, ], aes(x = x, y = y, xend = xend, yend = yend),
                                    curvature = pairs$curv[k], colour = ink, linewidth = 0.7,
                                    arrow = arrow(length = unit(0.11, "in"), ends = "both", type = "closed"))
map_plot <- map_plot +
  geom_point(data = nodes, aes(x, y, fill = Region), shape = 21, size = 13, colour = "#FFFFFF", stroke = 1.6) +
  geom_text(data = nodes, aes(x, y, label = lab), family = font_family, fontface = "bold", size = 6.2, colour = ink) +
  map_coord() + theme_void() +
  theme(panel.background = element_rect(fill = sea_fill, colour = soft, linewidth = 0.4))

draw_beta <- function() {
  grid.newpage(); grid.rect(gp = gpar(fill = "#FFFFFF", col = NA))
  pushViewport(viewport(x = 0.01, y = 0.58, width = 0.56, height = 0.82, just = c("left", "centre")))
  print(map_plot, newpage = FALSE); popViewport()
  for (i in seq_along(regions)) {
    kx <- 0.04 + (i - 1) * 0.165
    grid.circle(kx, 0.08, unit(9, "pt"), gp = gpar(fill = region_fill[regions[i]], col = NA))
    txt(i, kx, 0.08, 11, ink, "bold", "centre")
    txt(regions[i], kx + 0.017, 0.08, 12.5, soft)
  }
  x0 <- 0.61
  txt("Each pair of provinces is compared,", x0, 0.86, 15)
  txt("then all three are combined into one value", x0, 0.80, 15)
  grid.text(as.expression(substitute(beta[SIM] == E, list(E = beta_eq))), x0, 0.60, just = "left",
            gp = gpar(fontfamily = font_family, fontsize = 19, col = ink))
  keys <- list(list(quote(S[i]), "species in province i"),
               list(quote(S[T]), "species across all three provinces"),
               list(quote(b[ij]), "species in i but not in j"),
               list(quote(min*group("(", list(b[ij], b[ji]), ")")), "species replaced between i and j (one arrow)"),
               list(quote(sum(S[i], i, "") - S[T]), "species shared among provinces"))
  for (k in seq_along(keys)) {
    y <- 0.42 - (k - 1) * 0.058
    grid.text(as.expression(keys[[k]][[1]]), x0, y, just = "left", gp = gpar(fontfamily = font_family, fontsize = 13, col = ink))
    txt(keys[[k]][[2]], x0 + 0.115, y, 12.5, soft)
  }
  txt("0 = same species everywhere (homogeneous)", x0, 0.11, 12.5, face = "bold")
  txt("1 = no species shared (fully distinct)", x0, 0.06, 12.5, face = "bold")
}
save_all(draw_beta, "beta_sim_schematic", 13.33, 6)

cat("\nDone. Figures in", output_dir, "\n")
