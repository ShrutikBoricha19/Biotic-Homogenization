# =============================================================================
# STEP 4 ADDENDUM: SPATIOTEMPORAL MAPS OF THE THREE STUDY PROVINCES
#
# One map per time bin, in one figure: 3 maps on the top row, 2 on the bottom,
# the legend in the sixth place. Only the study provinces are highlighted -
# Basin and Range, Coastal Plain and Great Plains (with the Central Lowland
# west of the Mississippi, as in Step 5b and every later step) - and only
# the sites in those provinces are shown, coloured by province. Every other
# province is drawn as plain land.
#
# Sites: one point per site (locality / collection) and time bin, from the
# Step 5b master data, so the sites are exactly those used in Steps 6-8.
#
# Outputs (Outputs/4_maps/4_study_provinces/):
#   map_study_provinces_by_bin.png / .pdf    the 5-panel figure
#   study_province_sites_per_bin.csv         sites per province and time bin
#
# Inputs: Outputs/5b_great_plains/master_data_unique.csv (Step 5b),
#         Outputs/3d_resolved/time_bins.csv,
#         Outputs/3_spatial/spatial_layers.rds (Step 3, province outlines),
#         the Mississippi line cached by Step 5b in Outputs/maps/basemap
#         (downloaded once if it is missing).
# Run the whole file (Ctrl+Shift+S in RStudio) AFTER Step 5b (it needs the
# Great Plains assignment made there).
# Needs: dplyr, tidyr, ggplot2, sf, maps, patchwork.
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "sf", "maps", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)

suppressMessages(sf_use_s2(FALSE))

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

master_file <- file.path("Outputs", "5b_great_plains", "master_data_unique.csv")
bins_file   <- file.path("Outputs", "3d_resolved", "time_bins.csv")
layers_file <- file.path("Outputs", "3_spatial", "spatial_layers.rds")
basemap_dir <- file.path(work_dir, "Outputs", "maps", "basemap")      # Step 5b's download cache
output_dir  <- file.path(work_dir, "Outputs", "4_maps", "4_study_provinces")

regions <- c("Basin and Range", "Coastal Plain", "Great Plains")
western_central_lowland <- TRUE   # Great Plains outline includes the Central Lowland west of
                                  # the Mississippi (as Step 5b); FALSE = the Great Plains province only

rivers_url   <- paste0("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/",
                       "geojson/ne_10m_rivers_lake_centerlines.geojson")
rivers_local <- NA                # or the path of a local copy of that file
show_river   <- TRUE              # draw the Mississippi

# Colourblind-friendly (Okabe-Ito); pale tints for the provinces, full colours for the sites.
region_fill <- c("Basin and Range" = "#F6DCA8", "Coastal Plain" = "#C9E4F6", "Great Plains" = "#B9E2CF")
site_fill   <- c("Basin and Range" = "#E69F00", "Coastal Plain" = "#0072B2", "Great Plains" = "#009E73")
ink         <- "#1F2A44"          # text, outlines and site borders (dark navy)
ink_soft    <- "#3E4C6D"
land_fill   <- "#F4EFE4"          # land outside the study provinces
sea_fill    <- "#E4EEF7"
border_col  <- "#B9AE97"          # boundaries of the other provinces
river_col   <- "#2C6FB7"

map_crs <- "+proj=laea +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
margin_deg <- 1.5                 # space around the study provinces (degrees)
site_size  <- 2.2
fig_w <- 16; fig_h <- 9.4; fig_dpi <- 300; base_size <- 13

# =============================================================================
# HELPERS
# =============================================================================

read_text_csv <- function(path) {
  full <- file.path(work_dir, path)
  if (!file.exists(full)) stop("File not found:\n  ", full, "\nRun the earlier steps first (this addendum runs after Step 5b).")
  df <- read.csv(full, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
                 na.strings = c("", "NA"), encoding = "UTF-8")
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))
  cat(sprintf("  %-50s %7d rows\n", path, nrow(df)))
  df
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
fmt_age <- function(x) ifelse(x < 0.1, sprintf("%.4f", x), sprintf("%.2f", x))

# =============================================================================
# 1. READ INPUTS
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
md <- read_text_csv(master_file)
if (!"Spatial_Bin" %in% names(md)) stop("No Spatial_Bin column in ", master_file, ". Run Step 5b again.")
bins <- read_text_csv(bins_file) %>%
  transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma), Younger = as.numeric(Younger_Ma)) %>%
  arrange(Time_Bin)
stale <- setdiff(unique(as.integer(md$Time_Bin)), bins$Time_Bin)
if (length(stale)) stop("The master data holds time bin(s) ", paste(stale, collapse = ", "),
                        " that time_bins.csv does not. Rerun Steps 3c-5b first.")

layers_path <- file.path(work_dir, layers_file)
if (!file.exists(layers_path)) stop("File not found:\n  ", layers_path, "\nRun Step 3 again.")
prov <- repair_geom(readRDS(layers_path)$province) %>% st_transform(4326)
missing_units <- setdiff(regions, prov$Unit)
if (length(missing_units)) stop("Province(s) not in the Step 3 layers: ", paste(missing_units, collapse = ", "))

# Sites: one point per site and time bin, in the study provinces only.
sites <- md %>%
  mutate(Time_Bin = as.integer(Time_Bin), Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) %>%
  filter(Spatial_Bin %in% regions, !is.na(Latitude), !is.na(Longitude), Time_Bin %in% bins$Time_Bin) %>%
  distinct(Site_Key, Time_Bin, .keep_all = TRUE) %>%
  mutate(Region = factor(Spatial_Bin, levels = regions))
cat(sprintf("  %d site x bin points in the three provinces\n", nrow(sites)))

# =============================================================================
# 2. PROVINCE OUTLINES (Great Plains + Central Lowland west of the Mississippi)
# =============================================================================

cat("\n=== 2. PROVINCE OUTLINES ===\n")
geoms <- lapply(regions, function(r) union_units(prov, r)); names(geoms) <- regions
river <- NULL
if (western_central_lowland || show_river) {
  dir.create(basemap_dir, recursive = TRUE, showWarnings = FALSE)
  rivers_file <- if (!is.na(rivers_local)) rivers_local else file.path(basemap_dir, basename(rivers_url))
  if (!file.exists(rivers_file)) {
    cat("  Downloading Natural Earth rivers (one time, ~7 MB)...\n")
    ok <- tryCatch({ download.file(rivers_url, rivers_file, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (!ok) {
      if (file.exists(rivers_file)) file.remove(rivers_file)
      stop("Download failed. Download ne_10m_rivers_lake_centerlines.geojson from\n  ", rivers_url,
           "\nand set 'rivers_local' to its path (or set western_central_lowland and show_river to FALSE).")
    }
  }
  river <- st_read(rivers_file, quiet = TRUE)
  river <- river[river$name %in% "Mississippi", ]
  if (!nrow(river)) stop("No Mississippi in ", rivers_file)
}
if (western_central_lowland) {
  # River longitude at each latitude (westernmost point per 0.05-degree band, as Step 5b).
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
  if (!is.null(cl)) {
    geoms[["Great Plains"]] <- st_union(c(geoms[["Great Plains"]], safe_op(st_intersection, cl, west_poly)))
    cat("  Great Plains outline includes the Central Lowland west of the Mississippi\n")
  }
}
study <- st_sf(Region = factor(regions, levels = regions), geometry = do.call(c, unname(geoms)), crs = 4326)

# Map frame: the study provinces and their sites, plus a margin.
bb <- st_bbox(st_union(c(st_geometry(study), st_geometry(st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326)))))
frame_ll <- st_as_sfc(st_bbox(c(xmin = bb[["xmin"]] - margin_deg, xmax = bb[["xmax"]] + margin_deg,
                                ymin = bb[["ymin"]] - margin_deg, ymax = bb[["ymax"]] + margin_deg), crs = 4326))
clip_frame <- function(x) suppressWarnings(suppressMessages(
  st_intersection(st_set_agr(st_make_valid(x), "constant"), st_buffer(frame_ll, 3))))
countries <- st_as_sf(maps::map("world", regions = c("Canada", "USA", "Mexico"), fill = TRUE, plot = FALSE))
st_crs(countries) <- NA; st_crs(countries) <- 4326
countries_m <- clip_frame(countries) %>% st_transform(map_crs)
prov_m      <- clip_frame(prov) %>% st_transform(map_crs)
study_m     <- st_transform(study, map_crs)
river_m     <- if (show_river) suppressWarnings(st_transform(clip_frame(river), map_crs)) else NULL
lims <- st_bbox(st_transform(st_segmentize(frame_ll, 0.25), map_crs))
sites_m <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326) %>% st_transform(map_crs)

# =============================================================================
# 3. ONE MAP PER TIME BIN
# =============================================================================

cat("\n=== 3. MAPS ===\n")
counts <- sites %>% count(Time_Bin, Region, name = "n_sites", .drop = FALSE) %>%
  tidyr::complete(Time_Bin = bins$Time_Bin, Region, fill = list(n_sites = 0L))
abbr <- c("Basin and Range" = "BR", "Coastal Plain" = "CP", "Great Plains" = "GP")

bin_map <- function(b) {
  k <- which(bins$Time_Bin == b)
  pts <- sites_m[sites_m$Time_Bin == b, ]
  cn <- counts %>% filter(Time_Bin == b)
  sub <- paste(sprintf("%s %d", abbr[as.character(cn$Region)], cn$n_sites), collapse = "  \u00b7  ")
  p <- ggplot() +
    geom_sf(data = countries_m, fill = land_fill, colour = NA) +
    geom_sf(data = prov_m, fill = NA, colour = border_col, linewidth = 0.2) +
    geom_sf(data = study_m, aes(fill = Region), colour = ink_soft, linewidth = 0.35) +
    geom_sf(data = countries_m, fill = NA, colour = ink_soft, linewidth = 0.35)
  if (!is.null(river_m)) p <- p + geom_sf(data = river_m, colour = river_col, linewidth = 0.6)
  p +
    geom_sf(data = pts, aes(colour = Region), size = site_size, shape = 16) +
    geom_sf(data = pts, colour = ink, size = site_size, shape = 1, stroke = 0.35) +          # navy outline
    scale_fill_manual(values = region_fill, name = "Study province", drop = FALSE) +
    scale_colour_manual(values = site_fill, name = "Fossil sites", drop = FALSE) +
    guides(fill = guide_legend(order = 1, override.aes = list(colour = ink_soft)),
           colour = guide_legend(order = 2, override.aes = list(size = 4.5))) +
    coord_sf(crs = map_crs, xlim = lims[c("xmin", "xmax")], ylim = lims[c("ymin", "ymax")], expand = FALSE) +
    labs(title = sprintf("Bin %d  \u00b7  %s\u2013%s Ma", b, fmt_age(bins$Older[k]), fmt_age(bins$Younger[k])),
         subtitle = sprintf("%d sites:  %s", nrow(pts), sub)) +
    theme_void(base_size = base_size) +
    theme(panel.background = element_rect(fill = sea_fill, colour = ink_soft, linewidth = 0.4),
          plot.title = element_text(face = "bold", colour = ink, size = base_size * 1.15, margin = margin(b = 2)),
          plot.subtitle = element_text(colour = ink_soft, size = base_size * 0.85, margin = margin(b = 4)),
          plot.margin = margin(6, 8, 6, 8))
}
maps <- lapply(bins$Time_Bin, bin_map)

fig <- wrap_plots(c(maps, list(guide_area())), ncol = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.title = element_text(face = "bold", colour = ink, size = base_size),
        legend.text = element_text(colour = ink, size = base_size * 0.95),
        legend.key.size = unit(0.9, "cm"), legend.box = "vertical", legend.box.just = "left")
fig <- fig + plot_annotation(
  title = "Fossil sites in the three study provinces through time",
  subtitle = "Large-mammal localities and collections (FAUNMAP and PBDB) per 0.75-Myr time bin",
  caption = paste0("BR = Basin and Range, CP = Coastal Plain, GP = Great Plains",
                   if (western_central_lowland) " (with the Central Lowland west of the Mississippi)" else "",
                   if (show_river) "; blue line = Mississippi River" else "",
                   ". One point per site and time bin (sites as in Step 5b).\n",
                   "Province boundaries: Fenneman & Johnson (1946)",
                   if (show_river) "; river: Natural Earth." else "."),
  theme = theme(plot.title = element_text(face = "bold", colour = ink, size = base_size * 1.5),
                plot.subtitle = element_text(colour = ink_soft, size = base_size * 1.05),
                plot.caption = element_text(colour = ink_soft, size = base_size * 0.75, hjust = 0)))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "map_study_provinces_by_bin.png"), fig, width = fig_w, height = fig_h, dpi = fig_dpi, bg = "white")
ggsave(file.path(output_dir, "map_study_provinces_by_bin.pdf"), fig, width = fig_w, height = fig_h, bg = "white",
       device = if (capabilities("cairo")) cairo_pdf else pdf)
cat("  ", file.path(output_dir, "map_study_provinces_by_bin.png"), " (+ .pdf)\n", sep = "")

out <- counts %>% left_join(bins, by = "Time_Bin") %>%
  transmute(Time_Bin, Older_Ma = Older, Younger_Ma = Younger, Province = as.character(Region), n_sites) %>%
  arrange(Time_Bin, Province)
write.csv(out, file.path(output_dir, "study_province_sites_per_bin.csv"), row.names = FALSE)
cat("  ", file.path(output_dir, "study_province_sites_per_bin.csv"), "\n", sep = "")
cat("\n=== STEP 4 ADDENDUM COMPLETE ===\n")
