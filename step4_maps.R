# =============================================================================
# STEP 4: CONFERENCE MAPS - SPATIAL, TEMPORAL AND SPATIO-TEMPORAL
#
# Maps of the sites used in the analyses: every FAUNMAP locality and PBDB
# collection that has a time bin (Step 3c), a spatial unit (Step 3) and, after
# the sp./cf. resolution of Step 3d, at least one occurrence identified to
# species (sites left with only family-level records are not drawn; set
# 'only_sites_with_species' to FALSE to draw them).
# Hawaii lies outside the map frame and Alaska is not drawn as a unit; sites
# in Alaska, Hawaii and Canada north of 60 degrees N were removed in Step 3.
# Every map has a numbered key beside it (unit number - name - number of sites).
# Sites are drawn alike whatever their database (FAUNMAP or PBDB).
#
# Outputs (Outputs/4_maps/), each as .png (300 dpi) and .pdf:
#   1_spatial/
#     map_political    states and provinces          } all sites of all
#     map_division     physiographic divisions       } time bins
#     map_province     physiographic provinces       }
#   2_temporal/
#     map_bin1 ... map_binN   one map per time bin on the physiographic
#                             provinces, key with that bin's site counts
#     temporal_all_bins       all bins side by side (small multiples)
#   3_spatiotemporal/
#     map_spatiotemporal_province / _division / _political
#                             sites from all time bins, coloured by bin
#     sites_per_unit_and_bin  heat map: sites per province and time bin
#     sites_per_unit_and_bin_<level>.csv   the counts behind it
#
# Inputs: Outputs/3d_resolved/ (Step 3d): site_index.csv, time_bins.csv,
#           faunmap_fauna.csv, pbdb_occurrences.csv
#         Outputs/3_spatial/spatial_layers.rds (Step 3).
# Run the whole file (Ctrl+Shift+S in RStudio) after Steps 3, 3c and 3d.
# Needs: dplyr, tidyr, ggplot2, sf, maps, patchwork (ggrepel optional, for labels).
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "sf", "maps", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(ggplot2)
library(sf)

suppressMessages(sf_use_s2(FALSE))

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

site_index_file <- file.path("Outputs", "3d_resolved", "site_index.csv")
time_bins_file  <- file.path("Outputs", "3d_resolved", "time_bins.csv")
fauna_file      <- file.path("Outputs", "3d_resolved", "faunmap_fauna.csv")
pbdb_file       <- file.path("Outputs", "3d_resolved", "pbdb_occurrences.csv")
only_sites_with_species <- TRUE   # draw only sites with >= 1 occurrence identified to species
layers_file     <- file.path("Outputs", "3_spatial", "spatial_layers.rds")
spatial_file    <- file.path("Outputs", "3_spatial", "sites_spatial.csv")   # fallback for provinces
output_dir      <- file.path(work_dir, "Outputs", "4_maps")

temporal_level <- "province"     # units drawn under the temporal maps
spatiotemporal_levels <- c("province", "division", "political")

# Map frame (degrees). Hawaii (155-160 W) lies outside it.
map_crs <- "+proj=laea +lat_0=45 +lon_0=-98 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
map_lon <- c(-136, -52)
map_lat <- c(14, 61)

# Figure size (inches) and resolution.
map_width  <- 16
map_height <- 9
map_dpi    <- 300

# All sites are drawn alike, whatever their database (FAUNMAP or PBDB).
site_fill  <- "#B2182B"
point_size <- 2.2

unit_tints    <- c("#E9DDB9", "#CFE3C6", "#C8DCEC", "#DCD2EA", "#D2E9E2", "#F1D9C9")   # spatial maps
neutral_tints <- c("#F1ECE0", "#E3ECF1", "#EAE5F0", "#E6EFE3", "#F1E7E2", "#ECECEC")   # under coloured bins
empty_fill   <- "grey96"
no_data_fill <- "grey86"
sea_fill     <- "#EEF4F8"

level_titles <- c(political = "States and provinces",
                  division  = "Physiographic divisions (USA) and regions (Canada)",
                  province  = "Physiographic provinces (USA) and subregions (Canada)")
key_titles   <- c(political = "States and provinces", division = "Physiographic units",
                  province = "Provinces / subregions by division")
unit_cols    <- c(political = "State_Province", division = "Physio_Division", province = "Physio_Province")
level_short  <- c(political = "states and provinces", division = "physiographic divisions and regions",
                  province = "physiographic provinces and subregions")
st_titles    <- c(political = "Sites through time: states and provinces",
                  division  = "Sites through time: physiographic divisions",
                  province  = "Sites through time: physiographic provinces")

# =============================================================================
# HELPERS
# =============================================================================

read_text_csv <- function(file) {
  path <- file.path(work_dir, file)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun Steps 3, 3c and 3d first.")
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
                 na.strings = c("", "NA"), encoding = "UTF-8")
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))
  df[] <- lapply(df, fix_text)
  cat(sprintf("  %-45s %7d rows\n", file, nrow(df)))
  df
}

# Text written in Windows encoding (e.g. "Qu\xe9bec") is converted to UTF-8.
fix_text <- function(x) {
  bad <- !is.na(x) & !validUTF8(x)
  x[bad] <- iconv(x[bad], "latin1", "UTF-8")
  # "<U+00E9>" written by R in some locales -> the letter itself
  esc <- !is.na(x) & grepl("<U\\+[0-9A-Fa-f]{4,6}>", x)
  if (any(esc)) {
    x[esc] <- vapply(x[esc], function(s) {
      m <- gregexpr("<U\\+[0-9A-Fa-f]{4,6}>", s)[[1]]
      codes <- regmatches(s, list(m))[[1]]
      for (cd in unique(codes)) s <- gsub(cd, intToUtf8(strtoi(substr(cd, 4, nchar(cd) - 1), 16L)), s, fixed = TRUE)
      s
    }, character(1), USE.NAMES = FALSE)
  }
  Encoding(x) <- "UTF-8"
  x
}
# Names compared without accents or capitals ("Quebec" = "Québec").
fold <- function(x) tolower(iconv(x, "UTF-8", "ASCII//TRANSLIT"))

polygons_only <- function(x) {
  if (any(st_geometry_type(x) == "GEOMETRYCOLLECTION")) x <- st_collection_extract(x, "POLYGON", warn = FALSE)
  x <- x[st_geometry_type(x) %in% c("POLYGON", "MULTIPOLYGON"), ]
  st_cast(x, "MULTIPOLYGON", warn = FALSE)
}

# Repair invalid polygons (simplified outlines can cross themselves, which makes
# GEOS stop with "TopologyException: side location conflict").
repair_geom <- function(x) {
  x <- suppressWarnings(suppressMessages(polygons_only(st_make_valid(x))))
  bad <- !st_is_valid(x)
  if (any(bad, na.rm = TRUE)) {
    st_geometry(x)[which(bad)] <- suppressWarnings(st_buffer(st_geometry(x)[which(bad)], 0))
    x <- suppressWarnings(suppressMessages(polygons_only(st_make_valid(x))))
  }
  x[!st_is_empty(x), ]
}
# Run a geometry step; if GEOS fails, repair the input (zero buffer, then a
# 1-m precision grid) and try again.
geom_retry <- function(f, x) {
  tryCatch(f(x), error = function(e) {
    x2 <- repair_geom(suppressWarnings(st_buffer(repair_geom(x), 0)))
    tryCatch(f(x2), error = function(e2) {
      x3 <- repair_geom(st_set_precision(x2, 1))
      f(x3)
    })
  })
}

fmt_age <- function(x) ifelse(x < 0.1, sprintf("%.4f", x), sprintf("%.2f", x))   # 3.25, 0.0117

save_map <- function(p, dir, name, width = map_width, height = map_height) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  png_file <- file.path(dir, paste0(name, ".png"))
  ggsave(png_file, p, width = width, height = height, dpi = map_dpi, bg = "white")
  pdf_device <- if (capabilities("cairo")) cairo_pdf else pdf
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = width, height = height,
         device = pdf_device, bg = "white")
  cat("  ", png_file, " (+ .pdf)\n", sep = "")
}

# Numbered key drawn beside a map, grouped under bold headings.
key_panel <- function(key, title, subtitle = "Number on the map  -  name  (sites)") {
  lines <- bind_rows(lapply(unique(key$Group), function(g) {
    k <- key[key$Group == g, ]
    bind_rows(data.frame(type = "head", text = g, num = NA, n = NA),
              data.frame(type = "item", text = k$Unit, num = k$Number, n = k$n_sites))
  }))
  n_col <- ceiling(nrow(lines) / 40)
  per_col <- ceiling(nrow(lines) / n_col)
  lines$col <- pmin(n_col, ceiling(seq_len(nrow(lines)) / per_col))
  # a heading at the bottom of a column moves to the top of the next one
  for (i in which(lines$type == "head")) {
    if (i < nrow(lines) && lines$col[i + 1] > lines$col[i]) lines$col[i] <- lines$col[i + 1]
  }
  lines <- lines %>% group_by(col) %>% mutate(row = row_number()) %>% ungroup() %>%
    mutate(x = col - 1, y = -row)
  items <- filter(lines, type == "item")
  heads <- filter(lines, type == "head")
  ggplot() +
    geom_text(data = heads, aes(x = x, y = y, label = text), hjust = 0, fontface = "bold",
              size = 3.6, colour = "grey10") +
    geom_label(data = items, aes(x = x + 0.05, y = y, label = num,
                                 fill = ifelse(n > 0, "occupied", "empty")),
               size = 3, label.size = 0.2, label.padding = unit(0.12, "lines"),
               fontface = "bold", colour = "grey10") +
    geom_text(data = items, aes(x = x + 0.14, y = y,
                                label = ifelse(n > 0, sprintf("%s  (%d)", text, n), text),
                                colour = ifelse(n > 0, "occupied", "empty")),
              hjust = 0, size = 3.3) +
    scale_fill_manual(values = c(occupied = "white", empty = empty_fill), guide = "none") +
    scale_colour_manual(values = c(occupied = "grey10", empty = "grey55"), guide = "none") +
    scale_x_continuous(limits = c(-0.02, max(n_col, 2) + 0.15), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-(max(max(lines$row), 30) + 0.8), 0), expand = c(0, 0)) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10, margin = margin(b = 6)))
}

map_theme <- function() {
  theme_void(base_size = 13) +
    theme(panel.background = element_rect(fill = sea_fill, colour = "grey30", linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 17),
          plot.subtitle = element_text(colour = "grey30", size = 12, margin = margin(b = 6)),
          plot.caption = element_text(colour = "grey35", size = 9.5, hjust = 0),
          legend.position = "bottom", legend.text = element_text(size = 11),
          legend.title = element_text(size = 11, face = "bold"),
          legend.margin = margin(t = 4))
}

# =============================================================================
# 1. READ THE SITES, TIME BINS AND LAYERS
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
site_index <- read_text_csv(site_index_file)
time_bins  <- read_text_csv(time_bins_file)

# Keep sites that still have at least one occurrence identified to species.
if (only_sites_with_species) {
  fauna_occ <- read_text_csv(fauna_file)
  pbdb_occ  <- read_text_csv(pbdb_file)
  pb_name <- pbdb_occ[[intersect(c("accepted_name", "identified_name"), names(pbdb_occ))[1]]]
  is_species <- function(genus, species) {
    !is.na(genus) & trimws(genus) != "" & !is.na(species) &
      !grepl("^(sp|spp|indet)\\.?$", tolower(trimws(species))) & trimws(species) != ""
  }
  with_species <- unique(c(
    with(fauna_occ, paste(Site_Key, Stage_Number)[is_species(Genus, Species)]),
    with(pbdb_occ, paste(Site_Key, Stage_Number)[is_species(sub(" .*$", "", pb_name),
                                                            ifelse(grepl(" ", pb_name), sub("^\\S+\\s+", "", pb_name), NA))])))
  n_before <- nrow(site_index)
  site_index <- site_index[paste(site_index$Site_Key, site_index$Stage_Number) %in% with_species, ]
  cat(sprintf("  %d of %d site-bin records have an occurrence identified to species (the others are not drawn)\n",
              nrow(site_index), n_before))
}
layers_path <- file.path(work_dir, layers_file)
if (!file.exists(layers_path)) stop("File not found:\n  ", layers_path, "\nRun Step 3 (spatial binning) again.")
layers <- readRDS(layers_path)
for (nm in setdiff(names(layers), "settings")) {
  if (!is.null(layers[[nm]])) layers[[nm]] <- repair_geom(layers[[nm]])
}
work_crs <- layers$settings$work_crs
canada_north_limit <- layers$settings$canada_north_limit

# Older site_index.csv files have no spatial columns: take them from Step 3.
if (!all(unit_cols %in% names(site_index))) {
  sp <- read_text_csv(spatial_file)
  norm_key <- function(k) tolower(gsub("\\s+", " ", trimws(k)))
  sp <- sp[!duplicated(norm_key(sp$Site_Key)), ]
  m <- match(norm_key(site_index$Site_Key), norm_key(sp$Site_Key))
  for (col in c("Country", unit_cols)) site_index[[col]] <- sp[[col]][m]
  cat("  NOTE: spatial units taken from sites_spatial.csv (site_index.csv predates Step 3).\n")
}

# Spell unit names exactly as in the layers (protects against accents lost in a CSV).
for (lv in names(unit_cols)) {
  layer_units <- unique(layers[[lv]]$Unit)
  m <- match(fold(site_index[[unit_cols[[lv]]]]), fold(layer_units))
  site_index[[unit_cols[[lv]]]] <- ifelse(is.na(m), site_index[[unit_cols[[lv]]]], layer_units[m])
}

# Sites are localities: the analysis units of one FAUNMAP Machine Number (e.g.
# the ~800 levels of Vallecito Creek) count once per time bin, as in Steps 5-8.
site_index <- site_index %>%
  mutate(Site_Key_Unit = Site_Key,
         Site_Key = ifelse(Database == "FAUNMAP",
                           paste("FAUNMAP", vapply(strsplit(Site_Key, " | ", fixed = TRUE), `[`, "", 2), sep = " | "),
                           Site_Key)) %>%
  distinct(Site_Key, Stage_Number, .keep_all = TRUE)

time_bins <- time_bins %>%
  mutate(Bin_Number = as.integer(Bin_Number), Older_Ma = as.numeric(Older_Ma),
         Younger_Ma = as.numeric(Younger_Ma),
         Bin_Text = sprintf("Bin %d  (%s - %s Ma)", Bin_Number, fmt_age(Older_Ma), fmt_age(Younger_Ma))) %>%
  arrange(Bin_Number)
n_bins <- nrow(time_bins)
# Bin colours: old (dark purple) to young (orange), without the palest yellow.
bin_colours <- setNames(hcl.colors(n_bins + 1, "Plasma")[seq_len(n_bins)], time_bins$Bin_Text)

sites <- site_index %>%
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude),
         Stage_Number = as.integer(Stage_Number)) %>%
  filter(!is.na(Latitude), !is.na(Longitude), !is.na(Stage_Number)) %>%
  mutate(Bin_Text = factor(time_bins$Bin_Text[match(Stage_Number, time_bins$Bin_Number)],
                           levels = time_bins$Bin_Text))
cat(sprintf("  %d site-bin records with coordinates (%d distinct sites) in %d time bins\n",
            nrow(sites), n_distinct(sites$Site_Key), n_bins))

sites_sf <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(map_crs)

# =============================================================================
# 2. BASE LAYERS AND UNIT KEYS
# =============================================================================

cat("\n=== 2. PREPARING THE BASE MAP ===\n")
frame_ll <- st_as_sfc(st_bbox(c(xmin = map_lon[1], xmax = map_lon[2],
                                ymin = map_lat[1], ymax = map_lat[2]), crs = 4326))
countries <- st_as_sf(maps::map("world", regions = c("Canada", "USA", "Mexico"), fill = TRUE, plot = FALSE))
st_crs(countries) <- NA; st_crs(countries) <- 4326   # relabel the maps package's old datum
countries <- suppressWarnings(suppressMessages(
  countries %>% st_make_valid() %>% st_set_agr("constant") %>% st_intersection(frame_ll) %>%
    st_transform(map_crs)))
lims <- st_bbox(countries)

limit_line <- NULL
if (!is.na(canada_north_limit)) {
  limit_line <- st_sfc(st_linestring(cbind(seq(-141, -52, by = 0.5), canada_north_limit)), crs = 4326) %>%
    st_transform(map_crs)
}
study_area <- st_as_sfc(st_bbox(c(xmin = map_lon[1], xmax = map_lon[2], ymin = map_lat[1],
                                  ymax = min(map_lat[2], coalesce(canada_north_limit, 90))), crs = 4326))
study_area_work <- st_transform(study_area, work_crs)
study_area_map  <- st_set_agr(st_sf(geometry = st_transform(study_area, map_crs)), "constant")
country_order <- c("USA", "Canada", "Mexico")
use_repel <- requireNamespace("ggrepel", quietly = TRUE)

# Units of one level: merged by name across countries, numbered under headings
# (country, or parent division for provinces), with label positions.
prepare_level <- function(lv) {
  merge_units <- function(x, union = TRUE) {
    suppressWarnings(suppressMessages(
      x %>% group_by(Unit) %>%
        summarise(Countries = paste(intersect(country_order, unique(Country)), collapse = ", "),
                  .groups = "drop", do_union = union) %>%
        polygons_only()))
  }
  # If dissolving still fails, the pieces of a unit are kept side by side (same look).
  polys <- tryCatch(geom_retry(merge_units, layers[[lv]]),
                    error = function(e) merge_units(layers[[lv]], union = FALSE))
  polys <- repair_geom(polys)
  polys <- polys[lengths(st_intersects(polys, study_area_work)) > 0 | polys$Unit %in% sites[[unit_cols[[lv]]]], ]

  if (lv == "province") {
    parent_layer <- rbind(layers$us_division, layers$ca_division)
    pts <- suppressWarnings(suppressMessages(
      polys %>% st_set_agr("constant") %>% st_cast("POLYGON", warn = FALSE) %>%
        mutate(area = as.numeric(st_area(.))) %>% group_by(Unit) %>%
        slice_max(area, n = 1, with_ties = FALSE) %>% ungroup() %>%
        st_set_agr("constant")))
    pts <- tryCatch(suppressWarnings(st_point_on_surface(pts)),
                    error = function(e) suppressWarnings(st_centroid(pts)))
    near <- st_nearest_feature(pts, parent_layer)
    polys$Group <- parent_layer$Unit[near][match(polys$Unit, pts$Unit)]
    first_country <- tapply(match(sub(",.*", "", polys$Countries), country_order), polys$Group, min)
    group_levels <- names(first_country)[order(first_country, names(first_country))]
  } else {
    polys$Group <- sub(",.*", "", polys$Countries)
    group_levels <- intersect(country_order, polys$Group)
  }
  key <- st_drop_geometry(polys) %>%
    mutate(Group = factor(Group, levels = group_levels)) %>%
    arrange(Group, Unit) %>%
    mutate(Number = row_number(), Group = as.character(Group))

  draw <- polys %>% st_transform(map_crs) %>% left_join(select(key, Unit, Number), by = "Unit")
  touching <- st_is_within_distance(draw, draw, dist = 20000)
  tint <- integer(nrow(draw))
  for (i in order(lengths(touching), decreasing = TRUE)) {
    used <- tint[setdiff(touching[[i]], i)]
    tint[i] <- which(!seq_along(unit_tints) %in% used)[1]
    if (is.na(tint[i])) tint[i] <- 1L
  }
  draw$Tint <- tint

  draw <- repair_geom(draw)
  # Label inside the largest piece of each unit within the study area.
  clip <- function(x) suppressWarnings(suppressMessages(
    x %>% select(Number) %>% st_set_agr("constant") %>% st_intersection(study_area_map)))
  pieces <- tryCatch(geom_retry(clip, draw), error = function(e) select(draw, Number))
  pieces <- suppressWarnings(suppressMessages(
    pieces %>% st_set_agr("constant") %>% st_cast("POLYGON", warn = FALSE) %>%
      mutate(area = as.numeric(st_area(.))) %>%
      group_by(Number) %>% slice_max(area, n = 1, with_ties = FALSE) %>% ungroup() %>%
      st_set_agr("constant")))
  label_pts <- tryCatch(suppressWarnings(st_point_on_surface(pieces)),
                        error = function(e) suppressWarnings(st_centroid(pieces)))
  label_xy <- cbind(st_drop_geometry(label_pts)["Number"], st_coordinates(label_pts))
  cat(sprintf("  %-10s %3d units\n", lv, nrow(key)))
  list(level = lv, draw = draw, key = key, label_xy = label_xy)
}

needed_levels <- unique(c("political", "division", "province", temporal_level, spatiotemporal_levels))
preps <- setNames(lapply(needed_levels, prepare_level), needed_levels)

# One map: units (tinted where they hold sites), numbered labels, points, key.
make_map <- function(prep, pts, point_layers, title, subtitle, caption, key_title,
                     tints = unit_tints, key_subtitle = "Number on the map  -  name  (sites)") {
  ucol <- unit_cols[[prep$level]]
  counts <- st_drop_geometry(pts) %>% filter(!is.na(.data[[ucol]])) %>%
    distinct(Site_Key, Unit = .data[[ucol]]) %>% count(Unit, name = "n_sites")
  key <- prep$key %>% left_join(counts, by = "Unit") %>% mutate(n_sites = coalesce(n_sites, 0L))
  draw <- prep$draw %>% left_join(select(key, Unit, n_sites), by = "Unit")
  draw$Fill <- ifelse(draw$n_sites > 0, tints[draw$Tint], empty_fill)
  lab <- prep$label_xy %>% left_join(select(key, Number, n_sites), by = "Number")

  p <- ggplot() +
    geom_sf(data = countries, fill = no_data_fill, colour = NA) +
    geom_sf(data = draw, fill = draw$Fill, colour = "grey45", linewidth = 0.25) +
    geom_sf(data = countries, fill = NA, colour = "grey15", linewidth = 0.45)
  if (!is.null(limit_line)) p <- p + geom_sf(data = limit_line, colour = "grey25", linewidth = 0.4, linetype = "22")
  p <- p + point_layers
  if (use_repel) {
    p <- p + ggrepel::geom_label_repel(
      data = lab, aes(x = X, y = Y, label = Number), size = 3, fontface = "bold",
      colour = ifelse(lab$n_sites > 0, "grey5", "grey50"), fill = alpha("white", 0.85),
      label.size = 0.15, label.padding = unit(0.1, "lines"), min.segment.length = 0.2,
      segment.colour = "grey25", segment.size = 0.25, max.overlaps = Inf, seed = 1, box.padding = 0.15)
  } else {
    p <- p + geom_label(data = lab, aes(x = X, y = Y, label = Number), size = 3, fontface = "bold",
                        fill = alpha("white", 0.85), label.size = 0.15)
  }
  p <- p +
    coord_sf(crs = map_crs, xlim = lims[c("xmin", "xmax")], ylim = lims[c("ymin", "ymax")], expand = FALSE) +
    labs(title = title, subtitle = subtitle, caption = caption, x = NULL, y = NULL) +
    map_theme()
  patchwork::wrap_plots(p, key_panel(key, key_title, key_subtitle), widths = c(1.45, 1)) &
    theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(10, 12, 8, 12))
}

base_caption <- function(lv) {
  paste0("Sites in Alaska, Hawaii",
         if (!is.na(canada_north_limit)) paste0(" and Canada north of ", canada_north_limit, "\u00b0N (dashed line)"),
         " removed. Pale units hold no sites.",
         if (lv != "political") " Mexico: no physiographic layer (grey)." else "")
}

# Site points for the spatial and temporal maps: one colour for every site.
site_points <- function(pts, size = point_size) {
  list(geom_sf(data = pts, shape = 21, fill = site_fill, colour = "white", stroke = 0.3,
               size = size, alpha = 0.9))
}

# =============================================================================
# 3. SPATIAL MAPS (ALL TIME BINS TOGETHER, ONE MAP PER LEVEL)
# =============================================================================

cat("\n=== 3. SPATIAL MAPS ===\n")
spatial_pts <- sites_sf %>% distinct(Site_Key, .keep_all = TRUE)
for (lv in c("political", "division", "province")) {
  n_in <- n_distinct(spatial_pts$Site_Key[!is.na(spatial_pts[[unit_cols[[lv]]]])])
  p <- make_map(preps[[lv]], spatial_pts, site_points(spatial_pts),
                title = level_titles[[lv]],
                subtitle = sprintf("%d fossil sites, all time bins (%s - %s Ma)",
                                   n_in, fmt_age(max(time_bins$Older_Ma)), fmt_age(min(time_bins$Younger_Ma))),
                caption = base_caption(lv), key_title = key_titles[[lv]])
  save_map(p, file.path(output_dir, "1_spatial"), paste0("map_", lv))
}

# =============================================================================
# 4. TEMPORAL MAPS (ONE PER TIME BIN, AND ALL BINS SIDE BY SIDE)
# =============================================================================

cat("\n=== 4. TEMPORAL MAPS ===\n")
# remove per-bin maps of an earlier run (e.g. a bin that no longer exists)
unlink(list.files(file.path(output_dir, "2_temporal"), "^map_bin[0-9]+\\.(png|pdf)$", full.names = TRUE))
for (b in seq_len(n_bins)) {
  pts <- filter(sites_sf, Stage_Number == time_bins$Bin_Number[b])
  p <- make_map(preps[[temporal_level]], pts, site_points(pts),
                title = time_bins$Bin_Text[b],
                subtitle = sprintf("%d sites on the %s", n_distinct(pts$Site_Key),
                                   level_short[[temporal_level]]),
                caption = base_caption(temporal_level), key_title = key_titles[[temporal_level]],
                key_subtitle = "Number on the map  -  name  (sites in this bin)")
  save_map(p, file.path(output_dir, "2_temporal"), paste0("map_bin", time_bins$Bin_Number[b]))
}

# Small multiples: every bin on the same base map.
tp <- preps[[temporal_level]]
bin_n <- st_drop_geometry(sites_sf) %>% distinct(Bin_Text, Site_Key) %>% count(Bin_Text)
facet_labels <- setNames(sprintf("%s:  %d sites", time_bins$Bin_Text,
                                 coalesce(bin_n$n[match(time_bins$Bin_Text, bin_n$Bin_Text)], 0L)),
                         time_bins$Bin_Text)
p_all <- ggplot() +
  geom_sf(data = countries, fill = no_data_fill, colour = NA) +
  # one layer per tint (a per-unit fill vector does not work across facets)
  lapply(sort(unique(tp$draw$Tint)), function(t)
    geom_sf(data = tp$draw[tp$draw$Tint == t, ], fill = neutral_tints[t], colour = "grey55", linewidth = 0.15)) +
  geom_sf(data = countries, fill = NA, colour = "grey20", linewidth = 0.3)
if (!is.null(limit_line)) p_all <- p_all + geom_sf(data = limit_line, colour = "grey30", linewidth = 0.3, linetype = "22")
p_all <- p_all + site_points(sites_sf %>% distinct(Bin_Text, Site_Key, .keep_all = TRUE), size = 1.4) +
  facet_wrap(~ Bin_Text, ncol = ceiling(n_bins / 2), labeller = as_labeller(facet_labels)) +
  coord_sf(crs = map_crs, xlim = lims[c("xmin", "xmax")], ylim = lims[c("ymin", "ymax")], expand = FALSE) +
  labs(title = "Sites in each time bin",
       subtitle = paste("Fossil sites on the", level_short[[temporal_level]]),
       caption = sub(" Pale units hold no sites.", "", base_caption(temporal_level), fixed = TRUE)) +
  map_theme() +
  theme(strip.text = element_text(face = "bold", size = 12, hjust = 0, margin = margin(b = 4)),
        panel.spacing = unit(0.8, "lines"), plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 12, 8, 12))
save_map(p_all, file.path(output_dir, "2_temporal"), "temporal_all_bins")

# =============================================================================
# 5. SPATIO-TEMPORAL MAPS (ALL BINS, COLOURED BY BIN) AND SITES PER UNIT x BIN
# =============================================================================

cat("\n=== 5. SPATIO-TEMPORAL MAPS ===\n")
st_pts <- sites_sf %>% distinct(Site_Key, Stage_Number, .keep_all = TRUE) %>%
  arrange(desc(Stage_Number))          # youngest drawn first, oldest on top (fewer, easier to miss)
bin_counts <- st_drop_geometry(st_pts) %>% count(Bin_Text)
bin_labels <- setNames(sprintf("%s: %d", names(bin_colours),
                               coalesce(bin_counts$n[match(names(bin_colours), bin_counts$Bin_Text)], 0L)),
                       names(bin_colours))
bin_points <- list(
  geom_sf(data = st_pts, aes(fill = Bin_Text), shape = 21, colour = "grey15", stroke = 0.25,
          size = point_size, alpha = 0.9),
  scale_fill_manual(values = bin_colours, labels = bin_labels, name = "Time bin (sites)", drop = FALSE),
  guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 3, byrow = TRUE,
                             title.position = "top")))

for (lv in spatiotemporal_levels) {
  p <- make_map(preps[[lv]], st_pts, bin_points,
                title = st_titles[[lv]],
                subtitle = sprintf("%d sites in %d time bins of %s Myr (the last bin is shorter)",
                                   n_distinct(st_pts$Site_Key), n_bins,
                                   formatC(median(time_bins$Older_Ma - time_bins$Younger_Ma), format = "fg")),
                caption = base_caption(lv), key_title = key_titles[[lv]], tints = neutral_tints)
  save_map(p, file.path(output_dir, "3_spatiotemporal"), paste0("map_spatiotemporal_", lv))

  # Sites per unit and bin (table for every level; heat map for provinces).
  tab <- st_drop_geometry(st_pts) %>% filter(!is.na(.data[[unit_cols[[lv]]]])) %>%
    distinct(Unit = .data[[unit_cols[[lv]]]], Bin_Text, Site_Key) %>% count(Unit, Bin_Text) %>%
    left_join(select(preps[[lv]]$key, Unit, Number, Group), by = "Unit")
  wide <- tab %>% mutate(Bin_Text = paste0("Bin_", match(Bin_Text, time_bins$Bin_Text))) %>%
    tidyr::pivot_wider(names_from = Bin_Text, values_from = n, values_fill = 0) %>% arrange(Number)
  names(wide)[names(wide) == "Unit"] <- unit_cols[[lv]]
  write.csv(wide, file.path(output_dir, "3_spatiotemporal", paste0("sites_per_unit_and_bin_", lv, ".csv")),
            row.names = FALSE)

  if (lv == temporal_level) {
    full <- tidyr::expand_grid(Unit = unique(tab$Unit), Bin_Text = time_bins$Bin_Text) %>%
      left_join(select(tab, Unit, Bin_Text, n), by = c("Unit", "Bin_Text")) %>%
      left_join(select(preps[[lv]]$key, Unit, Number, Group), by = "Unit") %>%
      mutate(n = coalesce(n, 0L),
             Row = sprintf("%d  %s", Number, Unit),
             Bin_Text = factor(Bin_Text, levels = time_bins$Bin_Text),
             Group = factor(Group, levels = unique(preps[[lv]]$key$Group)),
             txt = ifelse(sqrt(n) >= 0.6 * sqrt(max(n)), "white", "grey10")) %>%
      mutate(Row = factor(Row, levels = rev(unique(Row[order(Number)]))))
    p_heat <- ggplot(full, aes(Bin_Text, Row, fill = n)) +
      geom_tile(colour = "white", linewidth = 1.2) +
      geom_text(aes(label = ifelse(n > 0, n, ""), colour = txt), size = 4, fontface = "bold") +
      scale_colour_identity() +
      scale_fill_gradientn(colours = c("#F4F4F2", "#CDE2FB", "#6DA7EC", "#256ABF", "#0D366B"),
                           name = "Sites", limits = c(0, NA), trans = "sqrt") +
      scale_x_discrete(position = "top",
                       labels = function(x) sub("^(Bin [0-9]+)  \\((.*)\\)$", "\\1\n\\2", x)) +
      facet_grid(Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
      labs(title = "Sites per province and time bin",
           subtitle = "Provinces (USA) and subregions (Canada) holding sites, grouped by division; numbers as on the maps",
           x = NULL, y = NULL) +
      theme_minimal(base_size = 14) +
      theme(panel.grid = element_blank(), axis.text.x.top = element_text(face = "bold", size = 9.5),
            axis.text.y = element_text(size = 11, colour = "grey10"),
            strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1, size = 11),
            strip.placement = "outside", panel.spacing = unit(0.4, "lines"),
            plot.title = element_text(face = "bold", size = 17), plot.title.position = "plot",
            plot.subtitle = element_text(colour = "grey30", size = 12),
            plot.background = element_rect(fill = "white", colour = NA), legend.position = "right")
    save_map(p_heat, file.path(output_dir, "3_spatiotemporal"), "sites_per_unit_and_bin",
             width = 13.33, height = max(5, 0.32 * n_distinct(full$Unit) + 2.5))
  }
}

cat("\n=== STEP 4 COMPLETE ===\n")
cat("Maps in:", output_dir, "\n")
