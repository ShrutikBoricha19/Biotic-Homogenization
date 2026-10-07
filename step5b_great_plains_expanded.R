# =============================================================================
# STEP 5b: CENTRAL LOWLAND WEST OF THE MISSISSIPPI JOINS THE GREAT PLAINS
#
# To strengthen the Great Plains sample, every site of the Central Lowland
# province that lies WEST of the Mississippi River - and all its faunal
# occurrences - is counted as a Great Plains site from here on, in every time
# bin. Earlier steps are not changed: this step reads the Step 5 master data
# and writes new copies with the updated spatial bin.
#
# West of the river: a site is west when its longitude is smaller than the
# Mississippi's longitude at the site's latitude (river line from Natural
# Earth 1:10m, downloaded once). North of the river's source (Lake Itasca,
# ~47.5 N) the source longitude is used. Sites without coordinates keep their
# province.
#
# Outputs (Outputs/5b_great_plains/):
#   master_data.csv, master_data_unique.csv   as Step 5, with Spatial_Bin
#       updated; the Step 5 province is kept in Spatial_Bin_Original and
#       Moved_To_Great_Plains marks the reassigned records
#   master_data_summary.csv     species, sites, occurrences per bin and province
#   reassigned_sites.csv        the Central Lowland sites that were moved
#   great_plains_sites_per_bin.csv   Great Plains sites per time bin: original,
#                                    added, total
#   map_great_plains_expanded.png / .pdf   spatio-temporal map with the river
#
# Inputs: Outputs/5_master_data/ (Step 5), Outputs/3d_resolved/time_bins.csv,
#         Outputs/3_spatial/spatial_layers.rds (Step 3, for the map).
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 5.
# Needs: dplyr, ggplot2, sf, maps, patchwork.
# =============================================================================

for (pkg in c("dplyr", "ggplot2", "sf", "maps", "patchwork")) {
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

master_dir   <- file.path("Outputs", "5_master_data")
bins_file    <- file.path("Outputs", "3d_resolved", "time_bins.csv")
layers_file  <- file.path("Outputs", "3_spatial", "spatial_layers.rds")
basemap_dir  <- file.path(work_dir, "Outputs", "maps", "basemap")     # download cache
output_dir   <- file.path(work_dir, "Outputs", "5b_great_plains")

source_province <- "Central Lowland"   # province whose western part is moved
target_province <- "Great Plains"      # province it joins

rivers_url <- paste0("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/",
                     "geojson/ne_10m_rivers_lake_centerlines.geojson")
rivers_local <- NA                     # or the path of a local copy of that file
river_name   <- "Mississippi"

# Map
map_crs <- "+proj=laea +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
map_lon <- c(-115, -79)               # frame of the map (degrees)
map_lat <- c(27, 51)
map_width <- 16; map_height <- 9; map_dpi <- 300
region_fills <- c("Great Plains" = "#CFE3C6",
                  "Central Lowland west of the Mississippi (added)" = "#9FD08F",
                  "Central Lowland east of the Mississippi (unchanged)" = "#E8E4D8")
other_fill   <- "#F4F2EC"
no_data_fill <- "grey88"
sea_fill     <- "#EEF4F8"
river_colour <- "#1F5FAD"

# =============================================================================
# HELPERS
# =============================================================================

read_text_csv <- function(path) {
  full <- file.path(work_dir, path)
  if (!file.exists(full)) stop("File not found:\n  ", full, "\nRun the earlier steps first.")
  df <- read.csv(full, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
                 na.strings = c("", "NA"), encoding = "UTF-8")
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))
  cat(sprintf("  %-45s %7d rows\n", path, nrow(df)))
  df
}
save_csv <- function(df, name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, paste0(name, ".csv"))
  ok <- tryCatch({ write.csv(df, path, row.names = FALSE, na = "", fileEncoding = "UTF-8"); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) stop("Could not write ", path, "\n  Close it (e.g. in Excel) and run again.")
  cat(sprintf("  %-30s %7d rows -> %s\n", name, nrow(df), path))
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
fmt_age <- function(x) ifelse(x < 0.1, sprintf("%.4f", x), sprintf("%.2f", x))

# =============================================================================
# 1. READ THE MASTER DATA
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
master    <- read_text_csv(file.path(master_dir, "master_data.csv"))
unique_md <- read_text_csv(file.path(master_dir, "master_data_unique.csv"))
time_bins <- read_text_csv(bins_file)

provinces_found <- sort(unique(na.omit(master$Spatial_Bin)))
for (p in c(source_province, target_province)) {
  if (!p %in% provinces_found) {
    stop("No sites in province '", p, "'. Provinces in the master data:\n  ",
         paste(provinces_found, collapse = "; "))
  }
}

# =============================================================================
# 2. THE MISSISSIPPI AND WHICH SITES LIE WEST OF IT
# =============================================================================

cat("\n=== 2. MISSISSIPPI RIVER ===\n")
dir.create(basemap_dir, recursive = TRUE, showWarnings = FALSE)
rivers_file <- if (!is.na(rivers_local)) rivers_local else file.path(basemap_dir, basename(rivers_url))
if (!file.exists(rivers_file)) {
  cat("  Downloading Natural Earth rivers (one time, ~7 MB)...\n")
  ok <- tryCatch({ download.file(rivers_url, rivers_file, mode = "wb", quiet = TRUE); TRUE },
                 error = function(e) FALSE)
  if (!ok) {
    if (file.exists(rivers_file)) file.remove(rivers_file)
    stop("Download failed. Download ne_10m_rivers_lake_centerlines.geojson from\n  ", rivers_url,
         "\nand set 'rivers_local' to its path.")
  }
}
rivers <- st_read(rivers_file, quiet = TRUE)
river <- rivers[rivers$name %in% river_name, ]
if (nrow(river) == 0) stop("No river called '", river_name, "' in ", rivers_file)
# Add points every 1 km along the river (in metres, then back to degrees).
aea <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +datum=WGS84 +units=m +no_defs"
river <- suppressWarnings(st_transform(st_segmentize(st_transform(river, aea), dfMaxLength = 1000), 4326))
rxy <- st_coordinates(river)[, c("X", "Y")]
cat(sprintf("  %s: %.2f to %.2f degrees N\n", river_name, min(rxy[, 2]), max(rxy[, 2])))

# River longitude at any latitude: the westernmost river point in each
# 0.05-degree band of latitude, interpolated between bands; beyond the ends
# of the river the end longitudes are used.
bands <- seq(floor(min(rxy[, 2]) * 20) / 20, ceiling(max(rxy[, 2]) * 20) / 20, by = 0.05)
band_lon <- vapply(bands, function(b) {
  inb <- abs(rxy[, 2] - b) <= 0.05
  if (any(inb)) min(rxy[inb, 1]) else NA_real_
}, numeric(1))
ok <- !is.na(band_lon)
river_lon_at <- approxfun(bands[ok], band_lon[ok], rule = 2)

prep <- function(df) {
  lat <- suppressWarnings(as.numeric(df$Latitude)); lon <- suppressWarnings(as.numeric(df$Longitude))
  west <- !is.na(lat) & !is.na(lon) & lon < river_lon_at(lat)
  move <- df$Spatial_Bin %in% source_province & west
  df %>% mutate(Spatial_Bin_Original = Spatial_Bin,
                West_Of_Mississippi = ifelse(is.na(lat) | is.na(lon), NA, west),
                Moved_To_Great_Plains = move,
                Spatial_Bin = ifelse(move, target_province, Spatial_Bin)) %>%
    relocate(Spatial_Bin_Original, Moved_To_Great_Plains, .after = Spatial_Bin)
}
master_out <- prep(master)
unique_out <- prep(unique_md)
# Sites are localities (Step 5): count Locality_Key where it exists.
if ("Locality_Key" %in% names(master_out)) master_out <- master_out %>% mutate(Site_Key_Unit = Site_Key, Site_Key = Locality_Key)

cl <- master_out %>% filter(Spatial_Bin_Original %in% source_province)
cat(sprintf("  %s: %d sites (%d occurrences); west of the river: %d sites (%d occurrences) -> %s\n",
            source_province, n_distinct(cl$Site_Key), nrow(cl),
            n_distinct(cl$Site_Key[cl$Moved_To_Great_Plains]), sum(cl$Moved_To_Great_Plains), target_province))
if (any(is.na(cl$West_Of_Mississippi))) {
  cat(sprintf("  NOTE: %d %s occurrences have no coordinates and stay in %s.\n",
              sum(is.na(cl$West_Of_Mississippi)), source_province, source_province))
}

# =============================================================================
# 3. TABLES
# =============================================================================

cat("\n=== 3. SAVING TABLES ===\n")
bins <- time_bins %>% transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma),
                                Younger = as.numeric(Younger_Ma)) %>% arrange(Time_Bin)

reassigned_sites <- master_out %>% filter(Moved_To_Great_Plains) %>%
  group_by(Site_Key, Time_Bin) %>%
  summarise(Site = first(Site), Latitude = first(Latitude), Longitude = first(Longitude),
            State_Province = first(State_Province), n_occurrences = n(),
            n_species = n_distinct(Species), .groups = "drop") %>%
  arrange(as.integer(Time_Bin), Site)

gp_per_bin <- master_out %>% filter(Spatial_Bin %in% target_province) %>%
  mutate(Time_Bin = as.integer(Time_Bin)) %>%
  group_by(Time_Bin) %>%
  summarise(sites_original = n_distinct(Site_Key[!Moved_To_Great_Plains]),
            sites_added = n_distinct(Site_Key[Moved_To_Great_Plains]),
            sites_total = n_distinct(Site_Key),
            species_original = n_distinct(Species[!Moved_To_Great_Plains]),
            species_total = n_distinct(Species), .groups = "drop") %>%
  right_join(select(bins, Time_Bin), by = "Time_Bin") %>%
  mutate(across(-Time_Bin, ~ coalesce(.x, 0L))) %>% arrange(Time_Bin)

summary_out <- master_out %>%
  group_by(Time_Bin, Time_Bin_Label, Spatial_Bin) %>%
  summarise(n_species = n_distinct(Species), n_sites = n_distinct(Site_Key),
            n_occurrences = n(), .groups = "drop") %>%
  arrange(as.integer(Time_Bin), desc(n_sites))

save_csv(master_out, "master_data")
save_csv(unique_out, "master_data_unique")
save_csv(summary_out, "master_data_summary")
save_csv(reassigned_sites, "reassigned_sites")
save_csv(gp_per_bin, "great_plains_sites_per_bin")
cat("\n  Great Plains sites per time bin:\n")
print(as.data.frame(gp_per_bin), row.names = FALSE)

# =============================================================================
# 4. MAP: THE EXPANDED GREAT PLAINS WITH THE MISSISSIPPI
# =============================================================================

cat("\n=== 4. MAP ===\n")
layers_path <- file.path(work_dir, layers_file)
if (!file.exists(layers_path)) stop("File not found:\n  ", layers_path, "\nRun Step 3 again.")
layers <- readRDS(layers_path)
prov <- repair_geom(layers$province) %>% st_transform(4326)

# Area west of the river, as a polygon (river line extended north and south).
lat_seq <- seq(map_lat[1] - 3, map_lat[2] + 3, by = 0.05)
west_poly <- st_sfc(st_polygon(list(rbind(
  cbind(river_lon_at(lat_seq), lat_seq),
  c(-150, max(lat_seq)), c(-150, min(lat_seq)),
  c(river_lon_at(min(lat_seq)), min(lat_seq))))), crs = 4326)

merge_unit <- function(u) {
  x <- prov[prov$Unit %in% u, ]
  if (!nrow(x)) return(NULL)
  suppressWarnings(suppressMessages(st_union(st_make_valid(x))))
}
gp_geom <- merge_unit(target_province)
cl_geom <- merge_unit(source_province)
safe_op <- function(f, a, b) tryCatch(suppressWarnings(suppressMessages(f(a, b))),
                                      error = function(e) suppressWarnings(suppressMessages(
                                        f(st_buffer(st_make_valid(a), 0), st_buffer(st_make_valid(b), 0)))))
regions <- st_sf(
  Region = names(region_fills),
  geometry = c(gp_geom, safe_op(st_intersection, cl_geom, west_poly), safe_op(st_difference, cl_geom, west_poly)),
  crs = 4326) %>%
  mutate(Region = factor(Region, levels = names(region_fills)))

frame_ll <- st_as_sfc(st_bbox(c(xmin = map_lon[1], xmax = map_lon[2], ymin = map_lat[1], ymax = map_lat[2]),
                              crs = 4326))
countries <- st_as_sf(maps::map("world", regions = c("Canada", "USA", "Mexico"), fill = TRUE, plot = FALSE))
st_crs(countries) <- NA; st_crs(countries) <- 4326
clip_frame <- function(x) suppressWarnings(suppressMessages(
  st_intersection(st_set_agr(st_make_valid(x), "constant"), st_buffer(frame_ll, 2))))
countries <- clip_frame(countries) %>% st_transform(map_crs)
prov_m    <- clip_frame(prov) %>% st_transform(map_crs)
regions_m <- regions %>% st_transform(map_crs)
river_m   <- suppressWarnings(st_transform(river, map_crs))
lims <- st_bbox(st_transform(st_segmentize(st_transform(frame_ll, aea), 20000), map_crs))

# Sites: one point per site and bin; Great Plains sites coloured by bin.
bin_text <- sprintf("Bin %d  (%s - %s Ma)", bins$Time_Bin, fmt_age(bins$Older), fmt_age(bins$Younger))
bin_colours <- setNames(hcl.colors(nrow(bins) + 1, "Plasma")[seq_len(nrow(bins))], bin_text)
pts <- master_out %>%
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude), Time_Bin = as.integer(Time_Bin)) %>%
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  distinct(Site_Key, Time_Bin, .keep_all = TRUE) %>%
  mutate(In_GP = Spatial_Bin %in% target_province,
         Bin_Text = factor(bin_text[match(Time_Bin, bins$Time_Bin)], levels = bin_text)) %>%
  arrange(In_GP, desc(Time_Bin)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% st_transform(map_crs)
gp_pts <- pts[pts$In_GP, ]
bin_n <- table(gp_pts$Bin_Text)
bin_labels <- setNames(sprintf("%s: %d", names(bin_colours), as.integer(bin_n[names(bin_colours)])), names(bin_colours))

# Labels for the river and the two provinces.
lab_pt <- function(g) {
  p <- tryCatch(st_point_on_surface(st_transform(g, map_crs)), error = function(e) st_centroid(st_transform(g, map_crs)))
  st_coordinates(p)[1, ]
}
labs_df <- rbind(
  data.frame(t(lab_pt(gp_geom)), lab = "GREAT PLAINS"),
  data.frame(t(lab_pt(regions$geometry[3])), lab = "CENTRAL LOWLAND\n(east, unchanged)"))
river_lab <- st_coordinates(st_transform(st_sfc(st_point(c(river_lon_at(36.5) + 0.35, 36.5)), crs = 4326), map_crs))

p_map <- ggplot() +
  geom_sf(data = countries, fill = no_data_fill, colour = NA) +
  geom_sf(data = prov_m, fill = other_fill, colour = "grey60", linewidth = 0.2) +
  geom_sf(data = regions_m, aes(fill = Region), colour = "grey35", linewidth = 0.35) +
  geom_sf(data = countries, fill = NA, colour = "grey20", linewidth = 0.4) +
  geom_sf(data = river_m, aes(linetype = "Mississippi River"), colour = river_colour, linewidth = 1.1) +
  geom_sf(data = pts[!pts$In_GP, ], colour = "grey55", size = 0.9, alpha = 0.6) +
  geom_sf(data = gp_pts, colour = "white", size = 2.8) +                       # white halo
  geom_sf(data = gp_pts, aes(colour = Bin_Text), size = 2.1) +
  geom_label(data = labs_df, aes(X, Y, label = lab), fontface = "bold", size = 4, colour = "grey20",
             lineheight = 0.9, fill = alpha("white", 0.75), label.size = 0, label.padding = unit(0.15, "lines")) +
  annotate("text", x = river_lab[1], y = river_lab[2], label = "Mississippi R.", colour = river_colour,
           fontface = "bold.italic", size = 4, hjust = 0, angle = 0) +
  scale_fill_manual(values = region_fills, name = NULL, drop = FALSE) +
  scale_colour_manual(values = bin_colours, labels = bin_labels, name = "Great Plains sites per time bin",
                      drop = FALSE) +
  scale_linetype_manual(values = c("Mississippi River" = "solid"), name = NULL) +
  guides(fill = guide_legend(order = 1, ncol = 1),
         linetype = guide_legend(order = 2, override.aes = list(colour = river_colour, linewidth = 1.1)),
         colour = guide_legend(order = 3, ncol = 1, title.position = "top",
                               override.aes = list(size = 4))) +
  coord_sf(crs = map_crs, xlim = lims[c("xmin", "xmax")], ylim = lims[c("ymin", "ymax")], expand = FALSE) +
  labs(title = "Great Plains + Central Lowland west of the Mississippi",
       subtitle = sprintf("%d Great Plains sites in all time bins, %d of them from the western Central Lowland",
                          n_distinct(gp_pts$Site_Key), n_distinct(gp_pts$Site_Key[gp_pts$Moved_To_Great_Plains])),
       caption = "Grey points: sites in other provinces. Province boundaries: Fenneman & Johnson (1946); river: Natural Earth.") +
  theme_void(base_size = 13) +
  theme(panel.background = element_rect(fill = sea_fill, colour = "grey30", linewidth = 0.4),
        plot.title = element_text(face = "bold", size = 17),
        plot.subtitle = element_text(colour = "grey30", size = 12, margin = margin(b = 6)),
        plot.caption = element_text(colour = "grey35", size = 9.5, hjust = 0),
        legend.position = "right", legend.box = "vertical", legend.text = element_text(size = 10.5),
        legend.box.just = "left", legend.spacing.y = unit(0.4, "lines"),
        legend.title = element_text(size = 11, face = "bold"))

# Side panel: Great Plains sites per bin, before and after.
tab <- gp_per_bin %>% mutate(Bin = bin_text[match(Time_Bin, bins$Time_Bin)])
n_r <- nrow(tab)
tab_long <- rbind(
  data.frame(row = 0, col = 1:4, txt = c("Time bin", "Great\nPlains", "Added", "Total"), bold = TRUE),
  do.call(rbind, lapply(seq_len(n_r), function(i) data.frame(
    row = i, col = 1:4, bold = c(FALSE, FALSE, FALSE, TRUE),
    txt = c(tab$Bin[i], tab$sites_original[i], paste0("+", tab$sites_added[i]), tab$sites_total[i])))),
  data.frame(row = n_r + 1, col = 1:4, bold = TRUE,
             txt = c("All bins", n_distinct(gp_pts$Site_Key[!gp_pts$Moved_To_Great_Plains]),
                     paste0("+", n_distinct(gp_pts$Site_Key[gp_pts$Moved_To_Great_Plains])),
                     n_distinct(gp_pts$Site_Key))))
x_pos <- c(0, 2.9, 3.75, 4.5)
p_tab <- ggplot(tab_long) +
  annotate("rect", xmin = -0.1, xmax = 4.9, ymin = -(seq_len(n_r)) - 0.5, ymax = -(seq_len(n_r)) + 0.5,
           fill = rep(c("#F6F6F4", "white"), length.out = n_r)) +
  annotate("segment", x = -0.1, xend = 4.9, y = -0.5, yend = -0.5, colour = "grey40") +
  annotate("segment", x = -0.1, xend = 4.9, y = -(n_r + 0.5), yend = -(n_r + 0.5), colour = "grey40") +
  geom_text(aes(x = x_pos[col], y = -row, label = txt, fontface = ifelse(bold, "bold", "plain"),
                hjust = ifelse(col == 1, 0, 0.5),
                colour = ifelse(col == 3 & row > 0, "#2E7D32", "grey10")),
            size = 3.6, lineheight = 0.9) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(-0.15, 4.95), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-(n_r + 1.7), 0.9), expand = c(0, 0)) +
  labs(title = "Great Plains sites per time bin",
       subtitle = "Added = Central Lowland sites west of the Mississippi") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10, margin = margin(b = 8)))

# Legends go under the table in the side panel.
p <- patchwork::wrap_plots(A = p_map, B = p_tab, C = patchwork::guide_area(),
                           design = "AB\nAC", widths = c(1.9, 1), heights = c(1, 1.25)) +
  patchwork::plot_layout(guides = "collect") &
  theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(10, 12, 8, 12))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
png_file <- file.path(output_dir, "map_great_plains_expanded.png")
ggsave(png_file, p, width = map_width, height = map_height, dpi = map_dpi, bg = "white")
ggsave(file.path(output_dir, "map_great_plains_expanded.pdf"), p, width = map_width, height = map_height,
       device = if (capabilities("cairo")) cairo_pdf else pdf, bg = "white")
cat("  ", png_file, " (+ .pdf)\n", sep = "")

cat("\n=== STEP 5b COMPLETE ===\n")
cat("Use the files in Outputs/5b_great_plains/ (not 5_master_data) from here on.\n")
