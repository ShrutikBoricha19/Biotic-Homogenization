# =============================================================================
# STEP 4: MAP ALL SITES FOR EACH GEOLOGICAL STAGE
#
# One map per stage (5 maps) of Canada, the USA and Mexico, showing every
# FAUNMAP and PBDB site in that stage. Both databases are drawn the same way.
#
# Each site is plotted ONCE per stage:
#   FAUNMAP site = SiteName      (all analysis units of a site -> one point)
#   PBDB site    = collection_name
# If a site has several coordinate records, the first one is used.
#
# Input:  Outputs/3_stages/site_index.csv  (from Step 3)
# Output: Outputs/maps/  (one PNG per stage + the plotted points as CSV)
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 3.
# Needs the packages dplyr, ggplot2 and maps.
# =============================================================================

for (pkg in c("dplyr", "ggplot2", "maps")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(ggplot2)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

site_index_file <- file.path("Outputs", "3_stages", "site_index.csv")

output_dir <- file.path(work_dir, "Outputs", "maps")

# Map extent (degrees) covering Canada, the USA (incl. Alaska) and Mexico.
map_xlim <- c(-170, -50)
map_ylim <- c(14, 84)

point_colour <- "#B2182B"
point_size   <- 1.6

stage_titles <- c(
  "Bin 1 - Zanclean (4.700-3.600 Ma)",
  "Bin 2 - Piacenzian (3.600-2.580 Ma)",
  "Bin 3 - Gelasian (2.580-1.800 Ma)",
  "Bin 4 - Calabrian (1.800-0.7741 Ma)",
  "Bin 5 - Chibanian (0.7741-0.129 Ma)"
)
stage_names <- c("Zanclean", "Piacenzian", "Gelasian", "Calabrian", "Chibanian")

# =============================================================================
# 1. READ THE SITE INDEX FROM STEP 3
# =============================================================================

cat("=== 1. READING SITE INDEX ===\n")

path <- file.path(work_dir, site_index_file)
if (!file.exists(path)) {
  stop("File not found:\n  ", path, "\nRun Step 3 first.")
}

site_index <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))

needed_cols <- c("Database", "SiteName", "Site_Key", "Latitude", "Longitude", "Stage_Number")
missing_cols <- setdiff(needed_cols, names(site_index))
if (length(missing_cols) > 0) {
  stop("site_index.csv is missing: ", paste(missing_cols, collapse = ", "))
}

cat(sprintf("  %d site-stage rows read from %s\n", nrow(site_index), site_index_file))

# =============================================================================
# 2. ONE POINT PER SITE PER STAGE
# =============================================================================

cat("\n=== 2. ONE POINT PER SITE ===\n")

map_points_all <- site_index %>%
  mutate(
    Latitude = suppressWarnings(as.numeric(Latitude)),
    Longitude = suppressWarnings(as.numeric(Longitude)),
    SiteName = trimws(SiteName),
    # Sites without a name fall back to their unique key.
    Map_Site = if_else(is.na(SiteName) | SiteName == "", Site_Key, SiteName),
    Has_Coordinates = !is.na(Latitude) & !is.na(Longitude)
  ) %>%
  arrange(Stage_Number, Database, Map_Site, desc(Has_Coordinates)) %>%
  group_by(Stage_Number, Database, Map_Site) %>%
  summarise(
    Latitude = first(Latitude),          # first record with coordinates
    Longitude = first(Longitude),
    Has_Coordinates = first(Has_Coordinates),
    n_index_rows = n(),
    .groups = "drop"
  ) %>%
  mutate(Stage = stage_names[Stage_Number])

map_points <- filter(map_points_all, Has_Coordinates)

no_coords <- filter(map_points_all, !Has_Coordinates)
outside <- filter(map_points, Longitude < map_xlim[1] | Longitude > map_xlim[2] |
                    Latitude < map_ylim[1] | Latitude > map_ylim[2])

counts <- map_points_all %>%
  group_by(Stage_Number, Stage) %>%
  summarise(
    FAUNMAP_sites = sum(Database == "FAUNMAP" & Has_Coordinates),
    PBDB_sites = sum(Database == "PBDB" & Has_Coordinates),
    Total_plotted = sum(Has_Coordinates),
    Without_coordinates = sum(!Has_Coordinates),
    .groups = "drop"
  )
print(as.data.frame(counts), row.names = FALSE)

if (nrow(outside) > 0) {
  cat(sprintf("\n  NOTE: %d sites fall outside the map extent and will not be visible.\n",
              nrow(outside)))
}

# =============================================================================
# 3. BASE MAP: CANADA, USA AND MEXICO
# =============================================================================

cat("\n=== 3. DRAWING MAPS ===\n")

base_map <- map_data("world", region = c("Canada", "USA", "Mexico"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (s in 1:5) {
  pts <- filter(map_points, Stage_Number == s)

  p <- ggplot() +
    geom_polygon(data = base_map, aes(x = long, y = lat, group = group),
                 fill = "grey92", colour = "grey55", linewidth = 0.25) +
    geom_point(data = pts, aes(x = Longitude, y = Latitude),
               colour = point_colour, size = point_size, alpha = 0.8) +
    coord_quickmap(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(
      title = stage_titles[s],
      subtitle = paste0(nrow(pts), " sites (FAUNMAP + PBDB)"),
      x = "Longitude", y = "Latitude"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.background = element_rect(fill = "aliceblue"),
      panel.grid = element_line(colour = "white", linewidth = 0.3),
      plot.title = element_text(face = "bold")
    )

  file <- file.path(output_dir, sprintf("map_stage%d_%s.png", s, stage_names[s]))
  ggsave(file, p, width = 8, height = 7, dpi = 300)
  cat(sprintf("  %-12s %5d sites -> %s\n", stage_names[s], nrow(pts), file))
}

# =============================================================================
# 4. SAVE THE PLOTTED POINTS
# =============================================================================

write.csv(select(map_points_all, Stage_Number, Stage, Database, Site = Map_Site,
                 Latitude, Longitude, Has_Coordinates, n_index_rows),
          file.path(output_dir, "map_points.csv"), row.names = FALSE, na = "")
cat("\n  Points (including sites without coordinates) -> ",
    file.path(output_dir, "map_points.csv"), "\n", sep = "")

cat("\n=== STEP 4 COMPLETE ===\n")
cat("Objects in your Environment: map_points, map_points_all, counts\n")
