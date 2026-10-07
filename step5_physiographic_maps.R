# =============================================================================
# STEP 5: PHYSIOGRAPHIC MAPS AND SITE-PHYSIOGRAPHY DATABASE (USA + CANADA)
#
# Physiographic layers:
#   USA     Fenneman & Johnson (1946), USGS physio_shp - your local file.
#           Levels: DIVISION (8) > PROVINCE (25) > SECTION (86).
#           Conterminous US only (no Alaska / Hawaii).
#   Canada  Natural Resources Canada, "Physiographic Regions of Canada"
#           (open.canada.ca record a3dfbaf4-1b20-4061-aa0a-e7a79953f52d),
#           downloaded once from NRCan's ArcGIS service and cached.
#           Levels: Region (7) > Subregion (21).
#   Mexico is not included (drawn plain grey).
#
# The two countries' levels are paired like this:
#   "division" level  = US DIVISION  + Canada Region
#   "province" level  = US PROVINCE  + Canada Subregion
#
# Outputs (Outputs/maps/physiographic/):
#   divisions/  5 stage maps at division level
#   provinces/  5 stage maps at province level
#   site_physio_database.csv / .rds
#       One row per site and stage, with its coordinates and its US division,
#       province and section, or Canadian region and subregion.
#       Key columns: Site_Key + Stage_Number (the same keys as Step 3).
#   occurrences_physio.csv / .rds
#       Every FAUNMAP fauna record and PBDB occurrence from Step 3c with the
#       physiographic units of its site attached - ready for analysis.
#   physio_units_key_divisions.csv / physio_units_key_provinces.csv
#       The numbered keys used on the maps.
#
# Reading the maps: each unit carries a NUMBER, listed in the key. Fill
# colours only separate neighbouring units (touching units never share a
# colour). Grey land = no physiographic layer (Alaska, Hawaii, Mexico).
#
# Input:  Outputs/3_stages/ (site_index.csv, faunmap_fauna.csv,
#         pbdb_occurrences.csv) from Step 3c
# Run the whole file (Ctrl+Shift+S in RStudio) after Steps 3 and 3c.
# Needs: dplyr, ggplot2, sf, maps (ggrepel optional, for tidier labels).
# =============================================================================

for (pkg in c("dplyr", "ggplot2", "sf", "maps")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(ggplot2)
library(sf)

sf_use_s2(FALSE)   # planar geometry; all overlay work is done in metres

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

stages_dir  <- file.path(work_dir, "Outputs", "3_stages")
output_dir  <- file.path(work_dir, "Outputs", "maps", "physiographic")
basemap_dir <- file.path(work_dir, "Outputs", "maps", "basemap")   # download cache

# --- USA (local shapefile) ---------------------------------------------------
# The folder holding the shapefile (or the .shp file itself). The folder is
# searched for the shapefile, so its exact name does not matter.
us_path <- "C:/Users/shrut/Downloads/physio_shp"
us_fields <- c(division = "DIVISION", province = "PROVINCE", section = "SECTION")

# --- Canada (NRCan, downloaded) -------------------------------------------------
canada_service <- "https://maps-cartes.services.geo.ca/server_serveur/rest/services/NRCan/phys_reg_en/MapServer"
canada_layers  <- c(division = 0, province = 1)      # 0 = Regions, 1 = Subregions
canada_fields  <- c(division = NA, province = NA)    # NA = detect the name column
# If the download does not work, download the dataset from open.canada.ca and
# give the files here, e.g. c(division = "C:/.../regions.shp", province = "C:/.../subregions.shp")
canada_local   <- c(division = NA, province = NA)

# --- Maps -------------------------------------------------------------------------
map_levels <- c("division", "province")   # one set of 5 maps per level
map_crs    <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
map_lon    <- c(-170, -50)                # extent drawn (degrees)
map_lat    <- c(14, 84)
simplify_m <- 1000                        # smooth outlines for drawing (~1 km)
show_state_lines <- TRUE                  # uses the file cached by Step 4

point_fill   <- "#B2182B"
point_size   <- 2
region_tints <- c("#E9DDB9", "#CFE3C6", "#C8DCEC", "#DCD2EA", "#D2E9E2", "#F1E3C4")
no_data_fill <- "grey88"
sea_fill     <- "#F3F7FA"

snap_km <- 10    # a site just outside every unit (e.g. on the coast) is linked to
                 # the nearest unit within this distance

stage_titles <- c(
  "Bin 1 - Zanclean (4.700-3.600 Ma)",
  "Bin 2 - Piacenzian (3.600-2.580 Ma)",
  "Bin 3 - Gelasian (2.580-1.800 Ma)",
  "Bin 4 - Calabrian (1.800-0.7741 Ma)",
  "Bin 5 - Chibanian (0.7741-0.129 Ma)"
)
stage_names <- c("Zanclean", "Piacenzian", "Gelasian", "Calabrian", "Chibanian")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(basemap_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================

# Read the US shapefile. `path` may be the .shp file or the folder holding it.
# If the exact file is not found, the folder (and its sub-folders) is searched
# for shapefiles and the one whose name contains "physio" is used.
read_us <- function(path) {
  folder <- if (dir.exists(path)) path else dirname(path)
  if (!dir.exists(folder)) stop("Folder not found: ", folder)

  if (!dir.exists(path) && file.exists(path)) {
    shp <- path
  } else {
    found <- list.files(folder, pattern = "\\.shp$", full.names = TRUE,
                        recursive = TRUE, all.files = TRUE, ignore.case = TRUE)
    cat("  Shapefiles found in", folder, ":\n")
    if (length(found) == 0) {
      cat("    (none)\n  Files in the folder:\n")
      cat(paste0("    ", list.files(folder, recursive = TRUE, all.files = TRUE)), sep = "\n")
      stop("No .shp file found in ", folder,
           ". If you downloaded a .zip, unzip it first.")
    }
    cat(paste0("    ", found), sep = "\n")
    pick <- grep("physio", basename(found), ignore.case = TRUE)
    shp <- found[if (length(pick)) pick[1] else 1]
  }
  cat("  Using", shp, "\n")

  x <- tryCatch(st_read(shp, quiet = TRUE), error = function(e) NULL)
  if (!is.null(x)) return(x)

  # A file literally named ".shp" opens on some systems but not others:
  # copy it to a temporary folder under the name "physio.*" and read that.
  stem <- sub("\\.shp$", "", shp, ignore.case = TRUE)
  tmp <- file.path(tempdir(), "physio_copy")
  dir.create(tmp, showWarnings = FALSE)
  for (ext in c("shp", "shx", "dbf", "prj", "cpg")) {
    src <- paste0(stem, ".", ext)
    if (file.exists(src)) file.copy(src, file.path(tmp, paste0("physio.", ext)), overwrite = TRUE)
  }
  cat("    (read via a temporary copy named physio.shp)\n")
  st_read(file.path(tmp, "physio.shp"), quiet = TRUE)
}

# Repair text read with the wrong encoding.
fix_text <- function(x) {
  x <- as.character(x)
  bad <- !is.na(x) & !validUTF8(x)
  x[bad] <- iconv(x[bad], "latin1", "UTF-8")
  Encoding(x) <- "UTF-8"
  x
}

# "INTERIOR PLAINS" -> "Interior Plains"; other names are left as written.
tidy_name <- function(x) {
  x <- gsub("\\s+", " ", trimws(fix_text(x)))
  upper <- !is.na(x) & x == toupper(x) & grepl("[A-Z]", x)
  if (any(upper)) x[upper] <- tools::toTitleCase(tolower(x[upper]))
  x[!is.na(x) & x == ""] <- NA
  x
}

# Print each attribute column with its number of distinct values and an example.
describe_fields <- function(df) {
  att <- st_drop_geometry(df)
  for (col in names(att)) {
    v <- att[[col]]
    cat(sprintf("      %-15s %-10s %4d distinct   e.g. %s\n", col, class(v)[1],
                n_distinct(v, na.rm = TRUE),
                substr(fix_text(as.character(v[!is.na(v)][1])), 1, 40)))
  }
}

find_field <- function(df, preferred, label) {
  cols <- names(st_drop_geometry(df))
  for (p in preferred[!is.na(preferred)]) {
    hit <- cols[toupper(cols) == toupper(p)]
    if (length(hit)) return(hit[1])
  }
  att <- st_drop_geometry(df)
  for (col in cols) {                      # fallback: first column of names
    v <- att[[col]]
    if (is.character(v) && n_distinct(v, na.rm = TRUE) >= 2) return(col)
  }
  stop(label, ": no name column found. Columns: ", paste(cols, collapse = ", "))
}

# Download one layer of an ArcGIS REST MapServer as GeoJSON (with paging).
arcgis_download <- function(service, layer, dest) {
  base <- paste0(service, "/", layer, "/query?where=1%3D1&outFields=*",
                 "&returnGeometry=true&outSR=4326",
                 "&maxAllowableOffset=0.002&geometryPrecision=5&f=geojson")
  fetch <- function(url) {
    tmp <- tempfile(fileext = ".geojson")
    suppressWarnings(download.file(url, tmp, mode = "wb", quiet = TRUE))
    st_read(tmp, quiet = TRUE)
  }
  pages <- list()
  offset <- 0
  repeat {
    url <- paste0(base, "&resultOffset=", offset, "&resultRecordCount=1000")
    page <- tryCatch(fetch(url), error = function(e) NULL)
    if (is.null(page) && offset == 0) page <- fetch(base)   # server without paging
    if (is.null(page) || nrow(page) == 0) break
    pages[[length(pages) + 1]] <- page
    if (nrow(page) < 1000) break
    offset <- offset + 1000
  }
  if (length(pages) == 0) stop("the service returned no features")
  out <- do.call(rbind, pages)
  st_write(out, dest, quiet = TRUE, delete_dsn = TRUE)
  out
}

# Keep only polygon parts, as MULTIPOLYGON.
polygons_only <- function(x) {
  if (any(st_geometry_type(x) == "GEOMETRYCOLLECTION")) {
    x <- st_collection_extract(x, "POLYGON", warn = FALSE)
  }
  x <- x[st_geometry_type(x) %in% c("POLYGON", "MULTIPOLYGON"), ]
  st_cast(x, "MULTIPOLYGON", warn = FALSE)
}

# One row per named unit, in the map projection (full detail, not simplified).
standardise <- function(layer, field, country) {
  layer %>%
    st_make_valid() %>%
    st_transform(map_crs) %>%
    transmute(Country = country, Unit = tidy_name(.data[[field]])) %>%
    filter(!is.na(Unit)) %>%
    polygons_only() %>%
    group_by(Country, Unit) %>%
    summarise(.groups = "drop") %>%
    st_make_valid() %>%
    polygons_only()
}

# Link each point to the unit it falls in (or the nearest within snap_km).
assign_units <- function(points, polys) {
  hit <- st_intersects(points, polys)
  idx <- vapply(hit, function(h) if (length(h)) h[1] else NA_integer_, integer(1))
  how <- ifelse(is.na(idx), "Outside all units", "Inside")
  out <- which(is.na(idx))
  if (length(out) > 0 && nrow(polys) > 0) {
    near <- st_nearest_feature(points[out, ], polys)
    d_km <- as.numeric(st_distance(points[out, ], polys[near, ], by_element = TRUE)) / 1000
    ok <- d_km <= snap_km
    idx[out[ok]] <- near[ok]
    how[out[ok]] <- sprintf("Nearest unit (%.1f km)", d_km[ok])
  }
  data.frame(Country = polys$Country[idx], Unit = polys$Unit[idx], How = how,
             stringsAsFactors = FALSE)
}

read_text_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
           colClasses = "character", na.strings = c("", "NA"))
}

# =============================================================================
# 1. LOAD THE US PHYSIOGRAPHIC SHAPEFILE
# =============================================================================

cat("=== 1. USA: FENNEMAN & JOHNSON (1946) ===\n")
us_raw <- read_us(us_path)

if (is.na(st_crs(us_raw))) {
  # physio.shp is often distributed without a .prj file. Its documented
  # projection is Albers Equal Area (29.5, 45.5, origin 23N 96W), NAD27.
  bb <- st_bbox(us_raw)
  if (all(abs(bb[c("xmin", "xmax")]) <= 180) && all(abs(bb[c("ymin", "ymax")]) <= 90)) {
    st_crs(us_raw) <- 4269
    cat("  No .prj file: coordinates are degrees -> NAD83 geographic assumed\n")
  } else {
    st_crs(us_raw) <- paste("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96",
                            "+x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs")
    cat("  No .prj file: USGS Albers Equal Area (NAD27) assumed, as documented\n")
  }
}

bb <- st_bbox(st_transform(us_raw, 4326))
cat(sprintf("  %d polygons | extent %.1f to %.1f lon, %.1f to %.1f lat\n",
            nrow(us_raw), bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))
if (bb["xmin"] < -130 || bb["xmax"] > -60 || bb["ymin"] < 20 || bb["ymax"] > 53) {
  warning("The US layer is not where the conterminous US should be - check st_crs().")
}
describe_fields(us_raw)

us_cols <- sapply(us_fields, function(f) {
  hit <- names(us_raw)[toupper(names(us_raw)) == toupper(f)]
  if (length(hit) == 0) stop("Column '", f, "' not found in the US shapefile.")
  hit[1]
})

us_units <- lapply(names(us_cols), function(lv) standardise(us_raw, us_cols[[lv]], "USA"))
names(us_units) <- names(us_cols)

cat("\n  US physiographic divisions and their provinces:\n")
us_tree <- st_drop_geometry(us_raw) %>%
  transmute(Division = tidy_name(.data[[us_cols[["division"]]]]),
            Province = tidy_name(.data[[us_cols[["province"]]]]),
            Section  = tidy_name(.data[[us_cols[["section"]]]])) %>%
  filter(!is.na(Division)) %>%
  distinct() %>%
  arrange(Division, Province, Section)
for (d in unique(us_tree$Division)) {
  cat("    ", d, "\n", sep = "")
  for (p in unique(us_tree$Province[us_tree$Division == d & !is.na(us_tree$Province)])) {
    cat("        - ", p, "\n", sep = "")
  }
}
cat(sprintf("  Totals: %d divisions, %d provinces, %d sections\n",
            n_distinct(us_tree$Division), n_distinct(us_tree$Province, na.rm = TRUE),
            n_distinct(us_tree$Section, na.rm = TRUE)))

# =============================================================================
# 2. LOAD THE CANADIAN PHYSIOGRAPHIC LAYERS
# =============================================================================

cat("\n=== 2. CANADA: NRCan PHYSIOGRAPHIC REGIONS ===\n")

canada_preferred <- list(
  division = c("REGION_EN", "REGION", "REGION_NAME", "ENGLISH_NAME", "NAME_EN", "NAME_E", "NAME"),
  province = c("SUBREGION_EN", "SUBREGION", "SUB_REGION", "SUBREGION_NAME",
               "ENGLISH_NAME", "NAME_EN", "NAME_E", "NAME")
)

ca_units <- list()
for (lv in names(canada_layers)) {
  label <- c(division = "Regions", province = "Subregions")[[lv]]
  cat("  ", label, "\n", sep = "")
  layer <- tryCatch({
    if (!is.na(canada_local[[lv]])) {
      cat("    Reading", canada_local[[lv]], "\n")
      st_read(canada_local[[lv]], quiet = TRUE)
    } else {
      cache <- file.path(basemap_dir, paste0("canada_physio_layer", canada_layers[[lv]], ".geojson"))
      if (!file.exists(cache)) {
        cat("    Downloading from NRCan (one time)...\n")
        arcgis_download(canada_service, canada_layers[[lv]], cache)
      }
      cat("    Reading", cache, "\n")
      st_read(cache, quiet = TRUE)
    }
  }, error = function(e) {
    cat("    NOT LOADED:", conditionMessage(e), "\n",
        "   Download 'Physiographic Regions of Canada' from open.canada.ca and set canada_local.\n")
    NULL
  })
  if (is.null(layer)) next

  describe_fields(layer)
  f <- find_field(layer, c(canada_fields[[lv]], canada_preferred[[lv]]), paste("Canada", label))
  ca_units[[lv]] <- standardise(layer, f, "Canada")
  cat(sprintf("    Using '%s': %d units\n", f, nrow(ca_units[[lv]])))
}

# Combined layers per level (full detail, used for linking sites).
physio_levels <- list(
  division = rbind(us_units$division, ca_units$division),
  province = rbind(us_units$province, ca_units$province),
  section  = us_units$section                       # US only
)

# =============================================================================
# 3. SITES AND THEIR PHYSIOGRAPHIC UNITS (THE DATABASE)
#
# One row per site and stage, keyed by Site_Key (as in Step 3):
#   FAUNMAP site = period + Machine Number + Analysis Unit
#   PBDB site    = collection_name
# =============================================================================

cat("\n=== 3. LINKING SITES TO PHYSIOGRAPHIC UNITS ===\n")

index_path <- file.path(stages_dir, "site_index.csv")
if (!file.exists(index_path)) stop("File not found:\n  ", index_path, "\nRun Steps 3 and 3c first.")

site_index <- read.csv(index_path, check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))

sites <- site_index %>%
  mutate(Latitude = suppressWarnings(as.numeric(Latitude)),
         Longitude = suppressWarnings(as.numeric(Longitude)),
         Has_Coordinates = !is.na(Latitude) & !is.na(Longitude)) %>%
  arrange(Stage_Number, Site_Key, desc(Has_Coordinates)) %>%
  distinct(Stage_Number, Site_Key, .keep_all = TRUE) %>%
  select(Stage_Number, Stage, Database, Dataset, SiteName, Site_Key,
         Latitude, Longitude, Has_Coordinates)

sites_sf <- sites %>%
  filter(Has_Coordinates) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(map_crs)

div  <- assign_units(sites_sf, physio_levels$division)
prov <- assign_units(sites_sf, physio_levels$province)
sect <- assign_units(sites_sf, physio_levels$section)

linked <- st_drop_geometry(sites_sf) %>%
  mutate(
    Physio_Country  = coalesce(prov$Country, div$Country),
    US_Division     = if_else(div$Country %in% "USA", div$Unit, NA_character_),
    US_Province     = if_else(prov$Country %in% "USA", prov$Unit, NA_character_),
    US_Section      = sect$Unit,
    CA_Region       = if_else(div$Country %in% "Canada", div$Unit, NA_character_),
    CA_Subregion    = if_else(prov$Country %in% "Canada", prov$Unit, NA_character_),
    # Harmonised columns: US division / Canadian region, and
    # US province / Canadian subregion.
    Physio_Division = div$Unit,
    Physio_Province = prov$Unit,
    Physio_Match    = prov$How
  )

site_physio_database <- sites %>%
  left_join(select(linked, Stage_Number, Site_Key, Physio_Country:Physio_Match),
            by = c("Stage_Number", "Site_Key")) %>%
  mutate(Physio_Match = if_else(Has_Coordinates, Physio_Match, "No coordinates")) %>%
  arrange(Stage_Number, Database, SiteName, Site_Key)

stopifnot(nrow(site_physio_database) == nrow(sites))

match_summary <- count(site_physio_database, Physio_Country, Physio_Match)
print(as.data.frame(match_summary), row.names = FALSE)
cat("  (sites outside all units are in areas without a layer: Mexico, Alaska, offshore)\n")

save_table <- function(df, name) {
  write.csv(df, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE, na = "")
  saveRDS(df, file.path(output_dir, paste0(name, ".rds")))
  cat(sprintf("  %-26s %7d rows -> %s.csv / .rds\n", name, nrow(df), file.path(output_dir, name)))
}

cat("\n  Saving the database:\n")
save_table(site_physio_database, "site_physio_database")

# Every fauna record / occurrence with the physiographic units of its site.
fauna_path <- file.path(stages_dir, "faunmap_fauna.csv")
pbdb_path  <- file.path(stages_dir, "pbdb_occurrences.csv")
occurrences_physio <- NULL
if (file.exists(fauna_path) && file.exists(pbdb_path)) {
  occ <- bind_rows(
    read_text_csv(fauna_path) %>% mutate(Database = "FAUNMAP"),
    read_text_csv(pbdb_path) %>% mutate(Database = "PBDB")
  ) %>%
    mutate(Stage_Number = as.integer(Stage_Number))

  physio_cols <- site_physio_database %>%
    select(Stage_Number, Site_Key, Physio_Country:Physio_Match)

  # Drop any same-named columns first so the join adds clean columns.
  occurrences_physio <- occ %>%
    select(-any_of(setdiff(names(physio_cols), c("Stage_Number", "Site_Key")))) %>%
    left_join(physio_cols, by = c("Stage_Number", "Site_Key")) %>%
    relocate(Database, Stage_Number, Site_Key, Physio_Country, Physio_Division,
             Physio_Province, .before = 1)

  stopifnot(nrow(occurrences_physio) == nrow(occ))
  save_table(occurrences_physio, "occurrences_physio")
} else {
  cat("  NOTE: faunmap_fauna.csv / pbdb_occurrences.csv not found - occurrence table skipped.\n")
}

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(site_physio_database, file.path(output_dir, "site_physio_database.xlsx"))
}

# =============================================================================
# 4. BASE LAYERS FOR THE MAPS
# =============================================================================

frame_ll <- st_as_sfc(st_bbox(c(xmin = map_lon[1], xmax = map_lon[2],
                                ymin = map_lat[1], ymax = map_lat[2]), crs = 4326))

countries <- st_as_sf(maps::map("world", regions = c("Canada", "USA", "Mexico"),
                                fill = TRUE, plot = FALSE))
st_crs(countries) <- NA       # the maps package labels its lon/lat data with an old
st_crs(countries) <- 4326     # ellipsoid; relabel it so it lines up with WGS84
countries <- countries %>% st_make_valid() %>% st_set_agr("constant")
countries <- suppressMessages(st_intersection(countries, frame_ll)) %>% st_transform(map_crs)
lims <- st_bbox(countries)

state_lines <- NULL
state_file <- file.path(basemap_dir, "ne_10m_admin_1_states_provinces_lines.geojson")
if (show_state_lines) {
  if (file.exists(state_file)) {
    state_lines <- st_read(state_file, quiet = TRUE) %>%
      filter(ADM0_A3 %in% c("CAN", "USA", "MEX")) %>%
      st_transform(map_crs)
  } else {
    cat("\n  NOTE: state/province lines skipped - run Step 4 once to download them.\n")
  }
}

# One point per site name per stage on the maps (as in Step 4).
map_points <- sites_sf %>%
  mutate(Map_Site = if_else(is.na(SiteName) | SiteName == "", Site_Key, SiteName)) %>%
  group_by(Stage_Number, Database, Map_Site) %>%
  slice(1) %>%
  ungroup()

use_repel <- requireNamespace("ggrepel", quietly = TRUE)

# =============================================================================
# 5. DRAW THE MAPS (ONE SET OF FIVE PER LEVEL)
# =============================================================================

country_code <- c(USA = "US", Canada = "CA")

for (lv in map_levels) {
  cat("\n=== 5. MAPS AT", toupper(lv), "LEVEL ===\n")
  level_dir <- file.path(output_dir, paste0(lv, "s"))
  dir.create(level_dir, showWarnings = FALSE)

  polys <- physio_levels[[lv]]
  if (simplify_m > 0) polys <- st_simplify(polys, dTolerance = simplify_m, preserveTopology = TRUE)

  # Number the units: US first, then Canada; same-named units share a number.
  key <- st_drop_geometry(polys) %>%
    group_by(Unit) %>%
    summarise(first = min(match(Country, c("USA", "Canada"))),
              Countries = paste(country_code[sort(unique(Country), decreasing = TRUE)],
                                collapse = ", "),
              .groups = "drop") %>%
    arrange(first, Unit) %>%
    mutate(Number = row_number(),
           Key_Label = sprintf("%2d  %s (%s)", Number, Unit, Countries)) %>%
    select(Number, Unit, Countries, Key_Label)
  write.csv(select(key, Number, Unit, Countries),
            file.path(output_dir, paste0("physio_units_key_", lv, "s.csv")), row.names = FALSE)

  polys <- polys %>%
    group_by(Unit) %>%
    summarise(.groups = "drop") %>%
    polygons_only() %>%
    left_join(key, by = "Unit")

  # Map colouring: touching units (within 20 km) get different tints.
  touching <- st_is_within_distance(polys, polys, dist = 20000)
  tint <- integer(nrow(polys))
  for (i in order(lengths(touching), decreasing = TRUE)) {
    used <- tint[setdiff(touching[[i]], i)]
    tint[i] <- which(!seq_along(region_tints) %in% used)[1]
    if (is.na(tint[i])) tint[i] <- 1L
  }
  polys$Tint <- factor(tint)

  # Number labels inside the largest piece of each unit.
  label_pts <- polys %>%
    select(Number) %>%
    st_set_agr("constant") %>%
    st_cast("POLYGON", warn = FALSE) %>%
    mutate(area = as.numeric(st_area(.))) %>%
    group_by(Number) %>%
    slice_max(area, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    st_set_agr("constant") %>%
    st_point_on_surface()
  label_xy <- cbind(st_drop_geometry(label_pts)["Number"], st_coordinates(label_pts))
  label_xy$Key <- key$Key_Label[match(label_xy$Number, key$Number)]

  key_cols <- if (nrow(key) > 30) 2 else 1

  for (s in 1:5) {
    pts <- filter(map_points, Stage_Number == s)

    p <- ggplot() +
      geom_sf(data = countries, fill = no_data_fill, colour = NA) +
      geom_sf(data = polys, aes(fill = Tint), colour = "grey35", linewidth = 0.2) +
      scale_fill_manual(values = region_tints, guide = "none")

    if (!is.null(state_lines)) {
      p <- p + geom_sf(data = state_lines, colour = "grey30", linewidth = 0.1, alpha = 0.35)
    }
    p <- p + geom_sf(data = countries, fill = NA, colour = "grey20", linewidth = 0.3)

    if (use_repel) {
      p <- p + ggrepel::geom_text_repel(
        data = label_xy, aes(x = X, y = Y, label = Number),
        size = 2.3, fontface = "bold", colour = "grey15",
        bg.color = "white", bg.r = 0.15, min.segment.length = 0.3,
        segment.colour = "grey30", segment.size = 0.2, max.overlaps = Inf, seed = 1
      )
    } else {
      p <- p + geom_text(data = label_xy, aes(x = X, y = Y, label = Number),
                         size = 2.3, fontface = "bold", colour = "grey15")
    }

    # Invisible layer that turns the numbered key into the legend.
    p <- p +
      geom_point(data = label_xy, aes(x = X, y = Y, shape = Key), alpha = 0) +
      scale_shape_manual(values = rep(32, nrow(key)), breaks = key$Key_Label,
                         name = paste("Physiographic", lv, "units")) +
      guides(shape = guide_legend(ncol = key_cols, override.aes = list(alpha = 0)))

    # Sites on top: solid points with a thin white ring.
    p <- p +
      geom_sf(data = pts, shape = 21, fill = point_fill, colour = "white",
              stroke = 0.3, size = point_size) +
      coord_sf(crs = map_crs, xlim = lims[c("xmin", "xmax")], ylim = lims[c("ymin", "ymax")],
               expand = FALSE) +
      labs(title = stage_titles[s],
           subtitle = paste0(nrow(pts), " sites (FAUNMAP + PBDB) on physiographic ", lv,
                             "s of the USA and Canada"),
           x = NULL, y = NULL) +
      theme_bw(base_size = 10) +
      theme(
        panel.background = element_rect(fill = sea_fill),
        panel.grid = element_line(colour = "white", linewidth = 0.3),
        plot.title = element_text(face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 6.5, family = "mono"),
        legend.title = element_text(size = 8, face = "bold"),
        legend.key.size = unit(0.3, "lines")
      )

    width <- if (key_cols == 2) 14 else 11
    file <- file.path(level_dir, sprintf("physio_%s_stage%d_%s.png", lv, s, stage_names[s]))
    ggsave(file, p, width = width, height = 7.5, dpi = 300)
    cat(sprintf("  %-12s %5d sites -> %s\n", stage_names[s], nrow(pts), file))
  }
}

cat("\n=== STEP 5 COMPLETE ===\n")
cat("Database: site_physio_database (one row per site and stage) and\n",
    "occurrences_physio (every fauna record / occurrence), in\n  ", output_dir, "\n", sep = "")
cat("Objects in your Environment: site_physio_database, occurrences_physio, us_tree\n")
