# =============================================================================
# STEP 3: SPATIAL BINNING OF FAUNMAP AND PBDB SITES
#
# Every FAUNMAP locality and PBDB collection is placed in a spatial unit at
# three levels, each saved in its own sub-folder:
#
#   political/  US state, Canadian province / territory, Mexican state
#               (Natural Earth 1:10m admin-1, downloaded once and cached)
#   division/   US physiographic DIVISION (Fenneman & Johnson 1946, physio_shp)
#               + Canadian physiographic Region (NRCan)
#   province/   US physiographic PROVINCE (physio_shp)
#               + Canadian physiographic Subregion (NRCan)
#               (the US SECTION is stored too, as US_Section)
#
# A site takes the unit its coordinates fall in; a site just off the coast
# takes the nearest unit within 'snap_km'. Mexico has no physiographic layer,
# so Mexican sites have a state but no physiographic division / province.
# Units of the same name in both countries (e.g. Interior Plains) are one unit.
#
# Removed before binning (listed in excluded_sites.csv):
#   - sites in Alaska and Hawaii
#   - Canadian sites north of 'canada_north_limit' (60 degrees N: north of
#     British Columbia, Alberta, Saskatchewan and Manitoba)
# Sites without coordinates are kept but get no spatial unit.
#
# Outputs (Outputs/3_spatial/):
#   <period>_localities.csv, <period>_fauna.csv, pbdb_clean.csv
#       The Step 2b files minus the removed sites, with the spatial columns
#       added (Country, State_Province, Physio_Division, Physio_Province,
#       US_Section, ...). Step 3c (time binning) reads these.
#   sites_spatial.csv      one row per site with all spatial units
#   excluded_sites.csv     removed sites and why
#   spatial_layers.rds     simplified spatial layers, used by the Step 4 maps
#   <level>/sites_<level>.csv          site -> unit, with the unit number
#   <level>/unit_summary_<level>.csv   sites, records and taxa per unit
#   <level>/map_<level>.png / .pdf     map with a numbered key (no Hawaii)
#
# Inputs: Outputs/2b_clean/ (Step 2b), the US physiographic shapefile
# (physio_shp, set 'us_path'), and two downloads that are cached in
# Outputs/maps/basemap/ (NRCan physiographic regions, Natural Earth admin-1).
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 2b; then Step 3c.
# Needs: dplyr, ggplot2, sf, maps, patchwork (ggrepel optional).
# =============================================================================

for (pkg in c("dplyr", "ggplot2", "sf", "maps", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(ggplot2)
library(sf)

suppressMessages(sf_use_s2(FALSE))   # planar geometry; overlays are done in metres

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

# Cleaned files from Step 2b.
periods <- c("Blancan", "Irvingtonian", "Rancholabrean")
loc_files   <- setNames(file.path("Outputs", "2b_clean", paste0(tolower(periods), "_localities.csv")), periods)
fauna_files <- setNames(file.path("Outputs", "2b_clean", paste0(tolower(periods), "_fauna.csv")), periods)
pbdb_file   <- file.path("Outputs", "2b_clean", "pbdb_clean.csv")

output_dir  <- file.path(work_dir, "Outputs", "3_spatial")
basemap_dir <- file.path(work_dir, "Outputs", "maps", "basemap")   # download cache (shared with Steps 4-5)

save_excel <- TRUE    # also save .xlsx copies of the main tables (needs writexl)

# --- USA physiography (local shapefile) --------------------------------------
# The folder holding physio.shp (or the .shp file itself).
us_path   <- "C:/Users/shrut/Downloads/physio_shp"
us_fields <- c(division = "DIVISION", province = "PROVINCE", section = "SECTION")

# --- Canada physiography (NRCan, downloaded once) -----------------------------
canada_service <- "https://maps-cartes.services.geo.ca/server_serveur/rest/services/NRCan/phys_reg_en/MapServer"
canada_layers  <- c(division = 0, province = 1)      # 0 = Regions, 1 = Subregions
canada_fields  <- c(division = NA, province = NA)    # NA = detect the name column
# If the download fails, download "Physiographic Regions of Canada" from
# open.canada.ca and give the files, e.g. c(division = "C:/.../regions.shp", province = "C:/.../subregions.shp")
canada_local   <- c(division = NA, province = NA)

# --- States and provinces (Natural Earth, downloaded once, ~40 MB) -------------
admin1_url   <- paste0("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/",
                       "geojson/ne_10m_admin_1_states_provinces.geojson")
admin1_local <- NA      # or the path of a local copy of that file

# --- Geographic limits ----------------------------------------------------------
remove_alaska      <- TRUE
remove_hawaii      <- TRUE
canada_north_limit <- 60      # degrees N; Canadian sites north of this are removed (NA = keep all)
snap_km            <- 25      # a site off the coast takes the nearest unit within this distance

# --- Maps -------------------------------------------------------------------------
map_crs  <- "+proj=laea +lat_0=45 +lon_0=-98 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
work_crs <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +datum=WGS84 +units=m +no_defs"
map_lon  <- c(-136, -52)      # extent drawn (degrees); Hawaii lies outside it
map_lat  <- c(14, 61)
simplify_m <- 1000            # smoother outlines for drawing (~1 km)
map_width  <- 16              # inches (map + key)
map_height <- 9
map_dpi    <- 300

dataset_colours <- c("FAUNMAP Blancan" = "#E69F00", "FAUNMAP Irvingtonian" = "#0072B2",
                     "FAUNMAP Rancholabrean" = "#009E73", "PBDB" = "#CC79A7")
dataset_shapes  <- c("FAUNMAP Blancan" = 21, "FAUNMAP Irvingtonian" = 22,
                     "FAUNMAP Rancholabrean" = 24, "PBDB" = 23)
point_size   <- 2.2
site_fill    <- "#B2182B"    # all sites drawn alike, whatever their database
unit_tints   <- c("#E9DDB9", "#CFE3C6", "#C8DCEC", "#DCD2EA", "#D2E9E2", "#F1D9C9")
empty_fill   <- "grey95"     # units without sites
no_data_fill <- "grey85"     # land without a layer (e.g. Mexico on physiographic maps)
sea_fill     <- "#EEF4F8"

level_titles <- c(political = "States and provinces",
                  division  = "Physiographic divisions (USA) and regions (Canada)",
                  province  = "Physiographic provinces (USA) and subregions (Canada)")

# =============================================================================
# HELPERS
# =============================================================================

clean_text <- function(x) { x <- trimws(as.character(x)); x[x == ""] <- NA_character_; x }
# Machine Number as a plain integer: "1234.00", "1,234.00", " 1234 " and 1234 all become "1234".
machine_key <- function(x) {
  x <- gsub("[,[:space:]]", "", clean_text(x))
  v <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(v) & v == round(v), format(round(v), scientific = FALSE, trim = TRUE), x)
}
# Analysis Unit compared ignoring capitals, extra spaces and a "." typed for a ","
# before a space ("CU 29.5. L6604" = "CU 29.5, L6604"; "Assemblage" = "assemblage ").
analysis_key <- function(x) {
  x <- tolower(gsub("\\s+", " ", clean_text(x)))
  gsub("[.,;] ", ", ", x)
}
to_num   <- function(x) suppressWarnings(as.numeric(clean_text(x)))
to_coord <- function(x, limit) { v <- to_num(x); v[!is.na(v) & abs(v) > limit] <- NA; v }

simplify_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))
column_options <- list(
  machine  = c("machinenumber", "machineno", "machine"),
  analysis = c("analysisunit"),
  site     = c("sitename", "collectionname"),
  lat      = c("latdd", "latitude", "lat"),
  lon      = c("longdd", "longitude", "lng", "long", "lon"),
  genus    = c("genus"),
  species  = c("species"),
  taxon    = c("acceptedname", "taxonname", "identifiedname")
)
find_col <- function(df, field, file, required = TRUE) {
  simple <- simplify_name(names(df))
  for (o in column_options[[field]]) {
    hit <- names(df)[simple == o]
    if (length(hit) >= 1) return(hit[1])
  }
  if (required) stop("No '", field, "' column found in ", file,
                     ".\nColumns in this file: ", paste(names(df), collapse = " | "))
  NA_character_
}

read_input <- function(file) {
  df <- read.csv(file.path(work_dir, file), check.names = FALSE, stringsAsFactors = FALSE,
                 colClasses = "character", na.strings = c("", "NA", "N/A"))
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))
  cat(sprintf("  %-45s %7d rows, %3d columns\n", file, nrow(df), ncol(df)))
  df
}

# Repair text read with the wrong encoding.
fix_text <- function(x) {
  x <- as.character(x)
  bad <- !is.na(x) & !validUTF8(x)
  x[bad] <- iconv(x[bad], "latin1", "UTF-8")
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

# Read the US shapefile; `path` may be the .shp file or a folder holding it.
read_us <- function(path) {
  folder <- if (dir.exists(path)) path else dirname(path)
  if (!dir.exists(folder)) stop("Folder not found: ", folder, "\nSet 'us_path' in SETTINGS.")
  if (!dir.exists(path) && file.exists(path)) {
    shp <- path
  } else {
    found <- list.files(folder, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE,
                        all.files = TRUE, ignore.case = TRUE)
    if (length(found) == 0) stop("No .shp file found in ", folder, ". If you downloaded a .zip, unzip it first.")
    pick <- grep("physio", basename(found), ignore.case = TRUE)
    shp <- found[if (length(pick)) pick[1] else 1]
  }
  cat("  Using", shp, "\n")
  x <- tryCatch(st_read(shp, quiet = TRUE), error = function(e) NULL)
  if (!is.null(x)) return(x)
  # A file literally named ".shp" opens on some systems but not others.
  stem <- sub("\\.shp$", "", shp, ignore.case = TRUE)
  tmp <- file.path(tempdir(), "physio_copy")
  dir.create(tmp, showWarnings = FALSE)
  for (ext in c("shp", "shx", "dbf", "prj", "cpg")) {
    src <- paste0(stem, ".", ext)
    if (file.exists(src)) file.copy(src, file.path(tmp, paste0("physio.", ext)), overwrite = TRUE)
  }
  st_read(file.path(tmp, "physio.shp"), quiet = TRUE)
}

find_field <- function(df, preferred, label) {
  cols <- names(st_drop_geometry(df))
  for (p in preferred[!is.na(preferred)]) {
    hit <- cols[toupper(cols) == toupper(p)]
    if (length(hit)) return(hit[1])
  }
  att <- st_drop_geometry(df)
  for (col in cols) if (is.character(att[[col]]) && n_distinct(att[[col]], na.rm = TRUE) >= 2) return(col)
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
  pages <- list(); offset <- 0
  repeat {
    page <- tryCatch(fetch(paste0(base, "&resultOffset=", offset, "&resultRecordCount=1000")),
                     error = function(e) NULL)
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
  if (any(st_geometry_type(x) == "GEOMETRYCOLLECTION")) x <- st_collection_extract(x, "POLYGON", warn = FALSE)
  x <- x[st_geometry_type(x) %in% c("POLYGON", "MULTIPOLYGON"), ]
  st_cast(x, "MULTIPOLYGON", warn = FALSE)
}

# One row per named unit (Country + Unit), in the working projection.
standardise <- function(layer, field, country) {
  suppressWarnings(suppressMessages(
    layer %>%
      st_make_valid() %>%
      st_transform(work_crs) %>%
      transmute(Country = country, Unit = tidy_name(.data[[field]])) %>%
      filter(!is.na(Unit)) %>%
      polygons_only() %>%
      group_by(Country, Unit) %>%
      summarise(.groups = "drop") %>%
      st_make_valid() %>%
      polygons_only()
  ))
}

# Link points to the unit they fall in, or the nearest unit within snap_km.
# Distance_km and Nearest_* describe the nearest unit for every point.
assign_units <- function(points, polys) {
  hit <- st_intersects(points, polys)
  idx <- vapply(hit, function(h) if (length(h)) h[1] else NA_integer_, integer(1))
  near <- idx
  d_km <- ifelse(is.na(idx), NA_real_, 0)
  out <- which(is.na(idx))
  if (length(out) > 0 && nrow(polys) > 0) {
    near[out] <- st_nearest_feature(points[out, ], polys)
    d_km[out] <- as.numeric(st_distance(points[out, ], polys[near[out], ], by_element = TRUE)) / 1000
  }
  inside_snap <- !is.na(near) & d_km <= snap_km
  data.frame(
    Country = ifelse(inside_snap, polys$Country[near], NA_character_),
    Unit    = ifelse(inside_snap, polys$Unit[near], NA_character_),
    Match   = case_when(!is.na(idx) ~ "Inside",
                        inside_snap ~ sprintf("Nearest unit (%.1f km)", d_km),
                        TRUE ~ sprintf("Outside all units (nearest %.0f km)", d_km)),
    Nearest_Country = polys$Country[near], Nearest_Unit = polys$Unit[near],
    Distance_km = round(d_km, 2), stringsAsFactors = FALSE)
}

use_excel <- save_excel && requireNamespace("writexl", quietly = TRUE)
save_out <- function(df, name, dir = output_dir, excel = FALSE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".csv"))
  write.csv(df, path, row.names = FALSE, na = "")
  if (excel && use_excel) writexl::write_xlsx(df, file.path(dir, paste0(name, ".xlsx")))
  cat(sprintf("  %-32s %7d rows -> %s\n", name, nrow(df), path))
}

# Same-named columns already in a file are kept as *_original.
add_columns <- function(df, extra) {
  clash <- intersect(names(df), names(extra))
  if (length(clash)) names(df)[names(df) %in% clash] <- paste0(clash, "_original")
  bind_cols(df, extra)
}

# =============================================================================
# 1. READ THE STEP 2b FILES
# =============================================================================

cat("=== 1. READING THE STEP 2b FILES ===\n")
if (!dir.exists(work_dir)) stop("The folder does not exist:\n  ", work_dir)
needed <- c(unname(loc_files), unname(fauna_files), pbdb_file)
missing <- needed[!file.exists(file.path(work_dir, needed))]
if (length(missing)) stop("Not found (run Steps 1, 2 and 2b first):\n  ", paste(missing, collapse = "\n  "))

loc_raw   <- lapply(loc_files, read_input)
fauna_raw <- lapply(fauna_files, read_input)
pbdb_raw  <- read_input(pbdb_file)

loc_cols <- lapply(periods, function(p) {
  sapply(c("machine", "analysis", "site", "lat", "lon"), function(f) find_col(loc_raw[[p]], f, loc_files[[p]]))
})
names(loc_cols) <- periods
pbdb_cols <- sapply(c("site", "lat", "lon"), function(f) find_col(pbdb_raw, f, pbdb_file))

# Coordinates of every locality row and PBDB occurrence.
loc_xy <- lapply(periods, function(p) {
  c <- loc_cols[[p]]
  data.frame(lat = to_coord(loc_raw[[p]][[c[["lat"]]]], 90), lon = to_coord(loc_raw[[p]][[c[["lon"]]]], 180))
})
names(loc_xy) <- periods
pbdb_xy <- data.frame(lat = to_coord(pbdb_raw[[pbdb_cols[["lat"]]]], 90),
                      lon = to_coord(pbdb_raw[[pbdb_cols[["lon"]]]], 180))

xy_key <- function(xy) ifelse(is.na(xy$lat) | is.na(xy$lon), NA_character_,
                              sprintf("%.6f|%.6f", xy$lat, xy$lon))
points <- bind_rows(c(loc_xy, list(pbdb_xy))) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  distinct() %>%
  mutate(.xy = xy_key(.))
cat(sprintf("  %d distinct site coordinates\n", nrow(points)))

points_sf <- st_as_sf(points, coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(work_crs)

# =============================================================================
# 2. SPATIAL LAYERS
# =============================================================================

cat("\n=== 2. SPATIAL LAYERS ===\n")
dir.create(basemap_dir, recursive = TRUE, showWarnings = FALSE)

# --- States and provinces ---
cat("  States and provinces (Natural Earth admin-1)\n")
admin1_file <- if (!is.na(admin1_local)) admin1_local else file.path(basemap_dir, basename(admin1_url))
if (!file.exists(admin1_file)) {
  cat("    Downloading (one time, ~40 MB)...\n")
  old_timeout <- getOption("timeout"); options(timeout = max(600, old_timeout))
  ok <- tryCatch({ download.file(admin1_url, admin1_file, mode = "wb", quiet = TRUE); TRUE },
                 error = function(e) FALSE)
  options(timeout = old_timeout)
  if (!ok) {
    if (file.exists(admin1_file)) file.remove(admin1_file)
    stop("Download failed. Download ne_10m_admin_1_states_provinces.geojson from\n  ", admin1_url,
         "\nand set 'admin1_local' to its path.")
  }
}
admin1_raw <- st_read(admin1_file, quiet = TRUE)
names(admin1_raw)[tolower(names(admin1_raw)) %in% c("adm0_a3", "name", "postal")] <-
  tolower(names(admin1_raw)[tolower(names(admin1_raw)) %in% c("adm0_a3", "name", "postal")])
country_names <- c(USA = "USA", CAN = "Canada", MEX = "Mexico")
political <- suppressWarnings(suppressMessages(
  admin1_raw %>%
    filter(adm0_a3 %in% names(country_names)) %>%
    st_make_valid() %>%
    st_transform(work_crs) %>%
    transmute(Country = unname(country_names[adm0_a3]), Unit = fix_text(name), Code = postal) %>%
    polygons_only()
))
state_codes <- st_drop_geometry(political) %>% distinct(Country, Unit, Code)
political <- select(political, Country, Unit)
cat(sprintf("    %d states / provinces / territories (USA, Canada, Mexico)\n", nrow(political)))

# --- USA physiography ---
cat("  USA physiography (Fenneman & Johnson 1946)\n")
us_raw <- read_us(us_path)
if (is.na(st_crs(us_raw))) {
  # physio.shp is often distributed without a .prj file. Its documented
  # projection is Albers Equal Area (29.5, 45.5, origin 23N 96W), NAD27.
  bb <- st_bbox(us_raw)
  if (all(abs(bb[c("xmin", "xmax")]) <= 180) && all(abs(bb[c("ymin", "ymax")]) <= 90)) {
    st_crs(us_raw) <- 4269
  } else {
    st_crs(us_raw) <- paste("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96",
                            "+x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs")
  }
  cat("    No .prj file: documented projection assumed\n")
}
us_cols <- sapply(us_fields, function(f) {
  hit <- names(us_raw)[toupper(names(us_raw)) == toupper(f)]
  if (length(hit) == 0) stop("Column '", f, "' not found in the US shapefile.")
  hit[1]
})
us_units <- lapply(names(us_cols), function(lv) standardise(us_raw, us_cols[[lv]], "USA"))
names(us_units) <- names(us_cols)
cat(sprintf("    %d divisions, %d provinces, %d sections\n",
            n_distinct(us_units$division$Unit), n_distinct(us_units$province$Unit),
            n_distinct(us_units$section$Unit)))

# --- Canada physiography ---
cat("  Canada physiography (NRCan)\n")
canada_preferred <- list(
  division = c("REGION_EN", "REGION", "REGION_NAME", "ENGLISH_NAME", "NAME_EN", "NAME_E", "NAME"),
  province = c("SUBREGION_EN", "SUBREGION", "SUB_REGION", "SUBREGION_NAME",
               "ENGLISH_NAME", "NAME_EN", "NAME_E", "NAME"))
ca_units <- list()
for (lv in names(canada_layers)) {
  label <- c(division = "Regions", province = "Subregions")[[lv]]
  layer <- tryCatch({
    if (!is.na(canada_local[[lv]])) {
      st_read(canada_local[[lv]], quiet = TRUE)
    } else {
      cache <- file.path(basemap_dir, paste0("canada_physio_layer", canada_layers[[lv]], ".geojson"))
      if (!file.exists(cache)) {
        cat("    Downloading", label, "from NRCan (one time)...\n")
        arcgis_download(canada_service, canada_layers[[lv]], cache)
      }
      st_read(cache, quiet = TRUE)
    }
  }, error = function(e) {
    cat("    ", label, " NOT LOADED: ", conditionMessage(e), "\n",
        "    Download 'Physiographic Regions of Canada' from open.canada.ca and set canada_local.\n", sep = "")
    NULL
  })
  if (is.null(layer)) next
  f <- find_field(layer, c(canada_fields[[lv]], canada_preferred[[lv]]), paste("Canada", label))
  ca_units[[lv]] <- standardise(layer, f, "Canada")
  cat(sprintf("    %-10s '%s': %d units\n", label, f, nrow(ca_units[[lv]])))
}

layers <- list(
  political = political,
  division  = rbind(us_units$division, ca_units$division),
  province  = rbind(us_units$province, ca_units$province)
)

# Simplified copies of the layers (without Alaska and Hawaii) for the Step 4 maps.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
map_layers <- lapply(c(layers, list(us_division = us_units$division, ca_division = ca_units$division)),
                     function(l) {
                       if (is.null(l)) return(NULL)
                       l <- l[!l$Unit %in% c("Alaska", "Hawaii"), ]
                       suppressWarnings(st_simplify(l, dTolerance = 500, preserveTopology = TRUE))
                     })
map_layers$settings <- list(work_crs = work_crs, canada_north_limit = canada_north_limit,
                            remove_alaska = remove_alaska, remove_hawaii = remove_hawaii)
saveRDS(map_layers, file.path(output_dir, "spatial_layers.rds"))
cat("  Map layers saved for Step 4:", file.path(output_dir, "spatial_layers.rds"), "\n")

# =============================================================================
# 3. ASSIGN EVERY SITE COORDINATE TO ITS UNITS
# =============================================================================

cat("\n=== 3. ASSIGNING SITES TO UNITS ===\n")

pol  <- assign_units(points_sf, layers$political)
div  <- assign_units(points_sf, layers$division)
prov <- assign_units(points_sf, layers$province)
sect <- assign_units(points_sf, us_units$section)

# Exclusions use the nearest state / province, whatever its distance.
point_units <- points %>%
  mutate(
    Country          = pol$Country,
    State_Province   = pol$Unit,
    State_Code       = state_codes$Code[match(paste(pol$Country, pol$Unit),
                                              paste(state_codes$Country, state_codes$Unit))],
    Physio_Division  = div$Unit,
    Physio_Province  = prov$Unit,
    US_Section       = sect$Unit,
    Political_Match  = pol$Match,
    Division_Match   = ifelse(pol$Country %in% "Mexico" & is.na(div$Unit), "No physiographic layer (Mexico)", div$Match),
    Province_Match   = ifelse(pol$Country %in% "Mexico" & is.na(prov$Unit), "No physiographic layer (Mexico)", prov$Match),
    Spatial_Exclusion = case_when(
      remove_alaska & pol$Nearest_Unit %in% "Alaska" ~ "Alaska",
      remove_hawaii & pol$Nearest_Unit %in% "Hawaii" ~ "Hawaii",
      !is.na(canada_north_limit) & pol$Nearest_Country %in% "Canada" & lat > canada_north_limit ~
        paste0("Canada north of ", canada_north_limit, " degrees N"),
      TRUE ~ NA_character_)
  )
cat(sprintf("  %d coordinates | %d inside a state/province | %d snapped to the coast | %d outside all\n",
            nrow(point_units), sum(pol$Match == "Inside"), sum(grepl("^Nearest", pol$Match)),
            sum(grepl("^Outside", pol$Match))))
print(as.data.frame(count(point_units, Spatial_Exclusion)), row.names = FALSE)

spatial_cols <- c("Country", "State_Province", "State_Code", "Physio_Division", "Physio_Province",
                  "US_Section", "Political_Match", "Division_Match", "Province_Match", "Spatial_Exclusion")

units_for <- function(xy) {
  out <- point_units[match(xy_key(xy), point_units$.xy), spatial_cols]
  out$Political_Match[is.na(xy$lat) | is.na(xy$lon)] <- "No coordinates"
  rownames(out) <- NULL
  out
}

# =============================================================================
# 4. ATTACH THE UNITS TO LOCALITIES, FAUNA AND PBDB; REMOVE EXCLUDED SITES
# =============================================================================

cat("\n=== 4. LOCALITIES, FAUNA AND PBDB ===\n")

loc_spatial <- lapply(periods, function(p) add_columns(loc_raw[[p]], units_for(loc_xy[[p]])))
names(loc_spatial) <- periods

# Fauna take the units of their locality (period + Machine Number + Analysis Unit).
fauna_spatial <- lapply(periods, function(p) {
  lc <- loc_cols[[p]]
  lookup <- loc_spatial[[p]] %>%
    mutate(.mk = machine_key(.data[[lc[["machine"]]]]), .ak = analysis_key(.data[[lc[["analysis"]]]])) %>%
    filter(!is.na(.mk), !is.na(.ak)) %>%
    distinct(.mk, .ak, .keep_all = TRUE) %>%
    select(.mk, .ak, all_of(spatial_cols))
  f <- fauna_raw[[p]]
  if (nrow(f) == 0) return(add_columns(f, lookup[0, spatial_cols]))
  fm <- find_col(f, "machine", fauna_files[[p]]); fa <- find_col(f, "analysis", fauna_files[[p]])
  keys <- data.frame(.mk = machine_key(f[[fm]]), .ak = analysis_key(f[[fa]]))
  extra <- left_join(keys, lookup, by = c(".mk", ".ak")) %>% select(all_of(spatial_cols))
  extra$Political_Match[is.na(extra$Political_Match)] <- "No matching locality"
  stopifnot(nrow(extra) == nrow(f))
  add_columns(f, extra)
})
names(fauna_spatial) <- periods

pbdb_spatial <- add_columns(pbdb_raw, units_for(pbdb_xy))

# Excluded sites (one row per site).
site_rows <- function(df, cols, period, dataset) {
  # Site keys as in Step 3c: FAUNMAP = period | Machine Number | Analysis Unit,
  # PBDB = "PBDB" | collection name.
  key <- if (is.na(period)) {
    paste("PBDB", clean_text(df[[cols[["site"]]]]), sep = " | ")
  } else {
    paste(period, machine_key(df[[cols[["machine"]]]]), clean_text(df[[cols[["analysis"]]]]), sep = " | ")
  }
  data.frame(Database = if (is.na(period)) "PBDB" else "FAUNMAP", Dataset = dataset, Site_Key = key,
             SiteName = clean_text(df[[cols[["site"]]]]),
             Latitude = to_coord(df[[cols[["lat"]]]], 90), Longitude = to_coord(df[[cols[["lon"]]]], 180),
             df[spatial_cols], stringsAsFactors = FALSE)
}
all_sites <- bind_rows(
  c(lapply(periods, function(p) site_rows(loc_spatial[[p]], loc_cols[[p]], p, paste("FAUNMAP", p))),
    list(site_rows(pbdb_spatial, pbdb_cols, NA, "PBDB")))
) %>% distinct()

# Records per site (fauna records / PBDB occurrences) and their taxa.
taxon_of <- function(df, file) {
  g <- find_col(df, "genus", file, FALSE); s <- find_col(df, "species", file, FALSE)
  t <- find_col(df, "taxon", file, FALSE)
  if (!is.na(t)) return(clean_text(df[[t]]))
  if (!is.na(g)) return(clean_text(trimws(paste(df[[g]], if (!is.na(s)) coalesce(df[[s]], "") else ""))))
  rep(NA_character_, nrow(df))
}
records <- bind_rows(
  c(lapply(periods, function(p) {
      f <- fauna_spatial[[p]]
      if (nrow(f) == 0) return(NULL)
      data.frame(Dataset = paste("FAUNMAP", p),
                 Site_Key = paste(p, machine_key(f[[find_col(f, "machine", "")]]),
                                  clean_text(f[[find_col(f, "analysis", "")]]), sep = " | "),
                 Taxon = taxon_of(f, fauna_files[[p]]), f[spatial_cols], stringsAsFactors = FALSE)
    }),
    list(data.frame(Dataset = "PBDB",
                    Site_Key = paste("PBDB", clean_text(pbdb_spatial[[pbdb_cols[["site"]]]]), sep = " | "),
                    Taxon = taxon_of(pbdb_spatial, pbdb_file), pbdb_spatial[spatial_cols],
                    stringsAsFactors = FALSE)))
)
# Fauna use the locality's Analysis Unit text, which may differ in case from
# the locality file: count records by the normalised key.
norm_site <- function(k) tolower(gsub("\\s+", " ", k))
all_sites$n_records <- as.integer(table(norm_site(records$Site_Key))[norm_site(all_sites$Site_Key)])
all_sites$n_records[is.na(all_sites$n_records)] <- 0L

excluded_sites <- filter(all_sites, !is.na(Spatial_Exclusion))
sites_spatial  <- filter(all_sites, is.na(Spatial_Exclusion)) %>% select(-Spatial_Exclusion)
records        <- filter(records, is.na(Spatial_Exclusion))

cat("  Sites removed:\n")
print(as.data.frame(count(excluded_sites, Dataset, Spatial_Exclusion)), row.names = FALSE)

keep <- function(df) filter(df, is.na(Spatial_Exclusion)) %>% select(-Spatial_Exclusion)
loc_out   <- lapply(loc_spatial, keep)
fauna_out <- lapply(fauna_spatial, keep)
pbdb_out  <- keep(pbdb_spatial)

for (p in periods) {
  cat(sprintf("  %-13s %6d localities kept (%d removed) | %6d fauna records kept (%d removed)\n", p,
              nrow(loc_out[[p]]), nrow(loc_spatial[[p]]) - nrow(loc_out[[p]]),
              nrow(fauna_out[[p]]), nrow(fauna_spatial[[p]]) - nrow(fauna_out[[p]])))
}
cat(sprintf("  %-13s %6d occurrences kept (%d removed)\n", "PBDB", nrow(pbdb_out),
            nrow(pbdb_spatial) - nrow(pbdb_out)))

# Check against the State/prov column of FAUNMAP (USA and Canada).
check <- bind_rows(lapply(periods, function(p) {
  col <- names(loc_out[[p]])[simplify_name(names(loc_out[[p]])) %in% c("stateprov", "stateprovince", "state")]
  if (!length(col)) return(NULL)
  data.frame(File = toupper(clean_text(loc_out[[p]][[col[1]]])), Assigned = loc_out[[p]]$State_Code,
             Country = loc_out[[p]]$Country)
})) %>% filter(!is.na(File), !is.na(Assigned), Country %in% c("USA", "Canada"))
if (nrow(check)) {
  cat(sprintf("  FAUNMAP State/prov agrees with the assigned state for %d of %d localities (%.1f%%)\n",
              sum(check$File == check$Assigned), nrow(check), 100 * mean(check$File == check$Assigned)))
}

# =============================================================================
# 5. SAVE THE TABLES
# =============================================================================

cat("\n=== 5. SAVING ===\n")
for (p in periods) {
  save_out(loc_out[[p]], paste0(tolower(p), "_localities"))
  save_out(fauna_out[[p]], paste0(tolower(p), "_fauna"))
}
save_out(pbdb_out, "pbdb_clean")
save_out(sites_spatial, "sites_spatial", excel = TRUE)
save_out(excluded_sites, "excluded_sites")

# =============================================================================
# 6. ONE SUB-FOLDER PER LEVEL: SITE TABLE, UNIT SUMMARY AND MAP
# =============================================================================

frame_ll <- st_as_sfc(st_bbox(c(xmin = map_lon[1], xmax = map_lon[2],
                                ymin = map_lat[1], ymax = map_lat[2]), crs = 4326))
countries <- st_as_sf(maps::map("world", regions = c("Canada", "USA", "Mexico"), fill = TRUE, plot = FALSE))
st_crs(countries) <- NA; st_crs(countries) <- 4326   # relabel the maps package's old datum
countries <- suppressWarnings(suppressMessages(
  countries %>% st_make_valid() %>% st_set_agr("constant") %>% st_intersection(frame_ll) %>%
    st_transform(map_crs)))
lims <- st_bbox(countries)

# 60 degrees N limit line across Canada.
limit_line <- NULL
if (!is.na(canada_north_limit)) {
  limit_line <- st_sfc(st_linestring(cbind(seq(-141, -52, by = 0.5), canada_north_limit)), crs = 4326) %>%
    st_transform(map_crs)
}

# Area used to decide which units are drawn: inside the map, south of the limit,
# excluding Alaska and Hawaii.
study_area <- st_as_sfc(st_bbox(c(xmin = map_lon[1], xmax = map_lon[2], ymin = map_lat[1],
                                  ymax = min(map_lat[2], coalesce(canada_north_limit, 90))), crs = 4326)) %>%
  st_transform(work_crs)
drop_units <- c(if (remove_alaska) "Alaska", if (remove_hawaii) "Hawaii")

use_repel <- requireNamespace("ggrepel", quietly = TRUE)
site_col  <- c(political = "State_Province", division = "Physio_Division", province = "Physio_Province")
match_col <- c(political = "Political_Match", division = "Division_Match", province = "Province_Match")
country_order <- c("USA", "Canada", "Mexico")

# Numbered key drawn as a panel beside the map, grouped under headings.
key_panel <- function(key, title) {
  lines <- bind_rows(lapply(unique(key$Group), function(g) {
    k <- key[key$Group == g, ]
    bind_rows(data.frame(type = "head", text = g, num = NA, n = NA),
              data.frame(type = "item", text = k$Unit, num = k$Number, n = k$n_sites))
  }))
  n_col <- ceiling(nrow(lines) / 40)
  per_col <- ceiling(nrow(lines) / n_col)
  # Start each column at a heading where possible (no orphaned headings at the bottom).
  lines$col <- pmin(n_col, ceiling(seq_len(nrow(lines)) / per_col))
  for (i in which(lines$type == "head")) if (i < nrow(lines) && lines$col[i + 1] > lines$col[i]) lines$col[i] <- lines$col[i + 1]
  lines <- lines %>% group_by(col) %>% mutate(row = row_number()) %>% ungroup()
  max_row <- max(lines$row)
  lines$x <- (lines$col - 1) * 1
  lines$y <- -lines$row
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
    scale_y_continuous(limits = c(-(max(max_row, 30) + 0.8), 0), expand = c(0, 0)) +
    labs(title = title, subtitle = "Number on the map  -  name  (sites)") +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10, margin = margin(b = 6)))
}

for (lv in names(layers)) {
  cat("\n=== 6. LEVEL:", toupper(lv), "===\n")
  lv_dir <- file.path(output_dir, lv)
  dir.create(lv_dir, recursive = TRUE, showWarnings = FALSE)

  polys <- layers[[lv]] %>% filter(!Unit %in% drop_units)
  # Merge same-named units across countries (e.g. Interior Plains).
  polys <- suppressWarnings(suppressMessages(
    polys %>% group_by(Unit) %>%
      summarise(Countries = paste(intersect(country_order, unique(Country)), collapse = ", "),
                .groups = "drop") %>%
      polygons_only()))

  # Heading for the key: country (political, division) or parent division (province).
  if (lv == "province") {
    parents <- suppressWarnings(suppressMessages(
      polys %>% st_set_agr("constant") %>% st_cast("POLYGON", warn = FALSE) %>%
        mutate(area = as.numeric(st_area(.))) %>% group_by(Unit) %>%
        slice_max(area, n = 1, with_ties = FALSE) %>% ungroup() %>%
        st_set_agr("constant") %>% st_point_on_surface()))
    pd <- assign_units(parents, rbind(us_units$division, ca_units$division))
    polys$Group <- pd$Nearest_Unit[match(polys$Unit, parents$Unit)]
  } else {
    polys$Group <- sub(",.*", "", polys$Countries)
  }

  # Sites of this level.
  lv_sites <- sites_spatial %>%
    transmute(Database, Dataset, Site_Key, SiteName, Latitude, Longitude, Country,
              Unit = .data[[site_col[[lv]]]], Match = .data[[match_col[[lv]]]], n_records)
  lv_records <- records %>% transmute(Dataset, Unit = .data[[site_col[[lv]]]], Taxon)

  # Units drawn: those reaching the study area, plus any unit holding a site.
  in_study <- lengths(st_intersects(polys, study_area)) > 0
  polys <- polys[in_study | polys$Unit %in% lv_sites$Unit, ]

  counts <- lv_sites %>% filter(!is.na(Unit)) %>% count(Unit, name = "n_sites")
  group_levels <- if (lv == "province") {
    gk <- st_drop_geometry(polys) %>% distinct(Group, Countries) %>%
      mutate(o = match(sub(",.*", "", Countries), country_order)) %>%
      group_by(Group) %>% summarise(o = min(o), .groups = "drop") %>% arrange(o, Group)
    gk$Group
  } else intersect(country_order, polys$Group)
  key <- st_drop_geometry(polys) %>%
    left_join(counts, by = "Unit") %>%
    mutate(n_sites = coalesce(n_sites, 0L), Group = factor(Group, levels = group_levels)) %>%
    arrange(Group, Unit) %>%
    mutate(Number = row_number(), Group = as.character(Group))

  # --- Tables ---
  lv_sites$Unit_Number <- key$Number[match(lv_sites$Unit, key$Unit)]
  lv_sites <- relocate(lv_sites, Unit_Number, .before = Unit)
  names(lv_sites)[names(lv_sites) == "Unit"] <- site_col[[lv]]
  save_out(lv_sites, paste0("sites_", lv), lv_dir)

  # Sites per unit and dataset, one column per dataset.
  in_unit <- filter(lv_sites, !is.na(.data[[site_col[[lv]]]]))
  per_ds <- data.frame(Unit = sort(unique(in_unit[[site_col[[lv]]]])), stringsAsFactors = FALSE)
  for (d in names(dataset_colours)) {
    n_d <- table(in_unit[[site_col[[lv]]]][in_unit$Dataset == d])
    per_ds[[paste0("n_sites_", gsub(" ", "_", d))]] <- as.integer(coalesce(as.numeric(n_d[per_ds$Unit]), 0))
  }
  rec <- lv_records %>% filter(!is.na(Unit)) %>% group_by(Unit) %>%
    summarise(n_records = n(), n_taxa = n_distinct(Taxon, na.rm = TRUE), .groups = "drop")
  unit_summary <- key %>%
    transmute(Unit_Number = Number, Unit, Countries, Group, n_sites) %>%
    left_join(per_ds, by = "Unit") %>% left_join(rec, by = "Unit") %>%
    mutate(across(where(is.numeric), ~ coalesce(.x, 0L)))
  names(unit_summary)[names(unit_summary) == "Unit"] <- site_col[[lv]]
  if (lv != "province") unit_summary$Group <- NULL else names(unit_summary)[names(unit_summary) == "Group"] <- "Physio_Division"
  save_out(unit_summary, paste0("unit_summary_", lv), lv_dir, excel = TRUE)
  cat(sprintf("  %d units drawn, %d hold sites | %d sites without a unit at this level\n",
              nrow(key), sum(key$n_sites > 0), sum(is.na(lv_sites$Unit_Number))))

  # --- Map ---
  draw <- polys %>% st_transform(map_crs)
  if (simplify_m > 0) draw <- suppressWarnings(st_simplify(draw, dTolerance = simplify_m, preserveTopology = TRUE))
  draw <- draw %>% left_join(select(key, Unit, Number, n_sites), by = "Unit")

  # Touching occupied units get different tints; empty units are pale grey.
  touching <- st_is_within_distance(draw, draw, dist = 20000)
  tint <- integer(nrow(draw))
  for (i in order(lengths(touching), decreasing = TRUE)) {
    used <- tint[setdiff(touching[[i]], i)]
    tint[i] <- which(!seq_along(unit_tints) %in% used)[1]
    if (is.na(tint[i])) tint[i] <- 1L
  }
  draw$Fill <- ifelse(draw$n_sites > 0, unit_tints[tint], empty_fill)

  label_pts <- suppressWarnings(suppressMessages(
    draw %>% select(Number, n_sites) %>% st_set_agr("constant") %>%
      st_intersection(st_transform(st_set_agr(st_sf(geometry = study_area), "constant"), map_crs)) %>%
      st_cast("POLYGON", warn = FALSE) %>%
      mutate(area = as.numeric(st_area(.))) %>% group_by(Number) %>%
      slice_max(area, n = 1, with_ties = FALSE) %>% ungroup() %>%
      st_set_agr("constant") %>% st_point_on_surface()))
  label_xy <- cbind(st_drop_geometry(label_pts)[c("Number", "n_sites")], st_coordinates(label_pts))

  pts <- lv_sites %>% filter(!is.na(Latitude), !is.na(Longitude)) %>%
    distinct(Latitude, Longitude) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% st_transform(map_crs)

  p_map <- ggplot() +
    geom_sf(data = countries, fill = no_data_fill, colour = NA) +
    geom_sf(data = draw, fill = draw$Fill, colour = "grey40", linewidth = 0.25) +
    geom_sf(data = countries, fill = NA, colour = "grey15", linewidth = 0.45)
  if (!is.null(limit_line)) {
    p_map <- p_map + geom_sf(data = limit_line, colour = "grey25", linewidth = 0.4, linetype = "22")
  }
  p_map <- p_map +
    geom_sf(data = pts, shape = 21, fill = site_fill, colour = "white", stroke = 0.3,
            size = point_size, alpha = 0.9)
  if (use_repel) {
    p_map <- p_map + ggrepel::geom_label_repel(
      data = label_xy, aes(x = X, y = Y, label = Number), size = 3, fontface = "bold",
      colour = ifelse(label_xy$n_sites > 0, "grey5", "grey50"), fill = alpha("white", 0.85),
      label.size = 0.15, label.padding = unit(0.1, "lines"), min.segment.length = 0.2,
      segment.colour = "grey25", segment.size = 0.25, max.overlaps = Inf, seed = 1,
      box.padding = 0.15)
  } else {
    p_map <- p_map + geom_label(data = label_xy, aes(x = X, y = Y, label = Number), size = 3,
                                fontface = "bold", fill = alpha("white", 0.85), label.size = 0.15)
  }
  caption <- paste0("Sites in Alaska, Hawaii",
                    if (!is.na(canada_north_limit)) paste0(" and Canada north of ", canada_north_limit, "\u00b0N (dashed line)"),
                    " removed. Pale units hold no sites.",
                    if (lv != "political") " Mexico: no physiographic layer (grey)." else "")
  p_map <- p_map +
    coord_sf(crs = map_crs, xlim = lims[c("xmin", "xmax")], ylim = lims[c("ymin", "ymax")], expand = FALSE) +
    labs(title = level_titles[[lv]],
         subtitle = sprintf("%d fossil sites in %d of %d units",
                            sum(!is.na(lv_sites$Unit_Number)), sum(key$n_sites > 0), nrow(key)),
         caption = caption, x = NULL, y = NULL) +
    guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1), nrow = 1),
           shape = guide_legend(nrow = 1)) +
    theme_void(base_size = 13) +
    theme(panel.background = element_rect(fill = sea_fill, colour = "grey30", linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 17),
          plot.subtitle = element_text(colour = "grey30", size = 12, margin = margin(b = 6)),
          plot.caption = element_text(colour = "grey35", size = 9.5, hjust = 0),
          legend.position = "bottom", legend.text = element_text(size = 11),
          legend.margin = margin(t = 4))

  key_title <- c(political = "States and provinces", division = "Physiographic units",
                 province = "Provinces / subregions by division")[[lv]]
  p <- patchwork::wrap_plots(p_map, key_panel(key, key_title), widths = c(1.45, 1)) &
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(10, 12, 8, 12))

  for (ext in c("png", "pdf")) {
    file <- file.path(lv_dir, paste0("map_", lv, ".", ext))
    if (ext == "png") ggsave(file, p, width = map_width, height = map_height, dpi = map_dpi, bg = "white")
    else ggsave(file, p, width = map_width, height = map_height, bg = "white")
    cat("  ", file, "\n", sep = "")
  }
}

cat("\n=== STEP 3 COMPLETE ===\n")
cat("Objects in your Environment: sites_spatial, excluded_sites, loc_out, fauna_out, pbdb_out\n")
cat("Next: Step 3c (time binning) reads the files in Outputs/3_spatial.\n")
