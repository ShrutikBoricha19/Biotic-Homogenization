# =============================================================================
# STEP 5: MAP ALL SITES ON THE PHYSIOGRAPHIC REGIONS OF NORTH AMERICA
#
# Same five stage maps as Step 4 (one point per site, FAUNMAP and PBDB not
# distinguished), drawn on the physiographic regions of:
#
#   USA     Fenneman & Johnson (1946), "Physiographic divisions of the
#           conterminous U.S." (USGS physio_shp) - your local copy.
#           Covers the conterminous US only (no Alaska / Hawaii).
#   Canada  Natural Resources Canada, "Physiographic Regions of Canada"
#           (open.canada.ca record a3dfbaf4-1b20-4061-aa0a-e7a79953f52d).
#           Downloaded automatically from NRCan's ArcGIS service.
#   Mexico  CONABIO, "Provincias Fisiograficas de Mexico", 1:4,000,000
#           (Cervantes-Zamora et al. 1990, after INEGI).
#           Downloaded automatically from CONABIO.
#
# Downloads are cached in Outputs/maps/basemap/ and reused. If a download
# fails (no internet, site down), download the data by hand and set
# canada_path / mexico_path below. Any country whose layer is missing is
# drawn plain grey; the maps are still produced.
#
# Reading the maps:
#   - Each physiographic unit has a NUMBER on the map, listed in the key.
#   - Fill colours only separate neighbouring units (touching units never
#     share a colour); they do not identify units by themselves.
#   - Grey land = no physiographic layer (e.g. Alaska, Hawaii).
#
# Also saved: the physiographic unit of every site (site_physio_units.csv).
#
# Input:  Outputs/3_stages/site_index.csv  (from Step 3)
# Output: Outputs/maps/physiographic/
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 3.
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

site_index_file <- file.path("Outputs", "3_stages", "site_index.csv")
output_dir      <- file.path(work_dir, "Outputs", "maps", "physiographic")
basemap_dir     <- file.path(work_dir, "Outputs", "maps", "basemap")   # download cache

# --- USA (local shapefile) ---------------------------------------------------
# A folder containing physio.shp, or the .shp file itself. Use forward
# slashes "/" in Windows paths.
us_path  <- "C:/Users/shrut/Downloads/physio_shp"
us_field <- "DIVISION"      # level to map: "DIVISION" (8), "PROVINCE" (25) or "SECTION" (86)

# --- Canada -------------------------------------------------------------------
# NULL = download from NRCan. Or a local file / folder / .gdb.
canada_path        <- NULL
canada_local_layer <- NULL  # layer name inside a local .gdb (NULL = choose automatically)
canada_service     <- "https://maps-cartes.services.geo.ca/server_serveur/rest/services/NRCan/phys_reg_en/MapServer"
canada_layer       <- 0     # 0 = Regions (7), 1 = Subregions (21), 2 = Divisions
canada_field       <- NULL  # NULL = detect the name column automatically

# --- Mexico -------------------------------------------------------------------
# NULL = download from CONABIO. Or a local folder / .shp.
mexico_path  <- NULL
mexico_url   <- "http://www.conabio.gob.mx/informacion/gis/maps/geo/rfisio4mgw.zip"
mexico_field <- NULL        # NULL = detect the province name column automatically

# --- Map appearance -------------------------------------------------------------
map_crs   <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
map_lon   <- c(-170, -50)   # extent drawn (degrees)
map_lat   <- c(14, 84)
simplify_m <- 1000          # smooth region outlines to ~1 km (0 = no smoothing)
show_state_lines <- TRUE    # state / province lines from Step 4's cached file

point_fill   <- "#B2182B"
point_size   <- 2
region_tints <- c("#E9DDB9", "#CFE3C6", "#C8DCEC", "#DCD2EA", "#D2E9E2", "#F1E3C4")
no_data_fill <- "grey88"
sea_fill     <- "#F3F7FA"

snap_km <- 10               # sites just off a unit's edge (coast) are assigned to the
                            # nearest unit within this distance

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

# Read a spatial file from a path that may be a .shp, a folder holding
# shapefiles, a .gdb, or a .geojson. `hint` picks among several files/layers.
read_spatial <- function(path, hint = NULL, layer = NULL) {
  if (!file.exists(path)) stop("Not found: ", path)

  if (dir.exists(path) && !grepl("\\.gdb/?$", path, ignore.case = TRUE)) {
    shp <- list.files(path, pattern = "\\.shp$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
    if (length(shp) == 0) stop("No .shp file found in folder: ", path)
    if (length(shp) > 1) {
      pick <- if (!is.null(hint)) grep(hint, basename(shp), ignore.case = TRUE) else integer(0)
      cat("    Several shapefiles found:", paste(basename(shp), collapse = ", "), "\n")
      shp <- shp[if (length(pick)) pick[1] else 1]
    }
    cat("    Reading", shp, "\n")
    return(st_read(shp, quiet = TRUE))
  }

  if (grepl("\\.gdb/?$", path, ignore.case = TRUE) && is.null(layer)) {
    lyr <- st_layers(path)$name
    cat("    Layers in geodatabase:", paste(lyr, collapse = ", "), "\n")
    pick <- if (!is.null(hint)) grep(hint, lyr, ignore.case = TRUE) else integer(0)
    layer <- lyr[if (length(pick)) pick[1] else 1]
  }

  cat("    Reading", path, if (!is.null(layer)) paste0("(layer: ", layer, ")"), "\n")
  if (is.null(layer)) st_read(path, quiet = TRUE) else st_read(path, layer = layer, quiet = TRUE)
}

# Repair text read with the wrong encoding (common in Spanish shapefiles).
fix_text <- function(x) {
  x <- as.character(x)
  bad <- !is.na(x) & !validUTF8(x)
  x[bad] <- iconv(x[bad], "latin1", "UTF-8")
  Encoding(x) <- "UTF-8"
  x
}

# "INTERIOR PLAINS" -> "Interior Plains"; other names are left as written.
tidy_name <- function(x) {
  x <- gsub("\\s+", " ", trimws(fix_text(as.character(x))))
  upper <- !is.na(x) & x == toupper(x) & grepl("[A-Z]", x)
  if (any(upper)) x[upper] <- tools::toTitleCase(tolower(x[upper]))
  x[x == ""] <- NA
  x
}

# Print every attribute column with its number of distinct values and an
# example, so the right name column can be chosen if detection is wrong.
describe_fields <- function(df) {
  att <- st_drop_geometry(df)
  for (col in names(att)) {
    v <- att[[col]]
    ex <- v[!is.na(v)][1]
    cat(sprintf("      %-15s %-10s %4d distinct   e.g. %s\n", col, class(v)[1],
                n_distinct(v, na.rm = TRUE), substr(fix_text(as.character(ex)), 1, 40)))
  }
}

pick_field <- function(df, preferred, override, label) {
  cols <- names(st_drop_geometry(df))
  if (!is.null(override)) {
    hit <- cols[toupper(cols) == toupper(override)]
    if (length(hit) == 0) {
      stop(label, ": column '", override, "' not found. Columns: ", paste(cols, collapse = ", "))
    }
    return(hit[1])
  }
  for (p in preferred) {
    hit <- cols[toupper(cols) == toupper(p)]
    if (length(hit)) return(hit[1])
  }
  # Fallback: first text column that looks like a set of names.
  att <- st_drop_geometry(df)
  for (col in cols) {
    v <- att[[col]]
    if (is.character(v) && n_distinct(v, na.rm = TRUE) >= 2 &&
        n_distinct(v, na.rm = TRUE) <= 300) return(col)
  }
  stop(label, ": could not find a name column. Set it in SETTINGS. Columns: ",
       paste(cols, collapse = ", "))
}

# Download one layer of an ArcGIS REST MapServer as GeoJSON (with paging).
arcgis_download <- function(service, layer, dest) {
  base <- paste0(service, "/", layer, "/query?where=1%3D1&outFields=*",
                 "&returnGeometry=true&outSR=4326",
                 "&maxAllowableOffset=0.005&geometryPrecision=5&f=geojson")
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
  if (length(pages) == 0) stop("The service returned no features.")

  out <- do.call(rbind, pages)
  st_write(out, dest, quiet = TRUE, delete_dsn = TRUE)
  out
}

# Keep only polygon parts, as MULTIPOLYGON (repairs can leave stray lines).
polygons_only <- function(x) {
  if (any(st_geometry_type(x) == "GEOMETRYCOLLECTION")) {
    x <- st_collection_extract(x, "POLYGON", warn = FALSE)
  }
  x <- x[st_geometry_type(x) %in% c("POLYGON", "MULTIPOLYGON"), ]
  st_cast(x, "MULTIPOLYGON", warn = FALSE)
}

# Bring a layer into the map projection with one row per named unit.
standardise <- function(layer, field, country) {
  layer %>%
    st_make_valid() %>%
    st_transform(map_crs) %>%
    transmute(Country = country, Unit = tidy_name(.data[[field]])) %>%
    filter(!is.na(Unit)) %>%
    polygons_only() %>%
    group_by(Country, Unit) %>%
    summarise(.groups = "drop") %>%          # dissolve pieces of each unit
    st_make_valid() %>%
    polygons_only()
}

# =============================================================================
# 1. SITES: ONE POINT PER SITE PER STAGE (same rule as Step 4)
# =============================================================================

cat("=== 1. READING SITES ===\n")

path <- file.path(work_dir, site_index_file)
if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun Step 3 first.")

site_index <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))

sites <- site_index %>%
  mutate(
    Latitude = suppressWarnings(as.numeric(Latitude)),
    Longitude = suppressWarnings(as.numeric(Longitude)),
    SiteName = trimws(SiteName),
    Site = if_else(is.na(SiteName) | SiteName == "", Site_Key, SiteName),
    Has_Coordinates = !is.na(Latitude) & !is.na(Longitude)
  ) %>%
  arrange(Stage_Number, Database, Site, desc(Has_Coordinates)) %>%
  group_by(Stage_Number, Database, Site) %>%
  summarise(Latitude = first(Latitude), Longitude = first(Longitude),
            Has_Coordinates = first(Has_Coordinates), .groups = "drop") %>%
  mutate(Stage = stage_names[Stage_Number])

cat(sprintf("  %d site-stage points, %d without coordinates (not mappable)\n",
            nrow(sites), sum(!sites$Has_Coordinates)))

sites_sf <- sites %>%
  filter(Has_Coordinates) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(map_crs)

# =============================================================================
# 2. PHYSIOGRAPHIC LAYERS
# =============================================================================

cat("\n=== 2. LOADING PHYSIOGRAPHIC LAYERS ===\n")

layers <- list()

# ---- USA -----------------------------------------------------------------------
cat("  USA (Fenneman & Johnson 1946)\n")
us <- tryCatch(read_spatial(us_path, hint = "physio"), error = function(e) {
  cat("    NOT LOADED:", conditionMessage(e), "\n"); NULL
})

if (!is.null(us)) {
  if (is.na(st_crs(us))) {
    # physio.shp is often distributed without a .prj file. Its documented
    # projection is Albers Equal Area (29.5, 45.5, origin 23N 96W), NAD27.
    bb <- st_bbox(us)
    if (all(abs(bb[c("xmin", "xmax")]) <= 180) && all(abs(bb[c("ymin", "ymax")]) <= 90)) {
      st_crs(us) <- 4269
      cat("    No .prj file: coordinates look like degrees -> NAD83 geographic assumed\n")
    } else {
      st_crs(us) <- paste("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96",
                          "+x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs")
      cat("    No .prj file: USGS Albers Equal Area (NAD27) assumed, as documented\n")
    }
  }
  # Sanity check: the conterminous US should land at about 125-66 W, 24-50 N.
  bb <- st_bbox(st_transform(us, 4326))
  cat(sprintf("    Extent: %.1f to %.1f lon, %.1f to %.1f lat\n",
              bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]))
  if (bb["xmin"] < -130 || bb["xmax"] > -60 || bb["ymin"] < 20 || bb["ymax"] > 53) {
    warning("The US layer is not where the conterminous US should be. ",
            "Its coordinate system is probably wrong - check st_crs().")
  }
  describe_fields(us)
  f <- pick_field(us, us_field, us_field, "USA")
  layers$USA <- standardise(us, f, "USA")
  cat(sprintf("    Using '%s': %d units\n", f, nrow(layers$USA)))
}

# ---- Canada --------------------------------------------------------------------
cat("  Canada (NRCan Physiographic Regions)\n")
canada <- tryCatch({
  if (!is.null(canada_path)) {
    read_spatial(canada_path, hint = c("region", "subregion", "division")[canada_layer + 1],
                 layer = canada_local_layer)
  } else {
    cache <- file.path(basemap_dir, paste0("canada_physio_layer", canada_layer, ".geojson"))
    if (!file.exists(cache)) {
      cat("    Downloading from NRCan (one time)...\n")
      arcgis_download(canada_service, canada_layer, cache)
    }
    cat("    Reading", cache, "\n")
    st_read(cache, quiet = TRUE)
  }
}, error = function(e) {
  cat("    NOT LOADED:", conditionMessage(e), "\n",
      "   Download 'Physiographic Regions of Canada' from open.canada.ca and set canada_path.\n")
  NULL
})

if (!is.null(canada)) {
  describe_fields(canada)
  f <- pick_field(canada,
                  c("REGION_EN", "REGION", "SUBREGION_EN", "SUBREGION", "DIVISION_EN",
                    "DIVISION", "NAME_EN", "ENGLISH_NAME", "NAME_E", "NAME"),
                  canada_field, "Canada")
  layers$Canada <- standardise(canada, f, "Canada")
  cat(sprintf("    Using '%s': %d units\n", f, nrow(layers$Canada)))
}

# ---- Mexico --------------------------------------------------------------------
cat("  Mexico (CONABIO Provincias Fisiograficas)\n")
mexico <- tryCatch({
  src <- mexico_path
  if (is.null(src)) {
    src <- file.path(basemap_dir, "mexico_rfisio4mgw")
    if (!dir.exists(src) || length(list.files(src, pattern = "\\.shp$", recursive = TRUE)) == 0) {
      cat("    Downloading from CONABIO (one time)...\n")
      zip <- file.path(basemap_dir, basename(mexico_url))
      old <- options(timeout = 600)
      on.exit(options(old), add = TRUE)
      ok <- tryCatch({
        suppressWarnings(download.file(mexico_url, zip, mode = "wb", quiet = TRUE))
        TRUE
      }, error = function(e) FALSE)
      if (!ok || !file.exists(zip) || file.size(zip) < 1000) {
        if (file.exists(zip)) file.remove(zip)   # so the next run tries again
        stop("download from CONABIO failed")
      }
      dir.create(src, showWarnings = FALSE)
      unzip(zip, exdir = src)
    }
  }
  read_spatial(src, hint = "fisio")
}, error = function(e) {
  cat("    NOT LOADED:", conditionMessage(e), "\n",
      "   Download 'Provincias Fisiograficas de Mexico' (rfisio4mgw) from CONABIO and set mexico_path.\n")
  NULL
})

if (!is.null(mexico)) {
  if (is.na(st_crs(mexico))) {
    st_crs(mexico) <- 4326
    cat("    No .prj file: WGS84 geographic assumed\n")
  }
  describe_fields(mexico)
  f <- pick_field(mexico,
                  c("PROVINCIA", "PROVINCIAS", "NOM_PROV", "NOMPROV", "PROV_FIS",
                    "PROVFIS", "NOMBRE", "NOM"),
                  mexico_field, "Mexico")
  layers$Mexico <- standardise(mexico, f, "Mexico")
  cat(sprintf("    Using '%s': %d units\n", f, nrow(layers$Mexico)))
}

if (length(layers) == 0) stop("No physiographic layer could be loaded.")

physio <- do.call(rbind, layers)
if (simplify_m > 0) physio <- st_simplify(physio, dTolerance = simplify_m, preserveTopology = TRUE)

# =============================================================================
# 3. NUMBER THE UNITS AND CHOOSE FILL COLOURS
#
# Units with the same name in different countries (e.g. Interior Plains in
# the US and Canada) share one number.
# =============================================================================

cat("\n=== 3. PREPARING THE KEY ===\n")

country_order <- c("USA", "Canada", "Mexico")
country_code  <- c(USA = "US", Canada = "CA", Mexico = "MX")

key <- st_drop_geometry(physio) %>%
  group_by(Unit) %>%
  summarise(
    first_country = min(match(Country, country_order)),
    Countries = paste(country_code[country_order[sort(unique(match(Country, country_order)))]],
                      collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(first_country, Unit) %>%
  mutate(Number = row_number(),
         Key_Label = sprintf("%2d  %s (%s)", Number, Unit, Countries)) %>%
  select(Number, Unit, Countries, Key_Label)

physio <- physio %>%
  group_by(Unit) %>%
  summarise(.groups = "drop") %>%            # merge same-named units across borders
  polygons_only() %>%
  left_join(key, by = "Unit")

# Map colouring: give touching units different tints (greedy, most
# connected unit first). Units within 20 km count as touching so that small
# gaps between the national datasets are ignored.
touching <- st_is_within_distance(physio, physio, dist = 20000)
tint <- integer(nrow(physio))
for (i in order(lengths(touching), decreasing = TRUE)) {
  used <- tint[setdiff(touching[[i]], i)]
  tint[i] <- which(!seq_along(region_tints) %in% used)[1]
  if (is.na(tint[i])) tint[i] <- 1L
}
physio$Tint <- factor(tint)

# Label position: a point inside the largest piece of each unit.
label_pts <- physio %>%
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

for (i in seq_len(nrow(key))) cat(" ", key$Key_Label[i], "\n")
write.csv(select(key, Number, Unit, Countries),
          file.path(output_dir, "physio_units_key.csv"), row.names = FALSE)

# =============================================================================
# 4. PHYSIOGRAPHIC UNIT OF EVERY SITE
# =============================================================================

cat("\n=== 4. ASSIGNING SITES TO PHYSIOGRAPHIC UNITS ===\n")

inside <- st_intersects(sites_sf, physio)
unit_idx <- vapply(inside, function(h) if (length(h)) h[1] else NA_integer_, integer(1))
match_type <- ifelse(is.na(unit_idx), "None", "Inside")

# Sites just outside every unit (usually on the coast): nearest unit within snap_km.
out <- which(is.na(unit_idx))
if (length(out) > 0) {
  near <- st_nearest_feature(sites_sf[out, ], physio)
  d_km <- as.numeric(st_distance(sites_sf[out, ], physio[near, ], by_element = TRUE)) / 1000
  ok <- d_km <= snap_km
  unit_idx[out[ok]] <- near[ok]
  match_type[out[ok]] <- sprintf("Nearest (%.1f km)", d_km[ok])
}

site_units <- sites_sf %>%
  st_drop_geometry() %>%
  mutate(Physio_Number = physio$Number[unit_idx],
         Physio_Unit = physio$Unit[unit_idx],
         Physio_Match = match_type) %>%
  select(Stage_Number, Stage, Database, Site, Latitude, Longitude,
         Physio_Number, Physio_Unit, Physio_Match)

cat(sprintf("  %d inside a unit | %d assigned to the nearest unit (<= %g km) | %d outside all units\n",
            sum(match_type == "Inside"), sum(grepl("^Nearest", match_type)),
            snap_km, sum(match_type == "None")))
if (any(match_type == "None")) {
  cat("  (sites outside all units are usually in areas with no layer, e.g. Alaska,\n",
      "   or in a country whose layer did not load)\n", sep = "")
}

write.csv(site_units, file.path(output_dir, "site_physio_units.csv"),
          row.names = FALSE, na = "")

# =============================================================================
# 5. BASE LAYERS: COUNTRY OUTLINES AND STATE / PROVINCE LINES
# =============================================================================

frame_ll <- st_as_sfc(st_bbox(c(xmin = map_lon[1], xmax = map_lon[2],
                                ymin = map_lat[1], ymax = map_lat[2]), crs = 4326))

countries <- st_as_sf(maps::map("world", regions = c("Canada", "USA", "Mexico"),
                                fill = TRUE, plot = FALSE))
st_crs(countries) <- NA       # the maps package labels its lon/lat data with an old
st_crs(countries) <- 4326     # ellipsoid; relabel it so it lines up with WGS84
countries <- countries %>%
  st_make_valid() %>%
  st_set_agr("constant")
countries <- suppressMessages(st_intersection(countries, frame_ll)) %>%
  st_transform(map_crs)

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

# =============================================================================
# 6. DRAW THE FIVE MAPS
# =============================================================================

cat("\n=== 5. DRAWING MAPS ===\n")

use_repel <- requireNamespace("ggrepel", quietly = TRUE)

for (s in 1:5) {
  pts <- filter(sites_sf, Stage_Number == s)

  p <- ggplot() +
    geom_sf(data = countries, fill = no_data_fill, colour = NA) +
    geom_sf(data = physio, aes(fill = Tint), colour = "grey35", linewidth = 0.2) +
    scale_fill_manual(values = region_tints, guide = "none")

  if (!is.null(state_lines)) {
    p <- p + geom_sf(data = state_lines, colour = "grey30", linewidth = 0.1, alpha = 0.35)
  }

  p <- p +
    geom_sf(data = countries, fill = NA, colour = "grey20", linewidth = 0.3)

  # Unit numbers (with a white halo when ggrepel is available).
  if (use_repel) {
    p <- p + ggrepel::geom_text_repel(
      data = label_xy, aes(x = X, y = Y, label = Number),
      size = 2.4, fontface = "bold", colour = "grey15",
      bg.color = "white", bg.r = 0.15, min.segment.length = 0.3,
      segment.colour = "grey30", segment.size = 0.2, seed = 1
    )
  } else {
    p <- p + geom_text(data = label_xy, aes(x = X, y = Y, label = Number),
                       size = 2.4, fontface = "bold", colour = "grey15")
  }

  # Invisible layer that turns the numbered key into the legend.
  p <- p +
    geom_point(data = transform(label_xy, Key = key$Key_Label[match(Number, key$Number)]),
               aes(x = X, y = Y, shape = Key), alpha = 0) +
    scale_shape_manual(values = rep(32, nrow(key)), breaks = key$Key_Label,
                       name = "Physiographic units") +
    guides(shape = guide_legend(ncol = 1, override.aes = list(alpha = 0)))

  # Sites on top: solid points with a thin white ring.
  p <- p +
    geom_sf(data = pts, shape = 21, fill = point_fill, colour = "white",
            stroke = 0.3, size = point_size) +
    coord_sf(crs = map_crs, xlim = lims[c("xmin", "xmax")], ylim = lims[c("ymin", "ymax")],
             expand = FALSE) +
    labs(title = stage_titles[s],
         subtitle = paste0(nrow(pts), " sites (FAUNMAP + PBDB)"),
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      panel.background = element_rect(fill = sea_fill),
      panel.grid = element_line(colour = "white", linewidth = 0.3),
      plot.title = element_text(face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 7, family = "mono"),
      legend.title = element_text(size = 8, face = "bold"),
      legend.key.size = unit(0.32, "lines")
    )

  file <- file.path(output_dir, sprintf("physio_map_stage%d_%s.png", s, stage_names[s]))
  ggsave(file, p, width = 11, height = 7.5, dpi = 300)
  cat(sprintf("  %-12s %5d sites -> %s\n", stage_names[s], nrow(pts), file))
}

cat("\n  Unit key           -> ", file.path(output_dir, "physio_units_key.csv"), "\n", sep = "")
cat("  Site units         -> ", file.path(output_dir, "site_physio_units.csv"), "\n", sep = "")

cat("\n=== STEP 5 COMPLETE ===\n")
cat("Layers drawn:", paste(names(layers), collapse = ", "), "\n")
cat("Objects in your Environment: physio, key, site_units, sites_sf\n")
