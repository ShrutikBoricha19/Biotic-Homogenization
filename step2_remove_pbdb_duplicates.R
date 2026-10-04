# =============================================================================
# STEP 2: REMOVE PBDB COLLECTIONS THAT DUPLICATE FAUNMAP LOCALITIES
#
# A PBDB collection is treated as a duplicate of a FAUNMAP locality if,
# in this order of priority:
#   1. its collection_name matches a FAUNMAP SiteName (ignoring case and
#      extra spaces)
#   2. its coordinates equal a FAUNMAP locality's (rounded to 6 decimals)
#   3. it lies within `duplicate_km` of a FAUNMAP locality
#
# Every occurrence (faunal record) of a duplicate collection is removed.
# All FAUNMAP Blancan + Irvingtonian localities are used as the reference.
# No five-fauna-per-site filter is applied.
#
# Run the whole file (Ctrl+Shift+S in RStudio). Step 1 does not need to be
# run first; this step only uses the PBDB file and the two locality files.
# =============================================================================

library(dplyr)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

pbdb_file             <- "PBDB_Genus - Species - CAN_MEX_USA.csv"
blancan_loc_file      <- "Blancan Localities Data (Updated).csv"
irvingtonian_loc_file <- "Irvingtonian Localities Data (Updated).csv"

duplicate_km <- 5     # PBDB collections this close to a FAUNMAP locality are removed

output_dir <- file.path(work_dir, "Outputs", "2_pbdb")

# =============================================================================
# 1. CHECK THAT THE FOLDER AND FILES EXIST
# =============================================================================

cat("=== 1. CHECKING INPUT FILES ===\n")

if (!dir.exists(work_dir)) {
  stop("The folder does not exist:\n  ", work_dir)
}

needed <- c(pbdb_file, blancan_loc_file, irvingtonian_loc_file)
found <- file.exists(file.path(work_dir, needed))

for (i in seq_along(needed)) {
  cat(if (found[i]) "  [found]   " else "  [MISSING] ", needed[i], "\n", sep = "")
}

if (!all(found)) {
  cat("\nCSV files actually present in the folder:\n")
  cat(paste0("  ", list.files(work_dir, pattern = "\\.csv$", ignore.case = TRUE)),
      sep = "\n")
  stop("Some input files were not found. Fix the file names in SETTINGS.")
}

# =============================================================================
# 2. READ THE FILES
# =============================================================================

cat("\n=== 2. READING FILES ===\n")

read_input <- function(file, header_word = NULL) {
  path <- file.path(work_dir, file)

  # PBDB downloads can start with several lines of metadata before the real
  # column header. If so, skip down to the line containing `header_word`.
  skip <- 0
  if (!is.null(header_word)) {
    first_lines <- readLines(path, n = 100, warn = FALSE)
    header_line <- which(grepl(header_word, first_lines, useBytes = TRUE))[1]
    if (is.na(header_line)) {
      stop("Could not find a '", header_word, "' column in ", file)
    }
    skip <- header_line - 1
    if (skip > 0) cat("  (skipping", skip, "metadata lines at the top of", file, ")\n")
  }

  df <- read.csv(
    path,
    skip = skip,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character",          # keep identifiers exactly as written
    na.strings = c("", "NA", "N/A")
  )

  # Remove hidden Excel byte-order marks and stray spaces from column names.
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))

  cat(sprintf("  %-45s %7d rows, %3d columns\n", file, nrow(df), ncol(df)))
  df
}

pbdb_raw             <- read_input(pbdb_file, header_word = "collection_name")
blancan_loc_raw      <- read_input(blancan_loc_file)
irvingtonian_loc_raw <- read_input(irvingtonian_loc_file)

# =============================================================================
# 3. FIND THE COLUMNS WE NEED IN EACH FILE
#
# Same rules as Step 1: case, spaces, dots and underscores are ignored and
# the first option listed wins when several are present (LATDD > Latitude).
# =============================================================================

cat("\n=== 3. IDENTIFYING COLUMNS ===\n")

simplify_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

column_options <- list(
  site = c("sitename", "collectionname"),
  lat  = c("latdd", "latitude", "lat"),
  lon  = c("longdd", "longitude", "lng", "long", "lon")
)

find_col <- function(df, field, file) {
  simple <- simplify_name(names(df))

  for (option in column_options[[field]]) {
    hit <- names(df)[simple == option]
    if (length(hit) > 1) {
      stop("Several columns in ", file, " are all called '", option, "': ",
           paste(hit, collapse = ", "))
    }
    if (length(hit) == 1) return(hit)
  }

  stop("No '", field, "' column found in ", file,
       ".\nColumns in this file: ", paste(names(df), collapse = " | "))
}

show_cols <- function(df, file) {
  cols <- sapply(c("site", "lat", "lon"), function(f) find_col(df, f, file))
  cat("  ", file, "\n", sep = "")
  for (f in names(cols)) cat(sprintf("      %-5s -> \"%s\"\n", f, cols[[f]]))
  cols
}

pbdb_cols             <- show_cols(pbdb_raw, pbdb_file)
blancan_loc_cols      <- show_cols(blancan_loc_raw, blancan_loc_file)
irvingtonian_loc_cols <- show_cols(irvingtonian_loc_raw, irvingtonian_loc_file)

# =============================================================================
# 4. PREPARE FAUNMAP REFERENCE LOCALITIES AND PBDB COLLECTIONS
# =============================================================================

cat("\n=== 4. PREPARING LOCALITIES AND COLLECTIONS ===\n")

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

# Lower case with single spaces, used only for comparing names.
name_key <- function(x) gsub("\\s+", " ", tolower(clean_text(x)))

to_coord <- function(x, limit) {
  v <- suppressWarnings(as.numeric(clean_text(x)))
  v[!is.na(v) & abs(v) > limit] <- NA     # impossible values -> missing
  v
}

prepare_faunmap <- function(df, cols, period) {
  data.frame(
    FAUNMAP_SiteName = clean_text(df[[cols[["site"]]]]),
    FAUNMAP_Period = period,
    Latitude = to_coord(df[[cols[["lat"]]]], 90),
    Longitude = to_coord(df[[cols[["lon"]]]], 180),
    stringsAsFactors = FALSE
  )
}

faunmap_sites <- bind_rows(
  prepare_faunmap(blancan_loc_raw, blancan_loc_cols, "Blancan"),
  prepare_faunmap(irvingtonian_loc_raw, irvingtonian_loc_cols, "Irvingtonian")
) %>%
  distinct()

faunmap_coords <- filter(faunmap_sites, !is.na(Latitude), !is.na(Longitude))
faunmap_coord_keys <- paste(round(faunmap_coords$Latitude, 6),
                            round(faunmap_coords$Longitude, 6))

cat(sprintf("  FAUNMAP localities: %d distinct | %d distinct names | %d with coordinates\n",
            nrow(faunmap_sites),
            n_distinct(name_key(faunmap_sites$FAUNMAP_SiteName), na.rm = TRUE),
            nrow(faunmap_coords)))

# PBDB: each collection is identified by its collection_name.
# Occurrences with no collection_name are kept as their own "collection"
# so they can still be checked by coordinates.
pbdb <- pbdb_raw %>%
  mutate(
    PBDB_Row = row_number(),
    Collection = clean_text(.data[[pbdb_cols[["site"]]]]),
    Collection = if_else(is.na(Collection),
                         paste0("(no collection_name) row ", PBDB_Row),
                         Collection),
    Latitude = to_coord(.data[[pbdb_cols[["lat"]]]], 90),
    Longitude = to_coord(.data[[pbdb_cols[["lon"]]]], 180)
  )

n_unnamed <- sum(grepl("^\\(no collection_name\\)", pbdb$Collection))

# One row per distinct collection + coordinate combination.
pbdb_locations <- distinct(pbdb, Collection, Latitude, Longitude)

n_multi <- sum(duplicated(pbdb_locations$Collection))

cat(sprintf("  PBDB occurrences: %d | collections: %d | locations checked: %d\n",
            nrow(pbdb), n_distinct(pbdb$Collection), nrow(pbdb_locations)))
cat(sprintf("  PBDB occurrences without a collection_name: %d\n", n_unnamed))
cat(sprintf("  PBDB locations without usable coordinates:  %d\n",
            sum(is.na(pbdb_locations$Latitude) | is.na(pbdb_locations$Longitude))))
if (n_multi > 0) {
  cat(sprintf(paste0(
    "  NOTE: %d collection names appear with more than one set of coordinates.\n",
    "        A collection is removed if ANY of its locations matches FAUNMAP.\n"),
    n_distinct(pbdb_locations$Collection[duplicated(pbdb_locations$Collection)])))
}

# =============================================================================
# 5. CHECK EACH PBDB LOCATION AGAINST FAUNMAP
# =============================================================================

cat("\n=== 5. CHECKING FOR DUPLICATES ===\n")

# Distance (m) from each point to its nearest FAUNMAP locality, using the
# haversine formula (same Earth radius as geosphere::distHaversine).
# Points are processed in chunks so memory use stays small.
nearest_point <- function(lat, lon, ref_lat, ref_lon, chunk = 1000) {
  rad <- pi / 180
  index <- rep(NA_integer_, length(lat))
  dist <- rep(NA_real_, length(lat))
  todo <- which(!is.na(lat) & !is.na(lon))

  for (rows in split(todo, ceiling(seq_along(todo) / chunk))) {
    dlat <- outer(lat[rows], ref_lat, "-") * rad
    dlon <- outer(lon[rows], ref_lon, "-") * rad
    a <- sin(dlat / 2)^2 +
      outer(cos(lat[rows] * rad), cos(ref_lat * rad)) * sin(dlon / 2)^2
    d <- 2 * 6378137 * asin(pmin(sqrt(a), 1))

    j <- max.col(-d, ties.method = "first")
    index[rows] <- j
    dist[rows] <- d[cbind(seq_along(rows), j)]
  }

  list(index = index, dist = dist)
}

nearest <- nearest_point(pbdb_locations$Latitude, pbdb_locations$Longitude,
                         faunmap_coords$Latitude, faunmap_coords$Longitude)

reasons <- c(
  "Name match",
  "Exact coordinates (rounded to 6 decimals)",
  paste("Within", duplicate_km, "km")
)

pbdb_locations <- pbdb_locations %>%
  mutate(
    Name_Matched_Site = faunmap_sites$FAUNMAP_SiteName[
      match(name_key(Collection), name_key(faunmap_sites$FAUNMAP_SiteName))
    ],
    Name_Match = !is.na(Name_Matched_Site),
    Exact_Coordinate_Match = !is.na(Latitude) & !is.na(Longitude) &
      paste(round(Latitude, 6), round(Longitude, 6)) %in% faunmap_coord_keys,
    Nearest_FAUNMAP_Site = faunmap_coords$FAUNMAP_SiteName[nearest$index],
    Nearest_FAUNMAP_Period = faunmap_coords$FAUNMAP_Period[nearest$index],
    Distance_km = round(nearest$dist / 1000, 3),
    Removal_Reason = case_when(
      Name_Match ~ reasons[1],
      Exact_Coordinate_Match ~ reasons[2],
      !is.na(Distance_km) & Distance_km <= duplicate_km ~ reasons[3],
      TRUE ~ NA_character_
    )
  )

# One row per removed collection, keeping its strongest reason.
pbdb_removed <- pbdb_locations %>%
  filter(!is.na(Removal_Reason)) %>%
  arrange(Collection, match(Removal_Reason, reasons), Distance_km) %>%
  distinct(Collection, .keep_all = TRUE) %>%
  left_join(count(pbdb, Collection, name = "Occurrences_Removed"), by = "Collection") %>%
  select(Collection, Latitude, Longitude, Removal_Reason, Name_Matched_Site, Nearest_FAUNMAP_Site,
         Nearest_FAUNMAP_Period, Distance_km, Occurrences_Removed)

cat("  Collections removed, by reason:\n")
for (r in reasons) {
  cat(sprintf("      %-45s %6d\n", r, sum(pbdb_removed$Removal_Reason == r)))
}

cat("\n  First removed collections (full list saved in pbdb_removed.csv):\n")
print(head(as.data.frame(pbdb_removed[, c("Collection", "Removal_Reason", "Name_Matched_Site",
                                          "Nearest_FAUNMAP_Site", "Distance_km")]), 10),
      row.names = FALSE)

# =============================================================================
# 6. REMOVE ALL OCCURRENCES OF DUPLICATE COLLECTIONS
# =============================================================================

cat("\n=== 6. REMOVING DUPLICATES ===\n")

pbdb_clean <- pbdb %>%
  filter(!Collection %in% pbdb_removed$Collection)

stopifnot(nrow(pbdb_clean) + sum(pbdb_removed$Occurrences_Removed) == nrow(pbdb))

cat(sprintf("  PBDB occurrences:  %7d original | %7d removed | %7d retained\n",
            nrow(pbdb), nrow(pbdb) - nrow(pbdb_clean), nrow(pbdb_clean)))
cat(sprintf("  PBDB collections:  %7d original | %7d removed | %7d retained\n",
            n_distinct(pbdb$Collection), nrow(pbdb_removed),
            n_distinct(pbdb_clean$Collection)))

# The saved file keeps all original PBDB columns plus standardised
# Latitude/Longitude, which later steps use.
pbdb_clean <- pbdb_clean %>%
  select(-PBDB_Row, -Collection)

# =============================================================================
# 7. SAVE
# =============================================================================

cat("\n=== 7. SAVING ===\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_csv <- function(df, name) {
  path <- file.path(output_dir, paste0(name, ".csv"))
  write.csv(df, path, row.names = FALSE, na = "")
  cat(sprintf("  %-16s %7d rows -> %s\n", name, nrow(df), path))
}

save_csv(pbdb_clean, "pbdb_clean")          # PBDB without FAUNMAP duplicates
save_csv(pbdb_removed, "pbdb_removed")      # one row per removed collection
save_csv(pbdb_locations, "pbdb_locations")  # every location and its check

cat("\n=== STEP 2 COMPLETE ===\n")
cat("Objects in your Environment: pbdb_clean, pbdb_removed, pbdb_locations,\n",
    "faunmap_sites\n", sep = "")
