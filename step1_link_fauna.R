# =============================================================================
# STEP 1: LINK FAUNMAP FAUNA TO THEIR LOCALITIES
#
# Each fauna record is matched to its locality using the pair
#   Machine Number + Analysis Unit
# and receives the locality's SiteName, latitude and longitude.
# If the locality files have no match, the SiteName is taken from the
# site-reference table instead (no coordinates are available then).
#
# Every fauna record is kept, matched or not.
# No five-fauna-per-site filter is applied.
#
# Run the whole file (Ctrl+Shift+S or Ctrl+Shift+Enter in RStudio).
# Progress is printed after every stage; the objects stay in the
# Environment pane so you can inspect them.
# =============================================================================

library(dplyr)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

blancan_fauna_file      <- "Blancan Fauna - Filtered.csv"
irvingtonian_fauna_file <- "Irvingtonian Fauna - Filtered.csv"
site_reference_file     <- "SiteName_Reference_Table.csv"
blancan_loc_file        <- "Blancan Localities Data (Updated).csv"
irvingtonian_loc_file   <- "Irvingtonian Localities Data (Updated).csv"

output_dir <- file.path(work_dir, "Outputs", "1_linked")

# =============================================================================
# 1. CHECK THAT THE FOLDER AND FILES EXIST
# =============================================================================

cat("=== 1. CHECKING INPUT FILES ===\n")

if (!dir.exists(work_dir)) {
  stop("The folder does not exist:\n  ", work_dir)
}

needed <- c(blancan_fauna_file, irvingtonian_fauna_file, site_reference_file,
            blancan_loc_file, irvingtonian_loc_file)
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

read_input <- function(file) {
  df <- read.csv(
    file.path(work_dir, file),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character",          # keep identifiers exactly as written
    na.strings = c("", "NA", "N/A")
  )

  # CSVs saved from Excel often carry a hidden byte-order mark on the first
  # column name, and column names may have stray spaces. Remove both.
  # Working on raw bytes avoids errors from the invisible characters.
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))

  cat(sprintf("  %-45s %7d rows, %3d columns\n", file, nrow(df), ncol(df)))
  df
}

blancan_fauna_raw      <- read_input(blancan_fauna_file)
irvingtonian_fauna_raw <- read_input(irvingtonian_fauna_file)
site_reference_raw     <- read_input(site_reference_file)
blancan_loc_raw        <- read_input(blancan_loc_file)
irvingtonian_loc_raw   <- read_input(irvingtonian_loc_file)

# =============================================================================
# 3. FIND THE COLUMNS WE NEED IN EACH FILE
#
# Column names are compared ignoring case, spaces, dots and underscores,
# so "Machine Number", "Machine.Number" and "machine_number" all match.
# =============================================================================

cat("\n=== 3. IDENTIFYING COLUMNS ===\n")

simplify_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

column_options <- list(
  machine  = c("machinenumber", "machineno", "machine"),
  analysis = c("analysisunit"),
  site     = c("sitename"),
  lat      = c("latdd", "latitude", "lat"),
  lon      = c("longdd", "longitude", "long", "lng", "lon")
)

find_col <- function(df, field, file, required = TRUE) {
  hit <- names(df)[simplify_name(names(df)) %in% column_options[[field]]]

  if (length(hit) == 1) return(hit)

  if (length(hit) > 1) {
    stop("More than one '", field, "' column in ", file, ": ",
         paste(hit, collapse = ", "))
  }
  if (required) {
    stop("No '", field, "' column found in ", file,
         ".\nColumns in this file: ", paste(names(df), collapse = " | "))
  }
  NA_character_
}

show_cols <- function(df, file, fields) {
  cols <- sapply(fields, function(f) find_col(df, f, file))
  cat("  ", file, "\n", sep = "")
  for (f in fields) cat(sprintf("      %-9s -> \"%s\"\n", f, cols[[f]]))
  cols
}

blancan_fauna_cols      <- show_cols(blancan_fauna_raw, blancan_fauna_file, c("machine", "analysis"))
irvingtonian_fauna_cols <- show_cols(irvingtonian_fauna_raw, irvingtonian_fauna_file, c("machine", "analysis"))
site_reference_cols     <- show_cols(site_reference_raw, site_reference_file, c("machine", "analysis", "site"))
blancan_loc_cols        <- show_cols(blancan_loc_raw, blancan_loc_file, c("machine", "analysis", "site", "lat", "lon"))
irvingtonian_loc_cols   <- show_cols(irvingtonian_loc_raw, irvingtonian_loc_file, c("machine", "analysis", "site", "lat", "lon"))

# =============================================================================
# 4. BUILD MATCHING KEYS AND LOOKUP TABLES
# =============================================================================

cat("\n=== 4. BUILDING LOOKUP TABLES ===\n")

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

# 1234.00 and 1234 become the same key; the original column is untouched.
machine_key <- function(x) sub("^([0-9]+)\\.0+$", "\\1", clean_text(x))

to_coord <- function(x, limit) {
  v <- suppressWarnings(as.numeric(clean_text(x)))
  v[!is.na(v) & abs(v) > limit] <- NA     # impossible values -> missing
  v
}

# Pairs that carry more than one SiteName are reported, not fatal:
# the first name is used and the conflicts are saved for checking.
name_conflicts <- list()

make_lookup <- function(df, cols, label, with_coords) {
  keyed <- data.frame(
    .mk = machine_key(df[[cols[["machine"]]]]),
    .ak = clean_text(df[[cols[["analysis"]]]]),
    SiteName = clean_text(df[[cols[["site"]]]]),
    stringsAsFactors = FALSE
  )
  if (with_coords) {
    keyed$Latitude  <- to_coord(df[[cols[["lat"]]]], 90)
    keyed$Longitude <- to_coord(df[[cols[["lon"]]]], 180)
  }

  dropped <- sum(is.na(keyed$.mk) | is.na(keyed$.ak))
  keyed <- filter(keyed, !is.na(.mk), !is.na(.ak))

  conflicts <- keyed %>%
    group_by(.mk, .ak) %>%
    filter(n_distinct(SiteName, na.rm = TRUE) > 1) %>%
    ungroup() %>%
    distinct(.mk, .ak, SiteName) %>%
    mutate(Source = label)
  if (nrow(conflicts) > 0) name_conflicts[[label]] <<- conflicts

  lookup <- keyed %>%
    group_by(.mk, .ak) %>%
    summarise(
      SiteName = SiteName[!is.na(SiteName)][1],
      Rows_For_Pair = n(),
      .groups = "drop"
    )

  if (with_coords) {
    # Coordinates come from the first row of the pair that has both values.
    coords <- keyed %>%
      filter(!is.na(Latitude), !is.na(Longitude)) %>%
      distinct(.mk, .ak, .keep_all = TRUE) %>%
      select(.mk, .ak, Latitude, Longitude)
    lookup <- left_join(lookup, coords, by = c(".mk", ".ak"))
  }

  cat(sprintf(
    "  %-28s %6d unique pairs | %4d rows without a full pair | %4d pairs with conflicting names\n",
    label, nrow(lookup), dropped, n_distinct(paste(conflicts$.mk, conflicts$.ak)) * (nrow(conflicts) > 0)
  ))
  lookup
}

reference_lookup <- make_lookup(site_reference_raw, site_reference_cols,
                                "Site reference table", with_coords = FALSE) %>%
  rename(Reference_SiteName = SiteName, Reference_Rows = Rows_For_Pair)

blancan_loc_lookup <- make_lookup(blancan_loc_raw, blancan_loc_cols,
                                  "Blancan localities", with_coords = TRUE) %>%
  rename(Locality_SiteName = SiteName, Locality_Rows = Rows_For_Pair)

irvingtonian_loc_lookup <- make_lookup(irvingtonian_loc_raw, irvingtonian_loc_cols,
                                       "Irvingtonian localities", with_coords = TRUE) %>%
  rename(Locality_SiteName = SiteName, Locality_Rows = Rows_For_Pair)

# =============================================================================
# 5. LINK EACH FAUNA FILE
# =============================================================================

cat("\n=== 5. LINKING FAUNA TO LOCALITIES ===\n")

link_fauna <- function(fauna, cols, loc_lookup, period) {
  keyed <- fauna %>%
    mutate(
      .mk = machine_key(.data[[cols[["machine"]]]]),
      .ak = clean_text(.data[[cols[["analysis"]]]])
    )

  linked <- keyed %>%
    left_join(reference_lookup, by = c(".mk", ".ak")) %>%
    left_join(loc_lookup, by = c(".mk", ".ak")) %>%
    mutate(
      FAUNMAP_Period = period,
      SiteName_Linked = coalesce(Locality_SiteName, Reference_SiteName),
      Match_Source = case_when(
        !is.na(Locality_SiteName) ~ "Locality file",
        !is.na(Reference_SiteName) ~ "Reference table only",
        TRUE ~ "No match"
      ),
      Has_Coordinates = !is.na(Latitude) & !is.na(Longitude)
    ) %>%
    select(-.mk, -.ak) %>%
    relocate(FAUNMAP_Period, SiteName_Linked, Latitude, Longitude,
             Match_Source, Has_Coordinates, .after = last_col())

  if (nrow(linked) != nrow(fauna)) {
    stop(period, ": linking changed the number of fauna records (",
         nrow(fauna), " -> ", nrow(linked), ").")
  }

  cat("  ", period, ":\n", sep = "")
  cat(sprintf("      fauna records ............... %6d\n", nrow(linked)))
  cat(sprintf("      matched to a locality ....... %6d\n", sum(linked$Match_Source == "Locality file")))
  cat(sprintf("      site name from reference only %6d\n", sum(linked$Match_Source == "Reference table only")))
  cat(sprintf("      not matched at all .......... %6d\n", sum(linked$Match_Source == "No match")))
  cat(sprintf("      with latitude/longitude ..... %6d\n", sum(linked$Has_Coordinates)))
  linked
}

blancan_fauna_linked <- link_fauna(blancan_fauna_raw, blancan_fauna_cols,
                                   blancan_loc_lookup, "Blancan")
irvingtonian_fauna_linked <- link_fauna(irvingtonian_fauna_raw, irvingtonian_fauna_cols,
                                        irvingtonian_loc_lookup, "Irvingtonian")

unmatched_fauna <- bind_rows(blancan_fauna_linked, irvingtonian_fauna_linked) %>%
  filter(Match_Source != "Locality file")

name_conflicts <- bind_rows(name_conflicts) %>%
  rename(any_of(c(Machine_Number = ".mk", Analysis_Unit = ".ak")))

# =============================================================================
# 6. SAVE
# =============================================================================

cat("\n=== 6. SAVING ===\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_csv <- function(df, name) {
  path <- file.path(output_dir, paste0(name, ".csv"))
  write.csv(df, path, row.names = FALSE, na = "")
  cat(sprintf("  %-22s %7d rows -> %s\n", name, nrow(df), path))
}

save_csv(blancan_fauna_linked, "blancan_fauna")
save_csv(irvingtonian_fauna_linked, "irvingtonian_fauna")
save_csv(unmatched_fauna, "unmatched_fauna")
if (nrow(name_conflicts) > 0) save_csv(name_conflicts, "site_name_conflicts")

cat("\n=== STEP 1 COMPLETE ===\n")
cat("Objects in your Environment: blancan_fauna_linked, irvingtonian_fauna_linked,\n",
    "unmatched_fauna, reference_lookup, blancan_loc_lookup, irvingtonian_loc_lookup\n", sep = "")
