# =============================================================================
# STEP 3: ASSIGN FAUNMAP AND PBDB RECORDS TO TIME BINS
#
# Each record is placed in ONE numbered time bin using the MIDPOINT of its
# age range. Bins are 'bin_width' (0.75 Myr) long, from 'bin_start' (4.00 Ma)
# to 'bin_end' (0.0117 Ma = end of the Late Pleistocene, 11,700 years ago);
# Bin 1 starts at the round age of 4 Ma and is therefore 0.80 Myr long:
#   Bin 1  4.00 - 3.20 Ma      Bin 4  1.70 - 0.95 Ma
#   Bin 2  3.20 - 2.45 Ma      Bin 5  0.95 - 0.20 Ma
#   Bin 3  2.45 - 1.70 Ma      Bin 6  0.20 - 0.0117 Ma (shorter: 0.19 Myr)
# A boundary age belongs to the younger bin (e.g. 3.20 Ma = Bin 2); both
# outer limits are included. The last bin is shorter than the others because
# the interval ends at 0.0117 Ma; set 'merge_short_last_bin' to TRUE to join
# it to the bin before it. The bins are saved to time_bins.csv, which every
# later step reads, so a change here carries through the whole pipeline.
#
# Before binning, sites in Alaska and Canadian sites north of 60 degrees N
# (north of British Columbia, Alberta, Saskatchewan and Manitoba: Yukon,
# Northwest Territories, Nunavut and the far north of Quebec) are removed.
#
# Records whose age range crosses a bin boundary are KEPT (and flagged).
# No five-fauna-per-site filter is applied.
# Latitude/longitude are carried into every output.
#
# Inputs (all from Step 2b, which cleans the Step 1 and Step 2 outputs):
#   Outputs/2b_clean/blancan_localities.csv, irvingtonian_localities.csv
#   Outputs/2b_clean/blancan_fauna.csv, irvingtonian_fauna.csv
#   Outputs/2b_clean/pbdb_clean.csv
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Steps 1, 2 and 2b.
# =============================================================================

for (pkg in c("dplyr", "sf", "maps")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

# Cleaned files from Step 2b (Canada/USA/Mexico only; no marine mammals or bats).
blancan_loc_file      <- file.path("Outputs", "2b_clean", "blancan_localities.csv")
irvingtonian_loc_file <- file.path("Outputs", "2b_clean", "irvingtonian_localities.csv")

blancan_fauna_file      <- file.path("Outputs", "2b_clean", "blancan_fauna.csv")
irvingtonian_fauna_file <- file.path("Outputs", "2b_clean", "irvingtonian_fauna.csv")
pbdb_file               <- file.path("Outputs", "2b_clean", "pbdb_clean.csv")

output_dir <- file.path(work_dir, "Outputs", "3_stages")

save_excel <- TRUE    # also save .xlsx copies (skipped if writexl is not installed)

# Time bins.
bin_width <- 0.75     # Myr
bin_grid  <- 3.95     # Ma, where the 0.75-Myr steps start (3.95, 3.20, 2.45 ...)
bin_start <- 4.00     # Ma, older limit of Bin 1 (extended by 0.05 Myr to the round age of 4 Ma)
bin_end   <- 0.0117   # Ma, end of the Late Pleistocene (11,700 years ago)
merge_short_last_bin <- FALSE   # TRUE: a last bin shorter than half the width joins the bin before it

# Geographic limits.
remove_alaska      <- TRUE
canada_north_limit <- 60      # degrees N; Canadian sites north of this are removed (NA = keep all)
coast_buffer_km    <- 25      # sites in the sea this close to Alaska count as Alaska

bin_bounds <- round(seq(bin_grid, bin_end, by = -bin_width), 6)
bin_bounds[1] <- bin_start            # Bin 1 runs from 4.00 Ma, so it is 0.80 Myr long
if (min(bin_bounds) > bin_end) bin_bounds <- c(bin_bounds, bin_end)
if (merge_short_last_bin && length(bin_bounds) > 2) {
  n_b <- length(bin_bounds)
  if (bin_bounds[n_b - 1] - bin_bounds[n_b] < bin_width / 2) bin_bounds <- bin_bounds[-(n_b - 1)]
}
n_bins      <- length(bin_bounds) - 1
stage_older <- bin_bounds[-length(bin_bounds)]
stage_young <- bin_bounds[-1]
fmt_age     <- function(x) trimws(formatC(x, format = "fg", digits = 4))
stage_names <- paste("Bin", seq_len(n_bins))
stage_labels <- sprintf("Bin %d (%s-%s Ma)", seq_len(n_bins), fmt_age(stage_older), fmt_age(stage_young))

# =============================================================================
# 1. CHECK THAT THE FOLDER AND FILES EXIST
# =============================================================================

cat("=== 1. CHECKING INPUT FILES ===\n")

if (!dir.exists(work_dir)) {
  stop("The folder does not exist:\n  ", work_dir)
}

needed <- c(blancan_loc_file, irvingtonian_loc_file,
            blancan_fauna_file, irvingtonian_fauna_file, pbdb_file)
found <- file.exists(file.path(work_dir, needed))

for (i in seq_along(needed)) {
  cat(if (found[i]) "  [found]   " else "  [MISSING] ", needed[i], "\n", sep = "")
}

if (!all(found)) {
  stop("Some input files were not found. Run Steps 1, 2 and 2b first, ",
       "or fix the file names in SETTINGS.")
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

  # Remove hidden Excel byte-order marks and stray spaces from column names.
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))

  cat(sprintf("  %-45s %7d rows, %3d columns\n", file, nrow(df), ncol(df)))
  df
}

blancan_loc_raw        <- read_input(blancan_loc_file)
irvingtonian_loc_raw   <- read_input(irvingtonian_loc_file)
blancan_fauna_raw      <- read_input(blancan_fauna_file)
irvingtonian_fauna_raw <- read_input(irvingtonian_fauna_file)
pbdb_raw               <- read_input(pbdb_file)

# =============================================================================
# 3. FIND THE COLUMNS WE NEED IN EACH FILE
#
# Same rules as Steps 1 and 2: case, spaces, dots and underscores are
# ignored, and the first option listed wins when several are present.
# =============================================================================

cat("\n=== 3. IDENTIFYING COLUMNS ===\n")

simplify_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

column_options <- list(
  machine    = c("machinenumber", "machineno", "machine"),
  analysis   = c("analysisunit"),
  site       = c("sitename", "collectionname"),
  lat        = c("latdd", "latitude", "lat"),
  lon        = c("longdd", "longitude", "lng", "long", "lon"),
  max_age    = c("maximumage", "maxage", "maxma"),
  min_age    = c("minimumage", "minage", "minma")
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

show_cols <- function(df, file, fields) {
  cols <- sapply(fields, function(f) find_col(df, f, file))
  cat("  ", file, "\n", sep = "")
  for (f in fields) cat(sprintf("      %-9s -> \"%s\"\n", f, cols[[f]]))
  cols
}

loc_fields <- c("machine", "analysis", "site", "lat", "lon", "max_age", "min_age")

blancan_loc_cols        <- show_cols(blancan_loc_raw, blancan_loc_file, loc_fields)
irvingtonian_loc_cols   <- show_cols(irvingtonian_loc_raw, irvingtonian_loc_file, loc_fields)
blancan_fauna_cols      <- show_cols(blancan_fauna_raw, blancan_fauna_file, c("machine", "analysis"))
irvingtonian_fauna_cols <- show_cols(irvingtonian_fauna_raw, irvingtonian_fauna_file, c("machine", "analysis"))
pbdb_cols               <- show_cols(pbdb_raw, pbdb_file, c("site", "lat", "lon", "max_age", "min_age"))

# =============================================================================
# 4. HELPER FUNCTIONS FOR AGES AND STAGES
# =============================================================================

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

# 1234.00 and 1234 become the same key (as in Step 1).
machine_key <- function(x) sub("^([0-9]+)\\.0+$", "\\1", clean_text(x))

to_num <- function(x) suppressWarnings(as.numeric(clean_text(x)))

to_coord <- function(x, limit) {
  v <- to_num(x)
  v[!is.na(v) & abs(v) > limit] <- NA     # impossible values -> missing
  v
}

# Bin k holds ages in (younger limit, older limit]; the youngest limit itself
# belongs to the last bin.
assign_stage <- function(age) {
  k <- rep(NA_integer_, length(age))
  for (i in seq_len(n_bins)) k[!is.na(age) & age <= stage_older[i] & age > stage_young[i]] <- i
  k[!is.na(age) & abs(age - bin_end) < 1e-9] <- n_bins
  k
}

# Adds midpoint, stage, NALMA label (for reference only), a flag for age
# ranges that cross a stage boundary (diagnostic only), and the reason a
# record could not be given a stage.
add_stage <- function(df, max_col, min_col) {
  df %>%
    mutate(
      Max_Ma = to_num(.data[[max_col]]),
      Min_Ma = to_num(.data[[min_col]]),
      Midpoint_Ma = if_else(
        !is.na(Max_Ma) & !is.na(Min_Ma) & Max_Ma >= Min_Ma,
        (Max_Ma + Min_Ma) / 2,
        NA_real_
      ),
      Stage_Number = assign_stage(Midpoint_Ma),
      Stage = stage_names[Stage_Number],
      Stage_Label = stage_labels[Stage_Number],
      NALMA = case_when(                       # approximate, for reference only
        is.na(Midpoint_Ma) ~ NA_character_,
        Midpoint_Ma > 1.800 & Midpoint_Ma <= 4.900 ~ "Blancan",
        Midpoint_Ma > 0.210 & Midpoint_Ma <= 1.800 ~ "Irvingtonian",
        Midpoint_Ma >= 0.0117 & Midpoint_Ma <= 0.210 ~ "Rancholabrean",
        TRUE ~ "Outside_Study_Interval"
      ),
      Crosses_Stage_Boundary = if_else(
        is.na(Midpoint_Ma),
        NA,
        coalesce(assign_stage(Max_Ma) != assign_stage(Min_Ma), TRUE)
      ),
      Exclusion_Reason = case_when(
        !is.na(Stage_Number) ~ NA_character_,
        is.na(Max_Ma) | is.na(Min_Ma) ~ "Missing or nonnumeric age",
        Max_Ma < Min_Ma ~ "Maximum age is younger than minimum age",
        Midpoint_Ma > bin_start ~ paste("Midpoint older than", bin_start, "Ma"),
        Midpoint_Ma < bin_end ~ paste("Midpoint younger than", bin_end, "Ma (Holocene)"),
        TRUE ~ "Could not assign midpoint"
      )
    )
}

# Geographic exclusions: Alaska, and Canada north of 'canada_north_limit'.
# Returns a reason, or NA for sites that stay (and for sites without
# coordinates, which cannot be checked).
library(sf)
suppressMessages(sf::sf_use_s2(FALSE))
world <- sf::st_as_sf(maps::map("world", regions = c("USA", "Canada"), fill = TRUE, plot = FALSE))
world <- sf::st_set_crs(sf::st_set_crs(world, NA), 4326)
world <- suppressMessages(suppressWarnings(sf::st_make_valid(world)))
world$Part <- ifelse(grepl("^Canada", world$ID), "Canada", "USA")
world <- suppressMessages(world %>% group_by(Part) %>% summarise(.groups = "drop"))
na_albers <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +datum=WGS84 +units=m"

geo_exclusion <- function(lat, lon) {
  out <- rep(NA_character_, length(lat))
  ok <- !is.na(lat) & !is.na(lon)
  if (!any(ok)) return(out)
  pts <- sf::st_as_sf(data.frame(lon = lon[ok], lat = lat[ok]), coords = c("lon", "lat"), crs = 4326)
  hit <- suppressMessages(sf::st_intersects(pts, world))
  part <- vapply(hit, function(h) if (length(h)) world$Part[h[1]] else NA_character_, character(1))
  miss <- which(is.na(part))
  if (length(miss)) {                         # sites in the sea: nearest land within the buffer
    d <- sf::st_distance(sf::st_transform(pts[miss, ], na_albers), sf::st_transform(world, na_albers))
    near <- apply(d, 1, which.min)
    part[miss] <- ifelse(as.numeric(d[cbind(seq_along(miss), near)]) / 1000 <= coast_buffer_km,
                         world$Part[near], "Sea / elsewhere")
  }
  # Mainland Alaska is part of the "USA" outline; the contiguous USA does not
  # reach 50 degrees N, so any US location north of 50 degrees N is Alaska.
  part[part %in% "USA"] <- ifelse(lat[ok][part %in% "USA"] > 50, "Alaska", "Other USA")
  reason <- rep(NA_character_, length(part))
  if (remove_alaska) reason[part %in% "Alaska"] <- "Alaska"
  if (!is.na(canada_north_limit)) {
    north <- lat[ok] > canada_north_limit & !part %in% c("Alaska", "Other USA")
    reason[is.na(reason) & north] <- paste0("Canada north of ", canada_north_limit, " degrees N")
  }
  out[ok] <- reason
  out
}

# =============================================================================
# 5. FAUNMAP LOCALITIES: ASSIGN STAGES
# =============================================================================

cat("\n=== 5. FAUNMAP LOCALITIES ===\n")

prepare_localities <- function(df, cols, period) {
  df %>%
    mutate(
      FAUNMAP_Period = period,
      Row = row_number(),
      Machine_Key = machine_key(.data[[cols[["machine"]]]]),
      Analysis_Key = clean_text(.data[[cols[["analysis"]]]]),
      SiteName_Std = clean_text(.data[[cols[["site"]]]]),
      Latitude = to_coord(.data[[cols[["lat"]]]], 90),
      Longitude = to_coord(.data[[cols[["lon"]]]], 180),
      Site_Key = if_else(
        !is.na(Machine_Key) & !is.na(Analysis_Key),
        paste(period, Machine_Key, Analysis_Key, sep = " | "),
        paste0(period, " | UNIDENTIFIED ROW ", Row)
      )
    ) %>%
    add_stage(cols[["max_age"]], cols[["min_age"]])
}

localities_all <- bind_rows(
  prepare_localities(blancan_loc_raw, blancan_loc_cols, "Blancan"),
  prepare_localities(irvingtonian_loc_raw, irvingtonian_loc_cols, "Irvingtonian")
)

# Geographic exclusions override the stage: such localities get no stage.
localities_all <- localities_all %>%
  mutate(Geo_Exclusion = geo_exclusion(Latitude, Longitude),
         Exclusion_Reason = coalesce(Geo_Exclusion, Exclusion_Reason),
         across(c(Stage_Number, Stage, Stage_Label), ~ if_else(is.na(Geo_Exclusion), .x, .x[NA_integer_])))
cat(sprintf("  Localities removed: %d in Alaska, %d in Canada north of %s degrees N\n",
            sum(localities_all$Geo_Exclusion %in% "Alaska"),
            sum(grepl("^Canada north", localities_all$Geo_Exclusion)), canada_north_limit))

faunmap_localities <- filter(localities_all, !is.na(Stage_Number))

for (p in c("Blancan", "Irvingtonian")) {
  cat(sprintf("  %-13s %6d locality rows | %6d with a stage | %6d without\n", p,
              sum(localities_all$FAUNMAP_Period == p),
              sum(faunmap_localities$FAUNMAP_Period == p),
              sum(localities_all$FAUNMAP_Period == p & is.na(localities_all$Stage_Number))))
}

# =============================================================================
# 6. FAUNMAP FAUNA: TAKE THE STAGE OF THEIR LOCALITY
#
# Fauna are matched to localities by period + Machine Number + Analysis Unit,
# exactly as in Step 1. If a pair occurs on several locality rows, the
# first row's age is used (the count is reported below).
# =============================================================================

cat("\n=== 6. FAUNMAP FAUNA ===\n")

locality_stage_lookup <- localities_all %>%
  filter(!is.na(Machine_Key), !is.na(Analysis_Key)) %>%
  group_by(FAUNMAP_Period, .mk = Machine_Key, .ak = Analysis_Key) %>%
  summarise(
    Site_Key = first(Site_Key),
    Pair_Locality_Rows = n(),
    Pair_Distinct_Ages = n_distinct(paste(Max_Ma, Min_Ma)),
    across(c(Max_Ma, Min_Ma, Midpoint_Ma, Stage_Number, Stage, Stage_Label,
             NALMA, Crosses_Stage_Boundary, Exclusion_Reason), first),
    .groups = "drop"
  )

n_multi_age <- sum(locality_stage_lookup$Pair_Distinct_Ages > 1)
cat(sprintf("  Machine Number + Analysis Unit pairs with more than one age range: %d\n",
            n_multi_age))
if (n_multi_age > 0) {
  cat("  (fauna of these pairs use the first row's age; see Pair_Distinct_Ages)\n")
}

link_stage <- function(fauna, cols, period) {
  keyed <- fauna %>%
    mutate(
      .mk = machine_key(.data[[cols[["machine"]]]]),
      .ak = clean_text(.data[[cols[["analysis"]]]])
    )

  lookup <- locality_stage_lookup %>%
    filter(FAUNMAP_Period == period) %>%
    select(-FAUNMAP_Period)

  # Keep any same-named column already in the file as *_original.
  clash <- setdiff(intersect(names(fauna), names(lookup)), c(".mk", ".ak"))
  if (length(clash) > 0) {
    names(keyed)[names(keyed) %in% clash] <- paste0(names(keyed)[names(keyed) %in% clash], "_original")
  }

  linked <- keyed %>%
    left_join(lookup, by = c(".mk", ".ak")) %>%
    mutate(
      Fauna_Row = row_number(),
      Exclusion_Reason = if_else(is.na(Site_Key), "No matching FAUNMAP locality",
                                 Exclusion_Reason)
    ) %>%
    select(-.mk, -.ak)

  if (nrow(linked) != nrow(fauna)) {
    stop(period, ": stage linking changed the number of fauna records.")
  }

  cat(sprintf("  %-13s %6d fauna | %6d with a stage | %6d with a stage AND coordinates\n",
              period, nrow(linked), sum(!is.na(linked$Stage_Number)),
              sum(!is.na(linked$Stage_Number) & !is.na(to_num(linked$Latitude)) &
                    !is.na(to_num(linked$Longitude)))))
  linked
}

fauna_all <- bind_rows(
  link_stage(blancan_fauna_raw, blancan_fauna_cols, "Blancan"),
  link_stage(irvingtonian_fauna_raw, irvingtonian_fauna_cols, "Irvingtonian")
) %>%
  mutate(Latitude = to_num(Latitude), Longitude = to_num(Longitude))

faunmap_fauna <- filter(fauna_all, !is.na(Stage_Number))

# =============================================================================
# 7. PBDB OCCURRENCES: ASSIGN STAGES
# =============================================================================

cat("\n=== 7. PBDB OCCURRENCES ===\n")

pbdb_all <- pbdb_raw %>%
  mutate(
    Row = row_number(),
    Collection = clean_text(.data[[pbdb_cols[["site"]]]]),
    Site_Key = if_else(is.na(Collection),
                       paste0("PBDB | UNNAMED ROW ", Row),
                       paste("PBDB", Collection, sep = " | ")),
    Latitude = to_coord(.data[[pbdb_cols[["lat"]]]], 90),
    Longitude = to_coord(.data[[pbdb_cols[["lon"]]]], 180)
  ) %>%
  add_stage(pbdb_cols[["max_age"]], pbdb_cols[["min_age"]])

pbdb_geo <- pbdb_all %>% distinct(Site_Key, Latitude, Longitude) %>%
  mutate(Geo_Exclusion = geo_exclusion(Latitude, Longitude)) %>%
  group_by(Site_Key) %>% summarise(Geo_Exclusion = first(na.omit(Geo_Exclusion))[1], .groups = "drop")
pbdb_all <- pbdb_all %>%
  left_join(pbdb_geo, by = "Site_Key") %>%
  mutate(Exclusion_Reason = coalesce(Geo_Exclusion, Exclusion_Reason),
         across(c(Stage_Number, Stage, Stage_Label), ~ if_else(is.na(Geo_Exclusion), .x, .x[NA_integer_])))
cat(sprintf("  PBDB collections removed: %d in Alaska, %d in Canada north of %s degrees N\n",
            n_distinct(pbdb_all$Site_Key[pbdb_all$Geo_Exclusion %in% "Alaska"]),
            n_distinct(pbdb_all$Site_Key[grepl("^Canada north", pbdb_all$Geo_Exclusion)]),
            canada_north_limit))

pbdb_occurrences <- filter(pbdb_all, !is.na(Stage_Number))

cat(sprintf("  PBDB          %6d occurrences | %6d with a stage | %6d without\n",
            nrow(pbdb_all), nrow(pbdb_occurrences),
            nrow(pbdb_all) - nrow(pbdb_occurrences)))

# =============================================================================
# 8. COMBINED SITE INDEX (ONE ROW PER SITE AND STAGE, WITH COORDINATES)
#
#   FAUNMAP site = period + Machine Number + Analysis Unit
#   PBDB site    = collection_name
# =============================================================================

cat("\n=== 8. BUILDING THE SITE INDEX ===\n")

index_cols <- c("Database", "Dataset", "SiteName", "Site_Key", "Latitude", "Longitude",
                "Stage_Number", "Stage", "Stage_Label", "NALMA", "Midpoint_Ma",
                "Crosses_Stage_Boundary")

site_index <- bind_rows(
  faunmap_localities %>%
    mutate(Database = "FAUNMAP", Dataset = paste("FAUNMAP", FAUNMAP_Period),
           SiteName = SiteName_Std),
  pbdb_occurrences %>%
    mutate(Database = "PBDB", Dataset = "PBDB", SiteName = Collection)
) %>%
  select(all_of(index_cols)) %>%
  distinct() %>%
  arrange(Stage_Number, Database, SiteName, Site_Key)

cat(sprintf("  %d site-stage rows | %d distinct sites | %d rows without coordinates\n",
            nrow(site_index), n_distinct(site_index$Site_Key),
            sum(is.na(site_index$Latitude) | is.na(site_index$Longitude))))

# =============================================================================
# 9. SUMMARIES AND EXCLUSIONS
# =============================================================================

cat("\n=== 9. SUMMARY BY STAGE ===\n")

fauna_counts <- bind_rows(
  faunmap_fauna %>% transmute(Dataset = paste("FAUNMAP", FAUNMAP_Period), Stage_Number),
  pbdb_occurrences %>% transmute(Dataset = "PBDB", Stage_Number)
) %>%
  count(Dataset, Stage_Number, name = "n_fauna_records")

stage_summary <- site_index %>%
  group_by(Dataset, Stage_Number, Stage_Label) %>%
  summarise(
    n_sites = n_distinct(Site_Key),
    n_sites_without_coordinates =
      n_distinct(Site_Key[is.na(Latitude) | is.na(Longitude)]),
    n_sites_crossing_stage_boundary =
      n_distinct(Site_Key[Crosses_Stage_Boundary %in% TRUE]),
    oldest_midpoint_Ma = max(Midpoint_Ma),
    youngest_midpoint_Ma = min(Midpoint_Ma),
    .groups = "drop"
  ) %>%
  left_join(fauna_counts, by = c("Dataset", "Stage_Number")) %>%
  mutate(n_fauna_records = coalesce(n_fauna_records, 0L)) %>%
  relocate(n_fauna_records, .after = n_sites) %>%
  arrange(Stage_Number, Dataset)

print(as.data.frame(stage_summary %>% select(-Stage_Label)), row.names = FALSE)

multi_stage_sites <- site_index %>%
  group_by(Database, Site_Key, SiteName) %>%
  summarise(
    n_stages = n_distinct(Stage_Number),
    stages = paste(stage_names[sort(unique(Stage_Number))], collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_stages > 1)

cat("\n  Sites appearing in more than one stage:", nrow(multi_stage_sites), "\n")

excluded <- bind_rows(
  localities_all %>%
    filter(is.na(Stage_Number)) %>%
    transmute(Dataset = paste("FAUNMAP", FAUNMAP_Period, "localities"), Row,
              SiteName = SiteName_Std, Max_Ma, Min_Ma, Midpoint_Ma, Exclusion_Reason),
  fauna_all %>%
    filter(is.na(Stage_Number)) %>%
    transmute(Dataset = paste("FAUNMAP", FAUNMAP_Period, "fauna"), Row = Fauna_Row,
              SiteName = if ("SiteName_Linked" %in% names(fauna_all)) SiteName_Linked else NA_character_,
              Max_Ma, Min_Ma, Midpoint_Ma, Exclusion_Reason),
  pbdb_all %>%
    filter(is.na(Stage_Number)) %>%
    transmute(Dataset = "PBDB occurrences", Row,
              SiteName = Collection, Max_Ma, Min_Ma, Midpoint_Ma, Exclusion_Reason)
)

cat("\n  Records without a stage, by reason:\n")
print(as.data.frame(count(excluded, Dataset, Exclusion_Reason)), row.names = FALSE)

# Every record must be either retained or listed as excluded.
stopifnot(
  nrow(faunmap_localities) + sum(grepl("localities$", excluded$Dataset)) == nrow(localities_all),
  nrow(faunmap_fauna) + sum(grepl("fauna$", excluded$Dataset)) == nrow(fauna_all),
  nrow(pbdb_occurrences) + sum(excluded$Dataset == "PBDB occurrences") == nrow(pbdb_all)
)
cat("\n  Check passed: every record is either in a stage or listed in excluded.csv\n")

# =============================================================================
# 10. SAVE
# =============================================================================

cat("\n=== 10. SAVING ===\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

use_excel <- save_excel && requireNamespace("writexl", quietly = TRUE)

save_out <- function(df, name) {
  path <- file.path(output_dir, paste0(name, ".csv"))
  write.csv(df, path, row.names = FALSE, na = "")
  if (use_excel) writexl::write_xlsx(df, file.path(output_dir, paste0(name, ".xlsx")))
  cat(sprintf("  %-20s %7d rows -> %s\n", name, nrow(df), path))
}

# Exclusion_Reason is blank for every retained record, so it is dropped here.
save_out(select(faunmap_localities, -Exclusion_Reason), "faunmap_localities")
save_out(select(faunmap_fauna, -Exclusion_Reason), "faunmap_fauna")
save_out(select(pbdb_occurrences, -Exclusion_Reason), "pbdb_occurrences")
save_out(site_index, "site_index")
save_out(stage_summary, "stage_summary")
save_out(multi_stage_sites, "multi_stage_sites")
save_out(excluded, "excluded")

time_bins <- data.frame(Bin_Number = seq_len(n_bins), Bin_Name = stage_names, Bin_Label = stage_labels,
                        Older_Ma = stage_older, Younger_Ma = stage_young,
                        Width_Myr = round(stage_older - stage_young, 4))
save_out(time_bins, "time_bins")
print(time_bins, row.names = FALSE)

cat("\n=== STEP 3 COMPLETE ===\n")
cat("Objects in your Environment: faunmap_localities, faunmap_fauna, pbdb_occurrences,\n",
    "site_index, stage_summary, multi_stage_sites, excluded\n", sep = "")
