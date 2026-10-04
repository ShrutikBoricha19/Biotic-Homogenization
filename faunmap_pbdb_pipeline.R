# =============================================================================
# FAUNMAP + PBDB: LINK, DE-DUPLICATE AND BIN BY GEOLOGICAL STAGE
#
# Step 1  Link FAUNMAP fauna to their localities (Machine Number + Analysis
#         Unit), adding site name, latitude/longitude and age information.
# Step 2  Remove PBDB collections that duplicate a FAUNMAP locality
#         (name match > identical coordinates > within `duplicate_km`).
# Step 3  Assign FAUNMAP localities, FAUNMAP fauna and PBDB occurrences to a
#         geological stage using the MIDPOINT of each age range, and build a
#         combined site index with latitude/longitude for every site.
#
# No sampling-bias or minimum-fauna-per-site filter is applied anywhere.
# Records whose age range crosses a stage boundary are kept (and flagged).
#
# Outputs (in <work_dir>/Outputs):
#   1_linked/   blancan_fauna, irvingtonian_fauna, unmatched_fauna
#   2_pbdb/     pbdb_clean, pbdb_removed, pbdb_locations
#   3_stages/   faunmap_localities, faunmap_fauna, pbdb_occurrences,
#               site_index, stage_summary, multi_stage_sites, excluded
# =============================================================================

library(dplyr)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

input_files <- c(
  pbdb               = "PBDB_Genus - Species - CAN_MEX_USA.csv",
  blancan_fauna      = "Blancan Fauna - Filtered.csv",
  irvingtonian_fauna = "Irvingtonian Fauna - Filtered.csv",
  site_reference     = "SiteName_Reference_Table.csv",
  blancan_loc        = "Blancan Localities Data (Updated).csv",
  irvingtonian_loc   = "Irvingtonian Localities Data (Updated).csv"
)

output_dir   <- file.path(work_dir, "Outputs")
duplicate_km <- 5     # PBDB collection this close to a FAUNMAP locality = duplicate
write_excel  <- TRUE  # also save .xlsx copies (skipped if writexl is missing)

# Accepted spellings of each FAUNMAP column across the input files.
col_names <- list(
  machine  = c("Machine Number", "Machine.Number", "MachineNumber"),
  analysis = c("Analysis Unit", "Analysis unit", "AnalysisUnit", "Analysis.Unit"),
  site     = c("SiteName", "Site Name", "Site.Name"),
  lat      = c("LATDD", "Latitude"),
  lon      = c("LONGDD", "Longitude"),
  max_age  = c("MaximumAge", "Maximum Age", "Maximum.Age"),
  min_age  = c("MinimumAge", "Minimum Age", "Minimum.Age")
)

# PBDB columns used (standard PBDB download names).
pbdb_cols <- c("collection_name", "lat", "lng", "max_ma", "min_ma")

stage_names <- c("Zanclean", "Piacenzian", "Gelasian", "Calabrian", "Chibanian")
stage_labels <- c(
  "Bin 1 - Zanclean (4.700-3.600 Ma)",
  "Bin 2 - Piacenzian (3.600-2.580 Ma)",
  "Bin 3 - Gelasian (2.580-1.800 Ma)",
  "Bin 4 - Calabrian (1.800-0.7741 Ma)",
  "Bin 5 - Chibanian (0.7741-0.129 Ma)"
)

# =============================================================================
# HELPERS
# =============================================================================

read_input <- function(file) {
  path <- file.path(work_dir, file)
  if (!file.exists(path)) stop("Input file not found: ", path)

  # Everything is read as text so identifiers are not altered on import.
  read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    na.strings = c("", "NA", "N/A")
  )
}

find_col <- function(df, field, file) {
  matches <- intersect(col_names[[field]], names(df))
  if (length(matches) != 1) {
    stop(
      "Could not identify exactly one '", field, "' column in ", file,
      ".\nAvailable columns: ", paste(names(df), collapse = ", ")
    )
  }
  matches
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

# 1234.00 and 1234 become the same key; the original column is untouched.
machine_key <- function(x) sub("^([0-9]+)\\.0+$", "\\1", clean_text(x))

name_key <- function(x) gsub("\\s+", " ", tolower(clean_text(x)))

to_num <- function(x) suppressWarnings(as.numeric(clean_text(x)))

# Out-of-range coordinates are treated as missing.
to_lat <- function(x) { v <- to_num(x); v[!is.na(v) & abs(v) > 90]  <- NA; v }
to_lon <- function(x) { v <- to_num(x); v[!is.na(v) & abs(v) > 180] <- NA; v }

first_non_na <- function(x) x[!is.na(x)][1]

# A boundary age belongs to the younger stage (3.600 Ma = Piacenzian).
# Both outer limits (4.700 and 0.129 Ma) are included.
assign_stage <- function(age) {
  case_when(
    age > 3.600  & age <= 4.700  ~ 1L,
    age > 2.580  & age <= 3.600  ~ 2L,
    age > 1.800  & age <= 2.580  ~ 3L,
    age > 0.7741 & age <= 1.800  ~ 4L,
    age >= 0.129 & age <= 0.7741 ~ 5L,
    TRUE ~ NA_integer_
  )
}

# Adds midpoint, stage, NALMA label (reference only), a boundary-crossing
# flag (diagnostic only) and, for records with no stage, the reason why.
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
      NALMA = case_when(
        is.na(Midpoint_Ma) ~ NA_character_,
        Midpoint_Ma > 1.800 & Midpoint_Ma <= 4.700 ~ "Blancan",
        Midpoint_Ma >= 0.129 & Midpoint_Ma <= 1.800 ~ "Irvingtonian",
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
        Midpoint_Ma > 4.700 ~ "Midpoint older than 4.700 Ma",
        Midpoint_Ma < 0.129 ~ "Midpoint younger than 0.129 Ma",
        TRUE ~ "Could not assign midpoint"
      )
    )
}

# Great-circle distance (m) from each point to its nearest reference point.
# Same Earth radius as geosphere::distHaversine. Processed in chunks so the
# distance matrix stays small.
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

save_out <- function(df, subdir, name) {
  dir <- file.path(output_dir, subdir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(df, file.path(dir, paste0(name, ".csv")), row.names = FALSE, na = "")
  if (write_excel && requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(df, file.path(dir, paste0(name, ".xlsx")))
  }
  cat("  Saved", file.path(subdir, name), "(", nrow(df), "rows )\n")
}

# =============================================================================
# LOAD INPUTS AND PREPARE FAUNMAP LOCALITIES
# =============================================================================

cat("=== LOADING INPUT FILES ===\n")
raw <- lapply(input_files, read_input)
for (nm in names(raw)) cat(sprintf("  %-20s %7d rows\n", nm, nrow(raw[[nm]])))

prep_localities <- function(df, period, file) {
  col <- function(field) find_col(df, field, file)

  df %>%
    mutate(
      FAUNMAP_Period = period,
      Row = row_number(),
      Machine_Key = machine_key(.data[[col("machine")]]),
      Analysis_Key = clean_text(.data[[col("analysis")]]),
      SiteName_Std = clean_text(.data[[col("site")]]),
      Latitude = to_lat(.data[[col("lat")]]),
      Longitude = to_lon(.data[[col("lon")]])
    ) %>%
    add_stage(col("max_age"), col("min_age"))
}

localities <- bind_rows(
  prep_localities(raw$blancan_loc, "Blancan", input_files[["blancan_loc"]]),
  prep_localities(
    raw$irvingtonian_loc, "Irvingtonian", input_files[["irvingtonian_loc"]]
  )
)

# =============================================================================
# STEP 1: LINK FAUNMAP FAUNA TO LOCALITIES
# =============================================================================

cat("\n=== STEP 1: LINKING FAUNMAP FAUNA TO LOCALITIES ===\n")

stop_on_name_conflicts <- function(df, label) {
  conflicts <- df %>%
    summarise(n_names = n_distinct(SiteName_Std, na.rm = TRUE), .groups = "drop") %>%
    filter(n_names > 1)

  if (nrow(conflicts) > 0) {
    print(head(conflicts, 20))
    stop("The same Machine Number + Analysis Unit pair has several site names in ", label)
  }
}

# One row per Machine Number + Analysis Unit pair in the site-reference table.
reference_lookup <- local({
  df <- raw$site_reference
  file <- input_files[["site_reference"]]

  keyed <- df %>%
    transmute(
      .mk = machine_key(.data[[find_col(df, "machine", file)]]),
      .ak = clean_text(.data[[find_col(df, "analysis", file)]]),
      SiteName_Std = clean_text(.data[[find_col(df, "site", file)]])
    ) %>%
    filter(!is.na(.mk), !is.na(.ak)) %>%
    group_by(.mk, .ak)

  stop_on_name_conflicts(keyed, file)
  summarise(keyed, Reference_SiteName = first_non_na(SiteName_Std), .groups = "drop")
})

# One row per period + pair in the locality files: site name, coordinates
# and the stage information of the first locality row for that pair.
locality_lookup <- local({
  keyed <- localities %>%
    filter(!is.na(Machine_Key), !is.na(Analysis_Key)) %>%
    group_by(FAUNMAP_Period, .mk = Machine_Key, .ak = Analysis_Key)

  stop_on_name_conflicts(keyed, "the locality files")

  keyed %>%
    summarise(
      Locality_SiteName = first_non_na(SiteName_Std),
      coord_row = which(!is.na(Latitude) & !is.na(Longitude))[1],
      Latitude = Latitude[coord_row],
      Longitude = Longitude[coord_row],
      Locality_Rows = n(),
      across(
        c(Max_Ma, Min_Ma, Midpoint_Ma, Stage_Number, Stage, Stage_Label,
          NALMA, Crosses_Stage_Boundary, Exclusion_Reason),
        first
      ),
      .groups = "drop"
    ) %>%
    select(-coord_row)
})

repeated_pairs <- sum(locality_lookup$Locality_Rows > 1)
if (repeated_pairs > 0) {
  warning(
    repeated_pairs, " Machine Number + Analysis Unit pairs occur on more than ",
    "one locality row; fauna take the age/stage of the first row. ",
    "See Locality_Rows in the linked fauna files."
  )
}

link_fauna <- function(fauna, period, file) {
  keyed <- fauna %>%
    mutate(
      .mk = machine_key(.data[[find_col(fauna, "machine", file)]]),
      .ak = clean_text(.data[[find_col(fauna, "analysis", file)]])
    )

  lookup <- locality_lookup %>%
    filter(FAUNMAP_Period == period) %>%
    select(-FAUNMAP_Period)

  # Columns added below replace any same-named columns in the fauna file.
  added <- c(names(reference_lookup), names(lookup),
             "FAUNMAP_Period", "SiteName_Linked", "Locality_Matched")
  keyed <- keyed[setdiff(names(keyed), setdiff(added, c(".mk", ".ak")))]

  linked <- keyed %>%
    left_join(reference_lookup, by = c(".mk", ".ak")) %>%
    left_join(lookup, by = c(".mk", ".ak")) %>%
    mutate(
      FAUNMAP_Period = period,
      SiteName_Linked = coalesce(Locality_SiteName, Reference_SiteName),
      Locality_Matched = !is.na(Locality_SiteName),
      Exclusion_Reason = if_else(
        Locality_Matched, Exclusion_Reason, "No matching FAUNMAP locality"
      )
    )

  if (nrow(linked) != nrow(fauna)) {
    stop(period, ": linking changed the number of fauna records.")
  }

  cat(sprintf(
    "  %-13s %6d fauna | %6d matched to a locality | %6d with a site name | %6d with coordinates\n",
    period, nrow(linked), sum(linked$Locality_Matched),
    sum(!is.na(linked$SiteName_Linked)),
    sum(!is.na(linked$Latitude) & !is.na(linked$Longitude))
  ))

  linked %>%
    select(-.mk, -.ak) %>%
    relocate(FAUNMAP_Period, SiteName_Linked, Locality_Matched, Latitude, Longitude,
             .after = all_of(names(fauna)[ncol(fauna)]))
}

blancan_fauna <- link_fauna(
  raw$blancan_fauna, "Blancan", input_files[["blancan_fauna"]]
)
irvingtonian_fauna <- link_fauna(
  raw$irvingtonian_fauna, "Irvingtonian", input_files[["irvingtonian_fauna"]]
)
faunmap_fauna_all <- bind_rows(blancan_fauna, irvingtonian_fauna)

save_out(blancan_fauna, "1_linked", "blancan_fauna")
save_out(irvingtonian_fauna, "1_linked", "irvingtonian_fauna")
save_out(filter(faunmap_fauna_all, !Locality_Matched), "1_linked", "unmatched_fauna")

# =============================================================================
# STEP 2: REMOVE PBDB COLLECTIONS THAT DUPLICATE FAUNMAP LOCALITIES
# =============================================================================

cat("\n=== STEP 2: REMOVING PBDB DUPLICATES OF FAUNMAP LOCALITIES ===\n")

missing_pbdb <- setdiff(pbdb_cols, names(raw$pbdb))
if (length(missing_pbdb) > 0) {
  stop("The PBDB file is missing: ", paste(missing_pbdb, collapse = ", "))
}

pbdb <- raw$pbdb %>%
  mutate(
    collection_name = clean_text(collection_name),
    Latitude = to_lat(lat),
    Longitude = to_lon(lng)
  )

if (anyNA(pbdb$collection_name)) {
  stop("Some PBDB records have no collection_name. Resolve these first.")
}

# All FAUNMAP localities (any age) are used as the duplicate reference.
faunmap_names <- unique(na.omit(name_key(localities$SiteName_Std)))
faunmap_coords <- localities %>%
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  distinct(SiteName_Std, FAUNMAP_Period, Latitude, Longitude)
faunmap_coord_keys <- paste(
  round(faunmap_coords$Latitude, 6), round(faunmap_coords$Longitude, 6)
)

# One row per distinct collection_name + coordinate combination.
pbdb_locations <- pbdb %>%
  distinct(collection_name, Latitude, Longitude)

n_multi <- sum(duplicated(pbdb_locations$collection_name))
if (n_multi > 0) {
  warning(
    n_multi, " extra coordinate rows found for reused PBDB collection names. ",
    "A collection is removed if ANY of its locations matches a FAUNMAP locality."
  )
}

nearest <- nearest_point(
  pbdb_locations$Latitude, pbdb_locations$Longitude,
  faunmap_coords$Latitude, faunmap_coords$Longitude
)

reasons <- c(
  "Name match",
  "Exact coordinates (rounded to 6 decimals)",
  paste("Within", duplicate_km, "km")
)

pbdb_locations <- pbdb_locations %>%
  mutate(
    Name_Match = name_key(collection_name) %in% faunmap_names,
    Exact_Coordinate_Match = !is.na(Latitude) & !is.na(Longitude) &
      paste(round(Latitude, 6), round(Longitude, 6)) %in% faunmap_coord_keys,
    Nearest_FAUNMAP_Site = faunmap_coords$SiteName_Std[nearest$index],
    Nearest_FAUNMAP_Period = faunmap_coords$FAUNMAP_Period[nearest$index],
    Distance_km = nearest$dist / 1000,
    Removal_Reason = case_when(
      Name_Match ~ reasons[1],
      Exact_Coordinate_Match ~ reasons[2],
      !is.na(Distance_km) & Distance_km <= duplicate_km ~ reasons[3],
      TRUE ~ NA_character_
    )
  )

# Strongest reason per collection: name > exact coordinates > distance.
pbdb_removed <- pbdb_locations %>%
  filter(!is.na(Removal_Reason)) %>%
  arrange(collection_name, match(Removal_Reason, reasons), Distance_km) %>%
  distinct(collection_name, .keep_all = TRUE) %>%
  select(collection_name, Latitude, Longitude, Removal_Reason,
         Nearest_FAUNMAP_Site, Nearest_FAUNMAP_Period, Distance_km)

# Remove every occurrence belonging to a removed collection.
pbdb_clean <- filter(pbdb, !collection_name %in% pbdb_removed$collection_name)

cat("  PBDB occurrences:", nrow(pbdb), "| retained:", nrow(pbdb_clean),
    "| removed:", nrow(pbdb) - nrow(pbdb_clean), "\n")
cat("  PBDB collections removed:", nrow(pbdb_removed), "\n")
print(table(factor(pbdb_removed$Removal_Reason, levels = reasons)))

save_out(select(pbdb_clean, -Latitude, -Longitude), "2_pbdb", "pbdb_clean")
save_out(pbdb_removed, "2_pbdb", "pbdb_removed")
save_out(pbdb_locations, "2_pbdb", "pbdb_locations")

# =============================================================================
# STEP 3: GEOLOGICAL-STAGE BINNING (MIDPOINT RULE)
# =============================================================================

cat("\n=== STEP 3: ASSIGNING GEOLOGICAL STAGES ===\n")

pbdb_all <- pbdb_clean %>%
  mutate(Row = row_number()) %>%
  add_stage("max_ma", "min_ma")

# Retained records; Exclusion_Reason is blank for all of them.
faunmap_localities <- filter(localities, !is.na(Stage_Number)) %>% select(-Exclusion_Reason)
faunmap_fauna <- filter(faunmap_fauna_all, !is.na(Stage_Number)) %>% select(-Exclusion_Reason)
pbdb_occurrences <- filter(pbdb_all, !is.na(Stage_Number)) %>% select(-Exclusion_Reason)

# Combined site index: one row per site and stage, with coordinates.
#   FAUNMAP site = period + Machine Number + Analysis Unit
#   PBDB site    = collection_name
index_cols <- c(
  "Database", "Dataset", "SiteName", "Site_Key", "Machine_Key", "Analysis_Key",
  "Latitude", "Longitude", "Stage_Number", "Stage", "Stage_Label", "NALMA",
  "Midpoint_Ma", "Crosses_Stage_Boundary"
)

site_index <- bind_rows(
  faunmap_localities %>%
    mutate(
      Database = "FAUNMAP",
      Dataset = paste("FAUNMAP", FAUNMAP_Period),
      SiteName = SiteName_Std,
      Site_Key = if_else(
        !is.na(Machine_Key) & !is.na(Analysis_Key),
        paste(FAUNMAP_Period, Machine_Key, Analysis_Key, sep = " | "),
        paste0(FAUNMAP_Period, " | UNIDENTIFIED ROW ", Row)
      )
    ),
  pbdb_occurrences %>%
    mutate(
      Database = "PBDB",
      Dataset = "PBDB",
      SiteName = collection_name,
      Site_Key = paste("PBDB", collection_name, sep = " | "),
      Machine_Key = NA_character_,
      Analysis_Key = NA_character_
    )
) %>%
  select(all_of(index_cols)) %>%
  distinct() %>%
  arrange(Stage_Number, Database, SiteName, Site_Key)

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
  left_join(
    bind_rows(
      faunmap_fauna %>% mutate(Dataset = paste("FAUNMAP", FAUNMAP_Period)),
      pbdb_occurrences %>% mutate(Dataset = "PBDB")
    ) %>%
      count(Dataset, Stage_Number, name = "n_fauna_records"),
    by = c("Dataset", "Stage_Number")
  ) %>%
  mutate(n_fauna_records = coalesce(n_fauna_records, 0L)) %>%
  relocate(n_fauna_records, .after = n_sites) %>%
  arrange(Stage_Number, Dataset)

multi_stage_sites <- site_index %>%
  group_by(Database, Site_Key, SiteName) %>%
  summarise(
    n_stages = n_distinct(Stage_Number),
    stages = paste(stage_names[sort(unique(Stage_Number))], collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_stages > 1)

excluded <- bind_rows(
  localities %>%
    filter(is.na(Stage_Number)) %>%
    transmute(Dataset = paste("FAUNMAP", FAUNMAP_Period, "localities"), Row,
              SiteName = SiteName_Std, Max_Ma, Min_Ma, Midpoint_Ma, Exclusion_Reason),
  faunmap_fauna_all %>%
    mutate(Row = ave(seq_along(FAUNMAP_Period), FAUNMAP_Period, FUN = seq_along)) %>%
    filter(is.na(Stage_Number)) %>%
    transmute(Dataset = paste("FAUNMAP", FAUNMAP_Period, "fauna"), Row,
              SiteName = SiteName_Linked, Max_Ma, Min_Ma, Midpoint_Ma, Exclusion_Reason),
  pbdb_all %>%
    filter(is.na(Stage_Number)) %>%
    transmute(Dataset = "PBDB occurrences", Row,
              SiteName = collection_name, Max_Ma, Min_Ma, Midpoint_Ma, Exclusion_Reason)
)

# Every record must be either retained or listed as excluded.
stopifnot(
  nrow(faunmap_localities) + sum(grepl("localities$", excluded$Dataset)) == nrow(localities),
  nrow(faunmap_fauna) + sum(grepl("fauna$", excluded$Dataset)) == nrow(faunmap_fauna_all),
  nrow(pbdb_occurrences) + sum(excluded$Dataset == "PBDB occurrences") == nrow(pbdb_all)
)

cat("  FAUNMAP localities with a stage:", nrow(faunmap_localities), "of", nrow(localities), "\n")
cat("  FAUNMAP fauna with a stage:     ", nrow(faunmap_fauna), "of", nrow(faunmap_fauna_all), "\n")
cat("  PBDB occurrences with a stage:  ", nrow(pbdb_occurrences), "of", nrow(pbdb_all), "\n")
cat("  Sites in more than one stage:   ", nrow(multi_stage_sites), "\n\n")
print(as.data.frame(stage_summary))
cat("\nExclusion reasons:\n")
print(as.data.frame(count(excluded, Dataset, Exclusion_Reason)))

save_out(faunmap_localities, "3_stages", "faunmap_localities")
save_out(faunmap_fauna, "3_stages", "faunmap_fauna")
save_out(pbdb_occurrences, "3_stages", "pbdb_occurrences")
save_out(site_index, "3_stages", "site_index")
save_out(stage_summary, "3_stages", "stage_summary")
save_out(multi_stage_sites, "3_stages", "multi_stage_sites")
save_out(excluded, "3_stages", "excluded")

cat("\n=== COMPLETE ===\nOutputs saved in:", output_dir, "\n")
cat("No five-fauna-per-site filter was applied.\n")
