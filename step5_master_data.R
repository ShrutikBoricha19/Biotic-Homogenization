# =============================================================================
# STEP 5: MASTER DATA FILE (FAUNMAP + PBDB)
#
# One table of every occurrence identified to species, from both databases,
# after the sp./cf. resolution of Step 3d. Columns:
#
#   Species          corrected "Genus species" name (Step 3d)
#   Genus
#   Site             site name (FAUNMAP SiteName / PBDB collection name)
#   Site_Key         unique site identifier (as in Steps 3-4)
#   Latitude, Longitude
#   Max_Ma, Min_Ma, Midpoint_Ma   age range of the site and its midpoint
#   Spatial_Bin      physiographic province (USA) / subregion (Canada)
#   Physio_Division  physiographic division (USA) / region (Canada)
#   State_Province, Country
#   Time_Bin         time bin number;  Time_Bin_Label  e.g. "Bin 1 (4-3.25 Ma)"
#   Original_Name    name before Step 3d (differs only for resolved sp./cf.)
#   Name_Status      "as identified" or "resolved from next bin (...)"
#   Source           FAUNMAP or PBDB (for tracing a record back only)
#
# Occurrences identified above genus level (no genus) are not included.
# Mexican sites have no physiographic province, so their Spatial_Bin is empty.
#
# Outputs (Outputs/5_master_data/):
#   master_data.csv / .xlsx           every occurrence (one row each)
#   master_data_unique.csv / .xlsx    one row per species x site x time bin
#   master_data_summary.csv           species, sites and occurrences per
#                                     spatial bin and time bin
#
# Input: Outputs/3d_resolved/ (Step 3d)
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 3d.
# Needs: dplyr (writexl optional, for the .xlsx copies).
# =============================================================================

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is not installed. Run: install.packages(\"dplyr\")")
}
library(dplyr)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

input_dir  <- file.path("Outputs", "3d_resolved")
output_dir <- file.path(work_dir, "Outputs", "5_master_data")
save_excel <- TRUE

# =============================================================================
# HELPERS
# =============================================================================

read_text_csv <- function(name) {
  path <- file.path(work_dir, input_dir, name)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun Step 3d first.")
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
                 na.strings = c("", "NA"), encoding = "UTF-8")
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))
  cat(sprintf("  %-26s %7d rows\n", name, nrow(df)))
  df
}
col <- function(df, name) if (name %in% names(df)) df[[name]] else rep(NA_character_, nrow(df))
clean <- function(x) { x <- trimws(gsub("\\s+", " ", as.character(x))); x[x == ""] <- NA; x }
num <- function(x) suppressWarnings(as.numeric(x))
is_open <- function(x) is.na(x) | grepl("^(sp|spp|indet|indeterminate)\\.?$", tolower(x))

use_excel <- save_excel && requireNamespace("writexl", quietly = TRUE)
save_out <- function(df, name, excel = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, paste0(name, ".csv"))
  ok <- tryCatch({ write.csv(df, path, row.names = FALSE, na = "", fileEncoding = "UTF-8"); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) stop("Could not write ", path, "\n  Close it (e.g. in Excel) and run again.")
  if (excel && use_excel) writexl::write_xlsx(df, file.path(output_dir, paste0(name, ".xlsx")))
  cat(sprintf("  %-26s %7d rows -> %s\n", name, nrow(df), path))
}

# =============================================================================
# 1. READ THE STEP 3d OCCURRENCES
# =============================================================================

cat("=== 1. READING STEP 3d OUTPUTS ===\n")
fauna     <- read_text_csv("faunmap_fauna.csv")
pbdb      <- read_text_csv("pbdb_occurrences.csv")
time_bins <- read_text_csv("time_bins.csv")

# =============================================================================
# 2. ONE TABLE FOR BOTH DATABASES
# =============================================================================

cat("\n=== 2. BUILDING THE MASTER TABLE ===\n")

status <- function(res) ifelse(is.na(res) | res == "unchanged", "as identified",
                               sub("^resolved \\((.*)\\)$", "resolved from next bin (\\1)", res))

fa <- fauna %>% transmute(
  Genus = clean(Genus), Epithet = tolower(clean(Species)),
  Site = clean(coalesce(col(fauna, "SiteName_Linked"), col(fauna, "SiteName"))),
  Site_Key, Latitude = num(Latitude), Longitude = num(Longitude),
  Max_Ma = num(Max_Ma), Min_Ma = num(Min_Ma), Midpoint_Ma = num(Midpoint_Ma),
  Spatial_Bin = clean(col(fauna, "Physio_Province")), Physio_Division = clean(col(fauna, "Physio_Division")),
  State_Province = clean(col(fauna, "State_Province")), Country = clean(col(fauna, "Country")),
  Time_Bin = as.integer(Stage_Number),
  Original_Name = clean(trimws(paste(coalesce(col(fauna, "IDConfidence_original"), ""),
                                     coalesce(col(fauna, "Genus_original"), Genus),
                                     coalesce(col(fauna, "Species_original"), "")))),
  Name_Status = status(col(fauna, "Resolution")),
  Source = "FAUNMAP")

name_col <- intersect(c("accepted_name", "identified_name"), names(pbdb))[1]
pb_name <- clean(pbdb[[name_col]])
pb <- pbdb %>% transmute(
  Genus = clean(sub(" .*$", "", pb_name)),
  Epithet = tolower(ifelse(grepl(" ", coalesce(pb_name, "")), sub("^\\S+\\s+", "", pb_name), NA)),
  Site = clean(coalesce(col(pbdb, "Collection"), col(pbdb, "collection_name"))),
  Site_Key, Latitude = num(Latitude), Longitude = num(Longitude),
  Max_Ma = num(Max_Ma), Min_Ma = num(Min_Ma), Midpoint_Ma = num(Midpoint_Ma),
  Spatial_Bin = clean(col(pbdb, "Physio_Province")), Physio_Division = clean(col(pbdb, "Physio_Division")),
  State_Province = clean(col(pbdb, "State_Province")), Country = clean(col(pbdb, "Country")),
  Time_Bin = as.integer(Stage_Number),
  Original_Name = clean(coalesce(col(pbdb, paste0(name_col, "_original")), pb_name)),
  Name_Status = status(col(pbdb, "Resolution")),
  Source = "PBDB")

all_occ <- bind_rows(fa, pb)
n_all <- nrow(all_occ)

master_data <- all_occ %>%
  filter(!is.na(Genus), !is_open(Epithet)) %>%
  mutate(Genus = paste0(toupper(substr(Genus, 1, 1)), tolower(substring(Genus, 2))),
         Epithet = sub("\\s.*$", "", Epithet),          # first word of the epithet
         Species = paste(Genus, Epithet),
         Time_Bin_Label = time_bins$Bin_Label[match(Time_Bin, as.integer(time_bins$Bin_Number))]) %>%
  select(Species, Genus, Site, Site_Key, Latitude, Longitude, Max_Ma, Min_Ma, Midpoint_Ma,
         Spatial_Bin, Physio_Division, State_Province, Country, Time_Bin, Time_Bin_Label,
         Original_Name, Name_Status, Source) %>%
  arrange(Time_Bin, Spatial_Bin, Site, Species)

cat(sprintf("  %d occurrences in Step 3d | %d identified to species (kept) | %d above genus level (left out)\n",
            n_all, nrow(master_data), sum(is.na(all_occ$Genus))))
if (any(!is.na(all_occ$Genus) & is_open(all_occ$Epithet))) {
  cat(sprintf("  NOTE: %d genus-level records (sp.) were still present and are left out.\n",
              sum(!is.na(all_occ$Genus) & is_open(all_occ$Epithet))))
}
cat(sprintf("  %d species | %d sites | names resolved in Step 3d: %d occurrences\n",
            n_distinct(master_data$Species), n_distinct(master_data$Site_Key),
            sum(master_data$Name_Status != "as identified")))

# One row per species x site x time bin (repeated records of a species at a
# site collapse into one).
master_data_unique <- master_data %>%
  group_by(Species, Site_Key, Time_Bin) %>%
  summarise(n_occurrences = n(),
            Name_Status = if (all(Name_Status == "as identified")) "as identified" else "includes resolved names",
            across(c(Genus, Site, Latitude, Longitude, Max_Ma, Min_Ma, Midpoint_Ma, Spatial_Bin,
                     Physio_Division, State_Province, Country, Time_Bin_Label, Source), first),
            .groups = "drop") %>%
  select(Species, Genus, Site, Site_Key, Latitude, Longitude, Max_Ma, Min_Ma, Midpoint_Ma,
         Spatial_Bin, Physio_Division, State_Province, Country, Time_Bin, Time_Bin_Label,
         n_occurrences, Name_Status, Source) %>%
  arrange(Time_Bin, Spatial_Bin, Site, Species)

master_data_summary <- master_data %>%
  group_by(Time_Bin, Time_Bin_Label, Spatial_Bin) %>%
  summarise(n_species = n_distinct(Species), n_sites = n_distinct(Site_Key),
            n_occurrences = n(), .groups = "drop") %>%
  arrange(Time_Bin, desc(n_sites))

# Checks: every row has a species name, a site, a time bin and coordinates.
stopifnot(!any(is.na(master_data$Species)), !any(is.na(master_data$Site_Key)),
          !any(is.na(master_data$Time_Bin)))
cat(sprintf("  Rows without coordinates: %d | without a spatial bin: %d (e.g. Mexico)\n",
            sum(is.na(master_data$Latitude) | is.na(master_data$Longitude)),
            sum(is.na(master_data$Spatial_Bin))))

# =============================================================================
# 3. SAVE
# =============================================================================

cat("\n=== 3. SAVING ===\n")
save_out(master_data, "master_data")
save_out(master_data_unique, "master_data_unique")
save_out(master_data_summary, "master_data_summary")

cat("\n  Species and sites per time bin:\n")
print(as.data.frame(master_data %>% group_by(Time_Bin) %>%
                      summarise(species = n_distinct(Species), sites = n_distinct(Site_Key),
                                occurrences = n(), .groups = "drop")), row.names = FALSE)

cat("\n=== STEP 5 COMPLETE ===\n")
cat("Objects in your Environment: master_data, master_data_unique, master_data_summary\n")
