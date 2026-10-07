# =============================================================================
# STEP 1: FAUNMAP LOCALITIES AND THEIR FAUNA
#
# Localities: the three FAUNMAP site workbooks
#   Blancan_Sites.xlsx, Irvingtonian_Sites.xlsx, Rancholabrean_Sites.xlsx
# Fauna: faunalf.csv, linked to the localities by
#   Machine Number + Analysis Unit
# Each fauna record receives its locality's period, SiteName, latitude and
# longitude.
#
# Localities with NO age at all are KEPT: they receive the minimum and maximum
# age of their land-mammal age (NALMA; 'nalma_ages' below), and the column
# Age_Source says so. Their midpoint is therefore the NALMA midpoint.
#
# Localities with an improper age are REMOVED (with all their fauna):
#   - only one of minimum / maximum age given (or not a number)
#   - maximum age younger than minimum age, or a negative age
#   - age inconsistent with the locality's land-mammal age (NALMA): the
#     midpoint lies more than 'nalma_tolerance' Myr outside the NALMA window
#     (e.g. an Irvingtonian site with a maximum age of 22 Ma)
# Every removed locality is listed with its reason in removed_localities.csv.
#
# Outputs (Outputs/1_linked/), used by Steps 2, 2b and 3:
#   blancan_localities.csv, irvingtonian_localities.csv, rancholabrean_localities.csv
#   blancan_fauna.csv, irvingtonian_fauna.csv, rancholabrean_fauna.csv
#   removed_localities.csv   localities removed for their age, and why
#   unmatched_fauna.csv      faunalf.csv records not linked to a kept locality
#   link_summary.csv         counts per period
#
# Run the whole file (Ctrl+Shift+S in RStudio).
# Needs: dplyr, readxl.
# =============================================================================

for (pkg in c("dplyr", "readxl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

locality_files <- c(
  Blancan       = "Blancan_Sites.xlsx",
  Irvingtonian  = "Irvingtonian_Sites.xlsx",
  Rancholabrean = "Rancholabrean_Sites.xlsx"
)
fauna_file <- "faunalf.csv"

output_dir <- file.path(work_dir, "Outputs", "1_linked")

# Ages given to localities with no age at all (Ma). Bell et al. (2004):
# Blancan 4.9-1.35, Irvingtonian 1.35-0.21, Rancholabrean 0.21-0.0117.
fill_no_age_from_nalma <- TRUE
nalma_ages <- list(
  Blancan       = c(older = 4.90, younger = 1.35),
  Irvingtonian  = c(older = 1.35, younger = 0.21),
  Rancholabrean = c(older = 0.21, younger = 0.0117)
)

# Age checks. NALMA windows (Ma) used to catch impossible ages; a locality is
# removed when its midpoint lies more than 'nalma_tolerance' outside its window.
check_nalma_ages <- TRUE
nalma_windows <- list(
  Blancan       = c(older = 4.90, younger = 1.35),
  Irvingtonian  = c(older = 1.90, younger = 0.21),
  Rancholabrean = c(older = 0.30, younger = 0.0117)
)
nalma_tolerance <- 0.5    # Myr

# =============================================================================
# HELPERS
# =============================================================================

simplify_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))
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
to_num <- function(x) suppressWarnings(as.numeric(clean_text(x)))
to_coord <- function(x, limit) { v <- to_num(x); v[!is.na(v) & abs(v) > limit] <- NA; v }

# First column matching one of the options (case, spaces, dots, underscores ignored).
find_col <- function(df, options, label, file, required = TRUE) {
  simple <- simplify_name(names(df))
  for (o in options) {
    hit <- names(df)[simple == o]
    if (length(hit) >= 1) return(hit[1])
  }
  if (required) stop("No '", label, "' column found in ", file,
                     ".\nColumns: ", paste(names(df), collapse = " | "))
  NA_character_
}

# Finds a file, allowing spaces or underscores in its name.
locate <- function(file) {
  path <- file.path(work_dir, file)
  if (file.exists(path)) return(path)
  pattern <- paste0("^", gsub("[ _]", "[ _]", gsub("\\.", "\\\\.", file)), "$")
  hit <- list.files(work_dir, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
  if (length(hit)) return(hit[1])
  stop("File not found:\n  ", path)
}

save_csv <- function(df, name) {
  path <- file.path(output_dir, paste0(name, ".csv"))
  ok <- tryCatch({ write.csv(df, path, row.names = FALSE, na = ""); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) {
    path <- sub("\\.csv$", "_new.csv", path)
    write.csv(df, path, row.names = FALSE, na = "")
    cat("  NOTE: file was open elsewhere - saved as", basename(path), "\n")
  }
  cat(sprintf("  %-28s %7d rows -> %s\n", name, nrow(df), path))
}

# =============================================================================
# 1. READ THE LOCALITY WORKBOOKS
# =============================================================================

cat("=== 1. READING FAUNMAP LOCALITIES ===\n")
if (!dir.exists(work_dir)) stop("The folder does not exist:\n  ", work_dir)

loc_raw <- lapply(names(locality_files), function(p) {
  path <- locate(locality_files[[p]])
  df <- as.data.frame(readxl::read_excel(path, col_types = "text"), stringsAsFactors = FALSE)
  names(df) <- trimws(names(df))
  cat(sprintf("  %-28s %6d rows, %3d columns\n", basename(path), nrow(df), ncol(df)))
  df
})
names(loc_raw) <- names(locality_files)

# =============================================================================
# 2. CHECK AGES; REMOVE LOCALITIES WITHOUT A PROPER AGE
# =============================================================================

cat("\n=== 2. CHECKING AGES ===\n")

loc_checked <- lapply(names(loc_raw), function(p) {
  df <- loc_raw[[p]]; f <- locality_files[[p]]
  c_mach <- find_col(df, c("machinenumber", "machineno", "machine"), "Machine Number", f)
  c_anal <- find_col(df, c("analysisunit"), "Analysis Unit", f)
  c_site <- find_col(df, c("sitename"), "SiteName", f)
  c_min  <- find_col(df, c("minimumage", "minage"), "MinimumAge", f)
  c_max  <- find_col(df, c("maximumage", "maxage"), "MaximumAge", f)
  c_lat  <- find_col(df, c("latdd", "latitude", "lat"), "latitude", f)
  c_lon  <- find_col(df, c("longdd", "longitude", "long", "lng", "lon"), "longitude", f)
  win <- nalma_windows[[p]]
  fill <- nalma_ages[[p]]
  out <- df %>% mutate(
    FAUNMAP_Period = p,
    .mk = machine_key(.data[[c_mach]]), .ak = analysis_key(.data[[c_anal]]),
    .site = clean_text(.data[[c_site]]),
    .min = to_num(.data[[c_min]]), .max = to_num(.data[[c_max]]),
    .no_age = is.na(.min) & is.na(.max) & fill_no_age_from_nalma,
    Age_Source = ifelse(.no_age, sprintf("NALMA %s (%s-%s Ma): no age in file", p,
                                         fill[["older"]], fill[["younger"]]), "file"),
    .min = ifelse(.no_age, fill[["younger"]], .min),
    .max = ifelse(.no_age, fill[["older"]], .max),
    .mid = (.min + .max) / 2,
    .lat = to_coord(.data[[c_lat]], 90), .lon = to_coord(.data[[c_lon]], 180),
    Age_Problem = case_when(
      is.na(.min) & is.na(.max)                   ~ "no age",
      is.na(.min) | is.na(.max)                   ~ "minimum or maximum age missing",
      .min < 0 | .max < 0                         ~ "negative age",
      .max < .min                                 ~ "maximum age younger than minimum age",
      check_nalma_ages & .mid > win[["older"]] + nalma_tolerance ~
        sprintf("age too old for %s (midpoint %.3g Ma)", p, .mid),
      check_nalma_ages & .mid < win[["younger"]] - nalma_tolerance ~
        sprintf("age too young for %s (midpoint %.3g Ma)", p, .mid),
      TRUE ~ NA_character_),
    Missing_Pair = is.na(.mk) | is.na(.ak))
  # Write the NALMA ages into the original age columns (used by Steps 2 and 3).
  out[[c_min]] <- ifelse(out$.no_age, as.character(out$.min), out[[c_min]])
  out[[c_max]] <- ifelse(out$.no_age, as.character(out$.max), out[[c_max]])
  out
})
names(loc_checked) <- names(loc_raw)

removed_localities <- bind_rows(lapply(loc_checked, function(df) {
  df %>% filter(!is.na(Age_Problem)) %>%
    transmute(FAUNMAP_Period, Machine_Number = .mk, Analysis_Unit = .ak, SiteName = .site,
              MinimumAge = .min, MaximumAge = .max, Reason = sub(" \\(midpoint.*$", "", Age_Problem),
              Detail = Age_Problem)
}))

for (p in names(loc_checked)) {
  df <- loc_checked[[p]]
  cat(sprintf("  %-14s %5d localities | %5d given NALMA ages (no age) | %5d removed for their age | %5d kept\n",
              p, nrow(df), sum(df$.no_age), sum(!is.na(df$Age_Problem)), sum(is.na(df$Age_Problem))))
}
cat("\n  Removed localities by reason:\n")
print(as.data.frame(count(removed_localities, FAUNMAP_Period, Reason)), row.names = FALSE)

loc_kept <- lapply(loc_checked, function(df) filter(df, is.na(Age_Problem)))

# Lookup: one row per Machine Number + Analysis Unit (pairs are unique across
# the three workbooks; duplicates within a workbook use the first row).
loc_lookup <- bind_rows(lapply(loc_kept, function(df) {
  df %>% filter(!Missing_Pair) %>%
    transmute(.mk, .ak, FAUNMAP_Period, SiteName_Linked = .site, Latitude = .lat, Longitude = .lon)
}))
dup_pairs <- loc_lookup %>% count(.mk, .ak) %>% filter(n > 1)
if (nrow(dup_pairs) > 0) {
  cat(sprintf("\n  NOTE: %d Machine Number + Analysis Unit pairs occur on more than one kept locality row;\n",
              nrow(dup_pairs)), "        the first row is used for their fauna.\n")
}
loc_lookup <- distinct(loc_lookup, .mk, .ak, .keep_all = TRUE)

all_pairs <- bind_rows(lapply(loc_checked, function(df) {
  df %>% filter(!Missing_Pair) %>% transmute(.mk, .ak, .removed_period = FAUNMAP_Period,
                                             .age_problem = Age_Problem)
})) %>% distinct(.mk, .ak, .keep_all = TRUE)

# =============================================================================
# 3. LINK FAUNA (faunalf.csv) TO THE LOCALITIES
# =============================================================================

cat("\n=== 3. LINKING FAUNA ===\n")

fauna_path <- locate(fauna_file)
fauna_raw <- read.csv(fauna_path, check.names = FALSE, stringsAsFactors = FALSE,
                      colClasses = "character", na.strings = c("", "NA", "N/A"))
names(fauna_raw) <- trimws(sub("^[^A-Za-z0-9]+", "", names(fauna_raw), useBytes = TRUE))
cat(sprintf("  %-28s %6d rows, %3d columns\n", basename(fauna_path), nrow(fauna_raw), ncol(fauna_raw)))

f_mach <- find_col(fauna_raw, c("machinenumber", "machineno", "machine"), "Machine Number", fauna_file)
f_anal <- find_col(fauna_raw, c("analysisunit"), "Analysis Unit", fauna_file)
cat(sprintf("  Linking on '%s' + '%s'\n", f_mach, f_anal))

# Columns this step adds; same-named columns already in faunalf.csv are kept as *_original.
added <- c("FAUNMAP_Period", "SiteName_Linked", "Latitude", "Longitude", "Match_Source", "Has_Coordinates")
fauna <- fauna_raw
clash <- intersect(names(fauna), added)
if (length(clash)) names(fauna)[names(fauna) %in% clash] <- paste0(clash, "_original")

fauna_linked <- fauna %>%
  mutate(.mk = machine_key(.data[[f_mach]]), .ak = analysis_key(.data[[f_anal]])) %>%
  left_join(loc_lookup, by = c(".mk", ".ak")) %>%
  left_join(all_pairs, by = c(".mk", ".ak")) %>%
  mutate(Match_Source = case_when(
           !is.na(FAUNMAP_Period) ~ "Locality file",
           !is.na(.age_problem)   ~ paste0("Locality removed (", .removed_period, ": ",
                                           sub(" \\(midpoint.*$", "", .age_problem), ")"),
           is.na(.mk) | is.na(.ak) ~ "No Machine Number / Analysis Unit",
           TRUE                   ~ "Not in the three locality workbooks"),
         Has_Coordinates = !is.na(Latitude) & !is.na(Longitude))
stopifnot(nrow(fauna_linked) == nrow(fauna_raw))

linked_by_period <- lapply(names(locality_files), function(p) {
  fauna_linked %>% filter(FAUNMAP_Period %in% p) %>%
    select(-.mk, -.ak, -.removed_period, -.age_problem)
})
names(linked_by_period) <- names(locality_files)
unmatched_fauna <- fauna_linked %>% filter(is.na(FAUNMAP_Period)) %>%
  select(-.mk, -.ak, -.removed_period, -.age_problem)

for (p in names(linked_by_period)) {
  cat(sprintf("  %-14s %7d fauna records linked (%d with coordinates)\n", p,
              nrow(linked_by_period[[p]]), sum(linked_by_period[[p]]$Has_Coordinates)))
}
cat(sprintf("  %-14s %7d fauna records not linked:\n", "Unmatched", nrow(unmatched_fauna)))
print(as.data.frame(count(unmatched_fauna, Match_Source, sort = TRUE)), row.names = FALSE)

# Kept localities with no fauna at all (Step 2b removes these).
no_fauna <- sapply(names(loc_kept), function(p) {
  k <- loc_kept[[p]]
  sum(!paste(k$.mk, k$.ak) %in% paste(linked_by_period[[p]][[f_mach]] %>% machine_key(),
                                      analysis_key(linked_by_period[[p]][[f_anal]])))
})

# =============================================================================
# 4. SAVE
# =============================================================================

cat("\n=== 4. SAVING ===\n")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (p in names(locality_files)) {
  stub <- tolower(p)
  save_csv(loc_kept[[p]] %>% select(-starts_with("."), -Age_Problem, -Missing_Pair),
           paste0(stub, "_localities"))
  save_csv(linked_by_period[[p]], paste0(stub, "_fauna"))
}
save_csv(removed_localities, "removed_localities")
save_csv(unmatched_fauna, "unmatched_fauna")

link_summary <- data.frame(
  Period = names(locality_files),
  localities_in_file = sapply(loc_checked, nrow),
  removed_for_age = sapply(loc_checked, function(d) sum(!is.na(d$Age_Problem))),
  localities_kept = sapply(loc_kept, nrow),
  kept_with_NALMA_ages = sapply(loc_kept, function(d) sum(d$.no_age)),
  kept_localities_without_fauna = no_fauna,
  fauna_records_linked = sapply(linked_by_period, nrow),
  row.names = NULL)
save_csv(link_summary, "link_summary")
print(link_summary, row.names = FALSE)

cat("\n=== STEP 1 COMPLETE ===\n")
cat("Objects in your Environment: loc_kept, linked_by_period, removed_localities, unmatched_fauna\n")
