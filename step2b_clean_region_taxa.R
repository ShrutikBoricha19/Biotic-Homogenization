# =============================================================================
# STEP 2b: KEEP ONLY CANADA, USA AND MEXICO; REMOVE MARINE MAMMALS AND BATS
#
# Runs after Steps 1 and 2 and before Step 3. Three clean-ups, in this order:
#
#   1. Region: every FAUNMAP locality and PBDB collection outside Canada, the
#      USA and Mexico is removed, with all its occurrences.
#        - Sites with coordinates: tested against the country outlines (maps
#          package). A site in the sea within 'coast_buffer_km' of one of the
#          three countries is kept (coastal sites with rounded coordinates).
#        - Sites without coordinates: the country column is used if the file
#          has one (e.g. PBDB 'cc'). If there is neither, the site is kept or
#          removed according to 'keep_unverified' and listed as unverified.
#   2. Taxa: every occurrence of a marine mammal or a bat (Chiroptera) is
#      removed. A record is removed when its order, family or genus (whichever
#      the file has) is in the lists in SETTINGS:
#        orders   Chiroptera; Cetacea, Sirenia, Desmostylia, Pinnipedia
#        families seals, sea lions, walruses, whales, dolphins, sea cows, bats
#                 (so pinnipeds filed under Carnivora and whales filed under
#                 Artiodactyla are caught too)
#        genera   sea otters (Enhydra, Enhydritherium) and the polar bear
#                 (Ursus maritimus)
#      Missing FAUNMAP orders are filled from faunalf.csv by genus, as in Step 6.
#   3. Empty sites: a FAUNMAP locality or PBDB collection left with no
#      occurrences is removed (with 'drop_localities_without_fauna', this
#      includes FAUNMAP localities that never had fauna records).
#
# Inputs:
#   Outputs/1_linked/<period>_localities.csv and <period>_fauna.csv   (Step 1;
#   period = blancan, irvingtonian, rancholabrean)
#   Outputs/2_pbdb/pbdb_clean.csv                                 (Step 2)
#   faunalf.csv (optional, to fill a missing Order by genus)
#
# Outputs (Outputs/2b_clean/) - same columns as the input files:
#   <period>_localities.csv   used by Step 3
#   <period>_fauna.csv        used by Steps 3 and 6
#   pbdb_clean.csv                                        used by Step 3
#   removed_sites.csv     every removed locality/collection and why
#   removed_records.csv   every removed occurrence and why
#   cleaning_summary.csv  counts before and after each clean-up
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Steps 1 and 2.
# Needs: dplyr, sf, maps.
# =============================================================================

for (pkg in c("dplyr", "sf", "maps")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(sf)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

periods <- c("Blancan", "Irvingtonian", "Rancholabrean")
loc_files   <- setNames(file.path("Outputs", "1_linked", paste0(tolower(periods), "_localities.csv")), periods)
fauna_files <- setNames(file.path("Outputs", "1_linked", paste0(tolower(periods), "_fauna.csv")), periods)
pbdb_file               <- file.path("Outputs", "2_pbdb", "pbdb_clean.csv")
order_lookup_file       <- "faunalf.csv"    # optional

output_dir <- file.path(work_dir, "Outputs", "2b_clean")

coast_buffer_km <- 25      # sea sites this close to Canada/USA/Mexico are kept
keep_unverified <- TRUE    # sites with neither coordinates nor a country: keep?
drop_localities_without_fauna <- TRUE   # FAUNMAP localities with no fauna left are removed

marine_orders <- c("Cetacea", "Sirenia", "Desmostylia", "Pinnipedia")
marine_families <- c(
  # pinnipeds
  "Otariidae", "Phocidae", "Odobenidae", "Desmatophocidae", "Enaliarctidae",
  # whales and dolphins
  "Balaenidae", "Balaenopteridae", "Eschrichtiidae", "Cetotheriidae", "Herpetocetidae",
  "Delphinidae", "Phocoenidae", "Monodontidae", "Physeteridae", "Kogiidae", "Ziphiidae",
  "Pontoporiidae", "Iniidae", "Platanistidae", "Lipotidae", "Kentriodontidae", "Albireonidae",
  # sea cows and desmostylians
  "Dugongidae", "Trichechidae", "Desmostylidae", "Paleoparadoxiidae"
)
bat_orders   <- c("Chiroptera")
bat_families <- c(
  "Vespertilionidae", "Molossidae", "Phyllostomidae", "Mormoopidae", "Natalidae",
  "Noctilionidae", "Emballonuridae", "Thyropteridae", "Furipteridae"
)
marine_genera  <- c("Enhydra", "Enhydritherium")   # sea otters
marine_species <- c("Ursus maritimus")              # polar bear

countries_kept <- c("Canada", "USA", "Mexico")
country_codes  <- c("CA", "CAN", "CANADA", "US", "USA", "UNITED STATES",
                    "UNITED STATES OF AMERICA", "MX", "MEX", "MEXICO", "MÉXICO")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================

read_input <- function(file, required = TRUE) {
  path <- file.path(work_dir, file)
  if (!file.exists(path)) {
    if (required) stop("File not found:\n  ", path, "\nRun Steps 1 and 2 first.")
    return(NULL)
  }
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                 colClasses = "character", na.strings = c("", "NA", "N/A"))
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))
  cat(sprintf("  %-50s %7d rows\n", file, nrow(df)))
  df
}

save_with_fallback <- function(writer, path) {
  ok <- tryCatch({ writer(path); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
  if (ok) return(invisible(path))
  alt <- sub("(\\.[A-Za-z]+)$", "_new\\1", path)
  writer(alt)
  cat("  NOTE: could not overwrite", basename(path), "- saved as", basename(alt), "\n")
  invisible(alt)
}
save_csv <- function(df, name) {
  save_with_fallback(function(f) write.csv(df, f, row.names = FALSE, na = ""),
                     file.path(output_dir, paste0(name, ".csv")))
}

simplify_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# First column whose simplified name matches one of the options (or NA).
find_col <- function(df, options) {
  simple <- simplify_name(names(df))
  for (o in options) {
    hit <- names(df)[simple == o]
    if (length(hit) >= 1) return(hit[1])
  }
  NA_character_
}

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
to_coord <- function(x, limit) {
  v <- suppressWarnings(as.numeric(clean_text(x)))
  v[!is.na(v) & abs(v) > limit] <- NA
  v
}
title_word <- function(x) {
  x <- clean_text(x)
  as.character(ifelse(is.na(x), NA_character_,
                      paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))))
}
# ifelse() that always returns text (also for empty tables).
chr_ifelse <- function(test, yes, no) as.character(ifelse(test, yes, no))
col_or_na <- function(df, col) if (is.na(col)) rep(NA_character_, nrow(df)) else df[[col]]

lat_opts     <- c("latdd", "latitude", "lat")
lon_opts     <- c("longdd", "longitude", "lng", "long", "lon")
country_opts <- c("country", "cc", "countrycode", "nation", "cc2")
machine_opts <- c("machinenumber", "machineno", "machine")
analysis_opts <- c("analysisunit")
site_opts    <- c("sitename", "collectionname")

# =============================================================================
# 1. READ FILES
# =============================================================================

cat("=== 1. READING FILES ===\n")
loc_raw   <- lapply(loc_files, read_input)
fauna_raw <- lapply(fauna_files, read_input)
for (p in names(fauna_raw)) {
  if (nrow(fauna_raw[[p]]) == 0) {
    cat(sprintf(paste0("  NOTE: %s has no fauna records - no faunalf.csv record matched a %s locality\n",
                       "        in Step 1 (check Step 1's unmatched_fauna.csv). Its localities will be\n",
                       "        removed as having no fauna.\n"), basename(fauna_files[[p]]), p))
  }
}
pbdb_raw <- read_input(pbdb_file)
order_lookup_raw <- read_input(order_lookup_file, required = FALSE)

# =============================================================================
# 2. REGION: INSIDE CANADA, USA OR MEXICO?
# =============================================================================

cat("\n=== 2. REGION CHECK ===\n")

suppressMessages(sf::sf_use_s2(FALSE))
countries <- sf::st_as_sf(maps::map("world", regions = countries_kept, fill = TRUE, plot = FALSE))
countries <- sf::st_set_crs(sf::st_set_crs(countries, NA), 4326)   # maps data are plain longitude/latitude
countries <- suppressMessages(suppressWarnings(sf::st_make_valid(countries)))
countries$Country <- sub(":.*$", "", countries$ID)
countries <- suppressMessages(countries %>% group_by(Country) %>% summarise(.groups = "drop"))
na_albers <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +datum=WGS84 +units=m"
countries_m <- sf::st_transform(countries, na_albers)

# Region of a set of points; returns country name, "outside" or NA (no coordinates).
region_of <- function(lat, lon) {
  out <- rep(NA_character_, length(lat))
  ok <- !is.na(lat) & !is.na(lon)
  if (!any(ok)) return(out)
  pts <- sf::st_as_sf(data.frame(lon = lon[ok], lat = lat[ok]), coords = c("lon", "lat"), crs = 4326)
  hit <- suppressMessages(sf::st_intersects(pts, countries))
  inside <- vapply(hit, function(h) if (length(h)) countries$Country[h[1]] else NA_character_, character(1))
  res <- inside
  miss <- which(is.na(inside))
  if (length(miss)) {
    d <- sf::st_distance(sf::st_transform(pts[miss, ], na_albers), countries_m)
    near <- apply(d, 1, which.min); dist_km <- as.numeric(d[cbind(seq_along(miss), near)]) / 1000
    res[miss] <- ifelse(dist_km <= coast_buffer_km, paste0(countries$Country[near], " (coast)"), "outside")
  }
  out[ok] <- res
  out
}

country_from_field <- function(x) {
  x <- toupper(trimws(x))
  ifelse(is.na(x), NA_character_, ifelse(x %in% country_codes, "in region (country column)", "outside"))
}

decide <- function(region_xy, region_field) {
  status <- case_when(
    !is.na(region_xy) & region_xy != "outside" ~ "keep",
    !is.na(region_xy) & region_xy == "outside" ~ "outside",
    !is.na(region_field) & region_field != "outside" ~ "keep",
    !is.na(region_field) & region_field == "outside" ~ "outside",
    TRUE ~ "unverified")
  status
}

# FAUNMAP localities.
loc <- lapply(names(loc_raw), function(p) {
  df <- loc_raw[[p]]
  c_lat <- find_col(df, lat_opts); c_lon <- find_col(df, lon_opts); c_cty <- find_col(df, country_opts)
  df %>% mutate(
    .period = p, .row = row_number(),
    .mk = machine_key(col_or_na(df, find_col(df, machine_opts))),
    .ak = analysis_key(col_or_na(df, find_col(df, analysis_opts))),
    .site = clean_text(col_or_na(df, find_col(df, site_opts))),
    .lat = to_coord(col_or_na(df, c_lat), 90), .lon = to_coord(col_or_na(df, c_lon), 180),
    .region_xy = region_of(.lat, .lon),
    .region_field = country_from_field(col_or_na(df, c_cty)),
    .status = decide(.region_xy, .region_field))
})
names(loc) <- names(loc_raw)

# PBDB collections (one decision per collection).
p_site <- find_col(pbdb_raw, c("collectionname", "sitename"))
p_coll <- find_col(pbdb_raw, c("collectionno", "collectionnumber"))
p_lat <- find_col(pbdb_raw, lat_opts); p_lon <- find_col(pbdb_raw, lon_opts)
p_cty <- find_col(pbdb_raw, country_opts)
pbdb <- pbdb_raw %>% mutate(
  .row = row_number(),
  .site = clean_text(col_or_na(pbdb_raw, p_site)),
  .coll = coalesce(clean_text(col_or_na(pbdb_raw, p_coll)), .site),
  .lat = to_coord(col_or_na(pbdb_raw, p_lat), 90), .lon = to_coord(col_or_na(pbdb_raw, p_lon), 180),
  .cty = clean_text(col_or_na(pbdb_raw, p_cty)))
first_value <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else x[NA_integer_][1] }
pbdb_sites <- pbdb %>% group_by(.coll) %>%
  summarise(.site = first_value(.site), .lat = first_value(.lat), .lon = first_value(.lon),
            .cty = first_value(.cty), .groups = "drop")
pbdb_sites <- pbdb_sites %>%
  mutate(.region_xy = region_of(.lat, .lon), .region_field = country_from_field(.cty),
         .status = decide(.region_xy, .region_field))
pbdb <- left_join(pbdb, select(pbdb_sites, .coll, .status, .region_xy), by = ".coll")

report_region <- function(label, status) {
  cat(sprintf("  %-22s %6d sites: %6d inside | %5d outside | %5d unverified (%s)\n", label,
              length(status), sum(status == "keep"), sum(status == "outside"), sum(status == "unverified"),
              if (keep_unverified) "kept" else "removed"))
}
for (p in names(loc)) report_region(paste("FAUNMAP", p), loc[[p]]$.status)
report_region("PBDB collections", pbdb_sites$.status)
if (!is.na(p_cty)) {
  disagree <- sum(!is.na(pbdb_sites$.region_xy) & !is.na(pbdb_sites$.region_field) &
                    (pbdb_sites$.region_xy == "outside") != (pbdb_sites$.region_field == "outside"))
  if (disagree > 0) cat(sprintf("  NOTE: %d PBDB collections where coordinates and '%s' disagree (coordinates used)\n",
                                disagree, p_cty))
}

region_keep <- function(status) status == "keep" | (status == "unverified" & keep_unverified)

# =============================================================================
# 3. TAXA: REMOVE MARINE MAMMALS AND BATS
# =============================================================================

cat("\n=== 3. MARINE MAMMALS AND BATS ===\n")

genus_order <- NULL
if (!is.null(order_lookup_raw)) {
  lg <- find_col(order_lookup_raw, c("genus", "genusname")); lo <- find_col(order_lookup_raw, c("order", "ordername"))
  if (!is.na(lg) && !is.na(lo)) {
    genus_order <- order_lookup_raw %>%
      transmute(.genus = title_word(.data[[lg]]), .lk_order = title_word(.data[[lo]])) %>%
      filter(!is.na(.genus), !is.na(.lk_order)) %>% count(.genus, .lk_order) %>%
      arrange(.genus, desc(n)) %>% distinct(.genus, .keep_all = TRUE) %>% select(-n)
  }
}

lc <- function(x) tolower(trimws(x))
# Adds .taxon_reason (NA = keep) to a data frame of occurrences.
flag_taxa <- function(df) {
  c_order <- find_col(df, c("order", "ordername"))
  c_family <- find_col(df, c("family", "familyname"))
  c_genus <- find_col(df, c("genus", "genusname"))
  c_species <- find_col(df, c("species", "speciesname", "specificepithet"))
  c_name <- find_col(df, c("acceptedname", "identifiedname", "taxonname"))
  out <- df %>% mutate(
    .order = title_word(col_or_na(df, c_order)),
    .family = title_word(col_or_na(df, c_family)),
    .genus = title_word(col_or_na(df, c_genus)),
    .binomial = coalesce(
      as.character(ifelse(!is.na(.genus) & !is.na(clean_text(col_or_na(df, c_species))),
                          paste(.genus, sub("^.* ", "", lc(col_or_na(df, c_species)))), NA_character_)),
      as.character(clean_text(col_or_na(df, c_name)))))
  if (is.na(c_genus)) {
    out$.genus <- title_word(sub(" .*$", "", out$.binomial))
  }
  if (!is.null(genus_order)) {
    out <- out %>% left_join(genus_order, by = ".genus") %>%
      mutate(.order = coalesce(.order, .lk_order)) %>% select(-.lk_order)
  }
  out %>% mutate(.taxon_reason = case_when(
    lc(.order) %in% lc(bat_orders) | lc(.family) %in% lc(bat_families) ~ "bat (Chiroptera)",
    lc(.order) %in% lc(marine_orders) ~ paste0("marine mammal (order ", .order, ")"),
    lc(.family) %in% lc(marine_families) ~ paste0("marine mammal (family ", .family, ")"),
    lc(.genus) %in% lc(marine_genera) ~ paste0("marine mammal (genus ", .genus, ")"),
    lc(.binomial) %in% lc(marine_species) ~ paste0("marine mammal (", .binomial, ")"),
    TRUE ~ NA_character_))
}

# =============================================================================
# 4. APPLY THE CLEAN-UPS
# =============================================================================

cat("\n=== 4. APPLYING ===\n")

removed_records <- list(); removed_sites <- list(); summary_rows <- list()
clean_fauna <- list(); clean_loc <- list()

for (p in names(loc)) {
  L <- loc[[p]]
  F <- fauna_raw[[p]]
  f_mk <- machine_key(col_or_na(F, find_col(F, machine_opts)))
  f_ak <- analysis_key(col_or_na(F, find_col(F, analysis_opts)))
  F <- F %>% mutate(.row = row_number(), .mk = f_mk, .ak = f_ak)

  # Region status of each fauna record: from its locality, else its own coordinates.
  loc_status <- L %>% filter(!is.na(.mk), !is.na(.ak)) %>%
    group_by(.mk, .ak) %>% summarise(.loc_status = if (any(region_keep(.status))) "keep" else first(.status),
                                     .groups = "drop")
  F <- F %>% left_join(loc_status, by = c(".mk", ".ak"))
  if (any(is.na(F$.loc_status))) {
    own <- region_of(to_coord(col_or_na(F, find_col(F, lat_opts)), 90),
                     to_coord(col_or_na(F, find_col(F, lon_opts)), 180))
    F$.loc_status <- coalesce(F$.loc_status, decide(own, rep(NA_character_, nrow(F))))
  }
  F <- flag_taxa(F)
  F <- F %>% mutate(.reason = case_when(
    !region_keep(.loc_status) ~ chr_ifelse(.loc_status == "outside", "site outside Canada/USA/Mexico",
                                       "site location unverified"),
    !is.na(.taxon_reason) ~ .taxon_reason,
    TRUE ~ NA_character_))

  kept_F <- filter(F, is.na(.reason))
  removed_records[[p]] <- F %>% filter(!is.na(.reason)) %>%
    transmute(Database = "FAUNMAP", Period = p, Site = paste(.mk, .ak, sep = " / "),
              Genus = .genus, Order = .order, Family = .family, Reason = .reason)

  # Localities: region, then "no fauna left".
  has_fauna_before <- paste(L$.mk, L$.ak) %in% paste(F$.mk, F$.ak)
  has_fauna_after  <- paste(L$.mk, L$.ak) %in% paste(kept_F$.mk, kept_F$.ak)
  L <- L %>% mutate(.loc_reason = case_when(
    !region_keep(.status) ~ chr_ifelse(.status == "outside", "outside Canada/USA/Mexico", "location unverified"),
    drop_localities_without_fauna & !has_fauna_before ~ "no fauna records",
    drop_localities_without_fauna & !has_fauna_after ~ "no species left after removing marine mammals/bats",
    TRUE ~ NA_character_))
  kept_L <- filter(L, is.na(.loc_reason))
  removed_sites[[p]] <- L %>% filter(!is.na(.loc_reason)) %>%
    transmute(Database = "FAUNMAP", Period = p, Site = .site, Key = paste(.mk, .ak, sep = " / "),
              Latitude = .lat, Longitude = .lon, Region = .region_xy, Reason = .loc_reason)

  clean_fauna[[p]] <- kept_F %>% select(-starts_with("."))
  clean_loc[[p]] <- kept_L %>% select(-starts_with("."))
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Database = paste("FAUNMAP", p),
    sites_before = nrow(L), sites_outside_region = sum(!region_keep(L$.status)),
    sites_without_fauna = sum(L$.loc_reason %in% "no fauna records"),
    sites_emptied_by_taxa = sum(L$.loc_reason %in% "no species left after removing marine mammals/bats"),
    sites_after = nrow(kept_L),
    records_before = nrow(F), records_removed_region = sum(!region_keep(F$.loc_status)),
    records_removed_taxa = sum(region_keep(F$.loc_status) & !is.na(F$.taxon_reason)),
    records_after = nrow(kept_F))
}

pbdb <- flag_taxa(pbdb) %>% mutate(.reason = case_when(
  !region_keep(.status) ~ chr_ifelse(.status == "outside", "site outside Canada/USA/Mexico", "site location unverified"),
  !is.na(.taxon_reason) ~ .taxon_reason,
  TRUE ~ NA_character_))
kept_P <- filter(pbdb, is.na(.reason))
removed_records[["PBDB"]] <- pbdb %>% filter(!is.na(.reason)) %>%
  transmute(Database = "PBDB", Period = NA_character_, Site = .site, Genus = .genus, Order = .order,
            Family = .family, Reason = .reason)
coll_after <- unique(kept_P$.coll)
pbdb_sites <- pbdb_sites %>% mutate(.site_reason = case_when(
  !region_keep(.status) ~ chr_ifelse(.status == "outside", "outside Canada/USA/Mexico", "location unverified"),
  !.coll %in% coll_after ~ "no species left after removing marine mammals/bats",
  TRUE ~ NA_character_))
removed_sites[["PBDB"]] <- pbdb_sites %>% filter(!is.na(.site_reason)) %>%
  transmute(Database = "PBDB", Period = NA_character_, Site = .site, Key = .coll,
            Latitude = .lat, Longitude = .lon, Region = .region_xy, Reason = .site_reason)
summary_rows[[length(summary_rows) + 1]] <- data.frame(
  Database = "PBDB",
  sites_before = nrow(pbdb_sites), sites_outside_region = sum(!region_keep(pbdb_sites$.status)),
  sites_without_fauna = 0L,
  sites_emptied_by_taxa = sum(pbdb_sites$.site_reason %in% "no species left after removing marine mammals/bats"),
  sites_after = length(coll_after),
  records_before = nrow(pbdb), records_removed_region = sum(!region_keep(pbdb$.status)),
  records_removed_taxa = sum(region_keep(pbdb$.status) & !is.na(pbdb$.taxon_reason)),
  records_after = nrow(kept_P))

# =============================================================================
# 5. SAVE
# =============================================================================

cat("\n=== 5. SAVING ===\n")
for (p in periods) {
  save_csv(clean_loc[[p]], paste0(tolower(p), "_localities"))
  save_csv(clean_fauna[[p]], paste0(tolower(p), "_fauna"))
}
save_csv(kept_P %>% select(-starts_with(".")), "pbdb_clean")

removed_records <- bind_rows(removed_records)
removed_sites <- bind_rows(removed_sites)
summary_tab <- bind_rows(summary_rows)
save_csv(removed_records, "removed_records")
save_csv(removed_sites, "removed_sites")
save_csv(summary_tab, "cleaning_summary")

print(summary_tab, row.names = FALSE)
cat("\n  Removed occurrences by reason:\n")
print(as.data.frame(count(removed_records, Database, Reason, sort = TRUE)), row.names = FALSE)
cat("\n  Removed sites by reason:\n")
print(as.data.frame(count(removed_sites, Database, Reason, sort = TRUE)), row.names = FALSE)

cat("\n=== STEP 2b COMPLETE ===\n")
cat("Cleaned files in:", output_dir, "\n")
cat("Steps 3 and 6 read their FAUNMAP and PBDB inputs from this folder.\n")
