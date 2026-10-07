# =============================================================================
# STEP 3d: RESOLVE "sp." AND "cf." IDENTIFICATIONS FROM THE NEXT TIME BIN
#
# Occurrences identified only to genus ("Equus sp.") or with "cf."
# ("cf. Lepus californicus", "cf. Titanotylopus sp.") are given a species
# name using the species of the same genus found in the NEXT time bin:
#
#   - Bin 1 records are resolved with Bin 2, Bin 2 with Bin 3, ... Bin 5 with
#     Bin 6. A bin is never resolved with an older bin (Bin 2 never uses Bin 1).
#   - Only firmly identified species (no "sp.", no "cf.") of the next bin are
#     used as candidates; names given by this script are not used again.
#   - One candidate species in the genus -> that name.
#     Several candidates -> one of them at random (reproducible: random_seed).
#   - "cf. X": if species X itself is a candidate it is taken (set
#     'cf_prefer_named' to FALSE to treat cf. like sp.).
#   - No candidate (no species of the genus in the next bin), and every
#     sp./cf. record of the last bin (no younger bin): the record is REMOVED.
#
# Where candidates are looked for ('match_scope'):
#   "province" (default) - the same physiographic province (Physio_Province,
#                          Step 3) in the next bin. A record without a province
#                          (e.g. Mexico) cannot be resolved and is removed.
#   "all"                - anywhere in the study area in the next bin.
# Matching within the province avoids giving a record a species that was
# never found in its region, which would add shared species between regions
# and bias the analyses towards homogenization.
#
# FAUNMAP and PBDB records are treated alike. Records identified above genus
# level (no genus) are left unchanged. Other uncertainty flags in FAUNMAP's
# IDConfidence column ("?", "aff.") are left unchanged; add them to
# 'uncertain_flags' to resolve them in the same way.
#
# Inputs (Step 3c, Outputs/3_stages/): faunmap_fauna.csv, pbdb_occurrences.csv,
#   site_index.csv, time_bins.csv
# Outputs (Outputs/3d_resolved/):
#   faunmap_fauna.csv, pbdb_occurrences.csv   records after resolution (same
#       columns as Step 3c, names updated; the original names are kept in
#       *_original columns and Resolution says what happened to each record)
#   site_index.csv, time_bins.csv             sites that still have records
#   resolution_log.csv      every sp./cf. record: candidates and result
#   removed_records.csv     records removed, with the reason
#   emptied_sites.csv       sites left without any occurrence, and why
#   resolution_summary.csv  counts per time bin
#   resolution_summary.png/.pdf
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 3c.
# Needs: dplyr, tidyr, ggplot2.
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(ggplot2)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

input_dir  <- file.path("Outputs", "3_stages")
output_dir <- file.path(work_dir, "Outputs", "3d_resolved")

match_scope     <- "province"   # "province" or "all" (see above)
cf_prefer_named <- TRUE         # "cf. X": take X when X is among the candidates
uncertain_flags <- c("cf")      # IDConfidence values resolved like sp. (e.g. c("cf", "aff", "?"))
random_seed     <- 2024         # makes the random choices repeatable

# =============================================================================
# HELPERS
# =============================================================================

read_text_csv <- function(name) {
  path <- file.path(work_dir, input_dir, name)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun Step 3c first.")
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
                 na.strings = c("", "NA"), encoding = "UTF-8")
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))
  cat(sprintf("  %-28s %7d rows\n", name, nrow(df)))
  df
}
save_csv <- function(df, name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, paste0(name, ".csv"))
  ok <- tryCatch({ write.csv(df, path, row.names = FALSE, na = ""); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) stop("Could not write ", path, "\n  Close it (e.g. in Excel) and run again.")
  cat(sprintf("  %-28s %7d rows -> %s\n", name, nrow(df), path))
}

clean <- function(x) { x <- trimws(gsub("\\s+", " ", as.character(x))); x[x == ""] <- NA; x }
cap1  <- function(x) ifelse(is.na(x), NA, paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2))))
flag_of <- function(x) gsub("[^a-z?]", "", tolower(coalesce(x, "")))   # " CF." -> "cf"
is_sp_word <- function(x) is.na(x) | grepl("^(sp|spp|indet|indeterminate)\\.?$", tolower(x))

# =============================================================================
# 1. READ THE STEP 3c OCCURRENCES
# =============================================================================

cat("=== 1. READING STEP 3c OUTPUTS ===\n")
fauna      <- read_text_csv("faunmap_fauna.csv")
pbdb       <- read_text_csv("pbdb_occurrences.csv")
site_index <- read_text_csv("site_index.csv")
time_bins  <- read_text_csv("time_bins.csv")
last_bin   <- max(as.integer(time_bins$Bin_Number))

if (match_scope == "province" && !"Physio_Province" %in% c(names(fauna), names(pbdb))) {
  stop("No Physio_Province column: run Step 3 (spatial binning) and Step 3c first, or set match_scope <- \"all\".")
}

# One table of all occurrences with genus, species epithet and uncertainty.
fa <- fauna %>% transmute(
  Source = "FAUNMAP", Row = row_number(), Site_Key,
  Bin = as.integer(Stage_Number), Province = clean(Physio_Province),
  Genus = cap1(sub("^(cf|aff)\\.?\\s+", "", clean(Genus), ignore.case = TRUE)),
  Epithet_raw = clean(Species),
  Flag = flag_of(IDConfidence))
# A "cf." written inside the species field counts as well.
fa <- fa %>% mutate(Flag = ifelse(grepl("(^|\\s)cf\\.?(\\s|$)", tolower(coalesce(Epithet_raw, ""))), "cf", Flag),
                    Epithet_raw = clean(sub("(^|\\s)cf\\.?(\\s|$)", " ", Epithet_raw, ignore.case = TRUE)))

pb_name <- clean(pbdb[[intersect(c("accepted_name", "identified_name"), names(pbdb))[1]]])
pb <- data.frame(
  Source = "PBDB", Row = seq_len(nrow(pbdb)), Site_Key = pbdb$Site_Key,
  Bin = as.integer(pbdb$Stage_Number),
  Province = if ("Physio_Province" %in% names(pbdb)) clean(pbdb$Physio_Province) else NA_character_,
  Name = pb_name, stringsAsFactors = FALSE) %>%
  mutate(Flag = ifelse(grepl("(^|\\s)cf\\.?(\\s|$)", tolower(coalesce(Name, ""))), "cf", ""),
         Name = clean(gsub("(^|\\s)(cf|aff)\\.?(\\s|$)", " ", Name, ignore.case = TRUE)),
         Genus = cap1(sub(" .*$", "", Name)),
         Epithet_raw = ifelse(grepl(" ", coalesce(Name, "")), sub("^\\S+\\s+", "", Name), NA)) %>%
  select(-Name)
if ("genus" %in% names(pbdb)) pb$Genus <- coalesce(pb$Genus, cap1(clean(pbdb$genus)))

occ <- bind_rows(fa, pb) %>%
  mutate(Epithet_raw = ifelse(is.na(Epithet_raw), NA, sub("\\s.*$", "", Epithet_raw)),   # keep first word
         Is_sp = !is.na(Genus) & is_sp_word(Epithet_raw),
         Is_cf = !is.na(Genus) & Flag %in% flag_of(uncertain_flags),
         Epithet = ifelse(Is_sp, NA, tolower(Epithet_raw)),
         Type = case_when(is.na(Genus) ~ "above genus (unchanged)",
                          Is_cf & !is.na(Epithet) ~ "cf. species",
                          Is_cf ~ "cf. sp.",
                          Is_sp ~ "sp.",
                          TRUE ~ "determinate"),
         Scope = if (match_scope == "province") Province else "all")

cat("\n  Occurrences by identification:\n")
print(as.data.frame(count(occ, Type)), row.names = FALSE)
other_flags <- setdiff(unique(fa$Flag[fa$Flag != ""]), flag_of(uncertain_flags))
if (length(other_flags)) {
  cat(sprintf("  Left unchanged: %d FAUNMAP records flagged %s in IDConfidence\n",
              sum(fa$Flag %in% other_flags), paste(sQuote(other_flags, FALSE), collapse = ", ")))
}

# =============================================================================
# 2. RESOLVE BIN BY BIN WITH THE NEXT BIN
# =============================================================================

cat(sprintf("\n=== 2. RESOLVING (candidates from the next bin, scope: %s) ===\n", match_scope))

# Candidate species: firmly identified species of each bin, genus and scope.
reference <- occ %>% filter(Type == "determinate") %>%
  distinct(Bin, Scope, Genus, Epithet)

set.seed(random_seed)
todo <- occ %>% filter(Type %in% c("sp.", "cf. sp.", "cf. species")) %>% arrange(Bin, Source, Row)
log_rows <- vector("list", nrow(todo))
for (i in seq_len(nrow(todo))) {
  r <- todo[i, ]
  cand <- if (r$Bin >= last_bin || is.na(r$Scope)) character(0) else
    sort(reference$Epithet[reference$Bin == r$Bin + 1 & reference$Scope %in% r$Scope &
                             reference$Genus == r$Genus])
  assigned <- NA_character_
  how <- NA_character_
  if (length(cand) > 0) {
    if (cf_prefer_named && r$Type == "cf. species" && r$Epithet %in% cand) {
      assigned <- r$Epithet; how <- "cf. species found in next bin"
    } else if (length(cand) == 1) {
      assigned <- cand; how <- "only species of the genus in next bin"
    } else {
      assigned <- cand[sample.int(length(cand), 1)]; how <- sprintf("random choice of %d species", length(cand))
    }
  }
  reason <- if (!is.na(assigned)) NA_character_ else if (r$Bin >= last_bin) "last time bin (no younger bin)" else
    if (is.na(r$Scope)) "no province (cannot match)" else "no species of the genus in next bin"
  log_rows[[i]] <- data.frame(
    Source = r$Source, Row = r$Row, Site_Key = r$Site_Key, Bin = r$Bin, Province = r$Province,
    Original = trimws(paste(if (r$Is_cf) "cf." else "", r$Genus, coalesce(r$Epithet, "sp."))),
    Type = r$Type, Reference_Bin = if (r$Bin < last_bin) r$Bin + 1L else NA_integer_,
    Candidates = paste(cand, collapse = "; "), n_candidates = length(cand),
    Assigned = ifelse(is.na(assigned), NA, paste(r$Genus, assigned)),
    Result = ifelse(is.na(assigned), "removed", "resolved"),
    How = coalesce(how, reason), stringsAsFactors = FALSE)
}
resolution_log <- bind_rows(log_rows)

per_bin <- resolution_log %>% count(Bin, Result) %>%
  tidyr::pivot_wider(names_from = Result, values_from = n, values_fill = 0)
cat("\n  sp./cf. records per bin:\n")
print(as.data.frame(per_bin), row.names = FALSE)
cat("\n  How records were resolved or why they were removed:\n")
print(as.data.frame(count(resolution_log, Result, How, sort = TRUE)), row.names = FALSE)

# =============================================================================
# 3. APPLY THE RESULTS TO THE OCCURRENCE TABLES
# =============================================================================

cat("\n=== 3. UPDATING THE OCCURRENCE TABLES ===\n")

status_of <- function(src, n) {
  st <- rep("unchanged", n)
  lg <- resolution_log[resolution_log$Source == src, ]
  st[lg$Row] <- ifelse(lg$Result == "resolved", paste0("resolved (", lg$Type, " -> ", lg$How, ")"),
                       paste0("removed (", lg$How, ")"))
  st
}

# FAUNMAP: Species gets the assigned epithet; the cf. flag is cleared.
fauna_out <- fauna
fauna_out$Genus_original <- fauna$Genus
fauna_out$Species_original <- fauna$Species
fauna_out$IDConfidence_original <- fauna$IDConfidence
lg <- resolution_log %>% filter(Source == "FAUNMAP", Result == "resolved")
fauna_out$Genus[lg$Row]        <- sub(" .*$", "", lg$Assigned)
fauna_out$Species[lg$Row]      <- sub("^\\S+\\s+", "", lg$Assigned)
fauna_out$IDConfidence[lg$Row] <- NA
fauna_out$Resolution <- status_of("FAUNMAP", nrow(fauna))

# PBDB: the name column gets the new binomial.
name_col <- intersect(c("accepted_name", "identified_name"), names(pbdb))[1]
pbdb_out <- pbdb
pbdb_out[[paste0(name_col, "_original")]] <- pbdb[[name_col]]
lg <- resolution_log %>% filter(Source == "PBDB", Result == "resolved")
pbdb_out[[name_col]][lg$Row] <- lg$Assigned
pbdb_out$Resolution <- status_of("PBDB", nrow(pbdb))

removed_records <- bind_rows(
  fauna_out %>% filter(startsWith(Resolution, "removed")) %>%
    transmute(Source = "FAUNMAP", Site_Key, Bin = Stage_Number, Province = Physio_Province,
              Name = trimws(paste(coalesce(IDConfidence, ""), Genus, coalesce(Species, ""))), Resolution),
  pbdb_out %>% filter(startsWith(Resolution, "removed")) %>%
    transmute(Source = "PBDB", Site_Key, Bin = Stage_Number,
              Province = if ("Physio_Province" %in% names(pbdb_out)) Physio_Province else NA_character_,
              Name = .data[[paste0(name_col, "_original")]], Resolution))
fauna_out <- filter(fauna_out, !startsWith(Resolution, "removed"))
pbdb_out  <- filter(pbdb_out, !startsWith(Resolution, "removed"))

# Sites that still have at least one record.
kept_sites <- unique(c(paste(fauna_out$Site_Key, fauna_out$Stage_Number),
                       paste(pbdb_out$Site_Key, pbdb_out$Stage_Number)))
site_index_out <- site_index %>% filter(paste(Site_Key, Stage_Number) %in% kept_sites)
cat(sprintf("  FAUNMAP %d -> %d records | PBDB %d -> %d records | sites %d -> %d\n",
            nrow(fauna), nrow(fauna_out), nrow(pbdb), nrow(pbdb_out),
            nrow(site_index), nrow(site_index_out)))

# Sites left without any occurrence, and why. Sites that had no records even
# before this step (e.g. a duplicate locality row) are listed separately.
site_key <- function(d) paste(d$Site_Key, d$Stage_Number)
had_records <- unique(c(site_key(fauna), site_key(pbdb)))
why <- resolution_log %>% filter(Result == "removed") %>%
  group_by(k = paste(Site_Key, Bin)) %>%
  summarise(Records_removed = n(), Reason = paste(sort(unique(How)), collapse = " + "), .groups = "drop")
emptied_sites <- site_index %>%
  mutate(k = site_key(site_index)) %>%
  filter(!k %in% kept_sites) %>%
  left_join(why, by = "k") %>%
  mutate(Reason = case_when(!k %in% had_records ~ "had no records before Step 3d",
                            TRUE ~ Reason)) %>%
  select(-k) %>%
  arrange(as.integer(Stage_Number), Reason)
cat(sprintf("  Sites left without any occurrence: %d (of %d)\n",
            sum(emptied_sites$Reason != "had no records before Step 3d"), nrow(site_index)))
print(as.data.frame(emptied_sites %>% count(Bin = as.integer(Stage_Number), Reason)), row.names = FALSE)

resolution_summary <- occ %>%
  filter(Type != "above genus (unchanged)") %>%
  mutate(Result = case_when(Type == "determinate" ~ "determinate",
                            TRUE ~ resolution_log$Result[match(paste(Source, Row),
                                                               paste(resolution_log$Source, resolution_log$Row))])) %>%
  count(Bin, Result) %>%
  tidyr::pivot_wider(names_from = Result, values_from = n, values_fill = 0) %>%
  arrange(Bin)
for (col in c("determinate", "resolved", "removed")) if (!col %in% names(resolution_summary)) resolution_summary[[col]] <- 0L
resolution_summary <- resolution_summary %>%
  mutate(sp_cf_records = resolved + removed,
         percent_resolved = round(100 * resolved / pmax(sp_cf_records, 1), 1)) %>%
  select(Bin, determinate, sp_cf_records, resolved, removed, percent_resolved)

# =============================================================================
# 4. SAVE
# =============================================================================

cat("\n=== 4. SAVING ===\n")
save_csv(fauna_out, "faunmap_fauna")
save_csv(pbdb_out, "pbdb_occurrences")
save_csv(site_index_out, "site_index")
save_csv(time_bins, "time_bins")
save_csv(resolution_log, "resolution_log")
save_csv(removed_records, "removed_records")
save_csv(emptied_sites, "emptied_sites")
save_csv(resolution_summary, "resolution_summary")
print(as.data.frame(resolution_summary), row.names = FALSE)

# Figure: species-level records per bin, by outcome.
fig <- resolution_summary %>%
  select(Bin, determinate, resolved, removed) %>%
  tidyr::pivot_longer(-Bin, names_to = "Outcome", values_to = "n") %>%
  mutate(Outcome = factor(Outcome, levels = c("removed", "resolved", "determinate"),
                          labels = c("sp./cf. removed (no congener in next bin)",
                                     "sp./cf. resolved from next bin",
                                     "Identified to species")),
         Bin = factor(paste("Bin", Bin), levels = paste("Bin", sort(unique(Bin)))))
p <- ggplot(fig, aes(Bin, n, fill = Outcome)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("#D9D9D9", "#E69F00", "#2A78D6"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Resolving sp. and cf. records with the next time bin",
       subtitle = sprintf("Candidates: species of the same genus in the next bin (%s). The last bin has no younger bin.",
                          if (match_scope == "province") "same physiographic province" else "whole study area"),
       x = NULL, y = "Occurrences") +
  theme_minimal(base_size = 15) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.subtitle = element_text(colour = "grey35", size = 11), legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA))
for (ext in c("png", "pdf")) {
  file <- file.path(output_dir, paste0("resolution_summary.", ext))
  if (ext == "png") ggsave(file, p, width = 11, height = 6.5, dpi = 300) else ggsave(file, p, width = 11, height = 6.5)
}
cat("  Figure:", file.path(output_dir, "resolution_summary.png"), "\n")

cat("\n=== STEP 3d COMPLETE ===\n")
cat("Objects in your Environment: fauna_out, pbdb_out, resolution_log, resolution_summary\n")
