# =============================================================================
# STEP 6: SITE x SPECIES MATRICES AND HEATMAPS BY GEOLOGICAL STAGE
#
# Builds one presence/absence matrix (sites x species) per geological stage
# from FAUNMAP and PBDB, and draws each matrix as a heatmap (PDF + PNG).
#
# Inputs:
#   Order Filtered Fauna.csv            FAUNMAP fauna (Machine Number,
#                                       Analysis Unit, Genus, Species, Order)
#   Outputs/3_stages/faunmap_localities.csv   FAUNMAP site, stage, lat/long (Step 3)
#   Outputs/3_stages/pbdb_occurrences.csv     PBDB occurrences with stage,
#                                             lat/long and taxonomy (Step 3)
#   Outputs/3_stages/site_index.csv           all sites per stage (for the
#                                             site accounting only)
#
# Rules:
#   - Sites: FAUNMAP SiteName, PBDB collection_name (one row per site per stage).
#   - FAUNMAP fauna get the stage, SiteName and coordinates of their locality,
#     matched on Machine Number + Analysis Unit (1234 and 1234.00 are equal).
#     If a pair has locality rows in several stages, its fauna are placed in
#     every one of those stages (as in the original script).
#   - Each genus gets one order (the order most of its records use).
#   - Only the orders in `orders_to_include` are kept.
#   - "sp." / "cf." / "aff." / "indet." identifications are resolved to a
#     species of the same genus recorded in the same stage; if none, the next
#     (younger) stage; otherwise the record is omitted. Every decision is saved.
#   - NO minimum-species (sampling) filter is applied.
#
# Outputs (Outputs/6_matrices/):
#   matrix_<stage>.xlsx / .csv      SiteName, Database, Latitude, Longitude, species...
#   heatmap_<stage>.pdf / .png      the matrix drawn as a heatmap
#   all_records_final.xlsx / .csv   every retained record after resolution
#   spcf_resolution_log.csv         what happened to each sp./cf. record
#   genus_order_conflicts.csv       genera recorded under more than one order
#   site_attrition.csv              every Step 3 site, and why it is or is not
#                                   in the matrices
#   summary.xlsx                    stage, order, resolution and site summaries
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 3.
# Needs: dplyr, tidyr, ggplot2, writexl.
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "writexl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

faunmap_fauna_file  <- "Order Filtered Fauna.csv"
faunmap_sites_file  <- file.path("Outputs", "3_stages", "faunmap_localities.csv")
pbdb_file           <- file.path("Outputs", "3_stages", "pbdb_occurrences.csv")
site_index_file     <- file.path("Outputs", "3_stages", "site_index.csv")

output_dir <- file.path(work_dir, "Outputs", "6_matrices")

orders_to_include <- c("Artiodactyla", "Perissodactyla", "Xenarthra", "Insectivora",
                       "Proboscidea", "Rodentia", "Lagomorpha")

# As in your original script, Insectivora records are grouped with Xenarthra.
# Set to NA to keep Insectivora as its own order.
merge_insectivora_into <- "Xenarthra"

# Order in which orders are stacked on the heatmaps (top to bottom).
order_levels <- c("Artiodactyla", "Perissodactyla", "Xenarthra", "Insectivora",
                  "Proboscidea", "Rodentia", "Lagomorpha")

stage_names <- c("Zanclean", "Piacenzian", "Gelasian", "Calabrian", "Chibanian")
stage_labels <- c(
  "Bin 1 - Zanclean (4.700-3.600 Ma)",
  "Bin 2 - Piacenzian (3.600-2.580 Ma)",
  "Bin 3 - Gelasian (2.580-1.800 Ma)",
  "Bin 4 - Calabrian (1.800-0.7741 Ma)",
  "Bin 5 - Chibanian (0.7741-0.129 Ma)"
)

present_colour <- "#2D7A3D"
absent_colour  <- "#F7F7F7"
png_dpi        <- 200

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================

read_input <- function(file) {
  path <- file.path(work_dir, file)
  if (!file.exists(path)) stop("File not found:\n  ", path)
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                 colClasses = "character", na.strings = c("", "NA", "N/A"))
  # Remove hidden byte-order marks / stray spaces and empty column names.
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))
  df <- df[, names(df) != "" & !is.na(names(df)), drop = FALSE]
  cat(sprintf("  %-45s %7d rows, %3d columns\n", file, nrow(df), ncol(df)))
  df
}

simplify_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# First column whose simplified name matches one of `options` (in order).
find_col <- function(df, options, label, required = TRUE) {
  simple <- simplify_name(names(df))
  for (o in options) {
    hit <- names(df)[simple == o]
    if (length(hit)) return(hit[1])
  }
  if (required) {
    stop("No ", label, " column found.\nColumns: ", paste(names(df), collapse = " | "))
  }
  NA_character_
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[!is.na(x) & x == ""] <- NA_character_
  x
}

# 1234.00 and 1234 become the same key.
machine_key <- function(x) sub("^([0-9]+)\\.0+$", "\\1", clean_text(x))

# "RODENTIA" / "rodentia" -> "Rodentia"
title_word <- function(x) {
  x <- clean_text(x)
  ifelse(is.na(x), NA_character_,
         paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2))))
}

# Species epithet: drops a repeated genus ("Equus simplicidens" -> "simplicidens"),
# and turns empty or genus-only entries into "sp.".
clean_species <- function(genus, species) {
  species <- gsub("\\s+", " ", clean_text(species))
  same_start <- !is.na(species) & !is.na(genus) &
    tolower(sub("\\s.*$", "", species)) == tolower(genus)
  species[same_start] <- clean_text(sub("^\\S+\\s*", "", species[same_start]))
  species <- tolower(species)
  species[is.na(species)] <- "sp."
  species
}

# Unresolved identifications: sp., spp., cf., aff., indet., "?" etc.
is_spcf <- function(x) {
  x <- tolower(trimws(as.character(x)))
  is.na(x) | x == "" |
    grepl("(^|\\s)(sp|spp|cf|aff|indet|indeterminate)\\.?(\\s|$)", x) |
    grepl("\\?", x)
}

save_both <- function(df, name) {
  write.csv(df, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE, na = "")
  write_xlsx(df, file.path(output_dir, paste0(name, ".xlsx")))
}

# =============================================================================
# 1. READ THE FILES
# =============================================================================

cat("=== 1. READING FILES ===\n")

needed <- c(faunmap_fauna_file, faunmap_sites_file, pbdb_file, site_index_file)
found <- file.exists(file.path(work_dir, needed))
for (i in seq_along(needed)) {
  cat(if (found[i]) "  [found]   " else "  [MISSING] ", needed[i], "\n", sep = "")
}
if (!found[1]) {
  cat("\n  Files in the folder whose name contains 'fauna':\n")
  cat(paste0("    ", list.files(work_dir, pattern = "fauna", ignore.case = TRUE)), sep = "\n")
}
if (!all(found)) stop("Some input files were not found (run Step 3 first for the Outputs files).")

faunmap_fauna_raw <- read_input(faunmap_fauna_file)
faunmap_sites_raw <- read_input(faunmap_sites_file)
pbdb_raw          <- read_input(pbdb_file)
site_index_raw    <- read_input(site_index_file)

# =============================================================================
# 2. FAUNMAP: LINK FAUNA TO SITES, STAGES AND COORDINATES
# =============================================================================

cat("\n=== 2. FAUNMAP FAUNA ===\n")

f_machine  <- find_col(faunmap_fauna_raw, c("machinenumber", "machineno", "machine"), "Machine Number")
f_analysis <- find_col(faunmap_fauna_raw, c("analysisunit"), "Analysis Unit")
f_genus    <- find_col(faunmap_fauna_raw, c("genus", "genusname"), "Genus")
f_species  <- find_col(faunmap_fauna_raw, c("species", "speciesname", "specificepithet"), "Species")
f_order    <- find_col(faunmap_fauna_raw, c("order", "ordername"), "Order")
f_period   <- find_col(faunmap_fauna_raw, c("faunmapperiod", "period", "nalma"), "period",
                       required = FALSE)

cat(sprintf("  Columns used: Machine='%s', Analysis Unit='%s', Genus='%s', Species='%s', Order='%s'%s\n",
            f_machine, f_analysis, f_genus, f_species, f_order,
            if (is.na(f_period)) "" else paste0(", Period='", f_period, "'")))

# One row per locality pair (period + Machine Number + Analysis Unit).
site_lookup <- faunmap_sites_raw %>%
  transmute(
    FAUNMAP_Period,
    .mk = Machine_Key,
    .ak = clean_text(Analysis_Key),
    SiteName = clean_text(SiteName_Std),
    Latitude = suppressWarnings(as.numeric(Latitude)),
    Longitude = suppressWarnings(as.numeric(Longitude)),
    Stage_Number = as.integer(Stage_Number)
  ) %>%
  filter(!is.na(.mk), !is.na(.ak), !is.na(SiteName), Stage_Number %in% 1:5)

join_by_period <- !is.na(f_period)
key_cols <- if (join_by_period) c("FAUNMAP_Period", ".mk", ".ak") else c(".mk", ".ak")

# A pair can have several locality rows (e.g. in both the Blancan and the
# Irvingtonian file, or with different ages). Every row is kept, so the
# fauna appear in every stage their locality appears in.
site_lookup <- distinct(site_lookup)
pair_stages <- site_lookup %>%
  group_by(across(all_of(key_cols))) %>%
  summarise(n_stages = n_distinct(Stage_Number), .groups = "drop")
n_multi <- sum(pair_stages$n_stages > 1)
cat(sprintf("  Machine Number + Analysis Unit pairs whose localities fall in more than one\n"),
    sprintf("  stage: %d (their fauna are placed in each of those stages)\n", n_multi), sep = "")

faunmap_records <- faunmap_fauna_raw %>%
  mutate(
    FAUNMAP_Period = if (join_by_period) title_word(.data[[f_period]]) else NA_character_,
    .mk = machine_key(.data[[f_machine]]),
    .ak = clean_text(.data[[f_analysis]]),
    Genus = title_word(.data[[f_genus]]),
    Order = title_word(.data[[f_order]])
  ) %>%
  mutate(Species = clean_species(Genus, .data[[f_species]]))

n_in <- nrow(faunmap_records)
faunmap_records <- faunmap_records %>%
  filter(!is.na(.mk), !is.na(.ak))
faunmap_linked <- faunmap_records %>%
  inner_join(site_lookup, by = key_cols, suffix = c(".fauna", ""),
             relationship = "many-to-many") %>%
  transmute(Database = "FAUNMAP", SiteName, Latitude, Longitude,
            Stage_Number, Genus, Species, Order)

n_unlinked <- nrow(anti_join(faunmap_records, site_lookup, by = key_cols))
cat(sprintf("  %d fauna records | %d with Machine Number + Analysis Unit\n",
            n_in, nrow(faunmap_records)))
cat(sprintf("  %d records match no staged locality (pair not in the locality files,\n", n_unlinked),
    "   or its locality has no stage 1-5)\n", sep = "")
cat(sprintf("  %d record-stage rows after linking\n", nrow(faunmap_linked)))

# =============================================================================
# 3. PBDB: OCCURRENCES ALREADY CARRY STAGE AND COORDINATES
# =============================================================================

cat("\n=== 3. PBDB OCCURRENCES ===\n")

p_order <- find_col(pbdb_raw, c("order", "ordername"), "order")
p_genus <- find_col(pbdb_raw, c("genus", "genusname"), "genus", required = FALSE)
p_ident <- find_col(pbdb_raw, c("identifiedname", "acceptedname", "taxonname"),
                    "identified_name / accepted_name")
p_species <- find_col(pbdb_raw, c("species", "speciesname"), "species", required = FALSE)

cat(sprintf("  Columns used: order='%s', identification='%s', genus=%s, species=%s\n",
            p_order, p_ident,
            if (is.na(p_genus)) "taken from the identification" else paste0("'", p_genus, "'"),
            if (is.na(p_species)) "taken from the identification" else paste0("'", p_species, "'")))

pbdb_records <- pbdb_raw %>%
  mutate(
    Ident = gsub("\\s+", " ", clean_text(.data[[p_ident]])),
    Genus = title_word(if (!is.na(p_genus)) .data[[p_genus]] else sub("\\s.*$", "", Ident)),
    # Text after the genus in the identification, e.g. "cf. simplicidens".
    Species_raw = if (!is.na(p_species)) .data[[p_species]] else
      ifelse(grepl("\\s", Ident), sub("^\\S+\\s+", "", Ident), NA_character_),
    Order = title_word(.data[[p_order]])
  ) %>%
  mutate(Species = clean_species(Genus, Species_raw)) %>%
  transmute(Database = "PBDB",
            SiteName = clean_text(collection_name),
            Latitude = suppressWarnings(as.numeric(Latitude)),
            Longitude = suppressWarnings(as.numeric(Longitude)),
            Stage_Number = as.integer(Stage_Number),
            Genus, Species, Order) %>%
  filter(!is.na(SiteName), Stage_Number %in% 1:5)

cat(sprintf("  %d occurrences | %d with a collection and stage | %d of these lack a genus\n",
            nrow(pbdb_raw), nrow(pbdb_records), sum(is.na(pbdb_records$Genus))))
cat("  (PBDB leaves genus empty when a fossil is identified only to family or higher)\n")

# =============================================================================
# 4. COMBINE AND KEEP THE SELECTED ORDERS
# =============================================================================

cat("\n=== 4. COMBINING AND FILTERING ORDERS ===\n")

all_raw <- bind_rows(faunmap_linked, pbdb_records)

# Site sets after each step, for the site accounting in section 7.
site_set <- function(df) distinct(df, Stage_Number, Database, SiteName)
sites_linked <- site_set(all_raw)

n_no_genus <- sum(is.na(all_raw$Genus))
all_raw <- filter(all_raw, !is.na(Genus))
sites_genus <- site_set(all_raw)
cat(sprintf("  Records without a genus removed: %d\n", n_no_genus))

# One order per genus. A genus recorded under different orders (e.g. FAUNMAP
# and PBDB disagree, or a typo) is given the order most of its records use;
# ties go to the order listed first in `order_levels`, then alphabetically.
genus_orders <- all_raw %>%
  filter(!is.na(Order)) %>%
  count(Genus, Order, name = "n_records") %>%
  group_by(Genus) %>%
  mutate(n_orders = n()) %>%
  arrange(Genus, desc(n_records), match(Order, order_levels), Order) %>%
  mutate(Assigned_Order = first(Order)) %>%
  ungroup()

order_conflicts <- genus_orders %>%
  filter(n_orders > 1) %>%
  select(Genus, Recorded_Order = Order, n_records, Assigned_Order)
write.csv(order_conflicts, file.path(output_dir, "genus_order_conflicts.csv"),
          row.names = FALSE, na = "")
cat(sprintf("  Genera recorded under more than one order: %d (assigned their majority order;\n",
            n_distinct(order_conflicts$Genus)),
    "   see genus_order_conflicts.csv)\n", sep = "")

all_raw <- all_raw %>%
  left_join(distinct(genus_orders, Genus, Assigned_Order), by = "Genus") %>%
  mutate(Order = coalesce(Assigned_Order, Order)) %>%
  select(-Assigned_Order)

order_counts <- count(all_raw, Order, sort = TRUE)
cat("  Orders present (records):\n")
print(as.data.frame(head(order_counts, 25)), row.names = FALSE)

all_raw <- all_raw %>%
  filter(Order %in% orders_to_include) %>%
  mutate(Order = if (!is.na(merge_insectivora_into))
           if_else(Order == "Insectivora", merge_insectivora_into, Order) else Order,
         GenusSpecies = paste(Genus, Species),
         is_spcf_flag = is_spcf(Species))

sites_orders <- site_set(all_raw)

cat(sprintf("\n  Records in the selected orders: %d (FAUNMAP %d, PBDB %d)\n",
            nrow(all_raw), sum(all_raw$Database == "FAUNMAP"), sum(all_raw$Database == "PBDB")))
cat(sprintf("  Sites: %d | genera: %d\n", n_distinct(all_raw$SiteName), n_distinct(all_raw$Genus)))

# =============================================================================
# 5. RESOLVE sp. / cf. IDENTIFICATIONS
#
# A sp./cf. record takes the alphabetically first species of the same genus
# recorded (as a firm identification) in the same stage; failing that, in
# the next younger stage; otherwise the record is omitted.
# =============================================================================

cat("\n=== 5. RESOLVING sp. / cf. IDENTIFICATIONS ===\n")

firm <- filter(all_raw, !is_spcf_flag)
spcf <- filter(all_raw, is_spcf_flag)

candidates <- firm %>%
  distinct(Genus, Stage_Number, Species) %>%
  group_by(Genus, Stage_Number) %>%
  summarise(Candidate = sort(Species)[1], .groups = "drop")

spcf_log <- spcf %>%
  left_join(rename(candidates, Same_Stage = Candidate), by = c("Genus", "Stage_Number")) %>%
  left_join(candidates %>% mutate(Stage_Number = Stage_Number - 1L) %>%
              rename(Next_Stage = Candidate),
            by = c("Genus", "Stage_Number")) %>%
  mutate(
    Resolved_Species = coalesce(Same_Stage, Next_Stage),
    Resolution = case_when(
      !is.na(Same_Stage) ~ "Resolved: same stage",
      !is.na(Next_Stage) ~ "Resolved: next stage",
      TRUE ~ "Omitted: no species of this genus in this or the next stage"
    )
  )

resolved <- spcf_log %>%
  filter(!is.na(Resolved_Species)) %>%
  mutate(Species = Resolved_Species, GenusSpecies = paste(Genus, Species),
         is_spcf_flag = FALSE) %>%
  select(all_of(names(firm)))

all_final <- bind_rows(firm, resolved) %>%
  distinct(Stage_Number, SiteName, GenusSpecies, .keep_all = TRUE) %>%
  select(-is_spcf_flag) %>%
  mutate(Stage = stage_names[Stage_Number]) %>%
  arrange(Stage_Number, SiteName, Order, GenusSpecies)

resolution_counts <- count(spcf_log, Resolution, name = "n_records")
cat(sprintf("  Firm identifications: %d | sp./cf. records: %d\n", nrow(firm), nrow(spcf)))
print(as.data.frame(resolution_counts), row.names = FALSE)
cat(sprintf("  Records retained (one per site, species and stage): %d\n", nrow(all_final)))
cat("  No minimum-species filter applied.\n")

save_both(all_final, "all_records_final")
write.csv(select(spcf_log, Database, SiteName, Stage_Number, Order, Genus,
                 Original_Species = Species, Resolution, Resolved_Species),
          file.path(output_dir, "spcf_resolution_log.csv"), row.names = FALSE, na = "")

# =============================================================================
# 6. MATRICES AND HEATMAPS, ONE PER STAGE
# =============================================================================

cat("\n=== 6. MATRICES AND HEATMAPS ===\n")

for (s in 1:5) {
  stage <- stage_names[s]
  bin_data <- filter(all_final, Stage_Number == s)
  cat(sprintf("  %-11s ", stage))
  if (nrow(bin_data) == 0) {
    cat("no records - skipped\n")
    next
  }

  # Site information: first available coordinates per site.
  site_info <- bin_data %>%
    group_by(SiteName) %>%
    summarise(Database = paste(sort(unique(Database)), collapse = "+"),
              Latitude = first(Latitude[!is.na(Latitude)]),
              Longitude = first(Longitude[!is.na(Longitude)]),
              n_species = n_distinct(GenusSpecies),
              .groups = "drop") %>%
    arrange(desc(n_species), SiteName)

  # Species order: by taxonomic order, then by number of sites, then name.
  species_info <- bin_data %>%
    group_by(GenusSpecies) %>%
    summarise(Order = first(Order), n_sites = n_distinct(SiteName), .groups = "drop") %>%
    mutate(Order = factor(Order, levels = order_levels)) %>%
    arrange(Order, desc(n_sites), GenusSpecies)

  # Presence/absence matrix with SiteName, Database, Latitude, Longitude first.
  matrix_wide <- bin_data %>%
    distinct(SiteName, GenusSpecies) %>%
    mutate(Presence = 1L) %>%
    pivot_wider(names_from = GenusSpecies, values_from = Presence, values_fill = 0L) %>%
    right_join(select(site_info, SiteName, Database, Latitude, Longitude), by = "SiteName") %>%
    select(SiteName, Database, Latitude, Longitude, all_of(species_info$GenusSpecies)) %>%
    arrange(match(SiteName, site_info$SiteName))
  save_both(matrix_wide, paste0("matrix_", stage))

  # Long form for the heatmap.
  heat <- matrix_wide %>%
    select(-Database, -Latitude, -Longitude) %>%
    pivot_longer(-SiteName, names_to = "GenusSpecies", values_to = "Presence") %>%
    left_join(select(species_info, GenusSpecies, Order), by = "GenusSpecies") %>%
    mutate(
      SiteName = factor(SiteName, levels = site_info$SiteName),
      # Reverse so the first order sits at the top of the plot.
      GenusSpecies = factor(GenusSpecies, levels = rev(species_info$GenusSpecies)),
      Presence = factor(Presence, levels = c(0, 1))
    )

  # Order blocks: separator lines and labels on the right.
  n_sp <- nrow(species_info)
  n_si <- nrow(site_info)
  blocks <- species_info %>%
    mutate(y = n_sp - row_number() + 1) %>%
    group_by(Order) %>%
    summarise(y_top = max(y), y_bottom = min(y), .groups = "drop") %>%
    mutate(y_mid = (y_top + y_bottom) / 2)
  separators <- blocks$y_bottom[blocks$y_bottom > 1] - 0.5
  right_margin <- 8 * max(nchar(as.character(blocks$Order))) + 30   # room for labels (pt)

  p <- ggplot(heat, aes(x = SiteName, y = GenusSpecies, fill = Presence)) +
    geom_tile(colour = "grey75", linewidth = 0.2) +
    geom_hline(yintercept = separators, colour = "black", linewidth = 0.6) +
    geom_text(data = blocks, aes(x = n_si + 0.5, y = y_mid, label = paste0("  ", Order)),
              inherit.aes = FALSE, hjust = 0, size = 3.6, fontface = "bold") +
    scale_fill_manual(values = c("0" = absent_colour, "1" = present_colour),
                      labels = c("0" = "Absent", "1" = "Present"),
                      name = NULL, drop = FALSE) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.5, n_si + 0.5), clip = "off") +   # labels sit outside
    labs(title = paste("Site by species matrix -", stage_labels[s]),
         subtitle = sprintf("%d sites x %d species (FAUNMAP + PBDB); no minimum-species filter",
                            n_si, n_sp),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(colour = "grey30", hjust = 0.5),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
      axis.text.y = element_text(face = "italic", size = 7),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      legend.position = "bottom",
      plot.margin = margin(10, right_margin, 10, 10)
    )

  width  <- min(max(12, n_si * 0.22 + 4), 50)
  height <- min(max(8, n_sp * 0.16 + 3), 50)
  ggsave(file.path(output_dir, paste0("heatmap_", stage, ".pdf")), p,
         width = width, height = height, limitsize = FALSE, bg = "white")
  ggsave(file.path(output_dir, paste0("heatmap_", stage, ".png")), p,
         width = width, height = height, dpi = png_dpi, limitsize = FALSE, bg = "white")

  cat(sprintf("%4d sites x %4d species -> matrix_%s.xlsx/.csv, heatmap_%s.pdf/.png\n",
              n_si, n_sp, stage, stage))
}

# =============================================================================
# 7. SITE ACCOUNTING: WHERE DID EACH STEP 3 SITE GO?
#
# Starts from every named site in each stage in Step 3 (the sites on the
# Step 4 maps) and records the first step at which a site dropped out.
# No minimum-species filter is involved; sites are lost only when they have
# no usable records left.
# =============================================================================

cat("\n=== 7. SITE ACCOUNTING ===\n")

index_sites <- site_index_raw %>%
  transmute(Stage_Number = as.integer(Stage_Number), Database,
            SiteName = clean_text(SiteName)) %>%
  filter(!is.na(SiteName), Stage_Number %in% 1:5) %>%
  distinct()

in_set <- function(df, set) {
  paste(df$Stage_Number, df$Database, df$SiteName) %in%
    paste(set$Stage_Number, set$Database, set$SiteName)
}

reasons_order <- c(
  "Retained in matrix",
  "No fauna records linked to this site",
  "Its fauna records have no genus",
  "No records in the selected orders",
  "Only sp./cf. records that could not be resolved"
)

site_attrition <- index_sites %>%
  mutate(
    Status = case_when(
      in_set(., site_set(all_final)) ~ reasons_order[1],
      !in_set(., sites_linked)       ~ reasons_order[2],
      !in_set(., sites_genus)        ~ reasons_order[3],
      !in_set(., sites_orders)       ~ reasons_order[4],
      TRUE                           ~ reasons_order[5]
    ),
    Stage = stage_names[Stage_Number]
  ) %>%
  arrange(Stage_Number, Database, Status, SiteName) %>%
  select(Stage_Number, Stage, Database, SiteName, Status)

attrition_summary <- site_attrition %>%
  count(Stage_Number, Stage, Database, Status) %>%
  mutate(Status = factor(Status, levels = reasons_order)) %>%
  tidyr::pivot_wider(names_from = Status, values_from = n, values_fill = 0,
                     names_sort = TRUE) %>%
  mutate(Sites_in_Step3 = rowSums(across(any_of(reasons_order)))) %>%
  relocate(Sites_in_Step3, .after = Database) %>%
  arrange(Stage_Number, Database)

# Sites in the matrices that are not in Step 3's index (should be none).
extra_sites <- site_set(all_final) %>%
  filter(!in_set(., index_sites))

print(as.data.frame(attrition_summary), row.names = FALSE)
if (nrow(extra_sites) > 0) {
  cat(sprintf("  NOTE: %d site-stage combinations are in the matrices but not in Step 3's index.\n",
              nrow(extra_sites)))
}
write.csv(site_attrition, file.path(output_dir, "site_attrition.csv"), row.names = FALSE, na = "")
cat("  Every site and its status: site_attrition.csv\n")

# =============================================================================
# 8. SUMMARIES
# =============================================================================

cat("\n=== 8. SUMMARIES ===\n")

stage_summary <- all_final %>%
  group_by(Stage_Number, Stage) %>%
  summarise(n_sites = n_distinct(SiteName),
            n_FAUNMAP_sites = n_distinct(SiteName[Database == "FAUNMAP"]),
            n_PBDB_sites = n_distinct(SiteName[Database == "PBDB"]),
            n_species = n_distinct(GenusSpecies),
            n_genera = n_distinct(Genus),
            n_records = n(), .groups = "drop")

order_summary <- all_final %>%
  group_by(Stage_Number, Stage, Order) %>%
  summarise(n_species = n_distinct(GenusSpecies), n_records = n(), .groups = "drop")

site_richness <- all_final %>%
  group_by(Stage_Number, Stage, SiteName, Database) %>%
  summarise(n_species = n_distinct(GenusSpecies), .groups = "drop") %>%
  arrange(Stage_Number, desc(n_species))

write_xlsx(list(StageSummary = stage_summary,
                OrderSummary = order_summary,
                Resolution = resolution_counts,
                SiteAccounting = attrition_summary,
                SiteRichness = site_richness),
           file.path(output_dir, "summary.xlsx"))

print(as.data.frame(stage_summary), row.names = FALSE)

cat("\n=== STEP 6 COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: all_final, spcf_log, stage_summary, order_summary\n")
