# =============================================================================
# STEP 9e: WHICH TAXA DROVE THE CHANGE? WINNERS, LOSERS, IMMIGRANTS, BODY SIZE
#
#   a) Occupancy (winners and losers): the share of a province's sites in
#      which each taxon occurs, per stage. Homogenization = fewer taxa known
#      from a single site and more taxa found at many sites; the taxa that
#      gain most are the "winners" that spread across the landscape.
#   b) Immigrants: taxa that entered North America during the Plio-
#      Pleistocene (from Eurasia via Beringia, or from South America in the
#      Great American Biotic Interchange), set by genus in 'immigrant_genera'
#      below - CHECK AND EDIT THIS LIST against the literature for your taxa.
#      If immigrants make up a larger share of widespread taxa than of all
#      taxa, the spread of immigrants is driving homogenization.
#   c) Body size (Rowan et al. 2024, Figs. 3b and 4b): Step 7's within-
#      province beta_SIM partitioned by body-size class (Eq. 4), and its
#      endemic and shared parts (Eqs. 5-6). Body masses come from the Smith et
#      al. table (aao5987, ln mass in grams), looked up only for the species
#      in the regional pools. Same site subsets as Step 7, so the grey line
#      equals Step 7.
#
# Inputs:  all_records_final.csv (Step 6), site_physio_database.csv (Step 5),
#          aao5987-smith-sm-137-190.xlsx (Smith et al.)
# Outputs (Outputs/9e_winners_immigrants_size/):
#   occupancy_classes_all_provinces.png/.pdf     a) taxa by occupancy class
#   occupancy_heatmap_<province>.png/.pdf        a) most widespread taxa
#   immigrants_all_provinces.png/.pdf            b) immigrant share
#   size_partition_<province>.png/.pdf           c) Fig. 3b style
#   size_endemic_shared_<province>.png/.pdf      c) Fig. 4b style
#   occupancy.csv, winners_losers.csv, immigrants.csv, species_size.csv,
#   size_partition.csv
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 7.
# Needs: dplyr, tidyr, ggplot2, readxl.
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "readxl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

records_file <- file.path("Outputs", "6_matrices", "all_records_final.csv")
physio_file  <- file.path("Outputs", "maps", "physiographic", "site_physio_database.csv")
output_dir   <- file.path(work_dir, "Outputs", "9e_winners_immigrants_size")

provinces <- c("Pacific Border", "Basin and Range", "Great Plains", "Coastal Plain")

# One colour per province (validated colour-blind-safe set, as in Step 7).
province_colours <- c(
  "Pacific Border"  = "#2a78d6",   # blue
  "Basin and Range" = "#eb6834",   # orange
  "Great Plains"    = "#1baf7a",   # aqua
  "Coastal Plain"   = "#4a3aa7"    # violet
)

stage_names <- c("Zanclean", "Piacenzian", "Gelasian", "Calabrian", "Chibanian")
stage_older <- c(4.700, 3.600, 2.580, 1.800, 0.7741)
stage_young <- c(3.600, 2.580, 1.800, 0.7741, 0.129)
stage_mid   <- (stage_older + stage_young) / 2

# a) Occupancy
occupancy_level  <- "species"   # "species" or "genus" (genus avoids chronospecies:
                                # e.g. Equus simplicidens -> Equus scotti counts as one taxon)
min_sites_stage  <- 3           # province x stage needs this many sites
widespread_share <- 0.25        # "widespread" = found at >= 25% of the province's sites
top_n_taxa       <- 20          # taxa shown in each occupancy heatmap

# b) Immigrant genera - a STARTING LIST. Check arrival times and add/remove
#    genera for your data before using these results.
immigrant_genera <- c(
  # From Eurasia (via Beringia)
  "Mammuthus" = "Eurasian", "Bison" = "Eurasian", "Ovibos" = "Eurasian",
  "Rangifer" = "Eurasian", "Alces" = "Eurasian", "Cervus" = "Eurasian",
  "Saiga" = "Eurasian", "Panthera" = "Eurasian", "Gulo" = "Eurasian",
  "Lemmus" = "Eurasian", "Dicrostonyx" = "Eurasian", "Microtus" = "Eurasian",
  # From South America (Great American Biotic Interchange)
  "Glossotherium" = "South American", "Paramylodon" = "South American",
  "Eremotherium" = "South American", "Nothrotheriops" = "South American",
  "Glyptotherium" = "South American", "Holmesina" = "South American",
  "Pampatherium" = "South American", "Dasypus" = "South American",
  "Erethizon" = "South American", "Coendou" = "South American",
  "Hydrochoerus" = "South American", "Neochoerus" = "South American",
  "Didelphis" = "South American"
)

# c) Body size
diet_file     <- file.path(work_dir, "aao5987-smith-sm-137-190.xlsx")   # Smith et al. table
downloads_dir <- "C:/Users/shrut/Downloads"
mass_genus_fallback <- FALSE    # TRUE: species missing from the table take the mean
                                # mass of their congeners in the table
size_breaks_kg <- c(1, 18, 80, 350, 1000)   # Rowan et al. classes, plus <1 kg for small mammals
size_labels    <- c("<1 kg", "1-18 kg", "18-80 kg", "80-350 kg", "350-1,000 kg", ">1,000 kg")
size_colours   <- setNames(c("#cfd8e3", "#a9b9cc", "#8399b4", "#5d7a9c", "#3a5a80", "#1f3b5e"),
                           size_labels)
# Keep these three identical to Step 7 so the overall values match Step 7.
n_units_sample <- 3
n_resamples    <- 999
step7_file     <- file.path("Outputs", "7_regional_pools", "beta_spatial.csv")

# Figure look (slide-ready).
base_size <- 18
slide_w   <- 13.33   # inches, 16:9
slide_h   <- 7.5
fig_dpi   <- 300
ink       <- "#1f1f1f"
ink_soft  <- "#5a5a5a"
grid_col  <- "#e6e6e3"
band_cols <- c("#f3f3f0", "#ffffff")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================

read_input <- function(file) {
  path <- file.path(work_dir, file)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun the earlier steps first.")
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                 colClasses = "character", na.strings = c("", "NA"))
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))
  cat(sprintf("  %-55s %7d rows\n", file, nrow(df)))
  df
}

# Save without stopping when a file is open in Excel / locked by OneDrive.
save_with_fallback <- function(writer, path) {
  ok <- tryCatch({ writer(path); TRUE }, error = function(e) FALSE,
                 warning = function(w) FALSE)
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
save_fig <- function(p, name, width = slide_w, height = slide_h) {
  pdf_device <- if (capabilities("cairo")) cairo_pdf else pdf
  save_with_fallback(function(f) ggsave(f, p, width = width, height = height,
                                         device = pdf_device, bg = "white"),
                     file.path(output_dir, paste0(name, ".pdf")))
  save_with_fallback(function(f) ggsave(f, p, width = width, height = height,
                                         dpi = fig_dpi, bg = "white"),
                     file.path(output_dir, paste0(name, ".png")))
}
file_stub <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# Multisite Simpson (turnover), nestedness and Sorensen dissimilarity, with the
# endemic and shared components (Baselga 2010; Rowan et al. 2024 Eqs. 1-3).
beta_multi_all <- function(pools) {
  pools <- pools[lengths(pools) > 0]
  k <- length(pools)
  if (k < 2) return(c(n = k, beta_SIM = NA, beta_SNE = NA, beta_SOR = NA,
                      beta_SIM_END = NA, beta_SIM_SH = NA))
  S_T <- length(unique(unlist(pools))); sum_S <- sum(lengths(pools))
  mins <- 0; maxs <- 0; a <- 0
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    b1 <- length(setdiff(pools[[i]], pools[[j]])); b2 <- length(setdiff(pools[[j]], pools[[i]]))
    mins <- mins + min(b1, b2); maxs <- maxs + max(b1, b2)
    a <- a + length(intersect(pools[[i]], pools[[j]]))
  }
  shared <- sum_S - S_T
  sim <- mins / (shared + mins)
  sor <- (mins + maxs) / (2 * shared + mins + maxs)
  c(n = k, beta_SIM = sim, beta_SNE = sor - sim, beta_SOR = sor,
    beta_SIM_END = mins / (mins + a), beta_SIM_SH = shared / sum_S)
}

# Pairwise Simpson dissimilarity.
beta_sim_pair <- function(A, B) {
  if (length(A) == 0 || length(B) == 0) return(NA_real_)
  b <- length(setdiff(A, B)); c <- length(setdiff(B, A))
  min(b, c) / (length(intersect(A, B)) + min(b, c))
}

# Shared slide theme.
theme_slide <- function() {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = ink),
      plot.title = element_text(face = "bold", size = base_size * 1.25, margin = margin(b = 4)),
      plot.subtitle = element_text(colour = ink_soft, size = base_size * 0.8, margin = margin(b = 12)),
      plot.caption = element_text(colour = ink_soft, size = base_size * 0.6, hjust = 0, lineheight = 1.1),
      plot.title.position = "plot", plot.caption.position = "plot",
      axis.title = element_text(colour = ink_soft, size = base_size * 0.85),
      axis.text = element_text(colour = ink, size = base_size * 0.72),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size * 0.9, hjust = 0),
      plot.margin = margin(18, 24, 14, 18),
      legend.position = "none"
    )
}

beta_lab <- quote(beta[SIM])

# Alternating background bands, one per stage (x axis in Ma, oldest left).
bands <- data.frame(xmin = stage_young, xmax = stage_older, Stage = stage_names,
                    lab = substr(stage_names, 1, 4),
                    fill = rep(band_cols, length.out = 5))
stage_bands <- function(y_lab = NULL, lab_size = base_size / 4.4) {
  out <- list(geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                                          fill = fill), inherit.aes = FALSE),
              scale_fill_identity())
  if (!is.null(y_lab)) {
    out <- c(out, list(geom_text(data = bands, aes(x = (xmin + xmax) / 2, y = y_lab, label = lab),
                                 size = lab_size, colour = ink_soft, inherit.aes = FALSE)))
  }
  out
}
age_axis <- function() {
  scale_x_reverse(limits = c(4.7, 0.129), breaks = c(4, 3, 2, 1), expand = c(0.01, 0))
}

# Line segments only between neighbouring stages that both have a value, so a
# missing stage leaves a gap.
gap_segments <- function(d, group, x = "Mid_Ma", y = "value") {
  d %>%
    arrange(.data[[group]], Stage_Number) %>%
    group_by(.data[[group]]) %>%
    mutate(x2 = lead(.data[[x]]), y2 = lead(.data[[y]]),
           adjacent = lead(Stage_Number) - Stage_Number == 1) %>%
    ungroup() %>%
    filter(!is.na(.data[[y]]), !is.na(y2), adjacent)
}

# =============================================================================
# 1. READ RECORDS AND LINK SITES TO PROVINCES (same rules as Steps 7 and 8)
# =============================================================================

cat("=== 1. READING RECORDS ===\n")
records <- read_input(records_file)
physio  <- read_input(physio_file)

site_province_all <- physio %>%
  transmute(Stage_Number = as.integer(Stage_Number), Database,
            SiteName = trimws(SiteName), Province = trimws(US_Province)) %>%
  filter(!is.na(SiteName), !is.na(Province))

site_province <- site_province_all %>%
  count(Stage_Number, Database, SiteName, Province) %>%
  group_by(Stage_Number, Database, SiteName) %>%
  arrange(desc(n), Province, .by_group = TRUE) %>%
  ungroup() %>%
  distinct(Stage_Number, Database, SiteName, .keep_all = TRUE) %>%
  select(Stage_Number, Database, SiteName, Province)

num_col <- function(df, col) {
  if (col %in% names(df)) suppressWarnings(as.numeric(df[[col]])) else rep(NA_real_, nrow(df))
}

recs <- records %>%
  mutate(Latitude = num_col(records, "Latitude"), Longitude = num_col(records, "Longitude")) %>%
  transmute(Stage_Number = as.integer(Stage_Number), Database, SiteName = trimws(SiteName),
            GenusSpecies, Genus = sub(" .*$", "", trimws(GenusSpecies)), Order,
            Latitude, Longitude) %>%
  filter(!is.na(GenusSpecies)) %>%
  inner_join(site_province, by = c("Stage_Number", "Database", "SiteName")) %>%
  mutate(Province = provinces[match(tolower(Province), tolower(provinces))]) %>%
  filter(!is.na(Province), Stage_Number %in% 1:5) %>%
  mutate(Site = paste(Database, SiteName, sep = " | "),
         Stage = factor(stage_names[Stage_Number], levels = stage_names))

cat(sprintf("  Species records in the four provinces: %d (%d sites x stages, %d species)\n",
            nrow(recs), n_distinct(paste(recs$Site, recs$Stage_Number)), n_distinct(recs$GenusSpecies)))

# Species list of each site (one province and stage).
site_lists <- function(d, taxon = "GenusSpecies") lapply(split(d[[taxon]], d$Site), unique)

recs <- mutate(recs, Taxon = if (occupancy_level == "genus") Genus else GenusSpecies)
recs$Immigrant <- unname(immigrant_genera[recs$Genus])
recs$Immigrant[is.na(recs$Immigrant)] <- "Native"

# =============================================================================
# 2. OCCUPANCY: WINNERS AND LOSERS (a)
# =============================================================================

cat(sprintf("\n=== 2. OCCUPANCY (%s level) ===\n", occupancy_level))

site_n <- recs %>% distinct(Province, Stage_Number, Site) %>% count(Province, Stage_Number, name = "n_sites")
occ <- recs %>%
  distinct(Province, Stage_Number, Site, Taxon, Genus, Immigrant) %>%
  group_by(Province, Stage_Number, Taxon) %>%
  summarise(Genus = first(Genus), Immigrant = first(Immigrant), n_sites_with = n(), .groups = "drop") %>%
  left_join(site_n, by = c("Province", "Stage_Number")) %>%
  filter(n_sites >= min_sites_stage) %>%
  mutate(occupancy = n_sites_with / n_sites,
         Stage = factor(stage_names[Stage_Number], levels = stage_names),
         Province = factor(Province, levels = provinces),
         Class = factor(case_when(
           n_sites_with == 1               ~ "One site only",
           occupancy < widespread_share    ~ "2+ sites, <25%",
           occupancy < 0.5                 ~ "25-50% of sites",
           TRUE                            ~ "50% of sites or more"),
           levels = c("One site only", "2+ sites, <25%", "25-50% of sites", "50% of sites or more")))
levels(occ$Class)[2] <- sprintf("2+ sites, <%.0f%%", 100 * widespread_share)
levels(occ$Class)[3] <- sprintf("%.0f-50%% of sites", 100 * widespread_share)
save_csv(mutate(occ, occupancy = round(occupancy, 3)), "occupancy")

class_share <- occ %>% count(Province, Stage, Stage_Number, Class) %>%
  group_by(Province, Stage) %>% mutate(share = n / sum(n)) %>% ungroup()
cat("  Share of taxa known from one site only (falling = homogenization):\n")
print(as.data.frame(class_share %>% filter(Class == "One site only") %>%
        transmute(Province, Stage, share = round(share, 3))), row.names = FALSE)

# Winners and losers between successive stages (absent = 0% occupancy).
wl <- bind_rows(lapply(provinces, function(p) bind_rows(lapply(1:4, function(s) {
  a <- filter(occ, Province == p, Stage_Number == s); b <- filter(occ, Province == p, Stage_Number == s + 1)
  if (nrow(a) == 0 || nrow(b) == 0) return(NULL)
  full_join(select(a, Taxon, occ_from = occupancy), select(b, Taxon, occ_to = occupancy), by = "Taxon") %>%
    mutate(across(c(occ_from, occ_to), ~ coalesce(.x, 0)), change = occ_to - occ_from,
           Province = p, Transition = paste(stage_names[s], "to", stage_names[s + 1])) %>%
    arrange(desc(change)) %>%
    { bind_rows(mutate(head(., 10), Type = "Winner"), mutate(head(arrange(., change), 10), Type = "Loser")) }
}))))
if (nrow(wl)) {
  wl <- wl %>% mutate(Immigrant = coalesce(unname(immigrant_genera[sub(" .*$", "", Taxon)]), "Native"),
                      across(c(occ_from, occ_to, change), ~ round(.x, 3))) %>%
    select(Province, Transition, Type, Taxon, Immigrant, occ_from, occ_to, change)
}
save_csv(wl, "winners_losers")

class_cols <- setNames(c("#cde2fb", "#8dbbf2", "#3987e5", "#0d366b"), levels(occ$Class))
p_class <- ggplot(class_share, aes(Stage, share, fill = Class)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.6) +
  facet_wrap(~ Province, ncol = 2, drop = FALSE) +
  scale_fill_manual(values = class_cols, name = NULL) +
  scale_x_discrete(labels = substr(stage_names, 1, 4), drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0)) +
  labs(title = "Are taxa becoming more widespread?",
       subtitle = sprintf("Share of %s by the proportion of the province's sites they occur at; shrinking light parts = homogenization",
                          ifelse(occupancy_level == "genus", "genera", "species")),
       x = NULL, y = "Share of taxa",
       caption = paste0("Stages with fewer than ", min_sites_stage, " sites in a province are left out.")) +
  theme_slide() +
  theme(legend.position = "top", legend.justification = "left",
        panel.grid.major.y = element_blank())
save_fig(p_class, "occupancy_classes_all_provinces", height = slide_h * 1.1)

for (p in provinces) {
  d <- filter(occ, Province == p)
  if (nrow(d) == 0) next
  top <- d %>% group_by(Taxon) %>% summarise(m = max(occupancy), first = min(Stage_Number), .groups = "drop") %>%
    arrange(desc(m)) %>% head(top_n_taxa) %>% arrange(first, desc(m))
  h <- d %>% filter(Taxon %in% top$Taxon) %>%
    select(Taxon, Immigrant, Stage, occupancy) %>%
    complete(Taxon, Stage, fill = list(occupancy = 0)) %>%
    group_by(Taxon) %>% mutate(Immigrant = first(na.omit(Immigrant))) %>% ungroup() %>%
    mutate(Label = ifelse(Immigrant != "Native", paste0(Taxon, " *"), Taxon))
  lab_levels <- rev(unique(h$Label[match(top$Taxon, h$Taxon)]))
  h <- mutate(h, Label = factor(Label, levels = lab_levels),
              txt = ifelse(occupancy > 0, sprintf("%.0f", 100 * occupancy), ""),
              txt_col = ifelse(occupancy > 0.5, "white", ink))
  present <- unique(as.character(d$Stage))
  h <- filter(h, Stage %in% present)
  pp <- ggplot(h, aes(Stage, Label, fill = occupancy)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = txt, colour = txt_col), size = base_size / 4.4, fontface = "bold") +
    scale_colour_identity() +
    scale_fill_gradientn(colours = c("#f4f4f2", "#cde2fb", "#6da7ec", "#256abf", "#0d366b"),
                         limits = c(0, 1), labels = function(x) paste0(100 * x, "%"),
                         name = "Sites\noccupied") +
    labs(title = paste(p, "- the most widespread taxa"),
         subtitle = sprintf("Percentage of the province's sites where each %s occurs, stage by stage",
                            ifelse(occupancy_level == "genus", "genus", "species")),
         x = NULL, y = NULL,
         caption = "* = immigrant genus (see immigrant_genera in the settings). Taxa ordered by first appearance.") +
    theme_slide() +
    theme(panel.grid = element_blank(), legend.position = "right",
          axis.text.y = element_text(size = base_size * 0.62, face = "italic"),
          legend.key.height = unit(1.2, "cm"))
  save_fig(pp, paste0("occupancy_heatmap_", file_stub(p)), height = slide_h * 1.1)
}
cat("  Occupancy figures saved\n")

# =============================================================================
# 3. IMMIGRANTS (b)
# =============================================================================

cat("\n=== 3. IMMIGRANTS ===\n")

imm <- bind_rows(
  occ %>% mutate(Group = "All taxa"),
  occ %>% filter(occupancy >= widespread_share) %>% mutate(Group = "Widespread taxa"),
  occ %>% filter(n_sites_with == 1) %>% mutate(Group = "One-site taxa")
) %>%
  group_by(Province, Stage_Number, Stage, Group) %>%
  summarise(n_taxa = n(), n_immigrant = sum(Immigrant != "Native"),
            n_Eurasian = sum(Immigrant == "Eurasian"), n_South_American = sum(Immigrant == "South American"),
            .groups = "drop") %>%
  mutate(share_immigrant = n_immigrant / n_taxa, Mid_Ma = stage_mid[Stage_Number],
         Group = factor(Group, levels = c("All taxa", "Widespread taxa", "One-site taxa")))
save_csv(mutate(imm, share_immigrant = round(share_immigrant, 3)), "immigrants")
print(as.data.frame(imm %>% filter(Group != "One-site taxa") %>%
        transmute(Province, Stage, Group, n_taxa, immigrant_share = round(share_immigrant, 3))),
      row.names = FALSE)

group_cols <- c("All taxa" = ink_soft, "Widespread taxa" = "#eb6834", "One-site taxa" = "#2a78d6")
imm_plot <- mutate(imm, value = share_immigrant, Series = interaction(Province, Group))
y_max <- max(0.1, imm_plot$value, na.rm = TRUE) * 1.2
p_imm <- ggplot(imm_plot) +
  stage_bands() +
  geom_segment(data = gap_segments(imm_plot, "Series"),
               aes(x = Mid_Ma, xend = x2, y = value, yend = y2, colour = Group), linewidth = 1.2) +
  geom_point(aes(Mid_Ma, value, colour = Group), shape = 21, fill = "white", size = 3.4, stroke = 1.6) +
  facet_wrap(~ Province, ncol = 2, drop = FALSE) +
  scale_colour_manual(values = group_cols, name = NULL) +
  age_axis() +
  scale_y_continuous(limits = c(0, y_max), labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Are immigrants the taxa that spread?",
       subtitle = "Share of taxa that are Plio-Pleistocene immigrants; orange above grey = immigrants over-represented among widespread taxa",
       x = "Age (Ma)", y = "Immigrant share",
       caption = sprintf(paste0("Widespread = found at %.0f%% or more of the province's sites. Immigrants set by genus ",
                                "(Eurasian or South American; edit immigrant_genera in the settings)."),
                         100 * widespread_share)) +
  theme_slide() +
  theme(legend.position = "top", legend.justification = "left",
        panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5))
save_fig(p_imm, "immigrants_all_provinces", height = slide_h * 1.1)
cat("  Immigrant figure saved\n")

# =============================================================================
# 4. BODY SIZE PARTITION (c)
# =============================================================================

cat("\n=== 4. BODY SIZE (Smith et al. masses) ===\n")

clean_name <- function(x) {
  x <- tolower(x)
  x <- gsub("\\b(cf|aff|nr|sp|spp|indet|n|gen|ex|gr)\\b\\.?", " ", x)
  x <- gsub("[^a-z ]", " ", x)
  x <- gsub("\\s+", " ", trimws(x))
  vapply(strsplit(x, " ", fixed = TRUE), function(v) paste(head(v, 2), collapse = " "), character(1))
}

if (!file.exists(diet_file)) {
  hit <- list.files(downloads_dir, pattern = "aao5987.*\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  if (length(hit) == 0) stop("Smith et al. table not found:\n  ", diet_file)
  diet_file <- hit[1]
}
raw <- as.data.frame(suppressMessages(readxl::read_excel(diet_file, col_names = FALSE, col_types = "text")))
header_row <- which(apply(raw, 1, function(r) any(trimws(r) == "Genus_species", na.rm = TRUE)))[1]
species_col <- which(trimws(unlist(raw[header_row, ])) == "Genus_species")[1]
hdr_rows <- max(1, header_row - 2):min(nrow(raw), header_row + 2)
header_text <- vapply(seq_along(raw), function(j) paste(na.omit(raw[hdr_rows, j]), collapse = " "), character(1))
mass_col <- grep("ln Mass", header_text, ignore.case = TRUE)[1]
if (is.na(mass_col)) stop("No 'ln Mass' column found in ", basename(diet_file))
smith <- raw[(header_row + 1):nrow(raw), c(species_col, mass_col)]
names(smith) <- c("Genus_species", "mass_cell")
smith <- smith %>%
  filter(!is.na(Genus_species)) %>%
  mutate(Name_key = clean_name(Genus_species),
         ln_mass_g = suppressWarnings(as.numeric(sub("^\\s*(-?[0-9]*\\.?[0-9]+).*$", "\\1", mass_cell))),
         mass_kg = exp(ln_mass_g) / 1000) %>%
  filter(!is.na(mass_kg))

pool_sp <- recs %>% distinct(GenusSpecies) %>% mutate(Name_key = clean_name(GenusSpecies),
                                                      Genus_key = sub(" .*$", "", Name_key))
by_sp <- smith %>% group_by(Name_key) %>% summarise(mass_kg_species = first(mass_kg), .groups = "drop")
by_gen <- smith %>% mutate(Genus_key = sub(" .*$", "", Name_key)) %>%
  filter(Genus_key %in% pool_sp$Genus_key) %>%
  group_by(Genus_key) %>% summarise(mass_kg_genus = exp(mean(log(mass_kg))), .groups = "drop")
species_size <- pool_sp %>%
  left_join(by_sp, by = "Name_key") %>% left_join(by_gen, by = "Genus_key") %>%
  mutate(mass_kg = coalesce(mass_kg_species, if (mass_genus_fallback) mass_kg_genus else NA_real_),
         Mass_source = case_when(!is.na(mass_kg_species) ~ "species (Smith et al.)",
                                 mass_genus_fallback & !is.na(mass_kg_genus) ~ "genus mean (Smith et al.)",
                                 TRUE ~ "not found"),
         Size_class = as.character(cut(mass_kg, c(0, size_breaks_kg, Inf), labels = size_labels, right = FALSE)),
         Size_class = coalesce(Size_class, "Unknown size")) %>%
  select(GenusSpecies, mass_kg, Mass_source, Size_class)
save_csv(mutate(species_size, mass_kg = signif(mass_kg, 4)), "species_size")
print(as.data.frame(count(species_size, Size_class)), row.names = FALSE)

size_cats <- c(size_labels, "Unknown size")
category_of <- setNames(species_size$Size_class, species_size$GenusSpecies)

# Additive partitions (Rowan et al. 2024 Eqs. 4-6): beta_SIM_f, END_f, SH_f.
beta_sim_partition <- function(pools, category_of, categories) {
  k <- length(pools)
  count <- function(sp) tabulate(match(category_of[sp], categories), nbins = length(categories))
  S_T <- length(unique(unlist(pools))); sum_S <- sum(lengths(pools))
  num <- numeric(length(categories)); min_sum <- 0; a_sum <- 0
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    b_ij <- setdiff(pools[[i]], pools[[j]]); b_ji <- setdiff(pools[[j]], pools[[i]])
    num <- num + if (length(b_ij) < length(b_ji)) count(b_ij) else
      if (length(b_ij) > length(b_ji)) count(b_ji) else (count(b_ij) + count(b_ji)) / 2
    min_sum <- min_sum + min(length(b_ij), length(b_ji))
    a_sum <- a_sum + length(intersect(pools[[i]], pools[[j]]))
  }
  S_if <- Reduce(`+`, lapply(pools, count)); S_Tf <- count(unique(unlist(pools)))
  c(num / (sum_S - S_T + min_sum), num / (min_sum + a_sum), (S_if - S_Tf) / sum_S)
}
unit_subsets <- function(k, n) {
  if (k < n || n < 2) return(list())
  if (choose(k, n) <= n_resamples) combn(k, n, simplify = FALSE) else
    replicate(n_resamples, sample.int(k, n), simplify = FALSE)
}

set.seed(2024)   # same seed and loop order as Step 7 -> same site subsets
nc <- length(size_cats)
size_part <- bind_rows(lapply(provinces, function(p) bind_rows(lapply(1:5, function(s) {
  d <- filter(recs, Province == p, Stage_Number == s)
  units <- lapply(split(d$GenusSpecies, d$Site), unique)
  subsets <- unit_subsets(length(units), n_units_sample)
  n_by <- tabulate(match(category_of[unique(d$GenusSpecies)], size_cats), nbins = nc)
  if (length(subsets) == 0) {
    return(data.frame(Province = p, Stage_Number = s, Size = size_cats, beta_SIM_f = NA, beta_END_f = NA,
                      beta_SH_f = NA, n_species = n_by, n_sites = length(units), beta_SIM = NA,
                      beta_SIM_lo95 = NA, beta_SIM_hi95 = NA))
  }
  parts <- sapply(subsets, function(ix) beta_sim_partition(units[ix], category_of, size_cats))
  tot <- colSums(parts[1:nc, , drop = FALSE])
  data.frame(Province = p, Stage_Number = s, Size = size_cats,
             beta_SIM_f = rowMeans(parts[1:nc, , drop = FALSE]),
             beta_END_f = rowMeans(parts[nc + 1:nc, , drop = FALSE]),
             beta_SH_f = rowMeans(parts[2 * nc + 1:nc, , drop = FALSE]),
             n_species = n_by, n_sites = length(units), beta_SIM = mean(tot),
             beta_SIM_lo95 = unname(quantile(tot, 0.025)), beta_SIM_hi95 = unname(quantile(tot, 0.975)))
})))) %>%
  mutate(Stage = factor(stage_names[Stage_Number], levels = stage_names), Mid_Ma = stage_mid[Stage_Number],
         Province = factor(Province, levels = provinces), Size = factor(Size, levels = size_cats))
save_csv(mutate(size_part, across(where(is.numeric), ~ round(.x, 4))), "size_partition")

s7_path <- file.path(work_dir, step7_file)
if (file.exists(s7_path)) {
  s7 <- read.csv(s7_path, stringsAsFactors = FALSE, check.names = FALSE)
  chk <- size_part %>% distinct(Province, Stage, beta_SIM) %>%
    mutate(Province = as.character(Province), Stage = as.character(Stage)) %>%
    inner_join(transmute(s7, Province, Stage, b7 = beta_SIM), by = c("Province", "Stage"))
  dif <- suppressWarnings(max(abs(chk$beta_SIM - chk$b7), na.rm = TRUE))
  cat(sprintf("  Check against Step 7: largest difference in overall beta_SIM = %s\n",
              if (is.finite(dif)) format(round(dif, 4), nsmall = 4) else "n/a"))
}

size_present <- size_part %>% group_by(Province, Size) %>%
  summarise(any = sum(n_species) > 0, .groups = "drop")

theme_panels <- function() {
  theme_slide() +
    theme(panel.grid.major.y = element_blank(),
          panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5),
          strip.background = element_rect(fill = "#ececea", colour = "#c9c9c4", linewidth = 0.5),
          strip.text = element_text(face = "bold", size = base_size * 0.8, hjust = 0,
                                    margin = margin(t = 8, b = 8, l = 6)))
}

for (p in provinces) {
  keep <- size_present %>% filter(Province == p, any, Size != "Unknown size") %>% pull(Size) %>% as.character()
  if (length(keep) == 0) next
  d <- size_part %>% filter(Province == p, Size %in% keep) %>%
    mutate(Size = factor(as.character(Size), levels = keep), bar_fill = size_colours[as.character(Size)])
  overall <- distinct(d, Size, Stage_Number, Mid_Ma, beta_SIM, beta_SIM_lo95, beta_SIM_hi95) %>%
    mutate(value = beta_SIM)
  unk <- size_part %>% filter(Province == p) %>% group_by(Stage_Number) %>%
    summarise(u = sum(beta_SIM_f[Size == "Unknown size"]), t = first(beta_SIM), .groups = "drop") %>%
    filter(!is.na(t), t > 0)
  unk_txt <- if (nrow(unk)) sprintf("%.0f%%", 100 * sum(unk$u) / sum(unk$t)) else "n/a"
  ncol_f <- if (length(keep) <= 3) length(keep) else 3

  p3 <- ggplot(d) +
    stage_bands(y_lab = 1.08, lab_size = base_size / 5) +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), colour = grid_col, linewidth = 0.4) +
    geom_ribbon(data = filter(overall, !is.na(beta_SIM_lo95)),
                aes(x = Mid_Ma, ymin = beta_SIM_lo95, ymax = beta_SIM_hi95, group = Size),
                fill = "#d6d6d2", alpha = 0.55) +
    geom_rect(aes(xmin = stage_young[Stage_Number] + 0.05, xmax = stage_older[Stage_Number] - 0.05,
                  ymin = 0, ymax = beta_SIM_f, fill = bar_fill), colour = ink, linewidth = 0.3, na.rm = TRUE) +
    geom_segment(data = gap_segments(overall, "Size"),
                 aes(x = Mid_Ma, xend = x2, y = value, yend = y2), colour = "#6b6b68", linewidth = 1) +
    geom_point(data = overall, aes(Mid_Ma, beta_SIM), shape = 21, size = 3.2, fill = "#9a9a96",
               colour = "white", stroke = 1.1, na.rm = TRUE) +
    facet_wrap(~ Size, ncol = ncol_f) +
    age_axis() +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(title = paste0(p, ": which body sizes make sites differ?"),
         subtitle = "Bars = each size class's share of multisite Simpson dissimilarity among sites (shares add up to the grey line)",
         x = "Age (Ma)", y = beta_lab,
         caption = paste0("After Rowan et al. (2024, Fig. 3b; Eq. 4). Masses: Smith et al. (aao5987). Grey = all species ",
                          "(band = 95% range of subsets of ", n_units_sample, " sites).\nSpecies of unknown size supply ",
                          unk_txt, " of beta_SIM and are not drawn.")) +
    theme_panels()
  save_fig(p3, paste0("size_partition_", file_stub(p)),
           height = if (length(keep) <= 3) slide_h * 0.8 else slide_h * 1.1)

  comp <- c("Shared species", "Endemic species")
  e <- d %>% select(Size, Stage_Number, bar_fill, beta_SH_f, beta_END_f) %>%
    pivot_longer(c(beta_SH_f, beta_END_f), names_to = "Component", values_to = "value") %>%
    mutate(Component = factor(ifelse(Component == "beta_SH_f", comp[1], comp[2]), levels = comp))
  ymax <- max(0.1, e$value, na.rm = TRUE) * 1.32
  grid_y <- pretty(c(0, ymax / 1.32), n = 3); grid_y <- grid_y[grid_y > 0 & grid_y < ymax]
  row_lab <- data.frame(Component = factor(comp, levels = comp), Size = factor(keep[1], levels = keep))
  p4 <- ggplot(e) +
    stage_bands() +
    geom_hline(yintercept = grid_y, colour = grid_col, linewidth = 0.4) +
    geom_rect(aes(xmin = stage_young[Stage_Number] + 0.05, xmax = stage_older[Stage_Number] - 0.05,
                  ymin = 0, ymax = value, fill = bar_fill), colour = ink, linewidth = 0.3, na.rm = TRUE) +
    geom_text(data = row_lab, aes(x = 4.6, y = ymax * 0.97, label = Component), hjust = 0, vjust = 1,
              size = base_size / 3.8, colour = ink) +
    facet_grid(Component ~ Size) +
    age_axis() +
    scale_y_continuous(limits = c(0, ymax), breaks = c(0, grid_y), expand = c(0, 0)) +
    labs(title = paste0(p, ": endemic versus shared species by body size"),
         subtitle = "Proportion of shared (top) and endemic (bottom) species among sites\nFalling endemic bars + rising shared bars through time = homogenization",
         x = "Age (Ma)", y = "Proportion",
         caption = paste0("After Rowan et al. (2024, Fig. 4b; Eqs. 5-6), mean of subsets of ", n_units_sample,
                          " sites. Bands = stages, Zanclean (left) to Chibanian (right). Species of unknown size not drawn.")) +
    theme_panels() +
    theme(strip.text.y = element_blank(), strip.background.y = element_blank(),
          strip.text.x = element_text(face = "bold", size = base_size * 0.66, hjust = 0,
                                      margin = margin(t = 8, b = 8, l = 6)),
          axis.text = element_text(colour = ink, size = base_size * 0.6),
          panel.spacing.x = unit(0.6, "lines"))
  save_fig(p4, paste0("size_endemic_shared_", file_stub(p)))
  cat(sprintf("  %-16s %d size classes\n", p, length(keep)))
}

cat("\n=== STEP 9e COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: occ, wl, imm, species_size, size_part\n")
