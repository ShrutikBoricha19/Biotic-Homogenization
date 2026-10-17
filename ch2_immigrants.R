# =============================================================================
# CHAPTER 2: DID IMMIGRANTS BUILD PROVINCES OR ERASE THEM?
#
# Tests whether South American (GABI) and Eurasian (Beringian) immigrants
# raised or lowered multisite beta_SIM among the three provinces (Basin and
# Range, Coastal Plain, Great Plains incl. the Central Lowland west of the
# Mississippi), and through what range dynamics.
#
#   H1 (Darroch et al. 2014)  immigrants occupy fewer provinces / cells than
#                             natives, so removing them LOWERS beta_SIM
#   H2 (Stigall 2012)         immigrants are more widespread, so removing them
#                             RAISES beta_SIM
#   H3 (late homogenizer)     a late Beringian grazer (Bison) carries the
#                             shared-grazer rise in Bin 5
#   H4 (diachrony)            first appearances are staggered in the direction
#                             of entry: Eurasian taxa earlier in the west (Basin
#                             and Range, Great Plains), South American taxa
#                             earlier in the Coastal Plain
#
# Analyses (each run at species AND genus level, as the taxonomic check):
#   1. Origin effects on beta_SIM: all taxa vs without South American vs
#      without Eurasian vs natives only, with the equal-species subsampling of
#      Chapter 1 (Step 7 addendum). Each removal is compared with a null model
#      that removes the same number of randomly chosen NATIVE taxa (Darroch et
#      al. 2014 used range-size comparisons; the random-removal null separates
#      "who" was removed from "how many").
#   2. Origin partition of beta_SIM (Rowan et al. 2024, Eq. 4 with origin as the
#      grouping), and diet x origin to ask how much of the grazer share is
#      immigrant-borne; shared/endemic proportions by origin (Eqs. 5-6).
#   3. Occupancy: provinces and grid cells occupied per taxon and bin, as a
#      share of the cells sampled in that bin; immigrants vs natives by
#      permutation of origin labels within bins (and within size classes), with
#      Cliff's delta, and a quasi-binomial GLM controlling for body mass, diet
#      and bin (occupancy depends on abundance and sampling as well as range).
#   4. First appearances per province with age uncertainty: site ages are drawn
#      uniformly within their Max-Min range (cf. the age resampling of PyRate,
#      Silvestro et al. 2014, used by Tseng et al. 2026), plus classical 95%
#      confidence intervals on the observed range (Strauss & Sadler 1989).
#      Province order is tested against the direction of entry (H4), and the
#      lag from first appearance to peak occupancy is reported.
#   5. Native displacement: change in occupancy of ecologically similar natives
#      (same diet group and size class) from the bin before an immigrant's
#      arrival to its arrival bin, relative to dissimilar natives
#      (difference-in-differences with a random-native null).
#   6. Leave-one-taxon-out contributions to beta_SIM and to the shared-grazer
#      proportion in every bin (H3: Bison in Bin 5).
#
# ORIGIN CODES: the script needs each genus coded Native / South American /
# Eurasian (and optionally Exclude). On the first run it writes
#   Outputs/ch2_immigrants/ch2_taxon_origins_template.csv
# pre-filled from a built-in starting list (see ORIGIN LOOKUP below) and
# "Native (assumed)" for every other genus. CHECK IT against Morgan (2005),
# Woodburne (2010) and Bell et al. (2004), correct it, and save it as
# taxon_origins.csv in work_dir; later runs use that file. Species-level rows
# (e.g. Ursus arctos) override their genus.
#
# SHORTCOMINGS kept in view (see also the README text printed at the end):
#   - few immigrants (Morgan 2005: ~16 South American genera reached temperate
#     North America; xenarthrans dominate, and rodents/notoungulates are not in
#     the large-mammal orders), so province-level tests have little power:
#     the script reports effect sizes and taxon-level permutation tests;
#   - occupancy mixes range, abundance and sampling; cell occupancy is scaled by
#     the cells sampled per bin, and the GLM adds body mass and diet;
#   - many FAUNMAP/PBDB sites have wide age ranges; first appearances use age
#     resampling and drop sites spanning more than max_site_age_span;
#   - origin coding and taxonomy are uncertain: edit the origin table, and
#     compare species- and genus-level results.
#
# Inputs:  Outputs/5b_great_plains/master_data_unique.csv (Step 5b),
#          Outputs/7_rowan_figures/species_traits.csv (Step 7),
#          Outputs/3d_resolved/time_bins.csv, taxon_origins.csv (optional).
# Outputs: Outputs/ch2_immigrants/ (tables and figures per taxonomic level).
# Run the whole file (Ctrl+Shift+S) after Step 7. Needs: dplyr, tidyr, ggplot2.
# Runtime: a few minutes with the default numbers of draws.
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2")) {
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

master_file  <- file.path("Outputs", "5b_great_plains", "master_data_unique.csv")
traits_file  <- file.path("Outputs", "7_rowan_figures", "species_traits.csv")
bins_file    <- file.path("Outputs", "3d_resolved", "time_bins.csv")
origins_file <- "taxon_origins.csv"                  # your checked origin table (in work_dir)
output_dir   <- file.path(work_dir, "Outputs", "ch2_immigrants")

regions <- c("Basin and Range", "Coastal Plain", "Great Plains")
west    <- c("Basin and Range", "Great Plains")      # entry side for Eurasian taxa (H4)
east    <- c("Coastal Plain")                        # entry side for South American taxa (H4)
levels_to_run <- c("species", "genus")

immigrant_since_ma <- 5      # an immigrant lineage counts as an immigrant if it reached North America
                             # after this age (older arrivals, e.g. Megalonyx ~9 Ma, count as established)
count_derived_as_immigrant <- TRUE   # lineages that evolved in North America from a window immigrant
                                     # (e.g. Paramylodon from Glossotherium) count with their origin
exclude_historic <- TRUE     # drop taxa coded Exclude (humans, domestic / historic introductions)

grid_deg        <- 1         # grid cell size for occupancy (degrees); try 2 as a check
n_draws         <- 999       # equal-species draws per bin and scenario (as Chapter 1)
n_null          <- 499       # random-native removals for the null model
n_draws_null    <- 99        # equal-species draws inside each null removal
n_perm          <- 9999      # label permutations for occupancy tests
n_age_draws     <- 999       # age resamples for first appearances
max_site_age_span <- 1.5     # Myr; sites with wider Max-Min ranges are left out of first appearances
random_seed     <- 2024

# Colours (colourblind-friendly, as Step 7)
origin_cols <- c("Native" = "#009E73", "South American" = "#E69F00", "Eurasian" = "#0072B2")
scen_cols   <- c("All taxa" = "#1F2A44", "Without South American" = "#E69F00",
                 "Without Eurasian" = "#56B4E9", "Natives only" = "#009E73")
prov_cols   <- c("Basin and Range" = "#E69F00", "Coastal Plain" = "#56B4E9", "Great Plains" = "#009E73")
ink <- "#1F2A44"; ink_soft <- "#3E4C6D"; band_fill <- "#F6EBD0"; strip_fill <- "#D6E6F5"
base_size <- 13; fig_dpi <- 300

# -----------------------------------------------------------------------------
# A sourced origin table, taxon_origins.csv (Source and Source_detail columns), is
# saved in the repository: copy it to work_dir and it is used instead of this list.
# ORIGIN LOOKUP: a STARTING POINT only, compiled from the GABI and Beringian
# literature (Morgan 2005; Woodburne 2010; Bell et al. 2004; Lister & Sher 2015;
# MacFadden & Hulbert 2009; Froese et al. 2017). Arrival ages (Ma, first
# appearance in North America) are approximate and NA where uncertain. Verify
# every row before using the results. Genera not listed default to Native.
# Derived = evolved in North America from an immigrant ancestor.
# -----------------------------------------------------------------------------
origin_lookup <- tribble(
  ~Taxon,             ~Origin,          ~Arrival_Ma, ~Derived, ~Note,
  # South American immigrants (GABI)
  "Glyptotherium",    "South American", 2.7,   FALSE, "Woodburne 2010 p. 250: arrived in North America ~2.7 Ma",
  "Holmesina",        "South American", 2.4,   FALSE, "Woodburne 2010 p. 251: entered the US ~2.4 Ma (GABI 1)",
  "Pampatherium",     "South American", 2.4,   FALSE, "Woodburne 2010 p. 251: entered the US from Mexico ~2.4 Ma",
  "Dasypus",          "South American", 2.4,   FALSE, "Woodburne 2010 p. 251: entered the US ~2.4 Ma (GABI 1)",
  "Pachyarmatherium", "South American", 2.4,   FALSE, "Woodburne 2010 p. 251: entered the US ~2.4 Ma (GABI 1)",
  "Glossotherium",    "South American", 3.0,   FALSE, "Woodburne 2010 p. 250: Arizona ~3.0 Ma; South American origin",
  "Paramylodon",      "South American", NA,    TRUE,  "North American descendant of the Glossotherium immigrant",
  "Eremotherium",     "South American", 2.4,   FALSE, "Woodburne 2010 p. 251: entered the US ~2.4 Ma (GABI 1)",
  "Nothrotheriops",   "South American", 1.3,   FALSE, "Woodburne 2010 p. 251: dispersal at 1.3 Ma (Leisey Shell Pit)",
  "Myrmecophaga",     "South American", 1.7,   FALSE, "Woodburne 2010 p. 251: Sonora ~1.6-1.8 Ma",
  "Megalonyx",        "South American", 7,     TRUE,  "Woodburne 2010 p. 248: from Pliometanastes ~7 Ma; established before the window",
  "Hydrochoerus",     "South American", NA,    FALSE, "rodent (not in the large-mammal orders)",
  "Neochoerus",       "South American", NA,    FALSE, "rodent (not in the large-mammal orders)",
  "Erethizon",        "South American", NA,    FALSE, "rodent (not in the large-mammal orders)",
  "Mixotoxodon",      "South American", NA,    FALSE, "notoungulate (not in the large-mammal orders)",
  # Eurasian / Beringian immigrants
  "Mammuthus",        "Eurasian",       1.4,   FALSE, "~1.5-1.3 Ma (Lister & Sher 2015); Florida 2.5 Ma record reworked (MacFadden & Hulbert 2009)",
  "Bison",            "Eurasian",       0.195, FALSE, "MIS 6, 195-135 ka (Froese et al. 2017)",
  "Ovibos",           "Eurasian",       NA,    FALSE, "verify arrival",
  "Praeovibos",       "Eurasian",       NA,    FALSE, "Irvingtonian; verify",
  "Soergelia",        "Eurasian",       NA,    FALSE, "Irvingtonian; verify",
  "Euceratherium",    "Eurasian",       NA,    FALSE, "shrub-ox; verify lineage and arrival",
  "Bootherium",       "Eurasian",       NA,    TRUE,  "North American ovibovine of Eurasian stock; verify",
  "Cervalces",        "Eurasian",       NA,    FALSE, "verify arrival",
  "Alces",            "Eurasian",       NA,    FALSE, "late Rancholabrean; verify",
  "Cervus",           "Eurasian",       NA,    FALSE, "Rancholabrean; verify",
  "Rangifer",         "Eurasian",       NA,    FALSE, "verify arrival",
  "Ovis",             "Eurasian",       NA,    FALSE, "Rancholabrean; verify",
  "Saiga",            "Eurasian",       NA,    FALSE, "Beringia, Rancholabrean",
  "Oreamnos",         "Eurasian",       NA,    TRUE,  "North American caprine of Old World stock; verify",
  "Panthera",         "Eurasian",       NA,    FALSE, "verify arrival (P. onca / P. atrox lineages)",
  "Gulo",             "Eurasian",       NA,    FALSE, "Irvingtonian; verify",
  "Xenocyon",         "Eurasian",       NA,    FALSE, "Irvingtonian; verify",
  "Cuon",             "Eurasian",       NA,    FALSE, "Rancholabrean; verify",
  "Ursus arctos",     "Eurasian",       NA,    FALSE, "late Rancholabrean (species row overrides genus)",
  "Canis lupus",      "Eurasian",       NA,    FALSE, "species row overrides the native genus Canis",
  "Vulpes lagopus",   "Eurasian",       NA,    FALSE, "Arctic fox; verify",
  "Alopex lagopus",   "Eurasian",       NA,    FALSE, "Arctic fox (older name); verify",
  # Humans and historic introductions: not part of the Pleistocene fauna
  "Homo",             "Exclude",        NA,    FALSE, "humans",
  "Bos",              "Exclude",        NA,    FALSE, "domestic cattle (historic)",
  "Sus",              "Exclude",        NA,    FALSE, "domestic pig (historic)",
  "Capra",            "Exclude",        NA,    FALSE, "domestic goat (historic)",
  "Felis catus",      "Exclude",        NA,    FALSE, "domestic cat (historic)",
  "Canis familiaris", "Exclude",        NA,    FALSE, "domestic dog",
  "Ovis aries",       "Exclude",        NA,    FALSE, "domestic sheep (historic)",
  "Equus asinus",     "Exclude",        NA,    FALSE, "donkey (historic)"
)

# =============================================================================
# HELPERS
# =============================================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_csv_wd <- function(f) {
  path <- if (file.exists(f)) f else file.path(work_dir, f)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun the earlier steps first.")
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
}
save_fig <- function(p, name, w, h, dir) {
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = fig_dpi, bg = "white")
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}
theme_ch2 <- function() {
  theme_classic(base_size = base_size) +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink_soft),
          axis.line = element_line(colour = ink), axis.ticks = element_line(colour = ink),
          strip.background = element_rect(fill = strip_fill, colour = ink_soft),
          strip.text = element_text(colour = ink), legend.position = "bottom",
          plot.caption = element_text(colour = ink_soft, hjust = 0, size = base_size * 0.7))
}
cliffs_delta <- function(x, y) {          # P(x > y) - P(x < y)
  if (!length(x) || !length(y)) return(NA_real_)
  mean(outer(x, y, ">")) - mean(outer(x, y, "<"))
}
two_sided_p <- function(null, obs) {
  null <- null[is.finite(null)]
  if (!length(null) || !is.finite(obs)) return(NA_real_)
  min(1, 2 * min((sum(null <= obs) + 1) / (length(null) + 1), (sum(null >= obs) + 1) / (length(null) + 1)))
}

# ---- Rowan et al. (2024) equations (as Step 7) --------------------------------
beta_sim <- function(x) {
  x <- x[, colSums(x) > 0, drop = FALSE]
  if (ncol(x) == 0) return(NA_real_)
  cmb <- combn(nrow(x), 2)
  mins <- vapply(seq_len(ncol(cmb)), function(k) {
    i <- cmb[1, k]; j <- cmb[2, k]
    min(sum(x[i, ] == 1 & x[j, ] == 0), sum(x[i, ] == 0 & x[j, ] == 1))
  }, numeric(1))
  shared <- sum(x) - ncol(x)
  if (shared + sum(mins) == 0) return(NA_real_)
  sum(mins) / (shared + sum(mins))
}
min_b <- function(x) {
  b_ij <- x[1, ] == 1 & x[2, ] == 0; b_ji <- x[1, ] == 0 & x[2, ] == 1
  sp <- if (sum(b_ij) < sum(b_ji)) colnames(x)[b_ij] else colnames(x)[b_ji]
  list(min_b = min(sum(b_ij), sum(b_ji)), species = sp)
}
endemic_n <- function(x) {
  cmb <- combn(nrow(x), 2)
  res <- lapply(seq_len(ncol(cmb)), function(i) min_b(x[cmb[, i], , drop = FALSE]))
  list(sum_b = sum(vapply(res, `[[`, numeric(1), "min_b")), species = unlist(lapply(res, `[[`, "species")))
}
shared_n <- function(x) sum(x) - sum(colSums(x) > 0)
simp_func <- function(x, grp) {                         # Eq. 4: beta_SIM split by group
  e <- endemic_n(x); f <- factor(grp[match(e$species, colnames(x))], levels = levels(grp))
  table(f) / (shared_n(x) + e$sum_b)
}
shared_func <- function(x, grp) {                       # Eq. 6: shared proportion by group
  g <- factor(grp[colnames(x)], levels = levels(grp))
  s_if <- rowSums(vapply(seq_len(nrow(x)), function(i) as.numeric(table(g[x[i, ] == 1])), numeric(nlevels(g))))
  s_tf <- as.numeric(table(g[colSums(x) > 0]))
  setNames((s_if - s_tf) / sum(x), levels(g))
}
pa_matrix <- function(d) {                              # d: Spatial_Bin, Taxon
  regs <- intersect(regions, unique(d$Spatial_Bin)); spp <- sort(unique(d$Taxon))
  x <- matrix(0L, length(regs), length(spp), dimnames = list(regs, spp))
  x[cbind(match(d$Spatial_Bin, regs), match(d$Taxon, spp))] <- 1L
  x
}
# Equal-species draws: every province cut to the poorest province's taxon count.
equal_draws <- function(x, n_rep, keep_picks = FALSE) {
  lists <- lapply(seq_len(nrow(x)), function(i) colnames(x)[x[i, ] == 1]); names(lists) <- rownames(x)
  n <- min(lengths(lists))
  if (n < 1 || length(lists) < 2) return(list(beta = rep(NA_real_, n_rep), picks = NULL, n = n))
  picks <- vector("list", n_rep); out <- numeric(n_rep)
  for (r in seq_len(n_rep)) {
    pk <- lapply(lists, function(s) if (length(s) > n) sample(s, n) else s)
    spp <- sort(unique(unlist(pk)))
    m <- t(vapply(pk, function(s) as.integer(spp %in% s), integer(length(spp)))); colnames(m) <- spp
    out[r] <- beta_sim(m)
    if (keep_picks) picks[[r]] <- m
  }
  list(beta = out, picks = if (keep_picks) picks else NULL, n = n)
}

# =============================================================================
# 1. INPUTS
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
md <- read_csv_wd(master_file)
traits <- read_csv_wd(traits_file)
time_bins <- read_csv_wd(bins_file) %>%
  transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma), Younger = as.numeric(Younger_Ma)) %>%
  arrange(Time_Bin) %>% mutate(Mid = (Older + Younger) / 2)
win_old <- max(time_bins$Older); win_young <- min(time_bins$Younger)

large <- traits$Species[as.logical(traits$Included) %in% TRUE]
occ_all <- md %>%
  mutate(Time_Bin = as.integer(Time_Bin), Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude),
         Max_Ma = as.numeric(Max_Ma), Min_Ma = as.numeric(Min_Ma)) %>%
  filter(Spatial_Bin %in% regions, Species %in% large, Time_Bin %in% time_bins$Time_Bin)
if (!"Genus" %in% names(occ_all)) occ_all$Genus <- sub(" .*$", "", occ_all$Species)
cat(sprintf("  %d large-mammal occurrences (species x site) in the three provinces\n", nrow(occ_all)))

# ---- Origin table ---------------------------------------------------------------
genera  <- sort(unique(occ_all$Genus))
species <- sort(unique(occ_all$Species))
orders  <- traits %>% transmute(Genus = sub(" .*$", "", Species), Order) %>% distinct(Genus, .keep_all = TRUE)
template <- data.frame(Taxon = genera, stringsAsFactors = FALSE) %>%
  left_join(origin_lookup, by = "Taxon") %>%
  mutate(Basis = ifelse(is.na(Origin), "assumed native (not in the built-in list) - check",
                        "built-in starting list - verify"),
         Origin = ifelse(is.na(Origin), "Native", Origin), Derived = ifelse(is.na(Derived), FALSE, Derived)) %>%
  left_join(orders, by = c("Taxon" = "Genus"))
sp_rows <- origin_lookup %>% filter(grepl(" ", Taxon), Taxon %in% species) %>%
  mutate(Basis = "built-in starting list - verify", Order = NA_character_)
template <- bind_rows(template, sp_rows) %>%
  select(Taxon, Order, Origin, Arrival_Ma, Derived, Basis, Note) %>% arrange(Origin, Taxon)
write.csv(template, file.path(output_dir, "ch2_taxon_origins_template.csv"), row.names = FALSE, na = "")

origins_path <- file.path(work_dir, origins_file)
if (file.exists(origins_path)) {
  origins <- read.csv(origins_path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  cat(sprintf("  Origin codes from %s (%d rows)\n", origins_file, nrow(origins)))
} else {
  origins <- template
  cat("  NOTE: no taxon_origins.csv yet - using the built-in starting list.\n",
      "       Check Outputs/ch2_immigrants/ch2_taxon_origins_template.csv, correct it and save it\n",
      "       as taxon_origins.csv in work_dir before trusting the results.\n")
}
origins <- origins %>% mutate(Derived = as.logical(Derived), Arrival_Ma = as.numeric(Arrival_Ma),
                              Origin = trimws(Origin))
bad <- setdiff(unique(origins$Origin), c("Native", "South American", "Eurasian", "Exclude"))
if (length(bad)) stop("Unknown Origin value(s) in the origin table: ", paste(bad, collapse = ", "),
                      "\nUse Native, South American, Eurasian or Exclude.")

# Group used in the analyses: window immigrants keep their origin; earlier arrivals count as Native.
assign_group <- function(o) {
  window <- is.na(o$Arrival_Ma) | o$Arrival_Ma <= immigrant_since_ma
  derived_ok <- count_derived_as_immigrant | !o$Derived
  ifelse(o$Origin == "Exclude", "Exclude",
         ifelse(o$Origin %in% c("South American", "Eurasian") & window & derived_ok, o$Origin, "Native"))
}
origins$Group <- assign_group(origins)
genus_group   <- setNames(origins$Group[!grepl(" ", origins$Taxon)], origins$Taxon[!grepl(" ", origins$Taxon)])
species_group <- setNames(origins$Group[grepl(" ", origins$Taxon)], origins$Taxon[grepl(" ", origins$Taxon)])
group_of_species <- function(sp) {
  g <- unname(species_group[sp]); gg <- unname(genus_group[sub(" .*$", "", sp)])
  ifelse(!is.na(g), g, ifelse(!is.na(gg), gg, "Native"))
}
occ_all$Group_sp <- group_of_species(occ_all$Species)
occ_all$Group_ge <- ifelse(is.na(genus_group[occ_all$Genus]), "Native", unname(genus_group[occ_all$Genus]))

excluded <- occ_all %>% filter(Group_sp == "Exclude") %>% distinct(Species)
if (nrow(excluded)) {
  cat(sprintf("  Excluded as humans / historic introductions: %s\n", paste(excluded$Species, collapse = ", ")))
  cat("  NOTE: these taxa are also in the Chapter 1 matrices; consider removing them there too.\n")
}
if (exclude_historic) occ_all <- occ_all %>% filter(Group_sp != "Exclude", Group_ge != "Exclude")
occ_all <- occ_all %>% mutate(Cell = paste(floor(Latitude / grid_deg), floor(Longitude / grid_deg), sep = "_"))

# Traits per species; per genus: diet by majority, mass by median.
tr_sp <- traits %>% filter(Species %in% occ_all$Species) %>%
  transmute(Taxon = Species, Diet = Diet_Group, Size = Size_Class, Mass_kg = as.numeric(Mass_kg))
tr_ge <- tr_sp %>% mutate(Genus = sub(" .*$", "", Taxon)) %>% group_by(Genus) %>%
  summarise(Diet = { d <- na.omit(Diet); if (length(d)) names(sort(table(d), decreasing = TRUE))[1] else NA_character_ },
            Mass_kg = if (all(is.na(Mass_kg))) NA_real_ else median(Mass_kg, na.rm = TRUE), .groups = "drop") %>%
  mutate(Size = as.character(cut(Mass_kg, c(0, 18, 80, 350, 1000, Inf), right = FALSE,
                                 labels = c("Size 1 (<18 kg)", "Size 2 (18-80 kg)", "Size 3 (80-350 kg)",
                                            "Size 4 (350-1,000 kg)", "Size 5 (>1,000 kg)")))) %>%
  rename(Taxon = Genus)

# =============================================================================
# 2. ANALYSES PER TAXONOMIC LEVEL
# =============================================================================

run_level <- function(level) {
  cat(sprintf("\n================ %s LEVEL ================\n", toupper(level)))
  out_dir <- file.path(output_dir, level); dir.create(out_dir, showWarnings = FALSE)
  set.seed(random_seed)
  occ <- occ_all %>% mutate(Taxon = if (level == "species") Species else Genus,
                            Group = if (level == "species") Group_sp else Group_ge)
  tr <- if (level == "species") tr_sp else tr_ge
  grp_tab <- occ %>% distinct(Taxon, Group)
  group_vec <- setNames(grp_tab$Group, grp_tab$Taxon)
  diet_vec  <- setNames(tr$Diet, tr$Taxon)
  writeLines(sprintf("Taxonomic level: %s", level), file.path(out_dir, "README.txt"))

  # ---- counts ----------------------------------------------------------------
  counts <- occ %>% distinct(Time_Bin, Spatial_Bin, Taxon, Group) %>%
    count(Time_Bin, Spatial_Bin, Group, name = "n_taxa") %>%
    complete(Time_Bin, Spatial_Bin, Group = names(origin_cols), fill = list(n_taxa = 0))
  write.csv(counts, file.path(out_dir, "ch2_taxa_per_bin_province_origin.csv"), row.names = FALSE)
  cat("  Immigrant taxa per bin (all provinces):\n")
  print(as.data.frame(occ %>% distinct(Time_Bin, Taxon, Group) %>% count(Time_Bin, Group) %>%
                        pivot_wider(names_from = Group, values_from = n, values_fill = 0)), row.names = FALSE)

  pa <- lapply(setNames(time_bins$Time_Bin, time_bins$Time_Bin), function(b)
    pa_matrix(occ %>% filter(Time_Bin == b) %>% distinct(Spatial_Bin, Taxon)))

  # ---- 2.1 origin effects on beta_SIM, with random-native null ---------------
  cat("\n  --- 2.1 Origin effects on beta_SIM ---\n")
  scenarios <- list("All taxa" = character(0), "Without South American" = "South American",
                    "Without Eurasian" = "Eurasian", "Natives only" = c("South American", "Eurasian"))
  scen_res <- list(); null_res <- list(); all_picks <- list()
  for (b in names(pa)) {
    x <- pa[[b]]
    if (nrow(x) < 2) next
    g <- group_vec[colnames(x)]
    full <- equal_draws(x, n_draws, keep_picks = TRUE)
    all_picks[[b]] <- full$picks
    for (s in names(scenarios)) {
      drop <- colnames(x)[g %in% scenarios[[s]]]
      xs <- x[, setdiff(colnames(x), drop), drop = FALSE]
      ed <- if (s == "All taxa" || length(drop) == 0) full else equal_draws(xs, n_draws)
      obs_delta <- mean(ed$beta) - mean(full$beta)
      null_delta <- rep(NA_real_, n_null)
      natives <- colnames(x)[g == "Native"]
      if (length(drop) > 0 && length(drop) <= length(natives) - 1) {
        null_delta <- vapply(seq_len(n_null), function(r) {
          xr <- x[, setdiff(colnames(x), sample(natives, length(drop))), drop = FALSE]
          mean(equal_draws(xr, n_draws_null)$beta) - mean(full$beta)
        }, numeric(1))
      }
      scen_res[[length(scen_res) + 1]] <- data.frame(
        Time_Bin = as.integer(b), Scenario = s, n_removed = length(drop), n_used = ed$n,
        beta_mean = mean(ed$beta), beta_lower95 = unname(quantile(ed$beta, 0.025, na.rm = TRUE)),
        beta_upper95 = unname(quantile(ed$beta, 0.975, na.rm = TRUE)),
        delta = obs_delta,
        null_mean = mean(null_delta, na.rm = TRUE),
        null_lower95 = if (all(is.na(null_delta))) NA else unname(quantile(null_delta, 0.025, na.rm = TRUE)),
        null_upper95 = if (all(is.na(null_delta))) NA else unname(quantile(null_delta, 0.975, na.rm = TRUE)),
        SES = if (all(is.na(null_delta)) || sd(null_delta, na.rm = TRUE) == 0) NA else
          (obs_delta - mean(null_delta, na.rm = TRUE)) / sd(null_delta, na.rm = TRUE),
        p_vs_random_natives = if (s == "All taxa" || length(drop) == 0) NA else two_sided_p(null_delta, obs_delta))
    }
  }
  scen_res <- bind_rows(scen_res) %>% left_join(time_bins, by = "Time_Bin") %>%
    group_by(Scenario) %>% mutate(p_holm = p.adjust(p_vs_random_natives, method = "holm")) %>% ungroup()
  write.csv(scen_res, file.path(out_dir, "ch2_beta_by_origin_scenario.csv"), row.names = FALSE)
  print(as.data.frame(scen_res %>% filter(Scenario != "All taxa") %>%
                        transmute(Bin = Time_Bin, Scenario, removed = n_removed, beta = round(beta_mean, 3),
                                  delta = sprintf("%+.3f", delta), null = sprintf("%+.3f", null_mean),
                                  SES = round(SES, 2), p = round(p_vs_random_natives, 3), p_holm = round(p_holm, 3))),
        row.names = FALSE)

  # ---- 2.2 origin partition of beta_SIM and diet x origin ---------------------
  cat("\n  --- 2.2 Origin partition (Rowan Eq. 4) ---\n")
  part <- list()
  for (b in names(all_picks)) {
    picks <- all_picks[[b]]
    og <- factor(group_vec, levels = names(origin_cols)); names(og) <- names(group_vec)
    dg <- factor(paste(ifelse(is.na(diet_vec[names(group_vec)]), "No diet", diet_vec[names(group_vec)]), group_vec, sep = " | "))
    names(dg) <- names(group_vec)
    eq_o <- rowMeans(vapply(picks, function(m) as.numeric(simp_func(m, og)), numeric(nlevels(og))))
    eq_d <- rowMeans(vapply(picks, function(m) as.numeric(simp_func(m, dg)), numeric(nlevels(dg))))
    x <- pa[[b]]
    sh <- shared_func(x, og)
    part[[length(part) + 1]] <- bind_rows(
      data.frame(Time_Bin = as.integer(b), Partition = "origin", Group = levels(og), beta_SIM_equal = eq_o,
                 beta_SIM_full = as.numeric(simp_func(x, og)), shared_prop_full = as.numeric(sh)),
      data.frame(Time_Bin = as.integer(b), Partition = "diet x origin", Group = levels(dg), beta_SIM_equal = eq_d,
                 beta_SIM_full = as.numeric(simp_func(x, dg)), shared_prop_full = NA))
  }
  part <- bind_rows(part)
  write.csv(part, file.path(out_dir, "ch2_beta_partition_by_origin.csv"), row.names = FALSE)
  grazer <- part %>% filter(Partition == "diet x origin", grepl("^Grazer", Group)) %>%
    group_by(Time_Bin) %>% mutate(share_of_grazer_beta = beta_SIM_equal / sum(beta_SIM_equal)) %>% ungroup()
  write.csv(grazer, file.path(out_dir, "ch2_grazer_beta_by_origin.csv"), row.names = FALSE)
  print(as.data.frame(grazer %>% transmute(Bin = Time_Bin, Group, beta = round(beta_SIM_equal, 3),
                                           share = round(share_of_grazer_beta, 2))), row.names = FALSE)

  # ---- 2.3 occupancy ---------------------------------------------------------
  cat("\n  --- 2.3 Occupancy: immigrants vs natives ---\n")
  sampled <- occ %>% distinct(Time_Bin, Spatial_Bin, Cell) %>% group_by(Time_Bin) %>%
    summarise(n_prov_sampled = n_distinct(Spatial_Bin), n_cells_sampled = n_distinct(Cell), .groups = "drop")
  occu <- occ %>% group_by(Time_Bin, Taxon, Group) %>%
    summarise(n_provinces = n_distinct(Spatial_Bin), n_cells = n_distinct(Cell), n_sites = n_distinct(Site_Key),
              .groups = "drop") %>%
    left_join(sampled, by = "Time_Bin") %>%
    mutate(prop_cells = n_cells / n_cells_sampled, prop_provinces = n_provinces / n_prov_sampled) %>%
    left_join(tr, by = "Taxon") %>%
    group_by(Taxon) %>% mutate(First_Bin = min(Time_Bin), Arrival_Bin = Time_Bin == First_Bin) %>% ungroup()
  write.csv(occu, file.path(out_dir, "ch2_occupancy_by_taxon_bin.csv"), row.names = FALSE)

  perm_test <- function(d, imm, stat_col, strata = c("Time_Bin")) {
    d <- d %>% filter(Group %in% c("Native", imm))
    if (sum(d$Group == imm) < 1 || sum(d$Group == "Native") < 1) return(NULL)
    stat <- function(lab) {
      per <- tapply(seq_len(nrow(d)), d$Time_Bin, function(ix) {
        a <- d[[stat_col]][ix][lab[ix] == imm]; n <- d[[stat_col]][ix][lab[ix] == "Native"]
        if (!length(a) || !length(n)) NA else mean(a) - mean(n)
      })
      w <- tapply(lab == imm, d$Time_Bin, sum)
      ok <- !is.na(per); sum(per[ok] * w[ok]) / sum(w[ok])        # weighted by immigrants per bin
    }
    obs <- stat(d$Group)
    key <- interaction(d[strata], drop = TRUE)
    null <- replicate(n_perm, { lab <- d$Group
      for (k in levels(key)) { ix <- which(key == k); lab[ix] <- lab[ix][sample.int(length(ix))] }
      stat(lab) })
    data.frame(Immigrant_group = imm, Statistic = stat_col, Strata = paste(strata, collapse = " x "),
               n_immigrant_taxon_bins = sum(d$Group == imm), mean_immigrant = mean(d[[stat_col]][d$Group == imm]),
               mean_native = mean(d[[stat_col]][d$Group == "Native"]), diff_weighted = obs,
               cliffs_delta = cliffs_delta(d[[stat_col]][d$Group == imm], d[[stat_col]][d$Group == "Native"]),
               p_perm = two_sided_p(null, obs))
  }
  occu_s <- occu %>% mutate(Size_strata = ifelse(is.na(Size), "No size", Size))
  occ_tests <- bind_rows(lapply(c("South American", "Eurasian"), function(imm) bind_rows(
    perm_test(occu_s, imm, "prop_cells"), perm_test(occu_s, imm, "prop_provinces"),
    perm_test(occu_s, imm, "prop_cells", c("Time_Bin", "Size_strata")),
    perm_test(occu_s %>% filter(Arrival_Bin | Group == "Native"), imm, "prop_cells") %>%
      mutate(Strata = "Time_Bin (immigrants in their first bin only)"))))
  write.csv(occ_tests, file.path(out_dir, "ch2_occupancy_tests.csv"), row.names = FALSE)
  if (nrow(occ_tests)) print(as.data.frame(occ_tests %>% transmute(Immigrant_group, Statistic, Strata,
    n = n_immigrant_taxon_bins, imm = round(mean_immigrant, 3), nat = round(mean_native, 3),
    diff = round(diff_weighted, 3), delta = round(cliffs_delta, 2), p = round(p_perm, 4))), row.names = FALSE)

  glm_d <- occu %>% filter(!is.na(Mass_kg), !is.na(Diet)) %>%
    mutate(Group = factor(Group, levels = names(origin_cols)), Bin = factor(Time_Bin), log10_mass = log10(Mass_kg))
  glm_out <- NULL
  if (nrow(glm_d) > 10 && n_distinct(glm_d$Group) > 1) {
    fit <- glm(cbind(n_cells, n_cells_sampled - n_cells) ~ Group + log10_mass + Diet + Bin,
               family = quasibinomial(), data = droplevels(glm_d))
    glm_out <- as.data.frame(summary(fit)$coefficients); glm_out$Term <- rownames(glm_out)
    write.csv(glm_out[, c("Term", names(glm_out)[1:4])], file.path(out_dir, "ch2_occupancy_glm.csv"), row.names = FALSE)
    cat("  Quasi-binomial GLM of cell occupancy (origin terms):\n")
    print(glm_out[grepl("^Group", glm_out$Term), c("Term", "Estimate", "Std. Error", "Pr(>|t|)")], row.names = FALSE)
  }

  # ---- 2.4 first appearances per province -------------------------------------
  cat("\n  --- 2.4 First appearances per province (age resampling) ---\n")
  fa_d <- occ %>% filter(Group != "Native") %>%
    distinct(Taxon, Group, Spatial_Bin, Site_Key, Max_Ma, Min_Ma) %>%
    filter(!is.na(Max_Ma), !is.na(Min_Ma), (Max_Ma - Min_Ma) <= max_site_age_span)
  fa_rows <- list(); order_rows <- list()
  for (tx in unique(fa_d$Taxon)) {
    dt <- fa_d %>% filter(Taxon == tx)
    provs <- intersect(regions, unique(dt$Spatial_Bin))
    draws <- sapply(provs, function(p) {
      dp <- dt %>% filter(Spatial_Bin == p)
      replicate(n_age_draws, max(runif(nrow(dp), dp$Min_Ma, dp$Max_Ma)))
    })
    draws <- matrix(draws, ncol = length(provs), dimnames = list(NULL, provs))
    for (p in provs) {
      dp <- dt %>% filter(Spatial_Bin == p); mids <- (dp$Max_Ma + dp$Min_Ma) / 2; H <- length(unique(mids))
      ext <- if (H >= 2) diff(range(mids)) * ((1 - 0.95)^(-1 / (H - 1)) - 1) else NA_real_
      fa_rows[[length(fa_rows) + 1]] <- data.frame(
        Taxon = tx, Group = dt$Group[1], Province = p, n_sites = nrow(dp),
        FA_median = median(draws[, p]), FA_lower95 = unname(quantile(draws[, p], 0.025)),
        FA_upper95 = unname(quantile(draws[, p], 0.975)),
        FA_observed_mid = max(mids), FA_StraussSadler95 = if (is.na(ext)) NA else max(mids) + ext)
    }
    pred <- if (dt$Group[1] == "Eurasian") list(early = west, late = east) else list(early = east, late = west)
    for (a in intersect(pred$early, provs)) for (l in intersect(pred$late, provs)) {
      order_rows[[length(order_rows) + 1]] <- data.frame(
        Taxon = tx, Group = dt$Group[1], Predicted_earlier = a, Predicted_later = l,
        prob_predicted_order = mean(draws[, a] > draws[, l]), median_lag_Myr = median(draws[, a] - draws[, l]))
    }
  }
  fa_tab <- bind_rows(fa_rows); ord_tab <- bind_rows(order_rows)
  if (nrow(fa_tab)) {
    # Records older than the literature arrival: possible misidentification, reworking
    # (cf. MacFadden & Hulbert 2009 for Florida Mammuthus) or a wrong site age. Check them.
    arr <- origins %>% filter(!is.na(Arrival_Ma)) %>% select(Taxon, Arrival_Ma)
    lit <- function(tx) { a <- arr$Arrival_Ma[match(tx, arr$Taxon)]
                          ifelse(is.na(a), arr$Arrival_Ma[match(sub(" .*$", "", tx), arr$Taxon)], a) }
    fa_tab <- fa_tab %>% mutate(Literature_Arrival_Ma = lit(Taxon),
                                Check_record = !is.na(Literature_Arrival_Ma) & FA_lower95 > Literature_Arrival_Ma + 0.1)
    if (any(fa_tab$Check_record)) {
      cat("  CHECK: first appearances older than the literature arrival (whole 95% interval):\n")
      print(as.data.frame(fa_tab %>% filter(Check_record) %>%
                            transmute(Taxon, Province, n_sites, FA_median = round(FA_median, 2),
                                      FA_lower95 = round(FA_lower95, 2), Literature_Arrival_Ma)), row.names = FALSE)
      flagged_sites <- fa_d %>% semi_join(fa_tab %>% filter(Check_record) %>% select(Taxon, Spatial_Bin = Province),
                                          by = c("Taxon", "Spatial_Bin")) %>%
        mutate(Literature_Arrival_Ma = lit(Taxon)) %>% filter(Max_Ma > Literature_Arrival_Ma + 0.1) %>%
        left_join(occ %>% distinct(Site_Key, Site, Source), by = "Site_Key") %>% arrange(Taxon, desc(Max_Ma))
      write.csv(flagged_sites, file.path(out_dir, "ch2_CHECK_records_older_than_arrival.csv"), row.names = FALSE)
    }
    peak <- occu %>% group_by(Taxon) %>% slice_max(prop_cells, n = 1, with_ties = FALSE) %>% ungroup() %>%
      left_join(time_bins %>% select(Time_Bin, Peak_Mid = Mid), by = "Time_Bin") %>% select(Taxon, Peak_Bin = Time_Bin, Peak_Mid)
    lag_tab <- fa_tab %>% group_by(Taxon, Group) %>% summarise(FA_continent = max(FA_median), .groups = "drop") %>%
      left_join(peak, by = "Taxon") %>% mutate(Lag_FA_to_peak_Myr = FA_continent - Peak_Mid)
    write.csv(fa_tab, file.path(out_dir, "ch2_first_appearances_by_province.csv"), row.names = FALSE)
    write.csv(ord_tab, file.path(out_dir, "ch2_first_appearance_order_tests.csv"), row.names = FALSE)
    write.csv(lag_tab, file.path(out_dir, "ch2_first_appearance_to_peak_lag.csv"), row.names = FALSE)
    if (nrow(ord_tab)) {
      summ <- ord_tab %>% group_by(Group) %>%
        summarise(n_pairs = n(), n_taxa = n_distinct(Taxon), n_supporting = sum(prob_predicted_order > 0.5),
                  mean_prob = mean(prob_predicted_order), .groups = "drop") %>%
        rowwise() %>% mutate(p_binomial = binom.test(n_supporting, n_pairs, 0.5)$p.value) %>% ungroup()
      write.csv(summ, file.path(out_dir, "ch2_first_appearance_order_summary.csv"), row.names = FALSE)
      cat("  H4 (direction of entry): taxon-province pairs supporting the predicted order\n")
      print(as.data.frame(summ), row.names = FALSE)
      cat("  (pairs within a taxon are not independent; treat the binomial p as descriptive)\n")
    }
  }

  # ---- 2.5 native displacement ---------------------------------------------
  cat("\n  --- 2.5 Native displacement (difference-in-differences) ---\n")
  occ_wide <- occu %>% select(Taxon, Time_Bin, prop_cells) %>%
    complete(Taxon, Time_Bin = time_bins$Time_Bin, fill = list(prop_cells = 0))
  nat <- tr %>% filter(Taxon %in% names(group_vec)[group_vec == "Native"], !is.na(Diet))
  disp <- list()
  imms <- occu %>% filter(Group != "Native") %>% distinct(Taxon, Group, First_Bin) %>% filter(First_Bin > min(time_bins$Time_Bin))
  for (i in seq_len(nrow(imms))) {
    tx <- imms$Taxon[i]; a <- imms$First_Bin[i]; ti <- tr %>% filter(Taxon == tx)
    if (!nrow(ti) || is.na(ti$Diet)) next
    sim <- nat$Taxon[nat$Diet == ti$Diet & (is.na(ti$Size) | nat$Size %in% ti$Size)]
    if (length(sim) < 1) sim <- nat$Taxon[nat$Diet == ti$Diet]
    others <- setdiff(nat$Taxon, sim)
    if (length(sim) < 1 || length(others) < 1) next
    chg <- function(set) {
      before <- occ_wide %>% filter(Taxon %in% set, Time_Bin == a - 1) %>% pull(prop_cells)
      after  <- occ_wide %>% filter(Taxon %in% set, Time_Bin == a) %>% pull(prop_cells)
      mean(after) - mean(before)
    }
    did <- chg(sim) - chg(others)
    null <- replicate(n_null, { s <- sample(nat$Taxon, length(sim)); chg(s) - chg(setdiff(nat$Taxon, s)) })
    disp[[length(disp) + 1]] <- data.frame(Immigrant = tx, Group = imms$Group[i], Arrival_Bin = a,
                                           Diet = ti$Diet, Size = ti$Size, n_similar_natives = length(sim),
                                           change_similar = chg(sim), change_other = chg(others),
                                           DiD = did, p_random_natives = two_sided_p(null, did))
  }
  disp <- bind_rows(disp)
  write.csv(disp, file.path(out_dir, "ch2_native_displacement.csv"), row.names = FALSE)
  if (nrow(disp)) print(as.data.frame(disp %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

  # ---- 2.6 leave-one-taxon-out contributions ----------------------------------
  cat("\n  --- 2.6 Leave-one-taxon-out contributions ---\n")
  loo <- list()
  for (b in names(pa)) {
    x <- pa[[b]]; if (nrow(x) < 2) next
    gz <- factor(ifelse(diet_vec[colnames(x)] %in% "Grazer", "Grazer", "Other"), levels = c("Grazer", "Other"))
    names(gz) <- colnames(x)
    b0 <- beta_sim(x); s0 <- shared_func(x, gz)[["Grazer"]]
    for (tx in colnames(x)) {
      xr <- x[, colnames(x) != tx, drop = FALSE]
      loo[[length(loo) + 1]] <- data.frame(
        Time_Bin = as.integer(b), Taxon = tx, Group = group_vec[[tx]], Diet = diet_vec[tx],
        n_provinces = sum(x[, tx]), beta_full = b0, delta_beta_if_removed = beta_sim(xr) - b0,
        shared_grazer_full = s0, delta_shared_grazer_if_removed = shared_func(xr, gz[colnames(xr)])[["Grazer"]] - s0)
    }
  }
  loo <- bind_rows(loo) %>% arrange(Time_Bin, delta_shared_grazer_if_removed)
  write.csv(loo, file.path(out_dir, "ch2_leave_one_taxon_out.csv"), row.names = FALSE)
  bis <- loo %>% filter(grepl("^Bison", Taxon))
  if (nrow(bis)) {
    cat("  H3 - Bison: change in shared-grazer proportion and beta_SIM if removed\n")
    print(as.data.frame(bis %>% transmute(Bin = Time_Bin, Taxon, provinces = n_provinces,
                                          d_shared_grazer = round(delta_shared_grazer_if_removed, 3),
                                          d_beta = round(delta_beta_if_removed, 3))), row.names = FALSE)
  }

  # ---- figures -----------------------------------------------------------------
  x_scale <- scale_x_reverse(breaks = c(3, 2, 1, 0))
  bands <- geom_rect(data = time_bins[seq(1, nrow(time_bins), 2), ], inherit.aes = FALSE,
                     aes(xmin = Younger, xmax = Older, ymin = -Inf, ymax = Inf), fill = band_fill)
  cap <- sprintf("%s level. Equal-species beta_SIM (%d draws); intervals = 95%% of draws.\nImmigrant groups from %s.",
                 tools::toTitleCase(level), n_draws,
                 if (file.exists(origins_path)) "taxon_origins.csv" else "the BUILT-IN STARTING ORIGIN CODES (unchecked)")

  d1 <- scen_res %>% mutate(Scenario = factor(Scenario, levels = names(scen_cols)))
  p1 <- ggplot(d1, aes(Mid, beta_mean, colour = Scenario)) + bands +
    geom_line(linewidth = 0.6, position = position_dodge(width = 0.12)) +
    geom_errorbar(aes(ymin = beta_lower95, ymax = beta_upper95), width = 0, linewidth = 0.5, position = position_dodge(width = 0.12)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.12)) +
    scale_colour_manual(values = scen_cols, name = NULL) + x_scale +
    labs(x = "Age (Ma)", y = expression(beta[SIM] ~ "(equal species)"), caption = cap) + theme_ch2()
  save_fig(p1, "ch2_fig1_beta_by_origin_scenario", 9, 5.6, out_dir)

  d2 <- scen_res %>% filter(Scenario != "All taxa", !is.na(null_mean)) %>%
    mutate(Scenario = factor(Scenario, levels = names(scen_cols)))
  if (nrow(d2)) {
    p2 <- ggplot(d2, aes(Mid)) + bands +
      geom_hline(yintercept = 0, colour = ink_soft, linewidth = 0.3) +
      geom_errorbar(aes(ymin = null_lower95, ymax = null_upper95), width = 0.12, colour = ink_soft, linewidth = 0.5) +
      geom_point(aes(y = null_mean), shape = 3, colour = ink_soft, size = 2) +
      geom_point(aes(y = delta, fill = Scenario), shape = 21, colour = ink, size = 3) +
      geom_text(aes(y = pmax(delta, null_upper95, na.rm = TRUE), label = ifelse(is.na(p_vs_random_natives), "",
                    sprintf("p = %.3f", p_vs_random_natives))), vjust = -0.8, size = 3, colour = ink_soft) +
      facet_wrap(~ Scenario, nrow = 1) + scale_fill_manual(values = scen_cols, guide = "none") + x_scale +
      labs(x = "Age (Ma)", y = expression(Delta * beta[SIM] ~ "after removal"),
           caption = paste0("Points: observed change. Grey: 95% of the change when the same number of random natives is removed (",
                            n_null, " removals).\np = two-sided, uncorrected (Holm values in the table).")) + theme_ch2()
    save_fig(p2, "ch2_fig2_removal_vs_random_natives", 12, 4.8, out_dir)
  }

  d3 <- occu %>% left_join(time_bins %>% select(Time_Bin, Mid), by = "Time_Bin") %>%
    mutate(Group = factor(Group, levels = names(origin_cols)), Bin = factor(paste("Bin", Time_Bin)))
  p3 <- ggplot(d3, aes(Group, prop_cells, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, colour = ink_soft, width = 0.6) +
    geom_jitter(aes(colour = Group), width = 0.15, height = 0, size = 1.4, alpha = 0.8) +
    facet_wrap(~ Bin, nrow = 1) + scale_fill_manual(values = origin_cols, drop = FALSE, name = NULL) +
    scale_colour_manual(values = origin_cols, drop = FALSE, guide = "none") +
    labs(x = NULL, y = sprintf("Share of sampled %g\u00b0 cells occupied", grid_deg),
         caption = "One point per taxon and bin.") +
    theme_ch2() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  save_fig(p3, "ch2_fig3_occupancy_by_origin", 12, 4.8, out_dir)

  if (nrow(fa_tab)) {
    d4 <- fa_tab %>% mutate(Province = factor(Province, levels = regions),
                            Taxon = reorder(Taxon, FA_median, max))
    p4 <- ggplot(d4, aes(FA_median, Taxon, colour = Province)) +
      geom_errorbarh(aes(xmin = FA_lower95, xmax = FA_upper95), height = 0, linewidth = 0.6,
                     position = position_dodge(width = 0.6)) +
      geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
      facet_grid(Group ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(9)) +
      geom_point(data = d4 %>% filter(Check_record), aes(x = FA_median), shape = 4, size = 4, colour = ink,
                 position = position_dodge(width = 0.6)) +
      scale_colour_manual(values = prov_cols, name = NULL) + scale_x_reverse() +
      labs(x = "First appearance (Ma)", y = NULL,
           caption = sprintf("Median and 95%% of %d age resamples (site ages drawn within Max-Min; sites spanning > %g Myr left out).\nx = older than the literature arrival: check the record (ch2_CHECK_records_older_than_arrival.csv).",
                             n_age_draws, max_site_age_span)) +
      theme_ch2() + theme(axis.text.y = element_text(face = "italic"), strip.text.y = element_text(angle = 0))
    save_fig(p4, "ch2_fig4_first_appearances", 9, 2 + 0.28 * n_distinct(d4$Taxon) * 1.5, out_dir)
  }

  d5 <- part %>% filter(Partition == "origin") %>% left_join(time_bins, by = "Time_Bin") %>%
    mutate(Group = factor(Group, levels = rev(names(origin_cols))))
  p5 <- ggplot(d5, aes(Mid, beta_SIM_equal, fill = Group)) + bands +
    geom_col(width = 0.18, colour = ink_soft, linewidth = 0.2) +
    scale_fill_manual(values = origin_cols, name = NULL) + x_scale +
    labs(x = "Age (Ma)", y = expression(beta[SIM] ~ "(equal species), split by origin"),
         caption = "Rowan et al. (2024) Eq. 4 with origin as the grouping; bars sum to the equal-species beta_SIM.") + theme_ch2()
  save_fig(p5, "ch2_fig5_beta_partition_by_origin", 9, 5.4, out_dir)

  invisible(list(scenarios = scen_res, partition = part, occupancy = occu, occ_tests = occ_tests,
                 first_appearances = fa_tab, displacement = disp, loo = loo))
}

results <- lapply(setNames(levels_to_run, levels_to_run), run_level)

# =============================================================================
# 3. SPECIES vs GENUS COMPARISON AND NOTES
# =============================================================================

if (all(c("species", "genus") %in% names(results))) {
  cmp <- inner_join(results$species$scenarios %>% select(Time_Bin, Scenario, delta_species = delta, p_species = p_vs_random_natives),
                    results$genus$scenarios %>% select(Time_Bin, Scenario, delta_genus = delta, p_genus = p_vs_random_natives),
                    by = c("Time_Bin", "Scenario")) %>%
    filter(Scenario != "All taxa") %>% mutate(same_sign = sign(delta_species) == sign(delta_genus))
  write.csv(cmp, file.path(output_dir, "ch2_species_vs_genus_removal_effects.csv"), row.names = FALSE)
  cat("\n=== SPECIES vs GENUS: removal effects with the same sign ===\n")
  print(as.data.frame(cmp %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
}

notes <- c(
  "CHAPTER 2 - HOW TO READ THE OUTPUTS",
  "",
  "ch2_beta_by_origin_scenario.csv: delta < 0 after removing an origin = that origin RAISED beta_SIM",
  "  (Darroch-type differentiation, H1); delta > 0 = it LOWERED beta_SIM (Stigall-type homogenization, H2).",
  "  p_vs_random_natives asks whether removing these taxa changes beta_SIM more than removing the same",
  "  number of random natives; SES = standardized effect size. Report effect sizes, not only p.",
  "ch2_grazer_beta_by_origin.csv: share of the grazer contribution to beta_SIM carried by each origin (H1).",
  "ch2_occupancy_tests.csv: diff/cliffs_delta < 0 = immigrants less widespread than natives (H1).",
  "ch2_first_appearance_order_tests.csv: prob_predicted_order = share of age resamples in which the",
  "  entry-side province records the taxon first (H4).",
  "ch2_leave_one_taxon_out.csv: delta_shared_grazer_if_removed < 0 = the taxon carried shared-grazer",
  "  structure (H3: Bison in Bin 5).",
  "",
  "SHORTCOMINGS: few immigrant taxa (low power); occupancy mixes range, abundance and sampling;",
  "wide site age ranges; origin codes need checking (taxon_origins.csv); species-genus agreement is",
  "reported in ch2_species_vs_genus_removal_effects.csv. Climate-envelope tests (Oscillayers) and",
  "Florida Museum / Mexican records are not included yet."
)
writeLines(notes, file.path(output_dir, "README_ch2.txt"))
cat("\n", paste(notes, collapse = "\n"), "\n")
cat("\n=== CHAPTER 2 COMPLETE: outputs in", output_dir, "===\n")
