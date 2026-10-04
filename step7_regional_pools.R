# =============================================================================
# STEP 7: REGIONAL SPECIES POOLS AND SIMPSON DISSIMILARITY THROUGH TIME
#
# For four physiographic provinces (Fenneman & Johnson 1946):
#   Pacific Border, Basin and Range, Great Plains, Coastal Plain
#
# 1. Regional species pool: every species recorded at any site in the province
#    during a geological stage (FAUNMAP + PBDB, from the Step 6 matrices).
# 2. Column graph of the number of sites per stage in each province.
# 3. Simpson dissimilarity (beta_SIM), formulas of Rowan et al. (2024, Nat.
#    Ecol. Evol.) after Baselga (2010):
#
#      pairwise:   beta_SIM = min(b, c) / (a + min(b, c))
#      multisite:  beta_SIM = sum_{i<j} min(b_ij, b_ji) /
#                  [ sum_i S_i - S_T + sum_{i<j} min(b_ij, b_ji) ]     (Eq. 1)
#
#    0 = same species (or one set nested in the other); 1 = no species shared.
#
#    A. HOMOGENIZATION (the Rowan et al. 2024 approach)
#       Rowan et al. computed ONE multisite beta_SIM per time bin, across the
#       spatial units (subregions) of that bin, and then read the TREND of
#       those values through time. Here: within each province and each stage,
#       multisite beta_SIM across the sites (or US physiographic sections, see
#       'spatial_unit') of that province. One value per province per stage.
#         falling values through time = sites becoming more alike (homogenization)
#         rising values through time  = sites becoming more different
#       Multisite beta_SIM rises with the number of units compared, so each
#       value is also computed on random subsets of 'n_units_sample' units
#       (the plotted value; bars = 95% range of the subsets). This makes stages
#       with many and few sites comparable. The full value (all units) and the
#       mean pairwise value are saved in the table as well.
#
#    B. TURNOVER OF THE POOL BETWEEN SUCCESSIVE STAGES
#       Pairwise beta_SIM between the pool of one stage and the pool of the
#       NEXT stage only (Zanclean-Piacenzian, Piacenzian-Gelasian, ...).
#       Non-adjacent stages are not compared. If a stage has no species in a
#       province, the two transitions touching it cannot be computed and the
#       line is broken there (not drawn across the gap).
#
# Inputs:
#   Outputs/6_matrices/all_records_final.csv              (Step 6)
#   Outputs/maps/physiographic/site_physio_database.csv   (Step 5)
#
# Outputs (Outputs/7_regional_pools/):
#   Figures (PNG 300 dpi + PDF, 16:9 slide format):
#     sites_per_stage_<province>.png/.pdf       one column graph per province
#     sites_per_stage_all_provinces.png/.pdf    the four together
#     beta_spatial_all_provinces.png/.pdf       A: within-stage beta_SIM through time
#     beta_spatial_<province>.png/.pdf          A: one province per slide
#     beta_consecutive_all_provinces.png/.pdf   B: successive-stage turnover
#   Tables:
#     regional_pools.xlsx            one sheet per province: species x stage (1/0)
#     regional_pools_long.csv        province, stage, species, order, n sites
#     pool_summary.csv               sites and species per province and stage
#     beta_spatial.csv               A: one row per province and stage
#     beta_spatial_trend.csv         A: change between successive stages (up/down)
#     beta_consecutive.csv           B: beta_SIM between successive stage pools
#     site_province_conflicts.csv    sites whose analysis units fall in >1 province
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Steps 5 and 6.
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

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

records_file <- file.path("Outputs", "6_matrices", "all_records_final.csv")
physio_file  <- file.path("Outputs", "maps", "physiographic", "site_physio_database.csv")
output_dir   <- file.path(work_dir, "Outputs", "7_regional_pools")

provinces <- c("Pacific Border", "Basin and Range", "Great Plains", "Coastal Plain")

# Analysis A (within-stage dissimilarity, Rowan et al. 2024 approach):
spatial_unit   <- "site"   # "site" = each site is a unit; "section" = US physiographic
                           # sections (sites pooled per section, like Rowan's subregions)
n_units_sample <- 3        # units drawn per subset (Rowan et al. compared 3 subregions)
n_resamples    <- 999      # random subsets per province and stage
set.seed(2024)             # same subsets every run

# One colour per province (validated colour-blind-safe set).
province_colours <- c(
  "Pacific Border"  = "#2a78d6",   # blue
  "Basin and Range" = "#eb6834",   # orange
  "Great Plains"    = "#1baf7a",   # aqua
  "Coastal Plain"   = "#4a3aa7"    # violet
)

stage_names <- c("Zanclean", "Piacenzian", "Gelasian", "Calabrian", "Chibanian")
stage_older <- c(4.700, 3.600, 2.580, 1.800, 0.7741)
stage_young <- c(3.600, 2.580, 1.800, 0.7741, 0.129)

# Figure look (slide-ready).
base_size  <- 18
slide_w    <- 13.33   # inches, 16:9
slide_h    <- 7.5
fig_dpi    <- 300
ink        <- "#1f1f1f"
ink_soft   <- "#5a5a5a"
grid_col   <- "#e6e6e3"
no_data    <- "#ececea"
heat_ramp  <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================

read_input <- function(file) {
  path <- file.path(work_dir, file)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun Steps 5 and 6 first.")
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

# Pairwise Simpson dissimilarity between two species sets.
beta_sim_pair <- function(A, B) {
  if (length(A) == 0 || length(B) == 0) return(c(a = NA, b = NA, c = NA, beta = NA))
  a <- length(intersect(A, B)); b <- length(setdiff(A, B)); c <- length(setdiff(B, A))
  c(a = a, b = b, c = c, beta = min(b, c) / (a + min(b, c)))
}

# Multisite Simpson dissimilarity (Rowan et al. 2024 Eqs. 1-3).
beta_sim_multi <- function(pools) {
  pools <- pools[lengths(pools) > 0]
  k <- length(pools)
  if (k < 2) return(c(n_units = k, beta_SIM = NA, beta_SIM_END = NA, beta_SIM_SH = NA))
  S_i <- lengths(pools)
  S_T <- length(unique(unlist(pools)))
  min_sum <- 0; a_sum <- 0
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    min_sum <- min_sum + min(length(setdiff(pools[[i]], pools[[j]])),
                             length(setdiff(pools[[j]], pools[[i]])))
    a_sum <- a_sum + length(intersect(pools[[i]], pools[[j]]))
  }
  shared <- sum(S_i) - S_T
  c(n_units = k,
    beta_SIM     = min_sum / (shared + min_sum),
    beta_SIM_END = min_sum / (min_sum + a_sum),
    beta_SIM_SH  = shared / sum(S_i))
}

# Mean of all pairwise beta_SIM values among a list of species sets.
beta_sim_mean_pair <- function(pools) {
  pools <- pools[lengths(pools) > 0]
  k <- length(pools)
  if (k < 2) return(NA_real_)
  v <- c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) v <- c(v, beta_sim_pair(pools[[i]], pools[[j]])[["beta"]])
  mean(v)
}

# Multisite beta_SIM on subsets of n units: every subset if there are few,
# otherwise n_resamples random subsets. Returns mean and 95% range.
beta_sim_resampled <- function(pools, n) {
  pools <- pools[lengths(pools) > 0]
  k <- length(pools)
  if (k < n || n < 2) return(c(mean = NA, lo = NA, hi = NA, n_subsets = 0))
  subsets <- if (choose(k, n) <= n_resamples) combn(k, n, simplify = FALSE) else
    replicate(n_resamples, sample.int(k, n), simplify = FALSE)
  v <- vapply(subsets, function(s) beta_sim_multi(pools[s])[["beta_SIM"]], numeric(1))
  c(mean = mean(v), lo = unname(quantile(v, 0.025)), hi = unname(quantile(v, 0.975)),
    n_subsets = length(v))
}

# Shared slide theme.
theme_slide <- function() {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = ink),
      plot.title = element_text(face = "bold", size = base_size * 1.25, margin = margin(b = 4)),
      plot.subtitle = element_text(colour = ink_soft, size = base_size * 0.85, margin = margin(b = 12)),
      plot.caption = element_text(colour = ink_soft, size = base_size * 0.65, hjust = 0),
      axis.title = element_text(colour = ink_soft, size = base_size * 0.9),
      axis.text = element_text(colour = ink, size = base_size * 0.8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size * 0.95, hjust = 0),
      plot.margin = margin(18, 24, 14, 18),
      legend.position = "none"
    )
}

beta_lab <- quote(beta[SIM])

# =============================================================================
# 1. READ AND LINK SITES TO PROVINCES
# =============================================================================

cat("=== 1. READING FILES ===\n")
records <- read_input(records_file)
physio  <- read_input(physio_file)

# Province of each site (per stage and database). A FAUNMAP site name can
# cover several analysis units; if they fall in different provinces, the
# province of most units is used and the site is listed in the conflicts file.
site_province_all <- physio %>%
  transmute(Stage_Number = as.integer(Stage_Number), Database,
            SiteName = trimws(SiteName), Province = trimws(US_Province),
            Section = if ("US_Section" %in% names(physio)) trimws(US_Section) else NA_character_) %>%
  filter(!is.na(SiteName), !is.na(Province))

site_province <- site_province_all %>%
  count(Stage_Number, Database, SiteName, Province) %>%
  group_by(Stage_Number, Database, SiteName) %>%
  mutate(n_provinces = n()) %>%
  arrange(desc(n), Province, .by_group = TRUE) %>%
  ungroup()

conflicts <- site_province %>% filter(n_provinces > 1)
site_province <- site_province %>%
  distinct(Stage_Number, Database, SiteName, .keep_all = TRUE) %>%
  select(Stage_Number, Database, SiteName, Province)

# Physiographic section of each site (majority within its province).
site_section <- site_province_all %>%
  filter(!is.na(Section)) %>%
  count(Stage_Number, Database, SiteName, Province, Section) %>%
  group_by(Stage_Number, Database, SiteName, Province) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-n)
site_province <- left_join(site_province, site_section,
                           by = c("Stage_Number", "Database", "SiteName", "Province"))
save_csv(conflicts, "site_province_conflicts")

# Match the requested provinces regardless of capitalisation.
available <- sort(unique(site_province$Province))
prov_match <- setNames(available[match(tolower(provinces), tolower(available))], provinces)
cat("\n  Provinces found in the site database:\n")
for (p in provinces) {
  cat(sprintf("    %-16s %s\n", p, if (is.na(prov_match[[p]])) "NOT FOUND" else "found"))
}
if (any(is.na(prov_match))) {
  cat("  Available US provinces:", paste(available, collapse = ", "), "\n")
}

recs <- records %>%
  transmute(Stage_Number = as.integer(Stage_Number), Database,
            SiteName = trimws(SiteName), GenusSpecies, Order) %>%
  filter(!is.na(GenusSpecies)) %>%
  inner_join(site_province, by = c("Stage_Number", "Database", "SiteName")) %>%
  mutate(Province = provinces[match(tolower(Province), tolower(provinces))]) %>%
  filter(!is.na(Province), Stage_Number %in% 1:5) %>%
  mutate(Stage = factor(stage_names[Stage_Number], levels = stage_names),
         Province = factor(Province, levels = provinces))

cat(sprintf("\n  Species records in the four provinces: %d (from %d sites)\n",
            nrow(recs), n_distinct(paste(recs$Database, recs$SiteName))))

# =============================================================================
# 2. REGIONAL SPECIES POOLS
# =============================================================================

cat("\n=== 2. REGIONAL SPECIES POOLS ===\n")

pools_long <- recs %>%
  group_by(Province, Stage, Stage_Number, GenusSpecies, Order) %>%
  summarise(n_sites_with_species = n_distinct(paste(Database, SiteName)), .groups = "drop") %>%
  arrange(Province, Stage_Number, Order, GenusSpecies)

# Sites per province and stage: sites contributing species to the pool, and
# all sites mapped in the province (some have no retained species).
sites_total <- site_province %>%
  mutate(Province = provinces[match(tolower(Province), tolower(provinces))]) %>%
  filter(!is.na(Province), Stage_Number %in% 1:5) %>%
  count(Province, Stage_Number, name = "n_sites_mapped")

pool_summary <- expand_grid(Province = factor(provinces, levels = provinces),
                            Stage_Number = 1:5) %>%
  left_join(recs %>%
              group_by(Province, Stage_Number) %>%
              summarise(n_sites = n_distinct(paste(Database, SiteName)),
                        n_FAUNMAP_sites = n_distinct(SiteName[Database == "FAUNMAP"]),
                        n_PBDB_sites = n_distinct(SiteName[Database == "PBDB"]),
                        n_species = n_distinct(GenusSpecies), .groups = "drop"),
            by = c("Province", "Stage_Number")) %>%
  left_join(mutate(sites_total, Province = factor(Province, levels = provinces)),
            by = c("Province", "Stage_Number")) %>%
  mutate(across(c(n_sites, n_FAUNMAP_sites, n_PBDB_sites, n_species, n_sites_mapped),
                ~ coalesce(.x, 0L)),
         Stage = factor(stage_names[Stage_Number], levels = stage_names)) %>%
  relocate(Stage, .after = Stage_Number)

print(as.data.frame(pool_summary %>% select(Province, Stage, n_sites, n_species)), row.names = FALSE)

save_csv(pools_long, "regional_pools_long")
save_csv(pool_summary, "pool_summary")

# One sheet per province: species x stage presence (1/0) and the number of sites.
pool_sheets <- lapply(provinces, function(p) {
  d <- filter(pools_long, Province == p)
  if (nrow(d) == 0) return(data.frame(Note = "No species recorded in this province"))
  d %>%
    mutate(Present = 1L) %>%
    select(Order, GenusSpecies, Stage, Present) %>%
    pivot_wider(names_from = Stage, values_from = Present, values_fill = 0L) %>%
    { for (s in stage_names) if (!s %in% names(.)) .[[s]] <- 0L; . } %>%
    select(Order, Species = GenusSpecies, all_of(stage_names)) %>%
    mutate(n_stages = rowSums(across(all_of(stage_names)))) %>%
    arrange(Order, Species)
})
names(pool_sheets) <- substr(provinces, 1, 31)
save_with_fallback(function(f) writexl::write_xlsx(pool_sheets, f),
                   file.path(output_dir, "regional_pools.xlsx"))

pool_sets <- lapply(provinces, function(p) {
  lapply(1:5, function(s) {
    unique(pools_long$GenusSpecies[pools_long$Province == p & pools_long$Stage_Number == s])
  })
})
names(pool_sets) <- provinces

# =============================================================================
# 3A. WITHIN-STAGE DISSIMILARITY ACROSS SITES (Rowan et al. 2024 approach)
# =============================================================================

unit_word <- if (spatial_unit == "section") "sections" else "sites"
cat(sprintf("\n=== 3A. WITHIN-STAGE beta_SIM ACROSS %s (Rowan et al. 2024 approach) ===\n",
            toupper(unit_word)))

recs_units <- recs %>%
  mutate(Unit = if (spatial_unit == "section") Section else paste(Database, SiteName, sep = " | "))
if (spatial_unit == "section") {
  n_no_section <- n_distinct(paste(recs_units$Database, recs_units$SiteName)[is.na(recs_units$Unit)])
  if (n_no_section > 0) cat(sprintf("  NOTE: %d sites have no US section and are left out of 3A\n", n_no_section))
}
recs_units <- filter(recs_units, !is.na(Unit))

beta_spatial <- bind_rows(lapply(provinces, function(p) bind_rows(lapply(1:5, function(s) {
  d <- filter(recs_units, Province == p, Stage_Number == s)
  units <- lapply(split(d$GenusSpecies, d$Unit), unique)
  full <- beta_sim_multi(units)
  rs <- beta_sim_resampled(units, n_units_sample)
  data.frame(Province = p, Stage_Number = s, Stage = stage_names[s],
             Mid_Ma = (stage_older[s] + stage_young[s]) / 2,
             n_units = length(units), n_species = n_distinct(d$GenusSpecies),
             beta_SIM = rs[["mean"]], beta_SIM_lo95 = rs[["lo"]], beta_SIM_hi95 = rs[["hi"]],
             n_subsets = rs[["n_subsets"]],
             beta_SIM_all_units = full[["beta_SIM"]],
             beta_SIM_END_all_units = full[["beta_SIM_END"]],
             beta_SIM_SH_all_units = full[["beta_SIM_SH"]],
             mean_pairwise_beta_SIM = beta_sim_mean_pair(units))
}))))
names(beta_spatial)[names(beta_spatial) == "n_units"] <- paste0("n_", unit_word)

# Change between successive stages that have a value. A change is called
# "clear" only when the 95% ranges of the two stages do not overlap.
beta_spatial_trend <- beta_spatial %>%
  filter(!is.na(beta_SIM)) %>%
  group_by(Province) %>%
  arrange(Stage_Number, .by_group = TRUE) %>%
  mutate(From = lag(Stage), From_n = lag(Stage_Number),
         b_from = lag(beta_SIM), lo_from = lag(beta_SIM_lo95), hi_from = lag(beta_SIM_hi95)) %>%
  ungroup() %>%
  filter(!is.na(From)) %>%
  transmute(
    Province = factor(Province, levels = provinces),
    Comparison = paste(From, "to", Stage),
    Skipped_stage = ifelse(Stage_Number - From_n > 1,
                           "yes - stage(s) in between had too few sites", ""),
    beta_SIM_from = round(b_from, 3), beta_SIM_to = round(beta_SIM, 3),
    Change = round(beta_SIM - b_from, 3),
    Direction = case_when(
      beta_SIM_lo95 > hi_from ~ "increase: sites MORE different",
      beta_SIM_hi95 < lo_from ~ "decrease: sites MORE alike (homogenization)",
      TRUE ~ ifelse(beta_SIM > b_from, "slight increase (95% ranges overlap)",
                    "slight decrease (95% ranges overlap)"))
  ) %>%
  arrange(Province)

save_csv(beta_spatial, "beta_spatial")
save_csv(beta_spatial_trend, "beta_spatial_trend")

cat(sprintf("  Multisite beta_SIM among %s within each stage (mean of subsets of %d %s):\n",
            unit_word, n_units_sample, unit_word))
print(as.data.frame(beta_spatial %>%
        transmute(Province, Stage, n = .data[[paste0("n_", unit_word)]],
                  beta_SIM = round(beta_SIM, 3), lo95 = round(beta_SIM_lo95, 3),
                  hi95 = round(beta_SIM_hi95, 3))), row.names = FALSE)
cat("\n  Trend (lower beta_SIM = sites more alike = homogenization):\n")
print(as.data.frame(select(beta_spatial_trend, Province, Comparison, beta_SIM_from,
                           beta_SIM_to, Direction)), row.names = FALSE)

# =============================================================================
# 3B. TURNOVER OF EACH POOL BETWEEN SUCCESSIVE STAGES
# =============================================================================

cat("\n=== 3B. POOL TURNOVER BETWEEN SUCCESSIVE STAGES ===\n")

beta_consecutive <- bind_rows(lapply(provinces, function(p) bind_rows(lapply(1:4, function(i) {
  A <- pool_sets[[p]][[i]]; B <- pool_sets[[p]][[i + 1]]
  r <- beta_sim_pair(A, B)
  empty <- stage_names[c(i, i + 1)][c(length(A), length(B)) == 0]
  data.frame(Province = p, Transition = paste(stage_names[i], "to", stage_names[i + 1]),
             Stage_A = stage_names[i], Stage_B = stage_names[i + 1],
             Stage_A_n = i, Boundary_Ma = stage_young[i],
             n_species_A = length(A), n_species_B = length(B),
             shared_a = r[["a"]], only_A_b = r[["b"]], only_B_c = r[["c"]],
             beta_SIM = r[["beta"]],
             Note = if (length(empty)) paste("no species in", paste(empty, collapse = " and ")) else "")
}))))

save_csv(beta_consecutive, "beta_consecutive")
print(as.data.frame(beta_consecutive %>% mutate(beta_SIM = round(beta_SIM, 3)) %>%
        select(Province, Transition, shared_a, beta_SIM, Note)), row.names = FALSE)

# =============================================================================
# 4. FIGURES
# =============================================================================

cat("\n=== 4. FIGURES ===\n")

stage_axis <- function() {
  paste0(stage_names, "\n", formatC(stage_older, format = "g"), "-",
         formatC(stage_young, format = "g"), " Ma")
}

# Alternating background bands, one per stage (x axis in Ma, oldest left).
bands <- data.frame(xmin = stage_young, xmax = stage_older, Stage = stage_names,
                    fill = rep(c("#f3f3f0", "#ffffff"), length.out = 5))

# Line segments only between neighbouring points that both have a value, so
# a missing stage leaves a visible gap instead of a line across it.
gap_segments <- function(d, x, order_col) {
  d %>%
    arrange(Province, .data[[order_col]]) %>%
    group_by(Province) %>%
    mutate(x2 = lead(.data[[x]]), y2 = lead(beta_SIM)) %>%
    ungroup() %>%
    filter(!is.na(beta_SIM), !is.na(y2))
}

# ---- 4a. Sites per stage: one column graph per province ----------------------
site_plot <- function(d, title, colour) {
  ymax <- max(1, d$n_sites) * 1.15
  ggplot(d, aes(x = Stage, y = n_sites)) +
    geom_col(fill = colour, width = 0.68) +
    geom_text(aes(label = n_sites), vjust = -0.45, size = base_size / 3.2,
              fontface = "bold", colour = ink) +
    scale_x_discrete(labels = setNames(stage_axis(), stage_names), drop = FALSE) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.02)),
                       breaks = scales::pretty_breaks(5)) +
    labs(title = title, x = NULL, y = "Number of sites",
         subtitle = "Sites contributing species to the regional pool (FAUNMAP + PBDB), oldest to youngest") +
    theme_slide()
}

for (p in provinces) {
  d <- filter(pool_summary, Province == p)
  save_fig(site_plot(d, paste(p, "- sites per geological stage"), province_colours[[p]]),
           paste0("sites_per_stage_", file_stub(p)))
}

p_all_sites <- ggplot(pool_summary, aes(x = Stage, y = n_sites, fill = Province)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n_sites), vjust = -0.4, size = base_size / 3.8,
            fontface = "bold", colour = ink) +
  facet_wrap(~ Province, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = province_colours) +
  scale_x_discrete(labels = substr(stage_names, 1, 4), drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), breaks = scales::pretty_breaks(4)) +
  labs(title = "Sites per geological stage in four physiographic provinces",
       subtitle = "Zanclean (4.7 Ma) to Chibanian (0.129 Ma); sites contributing species to each regional pool",
       x = NULL, y = "Number of sites") +
  theme_slide()
save_fig(p_all_sites, "sites_per_stage_all_provinces")
cat("  Site-count column graphs saved\n")

# ---- 4b. A: within-stage beta_SIM through time -------------------------------
spatial_caption <- paste0(
  "Each point: multisite Simpson dissimilarity among ", unit_word, " of the province in that stage ",
  "(mean of subsets of ", n_units_sample, " ", unit_word, "; bars = 95% range).\n",
  "Lower = ", unit_word, " share more species (more homogeneous). n = ", unit_word,
  " with species; no point = fewer than ", n_units_sample, ". Approach of Rowan et al. (2024).")

spatial_plot <- function(d, title, facet = FALSE) {
  n_col <- paste0("n_", unit_word)
  d <- d %>% mutate(Province = factor(Province, levels = provinces),
                    n_lab = paste0("n=", .data[[n_col]]),
                    n_col_ink = ifelse(is.na(beta_SIM), "#b0b0ab", ink_soft),
                    v_lab = ifelse(is.na(beta_SIM), "", sprintf("%.2f", beta_SIM)))
  sz <- if (facet) 0.8 else 1
  stage_text <- if (facet) substr(stage_names, 1, 4) else stage_names
  ggplot(d) +
    geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill)) +
    scale_fill_identity() +
    geom_hline(yintercept = seq(0, 1, 0.25), colour = grid_col, linewidth = 0.4) +
    geom_text(data = transform(bands, lab = stage_text),
              aes(x = (xmin + xmax) / 2, y = 1.1, label = lab),
              size = base_size / 3.9 * sz, colour = ink_soft) +
    geom_errorbar(aes(x = Mid_Ma, ymin = beta_SIM_lo95, ymax = beta_SIM_hi95, colour = Province),
                  width = 0.14, linewidth = 0.9, na.rm = TRUE) +
    geom_segment(data = gap_segments(d, "Mid_Ma", "Stage_Number"),
                 aes(x = Mid_Ma, xend = x2, y = beta_SIM, yend = y2, colour = Province),
                 linewidth = 1.3) +
    geom_point(aes(Mid_Ma, beta_SIM, colour = Province), size = 4.5 * sz, shape = 21,
               fill = "white", stroke = 1.8, na.rm = TRUE) +
    geom_text(aes(x = Mid_Ma - 0.12, y = beta_SIM, label = v_lab, colour = Province),
              hjust = 0, vjust = -0.6, size = base_size / 3.4 * sz, fontface = "bold", na.rm = TRUE) +
    geom_text(aes(x = Mid_Ma, y = -0.07, label = n_lab, colour = n_col_ink),
              size = base_size / 4.2 * sz) +
    scale_colour_manual(values = c(province_colours, setNames(c(ink_soft, "#b0b0ab"),
                                                              c(ink_soft, "#b0b0ab")))) +
    scale_x_reverse(limits = c(4.7, -0.15), breaks = c(4, 3, 2, 1, 0), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(-0.12, 1.16), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(title = title, x = "Age (Ma)", y = beta_lab, caption = spatial_caption) +
    theme_slide() +
    theme(plot.title.position = "plot", plot.caption.position = "plot") +
    if (facet) facet_wrap(~ Province, ncol = 2) else NULL
}

p_spatial <- spatial_plot(beta_spatial, "Are sites within each province becoming more alike?", facet = TRUE) +
  labs(subtitle = sprintf("Simpson dissimilarity among %s within each stage; falling line = homogenization",
                          unit_word))
save_fig(p_spatial, "beta_spatial_all_provinces", width = slide_w, height = slide_h * 1.25)

for (p in provinces) {
  d <- filter(beta_spatial, Province == p)
  pp <- spatial_plot(d, paste(p, "- dissimilarity among", unit_word, "through time")) +
    labs(subtitle = "Multisite Simpson dissimilarity within each stage; falling line = homogenization")
  save_fig(pp, paste0("beta_spatial_", file_stub(p)))
}
cat("  Within-stage beta_SIM figures saved\n")

# ---- 4c. B: successive-stage turnover through time ---------------------------
cons <- beta_consecutive %>% mutate(Province = factor(Province, levels = provinces))
label_pts <- cons %>% filter(!is.na(beta_SIM)) %>%
  group_by(Province) %>% slice_min(Boundary_Ma, n = 1) %>% ungroup()

# Spread end labels vertically so they never overlap (minimum gap in y units).
spread <- function(y, gap = 0.07) {
  if (length(y) == 0) return(y)
  o <- order(y); z <- y[o]
  for (k in seq_along(z)[-1]) z[k] <- max(z[k], z[k - 1] + gap)
  over <- max(0, z[length(z)] - 1.02); z <- z - over
  y[o] <- z; y
}
label_pts$label_y <- spread(label_pts$beta_SIM)
label_x <- 0.05

missing <- filter(cons, is.na(beta_SIM))
cons_caption <- paste0(
  "Each point compares a province's species pool with its pool in the NEXT stage; it sits on the boundary between them.\n",
  if (nrow(missing)) paste0("Gaps (no comparison possible): ",
                           paste(sprintf("%s %s", missing$Province, missing$Transition), collapse = "; "),
                           ".\n") else "",
  "Pairwise Simpson dissimilarity after Rowan et al. (2024).")

p_cons <- ggplot() +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill)) +
  scale_fill_identity() +
  geom_hline(yintercept = seq(0, 1, 0.25), colour = grid_col, linewidth = 0.4) +
  geom_text(data = bands, aes(x = (xmin + xmax) / 2, y = 1.07, label = Stage),
            size = base_size / 3.9, colour = ink_soft) +
  geom_segment(data = gap_segments(cons, "Boundary_Ma", "Stage_A_n"),
               aes(x = Boundary_Ma, xend = x2, y = beta_SIM, yend = y2, colour = Province),
               linewidth = 1.3) +
  geom_point(data = cons, aes(Boundary_Ma, beta_SIM, colour = Province),
             size = 4.2, shape = 21, fill = "white", stroke = 1.8, na.rm = TRUE) +
  geom_segment(data = label_pts, aes(x = Boundary_Ma - 0.06, xend = label_x + 0.03,
                                     y = beta_SIM, yend = label_y, colour = Province),
               linewidth = 0.5) +
  geom_text(data = label_pts, aes(x = label_x, y = label_y, label = Province, colour = Province),
            hjust = 0, size = base_size / 3.4, fontface = "bold") +
  scale_colour_manual(values = province_colours) +
  scale_x_reverse(limits = c(4.7, -1.45), breaks = c(4, 3, 2, 1, 0), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "Turnover of each regional pool between successive stages",
       subtitle = "Pairwise Simpson dissimilarity at each stage boundary (higher = more species replaced)",
       x = "Age (Ma)", y = beta_lab, caption = cons_caption) +
  theme_slide() +
  theme(plot.title.position = "plot", plot.caption.position = "plot")
save_fig(p_cons, "beta_consecutive_all_provinces", width = slide_w, height = slide_h * 1.05)
cat("  Successive-stage turnover figure saved\n")

cat("\n=== STEP 7 COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: pools_long, pool_summary, beta_spatial,",
    "beta_spatial_trend, beta_consecutive\n")
