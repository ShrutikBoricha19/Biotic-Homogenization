# =============================================================================
# STEP 9a: ARE THE FOUR PROVINCES BECOMING MORE ALIKE?
#          (between-province beta_SIM - the closest analogue of Rowan et al.)
#
# Rowan et al. (2024) compared three SUBREGIONS (Afar, Turkana, Southern Rift)
# in each time bin. Here the four physiographic provinces play that role:
# for each stage, the regional pools of the provinces are compared with each
# other.
#
#   a) Multisite beta_SIM across all provinces with species in that stage
#      (Rowan et al. Eq. 1), with its endemic and shared components.
#   b) The same with EQUAL SAMPLING: in each draw every province contributes
#      the same number of randomly chosen sites (the smallest number available
#      in that stage, at least 'min_sites_equal'), so a province with many
#      sites does not get a richer pool just because it was sampled more.
#      Mean and 95% range of 'n_resamples' draws.
#   c) Pairwise beta_SIM between every pair of provinces, stage by stage:
#      which provinces converged?
#
# Falling values through time = the provinces' faunas became more similar
# (homogenization at the scale of western/central North America).
#
# Inputs:  all_records_final.csv (Step 6), site_physio_database.csv (Step 5)
# Outputs (Outputs/9a_between_provinces/):
#   between_provinces_trend.png/.pdf    a) and b) through time
#   between_provinces_pairs.png/.pdf    c) one panel per pair of provinces
#   between_provinces.csv, between_provinces_pairs.csv
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Steps 5-7.
# Needs: dplyr, tidyr, ggplot2.
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

records_file <- file.path("Outputs", "6_matrices", "all_records_final.csv")
physio_file  <- file.path("Outputs", "maps", "physiographic", "site_physio_database.csv")
output_dir   <- file.path(work_dir, "Outputs", "9a_between_provinces")

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

# Equal sampling (b): every province contributes the same number of sites.
min_sites_equal <- 3      # provinces with fewer sites in a stage are left out of (b)
n_resamples     <- 999
set.seed(2024)

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

# =============================================================================
# 2. BETWEEN-PROVINCE DISSIMILARITY PER STAGE
# =============================================================================

cat("\n=== 2. BETWEEN-PROVINCE beta_SIM ===\n")

trend <- bind_rows(lapply(1:5, function(s) {
  d <- filter(recs, Stage_Number == s)
  pools <- lapply(setNames(provinces, provinces), function(p) unique(d$GenusSpecies[d$Province == p]))
  pools <- pools[lengths(pools) > 0]
  full <- beta_multi_all(pools)

  # (b) equal number of sites per province
  sites_per <- sapply(names(pools), function(p) n_distinct(d$Site[d$Province == p]))
  use <- names(sites_per)[sites_per >= min_sites_equal]
  eq <- c(mean = NA, lo = NA, hi = NA); n_eq <- NA
  if (length(use) >= 2) {
    n_eq <- min(sites_per[use])
    by_prov <- lapply(use, function(p) site_lists(filter(d, Province == p)))
    draws <- replicate(n_resamples, {
      pp <- lapply(by_prov, function(sl) unique(unlist(sl[sample.int(length(sl), n_eq)])))
      beta_multi_all(pp)[["beta_SIM"]]
    })
    eq <- c(mean = mean(draws), lo = unname(quantile(draws, 0.025)), hi = unname(quantile(draws, 0.975)))
  }
  data.frame(Stage_Number = s, Stage = stage_names[s], Mid_Ma = stage_mid[s],
             Provinces_with_species = paste(names(pools), collapse = "; "),
             n_provinces = length(pools),
             beta_SIM = full[["beta_SIM"]], beta_SIM_END = full[["beta_SIM_END"]],
             beta_SIM_SH = full[["beta_SIM_SH"]],
             Provinces_equal_sampling = paste(use, collapse = "; "),
             n_provinces_equal = length(use), sites_per_province_equal = n_eq,
             beta_SIM_equal = eq[["mean"]], beta_SIM_equal_lo95 = eq[["lo"]],
             beta_SIM_equal_hi95 = eq[["hi"]])
}))

pairs <- bind_rows(lapply(1:5, function(s) {
  d <- filter(recs, Stage_Number == s)
  cmb <- combn(provinces, 2)
  bind_rows(lapply(seq_len(ncol(cmb)), function(k) {
    A <- unique(d$GenusSpecies[d$Province == cmb[1, k]])
    B <- unique(d$GenusSpecies[d$Province == cmb[2, k]])
    data.frame(Stage_Number = s, Stage = stage_names[s], Mid_Ma = stage_mid[s],
               Province_A = cmb[1, k], Province_B = cmb[2, k],
               species_A = length(A), species_B = length(B),
               shared = length(intersect(A, B)), beta_SIM = beta_sim_pair(A, B))
  }))
})) %>%
  mutate(Pair = factor(paste(Province_A, "-", Province_B),
                       levels = unique(paste(Province_A, "-", Province_B))))

save_csv(mutate(trend, across(where(is.numeric), ~ round(.x, 4))), "between_provinces")
save_csv(mutate(pairs, beta_SIM = round(beta_SIM, 4)), "between_provinces_pairs")

print(as.data.frame(trend %>% transmute(Stage, n_provinces, beta_SIM = round(beta_SIM, 3),
                                        equal = round(beta_SIM_equal, 3),
                                        lo95 = round(beta_SIM_equal_lo95, 3),
                                        hi95 = round(beta_SIM_equal_hi95, 3),
                                        sites_each = sites_per_province_equal)), row.names = FALSE)

# =============================================================================
# 3. FIGURES
# =============================================================================

cat("\n=== 3. FIGURES ===\n")

series <- bind_rows(
  transmute(trend, Series = "All sites", Stage_Number, Mid_Ma, value = beta_SIM,
            lo = NA_real_, hi = NA_real_, n = n_provinces),
  transmute(trend, Series = "Equal sites per province", Stage_Number, Mid_Ma,
            value = beta_SIM_equal, lo = beta_SIM_equal_lo95, hi = beta_SIM_equal_hi95,
            n = n_provinces_equal)
) %>% mutate(Series = factor(Series, levels = c("All sites", "Equal sites per province")))
series_cols <- c("All sites" = ink, "Equal sites per province" = "#2a78d6")


p_trend <- ggplot(series) +
  stage_bands(y_lab = 1.07) +
  geom_ribbon(data = filter(series, !is.na(lo)),
              aes(x = Mid_Ma, ymin = lo, ymax = hi, group = Series), fill = "#2a78d6",
              alpha = 0.15, na.rm = TRUE) +
  geom_segment(data = gap_segments(series, "Series"),
               aes(x = Mid_Ma, xend = x2, y = value, yend = y2, colour = Series), linewidth = 1.3) +
  geom_point(aes(Mid_Ma, value, colour = Series), shape = 21, fill = "white", size = 4.2,
             stroke = 1.8, na.rm = TRUE) +
  geom_text(data = filter(series, Series == "All sites", !is.na(value)),
            aes(Mid_Ma, -0.05, label = paste0(n, " prov.")), size = base_size / 4.6,
            colour = ink_soft) +
  scale_colour_manual(values = series_cols, name = NULL) +
  age_axis() +
  scale_y_continuous(limits = c(-0.1, 1.12), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "Are the four provinces becoming more alike?",
       subtitle = "Multisite Simpson dissimilarity among the regional pools of the provinces; falling line = homogenization",
       x = "Age (Ma)", y = beta_lab,
       caption = paste0(
         "All sites: every site of each province pooled. Equal sites: each province contributes the same number ",
         "of random sites (at least ", min_sites_equal, "),\nmean and 95% range of ", n_resamples,
         " draws. 'prov.' = provinces with species in the stage. After Rowan et al. (2024, Eq. 1).")) +
  theme_slide() +
  theme(legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = base_size * 0.8))
save_fig(p_trend, "between_provinces_trend")

p_pairs <- ggplot(pairs) +
  stage_bands() +
  geom_segment(data = gap_segments(rename(pairs, value = beta_SIM), "Pair"),
               aes(x = Mid_Ma, xend = x2, y = value, yend = y2), colour = ink, linewidth = 1.1) +
  geom_point(aes(Mid_Ma, beta_SIM), shape = 21, fill = "white", colour = ink, size = 3.2,
             stroke = 1.5, na.rm = TRUE) +
  facet_wrap(~ Pair, ncol = 3) +
  age_axis() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "Which provinces converged?",
       subtitle = "Pairwise Simpson dissimilarity between two provinces' regional pools, stage by stage (falling = more alike)",
       x = "Age (Ma)", y = beta_lab,
       caption = "No point = one of the two provinces has no species in that stage. Bands = stages, Zanclean (left) to Chibanian (right).") +
  theme_slide() +
  theme(panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5),
        strip.text = element_text(face = "bold", size = base_size * 0.75, hjust = 0))
save_fig(p_pairs, "between_provinces_pairs", height = slide_h * 1.05)
cat("  Figures saved\n")

cat("\n=== STEP 9a COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: trend, pairs\n")
