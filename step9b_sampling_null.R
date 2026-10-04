# =============================================================================
# STEP 9b: IS THE TREND REAL OR A SAMPLING ARTEFACT?
#
#   a) Sampling correlations (Rowan et al. 2024, Extended Data Fig. 4): the
#      within-province beta_SIM of Step 7 plotted against the number of sites,
#      the number of species in the pool and the mean species per site. If
#      beta_SIM simply follows these, the "trend" may only reflect sampling.
#   b) Null model: in each province and stage, species are reshuffled among the
#      sites while every site keeps its number of species and every species
#      keeps its number of sites ("fixed-fixed" curveball algorithm). The real
#      beta_SIM is compared with 'n_null' shuffled versions:
#        SES = (observed - mean of nulls) / sd of nulls
#        SES > +1.96  sites more distinct than chance (real spatial structure)
#        SES near 0   sites look like random draws from one shared pool
#        SES < -1.96  sites more alike than chance
#      SES falling towards (or below) 0 through time = the spatial structure of
#      the fauna dissolved = homogenization that sampling cannot explain.
#   c) Trend test: linear regression of Step 7 beta_SIM on stage midpoint age
#      for each province (as Rowan et al. did; with 4-5 stages this is weak
#      supporting evidence only).
#
# Inputs:  all_records_final.csv (Step 6), site_physio_database.csv (Step 5),
#          Outputs/7_regional_pools/beta_spatial.csv (Step 7)
# Outputs (Outputs/9b_sampling_null/):
#   sampling_correlations.png/.pdf, null_model_ses.png/.pdf
#   sampling_measures.csv, null_model.csv, trend_test.csv
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 7.
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
output_dir   <- file.path(work_dir, "Outputs", "9b_sampling_null")

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

step7_file <- file.path("Outputs", "7_regional_pools", "beta_spatial.csv")

# Null model (b).
n_units_sample <- 3       # sites per subset, as in Step 7
n_subsets_null <- 200     # site subsets per matrix (fewer than Step 7, for speed)
n_null         <- 199     # shuffled matrices per province and stage
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
# 2. SAMPLING CORRELATIONS (a)
# =============================================================================

cat("\n=== 2. SAMPLING CORRELATIONS ===\n")

s7 <- read_input(step7_file) %>%
  transmute(Province, Stage, beta_SIM = as.numeric(beta_SIM))

sampling <- recs %>%
  group_by(Province, Stage_Number, Stage) %>%
  summarise(n_sites = n_distinct(Site), n_species = n_distinct(GenusSpecies),
            mean_species_per_site = n() / n_distinct(Site), .groups = "drop") %>%
  mutate(Stage = as.character(Stage)) %>%
  left_join(s7, by = c("Province", "Stage")) %>%
  mutate(Province = factor(Province, levels = provinces))
save_csv(mutate(sampling, across(where(is.numeric), ~ round(.x, 4))), "sampling_measures")

measures <- c(n_sites = "Number of sites", n_species = "Species in the regional pool",
              mean_species_per_site = "Mean species per site")
long <- sampling %>%
  filter(!is.na(beta_SIM)) %>%
  pivot_longer(all_of(names(measures)), names_to = "Measure", values_to = "x") %>%
  mutate(Measure = factor(measures[Measure], levels = measures))

fits <- long %>%
  group_by(Measure) %>%
  group_modify(function(d, k) {
    if (nrow(d) < 3 || length(unique(d$x)) < 2) return(data.frame(r2 = NA, p = NA, n = nrow(d)))
    m <- summary(lm(beta_SIM ~ x, data = d))
    data.frame(r2 = m$r.squared, p = coef(m)[2, 4], n = nrow(d))
  }) %>%
  ungroup() %>%
  mutate(label = ifelse(is.na(r2), "too few points",
                        sprintf("R-squared = %.2f, p %s, n = %d", r2,
                                ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)), n)))
print(as.data.frame(select(fits, Measure, r2, p, n)), row.names = FALSE)

p_samp <- ggplot(long, aes(x, beta_SIM)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = ink_soft, fill = "#e2e2de",
              linewidth = 0.9) +
  geom_point(aes(colour = Province), size = 4, alpha = 0.9) +
  geom_text(data = fits, aes(x = -Inf, y = 1.06, label = label), hjust = -0.05,
            size = base_size / 4.2, colour = ink) +
  facet_wrap(~ Measure, scales = "free_x") +
  scale_colour_manual(values = province_colours, name = NULL) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
  labs(title = "Does dissimilarity simply follow sampling?",
       subtitle = "Within-province beta_SIM (Step 7) against sampling measures; one point per province and stage",
       x = NULL, y = beta_lab,
       caption = paste0("A strong relationship means the beta_SIM trend could be a sampling artefact.\n",
                        "Line = linear fit of all points with 95% confidence band (after Rowan et al. 2024, Extended Data Fig. 4).")) +
  theme_slide() +
  theme(legend.position = "top", legend.justification = "left",
        panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5),
        panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.4))
save_fig(p_samp, "sampling_correlations")

# =============================================================================
# 3. NULL MODEL (b)
# =============================================================================

cat("\n=== 3. NULL MODEL (fixed-fixed curveball) ===\n")

# Curveball randomisation (Strona et al. 2014): keeps row and column totals.
curveball <- function(M, n_trades = 5 * nrow(M)) {
  rows <- lapply(seq_len(nrow(M)), function(i) which(M[i, ] == 1L))
  for (t in seq_len(n_trades)) {
    ij <- sample.int(length(rows), 2)
    a <- rows[[ij[1]]]; b <- rows[[ij[2]]]
    ua <- setdiff(a, b); ub <- setdiff(b, a)
    if (length(ua) > 0 && length(ub) > 0) {
      pool <- c(ua, ub)
      perm <- pool[sample.int(length(pool))]
      both <- intersect(a, b)
      rows[[ij[1]]] <- c(both, perm[seq_along(ua)])
      rows[[ij[2]]] <- c(both, perm[-seq_along(ua)])
    }
  }
  out <- matrix(0L, nrow(M), ncol(M))
  for (i in seq_along(rows)) out[i, rows[[i]]] <- 1L
  out
}

# Mean multisite beta_SIM over fixed subsets of rows (sites) of a 0/1 matrix.
mean_beta_subsets <- function(M, subs) {
  if (nrow(subs) == 3) {
    S <- rowSums(M); A <- tcrossprod(M)
    i <- subs[1, ]; j <- subs[2, ]; k <- subs[3, ]
    Aij <- A[cbind(i, j)]; Aik <- A[cbind(i, k)]; Ajk <- A[cbind(j, k)]
    mins <- pmin(S[i] - Aij, S[j] - Aij) + pmin(S[i] - Aik, S[k] - Aik) + pmin(S[j] - Ajk, S[k] - Ajk)
    tri <- rowSums(M[i, , drop = FALSE] * M[j, , drop = FALSE] * M[k, , drop = FALSE])
    S_T <- S[i] + S[j] + S[k] - Aij - Aik - Ajk + tri
    shared <- S[i] + S[j] + S[k] - S_T
    mean(mins / (shared + mins))
  } else {
    mean(apply(subs, 2, function(ix) {
      beta_multi_all(lapply(ix, function(r) which(M[r, ] == 1L)))[["beta_SIM"]]
    }))
  }
}

null_res <- bind_rows(lapply(provinces, function(p) bind_rows(lapply(1:5, function(s) {
  d <- filter(recs, Province == p, Stage_Number == s)
  sites <- sort(unique(d$Site)); k <- length(sites)
  base <- data.frame(Province = p, Stage_Number = s, Stage = stage_names[s], Mid_Ma = stage_mid[s],
                     n_sites = k, observed = NA_real_, null_mean = NA_real_, null_sd = NA_real_,
                     SES = NA_real_, p_more_distinct = NA_real_, p_more_alike = NA_real_)
  if (k < n_units_sample) return(base)
  spp <- sort(unique(d$GenusSpecies))
  M <- matrix(0L, k, length(spp))
  M[cbind(match(d$Site, sites), match(d$GenusSpecies, spp))] <- 1L
  subs <- if (choose(k, n_units_sample) <= n_subsets_null) combn(k, n_units_sample) else
    replicate(n_subsets_null, sample.int(k, n_units_sample))
  obs <- mean_beta_subsets(M, subs)
  nulls <- replicate(n_null, mean_beta_subsets(curveball(M), subs))
  base$observed <- obs; base$null_mean <- mean(nulls); base$null_sd <- sd(nulls)
  base$SES <- if (sd(nulls) > 0) (obs - mean(nulls)) / sd(nulls) else NA
  base$p_more_distinct <- (sum(nulls >= obs) + 1) / (n_null + 1)
  base$p_more_alike <- (sum(nulls <= obs) + 1) / (n_null + 1)
  cat(sprintf("  %-16s %-10s %3d sites  observed %.3f  null %.3f  SES %6.2f\n",
              p, stage_names[s], k, obs, mean(nulls), base$SES))
  base
}))))
save_csv(mutate(null_res, across(where(is.numeric), ~ round(.x, 4))), "null_model")

ses_plot_data <- null_res %>% mutate(Province = factor(Province, levels = provinces), value = SES)
y_lim <- range(c(-3, 3, ses_plot_data$SES), na.rm = TRUE) + c(-0.5, 0.8)

p_null <- ggplot(ses_plot_data) +
  stage_bands(y_lab = y_lim[2] - 0.35) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96, fill = "#e4e4df", alpha = 0.7) +
  geom_hline(yintercept = 0, colour = ink_soft, linewidth = 0.5, linetype = "dashed") +
  geom_segment(data = gap_segments(ses_plot_data, "Province"),
               aes(x = Mid_Ma, xend = x2, y = value, yend = y2, colour = Province), linewidth = 1.2) +
  geom_point(aes(Mid_Ma, SES, colour = Province), shape = 21, fill = "white", size = 3.6,
             stroke = 1.6, na.rm = TRUE) +
  facet_wrap(~ Province, ncol = 2) +
  scale_colour_manual(values = province_colours) +
  age_axis() +
  scale_y_continuous(limits = y_lim) +
  labs(title = "Are sites more distinct than chance would predict?",
       subtitle = "Standardized effect size of within-province beta_SIM against a null model; falling towards 0 = homogenization",
       x = "Age (Ma)", y = "SES of beta_SIM",
       caption = paste0("Null model: species reshuffled among sites keeping site richness and species frequencies (curveball, ",
                        n_null, " matrices; ", n_subsets_null, " subsets of ", n_units_sample,
                        " sites).\nGrey band = |SES| < 1.96 (not different from chance). Above: sites more distinct than random; ",
                        "below: more alike. No point = fewer than ", n_units_sample, " sites.")) +
  theme_slide() +
  theme(panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5))
save_fig(p_null, "null_model_ses", height = slide_h * 1.1)

# =============================================================================
# 4. TREND TEST (c)
# =============================================================================

cat("\n=== 4. TREND TEST (beta_SIM ~ age) ===\n")

s7_num <- s7 %>% mutate(Mid_Ma = stage_mid[match(Stage, stage_names)]) %>% filter(!is.na(beta_SIM))
trend_test <- bind_rows(lapply(provinces, function(p) {
  d <- filter(s7_num, Province == p)
  if (nrow(d) < 3) {
    return(data.frame(Province = p, n_stages = nrow(d), slope_per_Myr = NA, r2 = NA, p = NA,
                      spearman_rho = NA, Reading = "fewer than 3 stages"))
  }
  m <- summary(lm(beta_SIM ~ Mid_Ma, data = d))
  sl <- coef(m)[2, 1]; pv <- coef(m)[2, 4]
  data.frame(Province = p, n_stages = nrow(d), slope_per_Myr = round(sl, 4),
             r2 = round(m$r.squared, 3), p = round(pv, 4),
             spearman_rho = round(suppressWarnings(cor(d$Mid_Ma, d$beta_SIM, method = "spearman")), 3),
             Reading = ifelse(pv < 0.05 & sl > 0, "beta_SIM declines towards the present (homogenization)",
                       ifelse(pv < 0.05 & sl < 0, "beta_SIM rises towards the present",
                              "no significant linear trend")))
}))
save_csv(trend_test, "trend_test")
print(trend_test, row.names = FALSE)
cat("  (positive slope = higher beta_SIM in older stages, i.e. falling towards the present)\n")

cat("\n=== STEP 9b COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: sampling, fits, null_res, trend_test\n")
