# =============================================================================
# STEP 9d: MULTIVARIATE DISPERSION - IS THE CLOUD OF SITES SHRINKING?
#
# In each province, every site (of every stage) is placed in a principal
# coordinates analysis (PCoA) of pairwise Simpson dissimilarity. For each
# stage, the distance of each site to the centre (centroid) of its stage is
# measured (method of Anderson et al. 2006, as in vegan::betadisper, including
# negative eigenvalues and the small-sample bias adjustment).
#
#   Large distances = sites are spread out (different faunas).
#   Distances shrinking through time = sites cluster more tightly
#   = homogenization.
#
# Tests (permutation, 'n_perm' permutations):
#   - overall: do the stages differ in dispersion? (F test, labels permuted)
#   - successive stages only: is dispersion different between one stage and
#     the next? (difference in mean distance, labels permuted)
#
# No extra package is needed (the vegan method is built in).
#
# Inputs:  all_records_final.csv (Step 6), site_physio_database.csv (Step 5)
# Outputs (Outputs/9d_dispersion/):
#   dispersion_all_provinces.png/.pdf   distance to centroid by stage, 4 panels
#   ordination_<province>.png/.pdf      PCoA of the sites, coloured by stage
#   dispersion.csv, dispersion_tests.csv
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Steps 5 and 6.
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
output_dir   <- file.path(work_dir, "Outputs", "9d_dispersion")

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

min_sites_stage <- 3      # stages with fewer sites in a province are left out
n_perm          <- 999
set.seed(2024)

# Stages are ordered in time, so they get one hue from light (old) to dark (young).
stage_colours <- setNames(c("#9ec5f4", "#5f9ee9", "#2a78d6", "#184f95", "#0d2f5c"), stage_names)

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
# 2. PCoA AND DISTANCE TO STAGE CENTROID
# =============================================================================

cat("\n=== 2. DISPERSION ===\n")

# betadisper-style distances to group centroids from a dissimilarity matrix.
dispersion <- function(D, group) {
  n <- nrow(D)
  A <- -0.5 * D^2
  J <- diag(n) - matrix(1 / n, n, n)
  e <- eigen(J %*% A %*% J, symmetric = TRUE)
  tol <- max(abs(e$values)) * 1e-8
  pos <- e$values > tol; neg <- e$values < -tol
  Xp <- e$vectors[, pos, drop = FALSE] %*% diag(sqrt(e$values[pos]), sum(pos))
  Xn <- if (any(neg)) e$vectors[, neg, drop = FALSE] %*% diag(sqrt(-e$values[neg]), sum(neg)) else NULL
  z <- numeric(n)
  for (g in unique(group)) {
    ix <- which(group == g)
    dp <- rowSums(sweep(Xp[ix, , drop = FALSE], 2, colMeans(Xp[ix, , drop = FALSE]))^2)
    dn <- if (is.null(Xn)) 0 else rowSums(sweep(Xn[ix, , drop = FALSE], 2, colMeans(Xn[ix, , drop = FALSE]))^2)
    z[ix] <- sqrt(pmax(0, dp - dn)) * sqrt(length(ix) / (length(ix) - 1))   # bias adjustment
  }
  list(z = z, axes = Xp[, 1:min(2, ncol(Xp)), drop = FALSE],
       var_explained = e$values[pos][1:min(2, sum(pos))] / sum(e$values[pos]))
}

f_stat <- function(z, g) {
  g <- factor(g)
  if (nlevels(g) < 2) return(NA_real_)
  anova(lm(z ~ g))[["F value"]][1]
}

disp_sites <- list(); tests <- list()
for (p in provinces) {
  d <- filter(recs, Province == p)
  keep <- d %>% distinct(Stage_Number, Site) %>% count(Stage_Number) %>%
    filter(n >= min_sites_stage) %>% pull(Stage_Number)
  d <- filter(d, Stage_Number %in% keep)
  if (length(keep) < 1) { cat(sprintf("  %-16s too few sites\n", p)); next }
  d <- mutate(d, Unit = paste(Stage_Number, Site, sep = " :: "))
  units <- lapply(split(d$GenusSpecies, d$Unit), unique)
  unit_stage <- as.integer(sub(" :: .*$", "", names(units)))
  n <- length(units)
  D <- matrix(0, n, n)
  for (i in seq_len(n - 1)) for (j in (i + 1):n) D[i, j] <- D[j, i] <- beta_sim_pair(units[[i]], units[[j]])
  res <- dispersion(D, unit_stage)

  disp_sites[[p]] <- data.frame(Province = p, Unit = names(units), Stage_Number = unit_stage,
                                distance = res$z, PCoA1 = res$axes[, 1],
                                PCoA2 = if (ncol(res$axes) > 1) res$axes[, 2] else 0,
                                var1 = res$var_explained[1],
                                var2 = if (length(res$var_explained) > 1) res$var_explained[2] else 0)

  # Overall test.
  F_obs <- f_stat(res$z, unit_stage)
  F_null <- replicate(n_perm, f_stat(res$z, sample(unit_stage)))
  p_all <- if (is.na(F_obs)) NA else (sum(F_null >= F_obs) + 1) / (n_perm + 1)
  tests[[length(tests) + 1]] <- data.frame(Province = p, Comparison = "All stages",
                                           statistic = round(F_obs, 3), p = p_all)
  # Successive stages.
  ks <- sort(unique(unit_stage))
  for (a in seq_along(ks)[-length(ks)]) {
    s1 <- ks[a]; s2 <- ks[a + 1]
    if (s2 - s1 != 1) next
    ix <- unit_stage %in% c(s1, s2); zz <- res$z[ix]; gg <- unit_stage[ix]
    obs <- mean(zz[gg == s2]) - mean(zz[gg == s1])
    nul <- replicate(n_perm, { g2 <- sample(gg); mean(zz[g2 == s2]) - mean(zz[g2 == s1]) })
    tests[[length(tests) + 1]] <- data.frame(
      Province = p, Comparison = paste(stage_names[s1], "to", stage_names[s2]),
      statistic = round(obs, 3), p = (sum(abs(nul) >= abs(obs)) + 1) / (n_perm + 1))
  }
  cat(sprintf("  %-16s %3d sites in %d stages; overall p = %s\n", p, n, length(ks),
              ifelse(is.na(p_all), "n/a", sprintf("%.3f", p_all))))
}
disp_sites <- bind_rows(disp_sites) %>%
  mutate(Province = factor(Province, levels = provinces),
         Stage = factor(stage_names[Stage_Number], levels = stage_names),
         Mid_Ma = stage_mid[Stage_Number])
tests <- bind_rows(tests)

disp_summary <- disp_sites %>%
  group_by(Province, Stage_Number, Stage, Mid_Ma) %>%
  summarise(n_sites = n(), mean_distance = mean(distance), median_distance = median(distance),
            .groups = "drop")
save_csv(mutate(disp_summary, across(where(is.numeric), ~ round(.x, 4))), "dispersion")
save_csv(tests, "dispersion_tests")
cat("\n  Mean distance to stage centroid (falling = homogenization):\n")
print(as.data.frame(disp_summary %>% transmute(Province, Stage, n_sites,
                                               mean_distance = round(mean_distance, 3))), row.names = FALSE)
cat("\n  Tests ('statistic' = F for all stages; later minus earlier mean distance for successive stages):\n")
print(tests, row.names = FALSE)

# =============================================================================
# 3. FIGURES
# =============================================================================

cat("\n=== 3. FIGURES ===\n")

succ_lab <- tests %>%
  filter(Comparison != "All stages") %>%
  mutate(s1 = match(sub(" to .*$", "", Comparison), stage_names),
         x = stage_young[s1],
         label = ifelse(p < 0.001, "p<0.001", sprintf("p=%.2f", p)),
         Province = factor(Province, levels = provinces))
all_lab <- tests %>% filter(Comparison == "All stages") %>%
  mutate(label = sprintf("All stages: p %s", ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))),
         Province = factor(Province, levels = provinces))
y_top <- max(disp_sites$distance, na.rm = TRUE) * 1.38

p_disp <- ggplot(disp_sites) +
  stage_bands() +
  geom_boxplot(aes(x = Mid_Ma, y = distance, group = Stage), width = 0.35, fill = "white",
               colour = ink_soft, outlier.shape = NA, linewidth = 0.5) +
  geom_point(aes(x = Mid_Ma, y = distance, colour = Province),
             position = position_jitter(width = 0.08, height = 0, seed = 1),
             alpha = 0.55, size = 1.8) +
  geom_segment(data = gap_segments(rename(disp_summary, value = mean_distance), "Province"),
               aes(x = Mid_Ma, xend = x2, y = value, yend = y2, colour = Province), linewidth = 1.3) +
  geom_point(data = disp_summary, aes(Mid_Ma, mean_distance, colour = Province), shape = 23,
             fill = "white", size = 3.8, stroke = 1.6) +
  geom_text(data = succ_lab, aes(x = x, y = y_top * 0.84, label = label), size = base_size / 5.2,
            colour = ink_soft) +
  geom_text(data = all_lab, aes(x = 4.65, y = y_top * 0.985, label = label), hjust = 0, vjust = 1,
            size = base_size / 4.6, colour = ink) +
  facet_wrap(~ Province, ncol = 2, drop = FALSE) +
  scale_colour_manual(values = province_colours) +
  age_axis() +
  scale_y_continuous(limits = c(0, y_top), expand = c(0, 0)) +
  labs(title = "Are sites clustering more tightly through time?",
       subtitle = "Distance of each site to the centre of its stage in PCoA space; falling diamonds = homogenization",
       x = "Age (Ma)", y = "Distance to stage centroid",
       caption = paste0("Boxes = sites of each stage; diamonds = mean. p-values at stage boundaries compare successive ",
                        "stages only (", n_perm, " permutations).\nPCoA of pairwise Simpson dissimilarity per province ",
                        "(method of vegan::betadisper). Stages with fewer than ", min_sites_stage,
                        " sites are left out.")) +
  theme_slide() +
  theme(panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5))
save_fig(p_disp, "dispersion_all_provinces", height = slide_h * 1.1)

for (p in provinces) {
  d <- filter(disp_sites, Province == p)
  if (nrow(d) == 0) next
  cen <- d %>% group_by(Stage) %>% summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop")
  pp <- ggplot(d, aes(PCoA1, PCoA2)) +
    geom_hline(yintercept = 0, colour = grid_col) + geom_vline(xintercept = 0, colour = grid_col) +
    geom_point(aes(fill = Stage), shape = 21, colour = "white", size = 3.4, stroke = 0.6, alpha = 0.9) +
    geom_point(data = cen, aes(fill = Stage), shape = 23, colour = ink, size = 6, stroke = 1.2) +
    geom_label(data = cen, aes(label = Stage), nudge_y = 0.04 * diff(range(d$PCoA2)),
               vjust = 0, size = base_size / 4.4, label.size = 0, fill = alpha("white", 0.8)) +
    scale_fill_manual(values = stage_colours, drop = TRUE, name = NULL) +
    coord_equal() +
    labs(title = paste(p, "- sites in faunal space"),
         subtitle = "PCoA of pairwise Simpson dissimilarity; each point is a site, diamonds = centre of each stage",
         x = sprintf("PCoA 1 (%.0f%%)", 100 * d$var1[1]), y = sprintf("PCoA 2 (%.0f%%)", 100 * d$var2[1]),
         caption = "A tighter cloud of a stage's points = more similar sites. Old stages light, young stages dark.") +
    theme_slide() +
    theme(legend.position = "right", panel.grid.major.y = element_blank(),
          panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5))
  save_fig(pp, paste0("ordination_", file_stub(p)))
}
cat("  Figures saved\n")

cat("\n=== STEP 9d COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: disp_sites, disp_summary, tests\n")
