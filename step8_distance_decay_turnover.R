# =============================================================================
# STEP 8: DISTANCE DECAY, TURNOVER vs NESTEDNESS, AND SITES & SPECIES PER BIN
#
# Within each province (Basin and Range, Coastal Plain, Great Plains incl. the
# Central Lowland west of the Mississippi) and each time bin, using the sites
# of that province:
#
#   1. Distance decay of similarity. For every pair of sites: Simpson
#      similarity (1 - beta_sim = a / (a + min(b, c))) against the great-circle
#      distance between them (km). Line: least-squares fit; Mantel test
#      (Pearson, 'n_perm' permutations, vegan) for the correlation between
#      the similarity and distance matrices. A negative slope = sites further
#      apart share fewer species.
#
#   2. Baselga's partition of multiple-site dissimilarity (Baselga 2010,
#      Global Ecol. Biogeogr. 19:134-143):
#        beta_SOR = beta_SIM (turnover, replacement of species)
#                 + beta_SNE (nestedness, species loss / richness difference)
#      Multi-site values grow with the number of sites, so - as Baselga
#      recommends - each pool is subsampled to the same number of sites
#      ('n_sites_sample') 'n_resamples' times and the values are averaged
#      (pools with fewer sites use all their sites and are marked).
#
#   3. Number of sites and number of species per time bin, for the three
#      provinces in one figure.
#
# Species: the large mammals of Step 7 (species_traits.csv, Included = TRUE),
# so that the results match Rowan-style Figs. 2-4; set large_mammals_only to
# FALSE to use every species.
#
# Outputs (Outputs/8_spatial_beta/):
#   distance_decay.png/.pdf              3 provinces x time bins
#   turnover_nestedness.png/.pdf         stacked columns per province and bin
#   sites_species_per_bin.png/.pdf       sites and species, all provinces
#   distance_decay_stats.csv, turnover_nestedness.csv, sites_species_per_bin.csv
#   site_pairs.csv                       every site pair: distance and similarity
#
# Inputs: Outputs/5b_great_plains/master_data_unique.csv (Step 5b),
#         Outputs/3d_resolved/time_bins.csv, Outputs/7_rowan_figures/species_traits.csv (Step 7).
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 7.
# Needs: dplyr, tidyr, ggplot2, vegan (for the Mantel test).
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "vegan")) {
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

master_file <- file.path("Outputs", "5b_great_plains", "master_data_unique.csv")
bins_file   <- file.path("Outputs", "3d_resolved", "time_bins.csv")
traits_file <- file.path("Outputs", "7_rowan_figures", "species_traits.csv")
output_dir  <- file.path(work_dir, "Outputs", "8_spatial_beta")

provinces <- c("Basin and Range", "Coastal Plain", "Great Plains")
province_colours <- c("Basin and Range" = "#EB6834", "Coastal Plain" = "#4A3AA7", "Great Plains" = "#1BAF7A")

large_mammals_only <- TRUE   # same species as Step 7; FALSE = all species
min_sites          <- 3      # fewer sites in a pool: no distance decay / partition
n_perm             <- 999    # Mantel permutations
n_sites_sample     <- 5      # sites per subsample for Baselga's multi-site partition
n_resamples        <- 100
random_seed        <- 2024

turnover_colour   <- "#2A5C8A"
nestedness_colour <- "#9CC3E4"
base_size <- 13
fig_dpi   <- 300

# =============================================================================
# HELPERS
# =============================================================================

read_text_csv <- function(path) {
  full <- file.path(work_dir, path)
  if (!file.exists(full)) stop("File not found:\n  ", full, "\nRun the earlier steps first.")
  df <- read.csv(full, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
                 na.strings = c("", "NA"), encoding = "UTF-8")
  names(df) <- trimws(sub("^[^A-Za-z0-9]+", "", names(df), useBytes = TRUE))
  cat(sprintf("  %-50s %7d rows\n", path, nrow(df)))
  df
}
write_out <- function(df, name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, paste0(name, ".csv"))
  ok <- tryCatch({ write.csv(df, path, row.names = FALSE, na = "", fileEncoding = "UTF-8"); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) stop("Could not write ", path, "\n  Close it (e.g. in Excel) and run again.")
  cat(sprintf("  %-26s %7d rows -> %s\n", name, nrow(df), path))
}
save_fig <- function(p, name, width, height) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, paste0(name, ".png")), p, width = width, height = height, dpi = fig_dpi, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), p, width = width, height = height, bg = "white",
         device = if (capabilities("cairo")) cairo_pdf else pdf)
  cat(sprintf("  %s.png (+ .pdf)\n", file.path(output_dir, name)))
}
fmt_age <- function(x) ifelse(x < 0.1, sprintf("%.4f", x), sprintf("%.2f", x))

# Great-circle distance (km) between all pairs of points.
haversine_km <- function(lat, lon) {
  la <- lat * pi / 180; lo <- lon * pi / 180
  dlat <- outer(la, la, "-"); dlon <- outer(lo, lo, "-")
  a <- sin(dlat / 2)^2 + outer(cos(la), cos(la)) * sin(dlon / 2)^2
  2 * 6371.0088 * asin(pmin(sqrt(a), 1))   # sqrt(a) first, so the matrix shape is kept
}

# Baselga (2010) multiple-site Sorensen dissimilarity and its two components.
baselga_multi <- function(x) {
  x <- x[, colSums(x) > 0, drop = FALSE]
  s_i <- sum(x); s_t <- ncol(x)
  shared <- tcrossprod(x)                          # a_ij
  rich <- rowSums(x)
  b <- outer(rich, rep(1, nrow(x))) - shared       # b_ij = species in i not in j
  up <- upper.tri(b)
  mn <- sum(pmin(b, t(b))[up]); mx <- sum(pmax(b, t(b))[up])
  sor <- (mn + mx) / (2 * (s_i - s_t) + mn + mx)
  sim <- mn / (mn + s_i - s_t)
  c(beta_SOR = sor, beta_SIM = sim, beta_SNE = sor - sim)
}

theme_fig <- function() {
  theme_bw(base_size = base_size) +
    theme(strip.background = element_rect(fill = "#E3E3E3", colour = "grey35"),
          strip.text = element_text(face = "bold", colour = "grey10"),
          panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92"),
          plot.title = element_text(face = "bold", size = base_size * 1.3), plot.title.position = "plot",
          plot.subtitle = element_text(colour = "grey30", size = base_size * 0.82),
          plot.caption = element_text(colour = "grey35", size = base_size * 0.7, hjust = 0),
          plot.background = element_rect(fill = "white", colour = NA))
}

# =============================================================================
# 1. READ INPUTS
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
md <- read_text_csv(master_file) %>%
  filter(Spatial_Bin %in% provinces) %>%
  mutate(Time_Bin = as.integer(Time_Bin), Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude))
time_bins <- read_text_csv(bins_file) %>%
  transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma), Younger = as.numeric(Younger_Ma)) %>%
  arrange(Time_Bin) %>%
  mutate(Bin_Label = sprintf("Bin %d\n%s-%s Ma", Time_Bin, fmt_age(Older), fmt_age(Younger)),
         Bin_Short = sprintf("Bin %d (%s-%s Ma)", Time_Bin, fmt_age(Older), fmt_age(Younger)))

if (large_mammals_only) {
  tr <- read_text_csv(traits_file)
  keep <- tr$Species[toupper(tr$Included) == "TRUE"]
  md <- md %>% filter(Species %in% keep)
  cat(sprintf("  Large mammals only (Step 7): %d species\n", length(keep)))
}
md <- md %>% distinct(Spatial_Bin, Time_Bin, Site_Key, Species, .keep_all = TRUE)

# =============================================================================
# 2. SITE x SPECIES MATRICES, DISTANCE DECAY AND BASELGA PARTITION
# =============================================================================

cat("\n=== 2. ANALYSES PER PROVINCE AND TIME BIN ===\n")
set.seed(random_seed)
pairs_all <- list(); dd_stats <- list(); bp <- list(); counts <- list()
for (p in provinces) for (b in time_bins$Time_Bin) {
  d <- md %>% filter(Spatial_Bin == p, Time_Bin == b)
  sites <- d %>% distinct(Site_Key, .keep_all = TRUE) %>% arrange(Site_Key)
  counts[[length(counts) + 1]] <- data.frame(Province = p, Time_Bin = b, n_sites = nrow(sites),
                                             n_species = n_distinct(d$Species))
  if (nrow(sites) < min_sites) next
  spp <- sort(unique(d$Species))
  x <- matrix(0L, nrow(sites), length(spp), dimnames = list(sites$Site_Key, spp))
  x[cbind(match(d$Site_Key, sites$Site_Key), match(d$Species, spp))] <- 1L

  # --- distance decay (sites with coordinates) ---
  ok <- !is.na(sites$Latitude) & !is.na(sites$Longitude)
  if (sum(ok) >= min_sites) {
    xs <- x[ok, , drop = FALSE]
    dist_km <- haversine_km(sites$Latitude[ok], sites$Longitude[ok])
    a <- tcrossprod(xs); rich <- rowSums(xs)
    bmin <- pmin(outer(rich, rep(1, nrow(xs))) - a, t(outer(rich, rep(1, nrow(xs))) - a))
    sim <- ifelse(a + bmin > 0, a / (a + bmin), NA)       # Simpson similarity = 1 - beta_sim
    up <- which(upper.tri(dist_km), arr.ind = TRUE)
    pr <- data.frame(Province = p, Time_Bin = b, Site_1 = rownames(xs)[up[, 1]], Site_2 = rownames(xs)[up[, 2]],
                     Distance_km = dist_km[up], Similarity = sim[up])
    pairs_all[[length(pairs_all) + 1]] <- pr
    fit <- lm(Similarity ~ Distance_km, data = pr)
    man <- tryCatch(suppressMessages(vegan::mantel(as.dist(1 - sim), as.dist(dist_km), method = "pearson",
                                  permutations = n_perm, na.rm = TRUE)), error = function(e) NULL)
    dd_stats[[length(dd_stats) + 1]] <- data.frame(
      Province = p, Time_Bin = b, n_sites = nrow(xs), n_pairs = nrow(pr),
      mean_similarity = mean(pr$Similarity, na.rm = TRUE),
      intercept = coef(fit)[[1]], slope_per_1000km = coef(fit)[[2]] * 1000,
      r_squared = summary(fit)$r.squared,
      mantel_r = if (is.null(man)) NA else man$statistic,      # dissimilarity vs distance: positive = decay
      mantel_p = if (is.null(man)) NA else man$signif)
  }

  # --- Baselga partition, subsampled to n_sites_sample sites ---
  n <- nrow(x)
  if (n > n_sites_sample) {
    reps <- t(replicate(n_resamples, baselga_multi(x[sample.int(n, n_sites_sample), , drop = FALSE])))
    est <- colMeans(reps); sds <- apply(reps, 2, sd); used <- n_sites_sample; sub <- TRUE
  } else {
    est <- baselga_multi(x); sds <- c(beta_SOR = NA, beta_SIM = NA, beta_SNE = NA); used <- n; sub <- FALSE
  }
  allv <- baselga_multi(x)
  bp[[length(bp) + 1]] <- data.frame(
    Province = p, Time_Bin = b, n_sites = n, sites_per_sample = used, subsampled = sub,
    beta_SOR = est[["beta_SOR"]], beta_SIM = est[["beta_SIM"]], beta_SNE = est[["beta_SNE"]],
    sd_SOR = sds[["beta_SOR"]], sd_SIM = sds[["beta_SIM"]], sd_SNE = sds[["beta_SNE"]],
    all_sites_beta_SOR = allv[["beta_SOR"]], all_sites_beta_SIM = allv[["beta_SIM"]],
    all_sites_beta_SNE = allv[["beta_SNE"]],
    turnover_share = est[["beta_SIM"]] / est[["beta_SOR"]])
}
pairs_all <- bind_rows(pairs_all); dd_stats <- bind_rows(dd_stats)
bp <- bind_rows(bp); counts <- bind_rows(counts)

cat("  Distance decay:\n")
print(as.data.frame(dd_stats %>% transmute(Province, Time_Bin, n_sites, n_pairs,
                                           slope_per_1000km = round(slope_per_1000km, 3),
                                           mantel_r = round(mantel_r, 2), mantel_p)), row.names = FALSE)
cat("  Baselga partition (mean of subsamples):\n")
print(as.data.frame(bp %>% transmute(Province, Time_Bin, n_sites, sites_per_sample,
                                     beta_SOR = round(beta_SOR, 3), beta_SIM = round(beta_SIM, 3),
                                     beta_SNE = round(beta_SNE, 3), turnover_share = round(turnover_share, 2))),
      row.names = FALSE)

# =============================================================================
# 3. FIGURES
# =============================================================================

cat("\n=== 3. FIGURES ===\n")
lab_bin <- setNames(time_bins$Bin_Label, time_bins$Time_Bin)
as_facets <- function(df) df %>% mutate(Province = factor(Province, levels = provinces),
                                        Bin = factor(lab_bin[as.character(Time_Bin)], levels = time_bins$Bin_Label))

# ---- 3a. distance decay ----
all_cells <- as_facets(tidyr::expand_grid(Province = provinces, Time_Bin = time_bins$Time_Bin)) %>%
  left_join(as_facets(dd_stats) %>% select(Province, Bin, n_sites, mantel_r, mantel_p, slope_per_1000km),
            by = c("Province", "Bin")) %>%
  left_join(as_facets(counts) %>% select(Province, Bin, n_sites_all = n_sites), by = c("Province", "Bin")) %>%
  mutate(label = ifelse(is.na(n_sites),
                        sprintf("%d site%s\n(too few)", coalesce(n_sites_all, 0L), ifelse(coalesce(n_sites_all, 0L) == 1, "", "s")),
                        sprintf("Mantel r = %.2f%s\n%d sites", mantel_r,
                                ifelse(is.na(mantel_p), "", ifelse(mantel_p < 0.001, ", P < 0.001", sprintf(", P = %.3f", mantel_p))),
                                n_sites)),
         sig = !is.na(mantel_p) & mantel_p < 0.05)
pp <- as_facets(pairs_all)
p_dd <- ggplot(pp, aes(Distance_km, Similarity)) +
  geom_point(aes(colour = Province), alpha = 0.3, size = 1.3, stroke = 0) +
  geom_smooth(aes(colour = Province), method = "lm", formula = y ~ x, se = TRUE, fill = "grey60",
              linewidth = 0.9) +
  geom_text(data = all_cells, aes(x = Inf, y = Inf, label = label, fontface = ifelse(sig, "bold", "plain")),
            hjust = 1.05, vjust = 1.25, size = 2.9, colour = "grey15", lineheight = 0.95, inherit.aes = FALSE) +
  facet_grid(Province ~ Bin, scales = "free_x", drop = FALSE) +
  scale_colour_manual(values = province_colours, guide = "none") +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1.18)) +   # zoom (keeps the fitted bands intact)
  scale_x_continuous(labels = function(v) format(v, big.mark = ",", trim = TRUE)) +
  labs(title = "Distance decay of faunal similarity within each province",
       subtitle = "Each point is a pair of sites; line = least-squares fit with 95% band. Bold Mantel results: P < 0.05.",
       x = "Distance between sites (km)", y = expression("Simpson similarity (1 -" ~ beta[sim] * ")"),
       caption = sprintf("Mantel test: Pearson correlation between Simpson dissimilarity and distance, %d permutations; positive r = similarity decays with distance.", n_perm)) +
  theme_fig() +
  theme(axis.text.x = element_text(size = base_size * 0.62, angle = 30, hjust = 1))
save_fig(p_dd, "distance_decay", 16, 9)

# ---- 3b. turnover vs nestedness ----
bpl <- as_facets(bp) %>%
  select(Province, Bin, Time_Bin, Turnover = beta_SIM, Nestedness = beta_SNE, beta_SOR, sd_SOR, subsampled, n_sites) %>%
  pivot_longer(c(Turnover, Nestedness), names_to = "Component", values_to = "Value") %>%
  mutate(Component = factor(Component, levels = c("Nestedness", "Turnover")))
bp_tot <- as_facets(bp) %>% mutate(lab = ifelse(subsampled, sprintf("%d", n_sites), sprintf("%d*", n_sites)))
empty <- as_facets(counts) %>% anti_join(as_facets(bp) %>% select(Province, Bin), by = c("Province", "Bin"))
p_bp <- ggplot() +
  geom_col(data = bpl, aes(Bin, Value, fill = Component), width = 0.72, colour = "grey25", linewidth = 0.25) +
  geom_errorbar(data = bp_tot %>% filter(!is.na(sd_SOR)),
                aes(Bin, ymin = pmax(beta_SOR - sd_SOR, 0), ymax = pmin(beta_SOR + sd_SOR, 1)),
                width = 0.25, colour = "grey20", linewidth = 0.4) +
  geom_text(data = bp_tot, aes(Bin, pmin(beta_SOR + coalesce(sd_SOR, 0), 1) + 0.04, label = lab),
            size = 3.1, colour = "grey30") +
  geom_text(data = empty, aes(Bin, 0.05, label = "too few\nsites"), size = 2.8, colour = "grey50", lineheight = 0.9) +
  facet_wrap(~ Province, nrow = 1, drop = FALSE) +
  scale_fill_manual(values = c(Turnover = turnover_colour, Nestedness = nestedness_colour),
                    labels = c(Turnover = expression("Turnover (" * beta[SIM] * ")"),
                               Nestedness = expression("Nestedness (" * beta[SNE] * ")")),
                    breaks = c("Turnover", "Nestedness"), name = NULL) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_x_discrete(labels = function(v) sub("\n.*$", "", v)) +       # "Bin k" only; ages in the caption
  labs(title = "Turnover and nestedness of faunas within each province",
       subtitle = sprintf(paste0("Baselga's partition of multiple-site Sorensen dissimilarity (bar height = %s). ",
                                 "Mean of %d random sets of %d sites (error bar: SD); * = fewer sites, all used."),
                          "\u03b2SOR", n_resamples, n_sites_sample),
       x = NULL, y = "Multiple-site dissimilarity",
       caption = paste0("Numbers above bars: sites in the pool. Turnover = replacement of species between sites; nestedness = sites holding subsets of richer sites.\n",
                        paste(time_bins$Bin_Short, collapse = "   "))) +
  theme_fig() +
  theme(legend.position = "top", legend.text = element_text(size = base_size * 0.9),
        axis.text.x = element_text(size = base_size * 0.85), panel.grid.major.x = element_blank())
save_fig(p_bp, "turnover_nestedness", 16, 7)

# ---- 3c. sites and species per time bin, all provinces ----
cs <- counts %>%
  pivot_longer(c(n_sites, n_species), names_to = "Measure", values_to = "Count") %>%
  mutate(Measure = factor(Measure, levels = c("n_sites", "n_species"),
                          labels = c("Number of sites", "Number of species")),
         Province = factor(Province, levels = provinces),
         Bin = factor(time_bins$Bin_Label[match(Time_Bin, time_bins$Time_Bin)], levels = time_bins$Bin_Label))
p_cs <- ggplot(cs, aes(Bin, Count, fill = Province)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.35, size = 3.2, colour = "grey20") +
  facet_wrap(~ Measure, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = province_colours, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Sites and species in each province through time",
       subtitle = if (large_mammals_only) "Large mammals (as in Step 7); sites = localities and collections holding at least one of them"
                  else "All species; sites = localities and collections",
       x = NULL, y = NULL) +
  theme_fig() +
  theme(legend.position = "top", legend.text = element_text(size = base_size * 0.95),
        panel.grid.major.x = element_blank(), axis.text.x = element_text(size = base_size * 0.8, lineheight = 0.9))
save_fig(p_cs, "sites_species_per_bin", 13.33, 8)

# =============================================================================
# 4. TABLES
# =============================================================================

cat("\n=== 4. TABLES ===\n")
write_out(dd_stats %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), "distance_decay_stats")
write_out(bp %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), "turnover_nestedness")
write_out(counts, "sites_species_per_bin")
write_out(pairs_all %>% mutate(Distance_km = round(Distance_km, 2), Similarity = round(Similarity, 4)), "site_pairs")

cat("\n=== STEP 8 COMPLETE ===\n")
