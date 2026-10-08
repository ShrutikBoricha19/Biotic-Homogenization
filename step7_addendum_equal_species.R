# =============================================================================
# STEP 7 ADDENDUM: MULTISITE beta_SIM WITH EQUAL SPECIES NUMBERS PER PROVINCE
#
# The three provinces (Basin and Range, Coastal Plain, Great Plains incl. the
# Central Lowland west of the Mississippi) hold very different numbers of
# species in each time bin. Here every province is given the SAME number of
# species: in each time bin, n = the species count of the poorest province.
# n species are drawn at random from each province's own species list (the
# poorest province keeps all of its species), the multisite Simpson
# dissimilarity (Rowan et al. 2024, Eq. 1) of the three provinces is computed,
# and this is repeated 'n_resamples' times.
#
#   Points      mean beta_SIM of the random draws
#   Error bars  95% range of the draws (2.5th-97.5th percentiles)
#   p (above n) p value for the difference from the previous (older) time bin
#
# p values: the error bars only show how beta_SIM changes with WHICH species
# are drawn, not how it would change if a different set of sites had been
# found. The tests therefore use a bootstrap that includes both: in each of
# 'n_boot' iterations the sites of every province are resampled with
# replacement, and every province is then cut to the poorest province's
# species count as above. For successive bins A and B, p (two-sided) = 2 x the
# share of iterations in which beta_SIM(A) - beta_SIM(B) falls on the other
# side of zero; the p values are Holm-corrected over the successive pairs.
#
# Time bins are those of Step 3c (Bin 1 = 3.25-2.50 Ma); the age axis starts
# at the older limit of Bin 1.
#
# Outputs (Outputs/7_rowan_figures/):
#   Fig2_equal_species.png/.pdf       beta_SIM per bin with error bars
#   beta_sim_equal_species.csv        per bin: n species used, mean, SD, 95% range
#   beta_sim_equal_species_tests.csv  successive bins: difference, 95% CI, p, p (Holm)
#
# Inputs: Outputs/7_rowan_figures/regional_pa_matrices.rds and species_traits.csv (Step 7),
#         Outputs/5b_great_plains/master_data_unique.csv (Step 5b; sites for the tests),
#         Outputs/3d_resolved/time_bins.csv.
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 7.
# Needs: dplyr, ggplot2.
# =============================================================================

for (pkg in c("dplyr", "ggplot2")) {
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

pa_file     <- file.path("Outputs", "7_rowan_figures", "regional_pa_matrices.rds")
traits_file <- file.path("Outputs", "7_rowan_figures", "species_traits.csv")
master_file <- file.path("Outputs", "5b_great_plains", "master_data_unique.csv")
bins_file  <- file.path("Outputs", "3d_resolved", "time_bins.csv")
output_dir <- file.path(work_dir, "Outputs", "7_rowan_figures")

regions      <- c("Basin and Range", "Coastal Plain", "Great Plains")
n_resamples  <- 999
n_boot       <- 9999       # bootstrap iterations per bin for the p values (smallest p = 2/(n_boot+1))
alpha        <- 0.05
random_seed  <- 2024

band_fill   <- "#F6EBD0"    # shading of alternate time bins (pale sand, as Fig. 2)
# Colourblind-friendly colours, no greys or whites in the plots (as Step 7).
ink         <- "#1F2A44"    # text, axes, lines and outlines (dark navy)
ink_soft    <- "#3E4C6D"    # secondary text
point_fill  <- "#F0E442"    # beta_SIM points (Okabe-Ito yellow, navy outline)
page_fill   <- "white"      # figure background (change here, e.g. "#F7FAFD", for a tinted page)
fig_w <- 9; fig_h <- 5.6; fig_dpi <- 300; base_size <- 14

# =============================================================================
# HELPERS
# =============================================================================

# Multisite Simpson dissimilarity (Rowan et al. 2024, Eq. 1; Baselga 2010).
beta_sim <- function(x) {
  x <- x[, colSums(x) > 0, drop = FALSE]
  cmb <- combn(nrow(x), 2)
  mins <- vapply(seq_len(ncol(cmb)), function(k) {
    i <- cmb[1, k]; j <- cmb[2, k]
    min(sum(x[i, ] == 1 & x[j, ] == 0), sum(x[i, ] == 0 & x[j, ] == 1))
  }, numeric(1))
  shared <- sum(x) - ncol(x)                       # sum_i S_i - S_T
  sum(mins) / (shared + sum(mins))
}
fmt_age <- function(x) ifelse(x < 0.1, sprintf("%.4f", x), sprintf("%.2f", x))
fmt_p   <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
stars   <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

# beta_SIM of province species lists, each cut to n species at random.
beta_equal <- function(sp_lists, n) {
  picks <- lapply(sp_lists, function(s) if (length(s) > n) sample(s, n) else s)
  spp <- sort(unique(unlist(picks)))
  m <- t(vapply(picks, function(s) as.integer(spp %in% s), integer(length(spp))))
  beta_sim(m)
}
read_csv_utf8 <- function(f) {
  path <- file.path(work_dir, f)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun the earlier steps first.")
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
}

# =============================================================================
# 1. READ STEP 7 MATRICES
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
pa_path <- file.path(work_dir, pa_file)
if (!file.exists(pa_path)) stop("File not found:\n  ", pa_path, "\nRun Step 7 first.")
pa_list <- readRDS(pa_path)
time_bins <- read.csv(file.path(work_dir, bins_file), stringsAsFactors = FALSE) %>%
  transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma), Younger = as.numeric(Younger_Ma)) %>%
  arrange(desc(Older)) %>%
  mutate(Display_Bin = Time_Bin, Mid = (Older + Younger) / 2)
cat(sprintf("  %d time bins\n", nrow(time_bins)))

# =============================================================================
# 2. EQUAL-SPECIES RESAMPLING
# =============================================================================

cat(sprintf("\n=== 2. EQUAL SPECIES NUMBERS (%d random draws per bin) ===\n", n_resamples))
set.seed(random_seed)
res <- list()
for (k in seq_len(nrow(time_bins))) {
  b <- time_bins$Time_Bin[k]
  x <- pa_list[[as.character(b)]]
  if (is.null(x) || nrow(x) < 2) { cat(sprintf("  Bin %d: fewer than two provinces - skipped\n", time_bins$Display_Bin[k])); next }
  sp_lists <- lapply(seq_len(nrow(x)), function(i) colnames(x)[x[i, ] == 1])
  names(sp_lists) <- rownames(x)
  n_min <- min(lengths(sp_lists))
  draws <- replicate(n_resamples, {
    picks <- lapply(sp_lists, function(s) if (length(s) > n_min) sample(s, n_min) else s)
    spp <- sort(unique(unlist(picks)))
    m <- t(vapply(picks, function(s) as.integer(spp %in% s), integer(length(spp))))
    colnames(m) <- spp
    beta_sim(m)
  })
  res[[length(res) + 1]] <- data.frame(
    Display_Bin = time_bins$Display_Bin[k],
    Older = time_bins$Older[k], Younger = time_bins$Younger[k], Mid = time_bins$Mid[k],
    n_provinces = nrow(x), provinces = paste(rownames(x), collapse = "; "),
    species_per_province = paste(lengths(sp_lists), collapse = "/"),
    n_species_used = n_min,
    beta_SIM_mean = mean(draws), beta_SIM_sd = sd(draws),
    beta_SIM_lower95 = unname(quantile(draws, 0.025)), beta_SIM_upper95 = unname(quantile(draws, 0.975)))
}
res <- bind_rows(res)
print(as.data.frame(res %>% transmute(Bin = Display_Bin, Ages = sprintf("%s-%s Ma", fmt_age(Older), fmt_age(Younger)),
                                      species_per_province, n_species_used,
                                      mean = round(beta_SIM_mean, 3), lower95 = round(beta_SIM_lower95, 3),
                                      upper95 = round(beta_SIM_upper95, 3))), row.names = FALSE)
if (any(res$n_provinces < length(regions))) {
  cat("  NOTE: bins with fewer than three provinces use the provinces present.\n")
}

# =============================================================================
# 2b. p VALUES BETWEEN SUCCESSIVE BINS (sites resampled, then equal species)
# =============================================================================

cat(sprintf("\n=== 2b. TESTS BETWEEN SUCCESSIVE BINS (%d bootstrap iterations per bin) ===\n", n_boot))
traits <- read_csv_utf8(traits_file)
large  <- traits$Species[as.logical(traits$Included) %in% TRUE]
occ <- read_csv_utf8(master_file) %>%
  mutate(Time_Bin = as.integer(Time_Bin)) %>%
  filter(Spatial_Bin %in% regions, Species %in% large) %>%
  distinct(Time_Bin, Spatial_Bin, Site_Key, Species)

set.seed(random_seed)
boot <- list()
for (k in seq_len(nrow(res))) {
  b <- res$Display_Bin[k]
  d <- occ %>% filter(Time_Bin == b)
  regs <- intersect(regions, unique(d$Spatial_Bin))
  site_lists <- lapply(setNames(regs, regs), function(r) {
    dr <- d %>% filter(Spatial_Bin == r)
    split(dr$Species, dr$Site_Key)
  })
  boot[[as.character(b)]] <- replicate(n_boot, {
    sp_lists <- lapply(site_lists, function(sites) unique(unlist(sites[sample.int(length(sites), replace = TRUE)])))
    beta_equal(sp_lists, min(res$n_species_used[k], lengths(sp_lists)))
  })
}

tests <- bind_rows(lapply(seq_len(nrow(res) - 1), function(k) {
  a <- res$Display_Bin[k]; bb <- res$Display_Bin[k + 1]
  diff <- boot[[as.character(a)]] - boot[[as.character(bb)]]
  p <- min(1, 2 * min((sum(diff <= 0) + 1) / (n_boot + 1), (sum(diff >= 0) + 1) / (n_boot + 1)))
  data.frame(Bin_A = a, Bin_B = bb, diff_mean = mean(diff),
             diff_lower95 = unname(quantile(diff, 0.025)), diff_upper95 = unname(quantile(diff, 0.975)), p = p)
})) %>%
  mutate(p_holm = p.adjust(p, method = "holm"), Significant = p_holm < alpha)
print(as.data.frame(tests %>% transmute(Pair = sprintf("Bin %d vs Bin %d", Bin_A, Bin_B),
                                        Difference = sprintf("%+.3f", diff_mean),
                                        CI95 = sprintf("%+.3f to %+.3f", diff_lower95, diff_upper95),
                                        p = fmt_p(p), p_Holm = fmt_p(p_holm),
                                        Result = ifelse(Significant, "significant", "not significant"))),
      row.names = FALSE)
cat("  Difference = beta_SIM(older bin) - beta_SIM(younger bin).\n")

# =============================================================================
# 3. FIGURE (style of Fig. 2)
# =============================================================================

cat("\n=== 3. FIGURE ===\n")
x_max <- max(time_bins$Older)
bands <- time_bins %>% filter(Display_Bin %% 2 == 1)
# Column labels (bottom up): n, p against the previous (older) bin, bin name.
bin_labels <- res %>%
  left_join(tests %>% transmute(Display_Bin = Bin_B, p_holm, Significant), by = "Display_Bin") %>%
  transmute(Mid,
            bin_lab = paste0("Bin ", Display_Bin,
                             ifelse(n_provinces < length(regions), sprintf(" (%d provinces)", n_provinces), "")),
            p_lab = ifelse(is.na(p_holm), "", trimws(paste(fmt_p(p_holm), stars(p_holm)))),
            p_bold = Significant %in% TRUE,
            n_lab = sprintf("n = %d", n_species_used))

p <- ggplot(res, aes(Mid, beta_SIM_mean)) +
  geom_rect(data = bands, aes(xmin = Younger, xmax = Older, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = band_fill, colour = NA) +
  geom_errorbar(aes(ymin = beta_SIM_lower95, ymax = beta_SIM_upper95), width = 0.07, colour = ink, linewidth = 0.55) +
  geom_line(colour = ink, linewidth = 0.6) +
  geom_point(shape = 21, fill = point_fill, colour = ink, size = 3.6, stroke = 0.8) +
  geom_text(data = bin_labels, aes(Mid, 0.02, label = n_lab), vjust = 0, inherit.aes = FALSE,
            size = base_size / 4.3, colour = ink_soft) +
  geom_text(data = bin_labels, aes(Mid, 0.065, label = p_lab), vjust = 0, inherit.aes = FALSE,
            size = base_size / 4.3, colour = ink, fontface = ifelse(bin_labels$p_bold, "bold", "plain")) +
  geom_text(data = bin_labels, aes(Mid, 0.11, label = bin_lab), vjust = 0, inherit.aes = FALSE,
            size = base_size / 4.3, colour = ink_soft) +
  # right-hand guide, as in Fig. 2
  annotate("segment", x = -0.24, xend = -0.24, y = 0.53, yend = 0.98, colour = ink, linewidth = 0.5,
           arrow = arrow(length = unit(0.1, "in"), type = "open")) +
  annotate("segment", x = -0.24, xend = -0.24, y = 0.47, yend = 0.02, colour = ink, linewidth = 0.5,
           arrow = arrow(length = unit(0.1, "in"), type = "open")) +
  annotate("text", x = -0.34, y = 0.755, label = "Higher \u03b2 (provincialism)", angle = 90,
           size = base_size / 4.1, colour = ink) +
  annotate("text", x = -0.34, y = 0.245, label = "Lower \u03b2 (homogenization)", angle = 90,
           size = base_size / 4.1, colour = ink) +
  scale_x_reverse(breaks = c(x_max, seq(floor(x_max), 0, by = -1)), labels = function(v) sub("\\.?0+$", "", sprintf("%.2f", v))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(v) sprintf("%.2f", v)) +
  coord_cartesian(xlim = c(x_max, 0), ylim = c(0, 1), clip = "off") +
  labs(x = "Age (Ma)", y = expression(beta[SIM]),
       caption = sprintf(paste0("Each province subsampled to the species count of the poorest province (n). Points = mean of %d random draws;\n",
                                "error bars = 95%% range of the draws. p (above n) = difference from the previous, older bin\n",
                                "(bootstrap of sites, then equal species; %d iterations; Holm-corrected).\n",
                                "Bold with * p < 0.05, ** p < 0.01, *** p < 0.001 = significant."), n_resamples, n_boot)) +
  theme_classic(base_size = base_size) +
  theme(axis.line = element_line(colour = ink, linewidth = 0.4),
        axis.ticks = element_line(colour = ink, linewidth = 0.4),
        axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
        axis.title.y = element_text(size = base_size * 1.15),
        plot.caption = element_text(colour = ink_soft, size = base_size * 0.62, hjust = 0),
        plot.margin = margin(10, 40, 8, 10), plot.background = element_rect(fill = page_fill, colour = NA))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "Fig2_equal_species.png"), p, width = fig_w, height = fig_h, dpi = fig_dpi, bg = page_fill)
ggsave(file.path(output_dir, "Fig2_equal_species.pdf"), p, width = fig_w, height = fig_h, bg = page_fill,
       device = if (capabilities("cairo")) cairo_pdf else pdf)
cat("  ", file.path(output_dir, "Fig2_equal_species.png"), " (+ .pdf)\n", sep = "")

write.csv(res %>% mutate(across(starts_with("beta"), ~ round(.x, 4))),
          file.path(output_dir, "beta_sim_equal_species.csv"), row.names = FALSE)
cat("  ", file.path(output_dir, "beta_sim_equal_species.csv"), "\n", sep = "")
write.csv(tests %>% mutate(across(starts_with("diff"), ~ round(.x, 4)), across(c(p, p_holm), ~ signif(.x, 3))),
          file.path(output_dir, "beta_sim_equal_species_tests.csv"), row.names = FALSE)
cat("  ", file.path(output_dir, "beta_sim_equal_species_tests.csv"), "\n", sep = "")
cat("\n=== STEP 7 ADDENDUM COMPLETE ===\n")
