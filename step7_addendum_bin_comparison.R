# =============================================================================
# STEP 7 ADDENDUM 2: ARE THE TIME BINS SIGNIFICANTLY DIFFERENT?
#
# The equal-species figure (step7_addendum_equal_species.R) shows a mean
# beta_SIM per time bin. Its error bars only show how much beta_SIM changes
# with WHICH species are drawn; they do not show how much it would change if
# a different set of sites had been found. To test whether two bins really
# differ, both sources of uncertainty are needed:
#
#   1. Sites: in each iteration, the sites (localities/collections) of every
#      province are resampled with replacement (bootstrap), giving a new
#      species list per province.
#   2. Species: every province is then subsampled to the species count of
#      the poorest province (as in the equal-species figure) and the multisite
#      beta_SIM (Rowan et al. 2024, Eq. 1) is computed.
#
# Repeating this 'n_boot' times gives a distribution of beta_SIM for each bin.
# For every pair of bins A and B the difference beta_SIM(A) - beta_SIM(B) is
# computed for each iteration:
#   95% CI of the difference   2.5th-97.5th percentile of the differences
#   p (two-sided)              2 x the share of differences on the other side
#                              of zero (bootstrap percentile test)
#   p (Holm)                   p corrected for testing all pairs of bins
# A pair is called significantly different when p (Holm) < 0.05.
#
# 'common_n' (SETTINGS): FALSE = n is set separately in each bin (the poorest
# province of that bin; same as the equal-species figure). TRUE = one n for
# all bins (the smallest of those values), so that every bin is compared at
# the same number of species.
#
# Time bins are those of Step 3c (Bin 1 = 3.25-2.50 Ma).
#
# Outputs (Outputs/7_rowan_figures/bin_comparison/):
#   bin_bootstrap_summary.csv       per bin: n, bootstrap mean, SD, 95% CI
#   bin_pairwise_tests.csv          per pair: difference, 95% CI, p, p (Holm)
#   Fig_bin_comparison_trend.png/.pdf    beta_SIM per bin, bootstrap 95% CIs,
#                                        tests between successive bins
#   Fig_bin_comparison_matrix.png/.pdf   Holm-corrected p for all pairs
#
# Inputs: Outputs/5b_great_plains/master_data_unique.csv (Step 5b),
#         Outputs/7_rowan_figures/species_traits.csv (Step 7; large mammals),
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

master_file <- file.path("Outputs", "5b_great_plains", "master_data_unique.csv")
traits_file <- file.path("Outputs", "7_rowan_figures", "species_traits.csv")
bins_file   <- file.path("Outputs", "3d_resolved", "time_bins.csv")
output_dir  <- file.path(work_dir, "Outputs", "7_rowan_figures", "bin_comparison")

regions     <- c("Basin and Range", "Coastal Plain", "Great Plains")
n_boot      <- 9999       # bootstrap iterations per bin (smallest possible p = 2/(n_boot+1))
common_n    <- FALSE      # TRUE = the same species count in every bin
alpha       <- 0.05
random_seed <- 2024

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

# beta_SIM of province species lists, each cut to n species at random.
beta_equal <- function(sp_lists, n) {
  picks <- lapply(sp_lists, function(s) if (length(s) > n) sample(s, n) else s)
  spp <- sort(unique(unlist(picks)))
  m <- t(vapply(picks, function(s) as.integer(spp %in% s), integer(length(spp))))
  beta_sim(m)
}

fmt_age <- function(x) ifelse(x < 0.1, sprintf("%.4f", x), sprintf("%.2f", x))
fmt_p   <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
stars   <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "n.s.")))

read_csv_utf8 <- function(f) {
  path <- file.path(work_dir, f)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun the earlier steps first.")
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
}

# =============================================================================
# 1. READ INPUTS
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
md     <- read_csv_utf8(master_file)
traits <- read_csv_utf8(traits_file)
large  <- traits$Species[as.logical(traits$Included) %in% TRUE]
cat(sprintf("  %d large-mammal species (Step 7)\n", length(large)))

time_bins <- read_csv_utf8(bins_file) %>%
  transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma), Younger = as.numeric(Younger_Ma)) %>%
  arrange(desc(Older)) %>%
  mutate(Display_Bin = Time_Bin, Mid = (Older + Younger) / 2)

occ <- md %>%
  mutate(Time_Bin = as.integer(Time_Bin)) %>%
  filter(Spatial_Bin %in% regions, Species %in% large, Time_Bin %in% time_bins$Time_Bin) %>%
  distinct(Time_Bin, Spatial_Bin, Site_Key, Species)

# Site lists per bin and province: for each site, the species found there.
site_lists <- list()
for (b in time_bins$Time_Bin) {
  d <- occ %>% filter(Time_Bin == b)
  regs <- intersect(regions, unique(d$Spatial_Bin))
  site_lists[[as.character(b)]] <- lapply(setNames(regs, regs), function(r) {
    dr <- d %>% filter(Spatial_Bin == r)
    split(dr$Species, dr$Site_Key)
  })
}

obs <- bind_rows(lapply(seq_len(nrow(time_bins)), function(k) {
  sl <- site_lists[[as.character(time_bins$Time_Bin[k])]]
  rich <- vapply(sl, function(s) length(unique(unlist(s))), integer(1))
  data.frame(Display_Bin = time_bins$Display_Bin[k], Time_Bin = time_bins$Time_Bin[k],
             n_provinces = length(sl),
             sites_per_province = paste(vapply(sl, length, integer(1)), collapse = "/"),
             species_per_province = paste(rich, collapse = "/"),
             n_bin = if (length(rich)) min(rich) else NA_integer_)
}))
ok_bins <- obs$Time_Bin[obs$n_provinces >= 2]
if (any(obs$n_provinces < 2)) {
  cat(sprintf("  Bin(s) %s: fewer than two provinces - left out\n",
              paste(obs$Display_Bin[obs$n_provinces < 2], collapse = ", ")))
}
n_common <- min(obs$n_bin[obs$Time_Bin %in% ok_bins])
obs$n_target <- if (common_n) n_common else obs$n_bin
print(obs %>% select(Bin = Display_Bin, n_provinces, sites_per_province, species_per_province, n_used = n_target),
      row.names = FALSE)

# =============================================================================
# 2. BOOTSTRAP (sites resampled, then equal species numbers)
# =============================================================================

cat(sprintf("\n=== 2. BOOTSTRAP (%d iterations per bin; n = %s) ===\n", n_boot,
            if (common_n) sprintf("%d in every bin", n_common) else "poorest province of each bin"))
set.seed(random_seed)
boot <- list()
for (b in ok_bins) {
  sl <- site_lists[[as.character(b)]]
  n_target <- obs$n_target[obs$Time_Bin == b]
  boot[[as.character(b)]] <- replicate(n_boot, {
    sp_lists <- lapply(sl, function(sites) unique(unlist(sites[sample.int(length(sites), replace = TRUE)])))
    beta_equal(sp_lists, min(n_target, lengths(sp_lists)))
  })
  cat(sprintf("  Bin %d done\n", obs$Display_Bin[obs$Time_Bin == b]))
}

summ <- obs %>% filter(Time_Bin %in% ok_bins) %>%
  left_join(time_bins %>% select(Time_Bin, Older, Younger, Mid), by = "Time_Bin") %>%
  rowwise() %>%
  mutate(beta_SIM_mean = mean(boot[[as.character(Time_Bin)]]),
         beta_SIM_sd = sd(boot[[as.character(Time_Bin)]]),
         lower95 = unname(quantile(boot[[as.character(Time_Bin)]], 0.025)),
         upper95 = unname(quantile(boot[[as.character(Time_Bin)]], 0.975))) %>%
  ungroup()

# =============================================================================
# 3. PAIRWISE TESTS
# =============================================================================

cat("\n=== 3. PAIRWISE TESTS BETWEEN BINS ===\n")
pairs <- combn(summ$Display_Bin, 2)
tests <- bind_rows(lapply(seq_len(ncol(pairs)), function(k) {
  a <- pairs[1, k]; bb <- pairs[2, k]
  da <- boot[[as.character(summ$Time_Bin[summ$Display_Bin == a])]]
  db <- boot[[as.character(summ$Time_Bin[summ$Display_Bin == bb])]]
  diff <- da - db
  p <- min(1, 2 * min((sum(diff <= 0) + 1) / (n_boot + 1), (sum(diff >= 0) + 1) / (n_boot + 1)))
  data.frame(Bin_A = a, Bin_B = bb, Successive = (bb - a) == 1,
             diff_mean = mean(diff),
             diff_lower95 = unname(quantile(diff, 0.025)), diff_upper95 = unname(quantile(diff, 0.975)),
             p = p)
})) %>%
  mutate(p_holm = p.adjust(p, method = "holm"),
         Significant = p_holm < alpha,
         Result = ifelse(Significant,
                         ifelse(diff_mean > 0, "A higher than B", "A lower than B"),
                         "no significant difference"))

print(as.data.frame(tests %>% transmute(Pair = sprintf("Bin %d vs Bin %d", Bin_A, Bin_B),
                                        Difference = sprintf("%+.3f", diff_mean),
                                        CI95 = sprintf("%+.3f to %+.3f", diff_lower95, diff_upper95),
                                        p = fmt_p(p), p_Holm = fmt_p(p_holm), Result)), row.names = FALSE)
cat("  Difference = beta_SIM(A) - beta_SIM(B); A is the older bin.\n")

# =============================================================================
# 4. FIGURES
# =============================================================================

cat("\n=== 4. FIGURES ===\n")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, w = fig_w, h = fig_h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), p, width = w, height = h, dpi = fig_dpi, bg = page_fill)
  ggsave(file.path(output_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = page_fill,
         device = if (capabilities("cairo")) cairo_pdf else pdf)
  cat("  ", file.path(output_dir, paste0(name, ".png")), " (+ .pdf)\n", sep = "")
}

# 4a. Trend with bootstrap CIs and tests between successive bins -------------
x_max <- max(time_bins$Older)
bands <- time_bins %>% filter(Display_Bin %% 2 == 1)
bin_labels <- summ %>% transmute(Mid, lab = paste0(sprintf("Bin %d\nn = %d", Display_Bin, n_target),
                                                   ifelse(n_provinces < length(regions),
                                                          sprintf("\n(%d provinces)", n_provinces), "")))
succ <- tests %>% filter(Successive) %>%
  left_join(summ %>% select(Bin_A = Display_Bin, xa = Mid, ua = upper95), by = "Bin_A") %>%
  left_join(summ %>% select(Bin_B = Display_Bin, xb = Mid, ub = upper95), by = "Bin_B") %>%
  arrange(Bin_A) %>%
  mutate(y = ifelse(row_number() %% 2 == 1, 1.06, 1.15),       # alternate heights so brackets do not join
         xa = xa - 0.03 * sign(xa - xb), xb = xb + 0.03 * sign(xa - xb),
         lab = ifelse(Significant, stars(p_holm), "n.s."))

p1 <- ggplot(summ, aes(Mid, beta_SIM_mean)) +
  geom_rect(data = bands, aes(xmin = Younger, xmax = Older, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = band_fill, colour = NA) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.07, colour = ink, linewidth = 0.55) +
  geom_line(colour = ink, linewidth = 0.6) +
  geom_point(shape = 21, fill = point_fill, colour = ink, size = 3.6, stroke = 0.8) +
  geom_segment(data = succ, aes(x = xa, xend = xb, y = y, yend = y), inherit.aes = FALSE,
               colour = ink, linewidth = 0.45) +
  geom_segment(data = succ, aes(x = xa, xend = xa, y = y, yend = y - 0.02), inherit.aes = FALSE,
               colour = ink, linewidth = 0.45) +
  geom_segment(data = succ, aes(x = xb, xend = xb, y = y, yend = y - 0.02), inherit.aes = FALSE,
               colour = ink, linewidth = 0.45) +
  geom_text(data = succ, aes((xa + xb) / 2, y + 0.035, label = lab), inherit.aes = FALSE,
            size = base_size / 3.9, colour = ink, fontface = ifelse(succ$Significant, "bold", "plain")) +
  geom_text(data = bin_labels, aes(Mid, 0.02, label = lab), vjust = 0, inherit.aes = FALSE,
            size = base_size / 4.3, colour = ink_soft, lineheight = 0.9) +
  scale_x_reverse(breaks = c(x_max, seq(floor(x_max), 0, by = -1)),
                  labels = function(v) sub("\\.?0+$", "", sprintf("%.2f", v))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(v) sprintf("%.2f", v)) +
  coord_cartesian(xlim = c(x_max, 0), ylim = c(0, 1.2), clip = "off") +
  labs(x = "Age (Ma)", y = expression(beta[SIM]),
       caption = paste0(sprintf("Points = mean of %d bootstrap iterations (sites resampled, then each province cut to %s); ",
                                n_boot, if (common_n) sprintf("n = %d species", n_common) else "the poorest province's n species"),
                        "\nerror bars = 95% CI. Brackets: successive bins; * p < 0.05, ** p < 0.01, *** p < 0.001 ",
                        "(Holm-corrected over all pairs);\nn.s. = not significant.")) +
  theme_classic(base_size = base_size) +
  theme(axis.line = element_line(colour = ink, linewidth = 0.4),
        axis.ticks = element_line(colour = ink, linewidth = 0.4),
        axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
        axis.title.y = element_text(size = base_size * 1.15),
        plot.caption = element_text(colour = ink_soft, size = base_size * 0.62, hjust = 0),
        plot.margin = margin(10, 15, 8, 10), plot.background = element_rect(fill = page_fill, colour = NA))
save_fig(p1, "Fig_bin_comparison_trend")

# 4b. Matrix of Holm-corrected p for all pairs --------------------------------
bin_names <- sprintf("Bin %d\n%s-%s Ma", summ$Display_Bin, fmt_age(summ$Older), fmt_age(summ$Younger))
mat <- tests %>%
  mutate(row = factor(bin_names[match(Bin_B, summ$Display_Bin)], levels = rev(bin_names)),
         col = factor(bin_names[match(Bin_A, summ$Display_Bin)], levels = bin_names),
         p_class = cut(p_holm, c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("< 0.001", "0.001-0.01", "0.01-0.05", ">= 0.05 (n.s.)")),
         lab = sprintf("%+.2f\np = %s", diff_mean, fmt_p(p_holm)))

p2 <- ggplot(mat, aes(col, row)) +
  geom_tile(aes(fill = p_class), colour = page_fill, linewidth = 1.2) +
  geom_text(aes(label = lab, colour = p_class %in% c("< 0.001", "0.001-0.01")),
            size = base_size / 4.4, lineheight = 0.95) +
  scale_fill_manual(values = c("< 0.001" = "#08306B", "0.001-0.01" = "#2171B5",
                               "0.01-0.05" = "#9ECAE1", ">= 0.05 (n.s.)" = "#FBE7A8"),
                    drop = FALSE, name = "Holm-corrected p") +
  scale_colour_manual(values = c(`TRUE` = "#FDF3C4", `FALSE` = ink), guide = "none") +
  labs(x = "Older bin (A)", y = "Younger bin (B)",
       caption = "Cell text: difference in beta_SIM (A - B) and Holm-corrected p. Positive = older bin more provincial.") +
  theme_minimal(base_size = base_size) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = ink, size = base_size * 0.75),
        axis.title = element_text(colour = ink), legend.position = "right",
        plot.caption = element_text(colour = ink_soft, size = base_size * 0.62, hjust = 0),
        plot.background = element_rect(fill = page_fill, colour = NA))
save_fig(p2, "Fig_bin_comparison_matrix", w = 8.5, h = 6.2)

# =============================================================================
# 5. TABLES
# =============================================================================

write.csv(summ %>% transmute(Bin = Display_Bin, Older_Ma = Older, Younger_Ma = Younger,
                             n_provinces, sites_per_province, species_per_province, n_species_used = n_target,
                             beta_SIM_mean = round(beta_SIM_mean, 4), beta_SIM_sd = round(beta_SIM_sd, 4),
                             lower95 = round(lower95, 4), upper95 = round(upper95, 4)),
          file.path(output_dir, "bin_bootstrap_summary.csv"), row.names = FALSE)
write.csv(tests %>% mutate(across(c(diff_mean, diff_lower95, diff_upper95), ~ round(.x, 4)),
                           across(c(p, p_holm), ~ signif(.x, 3))),
          file.path(output_dir, "bin_pairwise_tests.csv"), row.names = FALSE)
cat("  ", file.path(output_dir, "bin_bootstrap_summary.csv"), "\n", sep = "")
cat("  ", file.path(output_dir, "bin_pairwise_tests.csv"), "\n", sep = "")
cat("\n=== STEP 7 ADDENDUM 2 COMPLETE ===\n")
