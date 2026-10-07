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
#
# Display: the first time bin of Step 3c (4.00-3.25 Ma; large mammals in one
# province only) is left out and the remaining bins are renumbered 1, 2, ...
# (Step 3c Bin 2 = Bin 1 here); the age axis starts at 3.25 Ma.
#
# Outputs (Outputs/7_rowan_figures/):
#   Fig2_equal_species.png/.pdf       beta_SIM per bin with error bars
#   beta_sim_equal_species.csv        per bin: n species used, mean, SD, 95% range
#
# Input: Outputs/7_rowan_figures/regional_pa_matrices.rds (Step 7),
#        Outputs/3d_resolved/time_bins.csv.
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

pa_file    <- file.path("Outputs", "7_rowan_figures", "regional_pa_matrices.rds")
bins_file  <- file.path("Outputs", "3d_resolved", "time_bins.csv")
output_dir <- file.path(work_dir, "Outputs", "7_rowan_figures")

regions      <- c("Basin and Range", "Coastal Plain", "Great Plains")
drop_bins    <- 1          # Step 3c bin(s) left out of this figure
n_resamples  <- 999
random_seed  <- 2024

band_fill <- "#E6EDF5"     # shading of alternate time bins (as Fig. 2)
fig_w <- 9; fig_h <- 5.2; fig_dpi <- 300; base_size <- 14

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

# =============================================================================
# 1. READ STEP 7 MATRICES
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
pa_path <- file.path(work_dir, pa_file)
if (!file.exists(pa_path)) stop("File not found:\n  ", pa_path, "\nRun Step 7 first.")
pa_list <- readRDS(pa_path)
time_bins <- read.csv(file.path(work_dir, bins_file), stringsAsFactors = FALSE) %>%
  transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma), Younger = as.numeric(Younger_Ma)) %>%
  filter(!Time_Bin %in% drop_bins) %>% arrange(desc(Older)) %>%
  mutate(Display_Bin = row_number(), Mid = (Older + Younger) / 2)
cat(sprintf("  %d time bins (Step 3c bin %s left out; the rest renumbered from 1)\n",
            nrow(time_bins), paste(drop_bins, collapse = ", ")))

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
    Display_Bin = time_bins$Display_Bin[k], Step3c_Bin = b,
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
# 3. FIGURE (style of Fig. 2)
# =============================================================================

cat("\n=== 3. FIGURE ===\n")
x_max <- max(time_bins$Older)
bands <- time_bins %>% filter(Display_Bin %% 2 == 1)
bin_labels <- res %>% transmute(Mid, lab = paste0(sprintf("Bin %d\nn = %d", Display_Bin, n_species_used),
                                                  ifelse(n_provinces < length(regions), sprintf("\n(%d provinces)", n_provinces), "")))

p <- ggplot(res, aes(Mid, beta_SIM_mean)) +
  geom_rect(data = bands, aes(xmin = Younger, xmax = Older, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = band_fill, colour = NA) +
  geom_errorbar(aes(ymin = beta_SIM_lower95, ymax = beta_SIM_upper95), width = 0.07, colour = "grey15", linewidth = 0.55) +
  geom_line(colour = "grey15", linewidth = 0.6) +
  geom_point(shape = 21, fill = "white", colour = "grey10", size = 3.6, stroke = 0.8) +
  geom_text(data = bin_labels, aes(Mid, 0.02, label = lab), vjust = 0, inherit.aes = FALSE, size = base_size / 4.3,
            colour = "grey30", lineheight = 0.9) +
  # right-hand guide, as in Fig. 2
  annotate("segment", x = -0.24, xend = -0.24, y = 0.53, yend = 0.98, colour = "grey15", linewidth = 0.5,
           arrow = arrow(length = unit(0.1, "in"), type = "open")) +
  annotate("segment", x = -0.24, xend = -0.24, y = 0.47, yend = 0.02, colour = "grey15", linewidth = 0.5,
           arrow = arrow(length = unit(0.1, "in"), type = "open")) +
  annotate("text", x = -0.34, y = 0.755, label = "Higher \u03b2 (provincialism)", angle = 90,
           size = base_size / 4.1, colour = "grey15") +
  annotate("text", x = -0.34, y = 0.245, label = "Lower \u03b2 (homogenization)", angle = 90,
           size = base_size / 4.1, colour = "grey15") +
  scale_x_reverse(breaks = c(x_max, seq(floor(x_max), 0, by = -1)), labels = function(v) sub("\\.?0+$", "", sprintf("%.2f", v))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(v) sprintf("%.2f", v)) +
  coord_cartesian(xlim = c(x_max, 0), ylim = c(0, 1), clip = "off") +
  labs(x = "Age (Ma)", y = expression(beta[SIM]),
       caption = sprintf(paste0("Each province subsampled to the species count of the poorest province (n).\n",
                                "Points = mean of %d random draws; error bars = 95%% range of the draws."), n_resamples)) +
  theme_classic(base_size = base_size) +
  theme(axis.line = element_line(colour = "grey25", linewidth = 0.4),
        axis.ticks = element_line(colour = "grey25", linewidth = 0.4),
        axis.text = element_text(colour = "grey15"), axis.title = element_text(colour = "grey15"),
        axis.title.y = element_text(size = base_size * 1.15),
        plot.caption = element_text(colour = "grey35", size = base_size * 0.62, hjust = 0),
        plot.margin = margin(10, 40, 8, 10), plot.background = element_rect(fill = "white", colour = NA))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "Fig2_equal_species.png"), p, width = fig_w, height = fig_h, dpi = fig_dpi, bg = "white")
ggsave(file.path(output_dir, "Fig2_equal_species.pdf"), p, width = fig_w, height = fig_h, bg = "white",
       device = if (capabilities("cairo")) cairo_pdf else pdf)
cat("  ", file.path(output_dir, "Fig2_equal_species.png"), " (+ .pdf)\n", sep = "")

write.csv(res %>% mutate(across(starts_with("beta"), ~ round(.x, 4))),
          file.path(output_dir, "beta_sim_equal_species.csv"), row.names = FALSE)
cat("  ", file.path(output_dir, "beta_sim_equal_species.csv"), "\n", sep = "")
cat("\n=== STEP 7 ADDENDUM COMPLETE ===\n")
