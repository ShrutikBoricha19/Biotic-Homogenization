# =============================================================================
# STEP 6: PRESENCE-ABSENCE MATRICES AND SPECIES ACCUMULATION CURVES
#
# For the three physiographic provinces - Basin and Range, Coastal Plain and
# Great Plains (Great Plains including the Central Lowland sites west of the
# Mississippi, Step 5b) - and every time bin:
#
#   1. a presence-absence matrix: sites (rows) x species (columns), 1 = the
#      species occurs at the site;
#   2. a species accumulation curve (vegan::specaccum, sites added in random
#      order, 'n_permutations' times; mean +/- 95% band);
#   3. an automatic check of whether the curve reaches an asymptote.
#
# The asymptote check (per province x time bin species pool):
#   Final slope      expected number of new species added by the last site
#                    of the curve (mean of the last increment). A flat curve
#                    has a slope near 0.
#   Clench model     S(n) = a n / (1 + b n) fitted to the curve; a / b is the
#                    estimated total number of species; Completeness = observed
#                    species / that total.
#   Chao2            incidence-based richness estimate (vegan::specpool);
#                    observed / Chao2 is a second completeness value.
#   Coverage         sample coverage for incidence data (Chao & Jost 2012):
#                    the share of the pool's incidences belonging to species
#                    already found.
#   Verdict          "Asymptote reached"      final slope < slope_threshold AND
#                                             Clench completeness >= completeness_threshold
#                    "Approaching asymptote"  one of the two
#                    "No asymptote"           neither
#                    "Too few sites"          fewer than min_sites sites
#   (Thresholds 0.1 and 70% follow Jimenez-Valverde & Hortal 2003.)
#
# Outputs:
#   Outputs/6_presence_absence/
#     <Province>/bin<k>.csv     one matrix per province and time bin
#                               (Site_Key, Site, Latitude, Longitude, species...)
#     presence_absence_matrices.rds   all matrices as an R list [[province]][[bin]]
#     matrix_summary.csv         sites, species and records per matrix
#   Outputs/6_species_accumulation/
#     accumulation_all_pools.png/.pdf      3 provinces x all bins
#     asymptote_verdicts.png/.pdf          verdict per province and bin
#     pools/<Province>_bin<k>.png          one plot per species pool
#     asymptote_diagnostics.csv            the numbers behind the verdicts
#     accumulation_curves.csv              the curves (mean and sd per site count)
#
# Input: Outputs/5b_great_plains/master_data_unique.csv (Step 5b),
#        Outputs/3d_resolved/time_bins.csv
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 5b.
# Needs: dplyr, tidyr, ggplot2, vegan.
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "vegan")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
suppressPackageStartupMessages(library(vegan))

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

master_file <- file.path("Outputs", "5b_great_plains", "master_data_unique.csv")
bins_file   <- file.path("Outputs", "3d_resolved", "time_bins.csv")
matrix_dir  <- file.path(work_dir, "Outputs", "6_presence_absence")
accum_dir   <- file.path(work_dir, "Outputs", "6_species_accumulation")

provinces <- c("Basin and Range", "Coastal Plain", "Great Plains")
province_colours <- c("Basin and Range" = "#EB6834", "Coastal Plain" = "#4A3AA7", "Great Plains" = "#1BAF7A")

n_permutations         <- 999
random_seed            <- 2024
min_sites              <- 5      # fewer sites: "Too few sites"
slope_threshold        <- 0.1    # new species per added site at the end of the curve
completeness_threshold <- 0.70   # observed / Clench asymptote

verdict_colours <- c("Asymptote reached" = "#2E7D32", "Approaching asymptote" = "#E69F00",
                     "No asymptote" = "#C62828", "Too few sites" = "#9E9E9E")

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
write_out <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ok <- tryCatch({ write.csv(df, path, row.names = FALSE, na = "", fileEncoding = "UTF-8"); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) stop("Could not write ", path, "\n  Close it (e.g. in Excel) and run again.")
}
save_fig <- function(p, path_no_ext, width, height) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height, bg = "white",
         device = if (capabilities("cairo")) cairo_pdf else pdf)
}
file_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
fmt_age <- function(x) ifelse(x < 0.1, sprintf("%.4f", x), sprintf("%.2f", x))

# Clench (Michaelis-Menten) model fitted to the mean accumulation curve.
fit_clench <- function(x, s) {
  if (length(x) < 3 || max(s) <= 0) return(NULL)
  # starting values from the linear form x/S = 1/a + (b/a) x
  lin <- tryCatch(coef(lm(I(x / s) ~ x)), error = function(e) c(NA, NA))
  a0 <- if (is.finite(lin[1]) && lin[1] > 0) 1 / lin[1] else s[1]
  b0 <- if (is.finite(lin[2]) && lin[2] > 0) lin[2] * a0 else a0 / max(s)
  fit <- tryCatch(nls(s ~ a * x / (1 + b * x), start = list(a = a0, b = b0),
                      control = nls.control(maxiter = 200, warnOnly = TRUE)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- coef(fit)
  if (any(!is.finite(cf)) || cf["a"] <= 0 || cf["b"] <= 0) return(NULL)
  list(a = unname(cf["a"]), b = unname(cf["b"]))
}

# Sample coverage for incidence data (Chao & Jost 2012, Eq. for incidence).
incidence_coverage <- function(m) {
  t_sites <- nrow(m); freq <- colSums(m); u <- sum(freq)
  q1 <- sum(freq == 1); q2 <- sum(freq == 2)
  if (u == 0 || t_sites < 2) return(NA_real_)
  if (q2 > 0) {
    1 - q1 / u * ((t_sites - 1) * q1 / ((t_sites - 1) * q1 + 2 * q2))
  } else {
    1 - q1 / u * ((t_sites - 1) * (q1 - 1) / ((t_sites - 1) * (q1 - 1) + 2))
  }
}

# =============================================================================
# 1. READ THE MASTER DATA
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
md <- read_text_csv(master_file)
time_bins <- read_text_csv(bins_file) %>%
  transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma), Younger = as.numeric(Younger_Ma)) %>%
  arrange(Time_Bin) %>%
  mutate(Bin_Text = sprintf("Bin %d\n%s-%s Ma", Time_Bin, fmt_age(Older), fmt_age(Younger)))

missing <- setdiff(provinces, unique(md$Spatial_Bin))
if (length(missing)) {
  stop("No records for: ", paste(missing, collapse = ", "), "\nProvinces in the data: ",
       paste(sort(unique(na.omit(md$Spatial_Bin))), collapse = "; "))
}
md <- md %>% filter(Spatial_Bin %in% provinces) %>%
  mutate(Time_Bin = as.integer(Time_Bin))
cat(sprintf("  %d species x site x bin records in the three provinces\n", nrow(md)))

# =============================================================================
# 2. PRESENCE-ABSENCE MATRICES
# =============================================================================

cat("\n=== 2. PRESENCE-ABSENCE MATRICES ===\n")
matrices <- list()
matrix_summary <- list()
for (p in provinces) {
  matrices[[p]] <- list()
  for (b in time_bins$Time_Bin) {
    d <- md %>% filter(Spatial_Bin == p, Time_Bin == b)
    sites <- d %>% distinct(Site_Key, .keep_all = TRUE) %>% arrange(Site_Key) %>%
      select(Site_Key, Site, Latitude, Longitude)
    spp <- sort(unique(d$Species))
    m <- matrix(0L, nrow = nrow(sites), ncol = length(spp), dimnames = list(sites$Site_Key, spp))
    if (nrow(d)) m[cbind(match(d$Site_Key, sites$Site_Key), match(d$Species, spp))] <- 1L
    matrices[[p]][[as.character(b)]] <- m
    write_out(cbind(sites, as.data.frame(m, check.names = FALSE)),
              file.path(matrix_dir, file_name(p), sprintf("bin%d.csv", b)))
    matrix_summary[[length(matrix_summary) + 1]] <- data.frame(
      Province = p, Time_Bin = b, n_sites = nrow(m), n_species = ncol(m), n_records = sum(m),
      mean_species_per_site = if (nrow(m)) round(mean(rowSums(m)), 2) else NA)
  }
}
matrix_summary <- bind_rows(matrix_summary)
saveRDS(matrices, file.path(matrix_dir, "presence_absence_matrices.rds"))
write_out(matrix_summary, file.path(matrix_dir, "matrix_summary.csv"))
print(as.data.frame(matrix_summary), row.names = FALSE)
cat("  Matrices saved in", matrix_dir, "\n")

# =============================================================================
# 3. SPECIES ACCUMULATION AND THE ASYMPTOTE CHECK
# =============================================================================

cat(sprintf("\n=== 3. SPECIES ACCUMULATION (%d random orders) ===\n", n_permutations))
set.seed(random_seed)
curves <- list(); fits <- list(); diag <- list()
for (p in provinces) for (b in time_bins$Time_Bin) {
  m <- matrices[[p]][[as.character(b)]]
  m <- m[, colSums(m) > 0, drop = FALSE]
  n <- nrow(m); s_obs <- ncol(m)
  row <- data.frame(Province = p, Time_Bin = b, n_sites = n, n_species = s_obs,
                    final_slope = NA_real_, clench_asymptote = NA_real_, clench_completeness = NA_real_,
                    clench_final_slope = NA_real_, chao2 = NA_real_, chao2_se = NA_real_,
                    chao2_completeness = NA_real_, coverage = NA_real_, uniques = sum(colSums(m) == 1),
                    duplicates = sum(colSums(m) == 2))
  if (n >= 2 && s_obs > 0) {
    # (permute reports when a small pool has fewer orders than n_permutations; not needed here)
    sa <- suppressMessages(specaccum(m, method = "random", permutations = n_permutations))
    cv <- data.frame(Province = p, Time_Bin = b, sites = sa$sites, richness = sa$richness,
                     sd = ifelse(is.na(sa$sd), 0, sa$sd))   # no spread among orders (e.g. 2 sites)
    curves[[length(curves) + 1]] <- cv
    row$final_slope <- sa$richness[n] - sa$richness[n - 1]
    cl <- fit_clench(sa$sites, sa$richness)
    if (!is.null(cl)) {
      row$clench_asymptote <- cl$a / cl$b
      row$clench_completeness <- s_obs / row$clench_asymptote
      row$clench_final_slope <- cl$a / (1 + cl$b * n)^2
      fits[[length(fits) + 1]] <- data.frame(Province = p, Time_Bin = b,
                                             sites = seq(1, n, length.out = 100)) %>%
        mutate(richness = cl$a * sites / (1 + cl$b * sites), asymptote = cl$a / cl$b)
    }
    sp <- tryCatch(specpool(m), error = function(e) NULL)
    if (!is.null(sp) && is.finite(sp$chao) && sp$chao > 0) {
      row$chao2 <- sp$chao; row$chao2_se <- sp$chao.se; row$chao2_completeness <- s_obs / sp$chao
    }
    row$coverage <- incidence_coverage(m)
  }
  diag[[length(diag) + 1]] <- row
}
curves <- bind_rows(curves); fits <- bind_rows(fits)

asymptote_diagnostics <- bind_rows(diag) %>%
  mutate(slope_ok = !is.na(final_slope) & final_slope < slope_threshold,
         completeness_ok = !is.na(clench_completeness) & clench_completeness >= completeness_threshold,
         Verdict = case_when(n_sites < min_sites ~ "Too few sites",
                             slope_ok & completeness_ok ~ "Asymptote reached",
                             slope_ok | completeness_ok ~ "Approaching asymptote",
                             TRUE ~ "No asymptote"),
         Reason = case_when(
           Verdict == "Too few sites" ~ sprintf("only %d site(s); at least %d needed", n_sites, min_sites),
           TRUE ~ paste0(sprintf("final slope %.2f (%s %.2f); ", final_slope,
                                 ifelse(slope_ok, "<", ">="), slope_threshold),
                         ifelse(is.na(clench_completeness), "Clench model could not be fitted",
                                sprintf("%.0f%% of the Clench asymptote (%s %.0f%%)", 100 * clench_completeness,
                                        ifelse(completeness_ok, ">=", "<"), 100 * completeness_threshold))))) %>%
  mutate(across(c(final_slope, clench_asymptote, clench_final_slope, chao2, chao2_se), ~ round(.x, 2)),
         across(c(clench_completeness, chao2_completeness, coverage), ~ round(.x, 3))) %>%
  select(Province, Time_Bin, n_sites, n_species, Verdict, Reason, final_slope, clench_asymptote,
         clench_completeness, clench_final_slope, chao2, chao2_se, chao2_completeness, coverage,
         uniques, duplicates)

write_out(asymptote_diagnostics, file.path(accum_dir, "asymptote_diagnostics.csv"))
write_out(curves %>% mutate(across(c(richness, sd), ~ round(.x, 3))),
          file.path(accum_dir, "accumulation_curves.csv"))
cat("\n  Asymptote check:\n")
print(as.data.frame(asymptote_diagnostics %>%
                      select(Province, Time_Bin, n_sites, n_species, final_slope,
                             clench_completeness, coverage, Verdict)), row.names = FALSE)

# =============================================================================
# 4. FIGURES
# =============================================================================

cat("\n=== 4. FIGURES ===\n")
lab_bin <- setNames(time_bins$Bin_Text, time_bins$Time_Bin)
facet_df <- function(df) df %>% mutate(Province = factor(Province, levels = provinces),
                                       Bin = factor(lab_bin[as.character(Time_Bin)], levels = time_bins$Bin_Text))
tags <- facet_df(asymptote_diagnostics) %>%
  mutate(tag = ifelse(Verdict == "Too few sites", sprintf("Too few sites (%d)", n_sites),
                      sprintf("%s\nslope %.2f | %s", Verdict, final_slope,
                              ifelse(is.na(clench_completeness), "no fit",
                                     sprintf("%.0f%% complete", 100 * clench_completeness)))))
asym_lines <- facet_df(fits %>% distinct(Province, Time_Bin, asymptote) %>%
                         left_join(asymptote_diagnostics %>% select(Province, Time_Bin, n_species),
                                   by = c("Province", "Time_Bin")) %>%
                         filter(asymptote <= 2 * n_species))      # only when it fits on the panel

p_all <- ggplot() +
  # every panel gets data, so pools without a curve (1 site) still draw
  geom_blank(data = facet_df(asymptote_diagnostics), aes(x = pmax(n_sites, 2), y = pmax(n_species, 1))) +
  geom_ribbon(data = facet_df(curves), aes(sites, ymin = pmax(richness - 1.96 * sd, 0),
                                           ymax = richness + 1.96 * sd, fill = Province), alpha = 0.22) +
  geom_line(data = facet_df(curves), aes(sites, richness, colour = Province), linewidth = 0.9) +
  geom_line(data = facet_df(fits), aes(sites, richness), linetype = "22", colour = "grey20", linewidth = 0.5) +
  geom_hline(data = asym_lines, aes(yintercept = asymptote), linetype = "dotted", colour = "grey30") +
  geom_label(data = tags, aes(x = -Inf, y = Inf, label = tag, colour = Verdict), hjust = -0.03, vjust = 1.08,
             size = 2.7, label.size = 0.25, fill = alpha("white", 0.85), lineheight = 0.95) +
  facet_grid(Province ~ Bin, scales = "free", drop = FALSE) +
  scale_colour_manual(values = c(province_colours, verdict_colours), guide = "none") +
  scale_fill_manual(values = province_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(title = "Species accumulation in each regional species pool",
       subtitle = sprintf(paste0("Mean of %d random site orders with 95%% band; dashed = Clench model, dotted = its asymptote. ",
                                 "Asymptote reached: final slope < %.1f and >= %.0f%% of the Clench asymptote."),
                          n_permutations, slope_threshold, 100 * completeness_threshold),
       x = "Number of sites", y = "Number of species") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95"), strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(colour = "grey30", size = 10))
save_fig(p_all, file.path(accum_dir, "accumulation_all_pools"), 16, 9)
cat("  accumulation_all_pools.png (+ .pdf)\n")

# Verdict table as a figure.
vt <- facet_df(asymptote_diagnostics) %>%
  mutate(label = ifelse(Verdict == "Too few sites", sprintf("%d site%s", n_sites, ifelse(n_sites == 1, "", "s")),
                        sprintf("%d sites, %d spp.\n%s complete", n_sites, n_species,
                                ifelse(is.na(clench_completeness), "-", sprintf("%.0f%%", 100 * clench_completeness)))),
         Verdict = factor(Verdict, levels = names(verdict_colours)))
p_vt <- ggplot(vt, aes(Bin, Province, fill = Verdict)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = label), colour = "white", fontface = "bold", size = 3.6, lineheight = 0.95) +
  scale_fill_manual(values = verdict_colours, name = NULL, drop = FALSE) +
  scale_y_discrete(limits = rev(provinces)) +
  scale_x_discrete(position = "top") +
  labs(title = "Has each species pool reached an asymptote?",
       subtitle = sprintf("Final slope < %.1f new species per site and >= %.0f%% of the Clench asymptote; at least %d sites",
                          slope_threshold, 100 * completeness_threshold, min_sites),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.text = element_text(face = "bold", colour = "grey15"),
        plot.title = element_text(face = "bold", size = 16), plot.title.position = "plot",
        plot.subtitle = element_text(colour = "grey30", size = 10.5), legend.position = "bottom")
save_fig(p_vt, file.path(accum_dir, "asymptote_verdicts"), 13.33, 5.5)
cat("  asymptote_verdicts.png (+ .pdf)\n")

# One plot per species pool.
for (i in seq_len(nrow(asymptote_diagnostics))) {
  r <- asymptote_diagnostics[i, ]
  cv <- curves %>% filter(Province == r$Province, Time_Bin == r$Time_Bin)
  if (!nrow(cv)) next
  ft <- fits %>% filter(Province == r$Province, Time_Bin == r$Time_Bin)
  pp <- ggplot(cv, aes(sites, richness)) +
    geom_ribbon(aes(ymin = pmax(richness - 1.96 * sd, 0), ymax = richness + 1.96 * sd),
                fill = province_colours[[r$Province]], alpha = 0.22) +
    geom_line(colour = province_colours[[r$Province]], linewidth = 1.1) +
    geom_point(colour = province_colours[[r$Province]], size = 1.3)
  if (nrow(ft)) {
    pp <- pp + geom_line(data = ft, linetype = "22", colour = "grey20")
    if (ft$asymptote[1] <= 2 * r$n_species) {
      pp <- pp + geom_hline(yintercept = ft$asymptote[1], linetype = "dotted", colour = "grey30") +
        annotate("text", x = 1, y = ft$asymptote[1], label = sprintf("Clench asymptote %.1f", ft$asymptote[1]),
                 hjust = 0, vjust = -0.5, size = 3.5, colour = "grey30")
    }
  }
  pp <- pp +
    labs(title = sprintf("%s - %s)", r$Province, gsub("\n", " (", lab_bin[as.character(r$Time_Bin)])),
         subtitle = sprintf("%s: %s\nChao2 %s | coverage %s | %d sites, %d species",
                            r$Verdict, r$Reason,
                            ifelse(is.na(r$chao2), "-", sprintf("%.1f", r$chao2)),
                            ifelse(is.na(r$coverage), "-", sprintf("%.0f%%", 100 * r$coverage)),
                            r$n_sites, r$n_species),
         x = "Number of sites", y = "Number of species") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = verdict_colours[[r$Verdict]], size = 10.5),
          panel.grid.minor = element_blank())
  dir.create(file.path(accum_dir, "pools"), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(accum_dir, "pools", sprintf("%s_bin%d.png", file_name(r$Province), r$Time_Bin)),
         pp, width = 8, height = 5.5, dpi = 300, bg = "white")
}
cat("  one plot per pool in", file.path(accum_dir, "pools"), "\n")

cat("\n=== STEP 6 COMPLETE ===\n")
cat("Objects in your Environment: matrices, asymptote_diagnostics, curves\n")
