# =============================================================================
# STEP 3b: WHICH TIME BINS GIVE THE BEST SPREAD OF SITES?
#
# For three physiographic provinces (Basin and Range, Coastal Plain, Great
# Plains; assigned to every site in Step 3), this script keeps the MIDPOINT
# age of every site (from Step 3c) and tests ways of cutting 4.00-0.0117 Ma
# into time bins:
#
#   - 0.75-Myr bins (current)   the Step 3c bins, read from time_bins.csv
#   - Geological stages         Zanclean (from 4 Ma), Piacenzian, Gelasian,
#                               Calabrian, Chibanian, Late Pleistocene
#   - Merged stages             every way of joining neighbouring stages
#                               (e.g. Zanc+Piac / Gela+Cala / Chib+Late)
#   - Epochs                    Pliocene / Early / Middle / Late Pleistocene
#   - Magnetic chrons           Gilbert / Gauss / Matuyama / Brunhes
#   - NALMA                     Blancan / Irvingtonian / Rancholabrean
#                               (boundaries 1.35 and 0.21 Ma, as in Step 1)
#   - Equal width               standard intervals of 250, 400, 500, 600, 750 and
#                               800 kyr and 1, 1.25 and 1.5 Myr, either aligned to
#                               round ages (e.g. 4.0, 3.5, 3.0 ... Ma, as the
#                               500-kyr bins of Rowan et al. 2024) or counted
#                               from 4 Ma
#   - Equal count               boundaries at quantiles of the site midpoints,
#                               so every bin holds about the same number of sites
#   - Optimized                 boundaries chosen to make the SMALLEST
#                               province x bin count as large as possible
#                               (exact dynamic programming; ties broken
#                               towards evenly filled bins)
#
# For every scheme it reports: number of bins, the smallest number of sites in
# any province x bin, how many province x bin cells have fewer than
# 'min_sites_ok' sites (too few for Step 7), and how evenly sites are spread
# (evenness 0-1, 1 = all bins equal). The recommended scheme is the one with
# at least 'min_bins' bins that has the fewest thin cells, then the largest
# minimum, then the highest evenness. The best scheme built from standard
# intervals only (no equal-count or optimized bins) is reported as well; set
# 'recommend_standard_only' to TRUE to recommend it.
#
# Note: many sites share the same midpoint (e.g. sites dated only to a NALMA
# or a magnetochron). A boundary cannot split sites with the same midpoint,
# so some imbalance can remain whatever the scheme.
#
# Counting unit: distinct sites (database + site name) per province and bin,
# as in the Step 7 sites-per-stage figure.
#
# Inputs:  Outputs/3_stages/site_index.csv (Step 3c; holds each site's
#          midpoint age and its Physio_Province from Step 3)
#          Outputs/3_stages/time_bins.csv (Step 3c; the current bins)
#          Outputs/6_matrices/all_records_final.csv (Step 6; only used when
#          only_sites_with_species = TRUE)
# Outputs (Outputs/3b_bin_explorer/):
#   midpoint_distribution.png/.pdf   where the site midpoints fall, per province
#   scheme_comparison.png/.pdf       sites per province x bin for each scheme
#   best_scheme_sites.png/.pdf       sites per bin for the recommended scheme
#   scheme_summary.csv               the comparison table
#   scheme_bins.csv                  bin boundaries of every scheme
#   site_bin_assignments.csv         the bin of every site under every scheme
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Steps 3 and 3c.
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

site_index_file <- file.path("Outputs", "3_stages", "site_index.csv")
time_bins_file  <- file.path("Outputs", "3_stages", "time_bins.csv")
records_file    <- file.path("Outputs", "6_matrices", "all_records_final.csv")
output_dir      <- file.path(work_dir, "Outputs", "3b_bin_explorer")

# Physiographic provinces (names as in the Physio_Province column of Step 3).
provinces <- c("Basin and Range", "Coastal Plain", "Great Plains")
province_colours <- c("Basin and Range" = "#eb6834", "Coastal Plain" = "#4a3aa7",
                      "Great Plains" = "#1baf7a")

only_sites_with_species <- FALSE  # TRUE: count only sites that have species in Step 6 (as Step 7)

age_old   <- 4.000                # study interval (Ma), as the Step 3c bins
age_young <- 0.0117
stage_bounds <- c(4.000, 3.600, 2.580, 1.800, 0.7741, 0.129, 0.0117)
stage_names  <- c("Zanclean", "Piacenzian", "Gelasian", "Calabrian", "Chibanian", "Late Pleistocene")

equal_widths   <- c(0.25, 0.4, 0.5, 0.6, 0.75, 0.8, 1.0, 1.25, 1.5)   # Myr (0.5 = 500 kyr)
width_anchors  <- c("round", "start")   # "round": boundaries at round ages (0.5, 1.0, 1.5 ... Ma,
                                        #   as in Rowan et al. 2024); "start": counted from 4 Ma
recommend_standard_only <- FALSE        # TRUE: recommend only standard intervals (current bins, stages, merged
                                        #   stages, epochs, chrons, NALMA, equal-width bins)
bin_numbers    <- 3:6                # bins tried for equal-count and optimized schemes
min_bin_width  <- 0.3                # Myr; optimized and equal-count bins are at least this long
min_sites_ok   <- 3                  # a province x bin needs at least this many sites (Step 7)
min_bins       <- 4                  # the recommended scheme must have at least this many bins

# Figure look (slide-ready).
base_size <- 18
slide_w   <- 13.33
slide_h   <- 7.5
fig_dpi   <- 300
ink       <- "#1f1f1f"
ink_soft  <- "#5a5a5a"
grid_col  <- "#e6e6e3"

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
save_with_fallback <- function(writer, path) {
  ok <- tryCatch({ writer(path); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
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
theme_slide <- function() {
  theme_minimal(base_size = base_size) +
    theme(text = element_text(colour = ink),
          plot.title = element_text(face = "bold", size = base_size * 1.25, margin = margin(b = 4)),
          plot.subtitle = element_text(colour = ink_soft, size = base_size * 0.8, margin = margin(b = 12)),
          plot.caption = element_text(colour = ink_soft, size = base_size * 0.6, hjust = 0, lineheight = 1.1),
          plot.title.position = "plot", plot.caption.position = "plot",
          axis.title = element_text(colour = ink_soft, size = base_size * 0.85),
          axis.text = element_text(colour = ink, size = base_size * 0.7),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.4),
          strip.text = element_text(face = "bold", size = base_size * 0.85, hjust = 0),
          plot.margin = margin(18, 24, 14, 18), legend.position = "none")
}
fmt_ma <- function(x) formatC(x, format = "fg", digits = 3)

# =============================================================================
# 1. SITES, MIDPOINTS AND PROVINCES
# =============================================================================

cat("=== 1. READING SITES ===\n")
site_index <- read_input(site_index_file)
time_bins  <- read_input(time_bins_file)
if (!"Physio_Province" %in% names(site_index)) {
  stop("site_index.csv has no Physio_Province column: it was written before the spatial step.\n",
       "  File date: ", format(file.mtime(file.path(work_dir, site_index_file)), "%Y-%m-%d %H:%M"), "\n",
       "  Run step3_spatial_binning.R (Step 3), then step3c_time_binning.R (Step 3c), then this script.\n",
       "  (Do not use the old step3_stage_binning.R - delete it from your folder.)")
}
current_bounds <- sort(unique(as.numeric(c(time_bins$Older_Ma, time_bins$Younger_Ma))), decreasing = TRUE)

found <- sort(unique(na.omit(trimws(site_index$Physio_Province))))
missing <- provinces[!tolower(provinces) %in% tolower(found)]
if (length(missing)) {
  stop("No sites in: ", paste(missing, collapse = ", "),
       "
Provinces in site_index.csv: ", paste(found, collapse = "; "))
}

sites <- site_index %>%
  transmute(Database, SiteName = trimws(SiteName), Site_Key, Stage_Number,
            Province = provinces[match(tolower(trimws(Physio_Province)), tolower(provinces))],
            Midpoint_Ma = suppressWarnings(as.numeric(Midpoint_Ma))) %>%
  filter(!is.na(Province), !is.na(Midpoint_Ma), Midpoint_Ma <= age_old, Midpoint_Ma >= age_young) %>%
  distinct(Site_Key, Stage_Number, .keep_all = TRUE)

if (only_sites_with_species && file.exists(file.path(work_dir, records_file))) {
  rec_sites <- read_input(records_file) %>% distinct(Database, SiteName = trimws(SiteName))
  sites <- semi_join(sites, rec_sites, by = c("Database", "SiteName"))
}
sites <- mutate(sites, Site = paste(Database, SiteName, sep = " | "),
                Province = factor(Province, levels = provinces))
cat(sprintf("  %d site records (%d distinct sites) with a midpoint in the %d provinces\n",
            nrow(sites), n_distinct(sites$Site), length(provinces)))
print(as.data.frame(count(distinct(sites, Province, Site), Province, name = "sites")), row.names = FALSE)
cat(sprintf("  %d distinct midpoint ages\n", n_distinct(sites$Midpoint_Ma)))

# Sites per province for a set of boundaries (old -> young).
count_bins <- function(bounds) {
  b <- sort(bounds, decreasing = TRUE)
  # a midpoint equal to a boundary goes to the younger bin, as in Step 3c
  bin <- findInterval(-sites$Midpoint_Ma, -b, left.open = FALSE, rightmost.closed = TRUE)
  bin[sites$Midpoint_Ma == age_young] <- length(b) - 1
  expand_grid(Province = factor(provinces, levels = provinces), Bin = seq_len(length(b) - 1)) %>%
    left_join(data.frame(Province = sites$Province, Bin = bin, Site = sites$Site) %>%
                distinct() %>% count(Province, Bin, name = "n_sites"),
              by = c("Province", "Bin")) %>%
    mutate(n_sites = coalesce(n_sites, 0L), Older = b[Bin], Younger = b[Bin + 1])
}

# =============================================================================
# 2. CANDIDATE SCHEMES
# =============================================================================

cat("\n=== 2. BUILDING SCHEMES ===\n")

schemes <- list(); family <- c()
add_scheme <- function(name, bounds, fam) {
  schemes[[name]] <<- sort(unique(round(bounds, 6)), decreasing = TRUE)
  family[name] <<- fam
}

# The current bins and standard geological schemes.
current_name <- sprintf("%g-Myr bins (current, Step 3c)", round(median(-diff(current_bounds)), 3))
add_scheme(current_name, current_bounds, "Current bins")
add_scheme("Geological stages", stage_bounds, "Stages")
add_scheme("Epochs: Pliocene / Early / Middle / Late Pleistocene",
           c(age_old, 2.58, 0.7741, 0.129, age_young), "Epochs")
add_scheme("Magnetic chrons: Gilbert / Gauss / Matuyama / Brunhes",
           c(age_old, 3.596, 2.581, 0.773, age_young), "Magnetic chrons")
add_scheme("NALMA (Blancan / Irvingtonian / Rancholabrean)", c(age_old, 1.35, 0.21, age_young), "NALMA")

# Every way of merging neighbouring stages (keeps the stage names meaningful).
abbr <- substr(stage_names, 1, 4)
inner <- stage_bounds[-c(1, length(stage_bounds))]
n_in <- length(inner)
for (mask in 0:(2^n_in - 2)) {                   # 2^n - 1 = keep all = the stages themselves
  keep <- as.logical(bitwAnd(mask, 2^(0:(n_in - 1))))
  b <- c(age_old, inner[keep], age_young)
  grp <- cumsum(c(TRUE, keep))
  lab <- paste(tapply(abbr, grp, paste, collapse = "+"), collapse = " / ")
  add_scheme(paste("Merged stages:", lab), b, "Merged stages")
}

# Equal-width bins: aligned to round ages, or counted from age_old. An edge bin
# shorter than half the width is merged into its neighbour.
merge_short_edges <- function(b, w) {
  if (length(b) > 3 && b[1] - b[2] < w / 2) b <- b[-2]
  n <- length(b)
  if (n > 3 && b[n - 1] - b[n] < w / 2) b <- b[-(n - 1)]
  b
}
for (w in equal_widths) {
  wl <- if (w < 1) sprintf("%g kyr", w * 1000) else sprintf("%g Myr", w)
  if ("round" %in% width_anchors) {
    inner_w <- round(seq(floor(age_old / w) * w, ceiling(age_young / w) * w, by = -w), 6)
    b <- c(age_old, inner_w[inner_w < age_old & inner_w > age_young], age_young)
    add_scheme(sprintf("%s bins (round ages)", wl), merge_short_edges(b, w), "Equal width")
  }
  if ("start" %in% width_anchors) {
    b <- seq(age_old, age_young, by = -w)
    if (b[length(b)] > age_young) b <- c(b, age_young)
    add_scheme(sprintf("%s bins (from %g Ma)", wl, age_old), merge_short_edges(b, w), "Equal width")
  }
}

# Possible cut points: halfway between neighbouring distinct midpoints.
u <- sort(unique(sites$Midpoint_Ma), decreasing = TRUE)
cuts <- round((u[-1] + u[-length(u)]) / 2, 4)
cuts <- cuts[cuts < age_old & cuts > age_young]
if (length(cuts) > 250) {                         # thin a very dense set of cut points
  cuts <- unique(round(cuts / 0.02) * 0.02)
}
cand <- c(age_old, sort(cuts, decreasing = TRUE), age_young)
cat(sprintf("  %d possible boundary positions\n", length(cand) - 2))

# Equal count: quantiles of the pooled site midpoints, snapped to cut points.
for (k in bin_numbers) {
  q <- quantile(sites$Midpoint_Ma, probs = rev(seq_len(k - 1) / k), names = FALSE)
  b <- c(age_old, vapply(q, function(x) cand[which.min(abs(cand - x))], numeric(1)), age_young)
  b <- unique(b)
  add_scheme(sprintf("Equal count, %d bins", k), b, "Equal count")
}

# Optimized: dynamic programming over cut points, maximizing the smallest
# province x bin count (ties: larger sum of log counts = more even).
n_c <- length(cand)
pos <- vapply(sites$Midpoint_Ma, function(m) max(which(cand >= m)), integer(1))   # interval index
pos[pos >= n_c] <- n_c - 1
cell <- array(NA_real_, c(n_c, n_c, 2))           # [start, end, (min, sumlog)]
for (a in 1:(n_c - 1)) {
  for (e in (a + 1):n_c) {
    if (cand[a] - cand[e] < min_bin_width - 1e-9) next
    inside <- pos >= a & pos < e
    cnt <- vapply(provinces, function(p) length(unique(sites$Site[inside & sites$Province == p])), numeric(1))
    cell[a, e, ] <- c(min(cnt), sum(log1p(cnt)))
  }
}
best_dp <- function(k) {
  fmin <- matrix(-Inf, k + 1, n_c); fsum <- matrix(-Inf, k + 1, n_c); from <- matrix(NA, k + 1, n_c)
  fmin[1, 1] <- Inf; fsum[1, 1] <- 0
  for (j in 1:k) for (e in 2:n_c) for (a in 1:(e - 1)) {
    if (!is.finite(fsum[j, a]) || is.na(cell[a, e, 1])) next
    m <- min(fmin[j, a], cell[a, e, 1]); s <- fsum[j, a] + cell[a, e, 2]
    if (m > fmin[j + 1, e] || (m == fmin[j + 1, e] && s > fsum[j + 1, e])) {
      fmin[j + 1, e] <- m; fsum[j + 1, e] <- s; from[j + 1, e] <- a
    }
  }
  if (!is.finite(fsum[k + 1, n_c])) return(NULL)
  path <- n_c; e <- n_c
  for (j in (k + 1):2) { e <- from[j, e]; path <- c(e, path) }
  cand[path]
}
for (k in bin_numbers) {
  b <- best_dp(k)
  if (!is.null(b)) add_scheme(sprintf("Optimized, %d bins", k), b, "Optimized")
}

# =============================================================================
# 3. COMPARE SCHEMES
# =============================================================================

cat("\n=== 3. COMPARING SCHEMES ===\n")

# Pielou evenness of the site counts over all bins (empty bins lower it).
evenness <- function(x) {
  n <- length(x)
  if (n < 2 || sum(x) == 0) return(0)
  p <- x / sum(x); p <- p[p > 0]
  -sum(p * log(p)) / log(n)
}

counts <- bind_rows(lapply(names(schemes), function(s) mutate(count_bins(schemes[[s]]), Scheme = s)))
summary_tab <- counts %>%
  group_by(Scheme) %>%
  summarise(n_bins = n_distinct(Bin),
            boundaries_Ma = paste(fmt_ma(sort(unique(c(Older, Younger)), decreasing = TRUE)), collapse = " | "),
            min_sites = min(n_sites),
            cells_below_ok = sum(n_sites < min_sites_ok),
            empty_cells = sum(n_sites == 0),
            max_sites = max(n_sites),
            evenness = mean(tapply(n_sites, Province, evenness)),
            .groups = "drop") %>%
  mutate(evenness = round(evenness, 3),
         Family = unname(family[Scheme]),
         Standard = !Family %in% c("Equal count", "Optimized"),
         Eligible = n_bins >= min_bins) %>%
  arrange(desc(Eligible), cells_below_ok, desc(min_sites), desc(evenness)) %>%
  relocate(Family, .after = Scheme)

best_standard <- summary_tab$Scheme[summary_tab$Eligible & summary_tab$Standard][1]
best_overall  <- summary_tab$Scheme[summary_tab$Eligible][1]
best <- if (recommend_standard_only) best_standard else best_overall
summary_tab <- mutate(summary_tab, Recommended = Scheme == best)
save_csv(summary_tab, "scheme_summary")

cat(sprintf("  %d schemes tested. Top 15 (ranked: >= %d bins, fewest cells with < %d sites,\n",
            nrow(summary_tab), min_bins, min_sites_ok))
cat("  then largest minimum, then evenness; full table in scheme_summary.csv):\n")
print(as.data.frame(head(select(summary_tab, Scheme, n_bins, min_sites, cells_below_ok,
                                empty_cells, max_sites, evenness), 15)), row.names = FALSE)
cat("\n  Best scheme in each family:\n")
print(as.data.frame(summary_tab %>% filter(Eligible) %>% group_by(Family) %>% slice(1) %>% ungroup() %>%
        select(Family, Scheme, n_bins, min_sites, cells_below_ok, evenness)), row.names = FALSE)
for (nm in unique(c(best_overall, best_standard))) {
  cat(sprintf("\n  %s: %s\n    Boundaries (Ma): %s\n",
              if (identical(nm, best_overall) && identical(nm, best_standard)) "Best overall and best standard scheme"
              else if (identical(nm, best_overall)) "Best overall" else "Best standard-interval scheme",
              nm, summary_tab$boundaries_Ma[summary_tab$Scheme == nm]))
}
cat(sprintf("\n  Recommended (%s): %s\n",
            if (recommend_standard_only) "standard intervals only" else "all schemes", best))

save_csv(counts %>% distinct(Scheme, Bin, Older, Younger) %>% arrange(Scheme, Bin), "scheme_bins")

assign_tab <- sites %>% select(Database, SiteName, Site_Key, Province, Midpoint_Ma)
for (s in names(schemes)) {
  b <- sort(schemes[[s]], decreasing = TRUE)
  bin <- findInterval(-sites$Midpoint_Ma, -b, rightmost.closed = TRUE)
  bin[sites$Midpoint_Ma == age_young] <- length(b) - 1
  assign_tab[[s]] <- sprintf("Bin %d (%s-%s Ma)", bin, fmt_ma(b[bin]), fmt_ma(b[bin + 1]))
}
save_csv(assign_tab, "site_bin_assignments")

# =============================================================================
# 4. FIGURES
# =============================================================================

cat("\n=== 4. FIGURES ===\n")

# 4a. Where do the site midpoints fall?
hist_w <- 0.05
p_hist <- ggplot(distinct(sites, Province, Site, Midpoint_Ma), aes(Midpoint_Ma, fill = Province)) +
  geom_vline(xintercept = current_bounds[-c(1, length(current_bounds))], colour = ink_soft,
             linetype = "dashed", linewidth = 0.5) +
  geom_histogram(binwidth = hist_w, boundary = age_young, colour = "white", linewidth = 0.2) +
  geom_text(data = data.frame(x = (current_bounds[-1] + current_bounds[-length(current_bounds)]) / 2,
                              lab = seq_len(length(current_bounds) - 1)),
            aes(x = x, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.4, size = base_size / 4.6, colour = ink_soft) +
  facet_wrap(~ Province, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = province_colours) +
  scale_x_reverse(breaks = 4:0, expand = c(0.01, 0)) +
  coord_cartesian(xlim = c(age_old, age_young)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3)), breaks = scales::pretty_breaks(3)) +
  labs(title = "Where do the site ages fall?",
       subtitle = sprintf("Number of sites per %.2f Myr of midpoint age; dashed lines and numbers = current Step 3c bins", hist_w),
       x = "Midpoint age (Ma)", y = "Sites",
       caption = "Tall single spikes = many sites with the same midpoint (e.g. dated only to a NALMA); a bin boundary cannot split them.") +
  theme_slide()
save_fig(p_hist, "midpoint_distribution", height = slide_h * 1.2)

# 4b. Scheme comparison (current, NALMA, equal width 1 Myr, best of each family, recommended).
first_of <- function(fam) summary_tab$Scheme[summary_tab$Family == fam & summary_tab$Eligible][1]
show <- unique(c(current_name,
                 "Geological stages",
                 "Magnetic chrons: Gilbert / Gauss / Matuyama / Brunhes",
                 first_of("Merged stages"),
                 if ("500 kyr bins (round ages)" %in% names(schemes)) "500 kyr bins (round ages)",
                 first_of("Equal width"),
                 first_of("Optimized"),
                 best))
show <- show[!is.na(show)]
wrap_lab <- function(x) gsub("(.{1,34})(\\s|$)", "\\1\n", x) %>% sub("\n$", "", .)
lab_of <- setNames(ifelse(show == best, paste0(wrap_lab(show), "\n(recommended)"), wrap_lab(show)), show)
cmp <- counts %>% filter(Scheme %in% show) %>%
  mutate(Scheme = factor(lab_of[Scheme], levels = lab_of),
         txt_col = ifelse(n_sites >= 40, "white", ink),
         flag = n_sites < min_sites_ok)
p_cmp <- ggplot(cmp) +
  geom_rect(aes(xmin = Younger, xmax = Older, ymin = 0, ymax = 1, fill = n_sites), colour = "white", linewidth = 1) +
  geom_rect(data = filter(cmp, flag), aes(xmin = Younger, xmax = Older, ymin = 0, ymax = 1),
            fill = NA, colour = "#e34948", linewidth = 0.9) +
  geom_text(aes(x = (Older + Younger) / 2, y = 0.5, label = n_sites, colour = txt_col,
                size = ifelse(Older - Younger < 0.6, base_size / 6.2, base_size / 4.6),
                angle = ifelse(Older - Younger < 0.3, 90, 0)),   # narrow cells: label turned sideways
            fontface = "bold") +
  scale_size_identity() +
  scale_colour_identity() +
  facet_grid(Scheme ~ Province, switch = "y") +
  scale_fill_gradientn(colours = c("#f4f4f2", "#cde2fb", "#6da7ec", "#256abf", "#0d366b"),
                       limits = c(0, NA), name = "Sites") +
  scale_x_reverse(limits = c(age_old, age_young), breaks = c(4, 3, 2, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
  labs(title = "How evenly do different time bins spread the sites?",
       subtitle = "Sites per province and bin for each binning scheme (same midpoint ages throughout)",
       x = "Age (Ma)", y = NULL,
       caption = sprintf("Red outline = fewer than %d sites (too few for the Step 7 analyses). Full comparison of all schemes in scheme_summary.csv.",
                         min_sites_ok)) +
  theme_slide() +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = base_size * 0.62, face = "plain"),
        strip.text.x = element_text(size = base_size * 0.8),
        strip.placement = "outside", panel.grid = element_blank(),
        panel.spacing.y = unit(0.35, "lines"), panel.spacing.x = unit(0.8, "lines"),
        legend.position = "right", axis.text.x = element_text(size = base_size * 0.6))
save_fig(p_cmp, "scheme_comparison", height = slide_h * max(1.15, 0.2 * length(show) + 0.2))

# 4c. The recommended scheme, as in the Step 7 sites-per-stage figure.
bb <- counts %>% filter(Scheme == best) %>%
  mutate(Label = factor(sprintf("%s-%s", fmt_ma(Older), fmt_ma(Younger)),
                        levels = unique(sprintf("%s-%s", fmt_ma(Older), fmt_ma(Younger)))))
p_best <- ggplot(bb, aes(Label, n_sites, fill = Province)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n_sites), vjust = -0.4, size = base_size / 3.8, fontface = "bold", colour = ink) +
  facet_wrap(~ Province, ncol = length(provinces)) +
  scale_fill_manual(values = province_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), breaks = scales::pretty_breaks(4)) +
  labs(title = paste("Sites per time bin -", best),
       subtitle = paste("Bin boundaries (Ma):", summary_tab$boundaries_Ma[summary_tab$Scheme == best]),
       x = "Time bin (Ma)", y = "Number of sites",
       caption = "Midpoint ages as in Step 3c; a site whose midpoint falls on a boundary goes to the younger bin.") +
  theme_slide() +
  theme(axis.text.x = element_text(size = base_size * 0.6))
save_fig(p_best, "best_scheme_sites", height = slide_h * 1.1)
cat("  Figures saved\n")

cat("\n=== STEP 3b COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: sites, schemes, counts, summary_tab\n")
