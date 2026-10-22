# =============================================================================
# STEP 7 ADDENDUM: MULTISITE beta_SIM WITH EQUAL AREA AND EQUAL SPECIES
#
# Chapter 1's equal-species subsampling (step7_addendum_equal_species.R) gives
# every province the same number of species, but not the same AREA or the same
# SPREAD of sites. Antell, Benson & Saupe (2024, Paleobiology 50:177-193) show
# that both matter: more sites, or more widely scattered sites, give more
# turnover and a higher beta even when species numbers are equal. This addendum
# standardizes both before the equal-species cut:
#
#   1. Sites -> equal-area grid cells (cell_km wide; Antell et al. 2024 used
#      cells ~100 km across). Each occupied cell is one "site".
#   2. Spatial subsampling within each province and bin (the divvy methods of
#      Antell et al. 2024):
#        clustr  - grow a nearest-neighbour cluster from a random seed cell until
#                  its maximum diameter would exceed clustr_km (suits the long,
#                  narrow Coastal Plain);
#        cookies - all cells within cookie_r_km of a random seed cell;
#      then draw the SAME number of cells (the quota) in every province and bin,
#      so area and extent are equal everywhere.
#   3. Species of the drawn cells pooled per province, then every province cut
#      to the species count of the poorest province (as Chapter 1).
#   4. Multisite beta_SIM (Rowan et al. 2024 Eq. 1), repeated n_iter times.
#
# Quota: by default the smallest number of occupied cells among the province-
# bins that have at least min_cells cells, so the same area is used everywhere.
# Province-bins with fewer than min_cells cells, or with no cluster/circle
# that holds the quota, are reported as INSUFFICIENT and that bin gets no point
# (the quota is not lowered for single bins; Antell et al. 2024).
#
#   Points      mean beta_SIM of the iterations
#   Error bars  95% range of the iterations (a quantile range, not a CI: it mixes
#               spatial, area and species-draw variation; Antell et al. 2024)
#   Brackets    p for the difference between successive bins: the iteration
#               values of the two bins are paired at random and the share of
#               differences on the other side of zero is doubled (two-sided);
#               Holm-corrected over the successive pairs.
#
# Outputs (Outputs/7_rowan_figures/):
#   Fig2_area_species_<method>.png/.pdf          same style as Fig2_equal_species
#   beta_sim_area_species_<method>.csv           values per bin
#   beta_sim_area_species_<method>_tests.csv     successive-bin tests
#   area_species_cells_per_province_bin.csv      occupied cells and sufficiency
#
# Inputs: Outputs/5b_great_plains/master_data_unique.csv (Step 5b),
#         Outputs/7_rowan_figures/species_traits.csv (Step 7),
#         Outputs/3d_resolved/time_bins.csv.
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 7.
# Needs: dplyr, ggplot2, sf.
# =============================================================================

for (pkg in c("dplyr", "ggplot2", "sf")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}
library(dplyr)
library(ggplot2)
suppressMessages(library(sf))

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

master_file <- file.path("Outputs", "5b_great_plains", "master_data_unique.csv")
traits_file <- file.path("Outputs", "7_rowan_figures", "species_traits.csv")
bins_file   <- file.path("Outputs", "3d_resolved", "time_bins.csv")
output_dir  <- file.path(work_dir, "Outputs", "7_rowan_figures")

regions        <- c("Basin and Range", "Coastal Plain", "Great Plains")
exclude_taxa   <- character(0)   # e.g. c("Homo", "Bos", "Sus") to drop humans and historic introductions
                                 # (left empty so the species match Step 7)
spatial_methods <- c("clustr", "cookies")   # one figure per method; clustr first (main)
cell_km        <- 100            # equal-area cell width (km)
clustr_km      <- 1500           # clustr: maximum diameter of a cluster (km)
cookie_r_km    <- 750            # cookies: radius of the circle (km; Antell et al. 2020 used 750 and 1,500)
min_cells      <- 5              # province-bins with fewer occupied cells are insufficient
quota_cells    <- NA             # NA = automatic (smallest cell count among sufficient province-bins);
                                 # or a number, e.g. 12 as recommended by Antell et al. (2024)
n_iter         <- 999            # iterations (spatial subsample + equal-species draw)
alpha          <- 0.05
random_seed    <- 2024
show_counts    <- FALSE          # TRUE = print cells and species used under each bin label

map_crs <- "+proj=laea +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

band_fill   <- "#F6EBD0"    # shading of alternate time bins (pale sand, as Fig. 2)
ink         <- "#1F2A44"    # text, axes, lines and outlines (dark navy)
ink_soft    <- "#3E4C6D"    # secondary text
point_fill  <- "#F0E442"    # beta_SIM points (Okabe-Ito yellow, navy outline)
page_fill   <- "white"
fig_w <- 9; fig_h <- 6; fig_dpi <- 300; base_size <- 14

# =============================================================================
# HELPERS
# =============================================================================

beta_sim <- function(x) {                      # multisite Simpson (Rowan et al. 2024 Eq. 1)
  x <- x[, colSums(x) > 0, drop = FALSE]
  cmb <- combn(nrow(x), 2)
  mins <- vapply(seq_len(ncol(cmb)), function(k) {
    i <- cmb[1, k]; j <- cmb[2, k]
    min(sum(x[i, ] == 1 & x[j, ] == 0), sum(x[i, ] == 0 & x[j, ] == 1))
  }, numeric(1))
  shared <- sum(x) - ncol(x)
  if (shared + sum(mins) == 0) return(NA_real_)
  sum(mins) / (shared + sum(mins))
}
lists_to_matrix <- function(lists) {
  spp <- sort(unique(unlist(lists)))
  m <- matrix(0L, length(lists), length(spp), dimnames = list(names(lists), spp))
  for (i in seq_along(lists)) m[i, spp %in% lists[[i]]] <- 1L
  m
}
gc_dist_km <- function(lon, lat) {             # great-circle distances (km)
  r <- pi / 180; la <- lat * r; lo <- lon * r
  d <- outer(sin(la), sin(la)) + outer(cos(la), cos(la)) * cos(outer(lo, lo, "-"))
  d[] <- pmin(1, pmax(-1, d))
  6371.0088 * acos(d)
}
# divvy's algorithms (Antell et al. 2024): candidate pools per seed cell.
cookie_pools <- function(D, n) {
  pools <- lapply(seq_len(nrow(D)), function(i) which(D[i, ] <= cookie_r_km))
  pools[lengths(pools) >= n]
}
clustr_pools <- function(D, n) {
  grow <- function(seed) {
    clust <- seed
    repeat {
      rest <- setdiff(seq_len(nrow(D)), clust); if (!length(rest)) break
      dd <- apply(D[clust, rest, drop = FALSE], 2, min); cand <- rest[dd == min(dd)]
      nxt <- if (length(cand) > 1) sample(cand, 1) else cand
      if (max(D[c(clust, nxt), c(clust, nxt)]) > clustr_km) break
      clust <- c(clust, nxt)
    }
    clust
  }
  pools <- lapply(seq_len(nrow(D)), grow)
  pools[lengths(pools) >= n]
}
fmt_age <- function(x) ifelse(x < 0.1, sprintf("%.4f", x), sprintf("%.2f", x))
fmt_p   <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
stars   <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
read_csv_utf8 <- function(f) {
  path <- file.path(work_dir, f)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun the earlier steps first.")
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
}

# =============================================================================
# 1. DATA: SITES -> EQUAL-AREA CELLS
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
time_bins <- read_csv_utf8(bins_file) %>%
  transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma), Younger = as.numeric(Younger_Ma)) %>%
  arrange(desc(Older)) %>% mutate(Display_Bin = Time_Bin, Mid = (Older + Younger) / 2)
traits <- read_csv_utf8(traits_file)
large  <- traits$Species[as.logical(traits$Included) %in% TRUE]
occ <- read_csv_utf8(master_file) %>%
  mutate(Time_Bin = as.integer(Time_Bin), Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) %>%
  filter(Spatial_Bin %in% regions, Species %in% large, Time_Bin %in% time_bins$Time_Bin,
         !is.na(Latitude), !is.na(Longitude),
         !(Species %in% exclude_taxa | sub(" .*$", "", Species) %in% exclude_taxa)) %>%
  distinct(Time_Bin, Spatial_Bin, Site_Key, Latitude, Longitude, Species)
xy <- st_coordinates(st_transform(st_as_sf(occ, coords = c("Longitude", "Latitude"), crs = 4326), map_crs))
occ$cell <- paste(floor(xy[, 1] / (cell_km * 1000)), floor(xy[, 2] / (cell_km * 1000)), sep = "_")
cell_xy <- occ %>% distinct(cell) %>%
  mutate(i = as.numeric(sub("_.*", "", cell)), j = as.numeric(sub(".*_", "", cell)))
cc <- st_coordinates(st_transform(st_as_sf(data.frame(X = (cell_xy$i + 0.5) * cell_km * 1000, Y = (cell_xy$j + 0.5) * cell_km * 1000),
                                           coords = c("X", "Y"), crs = map_crs), 4326))
cell_xy$lon <- cc[, 1]; cell_xy$lat <- cc[, 2]

cells_pb <- occ %>% group_by(Time_Bin, Spatial_Bin) %>%
  summarise(n_sites = n_distinct(Site_Key), n_cells = n_distinct(cell), n_species = n_distinct(Species), .groups = "drop") %>%
  tidyr::complete(Time_Bin = time_bins$Time_Bin, Spatial_Bin = regions, fill = list(n_sites = 0L, n_cells = 0L, n_species = 0L))
quota <- if (is.na(quota_cells)) min(cells_pb$n_cells[cells_pb$n_cells >= min_cells]) else quota_cells
cells_pb <- cells_pb %>% mutate(Enough_cells = n_cells >= max(min_cells, quota))
cat(sprintf("  %d-km cells; cell quota per province = %d (%s)\n", cell_km, quota,
            if (is.na(quota_cells)) "smallest count among province-bins with at least min_cells cells" else "set by hand"))
print(as.data.frame(cells_pb %>% transmute(Bin = Time_Bin, Province = Spatial_Bin, sites = n_sites, cells = n_cells,
                                           species = n_species, enough = Enough_cells)), row.names = FALSE)

# =============================================================================
# 2. ITERATIONS: SPATIAL SUBSAMPLE (AREA + EXTENT) -> EQUAL SPECIES -> beta_SIM
# =============================================================================

run_method <- function(method) {
  cat(sprintf("\n=== 2. %s: %d iterations per bin ===\n", toupper(method), n_iter))
  set.seed(random_seed)
  draws <- list(); status <- list()
  for (k in seq_len(nrow(time_bins))) {
    b <- time_bins$Time_Bin[k]
    pools <- list(); cell_species <- list(); why <- character(0)
    for (p in regions) {
      d <- occ %>% filter(Time_Bin == b, Spatial_Bin == p)
      cells <- unique(d$cell)
      if (length(cells) < max(min_cells, quota)) { why <- c(why, sprintf("%s: %d cells", p, length(cells))); next }
      ce <- cell_xy[match(cells, cell_xy$cell), ]
      D <- gc_dist_km(ce$lon, ce$lat)
      pl <- if (method == "clustr") clustr_pools(D, quota) else cookie_pools(D, quota)
      if (!length(pl)) { why <- c(why, sprintf("%s: no %s holds %d cells", p, method, quota)); next }
      pools[[p]] <- lapply(pl, function(ix) cells[ix])
      cell_species[[p]] <- split(d$Species, d$cell)
    }
    ok <- !length(why); why <- paste(why, collapse = "; ")   # list every province that fails
    status[[length(status) + 1]] <- data.frame(Display_Bin = b, Sufficient = ok, Reason = why)
    if (!ok) { cat(sprintf("  Bin %d: insufficient (%s)\n", b, why)); next }
    vals <- numeric(n_iter); n_used <- integer(n_iter)
    for (r in seq_len(n_iter)) {
      lists <- lapply(regions, function(p) {
        pool <- pools[[p]][[sample.int(length(pools[[p]]), 1)]]          # random seed cell -> its cluster / circle
        chosen <- pool[sample.int(length(pool), quota)]                    # rarefy to the cell quota
        unique(unlist(cell_species[[p]][chosen]))
      })
      names(lists) <- regions
      n <- min(lengths(lists)); n_used[r] <- n
      lists <- lapply(lists, function(s) if (length(s) > n) sample(s, n) else s)   # equal species
      vals[r] <- beta_sim(lists_to_matrix(lists))
    }
    draws[[as.character(b)]] <- vals
    cat(sprintf("  Bin %d: mean beta_SIM = %.3f (95%% range %.3f-%.3f); species per province after the area cut: median %d\n",
                b, mean(vals), quantile(vals, 0.025), quantile(vals, 0.975), as.integer(median(n_used))))
  }
  status <- bind_rows(status)
  res <- bind_rows(lapply(names(draws), function(b) {
    v <- draws[[b]]
    data.frame(Display_Bin = as.integer(b), beta_SIM_mean = mean(v), beta_SIM_sd = sd(v),
               beta_SIM_lower95 = unname(quantile(v, 0.025)), beta_SIM_upper95 = unname(quantile(v, 0.975)))
  })) %>% right_join(time_bins, by = "Display_Bin") %>% left_join(status, by = "Display_Bin") %>%
    mutate(Method = method, Cell_km = cell_km, Quota_cells = quota) %>% arrange(desc(Older))

  # successive-bin tests on the iteration distributions
  have <- res$Display_Bin[!is.na(res$beta_SIM_mean)]
  tests <- bind_rows(lapply(seq_len(nrow(res) - 1), function(k) {
    a <- res$Display_Bin[k]; bb <- res$Display_Bin[k + 1]
    if (!(a %in% have && bb %in% have)) return(NULL)
    diff <- draws[[as.character(a)]] - sample(draws[[as.character(bb)]])
    p <- min(1, 2 * min((sum(diff <= 0) + 1) / (n_iter + 1), (sum(diff >= 0) + 1) / (n_iter + 1)))
    data.frame(Bin_A = a, Bin_B = bb, diff_mean = mean(diff), diff_lower95 = unname(quantile(diff, 0.025)),
               diff_upper95 = unname(quantile(diff, 0.975)), p = p)
  }))
  if (nrow(tests)) tests <- tests %>% mutate(p_holm = p.adjust(p, method = "holm"), Significant = p_holm < alpha)
  list(res = res, tests = tests)
}

# =============================================================================
# 3. FIGURE (style of Fig. 2)
# =============================================================================

plot_method <- function(out, method) {
  res <- out$res; tests <- out$tests
  x_max <- max(time_bins$Older)
  bands <- time_bins %>% filter(Display_Bin %% 2 == 1)
  # the line only joins neighbouring bins that both have a value
  pts <- res %>% filter(!is.na(beta_SIM_mean)) %>% arrange(Display_Bin) %>%
    mutate(run = cumsum(c(1, diff(Display_Bin) != 1)))
  bin_labels <- res %>% transmute(Mid, lab = paste0(sprintf("Bin %d", Display_Bin),
                                                    ifelse(is.na(beta_SIM_mean), "\ninsufficient", "")))
  p <- ggplot(pts, aes(Mid, beta_SIM_mean)) +
    geom_rect(data = bands, aes(xmin = Younger, xmax = Older, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = band_fill, colour = NA) +
    geom_errorbar(aes(ymin = beta_SIM_lower95, ymax = beta_SIM_upper95), width = 0.07, colour = ink, linewidth = 0.55) +
    geom_line(aes(group = run), colour = ink, linewidth = 0.6) +
    geom_point(shape = 21, fill = point_fill, colour = ink, size = 3.6, stroke = 0.8)
  if (!is.null(tests) && nrow(tests)) {
    brk <- tests %>%
      left_join(res %>% select(Bin_A = Display_Bin, xa = Mid), by = "Bin_A") %>%
      left_join(res %>% select(Bin_B = Display_Bin, xb = Mid), by = "Bin_B") %>%
      mutate(y = ifelse(row_number() %% 2 == 1, 1.06, 1.15),
             xa = xa - 0.03 * sign(xa - xb), xb = xb + 0.03 * sign(xa - xb),
             lab = trimws(paste(fmt_p(p_holm), stars(p_holm))))
    p <- p +
      geom_segment(data = brk, aes(x = xa, xend = xb, y = y, yend = y), inherit.aes = FALSE, colour = ink, linewidth = 0.45) +
      geom_segment(data = brk, aes(x = xa, xend = xa, y = y, yend = y - 0.02), inherit.aes = FALSE, colour = ink, linewidth = 0.45) +
      geom_segment(data = brk, aes(x = xb, xend = xb, y = y, yend = y - 0.02), inherit.aes = FALSE, colour = ink, linewidth = 0.45) +
      geom_text(data = brk, aes((xa + xb) / 2, y + 0.035, label = lab), inherit.aes = FALSE,
                size = base_size / 4.4, colour = ink, fontface = ifelse(brk$Significant, "bold", "plain"))
  }
  if (show_counts) {
    cnt <- cells_pb %>% group_by(Display_Bin = Time_Bin) %>% summarise(c = paste(n_cells, collapse = "/"), .groups = "drop")
    bin_labels <- bin_labels %>% left_join(res %>% select(Mid, Display_Bin), by = "Mid") %>%
      left_join(cnt, by = "Display_Bin") %>% mutate(lab = paste0(lab, "\ncells ", c))
  }
  p <- p +
    geom_text(data = bin_labels, aes(Mid, 0.02, label = lab), vjust = 0, inherit.aes = FALSE, size = base_size / 4.3,
              colour = ink_soft, lineheight = 0.9) +
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
    coord_cartesian(xlim = c(x_max, 0), ylim = c(0, 1.2), clip = "off") +
    labs(x = "Age (Ma)", y = expression(beta[SIM]),
         caption = sprintf(paste0("Each province subsampled to the same area (%d grid cells of %d km, drawn within a %s),\n",
                                  "then to the species count of the poorest province (Antell, Benson & Saupe 2024).\n",
                                  "Points = mean of %d iterations; error bars = 95%% range. Lines join neighbouring bins only.\n",
                                  "Brackets: p for the difference between successive bins (iterations paired at random; Holm-corrected).\n",
                                  "Bold with * p < 0.05, ** p < 0.01, *** p < 0.001 = significant. Insufficient = a province had too few cells."),
                           quota, cell_km,
                           if (method == "clustr") sprintf("cluster no wider than %d km", clustr_km) else sprintf("circle of %d km radius", cookie_r_km),
                           n_iter)) +
    theme_classic(base_size = base_size) +
    theme(axis.line = element_line(colour = ink, linewidth = 0.4),
          axis.ticks = element_line(colour = ink, linewidth = 0.4),
          axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
          axis.title.y = element_text(size = base_size * 1.15),
          plot.caption = element_text(colour = ink_soft, size = base_size * 0.62, hjust = 0),
          plot.margin = margin(10, 40, 8, 10), plot.background = element_rect(fill = page_fill, colour = NA))
  name <- paste0("Fig2_area_species_", method)
  ggsave(file.path(output_dir, paste0(name, ".png")), p, width = fig_w, height = fig_h, dpi = fig_dpi, bg = page_fill)
  ggsave(file.path(output_dir, paste0(name, ".pdf")), p, width = fig_w, height = fig_h, bg = page_fill,
         device = if (capabilities("cairo")) cairo_pdf else pdf)
  cat("  ", file.path(output_dir, paste0(name, ".png")), " (+ .pdf)\n", sep = "")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cells_pb, file.path(output_dir, "area_species_cells_per_province_bin.csv"), row.names = FALSE)
for (m in spatial_methods) {
  out <- run_method(m)
  if (!is.null(out$tests) && nrow(out$tests)) {
    cat("  Successive bins (difference = older - younger):\n")
    print(as.data.frame(out$tests %>% transmute(Pair = sprintf("Bin %d vs Bin %d", Bin_A, Bin_B),
                                                Difference = sprintf("%+.3f", diff_mean), p = fmt_p(p), p_Holm = fmt_p(p_holm),
                                                Result = ifelse(Significant, "significant", "not significant"))), row.names = FALSE)
  }
  write.csv(out$res %>% mutate(across(starts_with("beta"), ~ round(.x, 4))),
            file.path(output_dir, paste0("beta_sim_area_species_", m, ".csv")), row.names = FALSE)
  if (!is.null(out$tests) && nrow(out$tests))
    write.csv(out$tests %>% mutate(across(starts_with("diff"), ~ round(.x, 4)), across(c(p, p_holm), ~ signif(.x, 3))),
              file.path(output_dir, paste0("beta_sim_area_species_", m, "_tests.csv")), row.names = FALSE)
  if (any(!is.na(out$res$beta_SIM_mean))) plot_method(out, m)
}
cat("\n=== STEP 7 ADDENDUM (AREA + SPECIES) COMPLETE ===\n")
