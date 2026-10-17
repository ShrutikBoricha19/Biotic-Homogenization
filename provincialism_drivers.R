# =============================================================================
# WHAT MAINTAINED NORTH AMERICAN PROVINCIALISM?  (analysis plan, Modules 0-6)
#
# Tests six candidate causes of provincialism among the Basin and Range,
# Coastal Plain and Great Plains (3.25 Ma - 11.7 ka) on the Chapter 1 data:
#   H0 sampling artefact      -> Module 0 (gate)
#   H1 environmental filtering, H2 barriers -> Modules 1-3
#   H3 immigration history    -> Module 5 = ch2_immigrants.R (run separately)
#   H4 in-place origination   -> Module 6
#   H5 stochastic assembly    -> Module 4
# Choose the modules to run in 'run_modules' (SETTINGS). Each writes tables
# and figures to Outputs/provincialism/<module>/; the synthesis collects the
# key numbers into one table laid out like "Reading the result" in the plan.
#
# Module 0  Sampling gate (Antell, Benson & Saupe 2024; divvy)
#   - sites -> equal-area grid cells; within every province and bin the cells
#     are subsampled to a fixed extent (circular 'cookies', or nearest-neighbour
#     'clustr' clusters of maximum diameter distMax) and rarefied to a fixed
#     cell quota; multisite beta_SIM (Rowan et al. 2024 Eq. 1) per subsample;
#     the age of every locality is redrawn within its Max-Min range in every
#     iteration, so one set of iterations carries spatial, area and age
#     uncertainty; spread reported as quantile ranges (not CIs);
#   - province-bins with fewer cells than the quota are reported as
#     insufficient (the quota is not lowered for single bins); lower quotas are
#     run as a separate, clearly labelled sensitivity check;
#   - Rowan et al.'s sensitivity checks: shifted bin starts (0.15-Myr steps),
#     identifications (all names vs names "as identified" only - the full
#     cf./aff. four-matrix check needs Step 3d variants, see notes), and
#     regressions of beta_SIM and its diet/size partitions on sampling tallies
#     (Holm-corrected);
#   - power analysis by simulation: known drops in beta_SIM are planted into
#     resampled data and the Chapter 1 site-bootstrap test is applied.
#   The in-script cookies/clustr replicate divvy's algorithms (findPool /
#   groupr) on a precomputed great-circle distance matrix, which is much faster
#   inside the age-resampling loop; set use_divvy = TRUE to call divvy itself.
#
# Module 1  Data-drawn bioregions
#   - bipartite cell x species network clustered with Infomap (igraph's
#     implementation; Vilhena & Antonelli 2015; Edler et al. 2017) and, as a
#     cross-check, UPGMA clustering of pairwise beta_SIM (Kreft & Jetz 2010);
#   - adjusted Rand index against the Fenneman provinces; enrichment of
#     cluster boundaries on the Rocky Mountain front and the Mississippi;
#   - writes per-bin occurrence files for the Infomap Bioregions application.
# Module 2  Distance decay (Nekola & White 1999; Soininen et al. 2007)
#   - pairwise beta_SIM (Baselga 2010) between cells; exponential decay of
#     similarity with great-circle distance (GLM, log link) per bin and per
#     functional group; cell-bootstrap intervals; Mantel tests (Darroch et
#     al. 2014).
# Module 3  Climate vs distance vs barriers (Ferrier et al. 2007, gdm package)
#   - GDM on Sorensen and on pairwise beta_SIM; predictors: geographic distance,
#     Rocky Mountain front and Mississippi barriers (0/1 per side, as the
#     nominal-variable method of Ferrier et al. 2007), and - when you supply
#     rasters - climate (Oscillayers / PaleoClim), open habitat (Andermann et
#     al. 2022) and topographic roughness (from a DEM); permutation
#     significance and backward elimination (gdm.varImp); [least-cost distance
#     (Ferrier et al. 2007, Fig. 3) is NOT implemented yet - see the notes];
#     variation partitioning with dbRDA (Borcard et al. 1992; Peres-Neto et al.
#     2006); separate grazer and browser models.
# Module 4  Raup-Crick (Chase et al. 2011): occupancy-weighted null, per cell
#   pair and province pair, per bin and functional group.
# Module 6  Origin of endemics: originated in place / immigrant endemic /
#   relict (extending Darroch et al. 2014), and PyRate input files per
#   province (Silvestro et al. 2014) for origination rates.
#
# MATERIALS (optional; the script runs without them and says what is missing):
#   climate_index.csv  - one row per raster: variable, age_ka, file
#                        (variables used: MAT, MAP, seasonality, open_habitat)
#   dem_file           - an elevation raster (GeoTIFF) for topographic roughness
#   taxon_origins.csv  - the checked origin table from ch2_immigrants.R
#
# Inputs:  Outputs/5b_great_plains/master_data_unique.csv (Step 5b),
#          Outputs/7_rowan_figures/species_traits.csv (Step 7),
#          Outputs/3d_resolved/time_bins.csv,
#          Outputs/3_spatial/spatial_layers.rds and Outputs/maps/basemap (maps).
# Needs: dplyr, tidyr, ggplot2, sf; igraph (Module 1), vegan (Modules 2-3),
#        gdm (Module 3), terra (rasters, optional), divvy (only if use_divvy).
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "sf")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}
library(dplyr)
library(tidyr)
library(ggplot2)
suppressMessages(library(sf))
suppressMessages(sf_use_s2(FALSE))

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"
run_modules <- c(0, 1, 2, 3, 4, 6)      # Module 5 is ch2_immigrants.R

master_file  <- file.path("Outputs", "5b_great_plains", "master_data_unique.csv")
traits_file  <- file.path("Outputs", "7_rowan_figures", "species_traits.csv")
bins_file    <- file.path("Outputs", "3d_resolved", "time_bins.csv")
layers_file  <- file.path("Outputs", "3_spatial", "spatial_layers.rds")
rivers_file  <- file.path("Outputs", "maps", "basemap", "ne_10m_rivers_lake_centerlines.geojson")
origins_file <- "taxon_origins.csv"                     # from ch2_immigrants.R (checked)
climate_index_file <- "climate_index.csv"               # optional (Module 3)
dem_file     <- NA                                      # optional elevation GeoTIFF (Module 3)
output_root  <- file.path(work_dir, "Outputs", "provincialism")

regions <- c("Basin and Range", "Coastal Plain", "Great Plains")
exclude_taxa <- c("Homo", "Bos", "Sus", "Capra", "Felis catus", "Canis familiaris", "Ovis aries",
                  "Equus asinus")                       # humans and historic introductions
map_crs <- "+proj=laea +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
random_seed <- 2024

# Module 0
cell_km        <- 100        # equal-area cell width (Antell et al. 2024 used ~100 km cells)
quota_main     <- 12         # cells per province and bin (Antell et al. 2024: >= 12)
quota_check    <- c(8, 5)    # lower quotas, reported only as a labelled sensitivity check
cookie_r_km    <- 750        # cookies radius (Antell et al. 2020: 750 and 1,500 km)
clustr_km      <- 1500       # clustr maximum diameter (divvy distMax; Womack et al. ~1,000-1,100 km)
n_age_reps     <- 100        # age redraws ...
n_spatial      <- 10         # ... x spatial subsamples each = 1,000 iterations
use_divvy      <- FALSE      # TRUE = call divvy::cookies / divvy::clustr (slower; same algorithms)
age_modes      <- c("redraw", "midpoint")  # "redraw" = ages drawn within Max-Min (main); "midpoint" = Chapter 1 binning,
                                           # run alongside to show how much of the uncertainty comes from ages
bin_shifts     <- c(0, 0.15, 0.30, 0.45, 0.60)   # Myr, younger shifts of the bin starts
n_draws_eq     <- 199        # equal-species draws for the sensitivity curves
power_drops    <- c(0, 0.05, 0.1, 0.2, 0.35, 0.5) # homogenization strengths planted
n_power_sims   <- 40
n_power_boot   <- 199

# Modules 1-4
grid_km_net    <- 200        # cell width for networks / pairwise analyses (sparser data: 2 degrees ~ 200 km)
min_cell_rich  <- 3          # minimum species per cell for pairwise analyses
n_boot_decay   <- 499
n_mantel       <- 999
n_rc           <- 999        # Raup-Crick null draws
gdm_perm       <- 50         # permutations for gdm.varImp (slow; raise for the final run)
min_cells_model <- 10        # minimum cells per bin for GDM / clustering
func_groups    <- list(All = NULL, Grazer = "Grazer", Browser = "Browser", `Mixed feeder` = "Mixed feeder",
                       `Size 3 (80-350 kg)` = "Size 3 (80-350 kg)", `Size 5 (>1,000 kg)` = "Size 5 (>1,000 kg)")

# Colours (colourblind-friendly, as Step 7)
ink <- "#1F2A44"; ink_soft <- "#3E4C6D"; band_fill <- "#F6EBD0"; strip_fill <- "#D6E6F5"
prov_cols <- c("Basin and Range" = "#E69F00", "Coastal Plain" = "#56B4E9", "Great Plains" = "#009E73")
okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#1F2A44")
base_size <- 13; fig_dpi <- 300

# =============================================================================
# SHARED HELPERS
# =============================================================================

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
mod_dir <- function(m) { d <- file.path(output_root, m); dir.create(d, showWarnings = FALSE); d }
read_wd <- function(f) {
  p <- if (file.exists(f)) f else file.path(work_dir, f)
  if (!file.exists(p)) stop("File not found:\n  ", p, "\nRun the earlier steps first.")
  read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
}
need <- function(pkg, module) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) message("  ", module, ": package '", pkg, "' not installed - install.packages(\"", pkg, "\")")
  ok
}
save_fig <- function(p, dir, name, w, h) {
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = fig_dpi, bg = "white")
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}
theme_p <- function() {
  theme_classic(base_size = base_size) +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink_soft),
          axis.line = element_line(colour = ink), axis.ticks = element_line(colour = ink),
          strip.background = element_rect(fill = strip_fill, colour = ink_soft), strip.text = element_text(colour = ink),
          legend.position = "bottom", plot.caption = element_text(colour = ink_soft, hjust = 0, size = base_size * 0.7))
}
two_sided_p <- function(null, obs) {
  null <- null[is.finite(null)]
  if (!length(null) || !is.finite(obs)) return(NA_real_)
  min(1, 2 * min((sum(null <= obs) + 1) / (length(null) + 1), (sum(null >= obs) + 1) / (length(null) + 1)))
}
qr95 <- function(x) { x <- x[is.finite(x)]; if (!length(x)) c(NA, NA, NA) else unname(quantile(x, c(0.025, 0.5, 0.975))) }

# Multisite beta_SIM (Rowan et al. 2024 Eq. 1 = betapart::beta.multi SIM).
beta_sim <- function(x) {
  x <- x[, colSums(x) > 0, drop = FALSE]
  if (nrow(x) < 2 || ncol(x) == 0) return(NA_real_)
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
beta_equal <- function(lists, n) {
  picks <- lapply(lists, function(s) if (length(s) > n) sample(s, n) else s)
  beta_sim(lists_to_matrix(picks))
}
# Pairwise beta_SIM (Baselga 2010) between the rows of a 0/1 matrix.
pair_beta_sim <- function(m) {
  a <- m %*% t(m); s <- rowSums(m)
  b <- matrix(s, nrow(m), nrow(m)) - a; c <- t(b)
  out <- pmin(b, c) / (a + pmin(b, c)); out[!is.finite(out)] <- NA
  dimnames(out) <- list(rownames(m), rownames(m)); out
}
gc_dist_km <- function(lon, lat) {                 # great-circle distances (km)
  r <- pi / 180; la <- lat * r; lo <- lon * r
  d <- outer(la, la, function(a, b) sin(a) * sin(b)) +
    outer(la, la, function(a, b) cos(a) * cos(b)) * cos(outer(lo, lo, "-"))
  d[] <- pmin(1, pmax(-1, d))             # keep the matrix dimensions
  6371.0088 * acos(d)
}
# Rowan et al. (2024) Eq. 4 partition (as Step 7).
endemic_species <- function(x) {
  cmb <- combn(nrow(x), 2)
  unlist(lapply(seq_len(ncol(cmb)), function(k) {
    y <- x[cmb[, k], , drop = FALSE]
    b_ij <- y[1, ] == 1 & y[2, ] == 0; b_ji <- y[1, ] == 0 & y[2, ] == 1
    if (sum(b_ij) < sum(b_ji)) colnames(y)[b_ij] else colnames(y)[b_ji]
  }))
}
simp_func <- function(x, grp) {
  e <- endemic_species(x); f <- factor(grp[match(e, colnames(x))], levels = levels(grp))
  table(f) / ((sum(x) - sum(colSums(x) > 0)) + length(e))
}

# =============================================================================
# DATA
# =============================================================================

cat("=== READING INPUTS ===\n")
setwd(work_dir)
md <- read_wd(master_file) %>%
  mutate(Time_Bin = as.integer(Time_Bin), Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude),
         Max_Ma = as.numeric(Max_Ma), Min_Ma = as.numeric(Min_Ma), Midpoint_Ma = as.numeric(Midpoint_Ma))
if (!"Genus" %in% names(md)) md$Genus <- sub(" .*$", "", md$Species)
traits <- read_wd(traits_file)
time_bins <- read_wd(bins_file) %>%
  transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma), Younger = as.numeric(Younger_Ma)) %>%
  arrange(Time_Bin) %>% mutate(Mid = (Older + Younger) / 2)
win_old <- max(time_bins$Older); win_young <- min(time_bins$Younger)
bin_breaks <- c(time_bins$Older, win_young)

large <- traits$Species[as.logical(traits$Included) %in% TRUE]
is_excluded <- function(sp) sp %in% exclude_taxa | sub(" .*$", "", sp) %in% exclude_taxa
occ_na <- md %>% filter(Species %in% large, !is_excluded(Species), !is.na(Latitude), !is.na(Longitude))
occ <- occ_na %>% filter(Spatial_Bin %in% regions, Time_Bin %in% time_bins$Time_Bin)
n_excl <- md %>% filter(Species %in% large, is_excluded(Species)) %>% distinct(Species)
cat(sprintf("  %d occurrences, %d species, %d sites in the three provinces%s\n", nrow(occ), n_distinct(occ$Species),
            n_distinct(occ$Site_Key),
            if (nrow(n_excl)) paste0("; excluded: ", paste(n_excl$Species, collapse = ", ")) else ""))

tr <- traits %>% transmute(Species, Diet = Diet_Group, Size = Size_Class)
diet_of <- setNames(tr$Diet, tr$Species); size_of <- setNames(tr$Size, tr$Species)
in_group <- function(sp, g) {
  if (is.null(g)) return(rep(TRUE, length(sp)))
  (diet_of[sp] %in% g) | (size_of[sp] %in% g)
}

# Sites with projected coordinates and equal-area cells.
if (!"West_Of_Mississippi" %in% names(occ_na)) occ_na$West_Of_Mississippi <- NA
sites <- occ_na %>% group_by(Site_Key) %>%
  summarise(Latitude = first(Latitude), Longitude = first(Longitude), Max_Ma = max(Max_Ma, na.rm = TRUE),
            Min_Ma = min(Min_Ma, na.rm = TRUE), Spatial_Bin = first(Spatial_Bin),
            West_Of_Mississippi = first(West_Of_Mississippi), .groups = "drop") %>%
  filter(is.finite(Max_Ma), is.finite(Min_Ma))
xy <- st_coordinates(st_transform(st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326), map_crs))
sites$X <- xy[, 1]; sites$Y <- xy[, 2]
cell_of <- function(X, Y, km) paste(floor(X / (km * 1000)), floor(Y / (km * 1000)), sep = "_")
cell_centroids <- function(cells, km) {
  ij <- do.call(rbind, strsplit(cells, "_")); storage.mode(ij) <- "numeric"
  pts <- st_as_sf(data.frame(cell = cells, X = (ij[, 1] + 0.5) * km * 1000, Y = (ij[, 2] + 0.5) * km * 1000),
                  coords = c("X", "Y"), crs = map_crs)
  ll <- st_coordinates(st_transform(pts, 4326))
  data.frame(cell = cells, lon = ll[, 1], lat = ll[, 2], X = (ij[, 1] + 0.5) * km * 1000, Y = (ij[, 2] + 0.5) * km * 1000)
}
west_of_miss <- function(x) if (is.logical(x)) x else toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES")

# Cell x species data for one bin (by binned midpoint), optional functional group.
cell_data <- function(b, km, group = NULL, data = occ) {
  d <- data %>% filter(Time_Bin == b, in_group(Species, group)) %>%
    left_join(sites %>% select(Site_Key, X, Y), by = "Site_Key") %>%
    mutate(cell = cell_of(X, Y, km))
  if (!nrow(d)) return(NULL)
  prov <- d %>% count(cell, Spatial_Bin) %>% group_by(cell) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(cell, Province = Spatial_Bin)
  side <- d %>% group_by(cell) %>% summarise(West_Miss = mean(west_of_miss(West_Of_Mississippi)) >= 0.5, .groups = "drop")
  pa <- d %>% distinct(cell, Species)
  m <- table(pa$cell, pa$Species); m <- unclass(m); m[m > 0] <- 1L
  storage.mode(m) <- "integer"
  keep <- rowSums(m) >= min_cell_rich
  info <- cell_centroids(rownames(m), km) %>% left_join(prov, by = "cell") %>% left_join(side, by = "cell") %>%
    mutate(richness = rowSums(m), West_Rockies = Province == "Basin and Range")
  list(m = m[keep, colSums(m[keep, , drop = FALSE]) > 0, drop = FALSE], info = info[keep, ])
}

# Study-province geometry for maps (as Step 10); NULL if the layers are missing.
study_geometry <- function() {
  if (!file.exists(layers_file)) return(NULL)
  out <- tryCatch({
    prov <- st_transform(st_make_valid(readRDS(layers_file)$province), 4326)
    un <- function(u) suppressMessages(st_union(st_make_valid(prov[prov$Unit %in% u, ])))
    g <- lapply(regions, un); names(g) <- regions
    river <- if (file.exists(rivers_file)) { r <- st_read(rivers_file, quiet = TRUE); r[r$name %in% "Mississippi", ] } else NULL
    if (!is.null(river) && "Central Lowland" %in% prov$Unit) {
      rxy <- st_coordinates(river)[, 1:2]
      bands <- seq(floor(min(rxy[, 2])), ceiling(max(rxy[, 2])), by = 0.05)
      bl <- vapply(bands, function(b) { k <- abs(rxy[, 2] - b) <= 0.05; if (any(k)) min(rxy[k, 1]) else NA }, 0)
      f <- approxfun(bands[!is.na(bl)], bl[!is.na(bl)], rule = 2); ls <- seq(20, 60, 0.05)
      west <- st_sfc(st_polygon(list(rbind(cbind(f(ls), ls), c(-150, 60), c(-150, 20), c(f(20), 20)))), crs = 4326)
      g[["Great Plains"]] <- suppressMessages(st_union(c(g[["Great Plains"]], st_intersection(un("Central Lowland"), west))))
    }
    list(study = st_sf(Province = regions, geometry = do.call(c, unname(g)), crs = 4326), river = river)
  }, error = function(e) { message("  maps: study geometry not available (", conditionMessage(e), ")"); NULL })
  out
}

summary_rows <- list()      # key results for the synthesis
add_summary <- function(module, hypothesis, metric, bin, value, note = "") {
  summary_rows[[length(summary_rows) + 1]] <<- data.frame(Module = module, Hypothesis = hypothesis, Metric = metric,
                                                         Time_Bin = bin, Value = value, Note = note)
}
x_scale <- scale_x_reverse(breaks = c(3, 2, 1, 0))
bands_layer <- function() geom_rect(data = time_bins[seq(1, nrow(time_bins), 2), ], inherit.aes = FALSE,
                                    aes(xmin = Younger, xmax = Older, ymin = -Inf, ymax = Inf), fill = band_fill)
set.seed(random_seed)

# =============================================================================
# MODULE 0: SAMPLING GATE (H0)
# =============================================================================

if (0 %in% run_modules) {
  cat("\n=== MODULE 0: SAMPLING GATE ===\n")
  dir0 <- mod_dir("M0_sampling")

  # In-script versions of divvy's cookies (findPool) and clustr (groupr).
  cookie_draw <- function(D, n, r_km, iter) {
    pools <- lapply(seq_len(nrow(D)), function(i) which(D[i, ] <= r_km))
    seeds <- which(lengths(pools) >= n)
    if (!length(seeds)) return(NULL)
    replicate(iter, { s <- seeds[sample.int(length(seeds), 1)]; p <- pools[[s]]
                      p[sample.int(length(p), n)] }, simplify = FALSE)
  }
  grow_tree <- function(seed, D, dmax) {
    clust <- seed
    repeat {
      rest <- setdiff(seq_len(nrow(D)), clust); if (!length(rest)) break
      dd <- apply(D[clust, rest, drop = FALSE], 2, min); cand <- rest[dd == min(dd)]
      nxt <- if (length(cand) > 1) sample(cand, 1) else cand
      if (max(D[c(clust, nxt), c(clust, nxt)]) > dmax) break
      clust <- c(clust, nxt)
    }
    clust
  }
  clustr_draw <- function(D, n, dmax, iter) {
    trees <- lapply(seq_len(nrow(D)), grow_tree, D = D, dmax = dmax)
    ok <- which(lengths(trees) >= n)
    if (!length(ok)) return(NULL)
    replicate(iter, { t <- trees[[ok[sample.int(length(ok), 1)]]]; t[sample.int(length(t), n)] }, simplify = FALSE)
  }
  divvy_draw <- function(cents, n, method, iter) {
    dat <- cents[, c("lon", "lat")]
    res <- tryCatch(if (method == "cookies") divvy::cookies(dat, xy = 1:2, iter = iter, nSite = n, r = cookie_r_km)
                    else divvy::clustr(dat, xy = 1:2, iter = iter, nSite = n, distMax = clustr_km),
                    error = function(e) NULL)
    if (is.null(res)) return(NULL)
    lapply(res, function(s) match(paste(s[[1]], s[[2]]), paste(dat$lon, dat$lat)))
  }

  quotas <- c(quota_main, quota_check)
  occ0 <- occ %>% select(Species, Site_Key, Spatial_Bin)
  site0 <- sites %>% filter(Site_Key %in% occ0$Site_Key) %>% mutate(cell = cell_of(X, Y, cell_km))
  cents_all <- cell_centroids(unique(site0$cell), cell_km)
  sp_by_site <- split(occ0$Species, occ0$Site_Key)
  res0 <- list(); suff <- list()
  for (am in age_modes) for (a in seq_len(n_age_reps)) {
    age <- if (am == "redraw") runif(nrow(site0), pmin(site0$Min_Ma, site0$Max_Ma), pmax(site0$Min_Ma, site0$Max_Ma)) else
      (site0$Min_Ma + site0$Max_Ma) / 2
    bin <- findInterval(-age, -bin_breaks)            # 1 = oldest bin
    bin[age > win_old | age < win_young] <- NA
    for (b in time_bins$Time_Bin) {
      sb <- site0[which(bin == b), ]
      for (meth in c("cookies", "clustr")) for (q in quotas) {
        draws <- list(); ok <- TRUE
        for (p in regions) {
          sp_cells <- sb %>% filter(Spatial_Bin == p)
          cells <- unique(sp_cells$cell)
          if (length(cells) < q) { ok <- FALSE; break }
          cents <- cents_all[match(cells, cents_all$cell), ]
          idx <- if (use_divvy) divvy_draw(cents, q, meth, n_spatial) else {
            D <- gc_dist_km(cents$lon, cents$lat)
            if (meth == "cookies") cookie_draw(D, q, cookie_r_km, n_spatial) else clustr_draw(D, q, clustr_km, n_spatial)
          }
          if (is.null(idx)) { ok <- FALSE; break }
          cell_sp <- lapply(split(sp_cells$Site_Key, sp_cells$cell), function(k) unique(unlist(sp_by_site[k])))
          draws[[p]] <- lapply(idx, function(ix) unique(unlist(cell_sp[cells[ix]])))
        }
        suff[[length(suff) + 1]] <- data.frame(Ages = am, age_rep = a, Time_Bin = b, Method = meth, Quota = q, Sufficient = ok)
        if (!ok) next
        bvals <- vapply(seq_len(n_spatial), function(k) beta_sim(lists_to_matrix(lapply(draws, `[[`, k))), numeric(1))
        res0[[length(res0) + 1]] <- data.frame(Ages = am, age_rep = a, Time_Bin = b, Method = meth, Quota = q, beta_SIM = bvals)
      }
    }
    if (a %% 20 == 0) cat(sprintf("  ages %s: %d / %d\n", am, a, n_age_reps))
  }
  res0 <- bind_rows(res0); suff <- bind_rows(suff)
  curve0 <- suff %>% group_by(Ages, Time_Bin, Method, Quota) %>% summarise(share_sufficient = mean(Sufficient), .groups = "drop") %>%
    left_join(if (nrow(res0)) res0 %>% group_by(Ages, Time_Bin, Method, Quota) %>%
                summarise(n_iter = n(), q025 = qr95(beta_SIM)[1], median = qr95(beta_SIM)[2], q975 = qr95(beta_SIM)[3],
                          .groups = "drop") else
                data.frame(Ages = character(), Time_Bin = integer(), Method = character(), Quota = numeric()),
              by = c("Ages", "Time_Bin", "Method", "Quota")) %>%
    left_join(time_bins, by = "Time_Bin") %>% arrange(Ages, Quota, Method, Time_Bin)
  write.csv(curve0, file.path(dir0, "M0_corrected_beta_sim.csv"), row.names = FALSE)
  cat("  Corrected beta_SIM (median and 95% quantile range of iterations; insufficient bins have no value):\n")
  print(as.data.frame(curve0 %>% transmute(Ages, Quota, Method, Bin = Time_Bin, sufficient = round(share_sufficient, 2),
                                           median = round(median, 3), q025 = round(q025, 3), q975 = round(q975, 3))),
        row.names = FALSE)
  for (i in which(curve0$Quota == quota_main))
    add_summary("M0", "H0", sprintf("corrected beta_SIM median (%s, quota %d, ages %s)", curve0$Method[i], quota_main, curve0$Ages[i]),
                curve0$Time_Bin[i], curve0$median[i],
                sprintf("sufficient in %.0f%% of age redraws", 100 * curve0$share_sufficient[i]))

  if (any(!is.na(curve0$median))) {
    d <- curve0 %>% filter(!is.na(median)) %>% mutate(Quota = factor(paste(Quota, "cells"), levels = paste(quotas, "cells")))
    p <- ggplot(d, aes(Mid, median, colour = Method)) + bands_layer() +
      geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, position = position_dodge(0.12)) +
      geom_line(position = position_dodge(0.12)) + geom_point(size = 2.4, position = position_dodge(0.12)) +
      facet_grid(Ages ~ Quota, labeller = labeller(Ages = c(redraw = "Ages redrawn (main)", midpoint = "Midpoint ages"))) +
      scale_colour_manual(values = c(cookies = "#0072B2", clustr = "#D55E00"), name = NULL) +
      x_scale + labs(x = "Age (Ma)", y = expression(beta[SIM] ~ "(spatially standardized)"),
                     caption = sprintf(paste0("%d-km cells; cookies r = %d km, clustr max. diameter = %d km; %d age redraws x %d subsamples.",
                                              "\nMedian and 95%% quantile range (not a CI). Bins without a point: a province had fewer cells than the quota (insufficient).",
                                              "\nOnly the %d-cell panel is the main result; lower quotas are a sensitivity check."),
                                       cell_km, cookie_r_km, clustr_km, n_age_reps, n_spatial, quota_main)) + theme_p()
    save_fig(p, dir0, "M0_fig1_corrected_beta_sim", 12, 7.5)
  }

  # --- Rowan et al. sensitivity checks -------------------------------------------
  eq_curve <- function(data, breaks) {
    nb <- length(breaks) - 1
    data$bin <- findInterval(-data$Midpoint_Ma, -breaks)
    bind_rows(lapply(seq_len(nb), function(k) {
      d <- data %>% filter(bin == k) %>% distinct(Spatial_Bin, Species)
      lists <- split(d$Species, d$Spatial_Bin)
      if (length(lists) < 2) return(NULL)
      n <- min(lengths(lists)); v <- replicate(n_draws_eq, beta_equal(lists, n))
      data.frame(bin = k, Older = breaks[k], Younger = breaks[k + 1], n_used = n, beta_mean = mean(v),
                 q025 = qr95(v)[1], q975 = qr95(v)[3])
    }))
  }
  shifts <- bind_rows(lapply(bin_shifts, function(s) {
    w <- time_bins$Older[1] - time_bins$Younger[1]
    br <- sort(unique(c(seq(win_old - s, win_young + 1e-6, by = -w), win_young)), decreasing = TRUE)
    eq_curve(occ, br) %>% mutate(Shift_Myr = s)
  })) %>% mutate(Mid = (Older + Younger) / 2)
  write.csv(shifts, file.path(dir0, "M0_sensitivity_bin_shifts.csv"), row.names = FALSE)
  p <- ggplot(shifts, aes(Mid, beta_mean, colour = factor(Shift_Myr))) +
    geom_line() + geom_point(size = 2) + scale_colour_manual(values = okabe, name = "Bin starts shifted (Myr)") +
    x_scale + labs(x = "Age (Ma)", y = expression(beta[SIM] ~ "(equal species)"),
                   caption = "Rowan et al. (2024) bin-placement check: bin starts shifted younger in 0.15-Myr steps.") + theme_p()
  save_fig(p, dir0, "M0_fig2_bin_shift_check", 9, 5.4)

  if ("Name_Status" %in% names(occ)) {
    nom <- bind_rows(eq_curve(occ, bin_breaks) %>% mutate(Names = "All names (main)"),
                     eq_curve(occ %>% filter(Name_Status == "as identified"), bin_breaks) %>%
                       mutate(Names = "Names as identified only")) %>% mutate(Mid = (Older + Younger) / 2)
    write.csv(nom, file.path(dir0, "M0_sensitivity_identifications.csv"), row.names = FALSE)
  }

  tallies <- occ %>% group_by(Time_Bin) %>%
    summarise(n_sites = n_distinct(Site_Key), n_species = n_distinct(Species), .groups = "drop") %>%
    left_join(occ %>% count(Time_Bin, Site_Key) %>% group_by(Time_Bin) %>% summarise(mean_sp_per_site = mean(n), .groups = "drop"),
              by = "Time_Bin")
  part_rows <- lapply(time_bins$Time_Bin, function(b) {
    d <- occ %>% filter(Time_Bin == b) %>% distinct(Spatial_Bin, Species)
    x <- lists_to_matrix(split(d$Species, d$Spatial_Bin)); if (nrow(x) < 2) return(NULL)
    dg <- factor(diet_of[colnames(x)]); names(dg) <- colnames(x)
    sg <- factor(size_of[colnames(x)]); names(sg) <- colnames(x)
    data.frame(Time_Bin = b, Response = c("beta_SIM", paste("Diet:", levels(dg)), paste("Size:", levels(sg))),
               value = c(beta_sim(x), as.numeric(simp_func(x, dg)), as.numeric(simp_func(x, sg))))
  })
  tal <- bind_rows(part_rows) %>% left_join(tallies, by = "Time_Bin")
  reg <- bind_rows(lapply(split(tal, tal$Response), function(d) bind_rows(lapply(c("n_sites", "n_species", "mean_sp_per_site"), function(v) {
    if (nrow(d) < 4 || sd(d[[v]]) == 0 || sd(d$value) == 0) return(NULL)
    f <- summary(lm(d$value ~ d[[v]]))
    data.frame(Response = d$Response[1], Tally = v, slope = coef(f)[2, 1], R2 = f$r.squared, p = coef(f)[2, 4], n_bins = nrow(d))
  })))) %>% mutate(p_holm = p.adjust(p, method = "holm"))
  write.csv(reg, file.path(dir0, "M0_sensitivity_sampling_tallies.csv"), row.names = FALSE)
  cat(sprintf("  Sampling-tally regressions: %d of %d significant after Holm correction (5 bins only: low power)\n",
              sum(reg$p_holm < 0.05, na.rm = TRUE), nrow(reg)))

  # --- power analysis: plant known drops, apply the Chapter 1 test --------------
  cat("  Power analysis ...\n")
  site_lists_bin <- function(b) {
    d <- occ %>% filter(Time_Bin == b)
    lapply(setNames(regions, regions), function(p) { dp <- d %>% filter(Spatial_Bin == p); split(dp$Species, dp$Site_Key) })
  }
  boot_beta <- function(sl, n_boot) {
    n_use <- min(sapply(sl, function(s) length(unique(unlist(s)))))
    replicate(n_boot, { lists <- lapply(sl, function(s) unique(unlist(s[sample.int(length(s), replace = TRUE)])))
                        beta_equal(lists, min(n_use, lengths(lists))) })
  }
  homogenize <- function(sl, h) {          # add species from other provinces into random sites
    allsp <- lapply(sl, function(s) unique(unlist(s)))
    for (p in names(sl)) {
      cand <- setdiff(unique(unlist(allsp[names(sl) != p])), allsp[[p]])
      add <- cand[runif(length(cand)) < h]
      for (s in add) { k <- sample.int(length(sl[[p]]), 1); sl[[p]][[k]] <- c(sl[[p]][[k]], s) }
    }
    sl
  }
  thin_to <- function(sl, sizes) lapply(setNames(names(sl), names(sl)), function(p)
    sl[[p]][sample.int(length(sl[[p]]), min(length(sl[[p]]), sizes[[p]]), replace = length(sl[[p]]) < sizes[[p]])])
  pow <- list()
  for (k in seq_len(nrow(time_bins) - 1)) {
    A <- site_lists_bin(time_bins$Time_Bin[k]); B <- site_lists_bin(time_bins$Time_Bin[k + 1])
    if (any(lengths(A) < 2) || any(lengths(B) < 2)) next
    sizes_B <- lengths(B)
    for (h in power_drops) {
      sims <- replicate(n_power_sims, {
        A1 <- lapply(A, function(s) s[sample.int(length(s), replace = TRUE)])
        hA <- homogenize(A, h)
        B1 <- thin_to(hA, sizes_B)
        la <- lapply(A, function(s) unique(unlist(s))); lb <- lapply(hA, function(s) unique(unlist(s)))
        true_drop <- beta_sim(lists_to_matrix(la)) - beta_sim(lists_to_matrix(lb))
        diff <- boot_beta(A1, n_power_boot) - boot_beta(B1, n_power_boot)
        c(true_drop = true_drop, p = min(1, 2 * min((sum(diff <= 0) + 1) / (n_power_boot + 1), (sum(diff >= 0) + 1) / (n_power_boot + 1))),
          right_direction = mean(diff) > 0)
      })
      pow[[length(pow) + 1]] <- data.frame(Pair = sprintf("Bin %d -> %d", time_bins$Time_Bin[k], time_bins$Time_Bin[k + 1]),
                                           strength = h, mean_planted_drop = mean(sims["true_drop", ]),
                                           power = mean(sims["p", ] < 0.05 & sims["right_direction", ] == 1))
    }
  }
  pow <- bind_rows(pow)
  write.csv(pow, file.path(dir0, "M0_power_analysis.csv"), row.names = FALSE)
  if (nrow(pow)) {
    p <- ggplot(pow, aes(mean_planted_drop, power, colour = Pair)) + geom_hline(yintercept = 0.8, linetype = 2, colour = ink_soft) +
      geom_line() + geom_point(size = 2) + scale_colour_manual(values = okabe, name = NULL) +
      labs(x = expression("Planted drop in " * beta[SIM]), y = "Power (share of simulations with p < 0.05)",
           caption = paste0("Earlier bin resampled; later bin = earlier bin homogenized by a known amount and thinned to the later bin's site numbers;\n",
                            "Chapter 1 site-bootstrap test (", n_power_boot, " iterations, ", n_power_sims, " simulations per point). Dashed: 80% power.")) + theme_p()
    save_fig(p, dir0, "M0_fig3_power", 9, 5.4)
    for (pr in unique(pow$Pair)) {
      d <- pow %>% filter(Pair == pr, power >= 0.8)
      add_summary("M0", "H0", paste("smallest drop detected with 80% power,", pr), NA,
                  if (nrow(d)) min(d$mean_planted_drop) else NA, if (!nrow(d)) "not reached in the simulated range" else "")
    }
  }
}

# =============================================================================
# MODULE 1: DATA-DRAWN BIOREGIONS
# =============================================================================

if (1 %in% run_modules && need("igraph", "Module 1")) {
  cat("\n=== MODULE 1: DATA-DRAWN BIOREGIONS ===\n")
  dir1 <- mod_dir("M1_bioregions")
  ari <- function(a, b) {
    t <- table(a, b); n <- sum(t); c2 <- function(x) x * (x - 1) / 2
    s <- sum(c2(t)); sa <- sum(c2(rowSums(t))); sb <- sum(c2(colSums(t))); e <- sa * sb / c2(n)
    if ((sa + sb) / 2 - e == 0) return(NA_real_)
    (s - e) / ((sa + sb) / 2 - e)
  }
  geo <- study_geometry()
  res1 <- list(); cells1 <- list()
  for (b in time_bins$Time_Bin) {
    cd <- cell_data(b, grid_km_net)
    if (is.null(cd) || nrow(cd$m) < min_cells_model) { cat(sprintf("  Bin %d: too few cells - skipped\n", b)); next }
    m <- cd$m; info <- cd$info
    # Infomap on the bipartite cell-species network
    el <- which(m == 1, arr.ind = TRUE)
    g <- igraph::graph_from_data_frame(data.frame(from = paste0("c:", rownames(m)[el[, 1]]), to = paste0("s:", colnames(m)[el[, 2]])),
                                       directed = FALSE)
    im <- igraph::cluster_infomap(g, nb.trials = 20)
    mem <- igraph::membership(im)
    info$Infomap <- as.integer(factor(mem[paste0("c:", info$cell)]))
    # UPGMA on pairwise beta_SIM (Kreft & Jetz 2010)
    D <- pair_beta_sim(m); D[is.na(D)] <- 1
    hc <- hclust(as.dist(D), method = "average")
    info$UPGMA_k3 <- cutree(hc, k = min(3, nrow(m) - 1))
    sil <- if (requireNamespace("cluster", quietly = TRUE)) sapply(2:min(8, nrow(m) - 1), function(k)
      mean(cluster::silhouette(cutree(hc, k), as.dist(D))[, 3])) else NA
    k_best <- if (all(is.na(sil))) 3 else (2:min(8, nrow(m) - 1))[which.max(sil)]
    info$UPGMA_best <- cutree(hc, k = k_best)
    # boundary enrichment on barriers: adjacent cells (sharing an edge) in different clusters
    ij <- do.call(rbind, strsplit(info$cell, "_")); storage.mode(ij) <- "numeric"
    adj <- which(abs(outer(ij[, 1], ij[, 1], "-")) + abs(outer(ij[, 2], ij[, 2], "-")) == 1 & upper.tri(diag(nrow(info))), arr.ind = TRUE)
    enrich <- function(cl, barrier) {
      if (!nrow(adj)) return(NA)
      cross <- barrier[adj[, 1]] != barrier[adj[, 2]]; bound <- cl[adj[, 1]] != cl[adj[, 2]]
      if (!any(bound) || mean(cross) == 0) return(NA)
      mean(cross[bound]) / mean(cross)
    }
    res1[[length(res1) + 1]] <- data.frame(
      Time_Bin = b, n_cells = nrow(m), n_species = ncol(m),
      Infomap_modules = length(unique(info$Infomap)), UPGMA_best_k = k_best,
      ARI_Infomap_vs_Fenneman = ari(info$Infomap, info$Province), ARI_UPGMA3_vs_Fenneman = ari(info$UPGMA_k3, info$Province),
      ARI_UPGMAbest_vs_Fenneman = ari(info$UPGMA_best, info$Province), ARI_Infomap_vs_UPGMAbest = ari(info$Infomap, info$UPGMA_best),
      Boundary_enrichment_Mississippi = enrich(info$Infomap, info$West_Miss),
      Boundary_enrichment_RockyFront = enrich(info$Infomap, info$West_Rockies))
    cells1[[length(cells1) + 1]] <- info %>% mutate(Time_Bin = b)
    # occurrence file for the Infomap Bioregions application
    occ %>% filter(Time_Bin == b) %>% distinct(Species, Longitude, Latitude) %>%
      write.csv(file.path(dir1, sprintf("M1_infomap_bioregions_input_bin%d.csv", b)), row.names = FALSE)
  }
  res1 <- bind_rows(res1); cells1 <- bind_rows(cells1)
  write.csv(res1, file.path(dir1, "M1_bioregion_summary.csv"), row.names = FALSE)
  write.csv(cells1, file.path(dir1, "M1_cells_with_clusters.csv"), row.names = FALSE)
  if (nrow(res1)) {
    print(as.data.frame(res1 %>% mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
    for (i in seq_len(nrow(res1))) {
      add_summary("M1", "H1/H2", "ARI Infomap vs Fenneman", res1$Time_Bin[i], res1$ARI_Infomap_vs_Fenneman[i])
      add_summary("M1", "H2", "boundary enrichment on the Mississippi (>1 = boundaries sit on it)", res1$Time_Bin[i], res1$Boundary_enrichment_Mississippi[i])
      add_summary("M1", "H2", "boundary enrichment on the Rocky Mountain front", res1$Time_Bin[i], res1$Boundary_enrichment_RockyFront[i])
    }
    d <- cells1 %>% mutate(Bin = paste("Bin", Time_Bin), Region = factor(Infomap))
    p <- ggplot() + {if (!is.null(geo)) geom_sf(data = geo$study, aes(fill = Province), alpha = 0.18, colour = ink_soft, linewidth = 0.3)} +
      {if (!is.null(geo) && !is.null(geo$river)) geom_sf(data = geo$river, colour = "#2C6FB7", linewidth = 0.4)} +
      geom_point(data = d, aes(lon, lat, colour = Region), size = 2.6) +
      facet_wrap(~ Bin, nrow = 2) + scale_fill_manual(values = prov_cols, name = "Fenneman province") +
      scale_colour_manual(values = rep(c(okabe, "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#AA4499", "#882255"), 3),
                          name = "Infomap bioregion") +
      coord_sf(xlim = range(d$lon) + c(-2, 2), ylim = range(d$lat) + c(-2, 2)) +
      labs(x = NULL, y = NULL, caption = sprintf("Cells of %d km; Infomap on the bipartite cell-species network (igraph). Check with the Infomap Bioregions application (input files written per bin).", grid_km_net)) +
      theme_p() + theme(legend.box = "vertical")
    save_fig(p, dir1, "M1_fig1_bioregion_maps", 13, 8)
  }
}

# =============================================================================
# MODULE 2: DISTANCE DECAY
# =============================================================================

if (2 %in% run_modules && need("vegan", "Module 2")) {
  cat("\n=== MODULE 2: DISTANCE DECAY ===\n")
  dir2 <- mod_dir("M2_distance_decay")
  fit_decay <- function(sim, dist) {
    ok <- is.finite(sim) & is.finite(dist); sim <- sim[ok]; dist <- dist[ok]
    st <- coef(lm(log(sim + 0.01) ~ dist))
    f <- tryCatch(glm(sim ~ dist, family = quasibinomial(link = "log"), start = c(min(st[1], -0.01), min(st[2], -1e-6)),
                      control = glm.control(maxit = 100)), error = function(e) NULL)
    if (is.null(f)) return(c(intercept = NA, slope = NA))
    c(intercept = unname(coef(f)[1]), slope = unname(coef(f)[2]))
  }
  res2 <- list(); pairs2 <- list()
  for (gname in names(func_groups)) for (b in time_bins$Time_Bin) {
    cd <- cell_data(b, grid_km_net, func_groups[[gname]])
    if (is.null(cd) || nrow(cd$m) < 8) next          # fewer cells give unstable slopes
    S <- 1 - pair_beta_sim(cd$m); Dk <- gc_dist_km(cd$info$lon, cd$info$lat) / 1000   # 1,000 km
    ut <- upper.tri(S); est <- fit_decay(S[ut], Dk[ut])
    boots <- replicate(n_boot_decay, { k <- sample.int(nrow(S), replace = TRUE); Sk <- S[k, k]; Dkk <- Dk[k, k]
                                       u <- upper.tri(Sk) & Dkk > 0; fit_decay(Sk[u], Dkk[u])["slope"] })
    mt <- tryCatch(vegan::mantel(as.dist(1 - S), as.dist(Dk), method = "spearman", permutations = n_mantel, na.rm = TRUE),
                   error = function(e) NULL)
    res2[[length(res2) + 1]] <- data.frame(Group = gname, Time_Bin = b, n_cells = nrow(S), n_pairs = sum(ut),
      mean_similarity = mean(S[ut], na.rm = TRUE), slope_per_1000km = est["slope"],
      slope_q025 = qr95(boots)[1], slope_q975 = qr95(boots)[3],
      halving_distance_km = if (is.finite(est["slope"]) && est["slope"] < 0) 1000 * log(2) / -est["slope"] else NA,
      mantel_r = if (is.null(mt)) NA else mt$statistic, mantel_p = if (is.null(mt)) NA else mt$signif)
    pairs2[[length(pairs2) + 1]] <- data.frame(Group = gname, Time_Bin = b, dist_km = 1000 * Dk[ut], similarity = S[ut])
  }
  res2 <- bind_rows(res2); pairs2 <- bind_rows(pairs2)
  write.csv(res2, file.path(dir2, "M2_distance_decay.csv"), row.names = FALSE)
  if (nrow(res2)) {
    print(as.data.frame(res2 %>% transmute(Group, Bin = Time_Bin, cells = n_cells, slope = round(slope_per_1000km, 3),
                                           q025 = round(slope_q025, 3), q975 = round(slope_q975, 3), mantel_r = round(mantel_r, 2),
                                           mantel_p = round(mantel_p, 3))), row.names = FALSE)
    for (i in which(res2$Group == "All")) add_summary("M2", "H1/H2", "decay slope per 1,000 km (all taxa)", res2$Time_Bin[i], res2$slope_per_1000km[i])
    d <- pairs2 %>% filter(Group == "All") %>% mutate(Bin = paste("Bin", Time_Bin))
    fits <- res2 %>% filter(Group == "All") %>% rowwise() %>%
      do(data.frame(Time_Bin = .$Time_Bin, dist_km = seq(0, max(d$dist_km), length.out = 100))) %>% ungroup()
    p <- ggplot(d, aes(dist_km, similarity)) + geom_point(alpha = 0.25, size = 0.8, colour = ink_soft) +
      geom_smooth(method = "glm", formula = y ~ x, method.args = list(family = quasibinomial(link = "log"),
                  start = c(-0.1, -1e-4)), se = FALSE, colour = "#D55E00") +
      facet_wrap(~ Bin, nrow = 1) + labs(x = "Distance between cells (km)", y = expression("Similarity (1 - " * beta[sim] * ")"),
                                         caption = "Exponential decay (GLM, log link). Slopes, intervals and Mantel tests in M2_distance_decay.csv.") + theme_p()
    save_fig(p, dir2, "M2_fig1_distance_decay", 13, 4.4)
    d2 <- res2 %>% left_join(time_bins, by = "Time_Bin")
    p <- ggplot(d2, aes(Mid, slope_per_1000km, colour = Group)) + bands_layer() + geom_hline(yintercept = 0, colour = ink_soft) +
      geom_errorbar(aes(ymin = slope_q025, ymax = slope_q975), width = 0, position = position_dodge(0.15)) +
      geom_point(size = 2.2, position = position_dodge(0.15)) + scale_colour_manual(values = okabe, name = NULL) + x_scale +
      labs(x = "Age (Ma)", y = "Decay slope (per 1,000 km)", caption = "More negative = faster loss of similarity with distance. 95% cell-bootstrap ranges.") + theme_p()
    save_fig(p, dir2, "M2_fig2_decay_slopes", 10, 5.4)
  }
}

# =============================================================================
# MODULE 3: CLIMATE vs DISTANCE vs BARRIERS (GDM)
# =============================================================================

if (3 %in% run_modules && need("gdm", "Module 3") && need("vegan", "Module 3")) {
  cat("\n=== MODULE 3: GDM - CLIMATE vs DISTANCE vs BARRIERS ===\n")
  dir3 <- mod_dir("M3_gdm")
  have_terra <- requireNamespace("terra", quietly = TRUE)
  clim_idx <- if (file.exists(file.path(work_dir, climate_index_file))) read_wd(climate_index_file) else NULL
  if (is.null(clim_idx)) cat("  No climate_index.csv: models use geographic distance and barriers only (H2 side of the test).\n",
                             "  Supply climate / open-habitat rasters to test H1 (see the notes at the end).\n")
  if (!is.null(clim_idx) && !have_terra) { cat("  terra is not installed - climate rasters skipped.\n"); clim_idx <- NULL }

  # Site climate: mean of the layers whose age lies within the site's age range (nearest layer if none).
  site_env <- NULL
  if (!is.null(clim_idx)) {
    pts <- terra::vect(sites, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
    site_env <- data.frame(Site_Key = sites$Site_Key)
    for (v in unique(clim_idx$variable)) {
      iv <- clim_idx %>% filter(variable == v) %>% mutate(age_Ma = as.numeric(age_ka) / 1000)
      vals <- sapply(seq_len(nrow(iv)), function(k) terra::extract(terra::rast(file.path(work_dir, iv$file[k])), pts)[, 2])
      vals <- matrix(vals, nrow = nrow(sites))
      site_env[[v]] <- vapply(seq_len(nrow(sites)), function(i) {
        inr <- iv$age_Ma >= min(sites$Min_Ma[i], sites$Max_Ma[i]) & iv$age_Ma <= max(sites$Min_Ma[i], sites$Max_Ma[i])
        k <- if (any(inr)) which(inr) else which.min(abs(iv$age_Ma - (sites$Min_Ma[i] + sites$Max_Ma[i]) / 2))
        mean(vals[i, k], na.rm = TRUE)
      }, numeric(1))
    }
    write.csv(site_env, file.path(dir3, "M3_site_climate.csv"), row.names = FALSE)
  }
  rough <- NULL
  if (!is.na(dem_file) && have_terra && file.exists(dem_file)) {
    dem <- terra::rast(dem_file)
    tri <- terra::terrain(dem, v = "TRI")
    rough <- data.frame(Site_Key = sites$Site_Key,
                        roughness = terra::extract(tri, terra::vect(sites, geom = c("Longitude", "Latitude"), crs = "EPSG:4326"))[, 2])
  }
  cell_preds <- function(cd, b) {
    keys <- occ %>% filter(Time_Bin == b) %>% left_join(sites %>% select(Site_Key, X, Y), by = "Site_Key") %>%
      mutate(cell = cell_of(X, Y, grid_km_net)) %>% distinct(cell, Site_Key)
    env <- cd$info %>% select(cell, lon, lat, West_Miss, West_Rockies) %>%
      mutate(Mississippi_side = as.numeric(West_Miss), Rocky_side = as.numeric(West_Rockies))
    add <- list(site_env, rough); add <- add[!vapply(add, is.null, TRUE)]
    for (a in add) {
      e <- keys %>% left_join(a, by = "Site_Key") %>% group_by(cell) %>% summarise(across(-Site_Key, ~ mean(.x, na.rm = TRUE)), .groups = "drop")
      env <- env %>% left_join(e, by = "cell")
    }
    env
  }
  res3 <- list(); imp3 <- list(); vp3 <- list()
  groups3 <- list(All = NULL, Grazer = "Grazer", Browser = "Browser")
  for (gname in names(groups3)) for (b in time_bins$Time_Bin) {
    cd <- cell_data(b, grid_km_net, groups3[[gname]])
    if (is.null(cd) || nrow(cd$m) < min_cells_model) next
    env <- cell_preds(cd, b)
    pred_cols <- setdiff(names(env), c("cell", "lon", "lat", "West_Miss", "West_Rockies"))
    pred_cols <- pred_cols[vapply(env[pred_cols], function(v) all(is.finite(v)) && sd(v) > 0, TRUE)]
    env_df <- data.frame(site = seq_len(nrow(env)), x = env$lon, y = env$lat, env[pred_cols])
    for (resp in c("Sorensen", "pairwise beta_SIM")) {
      bio <- if (resp == "Sorensen") data.frame(site = seq_len(nrow(cd$m)), cd$m, check.names = FALSE) else
        { D <- pair_beta_sim(cd$m); D[is.na(D)] <- 1; data.frame(site = seq_len(nrow(D)), D) }
      sp <- tryCatch(gdm::formatsitepair(bio, bioFormat = if (resp == "Sorensen") 1 else 3, dist = "bray",
                                         abundance = FALSE, siteColumn = "site", XColumn = "x", YColumn = "y", predData = env_df),
                     error = function(e) { message("  formatsitepair failed: ", conditionMessage(e)); NULL })
      if (is.null(sp)) next
      fit <- tryCatch(gdm::gdm(sp, geo = TRUE), error = function(e) NULL)
      if (is.null(fit)) next
      spl <- tryCatch(gdm::isplineExtract(fit), error = function(e) NULL)
      heights <- if (is.null(spl)) NULL else apply(spl$y, 2, max)
      res3[[length(res3) + 1]] <- data.frame(Group = gname, Time_Bin = b, Response = resp, n_cells = nrow(cd$m),
                                             deviance_explained = fit$explained, predictor = names(heights), spline_height = unname(heights))
      if (resp == "Sorensen" && gdm_perm > 0) {
        vi <- tryCatch(suppressWarnings(gdm::gdm.varImp(sp, geo = TRUE, nPerm = gdm_perm, parallel = FALSE, predSelect = TRUE)),
                       error = function(e) NULL)
        if (is.null(vi))            # backward elimination needs several predictors; fall back to the full model only
          vi <- tryCatch(suppressWarnings(gdm::gdm.varImp(sp, geo = TRUE, nPerm = gdm_perm, parallel = FALSE, predSelect = FALSE)),
                         error = function(e) { message("  gdm.varImp (", gname, ", Bin ", b, "): ", conditionMessage(e)); NULL })
        if (!is.null(vi)) {
          imp <- as.data.frame(vi[[2]]); sig <- as.data.frame(vi[[3]])
          imp3[[length(imp3) + 1]] <- data.frame(Group = gname, Time_Bin = b, predictor = rownames(imp),
                                                 importance_full_model = imp[, 1], p_full_model = sig[, 1])
        }
      }
    }
    # variation partitioning: environment vs space vs barriers (dbRDA on Sorensen)
    Ydist <- vegan::vegdist(cd$m, method = "bray", binary = TRUE)
    space <- tryCatch({ pc <- vegan::pcnm(dist(cbind(env$lon, env$lat))); as.data.frame(pc$vectors[, seq_len(min(3, ncol(pc$vectors))), drop = FALSE]) },
                      error = function(e) NULL)
    barr <- env[, c("Mississippi_side", "Rocky_side")]; barr <- barr[, vapply(barr, function(v) sd(v) > 0, TRUE), drop = FALSE]
    envv <- env[, setdiff(pred_cols, c("Mississippi_side", "Rocky_side")), drop = FALSE]
    sets <- list(Environment = envv, Space = space, Barriers = barr)
    sets <- sets[vapply(sets, function(s) !is.null(s) && ncol(s) > 0, TRUE)]
    if (length(sets) >= 2) {
      vp <- tryCatch(do.call(vegan::varpart, c(list(Ydist), unname(sets))), error = function(e) NULL)
      if (!is.null(vp)) {
        fr <- vp$part$indfract; lab <- rownames(fr)
        for (k in seq_along(sets)) lab <- gsub(paste0("X", k), names(sets)[k], lab)
        vp3[[length(vp3) + 1]] <- data.frame(Group = gname, Time_Bin = b, Sets = paste(names(sets), collapse = " | "),
                                             Fraction = lab, Adj_R2 = fr$Adj.R.square)
      }
    }
  }
  res3 <- bind_rows(res3); imp3 <- bind_rows(imp3); vp3 <- bind_rows(vp3)
  write.csv(res3, file.path(dir3, "M3_gdm_splines.csv"), row.names = FALSE)
  write.csv(imp3, file.path(dir3, "M3_gdm_importance_permutation.csv"), row.names = FALSE)
  write.csv(vp3, file.path(dir3, "M3_variation_partitioning.csv"), row.names = FALSE)
  if (nrow(res3)) {
    print(as.data.frame(res3 %>% group_by(Group, Time_Bin, Response) %>% summarise(dev_explained = round(first(deviance_explained), 1),
          top = predictor[which.max(spline_height)], .groups = "drop")), row.names = FALSE)
    s <- res3 %>% filter(Group == "All", Response == "Sorensen") %>% group_by(Time_Bin) %>%
      mutate(share = spline_height / sum(spline_height)) %>% ungroup()
    for (i in seq_len(nrow(s))) add_summary("M3", "H1/H2", paste("share of GDM spline height:", s$predictor[i]), s$Time_Bin[i], s$share[i],
                                            sprintf("deviance explained %.1f%%", s$deviance_explained[i]))
    p <- ggplot(res3 %>% left_join(time_bins, by = "Time_Bin") %>% mutate(Bin = factor(paste("Bin", Time_Bin))),
                aes(Bin, spline_height, fill = predictor)) + geom_col(colour = ink_soft, linewidth = 0.2) +
      facet_grid(Response ~ Group) + scale_fill_manual(values = okabe, name = NULL) +
      labs(x = NULL, y = "GDM I-spline height (turnover explained)",
           caption = paste0("Ferrier et al. (2007) GDM; barriers coded as 0/1 sides (Mississippi; Basin and Range = west of the Rocky Mountain front).",
                            if (is.null(clim_idx)) "\nNo climate rasters supplied: environment (H1) not yet in the models." else "")) + theme_p()
    save_fig(p, dir3, "M3_fig1_gdm_predictor_heights", 12, 7)
  }
}

# =============================================================================
# MODULE 4: RAUP-CRICK (H5)
# =============================================================================

if (4 %in% run_modules) {
  cat("\n=== MODULE 4: RAUP-CRICK ===\n")
  dir4 <- mod_dir("M4_raup_crick")
  rc_pairs <- function(m) {            # m: units x species 0/1; occupancy-weighted null (Chase et al. 2011)
    occw <- colSums(m) / sum(colSums(m)); spp <- colnames(m); rich <- rowSums(m)
    cache <- list(); out <- matrix(NA, nrow(m), nrow(m), dimnames = list(rownames(m), rownames(m)))
    for (i in 1:(nrow(m) - 1)) for (j in (i + 1):nrow(m)) {
      key <- paste(sort(c(rich[i], rich[j])), collapse = "_")
      if (is.null(cache[[key]])) cache[[key]] <- replicate(n_rc, length(intersect(sample(spp, rich[i], prob = occw),
                                                                               sample(spp, rich[j], prob = occw))))
      nul <- cache[[key]]; obs <- sum(m[i, ] & m[j, ])
      out[i, j] <- out[j, i] <- 2 * ((sum(nul > obs) + 0.5 * sum(nul == obs)) / n_rc) - 1
    }
    out
  }
  res4 <- list()
  for (gname in names(func_groups)) for (b in time_bins$Time_Bin) {
    cd <- cell_data(b, grid_km_net, func_groups[[gname]])
    if (!is.null(cd) && nrow(cd$m) >= 3) {
      rc <- rc_pairs(cd$m); v <- rc[upper.tri(rc)]
      res4[[length(res4) + 1]] <- data.frame(Group = gname, Time_Bin = b, Unit = "cell pairs", n_pairs = length(v),
        mean_RC = mean(v), share_divergent = mean(v > 0.95), share_stochastic = mean(abs(v) <= 0.95), share_convergent = mean(v < -0.95))
    }
    d <- occ %>% filter(Time_Bin == b, in_group(Species, func_groups[[gname]])) %>% distinct(Spatial_Bin, Species)
    if (n_distinct(d$Spatial_Bin) >= 2 && n_distinct(d$Species) >= 3) {
      rc <- rc_pairs(lists_to_matrix(split(d$Species, d$Spatial_Bin))); v <- rc[upper.tri(rc)]
      nm <- outer(rownames(rc), rownames(rc), paste, sep = " vs ")[upper.tri(rc)]
      res4[[length(res4) + 1]] <- data.frame(Group = gname, Time_Bin = b, Unit = paste("province pair:", nm), n_pairs = 1,
        mean_RC = v, share_divergent = as.numeric(v > 0.95), share_stochastic = as.numeric(abs(v) <= 0.95), share_convergent = as.numeric(v < -0.95))
    }
  }
  res4 <- bind_rows(res4)
  write.csv(res4, file.path(dir4, "M4_raup_crick.csv"), row.names = FALSE)
  if (nrow(res4)) {
    print(as.data.frame(res4 %>% filter(Group == "All") %>% mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
    s <- res4 %>% filter(Group == "All", Unit == "cell pairs")
    for (i in seq_len(nrow(s))) add_summary("M4", "H5", "share of cell pairs deterministic-divergent (RC > 0.95)", s$Time_Bin[i], s$share_divergent[i],
                                            sprintf("mean RC %.2f", s$mean_RC[i]))
    d <- res4 %>% filter(Unit == "cell pairs") %>% select(Group, Time_Bin, share_divergent, share_stochastic, share_convergent) %>%
      pivot_longer(starts_with("share"), names_to = "Class", values_to = "Share") %>%
      mutate(Class = factor(sub("share_", "", Class), levels = c("convergent", "stochastic", "divergent")), Bin = factor(paste("Bin", Time_Bin)))
    p <- ggplot(d, aes(Bin, Share, fill = Class)) + geom_col(colour = ink_soft, linewidth = 0.2) + facet_wrap(~ Group, nrow = 2) +
      scale_fill_manual(values = c(convergent = "#56B4E9", stochastic = "#F6EBD0", divergent = "#D55E00"), name = NULL) +
      labs(x = NULL, y = "Share of cell pairs", caption = "Raup-Crick (Chase et al. 2011), occupancy-weighted null; divergent > 0.95, convergent < -0.95.") + theme_p()
    save_fig(p, dir4, "M4_fig1_raup_crick", 12, 7)
  }
}

# =============================================================================
# MODULE 6: WHERE DID THE ENDEMICS COME FROM? (H4)
# =============================================================================

if (6 %in% run_modules) {
  cat("\n=== MODULE 6: ORIGIN OF ENDEMICS ===\n")
  dir6 <- mod_dir("M6_endemics")
  origin_tab <- if (file.exists(file.path(work_dir, origins_file))) read_wd(origins_file) else NULL
  if (is.null(origin_tab)) cat("  No taxon_origins.csv (from ch2_immigrants.R): immigrant endemics cannot be identified; all count as native.\n")
  origin_of <- function(sp) {
    if (is.null(origin_tab)) return(rep("Native", length(sp)))
    o <- setNames(origin_tab$Origin, origin_tab$Taxon)
    r <- ifelse(!is.na(o[sp]), o[sp], o[sub(" .*$", "", sp)]); ifelse(is.na(r), "Native", r)
  }
  pres <- occ %>% distinct(Time_Bin, Spatial_Bin, Species)
  first_na <- occ_na %>% filter(!is.na(Time_Bin)) %>% group_by(Species) %>% summarise(first_bin = min(Time_Bin), .groups = "drop")
  first_where <- occ_na %>% filter(!is.na(Time_Bin)) %>% inner_join(first_na, by = "Species") %>% filter(Time_Bin == first_bin) %>%
    group_by(Species) %>% summarise(first_provinces = paste(sort(unique(Spatial_Bin)), collapse = "; "), .groups = "drop")
  end6 <- pres %>% group_by(Time_Bin, Species) %>% filter(n() == 1) %>% ungroup() %>% rename(Home = Spatial_Bin) %>%
    left_join(first_na, by = "Species") %>% left_join(first_where, by = "Species") %>%
    rowwise() %>% mutate(
      earlier_multi = any(sapply(time_bins$Time_Bin[time_bins$Time_Bin < Time_Bin], function(e)
        sum(pres$Species == Species & pres$Time_Bin == e) >= 2)),
      Origin = origin_of(Species)) %>% ungroup() %>%
    mutate(Class = case_when(
      Origin %in% c("South American", "Eurasian") ~ "Immigrant endemic",
      earlier_multi ~ "Relict",
      first_bin == Time_Bin & first_provinces == Home & Time_Bin > min(time_bins$Time_Bin) ~ "Originated in place",
      first_bin == min(time_bins$Time_Bin) & Time_Bin == min(time_bins$Time_Bin) ~ "Present at window start (origin unknown)",
      TRUE ~ "Other (first recorded elsewhere or earlier)"))
  write.csv(end6, file.path(dir6, "M6_endemics_classified.csv"), row.names = FALSE)
  shares <- end6 %>% count(Time_Bin, Class) %>% group_by(Time_Bin) %>% mutate(share = n / sum(n)) %>% ungroup()
  write.csv(shares, file.path(dir6, "M6_endemic_class_shares.csv"), row.names = FALSE)
  print(as.data.frame(shares %>% mutate(share = round(share, 2))), row.names = FALSE)
  for (i in seq_len(nrow(shares))) add_summary("M6", "H4/H1/H3", paste("share of endemics:", shares$Class[i]), shares$Time_Bin[i], shares$share[i])
  p <- ggplot(shares %>% mutate(Bin = factor(paste("Bin", Time_Bin))), aes(Bin, share, fill = Class)) +
    geom_col(colour = ink_soft, linewidth = 0.2) +
    scale_fill_manual(values = c("Originated in place" = "#009E73", "Immigrant endemic" = "#E69F00", "Relict" = "#0072B2",
                                 "Present at window start (origin unknown)" = "#F6EBD0", "Other (first recorded elsewhere or earlier)" = "#CC79A7"), name = NULL) +
    labs(x = NULL, y = "Share of single-province species",
         caption = "Extends Darroch et al. (2014). First records use all North American occurrences in the master data; Bin 1 origins are unknown (left-censored).") +
    theme_p() + guides(fill = guide_legend(ncol = 2))
  save_fig(p, dir6, "M6_fig1_endemic_origins", 10, 6)

  # PyRate input per province, genus level (Silvestro et al. 2014; run PyRate separately)
  for (p in regions) {
    d <- occ %>% filter(Spatial_Bin == p, !is.na(Max_Ma), !is.na(Min_Ma)) %>%
      transmute(Species = Genus, Status = "extinct", min_age = pmin(Min_Ma, Max_Ma), max_age = pmax(Min_Ma, Max_Ma))
    write.table(d, file.path(dir6, paste0("M6_PyRate_input_", gsub(" ", "_", p), "_genus.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  cat("  PyRate inputs written (status set to 'extinct' for all genera - edit extant genera before running PyRate).\n")
}

# =============================================================================
# SYNTHESIS
# =============================================================================

synth <- bind_rows(summary_rows)
ch2_dir <- file.path(work_dir, "Outputs", "ch2_immigrants", "species")
if (file.exists(file.path(ch2_dir, "ch2_beta_by_origin_scenario.csv"))) {
  c2 <- read.csv(file.path(ch2_dir, "ch2_beta_by_origin_scenario.csv")) %>% filter(Scenario == "Natives only")
  synth <- bind_rows(synth, data.frame(Module = "M5", Hypothesis = "H3", Metric = "change in beta_SIM when immigrants are removed (ch2)",
                                       Time_Bin = c2$Time_Bin, Value = c2$delta, Note = sprintf("p vs random natives = %.3f", c2$p_vs_random_natives)))
}
if (nrow(synth)) {
  synth <- synth %>% arrange(Module, Metric, Time_Bin)
  write.csv(synth, file.path(output_root, "SYNTHESIS_key_results.csv"), row.names = FALSE)
  cat("\n=== SYNTHESIS: key results per module (see SYNTHESIS_key_results.csv) ===\n")
  print(as.data.frame(synth %>% mutate(Value = round(Value, 3))), row.names = FALSE)
}

notes <- c(
  "WHAT MAINTAINED PROVINCIALISM - NOTES",
  "Read the modules against the plan's 'Reading the result' table; report each cause's share per bin.",
  "Module 0 is the gate: use its corrected curve, not the Chapter 1 curve, when interpreting the others.",
  "Province-bins below the cell quota are reported as insufficient; lower quotas are a sensitivity check only.",
  "Identifications check: 'names as identified only' vs all names. The full Rowan four-matrix cf./aff. check",
  "  needs Step 3d to write the alternative matrices.",
  "Module 1 uses igraph's Infomap on the bipartite network; confirm with the Infomap Bioregions app (input files written).",
  "Module 3 without climate_index.csv tests only distance and barriers (H2). Oscillayers scales one spatial pattern",
  "  through time: treat early-bin climate results as provisional and compare with PaleoClim snapshots.",
  "Modern topography is assumed for roughness (relief of the West largely in place by 3.25 Ma).",
  "Not implemented yet: least-cost distance in GDM (Ferrier et al. 2007 Fig. 3), spatially constrained clustering",
  "  (Brook et al. 2025), and the PyRate runs themselves (input files are written).",
  "Module 5 = ch2_immigrants.R; its results are added to the synthesis if they exist.",
  "Module 6: Bin 1 endemics cannot be classified by origin (data start at 3.25 Ma)."
)
writeLines(notes, file.path(output_root, "README_provincialism.txt"))
cat("\n", paste(notes, collapse = "\n"), "\n\n=== DONE: outputs in", output_root, "===\n")
