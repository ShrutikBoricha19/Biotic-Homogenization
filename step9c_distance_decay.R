# =============================================================================
# STEP 9c: DISTANCE DECAY OF SIMILARITY
#
# For every pair of sites in a stage, Simpson similarity (1 - pairwise
# beta_SIM) is plotted against the distance between the two sites.
#   - Heterogeneous landscape: similarity drops quickly with distance (steep
#     negative slope); nearby sites share species, distant ones do not.
#   - Homogenization: the slope flattens and the line rises over time; distant
#     sites become as similar as nearby ones.
#
#   a) All sites of the four provinces together, one panel per stage, with a
#      linear fit and a Mantel test (permutation test of the correlation
#      between the similarity and distance matrices).
#   b) The slope (similarity change per 1000 km) and the mean similarity
#      through time, for all four provinces together and for each province
#      on its own (pairs within that province only).
#
# Each site is one point per stage (coordinates averaged over its records).
# Sites without coordinates are left out.
#
# Inputs:  all_records_final.csv (Step 6), site_physio_database.csv (Step 5)
# Outputs (Outputs/9c_distance_decay/):
#   distance_decay_by_stage.png/.pdf    a)
#   distance_decay_trend.png/.pdf       b)
#   distance_decay.csv                  slope, intercept, Mantel r and p
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
output_dir   <- file.path(work_dir, "Outputs", "9c_distance_decay")

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

min_sites_decay <- 5      # a province x stage needs at least this many sites for a slope
n_perm          <- 999    # Mantel test permutations
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
# 2. SITE PAIRS: SIMILARITY AND DISTANCE
# =============================================================================

cat("\n=== 2. SITE PAIRS ===\n")

haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}

sites <- recs %>%
  group_by(Stage_Number, Province, Site) %>%
  summarise(Latitude = mean(Latitude, na.rm = TRUE), Longitude = mean(Longitude, na.rm = TRUE),
            species = list(unique(GenusSpecies)), .groups = "drop") %>%
  mutate(has_xy = is.finite(Latitude) & is.finite(Longitude))
cat(sprintf("  %d site x stage units; %d without coordinates are left out\n",
            nrow(sites), sum(!sites$has_xy)))
sites <- filter(sites, has_xy)

# Similarity and distance matrices for one set of sites.
pair_matrices <- function(st) {
  n <- nrow(st)
  S <- matrix(NA_real_, n, n); D <- matrix(0, n, n)
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    S[i, j] <- S[j, i] <- 1 - beta_sim_pair(st$species[[i]], st$species[[j]])
    D[i, j] <- D[j, i] <- haversine_km(st$Latitude[i], st$Longitude[i], st$Latitude[j], st$Longitude[j])
  }
  list(S = S, D = D)
}

# Linear fit (similarity ~ distance in 1000 km) and Mantel test.
decay_stats <- function(st) {
  if (nrow(st) < 3) return(NULL)
  m <- pair_matrices(st)
  lt <- lower.tri(m$S)
  sim <- m$S[lt]; dist <- m$D[lt] / 1000
  if (length(unique(dist)) < 2) return(NULL)
  fit <- lm(sim ~ dist)
  r_obs <- suppressWarnings(cor(sim, dist))
  r_null <- replicate(n_perm, {
    o <- sample.int(nrow(st))
    suppressWarnings(cor(m$S[o, o][lt], dist))
  })
  list(pairs = data.frame(similarity = sim, distance_km = dist * 1000),
       stats = data.frame(n_sites = nrow(st), n_pairs = length(sim),
                          mean_similarity = mean(sim),
                          intercept = unname(coef(fit)[1]), slope_per_1000km = unname(coef(fit)[2]),
                          mantel_r = r_obs,
                          mantel_p = (sum(r_null <= r_obs, na.rm = TRUE) + 1) / (n_perm + 1)))
}

all_pairs <- list(); stats <- list()
for (s in 1:5) {
  st_all <- filter(sites, Stage_Number == s)
  res <- decay_stats(st_all)
  if (!is.null(res)) {
    # Label each pair: same province or different provinces.
    lt <- which(lower.tri(matrix(0, nrow(st_all), nrow(st_all))), arr.ind = TRUE)
    res$pairs$Pair_type <- ifelse(st_all$Province[lt[, 1]] == st_all$Province[lt[, 2]],
                                  "Same province", "Different provinces")
    all_pairs[[s]] <- mutate(res$pairs, Stage_Number = s)
    stats[[length(stats) + 1]] <- mutate(res$stats, Group = "All four provinces", Stage_Number = s)
  }
  for (p in provinces) {
    st_p <- filter(st_all, Province == p)
    if (nrow(st_p) < min_sites_decay) next
    res_p <- decay_stats(st_p)
    if (!is.null(res_p)) stats[[length(stats) + 1]] <- mutate(res_p$stats, Group = p, Stage_Number = s)
  }
  cat(sprintf("  %-10s %4d sites with coordinates\n", stage_names[s], nrow(st_all)))
}
all_pairs <- bind_rows(all_pairs) %>%
  mutate(Stage = factor(stage_names[Stage_Number], levels = stage_names))
stats <- bind_rows(stats) %>%
  mutate(Stage = stage_names[Stage_Number], Mid_Ma = stage_mid[Stage_Number]) %>%
  relocate(Group, Stage_Number, Stage)

save_csv(mutate(stats, across(where(is.numeric), ~ round(.x, 4))), "distance_decay")
cat("\n  All four provinces together:\n")
print(as.data.frame(stats %>% filter(Group == "All four provinces") %>%
        transmute(Stage, n_sites, mean_similarity = round(mean_similarity, 3),
                  slope_per_1000km = round(slope_per_1000km, 3), mantel_r = round(mantel_r, 3),
                  mantel_p = round(mantel_p, 3))), row.names = FALSE)

# =============================================================================
# 3. FIGURES
# =============================================================================

cat("\n=== 3. FIGURES ===\n")

lab_a <- stats %>% filter(Group == "All four provinces") %>%
  mutate(Stage = factor(Stage, levels = stage_names),
         label = sprintf("slope %.2f / 1000 km\nMantel r = %.2f, p %s", slope_per_1000km, mantel_r,
                         ifelse(mantel_p < 0.001, "< 0.001", sprintf("= %.3f", mantel_p))))
pair_cols <- c("Same province" = "#2a78d6", "Different provinces" = "#eb6834")

p_a <- ggplot(all_pairs, aes(distance_km, similarity)) +
  geom_point(aes(colour = Pair_type), alpha = 0.25, size = 1.3, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, colour = ink, linewidth = 1.1, se = FALSE) +
  geom_text(data = lab_a, aes(x = Inf, y = 1.17, label = label), hjust = 1.03, vjust = 1,
            size = base_size / 5, colour = ink, lineheight = 0.95) +
  facet_wrap(~ Stage, nrow = 1, drop = FALSE) +
  scale_colour_manual(values = pair_cols, name = NULL) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1.18)) +
  scale_x_continuous(labels = function(x) paste0(x / 1000, "k")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  labs(title = "Does similarity between sites fall with distance?",
       subtitle = "Every pair of sites in each stage (all four provinces); flatter, higher line through time = homogenization",
       x = "Distance between sites (km)", y = "Simpson similarity (1 - beta_SIM)",
       caption = paste0("Line = linear fit. Mantel test with ", n_perm,
                        " permutations (one-sided: similarity decreasing with distance).")) +
  theme_slide() +
  theme(legend.position = "top", legend.justification = "left",
        panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5),
        panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.4),
        axis.text = element_text(size = base_size * 0.6))
save_fig(p_a, "distance_decay_by_stage")

group_cols <- c("All four provinces" = ink, province_colours)
trend_long <- stats %>%
  select(Group, Stage_Number, Mid_Ma, mean_similarity, slope_per_1000km) %>%
  pivot_longer(c(mean_similarity, slope_per_1000km), names_to = "Measure", values_to = "value") %>%
  mutate(Measure = factor(ifelse(Measure == "mean_similarity", "Mean similarity between sites",
                                 "Slope: similarity per 1000 km"),
                          levels = c("Mean similarity between sites",
                                     "Slope: similarity per 1000 km")),
         Group = factor(Group, levels = names(group_cols)),
         Series = interaction(Group, Measure))

p_b <- ggplot(trend_long) +
  stage_bands() +
  geom_hline(data = data.frame(Measure = factor(levels(trend_long$Measure)[2],
                                                levels = levels(trend_long$Measure)), y = 0),
             aes(yintercept = y), colour = ink_soft, linetype = "dashed") +
  geom_segment(data = gap_segments(trend_long, "Series"),
               aes(x = Mid_Ma, xend = x2, y = value, yend = y2, colour = Group,
                   linewidth = Group == "All four provinces")) +
  geom_point(aes(Mid_Ma, value, colour = Group), shape = 21, fill = "white", size = 3.4, stroke = 1.5) +
  facet_wrap(~ Measure, scales = "free_y") +
  scale_colour_manual(values = group_cols, name = NULL) +
  scale_linewidth_manual(values = c(`TRUE` = 1.8, `FALSE` = 1), guide = "none") +
  age_axis() +
  labs(title = "Is the landscape becoming more uniform?",
       subtitle = "Rising mean similarity and a slope moving towards 0 through time = homogenization",
       x = "Age (Ma)", y = NULL,
       caption = paste0("Province lines use pairs of sites within that province only (at least ", min_sites_decay,
                        " sites); the black line uses all pairs across the four provinces.\n",
                        "Bands = stages, Zanclean (left) to Chibanian (right).")) +
  theme_slide() +
  theme(legend.position = "top", legend.justification = "left",
        panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5))
save_fig(p_b, "distance_decay_trend")
cat("  Figures saved\n")

cat("\n=== STEP 9c COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: sites, all_pairs, stats\n")
