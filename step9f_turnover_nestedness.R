# =============================================================================
# STEP 9f: TURNOVER VERSUS NESTEDNESS
#
# beta_SIM measures only one way in which sites can differ: species
# REPLACEMENT (turnover). Total dissimilarity (Sorensen, beta_SOR) also
# includes NESTEDNESS (beta_SNE): poorer sites holding a subset of the
# species of richer sites (Baselga 2010):
#
#   beta_SOR = [sum min(b_ij, b_ji) + sum max(b_ij, b_ji)] /
#              [2 (sum_i S_i - S_T) + sum min + sum max]
#   beta_SIM = sum min / [(sum_i S_i - S_T) + sum min]      (turnover)
#   beta_SNE = beta_SOR - beta_SIM                         (nestedness)
#
# How sites can become more alike:
#   - turnover falls: the same species occur everywhere (classic homogenization)
#   - nestedness rises while turnover falls: some sites lose species and become
#     subsets of richer sites (homogenization by local loss)
#
# Computed for each province and stage on the same subsets of sites as Step 7,
# so beta_SIM here equals Step 7.
#
# Inputs:  all_records_final.csv (Step 6), site_physio_database.csv (Step 5)
# Outputs (Outputs/9f_turnover_nestedness/):
#   turnover_nestedness_all_provinces.png/.pdf
#   turnover_nestedness.csv
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Steps 5-7.
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
output_dir   <- file.path(work_dir, "Outputs", "9f_turnover_nestedness")

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

# Keep these three identical to Step 7 so beta_SIM matches Step 7.
n_units_sample <- 3
n_resamples    <- 999
step7_file     <- file.path("Outputs", "7_regional_pools", "beta_spatial.csv")

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
# 2. TURNOVER AND NESTEDNESS PER PROVINCE AND STAGE
# =============================================================================

cat("\n=== 2. TURNOVER AND NESTEDNESS ===\n")

unit_subsets <- function(k, n) {
  if (k < n || n < 2) return(list())
  if (choose(k, n) <= n_resamples) combn(k, n, simplify = FALSE) else
    replicate(n_resamples, sample.int(k, n), simplify = FALSE)
}

set.seed(2024)   # same seed and loop order as Step 7 -> same site subsets
tn <- bind_rows(lapply(provinces, function(p) bind_rows(lapply(1:5, function(s) {
  d <- filter(recs, Province == p, Stage_Number == s)
  units <- lapply(split(d$GenusSpecies, d$Site), unique)
  subsets <- unit_subsets(length(units), n_units_sample)
  out <- data.frame(Province = p, Stage_Number = s, Stage = stage_names[s], Mid_Ma = stage_mid[s],
                    n_sites = length(units), beta_SIM = NA_real_, beta_SNE = NA_real_,
                    beta_SOR = NA_real_, beta_SOR_lo95 = NA_real_, beta_SOR_hi95 = NA_real_)
  if (length(subsets) == 0) return(out)
  v <- sapply(subsets, function(ix) beta_multi_all(units[ix])[c("beta_SIM", "beta_SNE", "beta_SOR")])
  out$beta_SIM <- mean(v["beta_SIM", ]); out$beta_SNE <- mean(v["beta_SNE", ])
  out$beta_SOR <- mean(v["beta_SOR", ])
  out$beta_SOR_lo95 <- unname(quantile(v["beta_SOR", ], 0.025))
  out$beta_SOR_hi95 <- unname(quantile(v["beta_SOR", ], 0.975))
  out
})))) %>%
  mutate(Province = factor(Province, levels = provinces),
         nestedness_share = beta_SNE / beta_SOR)
save_csv(mutate(tn, across(where(is.numeric), ~ round(.x, 4))), "turnover_nestedness")
print(as.data.frame(tn %>% transmute(Province, Stage, n_sites, turnover = round(beta_SIM, 3),
                                     nestedness = round(beta_SNE, 3), total = round(beta_SOR, 3))),
      row.names = FALSE)

s7_path <- file.path(work_dir, step7_file)
if (file.exists(s7_path)) {
  s7 <- read.csv(s7_path, stringsAsFactors = FALSE, check.names = FALSE)
  chk <- tn %>% mutate(Province = as.character(Province)) %>%
    inner_join(transmute(s7, Province, Stage, b7 = beta_SIM), by = c("Province", "Stage"))
  dif <- suppressWarnings(max(abs(chk$beta_SIM - chk$b7), na.rm = TRUE))
  cat(sprintf("\n  Check against Step 7: largest difference in beta_SIM = %s\n",
              if (is.finite(dif)) format(round(dif, 4), nsmall = 4) else "n/a"))
}

# =============================================================================
# 3. FIGURE
# =============================================================================

cat("\n=== 3. FIGURE ===\n")

comp_cols <- c("Turnover (beta_SIM)" = "#256abf", "Nestedness (beta_SNE)" = "#9ec5f4")
bars <- tn %>%
  select(Province, Stage_Number, beta_SIM, beta_SNE) %>%
  pivot_longer(c(beta_SIM, beta_SNE), names_to = "Component", values_to = "value") %>%
  mutate(Component = factor(ifelse(Component == "beta_SIM", names(comp_cols)[1], names(comp_cols)[2]),
                            levels = rev(names(comp_cols))),
         xmin = stage_young[Stage_Number] + 0.06, xmax = stage_older[Stage_Number] - 0.06) %>%
  filter(!is.na(value)) %>%
  arrange(Province, Stage_Number, Component) %>%
  group_by(Province, Stage_Number) %>%
  mutate(ymax = cumsum(value), ymin = ymax - value) %>%
  ungroup()
tot <- tn %>% filter(!is.na(beta_SOR)) %>% mutate(label = sprintf("%.2f", beta_SOR))
none <- tn %>% filter(is.na(beta_SOR)) %>% mutate(label = "too few\nsites")

p_tn <- ggplot(bars) +
  stage_bands() +
  geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), colour = grid_col, linewidth = 0.4) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Component),
            colour = "white", linewidth = 0.6) +
  geom_text(data = tot, aes(x = Mid_Ma, y = beta_SOR + 0.03, label = label), vjust = 0,
            size = base_size / 4.4, fontface = "bold", colour = ink) +
  geom_text(data = none, aes(x = Mid_Ma, y = 0.12, label = label), size = base_size / 5.5,
            colour = ink_soft, lineheight = 0.9) +
  facet_wrap(~ Province, ncol = 2, drop = FALSE) +
  scale_fill_manual(values = c(comp_cols, setNames(band_cols, band_cols)),
                    breaks = names(comp_cols), name = NULL) +
  age_axis() +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "How do sites differ: replacement or loss?",
       subtitle = "Multisite Sorensen dissimilarity among sites split into turnover and nestedness (number = total)",
       x = "Age (Ma)", y = "Dissimilarity",
       caption = paste0("Turnover = species replaced between sites (the beta_SIM of Step 7); nestedness = poorer sites hold a ",
                        "subset of richer sites' species.\nMean of subsets of ", n_units_sample,
                        " sites (Baselga 2010). Falling turnover = homogenization; rising nestedness = sites losing species.")) +
  theme_slide() +
  theme(legend.position = "top", legend.justification = "left",
        panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5),
        panel.grid.major.y = element_blank())
save_fig(p_tn, "turnover_nestedness_all_provinces", height = slide_h * 1.1)
cat("  Figure saved\n")

cat("\n=== STEP 9f COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: tn\n")
