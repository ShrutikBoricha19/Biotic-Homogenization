# =============================================================================
# STEP 7: REGIONAL SPECIES POOLS AND TEMPORAL SIMPSON DISSIMILARITY
#
# For four physiographic provinces (Fenneman & Johnson 1946):
#   Pacific Border, Basin and Range, Great Plains, Coastal Plain
#
# 1. Regional species pool: every species recorded at any site in the province
#    during a geological stage (FAUNMAP + PBDB, from the Step 6 matrices).
# 2. Column graph of the number of sites per stage in each province.
# 3. Simpson dissimilarity (beta_SIM) of each province's pool WITH ITSELF
#    THROUGH TIME (Zanclean v Piacenzian v Gelasian v Calabrian v Chibanian),
#    using the formulas of Rowan et al. (2024, Nat. Ecol. Evol.), after
#    Baselga (2010) and Baselga et al. (2007):
#
#      multisite:  beta_SIM = sum_{i<j} min(b_ij, b_ji) /
#                  [ sum_i S_i - S_T + sum_{i<j} min(b_ij, b_ji) ]     (Eq. 1)
#      pairwise:   beta_SIM = min(b, c) / (a + min(b, c))
#      endemic and shared components of the multisite value: Eqs. 2 and 3.
#
#    Here i, j are STAGES (not subregions): b_ij = species in stage i's pool
#    but not stage j's, S_i = pool size, S_T = all species across stages.
#    0 = identical pools (or one nested in the other); 1 = no species shared.
#
# Inputs:
#   Outputs/6_matrices/all_records_final.csv              (Step 6)
#   Outputs/maps/physiographic/site_physio_database.csv   (Step 5)
#
# Outputs (Outputs/7_regional_pools/):
#   Figures (PNG 300 dpi + PDF, 16:9 slide format):
#     sites_per_stage_<province>.png/.pdf      one column graph per province
#     sites_per_stage_all_provinces.png/.pdf   the four together
#     beta_pairwise_<province>.png/.pdf        stage-by-stage beta_SIM matrix
#     beta_pairwise_all_provinces.png/.pdf     the four together
#     beta_consecutive_all_provinces.png/.pdf  stage-to-stage turnover through time
#     beta_multisite_all_provinces.png/.pdf    overall beta_SIM across all stages
#   Tables:
#     regional_pools.xlsx            one sheet per province: species x stage (1/0)
#     regional_pools_long.csv        province, stage, species, order, n sites
#     pool_summary.csv               sites and species per province and stage
#     beta_pairwise.csv              a, b, c and beta_SIM for every stage pair
#     beta_consecutive.csv           beta_SIM between successive stages
#     beta_multisite.csv             multisite beta_SIM, endemic and shared parts
#     site_province_conflicts.csv    sites whose analysis units fall in >1 province
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Steps 5 and 6.
# Needs: dplyr, tidyr, ggplot2, writexl.
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "writexl")) {
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
output_dir   <- file.path(work_dir, "Outputs", "7_regional_pools")

provinces <- c("Pacific Border", "Basin and Range", "Great Plains", "Coastal Plain")

# One colour per province (validated colour-blind-safe set).
province_colours <- c(
  "Pacific Border"  = "#2a78d6",   # blue
  "Basin and Range" = "#eb6834",   # orange
  "Great Plains"    = "#1baf7a",   # aqua
  "Coastal Plain"   = "#4a3aa7"    # violet
)

stage_names <- c("Zanclean", "Piacenzian", "Gelasian", "Calabrian", "Chibanian")
stage_older <- c(4.700, 3.600, 2.580, 1.800, 0.7741)
stage_young <- c(3.600, 2.580, 1.800, 0.7741, 0.129)

# Figure look (slide-ready).
base_size  <- 18
slide_w    <- 13.33   # inches, 16:9
slide_h    <- 7.5
fig_dpi    <- 300
ink        <- "#1f1f1f"
ink_soft   <- "#5a5a5a"
grid_col   <- "#e6e6e3"
no_data    <- "#ececea"
heat_ramp  <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================

read_input <- function(file) {
  path <- file.path(work_dir, file)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun Steps 5 and 6 first.")
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

# Pairwise Simpson dissimilarity between two species sets.
beta_sim_pair <- function(A, B) {
  if (length(A) == 0 || length(B) == 0) return(c(a = NA, b = NA, c = NA, beta = NA))
  a <- length(intersect(A, B)); b <- length(setdiff(A, B)); c <- length(setdiff(B, A))
  c(a = a, b = b, c = c, beta = min(b, c) / (a + min(b, c)))
}

# Multisite Simpson dissimilarity (Rowan et al. 2024 Eqs. 1-3).
beta_sim_multi <- function(pools) {
  pools <- pools[lengths(pools) > 0]
  k <- length(pools)
  if (k < 2) return(c(n_stages = k, beta_SIM = NA, beta_SIM_END = NA, beta_SIM_SH = NA))
  S_i <- lengths(pools)
  S_T <- length(unique(unlist(pools)))
  min_sum <- 0; a_sum <- 0
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    min_sum <- min_sum + min(length(setdiff(pools[[i]], pools[[j]])),
                             length(setdiff(pools[[j]], pools[[i]])))
    a_sum <- a_sum + length(intersect(pools[[i]], pools[[j]]))
  }
  shared <- sum(S_i) - S_T
  c(n_stages = k,
    beta_SIM     = min_sum / (shared + min_sum),
    beta_SIM_END = min_sum / (min_sum + a_sum),
    beta_SIM_SH  = shared / sum(S_i))
}

# Shared slide theme.
theme_slide <- function() {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = ink),
      plot.title = element_text(face = "bold", size = base_size * 1.25, margin = margin(b = 4)),
      plot.subtitle = element_text(colour = ink_soft, size = base_size * 0.85, margin = margin(b = 12)),
      plot.caption = element_text(colour = ink_soft, size = base_size * 0.65, hjust = 0),
      axis.title = element_text(colour = ink_soft, size = base_size * 0.9),
      axis.text = element_text(colour = ink, size = base_size * 0.8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size * 0.95, hjust = 0),
      plot.margin = margin(18, 24, 14, 18),
      legend.position = "none"
    )
}

beta_lab <- quote(beta[SIM])

# =============================================================================
# 1. READ AND LINK SITES TO PROVINCES
# =============================================================================

cat("=== 1. READING FILES ===\n")
records <- read_input(records_file)
physio  <- read_input(physio_file)

# Province of each site (per stage and database). A FAUNMAP site name can
# cover several analysis units; if they fall in different provinces, the
# province of most units is used and the site is listed in the conflicts file.
site_province_all <- physio %>%
  transmute(Stage_Number = as.integer(Stage_Number), Database,
            SiteName = trimws(SiteName), Province = trimws(US_Province)) %>%
  filter(!is.na(SiteName), !is.na(Province))

site_province <- site_province_all %>%
  count(Stage_Number, Database, SiteName, Province) %>%
  group_by(Stage_Number, Database, SiteName) %>%
  mutate(n_provinces = n()) %>%
  arrange(desc(n), Province, .by_group = TRUE) %>%
  ungroup()

conflicts <- site_province %>% filter(n_provinces > 1)
site_province <- site_province %>%
  distinct(Stage_Number, Database, SiteName, .keep_all = TRUE) %>%
  select(Stage_Number, Database, SiteName, Province)
save_csv(conflicts, "site_province_conflicts")

# Match the requested provinces regardless of capitalisation.
available <- sort(unique(site_province$Province))
prov_match <- setNames(available[match(tolower(provinces), tolower(available))], provinces)
cat("\n  Provinces found in the site database:\n")
for (p in provinces) {
  cat(sprintf("    %-16s %s\n", p, if (is.na(prov_match[[p]])) "NOT FOUND" else "found"))
}
if (any(is.na(prov_match))) {
  cat("  Available US provinces:", paste(available, collapse = ", "), "\n")
}

recs <- records %>%
  transmute(Stage_Number = as.integer(Stage_Number), Database,
            SiteName = trimws(SiteName), GenusSpecies, Order) %>%
  inner_join(site_province, by = c("Stage_Number", "Database", "SiteName")) %>%
  mutate(Province = provinces[match(tolower(Province), tolower(provinces))]) %>%
  filter(!is.na(Province), Stage_Number %in% 1:5) %>%
  mutate(Stage = factor(stage_names[Stage_Number], levels = stage_names),
         Province = factor(Province, levels = provinces))

cat(sprintf("\n  Species records in the four provinces: %d (from %d sites)\n",
            nrow(recs), n_distinct(paste(recs$Database, recs$SiteName))))

# =============================================================================
# 2. REGIONAL SPECIES POOLS
# =============================================================================

cat("\n=== 2. REGIONAL SPECIES POOLS ===\n")

pools_long <- recs %>%
  group_by(Province, Stage, Stage_Number, GenusSpecies, Order) %>%
  summarise(n_sites_with_species = n_distinct(paste(Database, SiteName)), .groups = "drop") %>%
  arrange(Province, Stage_Number, Order, GenusSpecies)

# Sites per province and stage: sites contributing species to the pool, and
# all sites mapped in the province (some have no retained species).
sites_total <- site_province %>%
  mutate(Province = provinces[match(tolower(Province), tolower(provinces))]) %>%
  filter(!is.na(Province), Stage_Number %in% 1:5) %>%
  count(Province, Stage_Number, name = "n_sites_mapped")

pool_summary <- expand_grid(Province = factor(provinces, levels = provinces),
                            Stage_Number = 1:5) %>%
  left_join(recs %>%
              group_by(Province, Stage_Number) %>%
              summarise(n_sites = n_distinct(paste(Database, SiteName)),
                        n_FAUNMAP_sites = n_distinct(SiteName[Database == "FAUNMAP"]),
                        n_PBDB_sites = n_distinct(SiteName[Database == "PBDB"]),
                        n_species = n_distinct(GenusSpecies), .groups = "drop"),
            by = c("Province", "Stage_Number")) %>%
  left_join(mutate(sites_total, Province = factor(Province, levels = provinces)),
            by = c("Province", "Stage_Number")) %>%
  mutate(across(c(n_sites, n_FAUNMAP_sites, n_PBDB_sites, n_species, n_sites_mapped),
                ~ coalesce(.x, 0L)),
         Stage = factor(stage_names[Stage_Number], levels = stage_names)) %>%
  relocate(Stage, .after = Stage_Number)

print(as.data.frame(pool_summary %>% select(Province, Stage, n_sites, n_species)), row.names = FALSE)

save_csv(pools_long, "regional_pools_long")
save_csv(pool_summary, "pool_summary")

# One sheet per province: species x stage presence (1/0) and the number of sites.
pool_sheets <- lapply(provinces, function(p) {
  d <- filter(pools_long, Province == p)
  if (nrow(d) == 0) return(data.frame(Note = "No species recorded in this province"))
  d %>%
    mutate(Present = 1L) %>%
    select(Order, GenusSpecies, Stage, Present) %>%
    pivot_wider(names_from = Stage, values_from = Present, values_fill = 0L) %>%
    { for (s in stage_names) if (!s %in% names(.)) .[[s]] <- 0L; . } %>%
    select(Order, Species = GenusSpecies, all_of(stage_names)) %>%
    mutate(n_stages = rowSums(across(all_of(stage_names)))) %>%
    arrange(Order, Species)
})
names(pool_sheets) <- substr(provinces, 1, 31)
save_with_fallback(function(f) writexl::write_xlsx(pool_sheets, f),
                   file.path(output_dir, "regional_pools.xlsx"))

pool_sets <- lapply(provinces, function(p) {
  lapply(1:5, function(s) {
    unique(pools_long$GenusSpecies[pools_long$Province == p & pools_long$Stage_Number == s])
  })
})
names(pool_sets) <- provinces

# =============================================================================
# 3. SIMPSON DISSIMILARITY THROUGH TIME, WITHIN EACH PROVINCE
# =============================================================================

cat("\n=== 3. SIMPSON DISSIMILARITY (beta_SIM) ===\n")

beta_pairwise <- bind_rows(lapply(provinces, function(p) {
  bind_rows(lapply(1:4, function(i) bind_rows(lapply((i + 1):5, function(j) {
    r <- beta_sim_pair(pool_sets[[p]][[i]], pool_sets[[p]][[j]])
    data.frame(Province = p, Stage_A = stage_names[i], Stage_B = stage_names[j],
               shared_a = r[["a"]], only_A_b = r[["b"]], only_B_c = r[["c"]],
               beta_SIM = r[["beta"]])
  }))))
}))

beta_consecutive <- beta_pairwise %>%
  filter(match(Stage_B, stage_names) == match(Stage_A, stage_names) + 1) %>%
  mutate(Boundary_Ma = stage_young[match(Stage_A, stage_names)],
         Transition = paste(Stage_A, "to", Stage_B))

beta_multisite <- bind_rows(lapply(provinces, function(p) {
  r <- beta_sim_multi(pool_sets[[p]])
  data.frame(Province = p,
             Stages_with_species = paste(stage_names[lengths(pool_sets[[p]]) > 0], collapse = "; "),
             n_stages = r[["n_stages"]], beta_SIM = r[["beta_SIM"]],
             beta_SIM_END = r[["beta_SIM_END"]], beta_SIM_SH = r[["beta_SIM_SH"]],
             n_species_total = length(unique(unlist(pool_sets[[p]]))))
}))

save_csv(beta_pairwise, "beta_pairwise")
save_csv(beta_consecutive, "beta_consecutive")
save_csv(beta_multisite, "beta_multisite")

cat("  Multisite beta_SIM across stages (0 = same pool through time, 1 = complete turnover):\n")
print(as.data.frame(mutate(beta_multisite, across(where(is.numeric), ~ round(.x, 3))) %>%
        select(Province, n_stages, beta_SIM, beta_SIM_END, beta_SIM_SH)), row.names = FALSE)
cat("\n  Successive stages:\n")
print(as.data.frame(beta_consecutive %>% mutate(beta_SIM = round(beta_SIM, 3)) %>%
        select(Province, Transition, shared_a, beta_SIM)), row.names = FALSE)

# =============================================================================
# 4. FIGURES
# =============================================================================

cat("\n=== 4. FIGURES ===\n")

stage_axis <- function(n) {
  paste0(stage_names, "\n", formatC(stage_older, format = "g"), "-",
         formatC(stage_young, format = "g"), " Ma")
}

# ---- 4a. Sites per stage: one column graph per province ----------------------
site_plot <- function(d, title, colour, facet = FALSE) {
  ymax <- max(1, d$n_sites) * 1.15
  p <- ggplot(d, aes(x = Stage, y = n_sites)) +
    geom_col(fill = colour, width = 0.68) +
    geom_text(aes(label = n_sites), vjust = -0.45, size = base_size / 3.2,
              fontface = "bold", colour = ink) +
    scale_x_discrete(labels = setNames(stage_axis(), stage_names), drop = FALSE) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.02)),
                       breaks = scales::pretty_breaks(5)) +
    labs(x = NULL, y = "Number of sites") +
    theme_slide()
  if (facet) {
    p + facet_wrap(~ Province, ncol = 2, scales = "free_y") +
      labs(title = title,
           subtitle = "Sites contributing species to each regional pool (FAUNMAP + PBDB), oldest to youngest")
  } else {
    p + labs(title = title,
             subtitle = "Sites contributing species to the regional pool (FAUNMAP + PBDB), oldest to youngest")
  }
}

for (p in provinces) {
  d <- filter(pool_summary, Province == p)
  save_fig(site_plot(d, paste(p, "- sites per geological stage"), province_colours[[p]]),
           paste0("sites_per_stage_", file_stub(p)))
}

p_all_sites <- ggplot(pool_summary, aes(x = Stage, y = n_sites, fill = Province)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n_sites), vjust = -0.4, size = base_size / 3.8,
            fontface = "bold", colour = ink) +
  facet_wrap(~ Province, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = province_colours) +
  scale_x_discrete(labels = substr(stage_names, 1, 4), drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), breaks = scales::pretty_breaks(4)) +
  labs(title = "Sites per geological stage in four physiographic provinces",
       subtitle = "Zanclean (4.7 Ma) to Chibanian (0.129 Ma); sites contributing species to each regional pool",
       x = NULL, y = "Number of sites") +
  theme_slide()
save_fig(p_all_sites, "sites_per_stage_all_provinces")
cat("  Site-count column graphs saved\n")

# ---- 4b. Pairwise beta_SIM matrices ------------------------------------------
pair_grid <- function(p) {
  half <- beta_pairwise %>% filter(Province == p)
  # Lower triangle: row = later stage, column = earlier stage.
  expand_grid(Row = stage_names, Col = stage_names) %>%
    mutate(i = match(Col, stage_names), j = match(Row, stage_names)) %>%
    filter(j > i) %>%
    left_join(half, by = c("Col" = "Stage_A", "Row" = "Stage_B")) %>%
    mutate(Province = p)
}

pool_label <- function(p) {
  s <- filter(pool_summary, Province == p)
  setNames(sprintf("%s\n%d spp.\n%d sites", stage_names, s$n_species, s$n_sites), stage_names)
}

heat_plot <- function(g, title, subtitle, labels_x, labels_y, facet = FALSE, text_size = 6.5) {
  g <- g %>%
    mutate(Row = factor(Row, levels = rev(stage_names)),
           Col = factor(Col, levels = stage_names),
           lab = ifelse(is.na(beta_SIM), "no data", sprintf("%.2f", beta_SIM)),
           lab_col = ifelse(!is.na(beta_SIM) & beta_SIM > 0.5, "white", ink))
  p <- ggplot(g, aes(x = Col, y = Row)) +
    geom_tile(aes(fill = beta_SIM), colour = "white", linewidth = 1.5) +
    geom_text(aes(label = lab, colour = lab_col), size = text_size, fontface = "bold") +
    scale_colour_identity() +
    scale_fill_gradientn(colours = heat_ramp, limits = c(0, 1), na.value = no_data,
                         breaks = c(0, 0.25, 0.5, 0.75, 1), name = beta_lab) +
    scale_x_discrete(labels = labels_x, drop = TRUE, position = "bottom") +
    scale_y_discrete(labels = labels_y, drop = TRUE) +
    coord_equal() +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL,
         caption = paste0("0 = same species (or one pool nested in the other); 1 = no species shared.\n",
                          "Simpson dissimilarity after Rowan et al. (2024).")) +
    theme_slide() +
    theme(panel.grid = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title.position = "plot",
          plot.caption.position = "plot",
          axis.text = element_text(size = base_size * 0.75, lineheight = 0.9),
          legend.position = "right",
          legend.title = element_text(size = base_size, face = "bold"),
          legend.key.height = unit(1.6, "cm"),
          legend.key.width = unit(0.6, "cm"))
  if (facet) p <- p + facet_wrap(~ Province, ncol = 2)
  p
}

for (p in provinces) {
  m <- beta_multisite %>% filter(Province == p)
  sub <- if (is.na(m$beta_SIM)) "Fewer than two stages with species - no comparison possible" else
    sprintf("Pairwise Simpson dissimilarity between stage pools  |  multisite value (%d stages) = %.2f",
            m$n_stages, m$beta_SIM)
  lab <- pool_label(p)
  save_fig(heat_plot(pair_grid(p), paste(p, "- species pool turnover through time"), sub,
                     labels_x = lab, labels_y = lab),
           paste0("beta_pairwise_", file_stub(p)), width = slide_w, height = slide_h)
}

g_all <- bind_rows(lapply(provinces, pair_grid)) %>%
  mutate(Province = factor(Province, levels = provinces))
p_all_heat <- heat_plot(g_all, "Species pool turnover through time in four provinces",
                        "Pairwise Simpson dissimilarity between the stage pools of each province",
                        labels_x = setNames(substr(stage_names, 1, 4), stage_names),
                        labels_y = setNames(substr(stage_names, 1, 4), stage_names),
                        facet = TRUE, text_size = 4.6)
save_fig(p_all_heat, "beta_pairwise_all_provinces", width = slide_w, height = slide_h * 1.35)
cat("  Pairwise beta_SIM matrices saved\n")

# ---- 4c. Successive-stage turnover through time ------------------------------
bands <- data.frame(xmin = stage_young, xmax = stage_older, Stage = stage_names,
                    fill = rep(c("#f3f3f0", "#ffffff"), length.out = 5))
cons <- beta_consecutive %>% mutate(Province = factor(Province, levels = provinces))
label_pts <- cons %>% filter(!is.na(beta_SIM)) %>%
  group_by(Province) %>% slice_min(Boundary_Ma, n = 1) %>% ungroup()

# Spread end labels vertically so they never overlap (minimum gap in y units).
spread <- function(y, gap = 0.07) {
  o <- order(y); z <- y[o]
  for (k in seq_along(z)[-1]) z[k] <- max(z[k], z[k - 1] + gap)
  over <- max(0, z[length(z)] - 1.02); z <- z - over
  y[o] <- z; y
}
label_pts$label_y <- spread(label_pts$beta_SIM)
label_x <- 0.05

p_cons <- ggplot() +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = bands$fill) +
  geom_hline(yintercept = seq(0, 1, 0.25), colour = grid_col, linewidth = 0.4) +
  geom_text(data = bands, aes(x = (xmin + xmax) / 2, y = 1.07, label = Stage),
            size = base_size / 3.9, colour = ink_soft) +
  geom_line(data = filter(cons, !is.na(beta_SIM)),
            aes(Boundary_Ma, beta_SIM, colour = Province, group = Province), linewidth = 1.3) +
  geom_point(data = cons, aes(Boundary_Ma, beta_SIM, colour = Province),
             size = 4.2, shape = 21, fill = "white", stroke = 1.8, na.rm = TRUE) +
  geom_segment(data = label_pts, aes(x = Boundary_Ma - 0.06, xend = label_x + 0.03,
                                     y = beta_SIM, yend = label_y, colour = Province),
               linewidth = 0.5) +
  geom_text(data = label_pts, aes(x = label_x, y = label_y, label = Province, colour = Province),
            hjust = 0, size = base_size / 3.4, fontface = "bold") +
  scale_colour_manual(values = province_colours) +
  scale_x_reverse(limits = c(4.7, -1.45), breaks = c(4, 3, 2, 1, 0), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "Turnover of each regional pool between successive stages",
       subtitle = "Pairwise Simpson dissimilarity at each stage boundary (higher = more species replaced)",
       x = "Age (Ma)", y = beta_lab,
       caption = "Points sit on the boundary between the two stages compared. After Rowan et al. (2024).") +
  theme_slide() +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = grid_col))
save_fig(p_cons, "beta_consecutive_all_provinces")
cat("  Successive-stage turnover figure saved\n")

# ---- 4d. Multisite beta_SIM per province -------------------------------------
ms <- beta_multisite %>%
  mutate(Province = factor(Province, levels = provinces),
         lab = ifelse(is.na(beta_SIM), "n/a", sprintf("%.2f", beta_SIM)))
p_ms <- ggplot(ms, aes(x = Province, y = beta_SIM, fill = Province)) +
  geom_col(width = 0.62, na.rm = TRUE) +
  geom_text(aes(y = coalesce(beta_SIM, 0), label = lab), vjust = -0.5,
            size = base_size / 3, fontface = "bold", colour = ink) +
  geom_text(aes(y = 0, label = paste0(n_stages, " stages, ", n_species_total, " spp.")),
            vjust = 1.6, size = base_size / 4.2, colour = ink_soft) +
  scale_fill_manual(values = province_colours) +
  scale_y_continuous(limits = c(-0.08, 1.05), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "Overall turnover of each regional pool, Zanclean to Chibanian",
       subtitle = "Multisite Simpson dissimilarity across all stages with species (Rowan et al. 2024, Eq. 1)",
       x = NULL, y = beta_lab) +
  theme_slide()
save_fig(p_ms, "beta_multisite_all_provinces")
cat("  Multisite beta_SIM figure saved\n")

cat("\n=== STEP 7 COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: pools_long, pool_summary, beta_pairwise,",
    "beta_consecutive, beta_multisite\n")
