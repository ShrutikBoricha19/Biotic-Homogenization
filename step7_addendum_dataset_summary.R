# =============================================================================
# STEP 7 ADDENDUM: DATASET SUMMARY GRAPHIC
#
# One slide-sized image describing the large-mammal species in the analysis
# (every species in the province x species matrices of Step 7):
#   - total species, orders, families
#   - species per diet group and per body size class (with the Step 7 colours
#     and silhouettes)
#   - a few facts: heaviest and lightest species (exact body mass, Smith et
#     al. 2018), most widespread species, species found in every time bin,
#     most species-rich family, richest province x time bin
#
# Outputs (Outputs/7_rowan_figures/):
#   dataset_summary.png/.pdf    the graphic
#   dataset_summary_diet_size.png/.pdf   short version: only the diet and body size columns
#   dataset_summary_facts.csv   the numbers shown on it
#
# Inputs: Outputs/7_rowan_figures/regional_pa_matrices.rds and species_traits.csv
#         (Step 7); the silhouettes come from PhyloPic exactly as in Step 7 (same taxa, same
#         cache in Outputs/7_rowan_figures/phylopic/), or from step7_rowan_figures.R's built-in ones.
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 7.
# Needs: base R (grid); jsonlite and png for the PhyloPic silhouettes.
# =============================================================================

library(grid)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

pa_file     <- file.path("Outputs", "7_rowan_figures", "regional_pa_matrices.rds")
traits_file <- file.path("Outputs", "7_rowan_figures", "species_traits.csv")
step7_file  <- "step7_rowan_figures.R"      # source of the silhouettes (in work_dir)
output_dir  <- file.path(work_dir, "Outputs", "7_rowan_figures")

diet_groups <- c("Carnivore", "Omnivore", "Browser", "Mixed feeder", "Grazer")
size_groups <- c("Size 1 (<18 kg)", "Size 2 (18-80 kg)", "Size 3 (80-350 kg)",
                 "Size 4 (350-1,000 kg)", "Size 5 (>1,000 kg)")
diet_colours <- c("Carnivore" = "#D55E00", "Omnivore" = "#CC79A7", "Browser" = "#0072B2",
                  "Mixed feeder" = "#009E73", "Grazer" = "#E69F00")
size_colours <- setNames(c("#FDE725", "#5EC962", "#21918C", "#3B528B", "#440154"), size_groups)
ink       <- "#1F2A44"
ink_soft  <- "#3E4C6D"
rule_col  <- "#D6E6F5"
font      <- "sans"            # "sans" = Helvetica/Arial
fig_w <- 13.33; fig_h <- 7.5; fig_dpi <- 300
short_w <- 10; short_h <- 4.8              # size of the short (diet + body size) version
silhouette_source <- "phylopic"            # "phylopic" (taxa and cache as in Step 7) or "builtin"

# =============================================================================
# 1. DATA
# =============================================================================

setwd(work_dir)
if (!file.exists(pa_file) || !file.exists(traits_file)) stop("Run Step 7 first (", pa_file, ", ", traits_file, ").")
pa     <- readRDS(pa_file)
traits <- read.csv(traits_file, stringsAsFactors = FALSE, check.names = FALSE)

# Silhouettes: the path strings defined in Step 7.
s7 <- readLines(step7_file, warn = FALSE)
i0 <- grep("^diet_silhouette_paths <- c\\(", s7); i1 <- grep("^diet_silhouette <- function", s7)
if (!length(i0) || !length(i1)) stop("Silhouettes not found in ", step7_file)
eval(parse(text = s7[i0:(i1 - 1)]))
# PhyloPic: the taxa and download helpers defined in Step 7.
if (silhouette_source == "phylopic") {
  j0 <- grep("^phylopic_diet <- c\\(", s7); j1 <- grep("^phylopic_refresh <-", s7)
  k0 <- grep("^# PhyloPic silhouettes, straight from", s7); k1 <- grep("^# Silhouette grob from PhyloPic", s7)
  ok <- length(j0) && length(j1) && length(k0) && length(k1) &&
        all(vapply(c("jsonlite", "png"), requireNamespace, logical(1), quietly = TRUE))
  if (ok) {
    eval(parse(text = s7[j0:j1])); eval(parse(text = s7[k0:(k1 - 1)]))
  } else {
    message("  PhyloPic code or the jsonlite/png packages not found; using the built-in silhouettes.")
    silhouette_source <- "builtin"
  }
}

species <- sort(unique(unlist(lapply(pa, colnames))))
sp <- traits[match(species, traits$Species), ]
sp$Species <- species
n_species  <- length(species)
n_orders   <- length(unique(na.omit(sp$Order)))
n_families <- length(unique(na.omit(sp$Family)))

diet_n <- table(factor(sp$Diet_Group, levels = diet_groups))
size_n <- table(factor(sp$Size_Class, levels = size_groups))
n_no_diet <- sum(is.na(sp$Diet_Group) | !sp$Diet_Group %in% diet_groups)
n_no_mass <- sum(is.na(sp$Mass_kg))

# Facts
heavy <- sp[which.max(sp$Mass_kg), ]
light <- sp[which.min(sp$Mass_kg), ]
cells <- sapply(species, function(s) sum(vapply(pa, function(m) if (s %in% colnames(m)) sum(m[, s]) else 0, numeric(1))))
n_cells <- sum(vapply(pa, nrow, numeric(1)))
wide  <- names(cells)[cells == max(cells)]
in_all_bins <- species[vapply(species, function(s) all(vapply(pa, function(m) s %in% colnames(m) && any(m[, s] > 0), logical(1))), logical(1))]
fam_tab <- sort(table(sp$Family), decreasing = TRUE)
rich <- do.call(rbind, lapply(names(pa), function(b) data.frame(Bin = b, Province = rownames(pa[[b]]), S = rowSums(pa[[b]]))))
rich_top <- rich[which.max(rich$S), ]

fmt_kg <- function(x) paste(formatC(x, format = "f", digits = if (x < 100) 2 else 1, big.mark = ","), "kg")
lead_name <- function(x) if (length(x) == 1) x else sprintf("%s and %d others", x[1], length(x) - 1)

facts <- data.frame(
  Fact  = c("Species", "Orders", "Families", "Heaviest", "Lightest", "Most widespread",
            "In every time bin", "Most species-rich family", "Richest province x bin",
            "Without diet data", "Without body mass"),
  Value = c(n_species, n_orders, n_families,
            sprintf("%s (%s)", heavy$Species, fmt_kg(heavy$Mass_kg)),
            sprintf("%s (%s)", light$Species, fmt_kg(light$Mass_kg)),
            sprintf("%s (%d of %d province x bin cells)", paste(wide, collapse = "; "), max(cells), n_cells),
            sprintf("%d species: %s", length(in_all_bins), paste(in_all_bins, collapse = "; ")),
            sprintf("%s (%d species)", names(fam_tab)[1], fam_tab[[1]]),
            sprintf("%s, Bin %s (%d species)", rich_top$Province, rich_top$Bin, rich_top$S),
            n_no_diet, n_no_mass))
write.csv(facts, file.path(output_dir, "dataset_summary_facts.csv"), row.names = FALSE)
print(facts, right = FALSE)

# =============================================================================
# 2. GRAPHIC
# =============================================================================

silhouette <- function(path_string, x, y, h, col) {
  parts <- strsplit(path_string, ";", fixed = TRUE)[[1]]
  xy <- do.call(rbind, lapply(seq_along(parts), function(i) {
    v <- as.numeric(strsplit(trimws(parts[i]), " +")[[1]])
    data.frame(x = v[c(TRUE, FALSE)], y = v[c(FALSE, TRUE)], id = i)
  }))
  w <- diff(range(xy$x)); hh <- diff(range(xy$y))
  pushViewport(viewport(x, y, width = unit(h * w / hh, "in"), height = unit(h, "in")))
  grid.path((xy$x - min(xy$x)) / w, 1 - (xy$y - min(xy$y)) / hh, id = xy$id, rule = "evenodd",
            gp = gpar(fill = col, col = if (col == size_colours[[1]]) ink_soft else NA, lwd = 0.6))
  popViewport()
}
# PhyloPic silhouette centred at (x, y), h inches tall; FALSE if it could not be drawn.
phylopic_memo <- new.env()                    # one download attempt per taxon
phylopic_icon <- function(taxon, x, y, h, col) {
  if (is.na(taxon)) return(FALSE)
  if (!exists(taxon, envir = phylopic_memo, inherits = FALSE))
    assign(taxon, tryCatch(png::readPNG(phylopic_fetch(taxon)), error = function(e) {
      message("  PhyloPic: no silhouette for ", taxon, " (", conditionMessage(e), ")"); NULL }), envir = phylopic_memo)
  img <- get(taxon, envir = phylopic_memo)
  if (is.null(img)) return(FALSE)
  if (length(dim(img)) == 2) img <- array(rep(img, 2), c(dim(img), 2))
  alpha <- switch(as.character(dim(img)[3]), "2" = img[, , 2], "4" = img[, , 4], 1 - img[, , 1])
  tint <- function(cl) { out <- array(0, c(dim(img)[1:2], 4)); rgb <- col2rgb(cl) / 255
                         for (k in 1:3) out[, , k] <- rgb[k]; out[, , 4] <- alpha; out }
  w <- h * dim(img)[2] / dim(img)[1]
  if (w > 1.6 * h) { w <- 1.6 * h; h <- w * dim(img)[1] / dim(img)[2] }     # keep wide animals in their slot
  if (col == size_colours[[1]])                                              # thin outline for the pale yellow
    for (d in list(c(-0.7, 0), c(0.7, 0), c(0, -0.7), c(0, 0.7)))
      grid.raster(tint(ink_soft), x + unit(d[1], "pt"), y + unit(d[2], "pt"), unit(w, "in"), unit(h, "in"), interpolate = TRUE)
  grid.raster(tint(col), x, y, unit(w, "in"), unit(h, "in"), interpolate = TRUE)
  TRUE
}
icon <- function(taxon, path_string, x, y, h, col) {
  if (silhouette_source == "phylopic" && phylopic_icon(taxon, x, y, h, col)) return(invisible())
  silhouette(path_string, x, y, h, col)
}
txt <- function(label, x, y, size, col = ink, face = "plain", just = "left")
  grid.text(label, x, y, just = just, gp = gpar(fontfamily = font, fontsize = size, col = col, fontface = face))

# One column of labelled bars. x0, top: npc; col_w and bar_w: npc; step: npc between rows.
bar_block <- function(x0, title, groups, counts, colours, taxa, paths, labels,
                      top = 0.75, step = 0.105, col_w = 0.27, bar_w = 0.19, title_size = 17) {
  txt(title, x0, top, title_size, face = "bold")
  grid.lines(c(x0, x0 + col_w), c(top - step * 0.25, top - step * 0.25), gp = gpar(col = rule_col, lwd = 1.5))
  max_n <- max(counts, 1)
  bx <- unit(x0, "npc") + unit(0.6, "in")                     # bars start right of the silhouettes
  for (i in seq_along(groups)) {
    y <- top - step * 0.9 - (i - 1) * step
    icon(taxa[i], paths[[i]], unit(x0, "npc") + unit(0.25, "in"), unit(y, "npc"), 0.42, colours[[i]])
    grid.text(labels[i], bx, unit(y + step * 0.21, "npc"), just = "left",
              gp = gpar(fontfamily = font, fontsize = 12.5, col = ink_soft))
    len <- bar_w * counts[[i]] / max_n
    grid.rect(bx, unit(y - step * 0.17, "npc"), unit(len, "npc"), unit(step * 0.21, "npc"), just = c("left", "centre"),
              gp = gpar(fill = colours[[i]], col = if (colours[[i]] == size_colours[[1]]) ink_soft else NA, lwd = 0.6))
    grid.text(counts[[i]], bx + unit(len, "npc") + unit(0.07, "in"), unit(y - step * 0.17, "npc"), just = "left",
              gp = gpar(fontfamily = font, fontsize = 13, col = ink, fontface = "bold"))
  }
}

fact <- function(y, label, name, value, italic = TRUE, value_italic = FALSE) {
  txt(toupper(label), 0.70, y + 0.035, 10.5, ink_soft, "bold")
  txt(name, 0.70, y, 15, ink, if (italic) "italic" else "plain")
  txt(value, 0.70, y - 0.036, 12.5, ink_soft, if (value_italic) "italic" else "plain")
}

diet_taxa <- if (exists("phylopic_diet")) unname(phylopic_diet[diet_groups]) else rep(NA, 5)
size_taxa <- if (exists("phylopic_size")) unname(phylopic_size) else rep(NA, 5)
size_labels <- sub("^(Size \\d) \\((.*)\\)$", "\\1  \u00b7  \\2", size_groups)
diet_columns <- function(x_diet, x_size, ...) {
  bar_block(x_diet, "Diet", diet_groups, diet_n, diet_colours, diet_taxa,
            lapply(diet_groups, function(g) diet_silhouette_paths[[g]]), diet_groups, ...)
  bar_block(x_size, "Body size", size_groups, size_n, size_colours, size_taxa,
            as.list(size_silhouette_paths), size_labels, ...)
}

draw <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  # headline
  txt(n_species, 0.04, 0.89, 46, face = "bold")
  hx <- 0.04 + convertWidth(grobWidth(textGrob(n_species, gp = gpar(fontfamily = font, fontsize = 46, fontface = "bold"))),
                            "npc", valueOnly = TRUE) + 0.012
  txt("large-mammal species", hx, 0.905, 18)
  txt(sprintf("%d orders  \u00b7  %d families  \u00b7  3 provinces  \u00b7  5 time bins, 3.25 Ma to 11.7 ka",
              n_orders, n_families), hx, 0.868, 13, ink_soft)

  diet_columns(0.04, 0.36)

  txt("Facts", 0.70, 0.75, 17, face = "bold")
  grid.lines(c(0.70, 0.96), c(0.725, 0.725), gp = gpar(col = rule_col, lwd = 1.5))
  fact(0.645, "Heaviest", heavy$Species, sprintf("%s  \u00b7  %s", fmt_kg(heavy$Mass_kg), heavy$Family))
  fact(0.525, "Lightest", light$Species, sprintf("%s  \u00b7  %s", fmt_kg(light$Mass_kg), light$Family))
  fact(0.405, "Most widespread", lead_name(wide),
       sprintf("in %d of %d province \u00d7 time-bin cells", max(cells), n_cells))
  fact(0.285, "Found in every time bin", sprintf("%d species", length(in_all_bins)),
       if (length(in_all_bins)) paste(head(in_all_bins, 2), collapse = ", ") else "none",
       italic = FALSE, value_italic = length(in_all_bins) > 0)
  fact(0.165, "Most species-rich family", names(fam_tab)[1], sprintf("%d species", fam_tab[[1]]))

  foot <- sprintf(paste0("Species in the province \u00d7 species matrices (Step 7). Body mass and diet from Smith et al. (2018)%s.",
                         if (n_no_diet) " %d species without diet data are not in the diet counts." else "%s"),
                  if (n_no_mass) sprintf("; %d without body mass", n_no_mass) else "", if (n_no_diet) n_no_diet else "")
  txt(foot, 0.04, 0.045, 10.5, ink_soft)
}

draw_short <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  diet_columns(0.05, 0.53, top = 0.92, step = 0.172, col_w = 0.42, bar_w = 0.27, title_size = 18)
}

save_both <- function(fun, name, w, h) {
  png(file.path(output_dir, paste0(name, ".png")), width = w, height = h, units = "in", res = fig_dpi,
      type = if (capabilities("cairo")) "cairo" else "windows")
  fun(); invisible(dev.off())
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
  pdf_dev(file.path(output_dir, paste0(name, ".pdf")), width = w, height = h)
  fun(); invisible(dev.off())
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_both(draw, "dataset_summary", fig_w, fig_h)
save_both(draw_short, "dataset_summary_diet_size", short_w, short_h)
if (exists("phylopic_credits") && length(phylopic_credits))
  write.csv(do.call(rbind, phylopic_credits), file.path(output_dir, "phylopic_credits.csv"), row.names = FALSE)
cat("\nSaved dataset_summary, dataset_summary_diet_size (.png/.pdf) and dataset_summary_facts.csv to", output_dir, "\n")
