# =============================================================================
# STEP 8b (ADDENDUM TO STEP 8): ENDEMIC VS SHARED SPECIES BY DIET
#         (Rowan et al. 2024, Nat. Ecol. Evol., Fig. 4a, Eqs. 5-6)
#
# Rowan et al. split multisite beta_SIM into the species found at only one
# site (endemic) and the species found at two or more sites (shared), and
# partitioned both by diet:
#
#   Endemic:  beta_SIM_END_f = sum_{i<j} min(b_ij, b_ji)_f /
#                              [ sum_{i<j} min(b_ij, b_ji) + a_ij ]        (Eq. 5)
#   Shared:   beta_SIM_SH_f  = ( sum_i S_if - S_Tf ) / sum_i S_i            (Eq. 6)
#
# The values are calculated in Step 8 (same diets, same site subsets) and
# saved in diet_partition.csv; this script draws them. One figure per
# province: shared species (top row) and endemic species (bottom row) for
# each diet group, stage by stage. Diet groups with no species in a province
# are left out. Falling endemic bars + rising shared bars = homogenization.
#
# Input:   Outputs/8_diet_partitioning/diet_partition.csv   (Step 8)
# Outputs (Outputs/8_diet_partitioning/):
#   endemic_shared_<province>.png/.pdf   Fig. 4a-style figure, one per province
#   endemic_shared_summary.csv           the plotted values per diet and stage
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 8.
# Needs: dplyr, tidyr, ggplot2 (gtable and grid come with ggplot2).
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "gtable")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Run: install.packages(\"", pkg, "\")")
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

# =============================================================================
# SETTINGS
# =============================================================================

work_dir <- "C:/Users/shrut/OneDrive/Documents/Data D/Ph.D/Research/Dissertation_Chapter_1"

output_dir     <- file.path(work_dir, "Outputs", "8_diet_partitioning")
partition_file <- file.path(output_dir, "diet_partition.csv")

provinces <- c("Pacific Border", "Basin and Range", "Great Plains", "Coastal Plain")

# Must match Step 8 (only used in the figure caption).
n_units_sample <- 3

diet_colours <- c(
  "Carnivore"             = "#eb6834",   # orange
  "Omnivore"              = "#2a78d6",   # blue
  "Frugivore"             = "#e87ba4",   # magenta
  "Browser"               = "#4a3aa7",   # violet
  "Browser and frugivore" = "#4a3aa7",
  "Mixed feeder"          = "#1baf7a",   # aqua
  "Grazer"                = "#eda100"    # yellow
)

stage_names <- c("Zanclean", "Piacenzian", "Gelasian", "Calabrian", "Chibanian")
stage_older <- c(4.700, 3.600, 2.580, 1.800, 0.7741)
stage_young <- c(3.600, 2.580, 1.800, 0.7741, 0.129)

# Figure look (slide-ready).
base_size <- 18
slide_w   <- 13.33   # inches, 16:9
slide_h   <- 7.5
fig_dpi   <- 300
ink       <- "#1f1f1f"
ink_soft  <- "#5a5a5a"
grid_col  <- "#e6e6e3"
overall_line <- "#6b6b68"
overall_band <- "#d6d6d2"
strip_fill   <- "#ececea"
icon_height_in <- 0.46   # silhouette height in the panel titles (inches)

# =============================================================================
# HELPERS
# =============================================================================

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

# =============================================================================
# SILHOUETTES
# Game-icons.net (by Lorc, Delapouite and contributors), CC BY 3.0:
# feline, polar bear, squirrel, deer, camel, bison. Stored as polygon outlines
# (x y pairs on a 512 x 512 grid; ";" separates parts) so no download or
# extra package is needed.
# =============================================================================

silhouette_paths <- c(
  "Carnivore" = paste0(
    "110 65 98 66 83 72 72 80 64 89 57 104 52 121 50 138 51 150 55 163 62 178 73 195 84 210 ",
    "81 234 81 262 83 282 87 293 103 318 167 403 196 446 201 447 206 446 211 444 196 381 ",
    "184 342 177 327 163 300 201 311 220 313 243 314 267 313 288 311 306 307 322 301 357 ",
    "341 411 413 432 439 439 447 443 446 446 444 448 441 449 438 418 345 402 305 385 272 ",
    "390 256 398 240 406 229 416 218 433 204 444 210 451 212 457 211 464 209 473 202 481 ",
    "192 488 182 492 171 478 163 470 139 465 131 460 124 446 114 425 104 426 92 425 79 412 ",
    "91 395 109 355 159 335 183 295 179 243 177 191 179 131 184 123 183 112 178 101 169 92 ",
    "159 87 148 85 139 87 123 92 113 99 107 110 103 125 104 138 107 181 125 200 131 223 136 ",
    "240 136 247 134 257 128 278 110 289 97 297 85 304 73 307 65 299 67 285 72 240 100 233 ",
    "102 223 103 211 102 196 97 152 76 134 69 121 66 110 65;453 138 463 160 441 156;75 309 ",
    "57 363 46 389 31 420 20 442 24 446 29 447 47 432 69 412 89 390 113 361 90 331 75 ",
    "309;317 322 292 328 281 358 268 387 239 442 242 445 246 447 249 445 272 423 289 405 ",
    "305 386 334 342 317 322"),
  "Omnivore" = paste0(
    "238 81 212 83 191 86 164 94 144 102 130 109 103 127 80 145 62 164 48 184 37 203 30 223 ",
    "26 242 24 261 24 279 25 296 28 313 36 341 46 370 54 379 63 386 77 395 121 416 120 396 ",
    "138 403 137 388 158 388 125 351 130 342 134 331 140 295 200 308 256 301 259 330 264 ",
    "353 270 375 276 393 346 431 342 412 366 416 363 399 387 400 341 358 342 317 348 311 ",
    "355 307 363 305 374 302 385 302 394 304 393 313 395 323 408 354 419 369 441 390 443 ",
    "371 465 375 464 359 489 362 464 324 464 289 462 262 454 272 442 281 429 287 414 289 ",
    "394 288 376 285 363 278 351 266 366 256 371 262 379 267 389 270 398 271 414 271 422 ",
    "270 431 266 440 260 448 250 448 241 460 229 466 218 461 205 460 197 461 190 470 175 ",
    "470 168 469 162 463 159 453 156 436 171 417 168 402 169 387 171 368 177 341 165 336 ",
    "163 332 164 329 166 326 171 322 181 323 188 325 195 328 199 333 202 326 219 318 214 ",
    "311 207 307 199 304 189 304 179 307 169 311 160 318 151 328 146 339 145 350 149 369 ",
    "158 388 153 406 151 391 136 374 123 355 111 334 101 313 93 290 86 266 82 238 81;427 ",
    "183 432 184 436 187 438 192 439 198 438 203 436 208 432 211 427 212 422 211 418 208 ",
    "416 203 415 198 416 192 418 187 422 184 427 183;378 188 383 189 387 192 389 197 390 ",
    "203 389 208 387 213 383 216 378 218 373 216 370 213 367 208 366 203 367 197 370 192 ",
    "373 189 378 188;149 315 157 353 195 378 193 363 210 366 205 354 223 355 203 336 201 ",
    "329"),
  "Frugivore" = paste0(
    "206 24 183 25 159 30 135 38 111 49 89 61 68 76 49 92 32 109 55 108 82 110 101 114 119 ",
    "120 126 125 133 130 138 137 142 144 145 157 145 168 142 177 136 186 120 201 77 233 67 ",
    "243 57 254 49 267 43 282 38 299 36 318 37 334 40 349 45 365 51 380 59 394 69 408 80 ",
    "422 91 434 104 446 118 456 133 465 148 473 164 479 180 484 196 487 215 488 395 488 401 ",
    "487 405 484 406 480 404 475 392 466 369 457 346 451 369 425 384 402 392 382 394 366 ",
    "390 352 382 342 371 335 357 330 342 329 325 331 309 336 294 344 280 355 270 369 263 ",
    "386 260 405 257 401 254 389 255 374 261 357 266 349 280 333 289 326 299 320 311 315 ",
    "325 312 341 310 358 310 375 321 386 326 396 330 406 332 428 332 447 329 449 327 452 ",
    "322 453 315 451 308 449 305 441 302 409 299 390 279 385 269 385 262 390 253 403 242 ",
    "435 249 451 252 467 251 471 250 476 246 479 237 479 217 476 196 468 173 459 159 448 ",
    "147 433 139 424 136 414 135 412 127 415 100 414 93 410 92 402 93 394 96 386 101 375 ",
    "111 369 118 364 129 364 135 353 141 338 154 324 170 304 195 292 201 279 209 268 218 ",
    "257 229 236 253 219 279 204 306 190 342 184 366 186 334 191 305 198 281 208 260 219 ",
    "241 232 225 245 211 259 199 298 169 308 161 316 152 322 142 325 132 324 121 321 108 ",
    "313 86 301 68 287 53 271 41 254 32 235 27 206 24;433 170 438 171 442 174 445 178 446 ",
    "183 445 188 442 192 438 195 433 196 428 195 424 192 421 188 420 183 421 178 424 174 ",
    "428 171 433 170"),
  "Browser" = paste0(
    "139 20 123 25 129 42 135 56 145 71 155 81 166 88 178 92 194 94 211 94 251 89 272 90 ",
    "293 95 318 108 326 106 330 104 331 97 331 88 326 76 320 67 313 60 294 48 284 44 279 60 ",
    "290 65 299 71 307 78 312 85 290 75 276 71 258 67 223 63 212 59 205 55 200 48 195 39 ",
    "192 29 176 33 179 44 184 54 191 65 200 72 211 77 196 77 183 76 169 71 162 66 156 59 ",
    "147 43 139 20;229 24 229 39 233 52 255 54 250 47 247 40 245 32 245 25;349 37 338 50 ",
    "345 58 348 69 349 73 336 64 341 76 344 95 360 97 365 81 365 72 364 62 361 52 358 47 ",
    "349 37;268 107 250 108 244 109 240 113 236 120 248 132 256 137 267 142 288 147 305 148 ",
    "303 158 303 178 302 190 298 202 295 207 287 217 277 226 270 229 256 235 238 239 199 ",
    "241 175 244 150 250 127 258 111 268 96 281 86 294 80 309 87 320 102 320 102 381 84 410 ",
    "89 492 107 492 108 426 128 412 147 395 162 377 172 359 179 339 182 316 179 291 173 273 ",
    "188 267 196 289 198 303 199 316 198 329 195 341 192 352 185 369 224 367 258 363 253 ",
    "306 270 305 286 492 303 492 306 370 315 365 327 355 335 345 343 332 352 315 359 294 ",
    "367 258 371 223 372 190 348 182 423 182 431 165 424 165 415 162 417 154 432 154 347 ",
    "110 331 117 320 124 309 117 295 111 283 108 268 107;364 333 351 354 335 374 382 395 ",
    "361 443 376 450 409 385 365 348 364 333"),
  "Mixed feeder" = paste0(
    "421 27 403 30 387 37 379 31 371 30 368 31 363 37 362 42 364 48 368 53 375 59 382 105 ",
    "383 128 382 142 380 148 374 158 365 165 356 167 349 165 341 160 335 153 322 132 297 82 ",
    "282 60 273 51 264 44 253 39 241 38 234 39 224 44 217 52 205 79 199 89 191 96 185 99 ",
    "178 100 169 98 161 92 156 83 141 53 135 47 131 45 122 44 108 48 100 53 93 62 86 74 80 ",
    "89 74 108 69 131 60 132 49 137 43 143 35 154 28 170 21 190 17 216 16 236 16 269 23 353 ",
    "30 276 35 295 42 313 47 341 48 369 43 378 42 387 43 397 48 405 46 455 46 485 59 485 ",
    "146 485 145 479 139 471 128 464 111 455 111 407 116 398 118 388 116 378 111 370 112 ",
    "345 113 331 119 317 144 279 165 282 188 282 211 279 235 275 236 360 232 369 230 379 ",
    "232 389 236 398 238 414 239 436 237 485 295 485 294 477 288 471 280 465 261 455 259 ",
    "410 260 395 264 381 263 371 261 364 264 314 267 290 271 265 281 259 292 358 289 367 ",
    "288 377 291 387 296 395 298 406 300 428 303 484 365 485 362 479 355 473 324 459 320 ",
    "439 318 395 322 385 323 374 321 364 316 355 318 326 324 283 328 262 335 240 354 239 ",
    "370 236 383 231 395 226 405 219 412 212 419 204 427 185 430 175 432 153 434 107 435 ",
    "103 437 100 440 97 448 95 458 96 483 99 490 99 494 96 496 88 496 81 494 74 491 67 484 ",
    "61 480 58 472 55 457 53 433 32 427 28 422 27;428 46 431 51 430 57 420 57 409 52 418 48 ",
    "428 46;81 302 81 304 89 328 93 366 94 372 89 380 88 390 90 399 94 407 93 471 89 466 83 ",
    "461 76 457 68 454 67 404 71 395 72 386 70 377 66 369 66 331 69 316 81 302"),
  "Grazer" = paste0(
    "300 99 286 100 277 102 257 110 239 120 193 148 166 162 149 166 110 172 86 177 74 182 ",
    "64 189 56 200 54 207 48 235 37 263 24 274 21 277 20 282 20 299 23 307 26 312 41 302 41 ",
    "290 43 277 55 249 56 265 54 287 43 337 46 357 53 386 64 413 100 413 86 395 78 381 74 ",
    "367 73 352 77 337 94 334 102 330 106 327 127 382 136 397 146 413 187 413 158 383 149 ",
    "368 144 356 143 349 142 332 143 322 177 324 192 324 209 322 232 316 229 352 230 356 ",
    "236 369 243 380 267 413 306 413 277 386 271 378 267 369 265 357 265 344 274 320 294 ",
    "350 312 374 327 391 345 408 352 413 388 413 353 380 339 363 328 343 325 328 325 320 ",
    "354 327 360 327 366 326 375 347 383 362 394 375 403 382 416 366 424 352 428 336 429 ",
    "316 454 331 462 329 466 325 470 321 470 317 467 286 478 266 492 252 492 237 489 223 ",
    "480 199 467 176 464 173 453 170 417 164 392 148 349 117 335 109 314 100 300 99;373 134 ",
    "388 163 397 174 409 184 430 196 438 204 441 209 442 217 441 225 438 231 432 238 427 ",
    "242 421 245 416 244 399 238 384 231 373 223 364 214 357 202 354 195 353 185 356 168 ",
    "365 148;373 172 371 181 372 190 374 197 378 203 384 209 393 215 405 221 419 227 423 ",
    "219 424 215 423 213 393 193 382 183 373 172")
)

# Grob that draws one silhouette, right-aligned in its cell, aspect preserved.
silhouette_grob <- function(diet, colour, height_in = icon_height_in) {
  key <- if (diet == "Browser and frugivore") "Browser" else diet
  if (!key %in% names(silhouette_paths)) return(nullGrob())
  parts <- strsplit(silhouette_paths[[key]], ";", fixed = TRUE)[[1]]
  xy <- do.call(rbind, lapply(seq_along(parts), function(i) {
    v <- as.numeric(strsplit(trimws(parts[i]), " +")[[1]])
    data.frame(x = v[c(TRUE, FALSE)], y = v[c(FALSE, TRUE)], id = i)
  }))
  w <- diff(range(xy$x)); h <- diff(range(xy$y))
  vp <- viewport(x = unit(1, "npc") - unit(10, "pt"), y = unit(0.5, "npc"),
                 width = unit(height_in * w / h, "in"), height = unit(height_in, "in"),
                 just = c("right", "centre"))
  pathGrob((xy$x - min(xy$x)) / w, 1 - (xy$y - min(xy$y)) / h, id = xy$id,
           rule = "evenodd", gp = gpar(fill = colour, col = NA), vp = vp)
}

# =============================================================================
# 1. READ THE STEP 8 RESULTS
# =============================================================================

cat("=== 1. READING STEP 8 RESULTS ===\n")
if (!file.exists(partition_file)) stop("File not found:\n  ", partition_file, "\nRun Step 8 first.")
partition <- read.csv(partition_file, check.names = FALSE, stringsAsFactors = FALSE)
names(partition) <- trimws(sub("^[^A-Za-z0-9]+", "", names(partition), useBytes = TRUE))
if (!all(c("beta_END_f", "beta_SH_f") %in% names(partition))) {
  stop("diet_partition.csv has no endemic/shared columns - re-run the current Step 8 first.")
}
cat(sprintf("  %d rows read from %s\n", nrow(partition), basename(partition_file)))

n_col <- grep("^n_(sites|sections)$", names(partition), value = TRUE)[1]
unit_word <- sub("^n_", "", n_col)
categories <- unique(partition$Diet)
diet_levels <- setdiff(categories, "Unclassified")

partition <- partition %>%
  mutate(Province = factor(Province, levels = provinces),
         Stage = factor(Stage, levels = stage_names),
         Diet = factor(Diet, levels = categories))

# A diet group is drawn for a province only if it has species there.
diet_present <- partition %>%
  group_by(Province, Diet) %>%
  summarise(any_species = sum(n_species, na.rm = TRUE) > 0, .groups = "drop")

bands <- data.frame(xmin = stage_young, xmax = stage_older, Stage = stage_names,
                    lab = substr(stage_names, 1, 4),
                    fill = rep(c("#f3f3f0", "#ffffff"), length.out = 5))

theme_panels <- function() {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = ink),
      plot.title = element_text(face = "bold", size = base_size * 1.25, margin = margin(b = 4)),
      plot.subtitle = element_text(colour = ink_soft, size = base_size * 0.8, margin = margin(b = 12)),
      plot.caption = element_text(colour = ink_soft, size = base_size * 0.6, hjust = 0,
                                  lineheight = 1.1),
      plot.title.position = "plot", plot.caption.position = "plot",
      axis.title = element_text(colour = ink_soft, size = base_size * 0.85),
      axis.text = element_text(colour = ink, size = base_size * 0.7),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "#c9c9c4", fill = NA, linewidth = 0.5),
      panel.spacing.x = unit(1.1, "lines"), panel.spacing.y = unit(1.3, "lines"),
      strip.background = element_rect(fill = strip_fill, colour = "#c9c9c4", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = base_size * 0.95, hjust = 0,
                                margin = margin(t = 11, b = 11, l = 8)),
      plot.margin = margin(18, 24, 14, 18),
      legend.position = "none"
    )
}

# =============================================================================
# 2. ENDEMIC VS SHARED SPECIES BY DIET (Rowan et al. 2024 Fig. 4a, Eqs. 5-6)
#
#   Endemic (bottom row):  beta_SIM_END_f = sum_{i<j} min(b_ij, b_ji)_f /
#                                           [ sum_{i<j} min(b_ij, b_ji) + a_ij ]
#   Shared (top row):      beta_SIM_SH_f  = ( sum_i S_if - S_Tf ) / sum_i S_i
#
#   Endemic = species found at only one of the compared sites (more precisely,
#   the endemics of whichever site has fewer); shared = species found at two
#   or more sites. Within a row, the diet groups (plus unclassified species,
#   not drawn) add up to beta_SIM_END or beta_SIM_SH. Computed in Step 8 on the same site subsets.
#   Falling endemic bars + rising shared bars through time = homogenization.
# =============================================================================

cat("\n=== 2. ENDEMIC VS SHARED SPECIES ===\n")

es_table <- partition %>%
  filter(Diet != "Unclassified") %>%
  select(Province, Diet, Stage, beta_SH_f, beta_END_f) %>%
  pivot_longer(c(beta_SH_f, beta_END_f), names_to = "Component", values_to = "value") %>%
  mutate(Component = ifelse(Component == "beta_SH_f", "Shared", "Endemic"),
         value = round(value, 3)) %>%
  pivot_wider(names_from = Stage, values_from = value) %>%
  left_join(diet_present, by = c("Province", "Diet")) %>%
  filter(any_species) %>% select(-any_species) %>%
  arrange(Province, Component, Diet)
save_csv(es_table, "endemic_shared_summary")
print(as.data.frame(es_table), row.names = FALSE)

es_plot <- function(p) {
  keep <- diet_present %>% filter(Province == p, any_species, Diet != "Unclassified") %>%
    pull(Diet) %>% as.character()
  if (length(keep) == 0) return(NULL)
  comp_levels <- c("Shared species", "Endemic species")
  d <- partition %>%
    filter(Province == p, Diet %in% keep) %>%
    select(Diet, Stage_Number, beta_SIM, beta_SH_f, beta_END_f) %>%
    pivot_longer(c(beta_SH_f, beta_END_f), names_to = "Component", values_to = "value") %>%
    mutate(Component = factor(ifelse(Component == "beta_SH_f", comp_levels[1], comp_levels[2]),
                              levels = comp_levels),
           Diet = factor(as.character(Diet), levels = keep),
           bar_fill = diet_colours[as.character(Diet)])
  ymax <- max(0.1, d$value, na.rm = TRUE) * 1.32
  grid_y <- pretty(c(0, ymax / 1.32), n = 3)
  grid_y <- grid_y[grid_y > 0 & grid_y < ymax]
  # Stages with too few sites: darker band in every panel.
  no_value <- d %>% filter(is.na(beta_SIM)) %>% distinct(Stage_Number) %>%
    mutate(xmin = stage_young[Stage_Number], xmax = stage_older[Stage_Number])
  row_lab <- data.frame(Component = factor(comp_levels, levels = comp_levels),
                        Diet = factor(keep[1], levels = keep))
  n_col <- paste0("n_", unit_word)
  n_tab <- distinct(partition %>% filter(Province == p), Stage_Number, .data[[n_col]])
  inset <- 0.05

  ggplot(d) +
    geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill)) +
    geom_rect(data = no_value, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "#dcdcd7", inherit.aes = FALSE) +
    scale_fill_identity() +
    geom_hline(yintercept = grid_y, colour = grid_col, linewidth = 0.4) +
    geom_rect(aes(xmin = stage_young[Stage_Number] + inset, xmax = stage_older[Stage_Number] - inset,
                  ymin = 0, ymax = value, fill = bar_fill),
              colour = ink, linewidth = 0.3, na.rm = TRUE) +
    geom_text(data = row_lab, aes(x = 4.6, y = ymax * 0.97, label = Component),
              hjust = 0, vjust = 1, size = base_size / 3.6, colour = ink) +
    facet_grid(Component ~ Diet) +
    scale_x_reverse(limits = c(4.7, 0.129), breaks = c(4, 3, 2, 1), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, ymax), breaks = c(0, grid_y),
                       expand = c(0, 0)) +
    labs(title = paste0(p, ": endemic versus shared species by diet"),
         subtitle = paste0("Proportion of shared (top) and endemic (bottom) species among ", unit_word,
                           ", stage by stage\nFalling endemic bars + rising shared bars through time = homogenization"),
         x = "Age (Ma)", y = "Proportion",
         caption = paste0(
           "After Rowan et al. (2024, Fig. 4a; Eqs. 5-6), mean of subsets of ", n_units_sample, " ",
           unit_word, ". Background bands = stages, Zanclean (left) to Chibanian (right); darker band = fewer than ",
           n_units_sample, " ", unit_word, ".\n", tools::toTitleCase(unit_word), " per stage: ",
           paste(sprintf("%s %d", substr(stage_names, 1, 4), n_tab[[n_col]][order(n_tab$Stage_Number)]),
                 collapse = ", "),
           ". Unclassified species are not drawn. Diets: Smith et al. (aao5987). Silhouettes: game-icons.net (CC BY 3.0).")) +
    theme_panels() +
    theme(strip.text.y = element_blank(), strip.background.y = element_blank(),
          strip.text.x = element_text(face = "bold", size = base_size * 0.64, hjust = 0,
                                      margin = margin(t = 10, b = 10, l = 6)),
          axis.text = element_text(colour = ink, size = base_size * 0.62),
          panel.spacing.x = unit(0.6, "lines"), panel.spacing.y = unit(0.8, "lines"))
}

# Silhouettes in the column titles of a facet_grid figure.
add_silhouettes_grid <- function(plt, height_in) {
  g <- ggplotGrob(plt)
  lay <- ggplot_build(plt)$layout$layout
  for (col in sort(unique(lay$COL))) {
    diet <- as.character(lay$Diet[lay$COL == col][1])
    pos <- g$layout[g$layout$name %in% c(sprintf("strip-t-%d", col), sprintf("strip-t-%d-1", col)), ]
    if (nrow(pos) == 0) next
    pos <- pos[which.min(pos$t), ]
    g <- gtable::gtable_add_grob(g, silhouette_grob(diet, diet_colours[[diet]], height_in),
                                 t = pos$t, l = pos$l, b = pos$b, r = pos$r,
                                 z = Inf, clip = "off", name = paste0("silhouette-col-", col))
  }
  g
}

for (p in provinces) {
  plt <- es_plot(p)
  if (is.null(plt)) next
  n_panels <- length(unique(plt$data$Diet))
  w <- if (n_panels <= 3) slide_w * 0.75 else slide_w
  save_fig(add_silhouettes_grid(plt, 0.3), paste0("endemic_shared_", file_stub(p)),
           width = w, height = slide_h)
  cat(sprintf("  %-16s endemic/shared figure (%d diet groups)\n", p, n_panels))
}

cat("\n=== STEP 8b COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: partition, es_table\n")
