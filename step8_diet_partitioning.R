# =============================================================================
# STEP 8: ADDITIVE PARTITIONING OF beta_SIM BY DIETARY GROUP
#         (Rowan et al. 2024, Nat. Ecol. Evol., Fig. 3a and Eq. 4)
#
# Step 7 (analysis A) gives, for each province and stage, the multisite
# Simpson dissimilarity (beta_SIM) among the sites of that province. Here that
# value is split into the contributions of six dietary groups:
#
#   Carnivore, Omnivore, Frugivore, Browser, Mixed feeder, Grazer
#
#   beta_SIM_f = sum_{i<j} min(b_ij, b_ji)_f /
#                [ sum_i S_i - S_T + sum_{i<j} min(b_ij, b_ji) ]          (Eq. 4)
#
#   min(b_ij, b_ji)_f = species of diet f among the endemics of whichever unit
#   (i or j) has fewer endemics. When both have the same number, the two sets
#   are averaged. The denominator is the full one, so the groups ADD UP to the
#   overall value: sum_f beta_SIM_f = beta_SIM. A tall bar means that diet group
#   supplies many of the species that differ between sites.
#
# Diet of each species: every species in the four regional pools is looked up
# by name in Smith et al. (Science, aao5987, Table S7; file
# aao5987-smith-sm-137-190.xlsx); other rows of the table are not used. The
# species's row is searched (every cell of the row, not only one column) for
# the terms carnivore, omnivore, frugivore, browser, mixed feeder, grazer.
# Rules when a row holds several terms:
#   - "mixed feeder", or grazer + browser, or grazer + frugivore -> Mixed feeder
#   - a plant term (grazer/browser/frugivore) + carnivore/omnivore -> Omnivore
#   - grazer -> Grazer; browser -> Browser; frugivore -> Frugivore
#   - carnivore + omnivore -> the table's own "Recoded Diet" (carnivore or
#     omnivore), otherwise Omnivore
#   - carnivore -> Carnivore; omnivore -> Omnivore
#   - none of the terms (e.g. "herbivore", "insectivore") -> unclassified
# Pool species not in the table stay unclassified (set genus_fallback <- TRUE
# to give them the most common diet of their congeners instead). Unclassified
# species still count in the overall beta_SIM, so the grey line equals Step 7.
# You can set any species by hand in diet_overrides.csv (see SETTINGS).
#
# A diet group with no species in a province is left out of that province's
# figure.
#
# Inputs:
#   Outputs/6_matrices/all_records_final.csv              (Step 6)
#   Outputs/maps/physiographic/site_physio_database.csv   (Step 5)
#   aao5987-smith-sm-137-190.xlsx                         (Smith et al. diet table)
#   Outputs/7_regional_pools/beta_spatial.csv             (Step 7; used only as a check)
#
# Outputs (Outputs/8_diet_partitioning/):
#   diet_partition_<province>.png/.pdf   Fig. 3a-style figure, one per province
#   diet_partition.csv          beta_SIM_f for every province, stage and diet group
#   diet_trend_summary.csv      one row per province x diet: values per stage, change
#   species_diet.csv            every pool species with its diet group and the
#                               terms found in its row of the table
#   unclassified_species.csv    pool species without a group (fill in, save as
#                               diet_overrides.csv, run again)
#
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 7.
# Needs: dplyr, tidyr, ggplot2, readxl (gtable and grid come with ggplot2).
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "readxl", "gtable")) {
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

records_file <- file.path("Outputs", "6_matrices", "all_records_final.csv")
physio_file  <- file.path("Outputs", "maps", "physiographic", "site_physio_database.csv")
step7_file   <- file.path("Outputs", "7_regional_pools", "beta_spatial.csv")
output_dir   <- file.path(work_dir, "Outputs", "8_diet_partitioning")

# Smith et al. diet table. Put the file in work_dir, or give its full path here.
# If it is not found, the Downloads folder is searched as well.
diet_file <- file.path(work_dir, "aao5987-smith-sm-137-190.xlsx")
downloads_dir <- "C:/Users/shrut/Downloads"

# Optional: your own diet assignments (they win over the table). A CSV in
# work_dir with two columns: GenusSpecies, Diet_Group (one of diet_groups).
diet_overrides_file <- file.path(work_dir, "diet_overrides.csv")

genus_fallback <- FALSE    # TRUE: pool species missing from the table take the most
                           # common diet of their congeners in the table

provinces <- c("Pacific Border", "Basin and Range", "Great Plains", "Coastal Plain")

diet_groups <- c("Carnivore", "Omnivore", "Frugivore", "Browser", "Mixed feeder", "Grazer")

# Rowan et al. (2024) drew browsers and frugivores as one panel. TRUE copies that.
merge_browser_frugivore <- FALSE

# Keep these three identical to Step 7 so the overall values match Step 7.
spatial_unit   <- "site"   # "site" or "section"
n_units_sample <- 3
n_resamples    <- 999
set.seed(2024)

# One colour per diet group (validated colour-blind-safe set).
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

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================

read_input <- function(file) {
  path <- file.path(work_dir, file)
  if (!file.exists(path)) stop("File not found:\n  ", path, "\nRun Steps 5-7 first.")
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

# "Equus cf. simplicidens" -> "equus simplicidens"; "Equus sp." -> "equus".
clean_name <- function(x) {
  x <- tolower(x)
  x <- gsub("\\b(cf|aff|nr|sp|spp|indet|n|gen|ex|gr)\\b\\.?", " ", x)
  x <- gsub("[^a-z ]", " ", x)
  x <- gsub("\\s+", " ", trimws(x))
  w <- strsplit(x, " ", fixed = TRUE)
  vapply(w, function(v) paste(head(v, 2), collapse = " "), character(1))
}
genus_of <- function(x) sub(" .*$", "", x)

# Multisite Simpson dissimilarity (Rowan et al. 2024 Eq. 1).
beta_sim_multi <- function(pools) {
  pools <- pools[lengths(pools) > 0]
  k <- length(pools)
  if (k < 2) return(NA_real_)
  S_T <- length(unique(unlist(pools)))
  min_sum <- 0
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    min_sum <- min_sum + min(length(setdiff(pools[[i]], pools[[j]])),
                             length(setdiff(pools[[j]], pools[[i]])))
  }
  min_sum / (sum(lengths(pools)) - S_T + min_sum)
}

# Additive partition of multisite beta_SIM by category (Rowan et al. 2024 Eq. 4).
# Returns the share of each category; the shares add up to beta_SIM.
beta_sim_partition <- function(pools, category_of, categories) {
  k <- length(pools)
  count <- function(sp) tabulate(match(category_of[sp], categories), nbins = length(categories))
  S_T <- length(unique(unlist(pools)))
  num <- numeric(length(categories)); min_sum <- 0
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    b_ij <- setdiff(pools[[i]], pools[[j]]); b_ji <- setdiff(pools[[j]], pools[[i]])
    num <- num + if (length(b_ij) < length(b_ji)) count(b_ij) else
      if (length(b_ij) > length(b_ji)) count(b_ji) else (count(b_ij) + count(b_ji)) / 2
    min_sum <- min_sum + min(length(b_ij), length(b_ji))
  }
  setNames(num / (sum(lengths(pools)) - S_T + min_sum), categories)
}

# Subsets of n units, exactly as in Step 7 (all subsets if few, else random).
unit_subsets <- function(k, n) {
  if (k < n || n < 2) return(list())
  if (choose(k, n) <= n_resamples) combn(k, n, simplify = FALSE) else
    replicate(n_resamples, sample.int(k, n), simplify = FALSE)
}

beta_lab <- quote(beta[SIM])

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
silhouette_grob <- function(diet, colour) {
  key <- if (diet == "Browser and frugivore") "Browser" else diet
  if (!key %in% names(silhouette_paths)) return(nullGrob())
  parts <- strsplit(silhouette_paths[[key]], ";", fixed = TRUE)[[1]]
  xy <- do.call(rbind, lapply(seq_along(parts), function(i) {
    v <- as.numeric(strsplit(trimws(parts[i]), " +")[[1]])
    data.frame(x = v[c(TRUE, FALSE)], y = v[c(FALSE, TRUE)], id = i)
  }))
  w <- diff(range(xy$x)); h <- diff(range(xy$y))
  vp <- viewport(x = unit(1, "npc") - unit(10, "pt"), y = unit(0.5, "npc"),
                 width = unit(icon_height_in * w / h, "in"), height = unit(icon_height_in, "in"),
                 just = c("right", "centre"))
  pathGrob((xy$x - min(xy$x)) / w, 1 - (xy$y - min(xy$y)) / h, id = xy$id,
           rule = "evenodd", gp = gpar(fill = colour, col = NA), vp = vp)
}

# =============================================================================
# 1. READ RECORDS AND LINK SITES TO PROVINCES (same rules as Step 7)
# =============================================================================

cat("=== 1. READING RECORDS ===\n")
records <- read_input(records_file)
physio  <- read_input(physio_file)

site_province_all <- physio %>%
  transmute(Stage_Number = as.integer(Stage_Number), Database,
            SiteName = trimws(SiteName), Province = trimws(US_Province),
            Section = if ("US_Section" %in% names(physio)) trimws(US_Section) else NA_character_) %>%
  filter(!is.na(SiteName), !is.na(Province))

site_province <- site_province_all %>%
  count(Stage_Number, Database, SiteName, Province) %>%
  group_by(Stage_Number, Database, SiteName) %>%
  arrange(desc(n), Province, .by_group = TRUE) %>%
  ungroup() %>%
  distinct(Stage_Number, Database, SiteName, .keep_all = TRUE) %>%
  select(Stage_Number, Database, SiteName, Province)

site_section <- site_province_all %>%
  filter(!is.na(Section)) %>%
  count(Stage_Number, Database, SiteName, Province, Section) %>%
  group_by(Stage_Number, Database, SiteName, Province) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-n)
site_province <- left_join(site_province, site_section,
                           by = c("Stage_Number", "Database", "SiteName", "Province"))

recs <- records %>%
  transmute(Stage_Number = as.integer(Stage_Number), Database,
            SiteName = trimws(SiteName), GenusSpecies, Order) %>%
  filter(!is.na(GenusSpecies)) %>%
  inner_join(site_province, by = c("Stage_Number", "Database", "SiteName")) %>%
  mutate(Province = provinces[match(tolower(Province), tolower(provinces))]) %>%
  filter(!is.na(Province), Stage_Number %in% 1:5)

# The species of the four regional pools (all stages together).
species_diet <- recs %>%
  group_by(GenusSpecies) %>%
  summarise(Order = names(sort(table(Order), decreasing = TRUE))[1],
            n_records = n(), Provinces = paste(sort(unique(Province)), collapse = "; "),
            .groups = "drop") %>%
  mutate(Name_key = clean_name(GenusSpecies))

cat(sprintf("  Species records in the four provinces: %d (%d species in the regional pools)\n",
            nrow(recs), nrow(species_diet)))

# =============================================================================
# 2. LOOK UP THE POOL SPECIES IN THE SMITH ET AL. TABLE
# =============================================================================

cat("\n=== 2. DIET OF THE POOL SPECIES (Smith et al. table) ===\n")

if (!file.exists(diet_file)) {
  hit <- list.files(downloads_dir, pattern = "aao5987.*\\.xlsx$", full.names = TRUE,
                    ignore.case = TRUE)
  if (length(hit) == 0) {
    stop("Diet table not found:\n  ", diet_file,
         "\nPut aao5987-smith-sm-137-190.xlsx in your work folder or set diet_file.")
  }
  diet_file <- hit[1]
}
cat("  Reading", diet_file, "\n")

raw <- suppressMessages(readxl::read_excel(diet_file, col_names = FALSE, col_types = "text"))
raw <- as.data.frame(raw, stringsAsFactors = FALSE)

# Header row = the row with a cell "Genus_species"; data start below it.
header_row <- which(apply(raw, 1, function(r) any(trimws(r) == "Genus_species", na.rm = TRUE)))[1]
if (is.na(header_row)) stop("No 'Genus_species' column found in ", basename(diet_file))
species_col <- which(trimws(unlist(raw[header_row, ])) == "Genus_species")[1]

# Column titles span the rows just above and below "Genus_species" (e.g.
# "Recoded Diet: ..."). The long table caption at the top is not used.
hdr_rows <- max(1, header_row - 2):min(nrow(raw), header_row + 2)
header_text <- vapply(seq_along(raw), function(j) {
  paste(na.omit(raw[hdr_rows, j]), collapse = " ")
}, character(1))
recoded_col <- grep("Recoded Diet", header_text, ignore.case = TRUE)[1]
pbdb_col    <- grep("PBDB diet", header_text, ignore.case = TRUE)[1]
cat(sprintf("  Columns: species = '%s', PBDB diet = '%s', recoded diet = '%s'\n",
            header_text[species_col],
            if (is.na(pbdb_col)) "not found" else header_text[pbdb_col],
            if (is.na(recoded_col)) "not found" else gsub("\\s+", " ", header_text[recoded_col])))

smith <- raw[(header_row + 1):nrow(raw), , drop = FALSE]
smith <- smith[!is.na(smith[[species_col]]) & trimws(smith[[species_col]]) != "", , drop = FALSE]
smith_key <- clean_name(smith[[species_col]])

# Only the rows of species that are in the regional pools are used
# (plus their congeners when genus_fallback is TRUE).
keep <- smith_key %in% species_diet$Name_key |
  (genus_fallback & genus_of(smith_key) %in% genus_of(species_diet$Name_key))
smith <- smith[keep, , drop = FALSE]
smith_key <- smith_key[keep]

# Search every cell of each row for the six diet terms.
row_text <- tolower(apply(smith, 1, function(r) paste(na.omit(r), collapse = " | ")))
has <- function(term) grepl(paste0("\\b", term, "s?\\b"), row_text)
t_carn <- has("carnivore"); t_omni <- has("omnivore"); t_frug <- has("frugivore")
t_brow <- has("browser");   t_graz <- has("grazer");   t_mix  <- grepl("mixed[- ]?feed", row_text)
recoded <- if (is.na(recoded_col)) rep(NA_character_, nrow(smith)) else
  tolower(trimws(smith[[recoded_col]]))
plant <- t_graz | t_brow | t_frug

smith_diet <- data.frame(
  Genus_species = trimws(smith[[species_col]]),
  Name_key = smith_key,
  PBDB_diet = if (is.na(pbdb_col)) NA_character_ else smith[[pbdb_col]],
  Recoded_diet = if (is.na(recoded_col)) NA_character_ else smith[[recoded_col]],
  Terms_found = vapply(seq_len(nrow(smith)), function(i) {
    paste(c("carnivore", "omnivore", "frugivore", "browser", "grazer", "mixed feeder")[
      c(t_carn[i], t_omni[i], t_frug[i], t_brow[i], t_graz[i], t_mix[i])], collapse = ", ")
  }, character(1)),
  stringsAsFactors = FALSE
) %>%
  mutate(Diet_Group = case_when(
    t_mix | (t_graz & t_brow) | (t_graz & t_frug) ~ "Mixed feeder",
    plant & (t_carn | t_omni)                     ~ "Omnivore",
    t_graz                                        ~ "Grazer",
    t_brow                                        ~ "Browser",
    t_frug                                        ~ "Frugivore",
    t_carn & t_omni & recoded %in% "carnivore"    ~ "Carnivore",
    t_carn & t_omni                               ~ "Omnivore",
    t_carn                                        ~ "Carnivore",
    t_omni                                        ~ "Omnivore",
    TRUE                                          ~ NA_character_))

# a) species-level match
by_species <- smith_diet %>%
  group_by(Name_key) %>%
  summarise(Matched_name = first(Genus_species), Diet_species = first(Diet_Group),
            Smith_PBDB_diet = first(PBDB_diet), Smith_recoded = first(Recoded_diet),
            Terms_found = first(Terms_found), .groups = "drop")
species_diet <- left_join(species_diet, by_species, by = "Name_key")
cat(sprintf("  %d of %d pool species found in the table by name\n",
            sum(!is.na(species_diet$Matched_name)), nrow(species_diet)))

# b) genus fallback (only if genus_fallback is TRUE): most common group among
#    congeners in the table; ties stay unclassified
genus_diet <- smith_diet %>%
  filter(!is.na(Diet_Group)) %>%
  mutate(Genus = genus_of(Name_key)) %>%
  count(Genus, Diet_Group) %>%
  group_by(Genus) %>%
  mutate(n_genus = sum(n)) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  summarise(Diet_genus = if (n() > 1 && n[1] == n[2]) NA_character_ else Diet_Group[1],
            Genus_support = sprintf("%d of %d congeners", n[1], n_genus[1]), .groups = "drop")
species_diet <- species_diet %>%
  mutate(Genus = genus_of(Name_key)) %>%
  left_join(genus_diet, by = "Genus")

# c) your own assignments
overrides <- NULL
if (file.exists(diet_overrides_file)) {
  overrides <- read.csv(diet_overrides_file, stringsAsFactors = FALSE, check.names = FALSE)
  names(overrides) <- trimws(sub("^[^A-Za-z0-9]+", "", names(overrides), useBytes = TRUE))
  overrides <- overrides %>%
    transmute(GenusSpecies = trimws(GenusSpecies), Diet_override = trimws(Diet_Group)) %>%
    filter(Diet_override %in% diet_groups)
  cat(sprintf("  %d diet overrides read from %s\n", nrow(overrides), basename(diet_overrides_file)))
  species_diet <- left_join(species_diet, overrides, by = "GenusSpecies")
} else {
  species_diet$Diet_override <- NA_character_
}

species_diet <- species_diet %>%
  mutate(
    Diet_Group = case_when(
      !is.na(Diet_override)                             ~ Diet_override,
      !is.na(Matched_name)                              ~ Diet_species,
      genus_fallback & !is.na(Diet_genus)               ~ Diet_genus,
      TRUE                                              ~ NA_character_),
    Diet_Source = case_when(
      !is.na(Diet_override)                             ~ "your override",
      !is.na(Matched_name) & !is.na(Diet_species)       ~ "species (Smith et al.)",
      !is.na(Matched_name)                              ~ "species in table, no diet term",
      genus_fallback & !is.na(Diet_genus)               ~ paste0("genus (", Genus_support, ")"),
      TRUE                                              ~ "not in table")
  )
if (merge_browser_frugivore) {
  species_diet$Diet_Group[species_diet$Diet_Group %in% c("Browser", "Frugivore")] <-
    "Browser and frugivore"
}

diet_levels <- if (merge_browser_frugivore) {
  c("Carnivore", "Omnivore", "Browser and frugivore", "Mixed feeder", "Grazer")
} else diet_groups
categories <- c(diet_levels, "Unclassified")

species_diet <- species_diet %>%
  mutate(Diet_Group = factor(coalesce(Diet_Group, "Unclassified"), levels = categories)) %>%
  select(GenusSpecies, Order, Diet_Group, Diet_Source, Matched_name, Smith_PBDB_diet,
         Smith_recoded, Terms_found, n_records, Provinces) %>%
  arrange(Diet_Group, Order, GenusSpecies)

save_csv(species_diet, "species_diet")
unclassified <- filter(species_diet, Diet_Group == "Unclassified") %>%
  transmute(GenusSpecies, Order, Diet_Source, Smith_PBDB_diet, Smith_recoded, n_records,
            Provinces, Diet_Group = "")
save_csv(unclassified, "unclassified_species")

cat("  Pool species per diet group:\n")
print(as.data.frame(count(species_diet, Diet_Group, name = "n_species")), row.names = FALSE)
cat("\n  How the diet was found:\n")
print(as.data.frame(count(species_diet, Source = sub(" \\(.*", "", Diet_Source),
                          name = "n_species")), row.names = FALSE)
cat(sprintf(paste0("\n  %d pool species are unclassified -> unclassified_species.csv.\n",
                   "  To assign them, fill in its Diet_Group column, save it as diet_overrides.csv\n",
                   "  in your work folder and run this step again.\n"), nrow(unclassified)))

category_of <- setNames(as.character(species_diet$Diet_Group), species_diet$GenusSpecies)

# =============================================================================
# 4. PARTITION beta_SIM BY DIET (Rowan et al. 2024 Eq. 4)
# =============================================================================

unit_word <- if (spatial_unit == "section") "sections" else "sites"
cat(sprintf("\n=== 4. PARTITIONING beta_SIM AMONG %s BY DIET ===\n", toupper(unit_word)))

recs_units <- recs %>%
  mutate(Unit = if (spatial_unit == "section") Section else paste(Database, SiteName, sep = " | ")) %>%
  filter(!is.na(Unit))

# Same loop order and subsets as Step 7, so the overall values are identical.
partition <- bind_rows(lapply(provinces, function(p) bind_rows(lapply(1:5, function(s) {
  d <- filter(recs_units, Province == p, Stage_Number == s)
  units <- lapply(split(d$GenusSpecies, d$Unit), unique)
  subsets <- unit_subsets(length(units), n_units_sample)
  n_by_diet <- tabulate(match(category_of[unique(d$GenusSpecies)], categories),
                        nbins = length(categories))
  if (length(subsets) == 0) {
    return(data.frame(Province = p, Stage_Number = s, Diet = categories, beta_SIM_f = NA_real_,
                      beta_SIM_within_group = NA_real_, n_species = n_by_diet,
                      n_units = length(units), beta_SIM = NA_real_, beta_SIM_lo95 = NA_real_,
                      beta_SIM_hi95 = NA_real_))
  }
  parts <- sapply(subsets, function(ix) beta_sim_partition(units[ix], category_of, categories))
  totals <- colSums(parts)
  # Dissimilarity of each group on its own (same subsets; for reference).
  within <- vapply(categories, function(f) {
    v <- vapply(subsets, function(ix) {
      beta_sim_multi(lapply(units[ix], function(u) u[category_of[u] == f]))
    }, numeric(1))
    if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  }, numeric(1))
  data.frame(Province = p, Stage_Number = s, Diet = categories,
             beta_SIM_f = rowMeans(parts), beta_SIM_within_group = unname(within),
             n_species = n_by_diet, n_units = length(units),
             beta_SIM = mean(totals),
             beta_SIM_lo95 = unname(quantile(totals, 0.025)),
             beta_SIM_hi95 = unname(quantile(totals, 0.975)))
}))))

partition <- partition %>%
  mutate(Stage = factor(stage_names[Stage_Number], levels = stage_names),
         Mid_Ma = (stage_older[Stage_Number] + stage_young[Stage_Number]) / 2,
         Province = factor(Province, levels = provinces),
         Diet = factor(Diet, levels = categories)) %>%
  relocate(Stage, .after = Stage_Number)
names(partition)[names(partition) == "n_units"] <- paste0("n_", unit_word)

save_csv(mutate(partition, across(where(is.numeric), ~ round(.x, 4))), "diet_partition")

# Check: the overall values should equal Step 7 (beta_spatial.csv).
step7_path <- file.path(work_dir, step7_file)
if (file.exists(step7_path)) {
  s7 <- read.csv(step7_path, stringsAsFactors = FALSE, check.names = FALSE)
  chk <- partition %>% distinct(Province, Stage, beta_SIM) %>%
    mutate(Province = as.character(Province), Stage = as.character(Stage)) %>%
    inner_join(transmute(s7, Province, Stage, beta_step7 = beta_SIM), by = c("Province", "Stage"))
  dif <- suppressWarnings(max(abs(chk$beta_SIM - chk$beta_step7), na.rm = TRUE))
  cat(sprintf("  Check against Step 7: largest difference in overall beta_SIM = %s\n",
              if (is.finite(dif)) format(round(dif, 4), nsmall = 4) else "n/a"))
  if (is.finite(dif) && dif > 0.01) {
    cat("  NOTE: differs from Step 7 - make sure spatial_unit, n_units_sample,\n",
        "        n_resamples and the seed are the same in both scripts.\n")
  }
}

# Trend summary per province and diet group.
diet_present <- partition %>%
  group_by(Province, Diet) %>%
  summarise(any_species = sum(n_species) > 0, .groups = "drop")

trend_summary <- partition %>%
  filter(Diet != "Unclassified") %>%
  select(Province, Diet, Stage, beta_SIM_f) %>%
  mutate(beta_SIM_f = round(beta_SIM_f, 3)) %>%
  pivot_wider(names_from = Stage, values_from = beta_SIM_f) %>%
  left_join(diet_present, by = c("Province", "Diet")) %>%
  filter(any_species) %>%
  select(-any_species)

# Change from the first to the last stage that has a value.
first_last <- t(apply(as.matrix(trend_summary[stage_names]), 1, function(v) {
  ok <- which(!is.na(v))
  if (length(ok) < 2) return(c(NA, NA, NA))
  c(stage_names[min(ok)], stage_names[max(ok)], v[max(ok)] - v[min(ok)])
}))
trend_summary <- trend_summary %>%
  mutate(From = first_last[, 1], To = first_last[, 2],
         Change = round(as.numeric(first_last[, 3]), 3),
         Direction = case_when(
           is.na(Change)  ~ "fewer than 2 stages with a value",
           Change > 0.02  ~ "larger share of dissimilarity",
           Change < -0.02 ~ "smaller share of dissimilarity",
           TRUE           ~ "little change (< 0.02)"))

save_csv(trend_summary, "diet_trend_summary")
cat("  Share of beta_SIM per diet group and stage (bars in the figures):\n")
print(as.data.frame(select(trend_summary, Province, Diet, all_of(stage_names), Change)),
      row.names = FALSE)

# =============================================================================
# 5. FIGURES: ONE FIG. 3a-STYLE FIGURE PER PROVINCE
# =============================================================================

cat("\n=== 5. FIGURES ===\n")

bands <- data.frame(xmin = stage_young, xmax = stage_older, Stage = stage_names,
                    lab = substr(stage_names, 1, 4),
                    fill = rep(c("#f3f3f0", "#ffffff"), length.out = 5))

# Line segments only between neighbouring stages that both have a value.
gap_segments <- function(d) {
  d %>% arrange(Diet, Stage_Number) %>% group_by(Diet) %>%
    mutate(x2 = lead(Mid_Ma), y2 = lead(beta_SIM)) %>% ungroup() %>%
    filter(!is.na(beta_SIM), !is.na(y2))
}

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

diet_plot <- function(p) {
  keep <- diet_present %>% filter(Province == p, any_species, Diet != "Unclassified") %>%
    pull(Diet) %>% as.character()
  if (length(keep) == 0) return(NULL)
  d <- partition %>% filter(Province == p, Diet %in% keep) %>%
    mutate(Diet = factor(as.character(Diet), levels = keep),
           bar_fill = diet_colours[as.character(Diet)])
  overall <- distinct(d, Diet, Stage_Number, Mid_Ma, beta_SIM, beta_SIM_lo95, beta_SIM_hi95)
  n_col <- paste0("n_", unit_word)
  n_tab <- distinct(partition %>% filter(Province == p), Stage_Number, .data[[n_col]])
  unclass_share <- partition %>% filter(Province == p) %>%
    group_by(Stage_Number) %>%
    summarise(u = sum(beta_SIM_f[Diet == "Unclassified"]), t = first(beta_SIM), .groups = "drop") %>%
    filter(!is.na(t), t > 0)
  share_txt <- if (nrow(unclass_share)) sprintf("%.0f%%", 100 * sum(unclass_share$u) / sum(unclass_share$t)) else "n/a"
  inset <- 0.05
  ncol_f <- if (length(keep) <= 3) length(keep) else if (length(keep) == 4) 2 else 3

  ggplot(d) +
    geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill)) +
    scale_fill_identity() +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), colour = grid_col, linewidth = 0.4) +
    geom_text(data = bands, aes(x = (xmin + xmax) / 2, y = 1.09, label = lab),
              size = base_size / 4.6, colour = ink_soft) +
    # Grey: overall beta_SIM (all species), 95% range of the site subsets.
    geom_ribbon(data = overall, aes(x = Mid_Ma, ymin = beta_SIM_lo95, ymax = beta_SIM_hi95),
                fill = overall_band, alpha = 0.55) +
    # Coloured bars: this group's share of beta_SIM in each stage.
    geom_rect(aes(xmin = stage_young[Stage_Number] + inset, xmax = stage_older[Stage_Number] - inset,
                  ymin = 0, ymax = beta_SIM_f, fill = bar_fill),
              colour = NA, na.rm = TRUE) +
    geom_segment(data = gap_segments(overall),
                 aes(x = Mid_Ma, xend = x2, y = beta_SIM, yend = y2),
                 colour = overall_line, linewidth = 1) +
    geom_point(data = overall, aes(Mid_Ma, beta_SIM), shape = 21, size = 3.4,
               fill = "#9a9a96", colour = "white", stroke = 1.2, na.rm = TRUE) +
    facet_wrap(~ Diet, ncol = ncol_f) +
    scale_x_reverse(limits = c(4.7, 0.129), breaks = c(4, 3, 2, 1), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25),
                       labels = c("0", "0.25", "0.50", "0.75", "1.00"), expand = c(0, 0)) +
    labs(title = paste0(p, ": which diets make sites differ?"),
         subtitle = paste0("Coloured bars = each diet group's share of multisite Simpson dissimilarity among ",
                           unit_word, " (shares add up to the grey line)\nGrey = all species ",
                           "(band = 95% range of subsets of ", n_units_sample, " ", unit_word, ")"),
         x = "Age (Ma)", y = beta_lab,
         caption = paste0(
           "Additive partition of beta_SIM by diet after Rowan et al. (2024, Eq. 4). Diets from Smith et al. ",
           "(aao5987, Table S7). ", tools::toTitleCase(unit_word), " per stage: ",
           paste(sprintf("%s %d", substr(stage_names, 1, 4), n_tab[[n_col]][order(n_tab$Stage_Number)]),
                 collapse = ", "), ".\n",
           "Unclassified species (no diet term) supply ", share_txt,
           " of beta_SIM and are not drawn. Silhouettes: game-icons.net (CC BY 3.0).")) +
    theme_panels()
}

# Put each group's silhouette in the right-hand end of its panel title.
add_silhouettes <- function(plt) {
  g <- ggplotGrob(plt)
  lay <- ggplot_build(plt)$layout$layout
  for (r in seq_len(nrow(lay))) {
    diet <- as.character(lay$Diet[r])
    nm <- sprintf("strip-t-%d-%d", lay$COL[r], lay$ROW[r])
    pos <- g$layout[g$layout$name == nm, ]
    if (nrow(pos) != 1) next
    g <- gtable::gtable_add_grob(g, silhouette_grob(diet, diet_colours[[diet]]),
                                 t = pos$t, l = pos$l, b = pos$b, r = pos$r,
                                 z = Inf, clip = "off", name = paste0("silhouette-", r))
  }
  g
}

for (p in provinces) {
  plt <- diet_plot(p)
  if (is.null(plt)) {
    cat(sprintf("  %-16s no classified species - no figure\n", p))
    next
  }
  n_panels <- length(unique(plt$data$Diet))
  h <- if (n_panels <= 3) slide_h * 0.8 else slide_h * 1.15
  fig <- add_silhouettes(plt)
  save_fig(fig, paste0("diet_partition_", file_stub(p)), width = slide_w, height = h)
  missing <- setdiff(diet_levels, unique(as.character(plt$data$Diet)))
  cat(sprintf("  %-16s %d diet panels%s\n", p, n_panels,
              if (length(missing)) paste0(" (no species: ", paste(missing, collapse = ", "), ")") else ""))
}

cat("\n=== STEP 8 COMPLETE ===\n")
cat("Outputs in:", output_dir, "\n")
cat("Objects in your Environment: species_diet, partition, trend_summary\n")
