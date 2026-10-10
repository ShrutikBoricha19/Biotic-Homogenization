# =============================================================================
# STEP 7: ROWAN ET AL. (2024) FIGURES 2, 3 AND 4 WITH OUR DATA
#
# Rowan, J. et al. 2024. Long-term biotic homogenization in the East African
# Rift System over the last 6 million years of hominin evolution. Nature
# Ecology & Evolution, doi:10.1038/s41559-024-02462-0.
#
# The three regions here play the role of Rowan's three EARS subregions:
#   Basin and Range, Coastal Plain, Great Plains (with the Central Lowland
#   west of the Mississippi, Step 5b)
# and the time bins are ours (Step 3c; 0.75 Myr from 3.25 Ma, the last bin
# shorter). For every time bin the sites of each region are pooled into one
# regional species list (a region x species presence-absence matrix), exactly
# as in Rowan et al.'s code (Rowan_et_al_SI_Code_14MAY.R).
#
# Figures (as in the paper):
#   Fig. 2  Multisite Simpson dissimilarity beta_SIM of the three regions per
#           time bin (Eq. 1).
#   Fig. 3  beta_SIM additively partitioned into diet groups (a) and body size
#           classes (b) (Eq. 4); the equal-species beta_SIM per bin (every region
#           cut to the poorest region's species count, mean of 'equal_draws'
#           random draws, as in step7_addendum_equal_species.R; yellow points
#           joined by straight lines) in every panel. The bars use all species.
#   Fig. 4  Proportions of shared (Eq. 6, top rows) and endemic (Eq. 5, bottom
#           rows) species within diet groups (a) and size classes (b); the
#           equal-species beta_SIM curve (as Fig. 3) in every panel.
#
# Species included (as Rowan et al.): large mammals - species of
# 'large_mammal_orders' heavier than 'min_mass_kg' (1 kg). Species of these
# orders without a body mass are kept (they count in beta_SIM but cannot be
# given a size class).
#
# Traits (Smith et al. 2018, Science, aao5987, Table S7 - the fixed workbook):
#   Body mass: exp(ln Mass (g)) / 1000 kg. Size classes as Rowan et al.:
#     Size 1 (<18 kg), Size 2 (18-80 kg), Size 3 (80-350 kg),
#     Size 4 (350-1,000 kg), Size 5 (>1,000 kg).
#   Diet (Rowan et al.'s five groups), from the table's recoded diet and its
#   PBDB diet terms:
#     carnivore (and insectivore)         -> Carnivore
#     omnivore                            -> Omnivore
#     herbivore + grazer and browser/frugivore/folivore terms -> Mixed feeder
#     herbivore + grazer only             -> Grazer
#     herbivore + browser/frugivore/folivore only -> Browser (no frugivores in
#     these data, so the group is called Browser)
#     herbivore without these terms       -> unclassified
#   A species not in the table takes the mean ln mass and the most common diet
#   of its genus in the table (flagged in species_traits.csv). Anything can be
#   set by hand in species_traits_overrides.csv (see SETTINGS).
# Unclassified species still count in the overall beta_SIM; they are listed in
# unclassified_species.csv so that they can be filled in.
#
# Outputs (Outputs/7_rowan_figures/):
#   Fig2_beta_sim.png/.pdf, Fig3_partitioned_beta_sim.png/.pdf,
#   Fig4_shared_endemic.png/.pdf   (300 dpi PNG and vector PDF; curve = equal-species beta_SIM)
#   Fig3_partitioned_beta_sim_full_data, Fig4_shared_endemic_full_data
#                                  the same with the full-data beta_SIM curve (as Fig. 2)
#   beta_sim_by_bin.csv          beta_SIM, endemic and shared components per bin
#   partition_by_group.csv       Eqs. 4-6 for every diet group and size class
#   regression_beta_sim_age.csv  OLS of beta_SIM on bin midpoint (as Rowan et al.)
#   species_traits.csv           every species with order, mass, size class, diet
#   species_functional_groups.xlsx  the same as a workbook: all species, by diet
#                                group, by size class, diet x size counts,
#                                unclassified species, notes
#   unclassified_species.csv     large-mammal species without a diet or size class
#   regional_pa_matrices.rds     the region x species matrices per bin
#
# Inputs: Outputs/5b_great_plains/master_data_unique.csv (Step 5b),
#         Outputs/3d_resolved/time_bins.csv, faunmap_fauna.csv,
#         pbdb_occurrences.csv (Step 3d, for each genus's order),
#         the fixed Smith et al. workbook (aao5987-smith-sm-137-190_fixed.xlsx).
# Run the whole file (Ctrl+Shift+S in RStudio) after Step 5b.
# Needs: dplyr, tidyr, ggplot2, readxl, patchwork (grid/gtable come with ggplot2);
#        writexl for the Excel workbook.
# Silhouettes: from PhyloPic (phylopic.org API; silhouette_source = "phylopic", cached in
#   Outputs/7_rowan_figures/phylopic/, credits in phylopic_credits.csv), else the built-in diet icons (as in Step 8)
#   and size icons from game-icons.net (CC BY 3.0).
# =============================================================================

for (pkg in c("dplyr", "tidyr", "ggplot2", "readxl", "patchwork", "gtable")) {
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

master_file <- file.path("Outputs", "5b_great_plains", "master_data_unique.csv")
bins_file   <- file.path("Outputs", "3d_resolved", "time_bins.csv")
fauna_file  <- file.path("Outputs", "3d_resolved", "faunmap_fauna.csv")
pbdb_file   <- file.path("Outputs", "3d_resolved", "pbdb_occurrences.csv")
output_dir  <- file.path(work_dir, "Outputs", "7_rowan_figures")

# Smith et al. table (the fixed workbook). Searched in work_dir, then Downloads.
smith_file    <- file.path(work_dir, "aao5987-smith-sm-137-190_fixed.xlsx")
downloads_dir <- "C:/Users/shrut/Downloads"
smith_sheet   <- "Table S7 fixed"

# Optional hand assignments (win over the table): CSV in work_dir with columns
# Species, Order, Mass_kg, Diet_Group (leave a cell empty to keep the table value).
overrides_file <- file.path(work_dir, "species_traits_overrides.csv")

regions <- c("Basin and Range", "Coastal Plain", "Great Plains")

# Large mammals as in Rowan et al. (their orders, with the North American
# xenarthrans in place of the African hyraxes and aardvarks).
large_mammal_orders <- c("Artiodactyla", "Carnivora", "Perissodactyla", "Proboscidea", "Primates",
                         "Xenarthra", "Cingulata", "Pilosa", "Edentata")
min_mass_kg    <- 1
genus_fallback <- TRUE           # species missing from the table take their genus's values
insectivore_as <- "Carnivore"    # Rowan et al. have no insectivore group

diet_groups <- c("Carnivore", "Omnivore", "Browser", "Mixed feeder", "Grazer")
size_breaks <- c(0, 18, 80, 350, 1000, Inf)          # kg
size_groups <- c("Size 1 (<18 kg)", "Size 2 (18-80 kg)", "Size 3 (80-350 kg)",
                 "Size 4 (350-1,000 kg)", "Size 5 (>1,000 kg)")

# Colourblind-friendly colours, no greys or whites in the plots.
# Diet groups: Okabe-Ito palette (checked for protan, deutan and tritan vision;
# every panel is also labelled and carries an icon). Size classes: viridis
# (ordered light -> dark from small to large; readable in all colour-vision types).
diet_colours <- c("Carnivore" = "#D55E00", "Omnivore" = "#CC79A7", "Browser" = "#0072B2",
                  "Mixed feeder" = "#009E73", "Grazer" = "#E69F00")
size_colours <- setNames(c("#FDE725", "#5EC962", "#21918C", "#3B528B", "#440154"), size_groups)
ink          <- "#1F2A44"        # text, axes, lines and outlines (dark navy instead of grey/black)
ink_soft     <- "#3E4C6D"        # secondary text
point_fill   <- "#F0E442"        # overall beta_SIM points (Okabe-Ito yellow, navy outline)
strip_fill   <- "#D6E6F5"        # panel titles
band_fill    <- "#F6EBD0"        # shading of alternate time bins (pale sand)
equal_draws  <- 999            # Figs. 3-4 curve: random draws per bin, as step7_addendum_equal_species.R
equal_seed   <- 2024           #   (same seed and method, so the values match the addendum exactly)
page_fill    <- "white"          # figure background (change here, e.g. "#F7FAFD", for a tinted page)
panel_letters <- FALSE           # TRUE = "a"/"b" tags on Figs. 3 and 4

# Silhouettes on Figs. 3-4: "phylopic" downloads them from PhyloPic (phylopic.org; needs internet
# the first time and the packages jsonlite and png, then works offline from the cache folder);
# "builtin" uses the embedded silhouettes below. A taxon is looked up by its exact name, then by
# its genus; to force a particular picture, put its PhyloPic image UUID in phylopic_uuid.
# Any taxon PhyloPic cannot supply falls back to the built-in silhouette.
silhouette_source <- "phylopic"
phylopic_diet <- c("Carnivore" = "Smilodon fatalis", "Omnivore" = "Ursus americanus",
                   "Browser" = "Odocoileus virginianus", "Mixed feeder" = "Camelops hesternus",
                   "Grazer" = "Bison bison")
phylopic_size <- c("Vulpes vulpes", "Odocoileus virginianus", "Camelops hesternus",
                   "Bison bison", "Mammuthus columbi")         # Size 1 ... Size 5
phylopic_uuid <- c()        # optional, e.g. c("Bison bison" = "<image uuid from phylopic.org/images/...>")
phylopic_refresh <- FALSE   # TRUE = download again even if a silhouette is already cached

# Figure size (inches) and look.
fig2_w <- 9;  fig2_h <- 5.2
fig3_w <- 16; fig3_h <- 7.2
fig4_w <- 16; fig4_h <- 11
fig_dpi <- 300
base_size <- 14

# Silhouettes (x y pairs on a 512 x 512 grid; ";" separates parts), embedded so
# no download is needed. Diet: as in Step 8. Size classes: fox, deer, camel,
# bison, mammoth (game-icons.net, CC BY 3.0).
diet_silhouette_paths <- c(
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

size_silhouette_paths <- c(
  "Size 1" = paste0(
    "214 23 210 53 210 62 212 79 219 99 224 108 233 119 234 120 246 107 236 95 229 78 227 59 22",
    "8 43 248 51 264 61 279 74 287 84 292 92 284 95 270 104 253 119 240 137 212 180 232 176 230",
    " 208 250 206 272 209 289 217 300 225 304 215 306 204 302 168 292 163 280 155 266 143 265 1",
    "39 266 135 271 133 296 140 301 144 308 150 312 157 314 164 318 194 319 206 315 223 306 237",
    " 314 256 352 256 358 237 350 223 346 206 347 193 351 162 355 153 362 145 369 140 384 136 3",
    "96 133 400 136 400 141 398 144 382 156 372 163 362 168 359 204 361 215 366 225 376 217 393",
    " 209 415 206 434 208 432 176 454 180 425 137 412 119 395 104 381 95 374 92 378 84 386 74 4",
    "01 61 417 51 438 43 438 59 436 78 429 95 420 107 430 120 432 119 441 108 446 99 453 79 455",
    " 62 455 53 450 23 434 27 418 32 388 47 362 66 340 86 324 86 284 53 259 37 239 29 214 23;23",
    "9 228 244 238 246 246 245 253 241 259 234 265 214 276 158 299 128 314 113 324 99 335 86 34",
    "8 77 359 67 374 62 385 59 395 57 410 60 428 66 444 71 451 84 464 100 475 119 483 141 487 1",
    "52 488 177 488 190 486 222 478 316 449 345 442 364 440 381 440 396 445 409 454 421 467 436",
    " 462 434 442 429 423 423 406 414 391 404 378 389 368 374 360 360 356 345 354 329 354 311 3",
    "55 293 358 257 367 187 391 166 397 152 399 138 399 126 398 115 394 109 389 105 380 104 369",
    " 111 378 119 383 129 386 139 387 150 386 162 383 188 374 259 345 288 336 302 334 302 313 2",
    "99 296 294 278 285 260 279 252 271 244 262 238 251 232 239 229;328 230 338 230 341 231 343",
    " 235 342 240 338 242 327 242 324 241 322 238 322 234 324 231 328 230;389 388 398 397 406 4",
    "09 411 422 416 437 406 431 393 426 378 423 363 423 345 425 325 429 276 443 302 418 282 427",
    " 262 431 244 433 227 430 308 390 247 401 258 393 272 385 299 375 319 371 330 370 350 371 3",
    "70 376 380 381 389 388"),
  "Size 2" = paste0(
    "139 20 123 25 129 42 135 56 145 71 155 81 166 88 178 92 194 94 211 94 251 89 272 90 293 95",
    " 318 108 326 106 330 104 331 97 331 88 326 76 320 67 313 60 294 48 284 44 279 60 290 65 29",
    "9 71 307 78 312 85 290 75 276 71 258 67 223 63 212 59 205 55 200 48 195 39 192 29 176 33 1",
    "79 44 184 54 191 65 200 72 211 77 196 77 183 76 169 71 162 66 156 59 147 43 139 20;229 24 ",
    "229 39 233 52 255 54 250 47 247 40 245 32 245 25;349 37 338 50 345 58 348 69 349 73 336 64",
    " 341 76 344 95 360 97 365 81 365 72 364 62 361 52 358 47 349 37;268 107 250 108 244 109 24",
    "0 113 236 120 248 132 256 137 267 142 288 147 305 148 303 158 303 178 302 190 298 202 295 ",
    "207 287 217 277 226 270 229 256 235 238 239 199 241 175 244 150 250 127 258 111 268 96 281",
    " 86 294 80 309 87 320 102 320 102 381 84 410 89 492 107 492 108 426 128 412 147 395 162 37",
    "7 172 359 179 339 182 316 179 291 173 273 188 267 196 289 198 303 199 316 198 329 195 341 ",
    "192 352 185 369 224 367 258 363 253 306 270 305 286 492 303 492 306 370 315 365 327 355 33",
    "5 345 343 332 352 315 359 294 367 258 371 223 372 190 348 182 423 182 431 165 424 165 415 ",
    "162 417 154 432 154 347 110 331 117 320 124 309 117 295 111 283 108 268 107;364 333 351 35",
    "4 335 374 382 395 361 443 376 450 409 385 365 348 364 333"),
  "Size 3" = paste0(
    "421 27 403 30 387 37 379 31 371 30 368 31 363 37 362 42 364 48 368 53 375 59 382 105 383 1",
    "28 382 142 380 148 374 158 365 165 356 167 349 165 341 160 335 153 322 132 297 82 282 60 2",
    "73 51 264 44 253 39 241 38 234 39 224 44 217 52 205 79 199 89 191 96 185 99 178 100 169 98",
    " 161 92 156 83 141 53 135 47 131 45 122 44 108 48 100 53 93 62 86 74 80 89 74 108 69 131 6",
    "0 132 49 137 43 143 35 154 28 170 21 190 17 216 16 236 16 269 23 353 30 276 35 295 42 313 ",
    "47 341 48 369 43 378 42 387 43 397 48 405 46 455 46 485 59 485 146 485 145 479 139 471 128",
    " 464 111 455 111 407 116 398 118 388 116 378 111 370 112 345 113 331 119 317 144 279 165 2",
    "82 188 282 211 279 235 275 236 360 232 369 230 379 232 389 236 398 238 414 239 436 237 485",
    " 295 485 294 477 288 471 280 465 261 455 259 410 260 395 264 381 263 371 261 364 264 314 2",
    "67 290 271 265 281 259 292 358 289 367 288 377 291 387 296 395 298 406 300 428 303 484 365",
    " 485 362 479 355 473 324 459 320 439 318 395 322 385 323 374 321 364 316 355 318 326 324 2",
    "83 328 262 335 240 354 239 370 236 383 231 395 226 405 219 412 212 419 204 427 185 430 175",
    " 432 153 434 107 435 103 437 100 440 97 448 95 458 96 483 99 490 99 494 96 496 88 496 81 4",
    "94 74 491 67 484 61 480 58 472 55 457 53 433 32 427 28 422 27;428 46 431 51 430 57 420 57 ",
    "409 52 418 48 428 46;81 302 81 304 89 328 93 366 94 372 89 380 88 390 90 399 94 407 93 471",
    " 89 466 83 461 76 457 68 454 67 404 71 395 72 386 70 377 66 369 66 331 69 316 81 302"),
  "Size 4" = paste0(
    "300 99 286 100 277 102 257 110 239 120 193 148 166 162 149 166 110 172 86 177 74 182 64 18",
    "9 56 200 54 207 48 235 37 263 24 274 21 277 20 282 20 299 23 307 26 312 41 302 41 290 43 2",
    "77 55 249 56 265 54 287 43 337 46 357 53 386 64 413 100 413 86 395 78 381 74 367 73 352 77",
    " 337 94 334 102 330 106 327 127 382 136 397 146 413 187 413 158 383 149 368 144 356 143 34",
    "9 142 332 143 322 177 324 192 324 209 322 232 316 229 352 230 356 236 369 243 380 267 413 ",
    "306 413 277 386 271 378 267 369 265 357 265 344 274 320 294 350 312 374 327 391 345 408 35",
    "2 413 388 413 353 380 339 363 328 343 325 328 325 320 354 327 360 327 366 326 375 347 383 ",
    "362 394 375 403 382 416 366 424 352 428 336 429 316 454 331 462 329 466 325 470 321 470 31",
    "7 467 286 478 266 492 252 492 237 489 223 480 199 467 176 464 173 453 170 417 164 392 148 ",
    "349 117 335 109 314 100 300 99;373 134 388 163 397 174 409 184 430 196 438 204 441 209 442",
    " 217 441 225 438 231 432 238 427 242 421 245 416 244 399 238 384 231 373 223 364 214 357 2",
    "02 354 195 353 185 356 168 365 148;373 172 371 181 372 190 374 197 378 203 384 209 393 215",
    " 405 221 419 227 423 219 424 215 423 213 393 193 382 183 373 172"),
  "Size 5" = paste0(
    "300 112 278 125 264 132 250 137 236 140 226 141 212 136 196 132 159 125 163 133 163 141 13",
    "6 138 109 138 119 142 124 145 98 151 73 159 57 169 44 182 34 200 29 221 20 247 26 243 26 2",
    "65 27 293 21 328 28 316 28 343 49 286 61 362 59 373 53 391 60 390 60 401 50 432 60 427 60 ",
    "450 74 457 86 459 97 458 109 455 121 450 128 350 151 308 156 294 174 297 168 314 158 332 1",
    "92 328 207 324 222 319 218 298 236 296 265 452 277 456 293 458 308 457 325 454 322 399 329",
    " 402 321 386 321 383 329 389 319 346 328 348 317 313 316 292 353 279 342 265 335 250 334 2",
    "34 337 220 343 214 351 210 362 208 372 210 375 213 392 235 406 249 420 259 431 263 436 262",
    " 441 260 451 250 454 241 455 228 436 171 450 174 437 160 431 151 420 123 414 135 387 125 3",
    "64 136 344 120 313 138;382 171 409 186 405 190 401 192 395 192 391 191 387 188 384 184 383",
    " 177 382 171;489 187 483 219 474 246 463 264 457 270 451 275 444 279 435 280 428 280 420 2",
    "78 409 274 402 269 384 253 363 227 362 227 356 228 353 231 352 238 353 247 358 257 368 269",
    " 380 280 395 291 410 297 427 299 440 297 454 291 465 283 474 272 483 257 487 246 492 227 4",
    "92 214 492 201 489 187;450 313 430 317 410 315 408 331 405 342 401 352 396 358 391 361 380",
    " 364 364 365 347 365 345 366 343 371 342 384 344 397 348 404 350 407 357 409 370 405 394 3",
    "96 416 384 428 374 438 357 441 367 442 347 444 338 451 347;188 345 157 350 146 353 142 356",
    " 140 451 151 455 164 457 177 455 194 451")
)
diet_silhouette <- function(g) diet_silhouette_paths[[g]]


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
write_out <- function(df, name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, paste0(name, ".csv"))
  ok <- tryCatch({ write.csv(df, path, row.names = FALSE, na = "", fileEncoding = "UTF-8"); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) stop("Could not write ", path, "\n  Close it (e.g. in Excel) and run again.")
  cat(sprintf("  %-28s %6d rows -> %s\n", name, nrow(df), path))
}
save_fig <- function(p, name, width, height) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, paste0(name, ".png")), p, width = width, height = height, dpi = fig_dpi, bg = page_fill)
  ggsave(file.path(output_dir, paste0(name, ".pdf")), p, width = width, height = height, bg = page_fill,
         device = if (capabilities("cairo")) cairo_pdf else pdf)
  cat(sprintf("  %s.png (+ .pdf)\n", file.path(output_dir, name)))
}
mode_of <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA_character_ else names(sort(table(x), decreasing = TRUE))[1] }

# Diet group from the Smith table's recoded diet and PBDB diet.
diet_group_of <- function(recoded, pbdb) {
  r <- tolower(coalesce(recoded, "")); p <- tolower(coalesce(pbdb, ""))
  graze <- grepl("grazer", p); brow <- grepl("browser|frugivore|folivore", p)
  case_when(r == "carnivore" ~ "Carnivore",
            r == "insectivore" ~ insectivore_as,
            r == "omnivore" ~ "Omnivore",
            r == "herbivore" & graze & brow ~ "Mixed feeder",
            r == "herbivore" & graze ~ "Grazer",
            r == "herbivore" & brow ~ "Browser",
            TRUE ~ NA_character_)
}

# ---- Rowan et al. equations (ported from Rowan_et_al_SI_Code_14MAY.R) --------
# x: region x species 0/1 matrix (regions = rows).
min_b <- function(x) {                       # numerator of Eq. 1 within the sum
  b_ij <- x[1, ] == 1 & x[2, ] == 0
  b_ji <- x[1, ] == 0 & x[2, ] == 1
  if (sum(b_ij) < sum(b_ji)) sp <- colnames(x)[b_ij] else sp <- colnames(x)[b_ji]   # tie: b_ji, as in their code
  list(min_b = min(sum(b_ij), sum(b_ji)), species = sp)
}
endemic_n <- function(x) {                   # sum over pairs of min(b_ij, b_ji), and those species
  cmb <- combn(nrow(x), 2)
  res <- lapply(seq_len(ncol(cmb)), function(i) min_b(x[cmb[, i], , drop = FALSE]))
  list(sum_b = sum(vapply(res, `[[`, numeric(1), "min_b")), species = unlist(lapply(res, `[[`, "species")))
}
shared_n <- function(x) sum(x) - sum(colSums(x) > 0)   # sum_i S_i - S_T
a_sum <- function(x) {                       # sum over pairs of shared species a_ij
  cmb <- combn(nrow(x), 2)
  sum(vapply(seq_len(ncol(cmb)), function(i) sum(colSums(x[cmb[, i], , drop = FALSE]) == 2), numeric(1)))
}
beta_sim  <- function(x) { e <- endemic_n(x)$sum_b; e / (shared_n(x) + e) }          # Eq. 1
beta_end  <- function(x) { e <- endemic_n(x)$sum_b; e / (e + a_sum(x)) }             # Eq. 2
beta_sh   <- function(x) shared_n(x) / sum(x)                                        # Eq. 3
simp_func <- function(x, grp) {                                                      # Eq. 4
  e <- endemic_n(x); f <- factor(grp[match(e$species, colnames(x))], levels = levels(grp))
  table(f) / (shared_n(x) + e$sum_b)
}
endemic_func <- function(x, grp) {                                                   # Eq. 5
  e <- endemic_n(x); f <- factor(grp[match(e$species, colnames(x))], levels = levels(grp))
  table(f) / (e$sum_b + a_sum(x))
}
shared_func <- function(x, grp) {                                                    # Eq. 6
  s_if <- vapply(seq_len(nrow(x)), function(i) as.numeric(table(grp[x[i, ] == 1])), numeric(nlevels(grp)))
  s_if <- if (is.matrix(s_if)) rowSums(s_if) else sum(s_if)
  s_tf <- as.numeric(table(grp[colSums(x) > 0]))
  setNames((s_if - s_tf) / sum(x), levels(grp))
}

# Silhouette grob, drawn in the top-right corner of a panel.
icon_grob <- function(path_string, colour, height = 0.42, outline = NA) {
  parts <- strsplit(path_string, ";", fixed = TRUE)[[1]]
  xy <- do.call(rbind, lapply(seq_along(parts), function(i) {
    v <- as.numeric(strsplit(trimws(parts[i]), " +")[[1]])
    data.frame(x = v[c(TRUE, FALSE)], y = v[c(FALSE, TRUE)], id = i)
  }))
  w <- diff(range(xy$x)); h <- diff(range(xy$y))
  vp <- viewport(x = unit(1, "npc") - unit(6, "pt"), y = unit(1, "npc") - unit(4, "pt"),
                 width = unit(height * w / h, "in"), height = unit(height, "in"), just = c("right", "top"))
  pathGrob((xy$x - min(xy$x)) / w, 1 - (xy$y - min(xy$y)) / h, id = xy$id, rule = "evenodd",
           gp = gpar(fill = colour, col = outline, lwd = 0.6), vp = vp)
}
# PhyloPic silhouettes, straight from the PhyloPic API (no extra package beyond jsonlite/png).
phylopic_dir <- file.path(output_dir, "phylopic")
phylopic_api <- "https://api.phylopic.org"
phylopic_credits <- list()
phylopic_json <- function(path) {
  tmp <- tempfile(fileext = ".json")
  suppressWarnings(download.file(paste0(phylopic_api, path), tmp, quiet = TRUE, mode = "wb",
                                 headers = c(Accept = "application/vnd.phylopic.v2+json")))
  jsonlite::fromJSON(tmp, simplifyVector = FALSE)
}
phylopic_build <- local({ b <- NULL; function() { if (is.null(b)) b <<- phylopic_json("/")$build; b } })
phylopic_image_for_name <- function(name) {        # image record of a taxon's primary silhouette, or NULL
  q <- utils::URLencode(tolower(trimws(name)), reserved = TRUE)
  res <- phylopic_json(sprintf("/nodes?build=%s&filter_name=%s&page=0", phylopic_build(), q))
  items <- res$`_links`$items
  if (!length(items)) return(NULL)
  node_uuid <- sub("^/nodes/([^?]+).*$", "\\1", items[[1]]$href)
  node <- phylopic_json(sprintf("/nodes/%s?build=%s&embed_primaryImage=true", node_uuid, phylopic_build()))
  node$`_embedded`$primaryImage
}
phylopic_image_by_uuid <- function(uuid) phylopic_json(sprintf("/images/%s?build=%s", uuid, phylopic_build()))
phylopic_fetch <- function(taxon) {               # path of a cached PNG for the taxon, or NULL
  dir.create(phylopic_dir, recursive = TRUE, showWarnings = FALSE)
  key <- gsub("[^A-Za-z0-9]+", "_", taxon)
  png_path <- file.path(phylopic_dir, paste0(key, ".png")); meta_path <- file.path(phylopic_dir, paste0(key, ".csv"))
  if (!phylopic_refresh && file.exists(png_path) && file.exists(meta_path)) {
    phylopic_credits[[taxon]] <<- read.csv(meta_path, stringsAsFactors = FALSE)
    return(png_path)
  }
  img <- if (taxon %in% names(phylopic_uuid)) phylopic_image_by_uuid(phylopic_uuid[[taxon]]) else phylopic_image_for_name(taxon)
  if (is.null(img) && grepl(" ", taxon)) img <- phylopic_image_for_name(sub(" .*$", "", taxon))   # genus
  if (is.null(img)) stop("not found on PhyloPic")
  rasters <- img$`_links`$rasterFiles
  if (!length(rasters)) stop("no raster file")
  widths <- vapply(rasters, function(r) as.numeric(sub("x.*$", "", r$sizes)), numeric(1))
  pick <- rasters[[which.min(abs(widths - 512))]]
  suppressWarnings(download.file(pick$href, png_path, quiet = TRUE, mode = "wb"))
  meta <- data.frame(
    Taxon = taxon, Image_UUID = img$uuid,
    Contributor = if (!is.null(img$attribution)) img$attribution
                  else if (!is.null(img$`_links`$contributor$title)) img$`_links`$contributor$title else NA_character_,
    License = if (!is.null(img$`_links`$license$href)) img$`_links`$license$href else NA_character_,
    URL = paste0("https://www.phylopic.org/images/", img$uuid), stringsAsFactors = FALSE)
  write.csv(meta, meta_path, row.names = FALSE)
  phylopic_credits[[taxon]] <<- meta
  png_path
}
# Silhouette grob from PhyloPic (same corner and height as icon_grob); NULL if unavailable.
phylopic_grob <- function(taxon, colour, height = 0.42, outline = NA) {
  img <- tryCatch(png::readPNG(phylopic_fetch(taxon)),
                  error = function(e) { message("  PhyloPic: no silhouette for ", taxon, " (", conditionMessage(e), ")"); NULL })
  if (is.null(img)) return(NULL)
  if (length(dim(img)) == 2) img <- array(rep(img, 2), c(dim(img), 2))       # grey -> grey + alpha
  alpha <- switch(as.character(dim(img)[3]), "2" = img[, , 2], "4" = img[, , 4], 1 - img[, , 1])
  tint <- function(col) {                     # silhouette in one colour, keeping its transparency
    out <- array(0, c(dim(img)[1:2], 4)); rgb <- col2rgb(col) / 255
    for (k in 1:3) out[, , k] <- rgb[k]
    out[, , 4] <- alpha
    out
  }
  place <- function(arr, dx = 0, dy = 0)
    rasterGrob(arr, x = unit(1, "npc") - unit(6 - dx, "pt"), y = unit(1, "npc") - unit(4 - dy, "pt"),
               height = unit(height, "in"), just = c("right", "top"), interpolate = TRUE)
  grobs <- list()
  if (!is.na(outline)) {                      # thin outline for pale fills
    o <- tint(outline)
    grobs <- lapply(list(c(-0.7, 0), c(0.7, 0), c(0, -0.7), c(0, 0.7)), function(d) place(o, d[1], d[2]))
  }
  do.call(grobTree, c(grobs, list(place(tint(colour)))))
}
make_icon <- function(taxon, path_string, colour, outline = NA) {
  g <- if (silhouette_source == "phylopic") phylopic_grob(taxon, colour, outline = outline) else NULL
  if (is.null(g)) icon_grob(path_string, colour, outline = outline) else g
}
if (silhouette_source == "phylopic") {
  miss <- c("jsonlite", "png")[!vapply(c("jsonlite", "png"), requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    message("  For PhyloPic silhouettes install: ", paste(miss, collapse = ", "), "; using the built-in silhouettes.")
    silhouette_source <- "builtin"
  }
}

# Add one icon per facet column to the panels of a given facet row.
add_icons <- function(p, icons, row = 1, height = 0.42) {
  g <- ggplotGrob(p)
  # panels located by position (works for facet_wrap and facet_grid alike)
  pan <- g$layout[grepl("^panel", g$layout$name), ]
  pan <- pan[pan$t == sort(unique(pan$t))[row], ]
  pan <- pan[order(pan$l), ]
  for (k in seq_along(icons)) {
    if (k > nrow(pan) || is.null(icons[[k]])) next
    idx <- which(g$layout$name == pan$name[k])
    g <- gtable::gtable_add_grob(g, icons[[k]], t = g$layout$t[idx], l = g$layout$l[idx],
                                 b = g$layout$b[idx], r = g$layout$r[idx], z = Inf, clip = "off",
                                 name = paste0("icon-", row, "-", k))
  }
  g
}

# =============================================================================
# 1. READ INPUTS
# =============================================================================

cat("=== 1. READING INPUTS ===\n")
md <- read_text_csv(master_file)
time_bins <- read_text_csv(bins_file) %>%
  transmute(Time_Bin = as.integer(Bin_Number), Older = as.numeric(Older_Ma), Younger = as.numeric(Younger_Ma)) %>%
  arrange(Time_Bin) %>% mutate(Mid = (Older + Younger) / 2)

if (!file.exists(smith_file)) {
  alt <- list.files(downloads_dir, pattern = "aao5987.*fixed.*\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  if (length(alt)) smith_file <- alt[1] else stop("Smith et al. workbook not found:\n  ", smith_file)
}
smith <- readxl::read_excel(smith_file, sheet = smith_sheet)
cat(sprintf("  %-50s %7d rows\n", basename(smith_file), nrow(smith)))
names(smith)[grepl("^Recoded diet", names(smith))] <- "Recoded_diet"
smith <- smith %>% transmute(Order = trimws(Order), Family = trimws(Family),
                             Species = trimws(gsub("\\s+", " ", Genus_species)),
                             Smith_Recoded_Diet = Recoded_diet, Smith_PBDB_Diet = `PBDB diet`,
                             key = tolower(Species), Genus = sub(" .*$", "", Species),
                             ln_mass = suppressWarnings(as.numeric(`ln Mass (g)`)),
                             Diet_Group = diet_group_of(Recoded_diet, `PBDB diet`)) %>%
  filter(!is.na(Species))
smith_sp <- smith %>% filter(!duplicated(key))
smith_gen <- smith %>% group_by(gkey = tolower(Genus)) %>%
  summarise(Order_g = mode_of(Order), Family_g = mode_of(Family), ln_mass_g = mean(ln_mass, na.rm = TRUE),
            Diet_g = mode_of(Diet_Group), .groups = "drop") %>%
  mutate(ln_mass_g = ifelse(is.nan(ln_mass_g), NA, ln_mass_g))

# Order of every genus in our own records (FAUNMAP 'Order', PBDB 'order').
genus_order <- bind_rows(
  tryCatch(read_text_csv(fauna_file) %>% transmute(Genus, Order), error = function(e) NULL),
  tryCatch(read_text_csv(pbdb_file) %>% transmute(Genus = genus, Order = order), error = function(e) NULL)) %>%
  filter(!is.na(Genus), !is.na(Order)) %>%
  group_by(gkey = tolower(trimws(Genus))) %>% summarise(Order_data = mode_of(trimws(Order)), .groups = "drop")

md <- md %>% filter(Spatial_Bin %in% regions) %>% mutate(Time_Bin = as.integer(Time_Bin))

# =============================================================================
# 2. SPECIES TRAITS AND THE LARGE-MAMMAL FILTER
# =============================================================================

cat("\n=== 2. SPECIES TRAITS ===\n")
traits <- md %>% distinct(Species) %>%
  mutate(key = tolower(Species), gkey = tolower(sub(" .*$", "", Species))) %>%
  left_join(smith_sp %>% select(key, Order_s = Order, Family, ln_mass_s = ln_mass, Diet_s = Diet_Group,
                                Smith_Recoded_Diet, Smith_PBDB_Diet), by = "key") %>%
  left_join(smith_gen, by = "gkey") %>%
  left_join(genus_order, by = "gkey") %>%
  mutate(Order = coalesce(Order_s, Order_data, Order_g),
         Family = coalesce(Family, Family_g),
         ln_mass = coalesce(ln_mass_s, if (genus_fallback) ln_mass_g else NA_real_),
         Diet_Group = coalesce(Diet_s, if (genus_fallback) Diet_g else NA_character_),
         Mass_Source = case_when(!is.na(ln_mass_s) ~ "species (Smith et al.)",
                                 genus_fallback & !is.na(ln_mass_g) ~ "genus mean (Smith et al.)", TRUE ~ "none"),
         Diet_Source = case_when(!is.na(Diet_s) ~ "species (Smith et al.)",
                                 genus_fallback & !is.na(Diet_g) ~ "genus majority (Smith et al.)", TRUE ~ "none"),
         Mass_kg = exp(ln_mass) / 1000)

if (file.exists(overrides_file)) {
  ov <- read.csv(overrides_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  cat(sprintf("  Hand assignments from %s: %d species\n", basename(overrides_file), nrow(ov)))
  m <- match(tolower(traits$Species), tolower(trimws(ov$Species)))
  has <- !is.na(m)
  if ("Order" %in% names(ov)) traits$Order[has] <- coalesce(ov$Order[m[has]], traits$Order[has])
  if ("Mass_kg" %in% names(ov)) {
    new_mass <- suppressWarnings(as.numeric(ov$Mass_kg[m[has]]))
    traits$Mass_Source[has][!is.na(new_mass)] <- "hand assignment"
    traits$Mass_kg[has] <- coalesce(new_mass, traits$Mass_kg[has])
  }
  if ("Diet_Group" %in% names(ov)) {
    new_diet <- ov$Diet_Group[m[has]]
    traits$Diet_Source[has][!is.na(new_diet)] <- "hand assignment"
    traits$Diet_Group[has] <- coalesce(new_diet, traits$Diet_Group[has])
  }
}

traits <- traits %>%
  mutate(Size_Class = as.character(cut(Mass_kg, size_breaks, labels = size_groups, right = FALSE)),
         Large_Order = Order %in% large_mammal_orders,
         Included = Large_Order & (is.na(Mass_kg) | Mass_kg > min_mass_kg),
         Exclusion = case_when(Included ~ NA_character_,
                               is.na(Order) ~ "order unknown",
                               !Large_Order ~ paste("order", Order, "not a large-mammal order"),
                               TRUE ~ sprintf("body mass %.2f kg <= %g kg", Mass_kg, min_mass_kg))) %>%
  select(Species, Order, Family, Mass_kg, Size_Class, Diet_Group, Mass_Source, Diet_Source,
         Smith_Recoded_Diet, Smith_PBDB_Diet, Included, Exclusion) %>%
  arrange(desc(Included), Order, Species)

inc <- traits %>% filter(Included)
cat(sprintf("  %d species in the three regions | %d large mammals kept\n", nrow(traits), nrow(inc)))
print(as.data.frame(count(traits %>% filter(!Included), Exclusion, sort = TRUE)), row.names = FALSE)
cat(sprintf("  Large mammals with a diet group: %d | with a size class: %d\n",
            sum(!is.na(inc$Diet_Group)), sum(!is.na(inc$Size_Class))))
unclassified <- inc %>% filter(is.na(Diet_Group) | is.na(Size_Class)) %>%
  mutate(Missing = paste0(ifelse(is.na(Diet_Group), "diet ", ""), ifelse(is.na(Size_Class), "size", "")))

# =============================================================================
# 3. REGION x SPECIES MATRICES PER TIME BIN AND ROWAN'S METRICS
# =============================================================================

cat("\n=== 3. BETA DIVERSITY (Rowan et al. Eqs. 1-6) ===\n")
occ <- md %>% filter(Species %in% inc$Species) %>% distinct(Time_Bin, Spatial_Bin, Species)
diet_f <- factor(inc$Diet_Group, levels = diet_groups); names(diet_f) <- inc$Species
size_f <- factor(inc$Size_Class, levels = size_groups); names(size_f) <- inc$Species

pa_list <- list(); overall <- list(); parts <- list()
for (b in time_bins$Time_Bin) {
  d <- occ %>% filter(Time_Bin == b)
  regs <- intersect(regions, unique(d$Spatial_Bin))
  if (length(regs) < 2) {
    cat(sprintf("  Bin %d: fewer than two regions with large mammals - skipped\n", b)); next
  }
  spp <- sort(unique(d$Species))
  x <- matrix(0L, length(regs), length(spp), dimnames = list(regs, spp))
  x[cbind(match(d$Spatial_Bin, regs), match(d$Species, spp))] <- 1L
  pa_list[[as.character(b)]] <- x
  overall[[length(overall) + 1]] <- data.frame(
    Time_Bin = b, n_regions = length(regs), regions = paste(regs, collapse = "; "),
    n_species = ncol(x), species_per_region = paste(rowSums(x), collapse = "/"),
    beta_SIM = beta_sim(x), beta_SIM_END = beta_end(x), beta_SIM_SH = beta_sh(x))
  for (grp in c("Diet", "Size")) {
    f <- if (grp == "Diet") diet_f[spp] else size_f[spp]
    parts[[length(parts) + 1]] <- data.frame(
      Time_Bin = b, Group_Type = grp, Group = levels(f),
      beta_SIM_f = as.numeric(simp_func(x, f)),
      beta_SIM_END_f = as.numeric(endemic_func(x, f)),
      beta_SIM_SH_f = as.numeric(shared_func(x, f)),
      n_species = as.numeric(table(f)))
  }
}
two_region_bins <- vapply(pa_list, nrow, integer(1))
if (any(two_region_bins < length(regions))) {
  cat(sprintf("  NOTE: bin(s) %s have large mammals in only two of the three regions; beta_SIM uses those two.\n",
              paste(names(two_region_bins)[two_region_bins < length(regions)], collapse = ", ")))
}
beta_by_bin <- bind_rows(overall) %>% left_join(time_bins, by = "Time_Bin") %>% relocate(Older, Younger, Mid, .after = Time_Bin)
partition <- bind_rows(parts) %>% left_join(time_bins %>% select(Time_Bin, Mid), by = "Time_Bin")

# Check: the groups add up to the overall values (exactly, when every species has a group).
chk <- partition %>% group_by(Time_Bin, Group_Type) %>%
  summarise(sum_f = sum(beta_SIM_f), .groups = "drop") %>%
  left_join(beta_by_bin %>% select(Time_Bin, beta_SIM), by = "Time_Bin") %>%
  mutate(unassigned_share = beta_SIM - sum_f)
cat("  beta_SIM per bin:\n")
print(as.data.frame(beta_by_bin %>% transmute(Time_Bin, Mid, n_species, species_per_region,
                                              beta_SIM = round(beta_SIM, 3), beta_SIM_END = round(beta_SIM_END, 3),
                                              beta_SIM_SH = round(beta_SIM_SH, 3))), row.names = FALSE)
if (any(abs(chk$unassigned_share) > 1e-9)) {
  cat("  NOTE: part of beta_SIM comes from species without a diet group / size class (see unclassified_species.csv):\n")
  print(as.data.frame(chk %>% filter(abs(unassigned_share) > 1e-9) %>%
                        mutate(across(c(sum_f, beta_SIM, unassigned_share), ~ round(.x, 3)))), row.names = FALSE)
}

# OLS of beta_SIM on bin midpoint, as reported by Rowan et al.
fit <- lm(beta_SIM ~ Mid, data = beta_by_bin)
regression <- data.frame(n_bins = nrow(beta_by_bin), slope_per_Myr = coef(fit)[["Mid"]],
                         intercept = coef(fit)[[1]], r_squared = summary(fit)$r.squared,
                         p_value = if (nrow(beta_by_bin) > 2) summary(fit)$coefficients["Mid", 4] else NA)
cat(sprintf("  OLS beta_SIM ~ midpoint age: r2 = %.2f, P = %.3g (slope %.3f per Myr; positive = higher beta in the past)\n",
            regression$r_squared, regression$p_value, regression$slope_per_Myr))

# =============================================================================
# 4. FIGURES
# =============================================================================

# Equal-species beta_SIM (the curve of Figs. 3 and 4): in each bin every region
# is cut to the species count of the poorest region by random draws and the
# mean beta_SIM of 'equal_draws' draws is taken (as in the Step 7 addendum).
set.seed(equal_seed)
equal_beta <- bind_rows(lapply(names(pa_list), function(b) {
  x <- pa_list[[b]]
  sp_lists <- lapply(seq_len(nrow(x)), function(i) colnames(x)[x[i, ] == 1])
  n_min <- min(lengths(sp_lists))
  draws <- replicate(equal_draws, {
    picks <- lapply(sp_lists, function(s) if (length(s) > n_min) sample(s, n_min) else s)
    spp <- sort(unique(unlist(picks)))
    m <- t(vapply(picks, function(s) as.integer(spp %in% s), integer(length(spp))))
    colnames(m) <- spp
    beta_sim(m)
  })
  data.frame(Time_Bin = as.integer(b), n_species_equal = n_min, beta_SIM_equal_species = mean(draws))
}))
beta_by_bin <- beta_by_bin %>% left_join(equal_beta, by = "Time_Bin")
cat("  Equal-species beta_SIM (curve of Figs. 3-4):\n")
print(as.data.frame(beta_by_bin %>% transmute(Time_Bin, n_species_equal,
                                              beta_SIM_equal_species = round(beta_SIM_equal_species, 3))),
      row.names = FALSE)

cat("\n=== 4. FIGURES ===\n")
x_max <- max(time_bins$Older); x_min <- 0
bands <- time_bins %>% filter(Time_Bin %% 2 == 1)          # shade alternate bins, as in the paper
x_breaks <- seq(floor(x_max), 0, by = -1)
theme_rowan <- function() {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(colour = ink, linewidth = 0.4),
          axis.ticks = element_line(colour = ink, linewidth = 0.4),
          axis.text = element_text(colour = ink),
          axis.title = element_text(colour = ink),
          strip.background = element_rect(fill = strip_fill, colour = ink_soft, linewidth = 0.4),
          strip.text = element_text(colour = ink, size = base_size * 0.9),
          panel.spacing.x = unit(0.9, "lines"),
          plot.tag = element_text(face = "bold", size = base_size * 1.3, colour = ink),
          panel.background = element_rect(fill = page_fill, colour = NA),
          plot.background = element_rect(fill = page_fill, colour = NA))
}
band_layer <- geom_rect(data = bands, aes(xmin = Younger, xmax = Older, ymin = -Inf, ymax = Inf),
                        inherit.aes = FALSE, fill = band_fill, colour = NA)
x_scale <- scale_x_reverse(breaks = x_breaks)   # range set by coord_cartesian (keeps annotations outside it)
beta_lab <- expression(beta[SIM])

# ---- Fig. 2 ------------------------------------------------------------------
p2 <- ggplot(beta_by_bin, aes(Mid, beta_SIM)) +
  band_layer +
  geom_line(colour = ink, linewidth = 0.6) +
  geom_point(shape = 21, fill = point_fill, colour = ink, size = 3.6, stroke = 0.8) +
  scale_x_reverse(breaks = unique(c(x_max, x_breaks)),                # older limit of Bin 1 (3.25) marked too
                  labels = function(v) sub("\\.?0+$", "", sprintf("%.2f", v))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = function(v) sprintf("%.2f", v),
                     expand = expansion(mult = c(0.02, 0.04))) +
  # right-hand guide, as in the paper: up = provincialism, down = homogenization
  annotate("segment", x = -0.07 * x_max, xend = -0.07 * x_max, y = 0.53, yend = 0.98, colour = ink, linewidth = 0.5,
           arrow = arrow(length = unit(0.1, "in"), type = "open")) +
  annotate("segment", x = -0.07 * x_max, xend = -0.07 * x_max, y = 0.47, yend = 0.02, colour = ink, linewidth = 0.5,
           arrow = arrow(length = unit(0.1, "in"), type = "open")) +
  annotate("text", x = -0.10 * x_max, y = 0.755, label = "Higher \u03b2 (provincialism)", angle = 90,
           size = base_size / 4.1, colour = ink) +
  annotate("text", x = -0.10 * x_max, y = 0.245, label = "Lower \u03b2 (homogenization)", angle = 90,
           size = base_size / 4.1, colour = ink) +
  coord_cartesian(xlim = c(x_max, x_min), clip = "off", expand = TRUE) +
  labs(x = "Age (Ma)", y = beta_lab) +
  theme_rowan() +
  theme(plot.margin = margin(10, 40, 8, 10), axis.title.y = element_text(size = base_size * 1.15))
save_fig(p2, "Fig2_beta_sim", fig2_w, fig2_h)

# ---- Fig. 3 ------------------------------------------------------------------
# Two versions of Figs. 3 and 4: curve = equal-species beta_SIM (main files) or
# the full-data beta_SIM of Fig. 2 (files ending in _full_data).
overall_pts <- beta_by_bin %>% select(Mid, beta_SIM = beta_SIM_equal_species)   # equal-species curve
full_pts    <- beta_by_bin %>% select(Mid, beta_SIM)                            # all species (as Fig. 2)
fig3_panel <- function(type, groups, colours, tag, curve = overall_pts) {
  d <- partition %>% filter(Group_Type == type) %>% mutate(Group = factor(Group, levels = groups))
  ov <- tidyr::crossing(curve, Group = factor(groups, levels = groups))
  bw <- time_bins %>% transmute(Mid, w = (Older - Younger) * 0.8)
  d <- d %>% left_join(bw, by = "Mid") %>% filter(beta_SIM_f > 0)   # no bar where the group adds nothing
  ggplot() +
    band_layer +
    geom_rect(data = d, aes(xmin = Mid - w / 2, xmax = Mid + w / 2, ymin = 0, ymax = beta_SIM_f, fill = Group),
              colour = ink, linewidth = 0.25) +
    geom_line(data = ov, aes(Mid, beta_SIM), colour = ink, linewidth = 0.7) +   # observed values, no smoothing
    geom_point(data = ov, aes(Mid, beta_SIM), shape = 21, fill = point_fill, colour = ink, size = 2.3, stroke = 0.6) +
    facet_wrap(~ Group, nrow = 1, drop = FALSE) +
    scale_fill_manual(values = colours, guide = "none", drop = FALSE) +
    x_scale +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(v) sprintf("%.2f", v)) +
    coord_cartesian(xlim = c(x_max + 0.2, x_min - 0.2), ylim = c(0, 1.08), expand = FALSE) +
    labs(x = "Age (Ma)", y = beta_lab, tag = if (panel_letters) tag else NULL) +
    theme_rowan()
}
diet_icons <- lapply(diet_groups, function(g) make_icon(phylopic_diet[[g]], diet_silhouette(g), diet_colours[[g]]))
size_icons <- lapply(seq_along(size_groups), function(i)
  make_icon(phylopic_size[[i]], size_silhouette_paths[[i]], size_colours[[i]], outline = if (i == 1) ink_soft else NA))
if (length(phylopic_credits)) {
  write.csv(do.call(rbind, phylopic_credits), file.path(output_dir, "phylopic_credits.csv"), row.names = FALSE)
  cat("  PhyloPic silhouettes used; credits in phylopic_credits.csv (credit the contributors on your slide)\n")
}
g3a <- add_icons(fig3_panel("Diet", diet_groups, diet_colours, "a") + labs(x = NULL), diet_icons)
g3b <- add_icons(fig3_panel("Size", size_groups, size_colours, "b"), size_icons)
p3 <- patchwork::wrap_plots(patchwork::wrap_elements(g3a), patchwork::wrap_elements(g3b), ncol = 1)
save_fig(p3, "Fig3_partitioned_beta_sim", fig3_w, fig3_h)
g3a_full <- add_icons(fig3_panel("Diet", diet_groups, diet_colours, "a", full_pts) + labs(x = NULL), diet_icons)
g3b_full <- add_icons(fig3_panel("Size", size_groups, size_colours, "b", full_pts), size_icons)
save_fig(patchwork::wrap_plots(patchwork::wrap_elements(g3a_full), patchwork::wrap_elements(g3b_full), ncol = 1),
         "Fig3_partitioned_beta_sim_full_data", fig3_w, fig3_h)

# ---- Fig. 4 ------------------------------------------------------------------
fig4_panel <- function(type, groups, colours, tag, curve = overall_pts) {
  d <- partition %>% filter(Group_Type == type) %>%
    select(Mid, Group, Shared = beta_SIM_SH_f, Endemic = beta_SIM_END_f) %>%
    pivot_longer(c(Shared, Endemic), names_to = "Component", values_to = "Proportion") %>%
    mutate(Group = factor(Group, levels = groups), Component = factor(Component, levels = c("Shared", "Endemic"))) %>%
    left_join(time_bins %>% transmute(Mid, w = (Older - Younger) * 0.8), by = "Mid")
  y_top <- 1.2                                                       # room for beta_SIM (0-1) and the row labels
  d <- d %>% filter(Proportion > 0)                                  # no bar where the proportion is 0
  ov <- tidyr::crossing(curve, Group = factor(groups, levels = groups),
                        Component = factor(c("Shared", "Endemic"), levels = c("Shared", "Endemic")))
  lab <- data.frame(Group = factor(groups[1], levels = groups),
                    Component = factor(c("Shared", "Endemic"), levels = c("Shared", "Endemic")),
                    label = c("Shared species", "Endemic species"))
  ggplot(d) +
    band_layer +
    geom_rect(aes(xmin = Mid - w / 2, xmax = Mid + w / 2, ymin = 0, ymax = Proportion, fill = Group),
              colour = ink, linewidth = 0.25) +
    geom_line(data = ov, aes(Mid, beta_SIM), colour = ink, linewidth = 0.7) +   # equal-species beta_SIM, as Fig. 3
    geom_point(data = ov, aes(Mid, beta_SIM), shape = 21, fill = point_fill, colour = ink, size = 2.3, stroke = 0.6) +
    geom_text(data = lab, aes(x = x_max - 0.1, y = y_top * 0.98, label = label), hjust = 0, vjust = 1,
              size = base_size / 3.3, colour = ink) +
    facet_grid(Component ~ Group, drop = FALSE) +
    scale_fill_manual(values = colours, guide = "none", drop = FALSE) +
    x_scale +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(v) sprintf("%.2f", v)) +
    coord_cartesian(xlim = c(x_max + 0.2, x_min - 0.2), ylim = c(0, y_top), expand = FALSE) +
    labs(x = "Age (Ma)", y = expression("Proportion;" ~ beta[SIM] ~ "(points)"), tag = if (panel_letters) tag else NULL) +
    theme_rowan() +
    theme(strip.text.y = element_blank(), strip.background.y = element_blank(),
          panel.spacing.y = unit(0.6, "lines"))
}
g4a <- add_icons(fig4_panel("Diet", diet_groups, diet_colours, "a"), diet_icons, row = 1)
g4b <- add_icons(fig4_panel("Size", size_groups, size_colours, "b"), size_icons, row = 1)
p4 <- patchwork::wrap_plots(patchwork::wrap_elements(g4a), patchwork::wrap_elements(g4b), ncol = 1)
save_fig(p4, "Fig4_shared_endemic", fig4_w, fig4_h)
g4a_full <- add_icons(fig4_panel("Diet", diet_groups, diet_colours, "a", full_pts), diet_icons, row = 1)
g4b_full <- add_icons(fig4_panel("Size", size_groups, size_colours, "b", full_pts), size_icons, row = 1)
save_fig(patchwork::wrap_plots(patchwork::wrap_elements(g4a_full), patchwork::wrap_elements(g4b_full), ncol = 1),
         "Fig4_shared_endemic_full_data", fig4_w, fig4_h)

# =============================================================================
# 5. TABLES
# =============================================================================

cat("\n=== 5. SAVING TABLES ===\n")
write_out(beta_by_bin %>% mutate(across(starts_with("beta"), ~ round(.x, 4))), "beta_sim_by_bin")
write_out(partition %>% mutate(across(starts_with("beta"), ~ round(.x, 4))), "partition_by_group")
write_out(regression %>% mutate(across(everything(), ~ signif(.x, 4))), "regression_beta_sim_age")
write_out(traits %>% mutate(Mass_kg = round(Mass_kg, 3)), "species_traits")
write_out(unclassified %>% mutate(Mass_kg = round(Mass_kg, 3)), "unclassified_species")
saveRDS(pa_list, file.path(output_dir, "regional_pa_matrices.rds"))

# ---- Excel workbook: every species with its diet group and body-size class ----
# Where each species occurs (regions and time bins, from the master data).
occurs <- md %>% filter(Species %in% traits$Species) %>%
  group_by(Species) %>%
  summarise(Regions = paste(intersect(regions, unique(Spatial_Bin)), collapse = "; "),
            Time_Bins = paste(sort(unique(Time_Bin)), collapse = ", "),
            n_localities = n_distinct(Site_Key), .groups = "drop")
book_all <- traits %>% left_join(occurs, by = "Species") %>%
  transmute(Species, Order, Family, `Body mass (kg)` = round(Mass_kg, 2), `Size class` = Size_Class,
            `Diet group` = Diet_Group, `In the analyses (large mammal)` = ifelse(Included, "yes", "no"),
            `Reason left out` = Exclusion, `Mass source` = Mass_Source, `Diet source` = Diet_Source,
            `Smith et al. recoded diet` = Smith_Recoded_Diet, `Smith et al. PBDB diet` = Smith_PBDB_Diet,
            Regions, `Time bins` = Time_Bins, `Localities` = n_localities)
book_inc <- book_all %>% filter(`In the analyses (large mammal)` == "yes")
by_diet <- book_inc %>%
  mutate(`Diet group` = factor(coalesce(`Diet group`, "Unclassified"), levels = c(diet_groups, "Unclassified"))) %>%
  arrange(`Diet group`, Species) %>% mutate(`Diet group` = as.character(`Diet group`)) %>%
  select(`Diet group`, Species, Order, Family, `Body mass (kg)`, `Size class`, `Diet source`,
         `Smith et al. recoded diet`, `Smith et al. PBDB diet`, Regions, `Time bins`)
by_size <- book_inc %>%
  mutate(`Size class` = factor(coalesce(`Size class`, "No body mass"), levels = c(size_groups, "No body mass"))) %>%
  arrange(`Size class`, `Body mass (kg)`, Species) %>% mutate(`Size class` = as.character(`Size class`)) %>%
  select(`Size class`, Species, Order, Family, `Body mass (kg)`, `Mass source`, `Diet group`, Regions, `Time bins`)
cross <- book_inc %>%
  mutate(Diet = factor(coalesce(`Diet group`, "Unclassified"), levels = c(diet_groups, "Unclassified")),
         Size = factor(coalesce(`Size class`, "No body mass"), levels = c(size_groups, "No body mass"))) %>%
  count(Diet, Size) %>% tidyr::pivot_wider(names_from = Size, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>% rename(`Diet group` = Diet) %>%
  mutate(`Diet group` = as.character(`Diet group`))
cross <- bind_rows(cross, cross %>% summarise(across(where(is.numeric), sum)) %>% mutate(`Diet group` = "Total"))
notes <- data.frame(Notes = c(
  "Species of the three regions (Basin and Range, Coastal Plain, Great Plains incl. Central Lowland west of the Mississippi).",
  sprintf("'In the analyses' = large mammals as in Rowan et al. 2024: orders %s, body mass > %g kg.",
          paste(large_mammal_orders, collapse = ", "), min_mass_kg),
  "Body mass and diet: Smith et al. 2018 (Science, aao5987), Table S7 (fixed workbook). Mass = exp(ln mass) / 1000.",
  "Size classes (Rowan et al.): Size 1 <18 kg, Size 2 18-80 kg, Size 3 80-350 kg, Size 4 350-1,000 kg, Size 5 >1,000 kg.",
  paste0("Diet groups (Rowan et al.): carnivore (and insectivore) -> Carnivore; omnivore -> Omnivore; herbivore with grazer and ",
         "browser/frugivore/folivore terms -> Mixed feeder; grazer only -> Grazer; browser/frugivore/folivore only -> Browser (no frugivores in these data)."),
  "Source columns: 'species' = the species's own row in Smith et al.; 'genus' = mean mass / most common diet of its genus there; 'hand assignment' = species_traits_overrides.csv.",
  "Unclassified species can be given a diet group or mass in species_traits_overrides.csv (columns Species, Order, Mass_kg, Diet_Group); then run Step 7 again."))
book <- list(`All species` = book_all, `By diet group` = by_diet, `By size class` = by_size,
             `Diet x size (counts)` = cross, `Unclassified` = book_inc %>% filter(is.na(`Diet group`) | is.na(`Size class`)),
             Notes = notes)
book_path <- file.path(output_dir, "species_functional_groups.xlsx")
if (requireNamespace("writexl", quietly = TRUE)) {
  ok <- tryCatch({ writexl::write_xlsx(book, book_path); TRUE }, error = function(e) FALSE)
  if (ok) cat(sprintf("  %-28s %6d species -> %s\n", "species_functional_groups", nrow(book_all), book_path))
  else cat("  NOTE: could not write", book_path, "- close it in Excel and run again.\n")
} else {
  cat("  NOTE: install.packages(\"writexl\") to get species_functional_groups.xlsx (the CSV files have the same data).\n")
}

cat("\n=== STEP 7 COMPLETE ===\n")
cat("Objects in your Environment: beta_by_bin, partition, traits, pa_list\n")
