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
#         (Step 7). Silhouettes come from PhyloPic as in Step 7 (same taxa, same cache in
#         Outputs/7_rowan_figures/phylopic/), or from the built-in ones embedded below.
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
silhouette_source <- "phylopic"            # "phylopic" or "builtin" (the embedded silhouettes at the end of SETTINGS)
# PhyloPic taxa, as in Step 7 (same cache folder, so Step 7's downloads are reused).
phylopic_diet <- c("Carnivore" = "Smilodon fatalis", "Omnivore" = "Ursus americanus",
                   "Browser" = "Odocoileus virginianus", "Mixed feeder" = "Camelops hesternus",
                   "Grazer" = "Bison bison")
phylopic_size <- c("Vulpes vulpes", "Odocoileus virginianus", "Camelops hesternus",
                   "Bison bison", "Mammuthus columbi")         # Size 1 ... Size 5
phylopic_uuid <- c()        # optional, e.g. c("Bison bison" = "<image uuid from phylopic.org/images/...>")
phylopic_refresh <- FALSE   # TRUE = download again even if a silhouette is already cached

# Built-in silhouettes (x y pairs on a 512 x 512 grid; ";" separates parts), as in Step 7.
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

# =============================================================================
# 1. DATA
# =============================================================================

setwd(work_dir)
if (!file.exists(pa_file) || !file.exists(traits_file)) stop("Run Step 7 first (", pa_file, ", ", traits_file, ").")
pa     <- readRDS(pa_file)
traits <- read.csv(traits_file, stringsAsFactors = FALSE, check.names = FALSE)

if (silhouette_source == "phylopic" &&
    !all(vapply(c("jsonlite", "png"), requireNamespace, logical(1), quietly = TRUE))) {
  message("  For PhyloPic silhouettes install jsonlite and png; using the built-in silhouettes.")
  silhouette_source <- "builtin"
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
