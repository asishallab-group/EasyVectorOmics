source("r/tensoromics_functions.R")

#' Comprehensive test function for paralog analysis
test_paralog_functions <- function() {
  cat("=== Testing Mask Logic ===\n")

  n_paralogs <- 5L
  i_paralog <- 3L
  chunk_count <- mask_chunk_count(n_paralogs)
  cat("Chunk count for", n_paralogs, "paralogs:", chunk_count, "\n")

  bit_mask <- integer(chunk_count)
  bit_mask[i_paralog %/% 32 + 1L] <- bitwShiftL(1L, i_paralog %% 32 - 1)  # -1 because 1 obtains one bit already
  state <- mask_check_state(bit_mask, i_paralog)
  cat("Paralog", i_paralog, "active in mask:", state, "\n")

  cat("\n=== Testing Pattern Filtering ===\n")

  angles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  threshold <- 0.6

  dosage_mask <- filter_paralogs_by_pattern_dosage_effect(angles, threshold)
  cat("Dosage effect mask:", dosage_mask, "\n")

  subfunc_mask <- filter_paralogs_by_pattern_subfunctionalization(angles, threshold)
  cat("Subfunctionalization mask:", subfunc_mask, "\n")

  cat("\n=== Testing Work Array Size Calculation ===\n")

  max_subset_size <- 3L
  work_size_info <- calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, dosage_mask)
  cat("Work array size:", work_size_info$work_array_size, "\n")
  cat("Adjusted max subset size:", work_size_info$max_subset_size, "\n")

  cat("\n=== Testing Dosage Effect Detection ===\n")

  ancestor <- c(1.0, 1.0)
  paralogs <- matrix(c(1.1, 0.9, 1.2, 0.8, 1.3, 0.7, 1.4, 0.6, 1.5, 0.5), nrow = 2)
  dosage_result <- detect_dosage_effect(
    ancestor = ancestor,
    paralogs = paralogs,
    filtered_paralogs_mask = dosage_mask,
    max_subset_size = work_size_info$max_subset_size,
    gain_gamma = 0.1,
    max_angle = pi
  )
  cat("Dosage effect results:", dosage_result$n_results, "subsets\n")
  print(dosage_result$work_arr_paralog_subsets)

  cat("\n=== Testing Subfunctionalization Detection ===\n")

  norms <- sqrt(colSums(paralogs^2))
  sorted_perm <- order(norms)
  subfunc_result <- detect_subfunctionalization(
    ancestor = ancestor,
    paralogs = paralogs,
    rdi_threshold = 0.5,
    filtered_paralogs_mask = subfunc_mask,
    max_subset_size = work_size_info$max_subset_size,
    paralog_norms = norms,
    sorted_paralog_norms_perm = sorted_perm
  )
  cat("Subfunctionalization results:", subfunc_result$n_results, "subsets\n")
  print(subfunc_result$work_arr_paralog_subsets)

  cat("\n=== Testing Edge Cases ===\n")

  empty_mask <- tryCatch(mask_chunk_count(0L), error = function(e) e)
  cat("Empty paralog count test:", ifelse(inherits(empty_mask, "error"), "Error", "Success"), "\n")

  single_mask <- mask_chunk_count(1L)
  cat("Single paralog mask chunk count:", single_mask, "\n")
}

test_paralog_functions()
