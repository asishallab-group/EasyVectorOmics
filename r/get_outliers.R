# =====================
# Test cases for compute_family_scaling_r Fortran wrapper
# Each case checks a different combination of optional arguments and error handling

dyn.load("build/libtensor-omics.so")

# Common test data
n_genes <- as.integer(5)
n_families <- as.integer(2)
distances <- as.double(c(1, 2, 3, 4, 5))
gene_to_fam <- as.integer(c(1, 1, 2, 2, 2))
dscale <- double(n_families)
loess_x <- double(n_families)
loess_y <- double(n_families)
indices_used <- integer(n_families)
perm_tmp <- integer(n_genes)
stack_left_tmp <- integer(n_genes)
stack_right_tmp <- integer(n_genes)
workspace_weights <- double(n_families)
workspace_values <- matrix(0, 1, n_families)
error_code <- integer(1)


##
# [compute_family_scaling_r] Case 1: Basic LOESS scaling
# Expectation: dscale should be positive for both families (since both have >1 gene), error_code should be 0.
cat("\n[compute_family_scaling_r] Case 1: Basic LOESS scaling\n")
result1 <- .Fortran("compute_family_scaling_r",
  n_genes = n_genes,
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam,
  dscale = dscale,
  loess_x = loess_x,
  loess_y = loess_y,
  indices_used = indices_used,
  perm_tmp = perm_tmp,
  stack_left_tmp = stack_left_tmp,
  stack_right_tmp = stack_right_tmp,
  workspace_weights = workspace_weights,
  workspace_values = workspace_values,
  error_code = error_code
)
print(result1$dscale)
print(result1$error_code)

##
# [compute_family_scaling_r] Case 2: Invalid family indices
# Expectation: dscale should be all -1, error_code should be -2 (invalid family index in gene_to_fam).
cat("\n[compute_family_scaling_r] Case 2: Invalid family indices\n")
gene_to_fam_invalid <- as.integer(c(1, 3, 2, 2, 2))
result2 <- .Fortran("compute_family_scaling_r",
  n_genes = n_genes,
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam_invalid,
  dscale = dscale,
  loess_x = loess_x,
  loess_y = loess_y,
  indices_used = indices_used,
  perm_tmp = perm_tmp,
  stack_left_tmp = stack_left_tmp,
  stack_right_tmp = stack_right_tmp,
  workspace_weights = workspace_weights,
  workspace_values = workspace_values,
  error_code = error_code
)
print(result2$dscale)
print(result2$error_code)

##
# [compute_family_scaling_r] Case 3: Family with only one gene
# Expectation: dscale for the single-gene family should be 0, error_code should be 0.
cat("\n[compute_family_scaling_r] Case 3: Family with only one gene\n")
gene_to_fam_single <- as.integer(c(1, 2, 2, 2, 2))
result3 <- .Fortran("compute_family_scaling_r",
  n_genes = n_genes,
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam_single,
  dscale = dscale,
  loess_x = loess_x,
  loess_y = loess_y,
  indices_used = indices_used,
  perm_tmp = perm_tmp,
  stack_left_tmp = stack_left_tmp,
  stack_right_tmp = stack_right_tmp,
  workspace_weights = workspace_weights,
  workspace_values = workspace_values,
  error_code = error_code
)
print(result3$dscale)
print(result3$error_code)


# =====================
# Test cases for compute_rdi_r Fortran wrapper
# =====================

##
# [compute_rdi_r] Case 1: normal input
# Expectation: rdi = abs(distances) / dscale; should be c(0.5, 1, 0.75, 1, 1.25) for dscale = c(2,4).
cat("\n[compute_rdi_r] Case 1: normal input\n")
distances <- as.double(c(1, 2, 3, 4, 5))
gene_to_fam <- as.integer(c(1, 1, 2, 2, 2))
dscale <- as.double(c(2, 4))
rdi <- double(length(distances))
res1 <- .Fortran("compute_rdi_r",
  n_genes = as.integer(length(distances)),
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam,
  dscale = as.double(c(2, 4)),
  rdi = rdi
)
print(res1$rdi)

##
# [compute_rdi_r] Case 2: dscale with zeros
# Expectation: All rdi values should be 0 (since scaling is zero for all families).
cat("\n[compute_rdi_r] Case 2: dscale with zeros\n")
dscale <- as.double(c(0, 0))
rdi <- double(length(distances))
res2 <- .Fortran("compute_rdi_r",
  n_genes = as.integer(length(distances)),
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam,
  dscale = dscale,
  rdi = rdi
)
print(res2$rdi)

##
# [compute_rdi_r] Case 3: gene_to_fam out of range
# Expectation: rdi for genes with invalid family index should be -1, others normal.
cat("\n[compute_rdi_r] Case 3: gene_to_fam out of range\n")
gene_to_fam_bad <- as.integer(c(1, 3, 2, 2, 2)) # 3 doesn't exist in dscale
rdi <- double(length(distances))
res3 <- .Fortran("compute_rdi_r",
  n_genes = as.integer(length(distances)),
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam_bad,
  dscale = as.double(c(2, 4)),
  rdi = rdi
)
print(res3$rdi)

# [compute_rdi_r] Case 4: distances with NaN (skipped)
# distances_nan <- as.double(c(1, NaN, 3, 4, 5))
# rdi <- double(length(distances_nan))
# res4 <- .Fortran("compute_rdi_r",
#   n_genes = as.integer(length(distances_nan)),
#   n_families = n_families,
#   distances = distances_nan,
#   gene_to_fam = gene_to_fam,
#   dscale = as.double(c(2, 4)),
#   rdi = rdi
# )
# print(res4$rdi)

##
# [compute_rdi_r] Case 5: dscale shorter than max(gene_to_fam)
# Expectation: Should error or produce -1 for genes with out-of-bounds family index.
cat("\n[compute_rdi_r] Case 5: dscale shorter than max(gene_to_fam)\n")
dscale_short <- as.double(2) # only 1 value, but there are families 1 and 2
rdi <- double(length(distances))
tryCatch({
  res5 <- .Fortran("compute_rdi_r",
    n_genes = as.integer(length(distances)),
    n_families = n_families,
    distances = distances,
    gene_to_fam = gene_to_fam,
    dscale = dscale_short,
    rdi = rdi
  )
  print(res5$rdi)
}, error = function(e) cat("Unexpected error:", e$message, "\n"))

##
# [compute_rdi_r] Case 6: vectors with incompatible length
# Expectation: Should error due to mismatched vector lengths.
cat("\n[compute_rdi_r] Case 6: vectors with incompatible length\n")
distances_short <- as.double(c(1, 2, 3))
gene_to_fam_short <- as.integer(c(1, 1, 2))
rdi <- double(length(distances_short))
tryCatch({
  res6 <- .Fortran("compute_rdi_r",
    n_genes = as.integer(length(distances_short)),
    n_families = n_families,
    distances = distances_short,
    gene_to_fam = gene_to_fam_short,
    dscale = as.double(c(2, 4)),
    rdi = rdi
  )
  print(res6$rdi)
}, error = function(e) cat("Unexpected error:", e$message, "\n"))


# =====================
# End of compute_rdi_r tests
# =====================

# =====================
# Test cases for identify_outliers_r Fortran wrapper
# =====================

##
# [identify_outliers_r] Case 1: Simple RDI, percentile 50
# Expectation: Top 50% (highest 3) should be outliers (TRUE), others FALSE.
cat("\n[identify_outliers_r] Case 1: Simple RDI, percentile 50\n")
rdi <- as.double(c(0.1, 0.2, 0.3, 0.4, 0.5))
sorted_rdi <- double(length(rdi))
perm <- seq_along(rdi)
stack_left <- integer(length(rdi))
stack_right <- integer(length(rdi))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 50.0
res1 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(res1$is_outlier)
print(res1$threshold)

##
# [identify_outliers_r] Case 2: RDI with negative values (should be ignored)
# Expectation: Negative RDI values are ignored for outlier detection; only positive values considered.
cat("\n[identify_outliers_r] Case 2: RDI with negative values (should be ignored)\n")
rdi <- as.double(c(-1, 0.2, 0.3, 0.4, 0.5))
sorted_rdi <- double(length(rdi))
perm <- seq_along(rdi)
stack_left <- integer(length(rdi))
stack_right <- integer(length(rdi))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 80.0
res2 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(res2$is_outlier)
print(res2$threshold)

##
# [identify_outliers_r] Case 3: All RDI zeros
# Expectation: No outliers detected; all is_outlier should be FALSE.
cat("\n[identify_outliers_r] Case 3: All RDI zeros\n")
rdi <- as.double(rep(0, 5))
sorted_rdi <- double(length(rdi))
perm <- seq_along(rdi)
stack_left <- integer(length(rdi))
stack_right <- integer(length(rdi))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 90.0
res3 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(res3$is_outlier)
print(res3$threshold)

##
# [identify_outliers_r] Case 4: Percentile 0 (all outliers)
# Expectation: All genes should be outliers (all TRUE).
cat("\n[identify_outliers_r] Case 4: Percentile 0 (all outliers)\n")
rdi <- as.double(c(0.1, 0.2, 0.3, 0.4, 0.5))
sorted_rdi <- double(length(rdi))
perm <- seq_along(rdi)
stack_left <- integer(length(rdi))
stack_right <- integer(length(rdi))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 0.0
res5 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(res5$is_outlier)
print(res5$threshold)

##
# [identify_outliers_r] Case 5: Percentile 100 (1 outlier)
# Expectation: Only the highest RDI should be outlier (TRUE), rest FALSE.
cat("\n[identify_outliers_r] Case 5: Percentile 100 (1 outlier)\n")
rdi <- as.double(c(0.1, 0.2, 0.3, 0.4, 0.5))
sorted_rdi <- double(length(rdi))
perm <- seq_along(rdi)
stack_left <- integer(length(rdi))
stack_right <- integer(length(rdi))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 100.0
res6 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(res6$is_outlier)
print(res6$threshold)


# =====================
# End of identify_outliers_r tests
# =====================


# =====================
# Test cases for detect_outliers_r Fortran wrapper
# =====================

##
# [detect_outliers_r] Case 1: Typical input (LOESS only)
# Expectation: At least one gene should be detected as outlier (TRUE), error_code should be 0.
cat("\n[detect_outliers_r] Case 1: Typical input (LOESS only)\n")

# Re-initialize all relevant arrays for Case 1
n_genes <- as.integer(6)
n_families <- as.integer(2)
distances <- as.double(c(1, 2, 3, 4, 5, 6))
gene_to_fam <- as.integer(c(1, 1, 2, 2, 2, 2))
work_array <- double(n_genes)
perm <- seq_len(n_genes)
stack_left <- integer(n_genes)
stack_right <- integer(n_genes)
is_outlier <- logical(n_genes)
loess_x <- double(n_families)
loess_y <- double(n_families)
loess_n <- integer(n_families)
workspace_weights <- double(n_families)
workspace_values <- matrix(0, 1, n_families)
error_code <- integer(1)
percentile <- 80.0
res1 <- .Fortran("detect_outliers_r",
  n_genes = n_genes,
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam,
  work_array = work_array,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right,
  is_outlier = is_outlier,
  loess_x = loess_x,
  loess_y = loess_y,
  loess_n = loess_n,
  workspace_weights = workspace_weights,
  workspace_values = workspace_values,
  error_code = error_code,
  percentile = as.double(percentile)
)
print(res1$is_outlier)
print(res1$error_code)

##
# [detect_outliers_r] Case 2: Invalid gene_to_fam indices (should return error_code -2)
# Expectation: All is_outlier should be FALSE, error_code should be -2.
cat("\n[detect_outliers_r] Case 2: Invalid gene_to_fam indices (should return error_code -2)\n")

# Re-initialize all relevant arrays for Case 2
gene_to_fam_invalid <- as.integer(c(1, 3, 2, 2, 2, 2))
is_outlier <- rep(FALSE, n_genes)
work_array <- double(n_genes)
perm <- seq_len(n_genes)
stack_left <- integer(n_genes)
stack_right <- integer(n_genes)
loess_x <- double(n_families)
loess_y <- double(n_families)
loess_n <- integer(n_families)
workspace_weights <- double(n_families)
workspace_values <- matrix(0, 1, n_families)
error_code <- integer(1)

res2 <- .Fortran("detect_outliers_r",
  n_genes = n_genes,
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam_invalid,
  work_array = work_array,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right,
  is_outlier = is_outlier,
  loess_x = loess_x,
  loess_y = loess_y,
  loess_n = loess_n,
  workspace_weights = workspace_weights,
  workspace_values = workspace_values,
  error_code = error_code,
  percentile = as.double(percentile)
)
print(res2$is_outlier)
print(res2$error_code)

# =====================
# End of detect_outliers_r tests
# =====================
