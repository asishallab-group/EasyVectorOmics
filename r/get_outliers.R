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

# 1. Both is_ortholog and max_distance_bw_orths present
cat("\n[compute_family_scaling_r] Case 1: Both is_ortholog and max_distance_bw_orths present\n")
is_ortholog <- as.logical(c(TRUE, FALSE, TRUE, FALSE, TRUE))
max_distance_bw_orths <- as.double(c(10, 20))
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
  error_code = error_code,
  is_ortholog = is_ortholog,
  max_distance_bw_orths = max_distance_bw_orths
)
print(result1$dscale)
print(result1$error_code)

# 2. Both is_ortholog and max_distance_bw_orths empty (LOESS only)
cat("\n[compute_family_scaling_r] Case 2: Both is_ortholog and max_distance_bw_orths empty (LOESS only)\n")
is_ortholog <- rep(FALSE, n_genes)
max_distance_bw_orths <- rep(0, n_families)
result2 <- .Fortran("compute_family_scaling_r",
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
  error_code = error_code,
  is_ortholog = is_ortholog,
  max_distance_bw_orths = max_distance_bw_orths
)
print(result2$dscale)
print(result2$error_code)

# 3. Only is_ortholog present, max_distance_bw_orths empty (should error if any orthologs)
cat("\n[compute_family_scaling_r] Case 3: Only is_ortholog present, max_distance_bw_orths empty\n")
is_ortholog <- as.logical(c(TRUE, FALSE, TRUE, FALSE, TRUE))
max_distance_bw_orths <- rep(0, n_families)
result3 <- .Fortran("compute_family_scaling_r",
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
  error_code = error_code,
  is_ortholog = is_ortholog,
  max_distance_bw_orths = max_distance_bw_orths
)
print(result3$dscale)
print(result3$error_code)

# 4. Only max_distance_bw_orths present, is_ortholog empty (should use LOESS)
cat("\n[compute_family_scaling_r] Case 4: Only max_distance_bw_orths present, is_ortholog empty\n")
is_ortholog <- rep(FALSE, n_genes)
max_distance_bw_orths <- as.double(c(10, 20))
result4 <- .Fortran("compute_family_scaling_r",
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
  error_code = error_code,
  is_ortholog = is_ortholog,
  max_distance_bw_orths = max_distance_bw_orths
)
print(result4$dscale)
print(result4$error_code)

# 5. Invalid family indices (should return error_code -2)
cat("\n[compute_family_scaling_r] Case 5: Invalid family indices\n")
gene_to_fam_invalid <- as.integer(c(1, 3, 2, 2, 2))
is_ortholog <- as.logical(c(TRUE, FALSE, TRUE, FALSE, TRUE))
max_distance_bw_orths <- as.double(c(10, 20))
result5 <- .Fortran("compute_family_scaling_r",
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
  error_code = error_code,
  is_ortholog = is_ortholog,
  max_distance_bw_orths = max_distance_bw_orths
)
print(result5$dscale)
print(result5$error_code)

# 6. Family with only one gene (should set dscale=0 for that family)
cat("\n[compute_family_scaling_r] Case 6: Family with only one gene\n")
gene_to_fam_single <- as.integer(c(1, 2, 2, 2, 2))
is_ortholog <- as.logical(c(TRUE, FALSE, TRUE, FALSE, TRUE))
max_distance_bw_orths <- as.double(c(10, 20))
result6 <- .Fortran("compute_family_scaling_r",
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
  error_code = error_code,
  is_ortholog = is_ortholog,
  max_distance_bw_orths = max_distance_bw_orths
)
print(result6$dscale)
print(result6$error_code)

# =====================
# Test cases for compute_rdi_r Fortran wrapper
# =====================

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

# cat("\n[compute_rdi_r] Case 4: distances with NaN\n")
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

cat("\n[detect_outliers_r] Case 1: Typical input, both is_ortholog and max_distance_bw_orths present\n")

# Re-initialize all relevant arrays for Case 1
n_genes <- as.integer(6)
n_families <- as.integer(2)
distances <- as.double(c(1, 2, 3, 4, 5, 6))
gene_to_fam <- as.integer(c(1, 1, 2, 2, 2, 2))
is_ortholog <- as.logical(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
max_distance_bw_orths <- as.double(c(10, 20))
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
  is_ortholog = is_ortholog,
  percentile = as.double(percentile),
  max_distance_bw_orths = max_distance_bw_orths
)
print(res1$is_outlier)
print(res1$error_code)

cat("\n[detect_outliers_r] Case 2: Only LOESS fallback (no orthologs, no max_distance_bw_orths)\n")

# Re-initialize all relevant arrays for Case 2
is_ortholog <- rep(FALSE, n_genes)
max_distance_bw_orths <- double(n_families)
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
  is_ortholog = is_ortholog,
  percentile = as.double(percentile)
)
print(res2$is_outlier)
print(res2$error_code)

cat("\n[detect_outliers_r] Case 3: Invalid gene_to_fam indices (should return error_code -2)\n")

# Re-initialize all relevant arrays for Case 3
gene_to_fam_invalid <- as.integer(c(1, 3, 2, 2, 2, 2))
is_ortholog <- rep(FALSE, n_genes)
max_distance_bw_orths <- as.double(c(10, 20))
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

res3 <- .Fortran("detect_outliers_r",
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
  is_ortholog = is_ortholog,
  percentile = as.double(percentile),
  max_distance_bw_orths = max_distance_bw_orths
)
print(res3$is_outlier)
print(res3$error_code)

cat("\n[detect_outliers_r] Case 4: Only max_distance_bw_orths present (should use LOESS)\n")

# Re-initialize all relevant arrays for Case 4
is_ortholog <- rep(FALSE, n_genes)
max_distance_bw_orths <- as.double(c(10, 20))
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

res4 <- .Fortran("detect_outliers_r",
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
  max_distance_bw_orths = max_distance_bw_orths,
  percentile = as.double(percentile)
)
print(res4$is_outlier)
print(res4$error_code)

cat("\n[detect_outliers_r] Case 5: Only is_ortholog present (should fallback to LOESS if no orthologs)\n")

# Re-initialize all relevant arrays for Case 5
is_ortholog <- rep(FALSE, n_genes)
max_distance_bw_orths <- double(n_families)
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

res5 <- .Fortran("detect_outliers_r",
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
  is_ortholog = is_ortholog,
  percentile = as.double(percentile)
)
print(res5$is_outlier)
print(res5$error_code)

# =====================
# End of detect_outliers_r tests
# =====================
