
# =====================
# Test cases for compute_family_scaling_r (main/alloc) and compute_family_scaling_expert_r (expert/kernel) Fortran wrappers
# Each case checks a different combination of optional arguments and error handling

dyn.load("build/libtensor-omics.so")

# Common test data
n_genes <- as.integer(5)
n_families <- as.integer(2)
distances <- as.double(c(1, 2, 3, 4, 5))
gene_to_fam <- as.integer(c(1, 1, 2, 2, 2))

# =====================
# Tests for compute_family_scaling_r (main/alloc interface - user friendly)
# =====================

cat("=== Testing compute_family_scaling_r (main/alloc interface) ===\n")

##
# [compute_family_scaling_r] Case 1: Basic LOESS scaling
# Expectation: dscale should be positive for both families (since both have >1 gene), error_code should be 0.
cat("\n[compute_family_scaling_r] Case 1: Basic LOESS scaling\n")
dscale <- double(n_families)
loess_x <- double(n_families)
loess_y <- double(n_families)
indices_used <- integer(n_families)
error_code <- integer(1)
result1 <- .Fortran("compute_family_scaling_r",
  n_genes = n_genes,
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam,
  dscale = dscale,
  loess_x = loess_x,
  loess_y = loess_y,
  indices_used = indices_used,
  error_code = error_code
)
print(paste("dscale:", paste(result1$dscale, collapse=", ")))
print(paste("error_code:", result1$error_code))
print(paste("loess_x (medians):", paste(result1$loess_x, collapse=", ")))
print(paste("loess_y (stddevs):", paste(result1$loess_y, collapse=", ")))

##
# [compute_family_scaling_r] Case 2: Invalid family indices
# Expectation: dscale should be all -1, error_code should be -2 (invalid family index in gene_to_fam).
cat("\n[compute_family_scaling_r] Case 2: Invalid family indices\n")
gene_to_fam_invalid <- as.integer(c(1, 3, 2, 2, 2))
dscale <- double(n_families)
loess_x <- double(n_families)
loess_y <- double(n_families)
indices_used <- integer(n_families)
error_code <- integer(1)
result2 <- .Fortran("compute_family_scaling_r",
  n_genes = n_genes,
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam_invalid,
  dscale = dscale,
  loess_x = loess_x,
  loess_y = loess_y,
  indices_used = indices_used,
  error_code = error_code
)
print(paste("dscale:", paste(result2$dscale, collapse=", ")))
print(paste("error_code:", result2$error_code))

##
# [compute_family_scaling_r] Case 3: Family with only one gene
# Expectation: dscale for the single-gene family should be 0, error_code should be 0.
cat("\n[compute_family_scaling_r] Case 3: Family with only one gene\n")
gene_to_fam_single <- as.integer(c(1, 2, 2, 2, 2))
dscale <- double(n_families)
loess_x <- double(n_families)
loess_y <- double(n_families)
indices_used <- integer(n_families)
error_code <- integer(1)
result3 <- .Fortran("compute_family_scaling_r",
  n_genes = n_genes,
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam_single,
  dscale = dscale,
  loess_x = loess_x,
  loess_y = loess_y,
  indices_used = indices_used,
  error_code = error_code
)
print(paste("dscale:", paste(result3$dscale, collapse=", ")))
print(paste("error_code:", result3$error_code))

##
# [compute_family_scaling_r] Case 4: All distances zero
# Expectation: dscale should be 0 for all families, error_code should be 0.
cat("\n[compute_family_scaling_r] Case 4: All distances zero\n")
distances_zero <- as.double(c(0, 0, 0, 0, 0))
dscale <- double(n_families)
loess_x <- double(n_families)
loess_y <- double(n_families)
indices_used <- integer(n_families)
error_code <- integer(1)
result4 <- .Fortran("compute_family_scaling_r",
  n_genes = n_genes,
  n_families = n_families,
  distances = distances_zero,
  gene_to_fam = gene_to_fam,
  dscale = dscale,
  loess_x = loess_x,
  loess_y = loess_y,
  indices_used = indices_used,
  error_code = error_code
)
print(paste("dscale:", paste(result4$dscale, collapse=", ")))
print(paste("error_code:", result4$error_code))

##
# [compute_family_scaling_r] Case 5: Large dataset
# Expectation: Should handle larger datasets efficiently, error_code should be 0.
cat("\n[compute_family_scaling_r] Case 5: Large dataset\n")
n_genes_large <- as.integer(20)
n_families_large <- as.integer(4)
distances_large <- as.double(runif(n_genes_large, 0.5, 5.0))  # Random distances
gene_to_fam_large <- as.integer(sample(1:n_families_large, n_genes_large, replace=TRUE))
dscale <- double(n_families_large)
loess_x <- double(n_families_large)
loess_y <- double(n_families_large)
indices_used <- integer(n_families_large)
error_code <- integer(1)
result5 <- .Fortran("compute_family_scaling_r",
  n_genes = n_genes_large,
  n_families = n_families_large,
  distances = distances_large,
  gene_to_fam = gene_to_fam_large,
  dscale = dscale,
  loess_x = loess_x,
  loess_y = loess_y,
  indices_used = indices_used,
  error_code = error_code
)
print(paste("dscale:", paste(round(result5$dscale, 3), collapse=", ")))
print(paste("error_code:", result5$error_code))

##
# [compute_family_scaling_r] Case 6: Mixed family sizes
# Expectation: Should handle families of different sizes correctly, error_code should be 0.
cat("\n[compute_family_scaling_r] Case 6: Mixed family sizes\n")
n_genes_mixed <- as.integer(10)
n_families_mixed <- as.integer(3)
distances_mixed <- as.double(c(1.0, 1.1, 1.2, 1.3, 1.4, 2.0, 2.1, 3.0, 4.0, 5.0))
gene_to_fam_mixed <- as.integer(c(1, 1, 1, 1, 1, 2, 2, 3, 3, 3))  # Family 1: 5 genes, Family 2: 2 genes, Family 3: 3 genes
dscale <- double(n_families_mixed)
loess_x <- double(n_families_mixed)
loess_y <- double(n_families_mixed)
indices_used <- integer(n_families_mixed)
error_code <- integer(1)
result6 <- .Fortran("compute_family_scaling_r",
  n_genes = n_genes_mixed,
  n_families = n_families_mixed,
  distances = distances_mixed,
  gene_to_fam = gene_to_fam_mixed,
  dscale = dscale,
  loess_x = loess_x,
  loess_y = loess_y,
  indices_used = indices_used,
  error_code = error_code
)
print(paste("dscale:", paste(round(result6$dscale, 3), collapse=", ")))
print(paste("error_code:", result6$error_code))
print(paste("Family medians:", paste(round(result6$loess_x, 3), collapse=", ")))

# =====================
# Tests for compute_family_scaling_expert_r (expert/kernel interface - advanced users)
# =====================

cat("\n=== Testing compute_family_scaling_expert_r (expert/kernel interface) ===\n")

# For expert interface, user must provide all work arrays
perm_tmp <- integer(n_genes)
stack_left_tmp <- integer(n_genes)
stack_right_tmp <- integer(n_genes)
family_distances <- double(n_genes)

##
# [compute_family_scaling_expert_r] Case 1: Basic LOESS scaling with user-provided work arrays
# Expectation: Same results as main interface, but user controls memory allocation.
cat("\n[compute_family_scaling_expert_r] Case 1: Basic LOESS scaling with work arrays\n")
dscale <- double(n_families)
loess_x <- double(n_families)
loess_y <- double(n_families)
indices_used <- integer(n_families)
error_code <- integer(1)
result_expert1 <- .Fortran("compute_family_scaling_expert_r",
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
  family_distances = family_distances,
  error_code = error_code
)
print(paste("dscale:", paste(result_expert1$dscale, collapse=", ")))
print(paste("error_code:", result_expert1$error_code))

##
# [compute_family_scaling_expert_r] Case 2: Check that work arrays are modified
# Expectation: Work arrays should contain meaningful values after computation.
cat("\n[compute_family_scaling_expert_r] Case 2: Work arrays modification check\n")
# Reset work arrays to known values
perm_tmp <- integer(n_genes)
stack_left_tmp <- integer(n_genes)
stack_right_tmp <- integer(n_genes)
family_distances <- double(n_genes)
print(paste("family_distances before:", paste(family_distances, collapse=", ")))

dscale <- double(n_families)
loess_x <- double(n_families)
loess_y <- double(n_families)
indices_used <- integer(n_families)
error_code <- integer(1)
result_expert2 <- .Fortran("compute_family_scaling_expert_r",
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
  family_distances = family_distances,
  error_code = error_code
)
print(paste("family_distances after:", paste(round(result_expert2$family_distances, 3), collapse=", ")))
print(paste("perm_tmp:", paste(result_expert2$perm_tmp, collapse=", ")))

##
# [compute_family_scaling_expert_r] Case 3: Performance comparison with main interface
# Expectation: Expert interface should give identical results but with user-controlled memory.
cat("\n[compute_family_scaling_expert_r] Case 3: Results comparison with main interface\n")
# Use same test data as Case 1 of main interface
perm_tmp <- integer(n_genes)
stack_left_tmp <- integer(n_genes)
stack_right_tmp <- integer(n_genes)
family_distances <- double(n_genes)
dscale_expert <- double(n_families)
loess_x_expert <- double(n_families)
loess_y_expert <- double(n_families)
indices_used_expert <- integer(n_families)
error_code_expert <- integer(1)

result_expert3 <- .Fortran("compute_family_scaling_expert_r",
  n_genes = n_genes,
  n_families = n_families,
  distances = distances,
  gene_to_fam = gene_to_fam,
  dscale = dscale_expert,
  loess_x = loess_x_expert,
  loess_y = loess_y_expert,
  indices_used = indices_used_expert,
  perm_tmp = perm_tmp,
  stack_left_tmp = stack_left_tmp,
  stack_right_tmp = stack_right_tmp,
  family_distances = family_distances,
  error_code = error_code_expert
)

# Compare with result1 from main interface
print("Comparison with main interface:")
print(paste("Main dscale:", paste(result1$dscale, collapse=", ")))
print(paste("Expert dscale:", paste(result_expert3$dscale, collapse=", ")))
print(paste("Difference:", paste(round(abs(result1$dscale - result_expert3$dscale), 6), collapse=", ")))
print(paste("Results identical:", all(abs(result1$dscale - result_expert3$dscale) < 1e-10)))



# =====================
# Test cases for compute_rdi_r Fortran wrapper
# =====================

cat("\n=== Testing compute_rdi_r ===\n")

##
# [compute_rdi_r] Case 1: normal input
# Expectation: rdi = abs(distances) / dscale; should be c(0.5, 1, 0.75, 1, 1.25) for dscale = c(2,4).
cat("\n[compute_rdi_r] Case 1: normal input\n")
distances <- as.double(c(1, 2, 3, 4, 5))
gene_to_fam <- as.integer(c(1, 1, 2, 2, 2))
n_families_rdi <- as.integer(2)
dscale <- as.double(c(2, 4))
rdi <- double(length(distances))
sorted_rdi <- double(length(distances))
perm <- as.integer(seq_along(distances))
stack_left <- integer(length(distances))
stack_right <- integer(length(distances))
res1 <- .Fortran("compute_rdi_r",
  n_genes = as.integer(length(distances)),
  n_families = n_families_rdi,
  distances = distances,
  gene_to_fam = gene_to_fam,
  dscale = dscale,
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right
)
print(paste("RDI:", paste(round(res1$rdi, 3), collapse=", ")))
print(paste("Sorted RDI:", paste(round(res1$sorted_rdi, 3), collapse=", ")))

##
# [compute_rdi_r] Case 2: dscale with zeros
# Expectation: All rdi values should be 0 (since scaling is zero for all families).
cat("\n[compute_rdi_r] Case 2: dscale with zeros\n")
dscale_zero <- as.double(c(0, 0))
rdi <- double(length(distances))
sorted_rdi <- double(length(distances))
perm <- as.integer(seq_along(distances))
stack_left <- integer(length(distances))
stack_right <- integer(length(distances))
res2 <- .Fortran("compute_rdi_r",
  n_genes = as.integer(length(distances)),
  n_families = n_families_rdi,
  distances = distances,
  gene_to_fam = gene_to_fam,
  dscale = dscale_zero,
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right
)
print(paste("RDI with zero scaling:", paste(res2$rdi, collapse=", ")))

##
# [compute_rdi_r] Case 3: gene_to_fam out of range
# Expectation: rdi for genes with invalid family index should be 0, others normal.
cat("\n[compute_rdi_r] Case 3: gene_to_fam out of range\n")
gene_to_fam_bad <- as.integer(c(1, 3, 2, 2, 2)) # 3 doesn't exist in dscale
rdi <- double(length(distances))
sorted_rdi <- double(length(distances))
perm <- as.integer(seq_along(distances))
stack_left <- integer(length(distances))
stack_right <- integer(length(distances))
res3 <- .Fortran("compute_rdi_r",
  n_genes = as.integer(length(distances)),
  n_families = n_families_rdi,
  distances = distances,
  gene_to_fam = gene_to_fam_bad,
  dscale = as.double(c(2, 4)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right
)
print(paste("RDI with invalid family index:", paste(round(res3$rdi, 3), collapse=", ")))

##
# [compute_rdi_r] Case 4: High precision test
# Expectation: RDI calculations should be precise and match expected values.
cat("\n[compute_rdi_r] Case 4: High precision test\n")
distances_precise <- as.double(c(2.0, 4.0, 6.0))
gene_to_fam_precise <- as.integer(c(1, 1, 1))
n_families_precise <- as.integer(1)
dscale_precise <- as.double(c(2.0))  # All genes in family 1, scaling = 2.0
rdi <- double(length(distances_precise))
sorted_rdi <- double(length(distances_precise))
perm <- as.integer(seq_along(distances_precise))
stack_left <- integer(length(distances_precise))
stack_right <- integer(length(distances_precise))
res4 <- .Fortran("compute_rdi_r",
  n_genes = as.integer(length(distances_precise)),
  n_families = n_families_precise,
  distances = distances_precise,
  gene_to_fam = gene_to_fam_precise,
  dscale = dscale_precise,
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right
)
# Expected RDI: [2.0/2.0, 4.0/2.0, 6.0/2.0] = [1.0, 2.0, 3.0]
print(paste("Expected RDI: [1.0, 2.0, 3.0]"))
print(paste("Actual RDI:", paste(res4$rdi, collapse=", ")))
print(paste("High precision test passed:", all(abs(res4$rdi - c(1.0, 2.0, 3.0)) < 1e-10)))

##
# [compute_rdi_r] Case 5: Negative distances
# Expectation: Should handle negative distances appropriately.
cat("\n[compute_rdi_r] Case 5: Negative distances\n")
distances_negative <- as.double(c(-1, 2, -3, 4, 5))
rdi <- double(length(distances_negative))
sorted_rdi <- double(length(distances_negative))
perm <- as.integer(seq_along(distances_negative))
stack_left <- integer(length(distances_negative))
stack_right <- integer(length(distances_negative))
res5 <- .Fortran("compute_rdi_r",
  n_genes = as.integer(length(distances_negative)),
  n_families = n_families_rdi,
  distances = distances_negative,
  gene_to_fam = gene_to_fam,
  dscale = as.double(c(2, 4)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right
)
print(paste("RDI with negative distances:", paste(round(res5$rdi, 3), collapse=", ")))


# =====================
# Test cases for identify_outliers_r Fortran wrapper
# =====================

cat("\n=== Testing identify_outliers_r ===\n")

##
# [identify_outliers_r] Case 1: Simple RDI, percentile 50
# Expectation: Top 50% (highest 3) should be outliers (TRUE), others FALSE.
cat("\n[identify_outliers_r] Case 1: Simple RDI, percentile 50\n")
rdi <- as.double(c(0.3, 0.1, 0.5, 0.2, 0.4)) # intentionally unsorted
sorted_rdi <- sort(rdi[rdi >= 0])
if (length(sorted_rdi) < length(rdi)) sorted_rdi <- c(sorted_rdi, rep(0, length(rdi) - length(sorted_rdi)))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 50.0
res1 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(paste("RDI:", paste(rdi, collapse=", ")))
print(paste("Is outlier:", paste(res1$is_outlier, collapse=", ")))
print(paste("Threshold:", res1$threshold))

##
# [identify_outliers_r] Case 2: RDI with negative values (should be ignored)
# Expectation: Negative RDI values are ignored for outlier detection; only positive values considered.
cat("\n[identify_outliers_r] Case 2: RDI with negative values (should be ignored)\n")
rdi <- as.double(c(0.3, -1, 0.5, 0.2, 0.4)) # intentionally unsorted, with negative
sorted_rdi <- sort(rdi[rdi >= 0])
if (length(sorted_rdi) < length(rdi)) sorted_rdi <- c(sorted_rdi, rep(0, length(rdi) - length(sorted_rdi)))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 80.0
res2 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(paste("RDI with negative:", paste(rdi, collapse=", ")))
print(paste("Is outlier:", paste(res2$is_outlier, collapse=", ")))
print(paste("Threshold:", res2$threshold))

##
# [identify_outliers_r] Case 3: All RDI zeros
# Expectation: No outliers detected; all is_outlier should be FALSE.
cat("\n[identify_outliers_r] Case 3: All RDI zeros\n")
rdi <- as.double(c(0, 0, 0, 0, 0)) # all zero
sorted_rdi <- sort(rdi[rdi >= 0])
if (length(sorted_rdi) < length(rdi)) sorted_rdi <- c(sorted_rdi, rep(0, length(rdi) - length(sorted_rdi)))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 90.0
res3 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(paste("All zeros outliers:", paste(res3$is_outlier, collapse=", ")))
print(paste("Threshold:", res3$threshold))

##
# [identify_outliers_r] Case 4: Percentile 0 (all outliers)
# Expectation: All genes should be outliers (all TRUE).
cat("\n[identify_outliers_r] Case 4: Percentile 0 (all outliers)\n")
rdi <- as.double(c(0.3, 0.1, 0.5, 0.2, 0.4)) # intentionally unsorted
sorted_rdi <- sort(rdi[rdi >= 0])
if (length(sorted_rdi) < length(rdi)) sorted_rdi <- c(sorted_rdi, rep(0, length(rdi) - length(sorted_rdi)))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 0.0
res4 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(paste("0% percentile outliers:", paste(res4$is_outlier, collapse=", ")))
print(paste("All are outliers:", all(res4$is_outlier)))

##
# [identify_outliers_r] Case 5: Percentile 100 (1 outlier)
# Expectation: Only the highest RDI should be outlier (TRUE), rest FALSE.
cat("\n[identify_outliers_r] Case 5: Percentile 100 (1 outlier)\n")
rdi <- as.double(c(0.3, 0.1, 0.5, 0.2, 0.4)) # intentionally unsorted
sorted_rdi <- sort(rdi[rdi >= 0])
if (length(sorted_rdi) < length(rdi)) sorted_rdi <- c(sorted_rdi, rep(0, length(rdi) - length(sorted_rdi)))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 100.0
res5 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(paste("100% percentile outliers:", paste(res5$is_outlier, collapse=", ")))
print(paste("Only highest is outlier:", sum(res5$is_outlier) == 1))

##
# [identify_outliers_r] Case 6: All RDI negative (should be ignored)
# Expectation: All negative RDI values should be ignored for outlier detection; all is_outlier should be FALSE.
cat("\n[identify_outliers_r] Case 6: All RDI negative (should be ignored)\n")
rdi <- as.double(c(-0.1, -0.2, -0.3, -0.4, -0.5)) # all negative
sorted_rdi <- sort(rdi[rdi >= 0])
if (length(sorted_rdi) < length(rdi)) sorted_rdi <- c(sorted_rdi, rep(0, length(rdi) - length(sorted_rdi)))
is_outlier <- logical(length(rdi))
threshold <- double(1)
percentile <- 80.0
res6 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  is_outlier = is_outlier,
  threshold = threshold,
  percentile = as.double(percentile)
)
print(paste("All negative outliers:", paste(res6$is_outlier, collapse=", ")))
print(paste("No outliers detected:", !any(res6$is_outlier)))

##
# [identify_outliers_r] Case 7: Default percentile test
# Expectation: Should use 95% as default percentile.
cat("\n[identify_outliers_r] Case 7: Default percentile test\n")
rdi <- as.double(c(0.1, 0.2, 0.3, 0.4, 0.5)) # sorted for clarity
sorted_rdi <- sort(rdi[rdi >= 0])
if (length(sorted_rdi) < length(rdi)) sorted_rdi <- c(sorted_rdi, rep(0, length(rdi) - length(sorted_rdi)))
is_outlier <- logical(length(rdi))
threshold <- double(1)
# Don't specify percentile to test default
res7 <- .Fortran("identify_outliers_r",
  n_genes = as.integer(length(rdi)),
  rdi = rdi,
  sorted_rdi = sorted_rdi,
  is_outlier = is_outlier,
  threshold = threshold
)
print(paste("Default percentile outliers:", paste(res7$is_outlier, collapse=", ")))
print(paste("Should detect only highest:", sum(res7$is_outlier)))

# =====================
# Test cases for detect_outliers_r Fortran wrapper (comprehensive workflow)
# =====================

cat("\n=== Testing detect_outliers_r (complete outlier detection workflow) ===\n")

##
# [detect_outliers_r] Case 1: Typical input (LOESS only)
# Expectation: At least one gene should be detected as outlier (TRUE), error_code should be 0.
cat("\n[detect_outliers_r] Case 1: Typical input\n")

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
  error_code = error_code,
  percentile = as.double(percentile)
)
print(paste("Distances:", paste(distances, collapse=", ")))
print(paste("Gene to family:", paste(gene_to_fam, collapse=", ")))
print(paste("Is outlier:", paste(res1$is_outlier, collapse=", ")))
print(paste("Error code:", res1$error_code))
print(paste("Family medians:", paste(round(res1$loess_x, 3), collapse=", ")))
print(paste("Family stddevs:", paste(round(res1$loess_y, 3), collapse=", ")))

##
# [detect_outliers_r] Case 2: Invalid gene_to_fam indices (should return error_code -2)
# Expectation: All is_outlier should be FALSE, error_code should be -2.
cat("\n[detect_outliers_r] Case 2: Invalid gene_to_fam indices\n")

# Re-initialize all relevant arrays for Case 2
gene_to_fam_invalid <- as.integer(c(1, 3, 2, 2, 2, 2))  # family 3 doesn't exist
is_outlier <- rep(FALSE, n_genes)
work_array <- double(n_genes)
perm <- seq_len(n_genes)
stack_left <- integer(n_genes)
stack_right <- integer(n_genes)
loess_x <- double(n_families)
loess_y <- double(n_families)
loess_n <- integer(n_families)
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
  error_code = error_code,
  percentile = as.double(percentile)
)
print(paste("Invalid family - Is outlier:", paste(res2$is_outlier, collapse=", ")))
print(paste("Invalid family - Error code:", res2$error_code))

##
# [detect_outliers_r] Case 3: Default percentile test
# Expectation: Should use 95% as default percentile.
cat("\n[detect_outliers_r] Case 3: Default percentile test\n")
n_genes <- as.integer(8)
n_families <- as.integer(2)
distances <- as.double(c(1, 1.1, 1.2, 1.3, 10, 10.1, 10.2, 50))  # Last gene is clear outlier
gene_to_fam <- as.integer(c(1, 1, 1, 1, 2, 2, 2, 2))
work_array <- double(n_genes)
perm <- seq_len(n_genes)
stack_left <- integer(n_genes)
stack_right <- integer(n_genes)
is_outlier <- logical(n_genes)
loess_x <- double(n_families)
loess_y <- double(n_families)
loess_n <- integer(n_families)
error_code <- integer(1)

# Don't specify percentile to test default behavior
res3 <- .Fortran("detect_outliers_r",
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
  error_code = error_code
)
print(paste("Default percentile - Distances:", paste(distances, collapse=", ")))
print(paste("Default percentile - Is outlier:", paste(res3$is_outlier, collapse=", ")))
print(paste("Default percentile - Error code:", res3$error_code))
print(paste("Outlier detected for extreme value:", res3$is_outlier[8]))

##
# [detect_outliers_r] Case 4: Single gene families
# Expectation: Single gene families should have scaling 0, no outliers detected.
cat("\n[detect_outliers_r] Case 4: Single gene families\n")
n_genes <- as.integer(3)
n_families <- as.integer(3)
distances <- as.double(c(1, 10, 100))  # Each gene in different family
gene_to_fam <- as.integer(c(1, 2, 3))
work_array <- double(n_genes)
perm <- seq_len(n_genes)
stack_left <- integer(n_genes)
stack_right <- integer(n_genes)
is_outlier <- logical(n_genes)
loess_x <- double(n_families)
loess_y <- double(n_families)
loess_n <- integer(n_families)
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
  error_code = error_code,
  percentile = as.double(90.0)
)
print(paste("Single families - Is outlier:", paste(res4$is_outlier, collapse=", ")))
print(paste("Single families - Error code:", res4$error_code))
print(paste("Single families - Family medians:", paste(res4$loess_x, collapse=", ")))
print(paste("No outliers expected:", !any(res4$is_outlier)))

##
# [detect_outliers_r] Case 5: Extreme percentile values
# Expectation: 0% percentile should detect all as outliers, 100% only the highest.
cat("\n[detect_outliers_r] Case 5: Extreme percentile values\n")
n_genes <- as.integer(5)
n_families <- as.integer(1)
distances <- as.double(c(1, 2, 3, 4, 5))
gene_to_fam <- as.integer(c(1, 1, 1, 1, 1))

# Test 0% percentile
work_array <- double(n_genes)
perm <- seq_len(n_genes)
stack_left <- integer(n_genes)
stack_right <- integer(n_genes)
is_outlier <- logical(n_genes)
loess_x <- double(n_families)
loess_y <- double(n_families)
loess_n <- integer(n_families)
error_code <- integer(1)

res5a <- .Fortran("detect_outliers_r",
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
  error_code = error_code,
  percentile = as.double(0.0)
)
print(paste("0% percentile - All outliers:", all(res5a$is_outlier)))

# Test 100% percentile
work_array <- double(n_genes)
perm <- seq_len(n_genes)
stack_left <- integer(n_genes)
stack_right <- integer(n_genes)
is_outlier <- logical(n_genes)
loess_x <- double(n_families)
loess_y <- double(n_families)
loess_n <- integer(n_families)
error_code <- integer(1)

res5b <- .Fortran("detect_outliers_r",
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
  error_code = error_code,
  percentile = as.double(100.0)
)
print(paste("100% percentile - Only one outlier:", sum(res5b$is_outlier) == 1))
print(paste("100% percentile - Highest is outlier:", res5b$is_outlier[5]))

##
# [detect_outliers_r] Case 6: Large mixed dataset
# Expectation: Should handle large datasets with mixed family sizes efficiently.
cat("\n[detect_outliers_r] Case 6: Large mixed dataset\n")
n_genes_large <- as.integer(50)
n_families_large <- as.integer(5)
# Create synthetic data with some clear outliers
set.seed(123)  # For reproducibility
distances_large <- c(
  rnorm(10, 1, 0.1),    # Family 1: tight cluster around 1
  rnorm(10, 2, 0.1),    # Family 2: tight cluster around 2
  rnorm(10, 3, 0.1),    # Family 3: tight cluster around 3
  rnorm(10, 4, 0.1),    # Family 4: tight cluster around 4
  rnorm(8, 5, 0.1),     # Family 5: tight cluster around 5
  c(20, 25)             # Two clear outliers
)
gene_to_fam_large <- c(rep(1:5, each=10), c(1, 2))

work_array <- double(n_genes_large)
perm <- seq_len(n_genes_large)
stack_left <- integer(n_genes_large)
stack_right <- integer(n_genes_large)
is_outlier <- logical(n_genes_large)
loess_x <- double(n_families_large)
loess_y <- double(n_families_large)
loess_n <- integer(n_families_large)
error_code <- integer(1)

res6 <- .Fortran("detect_outliers_r",
  n_genes = n_genes_large,
  n_families = n_families_large,
  distances = distances_large,
  gene_to_fam = gene_to_fam_large,
  work_array = work_array,
  perm = perm,
  stack_left = stack_left,
  stack_right = stack_right,
  is_outlier = is_outlier,
  loess_x = loess_x,
  loess_y = loess_y,
  loess_n = loess_n,
  error_code = error_code,
  percentile = as.double(90.0)
)
print(paste("Large dataset - Error code:", res6$error_code))
print(paste("Large dataset - Outliers detected:", sum(res6$is_outlier)))
print(paste("Large dataset - Clear outliers detected:", res6$is_outlier[49] && res6$is_outlier[50]))

# =====================
# End of all tests
# =====================

cat("\n=== All tests completed ===\n")
cat("Summary:\n")
cat("- compute_family_scaling_r (main/alloc): User-friendly interface with automatic memory management\n")
cat("- compute_family_scaling_expert_r (expert/kernel): Advanced interface requiring user-provided work arrays\n")
cat("- compute_rdi_r: Relative Distance Index calculation\n")
cat("- identify_outliers_r: Outlier identification based on RDI percentiles\n")
cat("- detect_outliers_r: Complete workflow from distances to outlier detection\n")
