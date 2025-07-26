# =====================
# Test cases for which_r Fortran wrapper
# =====================

dyn.load("build/libtensor-omics.so")

cat("\n[which_r] Case 1: Simple mask\n")
mask <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
n <- as.integer(length(mask))
m_max <- as.integer(n)
idx_out <- integer(m_max)
m_out <- integer(1)
res1 <- .Fortran("which_r",
  mask = as.logical(mask),
  n = n,
  idx_out = idx_out,
  m_max = m_max,
  m_out = m_out
)
print(res1$idx_out[1:res1$m_out])
print(res1$m_out)
# Expected: idx_out = c(1, 3, 0, 0, 0), m_out = 2

cat("\n[which_r] Case 2: All FALSE\n")
mask <- rep(FALSE, 5)
res2 <- .Fortran("which_r",
  mask = as.logical(mask),
  n = as.integer(length(mask)),
  idx_out = integer(length(mask)),
  m_max = as.integer(length(mask)),
  m_out = integer(1)
)
print(res2$idx_out[1:res2$m_out])
print(res2$m_out)
# Expected: idx_out = integer(0), m_out = 0

cat("\n[which_r] Case 3: All TRUE\n")
mask <- rep(TRUE, 5)
res3 <- .Fortran("which_r",
  mask = as.logical(mask),
  n = as.integer(length(mask)),
  idx_out = integer(length(mask)),
  m_max = as.integer(length(mask)),
  m_out = integer(1)
)
print(res3$idx_out[1:res3$m_out])
print(res3$m_out)
# Expected: idx_out = c(1, 2, 3, 4, 5), m_out = 5

cat("\n[which_r] Case 4: Empty mask\n")
mask <- logical(0)
res4 <- .Fortran("which_r",
  mask = as.logical(mask),
  n = as.integer(length(mask)),
  idx_out = integer(0),
  m_max = as.integer(0),
  m_out = integer(1)
)
print(res4$idx_out)
print(res4$m_out)
# Expected: idx_out = integer(0), m_out = 0

# =====================
# End of which_r tests
# =====================
