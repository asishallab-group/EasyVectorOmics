dyn.load("build/tox_normalization.so") # Ajusta el path

# Test para enteros
test_sort_integer <- function() {
  data <- as.integer(c(10, 3, 7, 1))
  n <- as.integer(length(data))
  perm <- as.integer(1:n)
  stack_left <- integer(20)
  stack_right <- integer(20)
  res <- .Fortran("sort_array_integer_r",
                  data, perm, stack_left, stack_right)
  sorted <- data[res[[2]]]
  expected <- c(1, 3, 7, 10)
  print(sorted)
  stopifnot(all(sorted == expected))
  cat("PASS: test_sort_integer\n")
}

# test_sort_integer()