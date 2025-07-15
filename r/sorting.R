dyn.load("build/libtensor-omics.so")

test_sort_real <- function() {
    n <- 6
    array <- c(3.2, 1.5, 9.0, 2.1, 7.3, 4.4)
    perm <- as.integer(1:n)
    stack_left <- integer(n)
    stack_right <- integer(n)

    res <- .Fortran("sort_real_r",
                    as.double(array),
                    perm = as.integer(perm),
                    stack_left = as.integer(stack_left),
                    stack_right = as.integer(stack_right),
                    n = as.integer(n))

    sorted_indices <- res$perm
    print(array[sorted_indices])
}

test_sort_integer <- function() {
    n <- 6
    x <- c(5L, 3L, 8L, 1L, 4L, 2L)
    perm <- as.integer(1:n)
    stack_left <- integer(n)
    stack_right <- integer(n)

    res <- .Fortran("sort_integer_r",
                    array = as.integer(x),
                    perm = perm,
                    stack_left = stack_left,
                    stack_right = stack_right,
                    n = as.integer(n))

    sorted_x <- x[res$perm]
    print(sorted_x)
}

sort_test_char <- function() {
  strings <- c("test", "dog", "delta", "zeta", "alpha", "beta")
  n <- length(strings)
  strlen <- max(nchar(strings))

  # Create a character matrix of size (strlen x n), transposed for Fortran column-major order
  char_matrix <- matrix(" ", nrow = strlen, ncol = n)
  for (i in seq_len(n)) {
    # Pad each string to fixed length and split into individual characters
    padded <- sprintf("%-*s", strlen, strings[i])
    chars <- strsplit(padded, "")[[1]]
    char_matrix[, i] <- chars
  }

  # Flatten the matrix into a character vector
  char_vec <- as.vector(char_matrix)
  # Convert character vector to raw ASCII codes (needed for Fortran)
  char_raw <- as.raw(sapply(char_vec, charToRaw))

  # Prepare integer arrays
  perm <- as.integer(seq_len(n))
  stack_left <- integer(n)
  stack_right <- integer(n)

  # Call Fortran subroutine
  result <- .Fortran("sort_character_r",
    char_data = char_raw,
    perm = perm,
    stack_left = stack_left,
    stack_right = stack_right,
    n = as.integer(n),
    strlen = as.integer(strlen)
  )

  # Use the returned permutation to reorder the original strings
  sorted <- strings[result$perm]
  print(sorted)
}




test_sort_real()
test_sort_integer()
sort_test_char()