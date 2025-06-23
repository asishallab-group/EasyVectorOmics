dyn.load("build/tox_normalization.so") # Ajusta el path

# Test para enteros
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
  strings <- c("delta", "alpha", "beta")
  n <- length(strings)
  strlen <- max(nchar(strings))
  
  # Rellena a longitud fija
  padded <- format(strings, width = strlen, justify = "left")

  # Transforma en matriz de 1-char (n x strlen)
  char_matrix <- matrix("", nrow = n, ncol = strlen)
  for (i in 1:n) {
    char_matrix[i, ] <- strsplit(padded[i], "")[[1]]
  }

  # Transpone porque Fortran espera columna a columna
  flat_chars <- as.character(t(char_matrix))

  perm <- as.integer(seq_len(n))
  stack_left <- integer(n)
  stack_right <- integer(n)

  result <- .Fortran("sort_character_r",
    char_data = flat_chars,
    perm = perm,
    stack_left = stack_left,
    stack_right = stack_right,
    n = as.integer(n),
    strlen = as.integer(strlen)
  )

  sorted <- strings[result$perm]
  print(sorted)
}


test_sort_real()
test_sort_integer()
sort_test_char()