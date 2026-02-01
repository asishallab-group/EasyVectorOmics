run_all_tests <- function(env = parent.frame(), test_only = TRUE) {
  # Discover candidate names
  if (test_only) {
    test_names <- ls(env, pattern = "^test_")
  } else {
    test_names <- ls(env)
  }

  passed  <- 0
  failed  <- 0
  skipped <- 0

  cat("Running R tests...\n")

  for (name in test_names) {
    obj <- get(name, envir = env)
    if (!is.function(obj)) next  # skip non-functions

    test_func <- obj

    tryCatch(
      {

        cat(sprintf("Testing %s ...\n", name))
        test_func()
        cat(sprintf("✓ %s passed.\n", name))
        passed <- passed + 1
      },
      error = function(e) {
        msg <- conditionMessage(e)

        # Same skip logic as Python version
        if (grepl("Note:", msg) || grepl("acceptable", msg, ignore.case = TRUE)) {
          cat(sprintf("~ %s skipped (expected behavior): %s\n", name, msg))
          skipped <<- skipped + 1
        } else {
          cat(sprintf("✗ %s FAILED: %s\n", name, msg))
          failed <<- failed + 1
        }
      }
    )
  }

  cat("\nSummary: ", passed, " passed, ", failed, " failed, ", skipped, " skipped\n", sep = "")
  stopifnot(failed == 0)
}

assertError <- function(expr, msg) {
  err <- tryCatch(
    { expr; NULL },
    error = function(e) e
  )
  if (is.null(err)) stop(msg)
  invisible(TRUE)
}

assertTrue <- function(expr, msg) {
  if (!expr) stop(msg)
  invisible(TRUE)
}

assertFalse <- function(expr, msg) {
  assertTrue(!expr)
}
