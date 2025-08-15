# Load shared library
dyn.load("build/libtensor-omics.so")

omics_vector_RAP_projection <- function(vecs, vecs_selection_mask, axes_selection_mask) {
  n_selected_vecs <- sum(vecs_selection_mask == 1)
  n_selected_axes <- sum(axes_selection_mask == 1)
  res <- .Fortran("omics_vector_RAP_projection_r",
                  vecs = as.double(vecs),
                  n_axes = as.integer(nrow(vecs)),
                  n_vecs = as.integer(ncol(vecs)),
                  vecs_selection_mask = as.integer(vecs_selection_mask),
                  n_selected_vecs = n_selected_vecs,
                  axes_selection_mask = as.integer(axes_selection_mask),
                  n_selected_axes = n_selected_axes,
                  projections = matrix(as.double(1), nrow = n_selected_axes, ncol = n_selected_vecs))
  return(res$projections)
}

omics_field_RAP_projection <- function(vecs, vecs_selection_mask, axes_selection_mask) {
  n_selected_vecs <- sum(vecs_selection_mask == 1)
  n_selected_axes <- sum(axes_selection_mask == 1)
  res <- .Fortran("omics_field_RAP_projection_r",
                  vecs = as.double(vecs),
                  n_axes = as.integer(nrow(vecs) / 2),
                  n_vecs = as.integer(ncol(vecs)),
                  vecs_selection_mask = as.integer(vecs_selection_mask),
                  n_selected_vecs = n_selected_vecs,
                  axes_selection_mask = as.integer(axes_selection_mask),
                  n_selected_axes = n_selected_axes,
                  projections = matrix(as.double(1), nrow = n_selected_axes, ncol = n_selected_vecs))
  return(res$projections)
}

test_omics_vector_RAP_projection_call <- function() {
  n_axes <- 5
  n_selected_axes <- 2
  n_vecs <- 10
  n_selected_vecs <- 5

  vecs <- matrix(runif(n_axes * n_vecs), nrow = n_axes)

  vecs_selection_mask <- integer(n_vecs)
  vecs_selection_mask[sample(n_vecs, n_selected_vecs)] <- 1
  n_selected_vecs <- sum(vecs_selection_mask)

  axes_selection_mask <- integer(n_axes)
  axes_selection_mask[sample(n_axes, n_selected_axes)] <- 1

  projections <- omics_vector_RAP_projection(vecs, vecs_selection_mask, axes_selection_mask)

  for (i_vec in seq_len(n_selected_vecs)) {
    col_sum <- sum(projections[, i_vec])
    stopifnot(abs(col_sum) < 1e-6)
  }
}

test_omics_field_RAP_projection_call <- function() {
  n_axes <- 5
  n_selected_axes <- 2
  n_vecs <- 10
  n_selected_vecs <- 5

  vecs <- matrix(runif(2 * n_axes * n_vecs), nrow = 2 * n_axes)

  vecs_selection_mask <- integer(n_vecs)
  vecs_selection_mask[sample(n_vecs, n_selected_vecs)] <- 1
  n_selected_vecs <- sum(vecs_selection_mask)

  axes_selection_mask <- integer(n_axes)
  axes_selection_mask[sample(n_axes, n_selected_axes)] <- 1

  projections <- omics_field_RAP_projection(vecs, vecs_selection_mask, axes_selection_mask)

  for (i_vec in seq_len(n_selected_vecs)) {
    col_sum <- sum(projections[, i_vec])
    stopifnot(abs(col_sum) < 1e-6)
  }
}


test_omics_vector_RAP_projection_call()
test_omics_field_RAP_projection_call()
cat("Done\n")
