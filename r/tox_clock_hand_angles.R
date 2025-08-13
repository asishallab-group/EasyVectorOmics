# R interface and tests for clock hand angle calculations
# Uses Fortran wrappers for RAP projections and angle calculations

# Load the compiled Fortran library
dyn.load("build/libtensor-omics.so")  # Uncomment when library is built

# Constants
PI <- pi
TOL <- 1e-12

# ==================== WRAPPER FUNCTIONS ====================

#' Calculate signed clock hand angle between two normalized vectors
#' 
#' @param v1 First normalized vector in RAP space
#' @param v2 Second normalized vector in RAP space
#' @param selected_axes_for_signed Indices of 3 axes for directionality (ignored if dim <= 3)
#' @return Signed angle in radians [-π, π]
clock_hand_angle_between_vectors <- function(v1, v2, selected_axes_for_signed = c(1, 2, 3)) {
  n_dims <- length(v1)
  if (length(v2) != n_dims) {
    stop("Vectors must have the same dimension")
  }
  
  # Ensure selected_axes_for_signed has exactly 3 elements
  if (length(selected_axes_for_signed) != 3) {
    stop("selected_axes_for_signed must have exactly 3 elements")
  }
  
  # Call Fortran wrapper
  result <- .Fortran("clock_hand_angle_between_vectors_r",
                    v1 = as.double(v1),
                    v2 = as.double(v2),
                    n_dims = as.integer(n_dims),
                    signed_angle = as.double(0),
                    selected_axes_for_signed = as.integer(selected_axes_for_signed))
  
  return(result$signed_angle)
}

#' Calculate signed clock hand angles for multiple vector pairs
#' 
#' @param origins Matrix of origin vectors (n_dims x n_vecs)
#' @param targets Matrix of target vectors (n_dims x n_vecs)
#' @param vecs_selection_mask Logical vector indicating which pairs to compute
#' @param selected_axes_for_signed Indices of 3 axes for directionality
#' @return Vector of signed angles in radians [-π, π]
clock_hand_angles_for_shift_vectors <- function(origins, targets, 
                                               vecs_selection_mask = NULL,
                                               selected_axes_for_signed = c(1, 2, 3)) {
  if (!is.matrix(origins) || !is.matrix(targets)) {
    stop("origins and targets must be matrices")
  }
  
  if (!identical(dim(origins), dim(targets))) {
    stop("origins and targets must have the same dimensions")
  }
  
  n_dims <- nrow(origins)
  n_vecs <- ncol(origins)
  
  # Default selection mask (all TRUE)
  if (is.null(vecs_selection_mask)) {
    vecs_selection_mask <- rep(TRUE, n_vecs)
  }
  
  if (length(vecs_selection_mask) != n_vecs) {
    stop("vecs_selection_mask length must equal number of vector pairs")
  }
  
  n_selected_vecs <- sum(vecs_selection_mask)
  
  # Ensure selected_axes_for_signed has exactly 3 elements
  if (length(selected_axes_for_signed) != 3) {
    stop("selected_axes_for_signed must have exactly 3 elements")
  }
  
  # Call Fortran wrapper
  result <- .Fortran("clock_hand_angles_for_shift_vectors_r",
                    origins = as.double(origins),
                    targets = as.double(targets),
                    n_dims = as.integer(n_dims),
                    n_vecs = as.integer(n_vecs),
                    vecs_selection_mask = as.logical(vecs_selection_mask),
                    n_selected_vecs = as.integer(n_selected_vecs),
                    selected_axes_for_signed = as.integer(selected_axes_for_signed),
                    signed_angles = as.double(rep(0, n_selected_vecs)))
  
  return(result$signed_angles)
}

# ==================== TEST FUNCTIONS ====================

#' Assert that two values are approximately equal
assert_equal <- function(actual, expected, tolerance = TOL, message = "") {
  if (abs(actual - expected) > tolerance) {
    stop(sprintf("Test failed: %s. Expected %f, got %f (diff: %e)", 
                message, expected, actual, abs(actual - expected)))
  }
  cat(sprintf("✓ %s\n", message))
}

#' Assert that a condition is true
assert_true <- function(condition, message = "") {
  if (!condition) {
    stop(sprintf("Test failed: %s", message))
  }
  cat(sprintf("✓ %s\n", message))
}

# ==================== 2D TESTS ====================

test_identical_vectors_2d <- function() {
  v1 <- c(1.0, 0.0)
  v2 <- c(1.0, 0.0)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_equal(signed_angle, 0.0, TOL, "Identical 2D vectors should give 0 angle")
}

test_opposite_vectors_2d <- function() {
  v1 <- c(1.0, 0.0)
  v2 <- c(-1.0, 0.0)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_equal(abs(signed_angle), PI, TOL, "Opposite 2D vectors should give ±π")
}

test_perpendicular_vectors_2d <- function() {
  v1 <- c(1.0, 0.0)
  v2 <- c(0.0, 1.0)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_equal(abs(signed_angle), PI/2, TOL, "Perpendicular 2D vectors magnitude")
  assert_true(signed_angle > 0, "Counterclockwise rotation should be positive")
}

test_45_degree_rotation_2d <- function() {
  v1 <- c(1.0, 0.0)
  v2 <- c(sqrt(2)/2, sqrt(2)/2)  # 45 degrees
  expected <- PI/4
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_equal(signed_angle, expected, TOL, "45-degree counterclockwise rotation")
}

test_clockwise_vs_counterclockwise_2d <- function() {
  v1 <- c(1.0, 0.0)
  v2_ccw <- c(0.0, 1.0)   # 90° counterclockwise
  v2_cw <- c(0.0, -1.0)   # 90° clockwise
  
  angle_ccw <- clock_hand_angle_between_vectors(v1, v2_ccw)
  angle_cw <- clock_hand_angle_between_vectors(v1, v2_cw)
  
  assert_true(angle_ccw > 0, "Counterclockwise should be positive")
  assert_true(angle_cw < 0, "Clockwise should be negative")
  assert_equal(abs(angle_ccw), abs(angle_cw), TOL, "Magnitudes should be equal")
}

# ==================== 3D TESTS ====================

test_identical_vectors_3d <- function() {
  v1 <- c(1.0, 1.0, 1.0)
  v2 <- c(1.0, 1.0, 1.0)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_equal(signed_angle, 0.0, TOL, "Identical 3D vectors should give 0 angle")
}

test_perpendicular_vectors_3d <- function() {
  v1 <- c(1.0, 0.0, 0.0)
  v2 <- c(0.0, 1.0, 0.0)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_equal(abs(signed_angle), PI/2, TOL, "Perpendicular 3D vectors")
}

test_arbitrary_3d_rotation <- function() {
  v1 <- c(1.0, 2.0, 3.0)
  v2 <- c(2.0, 1.0, 3.0)
  
  # Normalize vectors
  v1 <- v1 / sqrt(sum(v1^2))
  v2 <- v2 / sqrt(sum(v2^2))
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  
  # Check magnitude is correct
  dot_product <- sum(v1 * v2)
  expected_magnitude <- acos(max(-1, min(1, dot_product)))
  assert_equal(abs(signed_angle), expected_magnitude, TOL, "3D arbitrary rotation magnitude")
}

# ==================== HIGH DIMENSIONAL TESTS ====================

test_high_dimensional_basic <- function() {
  v1 <- c(1.0, 0.0, 0.0, 0.0, 0.0)
  v2 <- c(0.0, 1.0, 0.0, 0.0, 0.0)
  selected_axes <- c(1, 2, 3)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, selected_axes)
  assert_equal(abs(signed_angle), PI/2, TOL, "High-dimensional perpendicular vectors")
}

test_high_dimensional_selected_axes <- function() {
  # Vectors that are perpendicular in dimensions 3 and 5
  v1 <- c(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
  v2 <- c(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
  selected_axes <- c(3, 5, 1)  # Use dimensions 3, 5, 1 for orientation
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2, selected_axes)
  assert_equal(abs(signed_angle), PI/2, TOL, "High-dimensional with selected axes")
}

# ==================== EDGE CASES ====================

test_denormalized_vectors <- function() {
  # Large magnitude vectors
  v1 <- c(100.0, 0.0)
  v2 <- c(0.0, 50.0)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_equal(abs(signed_angle), PI/2, TOL, "Denormalized vectors should work")
}

test_tiny_vectors_precision <- function() {
  tiny <- 1e-14
  v1 <- c(tiny, 0.0)
  v2 <- c(0.0, tiny)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_equal(abs(signed_angle), PI/2, 1e-10, "Tiny vectors precision")
}

test_huge_vectors_precision <- function() {
  huge_val <- 1e14
  v1 <- c(huge_val, 0.0)
  v2 <- c(0.0, huge_val)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_equal(abs(signed_angle), PI/2, TOL, "Huge vectors precision")
}

test_nearly_identical_vectors <- function() {
  epsilon <- 1e-15
  v1 <- c(1.0, 0.0)
  v2 <- c(1.0, epsilon)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_true(abs(signed_angle) < 1e-10, "Nearly identical vectors should have tiny angle")
}

test_nearly_opposite_vectors <- function() {
  epsilon <- 1e-15
  v1 <- c(1.0, 0.0)
  v2 <- c(-1.0, epsilon)
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_true(abs(abs(signed_angle) - PI) < 1e-10, "Nearly opposite vectors should be close to π")
}

test_mixed_positive_negative <- function() {
  v1 <- c(1.0, -2.0, 3.0)
  v2 <- c(-2.0, 1.0, -3.0)
  
  # Normalize
  v1 <- v1 / sqrt(sum(v1^2))
  v2 <- v2 / sqrt(sum(v2^2))
  
  signed_angle <- clock_hand_angle_between_vectors(v1, v2)
  assert_true(abs(signed_angle) >= 0 && abs(signed_angle) <= PI, "Mixed sign vectors in valid range")
}

# ==================== SHIFT VECTORS TESTS ====================

test_single_pair_shift_vectors <- function() {
  origins <- matrix(c(1.0, 0.0), nrow = 2, ncol = 1)
  targets <- matrix(c(0.0, 1.0), nrow = 2, ncol = 1)
  
  signed_angles <- clock_hand_angles_for_shift_vectors(origins, targets)
  assert_equal(abs(signed_angles[1]), PI/2, TOL, "Single pair shift vectors")
}

test_multiple_pairs_shift_vectors <- function() {
  # Three different rotations
  origins <- matrix(c(1.0, 0.0,   # 90° CCW
                     1.0, 0.0,   # 180°
                     1.0, 0.0),  # 90° CW
                   nrow = 2, ncol = 3)
  
  targets <- matrix(c(0.0, 1.0,   # 90° CCW
                     -1.0, 0.0,   # 180°
                     0.0, -1.0),  # 90° CW
                   nrow = 2, ncol = 3)
  
  signed_angles <- clock_hand_angles_for_shift_vectors(origins, targets)
  
  assert_equal(signed_angles[1], PI/2, TOL, "First rotation (90° CCW)")
  assert_equal(abs(signed_angles[2]), PI, TOL, "Second rotation (180°)")
  assert_equal(signed_angles[3], -PI/2, TOL, "Third rotation (90° CW)")
}

test_shift_vectors_with_selection_mask <- function() {
  # Four vectors, but only select 2nd and 4th
  origins <- matrix(c(1.0, 0.0,   # Not selected
                     1.0, 0.0,   # Selected (180°)
                     1.0, 0.0,   # Not selected
                     1.0, 0.0),  # Selected (45°)
                   nrow = 2, ncol = 4)
  
  targets <- matrix(c(0.0, 1.0,                              # Not selected
                     -1.0, 0.0,                             # Selected (180°)
                     0.0, -1.0,                             # Not selected
                     sqrt(2)/2, sqrt(2)/2),                 # Selected (45°)
                   nrow = 2, ncol = 4)
  
  vecs_selection_mask <- c(FALSE, TRUE, FALSE, TRUE)
  
  signed_angles <- clock_hand_angles_for_shift_vectors(origins, targets, vecs_selection_mask)
  
  assert_equal(abs(signed_angles[1]), PI, TOL, "Second vector (180°)")
  assert_equal(signed_angles[2], PI/4, TOL, "Fourth vector (45°)")
}

# ==================== CONSISTENCY TESTS ====================

test_consistency_between_functions <- function() {
  v1 <- c(1.0, 0.0)
  v2 <- c(0.0, 1.0)
  
  # Single function
  single_angle <- clock_hand_angle_between_vectors(v1, v2)
  
  # Batch function
  origins <- matrix(v1, nrow = 2, ncol = 1)
  targets <- matrix(v2, nrow = 2, ncol = 1)
  batch_angles <- clock_hand_angles_for_shift_vectors(origins, targets)
  
  assert_equal(single_angle, batch_angles[1], TOL, "Single vs batch consistency")
}

test_mathematical_properties <- function() {
  v1 <- c(1.0, 2.0)
  v2 <- c(3.0, 1.0)
  
  # Normalize
  v1 <- v1 / sqrt(sum(v1^2))
  v2 <- v2 / sqrt(sum(v2^2))
  
  angle_12 <- clock_hand_angle_between_vectors(v1, v2)
  angle_21 <- clock_hand_angle_between_vectors(v2, v1)
  
  assert_equal(angle_12, -angle_21, TOL, "Anti-commutativity of signed angles")
}

# ==================== TEST RUNNER ====================

#' Run all clock hand angle tests
run_all_clock_hand_angle_tests <- function() {
  cat("Running R clock hand angle tests...\n")
  cat("=====================================\n")
  
  test_functions <- list(
    # 2D tests
    test_identical_vectors_2d,
    test_opposite_vectors_2d,
    test_perpendicular_vectors_2d,
    test_45_degree_rotation_2d,
    test_clockwise_vs_counterclockwise_2d,
    
    # 3D tests
    test_identical_vectors_3d,
    test_perpendicular_vectors_3d,
    test_arbitrary_3d_rotation,
    
    # High dimensional tests
    test_high_dimensional_basic,
    test_high_dimensional_selected_axes,
    
    # Edge cases
    test_denormalized_vectors,
    test_tiny_vectors_precision,
    test_huge_vectors_precision,
    test_nearly_identical_vectors,
    test_nearly_opposite_vectors,
    test_mixed_positive_negative,
    
    # Shift vectors tests
    test_single_pair_shift_vectors,
    test_multiple_pairs_shift_vectors,
    test_shift_vectors_with_selection_mask,
    
    # Consistency tests
    test_consistency_between_functions,
    test_mathematical_properties
  )
  
  test_names <- c(
    "test_identical_vectors_2d",
    "test_opposite_vectors_2d", 
    "test_perpendicular_vectors_2d",
    "test_45_degree_rotation_2d",
    "test_clockwise_vs_counterclockwise_2d",
    "test_identical_vectors_3d",
    "test_perpendicular_vectors_3d",
    "test_arbitrary_3d_rotation",
    "test_high_dimensional_basic",
    "test_high_dimensional_selected_axes",
    "test_denormalized_vectors",
    "test_tiny_vectors_precision",
    "test_huge_vectors_precision",
    "test_nearly_identical_vectors",
    "test_nearly_opposite_vectors",
    "test_mixed_positive_negative",
    "test_single_pair_shift_vectors",
    "test_multiple_pairs_shift_vectors",
    "test_shift_vectors_with_selection_mask",
    "test_consistency_between_functions",
    "test_mathematical_properties"
  )
  
  passed <- 0
  failed <- 0
  
  for (i in seq_along(test_functions)) {
    cat(sprintf("\n--- %s ---\n", test_names[i]))
    tryCatch({
      test_functions[[i]]()
      passed <- passed + 1
    }, error = function(e) {
      cat(sprintf("✗ FAILED: %s\n", e$message))
      failed <- failed + 1
    })
  }
  
  cat("\n=====================================\n")
  cat(sprintf("Tests completed: %d passed, %d failed\n", passed, failed))
  
  if (failed == 0) {
    cat("🎉 All R clock hand angle tests passed!\n")
  } else {
    cat("❌ Some tests failed. Please check the output above.\n")
  }
  
  return(failed == 0)
}


# Run example if script is executed directly
if (!interactive()) {
  # Run tests when library is available
  run_all_clock_hand_angle_tests()
}
