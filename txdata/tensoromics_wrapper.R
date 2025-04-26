# Load libraries
dyn.load("libtensoromics_l.so")
dyn.load("libtensoromics_g.so")

# Wrapper functions

# Initialize the TensorOmics object (LFortran side)
init_tensoromics <- function(data_estimates) {
  n_cond <- length(data_estimates)
  n_genes <- sum(data_estimates)
  
  # Call LFortran function to initialize the TensorOmics object
  tom <- .C("init_tensoromics",
            as.integer(n_cond),
            as.integer(n_genes),
            PACKAGE = "tensoromics_l")
  
  # Return the initialized object in a list
  list(n_conditions = n_cond, n_genes = n_genes, tom = tom)
}

# Calculate memory requirements (LFortran side)
# Load libraries
dyn.load("./libtensoromics_l.so")
dyn.load("./libtensoromics_g.so")

# Wrapper functions
init_tensoromics <- function(data_estimates) {
  n_cond <- length(data_estimates)
  n_genes <- sum(data_estimates)
  
  handle <- .C("init_tensoromics",
               as.integer(n_cond),
               as.integer(n_genes),
               handle = raw(8),  # 64-bit pointer
               PACKAGE = "tensoromics_l")$handle
  
  list(handle = handle, n_conditions = n_cond, n_genes = n_genes)
}

calculate_memory <- function(tom) {
  .C("calculate_memory_requirements",
     as.integer(tom$n_conditions),
     as.integer(tom$n_genes),
     PACKAGE = "tensoromics_l")
}

update_tensoromics <- function(tom, patch_matrix) {
  indices <- integer(ncol(patch_matrix))
  flattened_patch <- as.single(c(patch_matrix))
  
  .C("update_tensoromics",
     tom$handle,
     flattened_patch,
     as.integer(ncol(patch_matrix)),
     indices = indices,
     PACKAGE = "tensoromics_l")$indices
}

save_tensoromics <- function(tom, filename) {
  .C("transfer_for_save",
     tom$handle,
     PACKAGE = "tensoromics_l")
  
  .C("save_tensoromics",
     as.character(filename),
     PACKAGE = "tensoromics_g")
}


# Test 1: Initialize TensorOmics object
test_init_tensoromics <- function() {
  # Test data estimates (n_conditions = 3, n_genes = 5)
  data_estimates <- c(3, 3, 3, 3, 3)  # 5 genes, each with 3 conditions
  
  tom <- init_tensoromics(data_estimates)
  
  # Check the dimensions
  cat("Test 1: Initialize TensorOmics\n")
  cat("Expected n_conditions: 3, n_genes: 15\n")
  cat("Returned n_conditions: ", tom$n_conditions, "\n")
  cat("Returned n_genes: ", tom$n_genes, "\n")
  stopifnot(tom$n_conditions == 3)
  stopifnot(tom$n_genes == 15)
}

# Test 2: Calculate memory requirements
test_calculate_memory <- function() {
  tom <- list(n_conditions = 3, n_genes = 15)
  
  # Call the memory calculation function via wrapper
  cat("Test 2: Calculate Memory Requirements\n")
  calculate_memory(tom)
}

# Test 3: Update TensorOmics object
test_update_tensoromics <- function() {
  # Initialize TensorOmics object
  data_estimates <- c(3, 3, 3, 3, 3)  # 5 genes, each with 3 conditions
  tom <- init_tensoromics(data_estimates)
  
  # Create a patch matrix (2 patches for 3 conditions)
  patch_matrix <- matrix(runif(3 * 2), nrow = 3, ncol = 2)
  
  # Call the update function via wrapper
  updated_indices <- update_tensoromics(tom, patch_matrix)
  
  cat("Test 3: Update TensorOmics\n")
  cat("Updated indices: ", updated_indices, "\n")
  stopifnot(length(updated_indices) == 2)
}

# Test 4: Save TensorOmics object
test_save_tensoromics <- function() {
  # Initialize TensorOmics object
  data_estimates <- c(3, 3, 3, 3, 3)  # 5 genes, each with 3 conditions
  tom <- init_tensoromics(data_estimates)
  
  # Call the save function via wrapper (using GFortran save)
  filename <- "test_tensoromics_data.bin"
  save_tensoromics(tom, filename)
  
  # Verify if the file is created (just checking existence in this case)
  cat("Test 4: Save TensorOmics\n")
  if (file.exists(filename)) {
    cat("File successfully saved: ", filename, "\n")
    file.remove(filename)  # Clean up the saved file
  } else {
    stop("Failed to save tensoromics data.")
  }
}

# Run the tests
cat("Running all tests...\n")
test_init_tensoromics()
test_calculate_memory()
test_update_tensoromics()
test_save_tensoromics()

cat("All tests passed!\n")
