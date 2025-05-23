require(RColorBrewer)
colors <- brewer.pal(7, "Dark2")
colors2 <- brewer.pal(9, "Set1")

pdf("figure_1g.pdf", width = 12, height = 12)
# par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))  # 4 panels (2x2) with minimal margins
layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(1.1, 1.3))
par(mar = c(1, 1, 1, 1))

# Function to draw the green amoeba-like contour
draw_contour <- function() {
  a <- 1.9  # length along E1 axis
  b <- 1.5  # width along E2 axis
  theta <- -2 * pi / 3  # rotation angle in radians (approx. -120°)

  t <- seq(0, 2 * pi, length.out = 200)
  r <- 1 + 0.15 * sin(3 * t)  # wavy radius for amoeba effect
  x <- a * r * cos(t)
  y <- b * r * sin(t)

  # Apply rotation
  rotation_matrix <- matrix(c(cos(theta), -sin(theta),
                              sin(theta), cos(theta)), ncol = 2)
  rotated <- rotation_matrix %*% rbind(x, y)

  # Plot contour
  lines(rotated[1, ] - 0.3, rotated[2, ], col = "#66A61E", lwd = 3)
  text(0, 1.8, "Photosynthesis", cex = 2.5, col = "#66A61E")
}

# Function to draw E1 and E2 axes with arrows and labels
draw_axis <- function() {
  theta <- -pi / 6  # 30 degrees in radians
  rotation_matrix <- matrix(c(cos(theta), -sin(theta),
                              sin(theta), cos(theta)), ncol = 2)

  e1 <- rotation_matrix %*% c(0, 2) * 0.8
  e2 <- rotation_matrix %*% c(2, 0) * 0.65

  # Arrows for E1 in both directions
  arrows(-e1[1] - 0.6, -e1[2] - 0.2, e1[1] - 0.6, e1[2] - 0.2, col = "black", lwd = 2, length = 0.2)
  arrows(e1[1] - 0.6, e1[2] - 0.2, -e1[1] - 0.6, -e1[2] - 0.2, col = "black", lwd = 2, length = 0.2)

  # Arrows for E2 in both directions
  arrows(-e2[1], -e2[2], e2[1], e2[2], col = "black", lwd = 2, length = 0.2)
  arrows(e2[1], e2[2], -e2[1], -e2[2], col = "black", lwd = 2, length = 0.2)

  # Axis labels
  text(0.4, -1.9, labels = expression(E[1]), pos = 3, cex = 2.3, srt = 30)
  text(0.9, 0.15, labels = expression(E[2]), pos = 4, cex = 2.3, srt = 30)
}

# Function to draw a single vector
draw_vector <- function(x0, y0, length, angle_deg, col = "black", lwd = 4, arrow_size = 0.1, lty = "solid") {
  angle_rad <- angle_deg * pi / 180
  x1 <- x0 + length * cos(angle_rad)
  y1 <- y0 + length * sin(angle_rad)
  arrows(x0, y0, x1, y1, length = arrow_size, col = "black", lwd = lwd, lty = lty)
}

# Draws a vector and its projections onto the E1 and E2 axes
draw_vector_with_projections <- function(x0, y0, length, angle_deg,
                                         e1_angle = 30, e2_angle = 120,
                                         e1_offset = c(0, 0),
                                         e2_offset = c(-0.6, -0.2),
                                         col = "gray60", lwd = 1.5, arrow_size = 0.1,
                                         proj_col = "black", proj_lty = "dotted") {
  angle_rad <- angle_deg * pi / 180
  e1_rad <- e1_angle * pi / 180
  e2_rad <- e2_angle * pi / 180

  vx <- length * cos(angle_rad)
  vy <- length * sin(angle_rad)
  x_tip <- x0 + vx
  y_tip <- y0 + vy

  e1 <- c(cos(e1_rad), sin(e1_rad))
  e2 <- c(cos(e2_rad), sin(e2_rad))

  project_to_axis <- function(xt, yt, base_dir, offset) {
    v <- c(xt, yt) - offset
    proj_len <- sum(v * base_dir)
    proj_point <- offset + proj_len * base_dir
    return(proj_point)
  }

  proj_e1_pt <- project_to_axis(x_tip, y_tip, e1, e1_offset)
  proj_e2_pt <- project_to_axis(x_tip, y_tip, e2, e2_offset)

  arrows(x0, y0, x_tip, y_tip, length = arrow_size, col = col, lwd = lwd)
  segments(x_tip, y_tip, proj_e1_pt[1], proj_e1_pt[2], col = proj_col, lty = proj_lty, lwd = 2)
  segments(x_tip, y_tip, proj_e2_pt[1], proj_e2_pt[2], col = proj_col, lty = proj_lty, lwd = 2)
  points(proj_e1_pt[1], proj_e1_pt[2], col = colors[[4]], pch = 16, cex = 1.5)
  points(proj_e2_pt[1], proj_e2_pt[2], col = colors[[4]], pch = 16, cex = 1.5)
}

# Draws a spline-smoothed semantic axis
draw_semantic_axis <- function() {
  x_points <- c(-.6, -0.65, -0.65, -0.5, -0.45, -0.4, 0)
  y_points <- c(-1.5, -1.1, -0.7, -.2, 0.5, 1, 1.4)
  t <- seq(0, 1, length.out = length(x_points))
  spline_x <- spline(t, x_points, n = 200)
  spline_y <- spline(t, y_points, n = 200)
  lines(spline_x$y, spline_y$y, col = colors[[4]], lwd = 5)
}

# Define parameters for 10 vectors to be drawn
x_vec <- c(-0.7, -1.2, -0.7,  0.0, -0.3, -0.8,  0.0, -0.15, -0.4,  0.58)
y_vec <- c( 0.9,  0.4, -0.3,  1.1,  0.5, -0.6, -0.8, -1.10, -0.9, -0.25)
len_vec <- c(0.25, 0.3, 0.2, 0.3, 0.4, 0.2, 0.4,  0.3,  0.1,  0.4)
angle_vec <- c(150, 170, 200,  50,  80, 220, 350, 270, 350,  70)

# --- Top-left Panel: Vector field with E1 and E2, with vector projections ---
plot(NA, xlim = c(-1.95, 1.95), ylim = c(-1.95, 1.95), asp = 1, axes = FALSE, xlab = "", ylab = "")

# Draw amoeba-like contour representing the biological response area
draw_contour()

# Draw principal axes E1 and E2 (semantic structure)
draw_axis()

# Draw all directional vectors with dotted projections onto E1 and E2
for (i in seq_along(x_vec)) {
  draw_vector_with_projections(x_vec[i], y_vec[i], len_vec[i], angle_vec[i])
}

# Draw the smoothed semantic axis that runs through the space
draw_semantic_axis()


# --- Bottom-left Panel: Same field, now highlighting semantic flow ---
plot(NA, xlim = c(-2, 2), ylim = c(-2, 2), asp = 1, axes = FALSE, xlab = "", ylab = "")

# Draw contour and axes again
draw_contour()
draw_axis()

# Highlight semantic axis direction with thick arrow
arrows(-0.35, -.1, -0.43, -.7, col = colors[[4]], lwd = 12, lty = 1, length = 0.2)

# --- Thick curved arrow (semantic axis trajectory) ---
# Step 1: Generate a symmetric bell-shaped curve (normal distribution)
curve_y <- seq(-1, 1, length.out = 100)
curve_x <- dnorm(curve_y, mean = 0, sd = 0.4)

# Step 2: Scale down to make the shape smaller
curve_x <- curve_x * 0.6
curve_y <- curve_y * 0.6
curve_x <- curve_x - mean(curve_x)  # Center the curve horizontally

# Step 3: Rotate the shape slightly to align with the semantic axis
theta <- pi / 20  # Small positive angle (≈9°)
rotation_matrix <- matrix(c(cos(theta), -sin(theta),
                            sin(theta),  cos(theta)), ncol = 2)
rotated_coords <- rotation_matrix %*% rbind(curve_x, curve_y)

# Step 4: Plot the rotated curve as a thick semantic arrow
lines(rotated_coords[1, ] - 0.4, rotated_coords[2, ] + 0.5, lwd = 12, col = colors[[4]])

# Draw the original vectors again, this time without projections
for (i in seq_along(x_vec)) {
  draw_vector(x_vec[i], y_vec[i], len_vec[i], angle_vec[i])
}

# Re-draw semantic axis curve to connect vector tips smoothly
draw_semantic_axis()

# --- Top-right Panel: Signal curves projected onto the semantic axis ---

# Create an empty plot with specified axis limits and no axis lines
plot(NA, xlim = c(0, 10), ylim = c(0, 1.2), type = "n",
     xlab = "", ylab = "", axes = FALSE, asp = 1)

# Draw axis arrows
arrows(0, -3, 8, -3, length = 0.1, lwd = 1.5)  # Horizontal axis: Semantic Axis (E1)
arrows(0, -3, 0, 4, length = 0.1, lwd = 1.5)   # Vertical axis: Signal Strength

# Add axis labels
text(7, -3.5, "Semantic Axis", cex = 2.5)
text(2, 4.5, "Centered Value", cex = 2.5)

# --- First Curve: Directional Spread (Dir-Spread) ---

# Generate a sequence of x values
x <- seq(1, 8, length.out = 200)

# Compute a right-skewed normal distribution
y1 <- dnorm(x, mean = 1.6, sd = 1.3)

# Normalize height and scale for display
y1 <- y1 / max(y1) * 3

# Label the curve
text(4.5, 2.5, "Dir-Spread", cex = 2.3, col = colors2[[1]])

# Draw the curve
lines(x, y1, lwd = 4, col = colors2[[1]])

# --- Second Curve: Significance Strength (Sign-Strn) ---

# Refine x-resolution for smoother drawing
x <- seq(1, 8, length.out = 500)

# Combine multiple Gaussian peaks (positive and negative) to simulate signal strength dynamics
y <- dnorm(x, mean = 1.8, sd = 1.2) - 
     0.5 * dnorm(x, mean = 4, sd = 1.5) + 
     0.4 * dnorm(x, mean = 6.5, sd = 0.7)

# Normalize height to max value and rescale
y <- y / max(y) * 2.7

# Draw the composite curve
lines(x, y, lwd = 4, col = colors2[[2]])

# Label the second curve
text(5, -1.5, "Sign-Strn", cex = 2.3, col = colors2[[2]])



dev.off()
