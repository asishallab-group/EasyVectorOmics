# Load required libraries
require(circlize)
require(RColorBrewer)
require(latex2exp)

# Function to add alpha transparency to colors
addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

# --- Configuration ---
# Colors
colors <- brewer.pal(7, "Dark2")
alpha.cols <- sapply(colors, addAlpha)
vec_colors <- c(colors[[2]], colors[[3]], colors[[6]], colors[[7]])
vec_alpha.cols <- c(alpha.cols[[2]], alpha.cols[[3]], alpha.cols[[6]], alpha.cols[[7]])

# Output file
# pdf('timeSeriesTissueClock_updated_with_ASV.pdf', width = 5, height = 5)
png('figure_1e_with_ASV.png', width = 800, height = 800)

# --- Base Plot Setup ---
par(mar = c(1, 1, 1, 1)) # Set margins for the main plot
plot(NA, xlim = c(0, 2.5), ylim = c(0, 2.5), axes = FALSE, xlab = '', ylab = '', asp = 1)

# Legend
legend(x = 0, y = 2.5,
       legend = c(TeX('$Viewing \\ angle \\ \\alpha = 30 \\degree :$'),
                  TeX('$deviation\\ from\\ diagonal\\ perspective$')),
       bty = "n", cex = 1.2, text.col = c(colors[4], colors[4]))

# --- Plot Elements ---
# Point O (Origin)
O <- c(1, 1)
points(O[1], O[2], bg = alpha.cols[1], col = colors[1], pch = 21, cex = 3)
text(O[1], O[2], label = TeX('$\\textit{O}$'), pos = 2, offset = 1, col = colors[1], cex = 2)

# Point D (Space origin)
D <- c(1.7, 0.5)
points(D[1], D[2], bg = alpha.cols[4], col = colors[4], pch = 21, cex = 2)
text(x = 1.5, y = 0.5, label = TeX('$\\alpha = 30^\\circ$'), cex = 1.3, col=colors[4])
text(x = 2.1, y = 0.5, label = TeX('Space origin and axes'), cex = 1)

# Anchor axis (AA)
lines(c(O[1], O[1]), c(O[2], O[2] + 1), col = "#bdbdbd", lwd = 2)
text(O[1]+0.3, O[2]+1.1, label = "Anchor Axis (AA), projected tissue X", adj = 1, col = "gray40", cex = 1)

# New axes at D
axis_length <- 0.2
lines(c(D[1], D[1]), c(D[2], D[2] + axis_length), col = "black", lwd = 1.3)
lines(c(D[1], D[1] - axis_length), c(D[2], D[2] - axis_length), col = "black", lwd = 1.3)
lines(c(D[1], D[1] + axis_length), c(D[2], D[2] - axis_length), col = "black", lwd = 1.3)
text(D[1], D[2] + axis_length + 0.01, label = TeX('\\textit{tissue X}'), pos = 3, cex = 0.7)
text(D[1]-0.2, D[2] + axis_length - .47, label = TeX('\\textit{tissue Y}'), pos = 3, cex = 0.7)
text(D[1]+0.25 , D[2] + axis_length - .47, label = TeX('\\textit{tissue Z}'), pos = 3, cex = 0.7)

# Dashed vector D
segments(O[1], O[2], D[1], D[2], col = colors[4], lwd = 1.5, lty = 2)

# Clock vectors t1–t4
clock_pts <- list(
  t1 = O + c(0.2, 0.65),
  t2 = O + c(0.65, 0.1),
  t3 = O + c(0.2, -0.6),
  t4 = O + c(-0.5, -0.3)
)

# Draw vectors and labels
for (i in 1:4) {
  pt <- clock_pts[[i]]
  arrows(O[1], O[2], pt[1], pt[2], col = vec_colors[i], lwd = 4, length = 0.13)
  points(pt[1], pt[2], pch = 21, bg = vec_alpha.cols[i], col = vec_colors[i], cex = 3)
  # *** MODIFIED: Added 'clock' superscript ***
  text(pt[1]+.02, pt[2]+.02, label = TeX(paste0('$t_', i, '^{{clock}}$')), pos = 3, col = vec_colors[i], cex = 1.5)
  text(pt[1]+.02, pt[2] + .1, label = TeX('$\\rightarrow$'), pos = 3, col = vec_colors[i], cex = 1.5)
}

# Phi arcs
phi_label_positions <- list(
  c(1.05, 1.33),  # phi_1
  c(1.25, 1.2),  # phi_2
  c(1.2, 0.7),  # phi_3
  c(0.8, 0.66)   # phi_4
)
anchor_vec_phi <- c(0, 1) # Reference for angle calculation (vertical Anchor Axis)
for (i in 1:4) {
  vec <- clock_pts[[i]] - O
  draw.sector(
    start.degree = atan2(anchor_vec_phi[2], anchor_vec_phi[1]) * 180 / pi,
    end.degree = atan2(vec[2], vec[1]) * 180 / pi,
    clock.wise = TRUE,
    center = O,
    border = vec_colors[i],
    col = NA,
    lwd = 2,
    rou1 = 0.35 + i * 0.03
  )
  label_pos <- phi_label_positions[[i]]
  text(label_pos[1], label_pos[2],
       label = TeX(paste0('$\\varphi_', i, '$')),
       col = vec_colors[i],
       cex = 1.8)
}

# --- *** NEW: Calculate Angles and ASV *** ---
# Calculate angles relative to the vertical anchor axis (0, 1) pointing up from O
anchor_vec_asv <- c(0, 1)
angles_rad <- sapply(clock_pts, function(pt) {
  vec <- pt - O
  # Calculate angle relative to the anchor vector (positive y-axis from O)
  angle_rel_anchor <- atan2(vec[2], vec[1]) - atan2(anchor_vec_asv[2], anchor_vec_asv[1])
  # Normalize angle to be within (-pi, pi]
  angle_rel_anchor <- atan2(sin(angle_rel_anchor), cos(angle_rel_anchor))
  return(angle_rel_anchor)
})

# Ensure angles are in a consistent range, e.g., [0, 2*pi) for easier diff calc
# angles_rad <- ifelse(angles_rad < 0, angles_rad + 2 * pi, angles_rad) # Optional: use if angle_diff needs [0, 2pi)

# Helper function to calculate shortest angle difference (handles wrap-around)
angle_diff <- function(a2, a1) {
  diff <- a2 - a1
  # Normalize to (-pi, pi] range
  while (diff <= -pi) diff <- diff + 2 * pi
  while (diff > pi) diff <- diff - 2 * pi
  return(diff)
}

# Calculate angular differences (ASV) in radians
# ASV(t_i) = angle(t_i) - angle(t_{i-1})
# Start with 0 for the first time point
asv_rad <- c(0, sapply(2:length(angles_rad), function(i) angle_diff(angles_rad[i], angles_rad[i-1])))

# Convert ASV to degrees for plotting
asv_deg <- asv_rad * 180 / pi

# --- Inset Plots ---

# Reset graphics parameters for insets (margins, position relative to axis)
par(mgp = c(2, 0.5, 0)) # control axis title, label, line positions
par(mar = c(3, 3, 1, 1)) # Smaller margins for insets: c(bottom, left, top, right)

# --- Inset 1: distance vs time ---
# Define plotting region: c(xmin, xmax, ymin, ymax) in normalized device coordinates (0-1)
par(fig = c(0.6, 0.98, 0.75, 0.98), new = TRUE) # Position: top-right

distances <- sapply(clock_pts, function(pt) sqrt(sum((pt - O)^2)))

plot(
  1:4,
  distances,
  type = "l",
  xaxt = "n",       # Turn off x-axis labels (will be added on bottom plot)
  ylab = "Shift",
  xlab = "",        # Turn off x-axis title
  main = "",        # Turn off main title
  col = "black",
  cex.axis = 0.8,
  cex.lab = 0.9     # Slightly smaller label size for insets
)
# No axis(1, ...) here
points(1:4, distances, pch = 16, col = vec_colors, cex = 1.5) # Smaller points for inset
box()

# --- *** NEW: Inset 2: ASV vs time *** ---
# Define plotting region below the first inset, aligning left/right edges
par(fig = c(0.6, 0.98, 0.52, 0.75), new = TRUE) # Position below first inset

plot(
  1:4,
  asv_deg,          # Use degrees for the y-axis
  type = "l",
  xaxt = "n",       # Turn off default x-axis
  ylab = "ASV (deg)", # Y-axis label
  xlab = "Time",    # X-axis label (only on the bottom plot)
  main = "",
  col = "black",
  cex.axis = 0.8,
  cex.lab = 0.9
)
# Add custom x-axis labels
axis(1, at = 1:4, labels = paste0("t", 1:4), cex.axis = 0.8)
asv_point_colors <- c("grey50", vec_colors[2:4]) # t2-t1 uses t2 color, t3-t2 uses t3 color, etc.
points(1:4, asv_deg, pch = 16, col = asv_point_colors, cex = 1.5)
box()

# --- Finalize ---
# Close the graphics device
dev.off()
