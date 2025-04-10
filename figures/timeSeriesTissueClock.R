require(circlize)
require(RColorBrewer)
require(latex2exp)

# Helper function
addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

# Colors
colors <- brewer.pal(5, "Dark2")
alpha.cols <- sapply(colors, addAlpha)

# Output PDF
pdf('timeSeriesTissueClock.pdf', width = 5, height = 5)

# Setup main plot
par(mar = c(1, 1, 1, 1))
plot(
  NA,
  xlim = c(0, 2.5),
  ylim = c(0, 2.5),
  axes = FALSE,
  xlab = '',
  ylab = '',
  asp = 1
)

# Center point O
O <- c(1, 1)
points(O[1], O[2], col = "red", pch = 20, cex = 3)
text(O[1], O[2], label = TeX('$\\bar{O}$'), pos = 4, offset = 1.25, col = "red", cex = 1.2)

# Triangle for tissue space
lines(c(1, 1), c(1, (1 + sqrt(2))), col = "black", lwd = 1.5)  # to tissue X
lines(c(1, 2), c(1, 0), col = "black", lwd = 1.5)              # to tissue Z
lines(c(1, 0), c(1, 0), col = "black", lwd = 1.5)              # to tissue Y

# Tissue labels
text(1, 2.2, label = TeX('\\textit{tissue} X'), pos = 4, offset = .2, cex = 1.0)
text(0.5, 0, label = TeX('\\textit{tissue} Y'), pos = 3, offset = .1, cex = 1.0)
text(2, 0, label = TeX('\\textit{tissue} Z'), pos = 3, offset = 1.25, cex = 1.0)

# Define 4 time points relative to O
clock_pts <- list(
  t1 = O + c(0.0, 0.75),
  t2 = O + c(0.6, 0.3),
  t3 = O + c(0.4, -0.4),
  t4 = O + c(-0.3, -0.5)
)

# Plot time point vectors from O
i <- 1
for (tp in names(clock_pts)) {
  pt <- clock_pts[[tp]]
  arrows(O[1], O[2], pt[1], pt[2],
         col = colors[i],
         lwd = 3,
         length = 0.1)
  points(pt[1], pt[2], pch = 21, bg = alpha.cols[[i]], col = colors[i], cex = 2)
  text(pt[1], pt[2], label = TeX(paste0('$t_', i, '$')), pos = 3, col = colors[i], cex = 1.2)
  i <- i + 1
}
# Vector from O to t1 and t2
v1 <- clock_pts$t1 - O
v2 <- clock_pts$t2 - O

# Angle in degrees between v1 and v2
angle_rad <- acos(sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2))))
angle_deg <- angle_rad * (180 / pi)

# Draw arc representing phi
draw.sector(
  start.degree = atan2(v1[2], v1[1]) * 180 / pi,
  end.degree   = atan2(v2[2], v2[1]) * 180 / pi,
  clock.wise = TRUE,
  center = O,
  border = colors[[5]],
  col = alpha.cols[[5]],
  lwd = 2,
  rou1 = 0.35
)

# Place φ label approximately between v1 and v2
mid_angle <- atan2(v1[2] + v2[2], v1[1] + v2[1])
label_pos <- O + 0.4 * c(cos(mid_angle), sin(mid_angle))
text(label_pos[1], label_pos[2],
     label = TeX('$\\varphi$'),
     col = colors[[5]],
     cex = 2)


# Mini plot in top right (distance vs time)
par(fig = c(0.65, 0.98, 0.70, 0.98), new = TRUE)
distances <- sapply(clock_pts, function(pt) sqrt(sum((pt - O)^2)))
plot(
  1:4,
  distances,
  type = "b",
  pch = 16,
  xaxt = "n",
  ylab = "Distance",
  xlab = "Time",
  main = "",
  col = "black",
  cex.axis = 0.75,
  cex.lab = 0.8
)
axis(1, at = 1:4, labels = paste0("t", 1:4), cex.axis = 0.75)
box()

dev.off()
