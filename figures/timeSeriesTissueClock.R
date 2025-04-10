require(circlize)
require(RColorBrewer)
require(latex2exp)

# Helper function
addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

# Colors
colors <- brewer.pal(7, "Dark2")
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

legend(
  x = -.2,
  y = 2.63,
  legend = c(TeX('$Viewing \\ angle \\ \\alpha = 30 \\degree :$'),
    TeX('$deviation\\ from\\ diagonal\\ perspective$')),
  bty = "n",
  cex = 0.7
)

# Center point O
O <- c(1, 1)
points(O[1], O[2], bg = alpha.cols[1], col = colors[1], pch = 21, cex = 2)
text(O[1], O[2], label = TeX('$\\textit{O}$'), pos = 2, offset = 1, col = colors[1], cex = 1.2)
text(O[1] + .04, O[2] + .08, label = TeX('$\\rightarrow$'), pos = 2, offset = 1, col = colors[1], cex = 1.2)

# Point D
D <- c(1.7, 0.5)
points(D[1], D[2], bg = alpha.cols[4], col = colors[4], pch = 21, cex = 2)
text(D[1], D[2], label = TeX('$\\textit{D}$'), pos = 4, offset = 1, col = colors[4], cex = 1.2)
text(D[1] - .04, D[2] + .08, label = TeX('$\\rightarrow$'), pos = 4, offset = 1, col = colors[4], cex = 1.2)


# Triangle for tissue space
lines(c(1, 1), c(1, (1 + sqrt(2))), col = "black", lwd = 1.5)  # to tissue X
lines(c(1, 2), c(1, 0), col = "black", lwd = 1.5)              # to tissue Z
lines(c(1, 0), c(1, 0), col = "black", lwd = 1.5)              # to tissue Y

# Tissue labels
text(1, 2.2, label = TeX('\\textit{tissue} X anchor'), pos = 4, offset = .2, cex = 1.0)
text(0.4, 0, label = TeX('\\textit{tissue} Y'), pos = 3, offset = 1, cex = 1.0)
text(2.1, 0, label = TeX('\\textit{tissue} Z'), pos = 3, offset = 1.25, cex = 1.0)

# Define 4 time points relative to O
clock_pts <- list(
  t1 = O + c(0.2, 0.65),
  t2 = O + c(0.65, 0.1),
  t3 = O + c(0.2, -0.6),
  t4 = O + c(-0.5, -0.3)
)

# Vector from O to t1 and t2
anch <- c(0, 1)
v1 <- clock_pts$t1 - O
v2 <- clock_pts$t2 - O

# Angle in degrees between v1 and v2
angle_rad <- acos(sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2))))
angle_deg <- angle_rad * (180 / pi)

# Draw arc representing phi
draw.sector(
  start.degree = atan2(anch[2], anch[1]) * 180 / pi,
  end.degree   = atan2(v1[2], v1[1]) * 180 / pi,
  clock.wise = TRUE,
  center = O,
  border = colors[[3]],
  col = alpha.cols[[3]],
  lwd = 2,
  rou1 = 0.4
)

draw.sector(
  start.degree = atan2(anch[2], anch[1]) * 180 / pi,
  end.degree   = atan2(v2[2], v2[1]) * 180 / pi,
  clock.wise = TRUE,
  center = O,
  border = colors[[3]],
  col = alpha.cols[[3]],
  lwd = 2,
  rou1 = 0.45
)

v3 <- clock_pts$t3 - O
draw.sector(
  start.degree = atan2(anch[2], anch[1]) * 180 / pi,
  end.degree   = atan2(v3[2], v3[1]) * 180 / pi,
  clock.wise = TRUE,
  center = O,
  border = colors[[3]],
  col = alpha.cols[[3]],
  lwd = 2,
  rou1 = 0.5
)

v4 <- clock_pts$t4 - O
draw.sector(
  start.degree = atan2(anch[2], anch[1]) * 180 / pi,
  end.degree   = atan2(v4[2], v4[1]) * 180 / pi,
  clock.wise = TRUE,
  center = O,
  border = colors[[3]],
  col = alpha.cols[[3]],
  lwd = 2,
  rou1 = 0.55
)

# Place φ label approximately between v1 and v2
mid_angle <- atan2(v1[2] + v2[2], v1[1] + v2[1])
label_pos <- O + 0.4 * c(cos(mid_angle), sin(mid_angle))
text(label_pos[1] + .2, label_pos[2] + .2,
     label = TeX('$\\varphi$'),
     col = colors[[3]],
     cex = 2)


# Plot time point vectors from O
vec_colors <- c(
  colors[[2]],
  colors[[5]],
  colors[[6]],
  colors[[7]]
)

vec_alpha.cols <- c(
  alpha.cols[[2]],
  alpha.cols[[5]],
  alpha.cols[[6]],
  alpha.cols[[7]]
)

# Draw vectors
i <- 1
for (tp in names(clock_pts)) {
  pt <- clock_pts[[tp]]
  arrows(O[1], O[2], pt[1], pt[2],
         col = vec_colors[i],
         lwd = 2,
         length = 0.13)
  points(pt[1], pt[2], pch = 21, bg = vec_alpha.cols[[i]], col = vec_colors[i], cex = 2)
  i <- i + 1
}
# Labels for vectors
text(clock_pts[[1]][1], clock_pts[[1]][2], label = TeX(paste0('$t_', 1, '$')), pos = 3, col = vec_colors[1], cex = 1.2)
text(clock_pts[[2]][1], clock_pts[[2]][2], label = TeX(paste0('$t_', 2, '$')), pos = 3, col = vec_colors[2], cex = 1.2)
text(clock_pts[[3]][1], clock_pts[[3]][2], label = TeX(paste0('$t_', 3, '$')), pos = 4, col = vec_colors[3], cex = 1.2)
text(clock_pts[[4]][1], clock_pts[[4]][2], label = TeX(paste0('$t_', 4, '$')), pos = 2, col = vec_colors[4], cex = 1.2)

text(clock_pts[[1]][1], clock_pts[[1]][2] + .1, label = TeX('$\\rightarrow$'), pos = 3, col = vec_colors[1], cex = 1.2)
text(clock_pts[[2]][1], clock_pts[[2]][2] + .1, label = TeX('$\\rightarrow$'), pos = 3, col = vec_colors[2], cex = 1.2)
text(clock_pts[[3]][1] - .04, clock_pts[[3]][2] + .09, label = TeX('$\\rightarrow$'), pos = 4, col = vec_colors[3], cex = 1.2)
text(clock_pts[[4]][1] + .05, clock_pts[[4]][2] + .09, label = TeX('$\\rightarrow$'), pos = 2, col = vec_colors[4], cex = 1.2)

# Vector O to D
arrows(
  O[1], O[2], D[1], D[2],
  col = colors[4],
  lwd = 2,
  length = 0.13
)

# Mini plot in top right (distance vs time)
par(fig = c(0.6, 1, 0.65, 1), new = TRUE)
distances <- sapply(clock_pts, function(pt) sqrt(sum((pt - O)^2)))
par(mgp = c(1, 0.5, 0))
par(mar = c(4,4,1,1))

plot(
  1:4,
  distances,
  type = "l",
  pch = 16,
  xaxt = "n",
  ylab = "Shift",
  xlab = "Time",
  main = "",
  col = "black",
  cex.axis = 0.5,
  cex.lab = 0.6
)
axis(1, at = 1:4, labels = paste0("t", 1:4), cex.axis = 0.5)
points(1:4, distances, pch = 16, cex = .7, col = "black")
box()

dev.off()
