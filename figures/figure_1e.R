require(circlize)
require(RColorBrewer)
require(latex2exp)

addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

# Colors
colors <- brewer.pal(7, "Dark2")
alpha.cols <- sapply(colors, addAlpha)

# Output
# pdf('timeSeriesTissueClock_updated.pdf', width = 5, height = 5)
png('figure_1e.png', width = 800, height = 800)

par(mar = c(1, 1, 1, 1))
plot(NA, xlim = c(0, 2.5), ylim = c(0, 2.5), axes = FALSE, xlab = '', ylab = '', asp = 1)

legend(x = 0, y = 2.5,
       legend = c(TeX('$Viewing \\ angle \\ \\alpha = 30 \\degree :$'),
                  TeX('$deviation\\ from\\ diagonal\\ perspective$')),
       bty = "n", cex = 1.2, text.col = c(colors[4], colors[4]))

# === Point O ===
O <- c(1, 1)
points(O[1], O[2], bg = alpha.cols[1], col = colors[1], pch = 21, cex = 3)
text(O[1], O[2], label = TeX('$\\textit{O}$'), pos = 2, offset = 1, col = colors[1], cex = 2)
# text(O[1] + .04, O[2] + .08, label = TeX('$\\rightarrow$'), pos = 2, offset = 1, col = colors[1], cex = 2)

# === Point D ===
D <- c(1.7, 0.5)
points(D[1], D[2], bg = alpha.cols[4], col = colors[4], pch = 21, cex = 2)

text(x = 1.5, y = 0.5, label = TeX('$\\alpha = 30^\\circ$'), cex = 1.3, col=colors[4])
text(x = 2.1, y = 0.5, label = TeX('Space origin and axes'), cex = 1)

# === Anchor axis (AA) ===
lines(c(O[1], O[1]), c(O[2], O[2] + 1), col = "#bdbdbd", lwd = 2)
text(O[1]+0.3, O[2]+1.1, label = "Anchor Axis (AA), projected tissue X", adj = 1, col = "gray40", cex = 1)

# === New axes at D ===
axis_length <- 0.2
lines(c(D[1], D[1]), c(D[2], D[2] + axis_length), col = "black", lwd = 1.3)
lines(c(D[1], D[1] - axis_length), c(D[2], D[2] - axis_length), col = "black", lwd = 1.3)
lines(c(D[1], D[1] + axis_length), c(D[2], D[2] - axis_length), col = "black", lwd = 1.3)

text(D[1], D[2] + axis_length + 0.01, label = TeX('\\textit{tissue X}'), pos = 3, cex = 0.7)
text(D[1]-0.2, D[2] + axis_length - .47, label = TeX('\\textit{tissue Y}'), pos = 3, cex = 0.7)
text(D[1]+0.25 , D[2] + axis_length - .47, label = TeX('\\textit{tissue Z}'), pos = 3, cex = 0.7)

# === Dashed vector D ===
segments(O[1], O[2], D[1], D[2], col = colors[4], lwd = 1.5, lty = 2)

# === Clock vectors t1–t4 ===
clock_pts <- list(
  t1 = O + c(0.2, 0.65),
  t2 = O + c(0.65, 0.1),
  t3 = O + c(0.2, -0.6),
  t4 = O + c(-0.5, -0.3)
)

vec_colors <- c(colors[[2]], colors[[3]], colors[[6]], colors[[7]])
vec_alpha.cols <- c(alpha.cols[[2]], alpha.cols[[3]], alpha.cols[[6]], alpha.cols[[7]])

for (i in 1:4) {
  pt <- clock_pts[[i]]
  arrows(O[1], O[2], pt[1], pt[2], col = vec_colors[i], lwd =4, length = 0.13)
  points(pt[1], pt[2], pch = 21, bg = vec_alpha.cols[i], col = vec_colors[i], cex = 3)
  text(pt[1]+.02, pt[2]+.02, label = TeX(paste0('$t_', i, '$')), pos = 3, col = vec_colors[i], cex = 1.8)
  text(pt[1]+.02, pt[2] + .1, label = TeX('$\\rightarrow$'), pos = 3, col = vec_colors[i], cex = 1.8)
}

# === Phi arcs ===
# Define manual positions for phi labels
phi_label_positions <- list(
  c(1.05, 1.33),  # phi_1
  c(1.25, 1.2),  # phi_2
  c(1.2, 0.7),  # phi_3
  c(0.8, 0.66)   # phi_4
)

# Phi arcs
anchor <- c(0, 1)
for (i in 1:4) {
  vec <- clock_pts[[i]] - O
  draw.sector(
    start.degree = atan2(anchor[2], anchor[1]) * 180 / pi,
    end.degree = atan2(vec[2], vec[1]) * 180 / pi,
    clock.wise = TRUE,
    center = O,
    border = vec_colors[i],
    col = NA,
    lwd = 2,
    rou1 = 0.35 + i * 0.03
  )
  
  # Use manual positions
  label_pos <- phi_label_positions[[i]]
  text(label_pos[1], label_pos[2], 
       label = TeX(paste0('$\\varphi_', i, '$')), 
       col = vec_colors[i], 
       cex = 1.8)
}


# === Mini plot: distance vs time ===
par(fig = c(0.6, 1, 0.65, 1), new = TRUE)
distances <- sapply(clock_pts, function(pt) sqrt(sum((pt - O)^2)))
par(mgp = c(2, 0.5, 0))
par(mar = c(4, 4, 1, 1))

plot(
  1:4,
  distances,
  type = "l",
  xaxt = "n",
  ylab = "Shift",
  xlab = "Time",
  main = "",
  col = "black",
  cex.axis = 0.8,
  cex.lab = 1
)
axis(1, at = 1:4, labels = paste0("t", 1:4), cex.axis = 0.8)
points(1:4, distances, pch = 16, col = vec_colors, cex = 2)
box()

dev.off()
