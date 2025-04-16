require(circlize)
require(RColorBrewer)
require(latex2exp)

addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

colors <- brewer.pal(7, "Dark2")
alpha.cols <- sapply(colors, addAlpha)

pdf('figure_1c.pdf', width=4, height=4)


plot(
  # c(1, .25, 1),
  # c(1, 2, .25),
  # col = c(colors[[4]], colors[[1]], colors[[2]]),
  # bg = c(alpha.cols[[4]], alpha.cols[[1]], alpha.cols[[2]]),
  NA,
  pch = c(20, 18, 18),
  lwd = 2.5,
  cex = 3,
  xlim = c(0, 2.5),
  ylim = c(0, 2.5),
  axes = FALSE,
  xlab = '',
  ylab = ''
)

D <- c(1, 1)
p <- c(1, .25)
o <- c(.25, 2)
points(D[1], D[2], pch = 21, bg = alpha.cols[[4]], col = colors[[4]], cex = 2)
points(o[1], o[2], pch = 21, bg = alpha.cols[[1]], col = colors[[1]], cex = 2)
points(p[1], p[2], pch = 21, bg = alpha.cols[[2]], col = colors[[2]], cex = 2)

draw.sector(
  127,
  270,
  clock.wise = FALSE,
  center = c(1, 1),
  border = colors[[3]],
  col = alpha.cols[[3]],
  lwd = 2,
  rou1 = .55
)
# Base vectors coming from (0,0)
arrows(1,
       1,
       1,
       .25,
       col = colors[[2]],
       lwd = 2,
       length = .1)
arrows(1,
       1,
       .25,
       2,
       col = colors[[1]],
       lwd = 2,
       length = .1)

# Coordinate System
lines(c(1, 1), c(1, (1 + sqrt(2))), col = "black", lwd = 1.5)
lines(c(1, 2), c(1, 0), col = "black", lwd = 1.5)
lines(c(1, 0), c(1, 0), col = "black", lwd = 1.5)

# Comment: \\vec{} is not supported, therefore each arrow is created and placed manually
# Label for D vector
text(
  .9,
  1,
  label = TeX('$\\textit{D}$'),
  pos = 4,
  offset = 1.25,
  col = colors[[4]],
  cex = .8
)
# Arrow for D vector
text(
  .9 - .04,
  1 + .11,
  label = TeX('$\\rightarrow$'),
  pos = 4,
  offset = 1.25,
  col = colors[[4]],
  cex = .8
)

# # Label for o point
# text(
#   .25,
#   1.85,
#   label = TeX('$\\textit{o}$'),
#   pos = 3,
#   offset = 1.25,
#   col = colors[[1]],
#   cex = .8
# )

# # Arrow for o label
# text(
#   .25,
#   1.85 + 0.1,
#   label = TeX('$\\rightarrow$'),
#   pos = 3,
#   offset = 1.25,
#   col = colors[[1]],
#   cex = .8
# )

# # Label for p point
# text(
#   1.15,
#   .26,
#   label = TeX('$\\textit{p}$'),
#   pos = 1,
#   offset = .755,
#   col = colors[[2]],
#   cex = .8
# )

# # Arrow for p point
# text(
#   1.15,
#   .26 + 0.1,
#   label = TeX('$\\rightarrow$'),
#   pos = 1,
#   offset = .755,
#   col = colors[[2]],
#   cex = .8
# )

# phi label
text(
  .65,
  1,
  label = TeX('$\\varphi$'),
  pos = 2,
  offset = 1.25,
  col = colors[[3]],
  cex = .8
)

# Label for base vector p
text(
  .95 +0.3,
  .45,
  label = TeX('$\\textit{{p}^{clock}}$'),
  pos = 1,
  offset = 1.5,
  col = colors[[2]],
  cex = .8
)

# Arrow for vec p
text(
  .95 - 0.04,
  .45 - .38,
  label = TeX('$\\rightarrow$'),
  pos = 4,
  offset = .5,
  col = colors[[2]],
  cex = .8
)

# Label for base vector o
text(
  .4-0.2 ,
  1.85 + 0.3,
  label = TeX('$\\textit{{o}^{clock}}$'),
  pos = 4,
  offset = 1.0,
  col = colors[[1]],
  cex = .8
)

# Arrow for vec o
text(
  .4 -0.27,
  1.85 + .38,
  label = TeX('$\\rightarrow$'),
  pos = 4,
  offset = 1.0,
  col = colors[[1]],
  cex = .8
)

# Axes Labels
text(
  1.89,
  0.08,
  label = TeX('\\textit{tissue} Z'),
  pos = 3,
  offset = 1.25,
  cex = .8
)
text(
  1,
  2.2,
  label = TeX('\\textit{tissue} X'),
  pos = 4,
  offset = .2,
  cex = .8
)
text(
  .15,
  .35,
  label = TeX('\\textit{tissue} Y'),
  pos = 3,
  offset = .15,
  cex = .8
)

# Defining p and o to calc distance
p <- c(1, 0.25)
o <- c(0.25, 2)
diff_vec <- p - o
unit_vec <- diff_vec / sqrt(sum(diff_vec^2))

# Unit vector
arrows(o[1] - .05,
       o[2],
       o[1] + unit_vec[1] - .05,
       o[2] + unit_vec[2],
       col = colors[[6]],
       lwd = 2.5,
       lty = 1,
       length = 0.1)

# Label for unit vector
text(o[1] - .05,
     o[2] - .4,
     label = TeX('$\\textit{{s}_{po}^{norm}}$'),
     pos = 1,
     col = colors[[6]],
     cex = .8)

# Arrow for unit vector
text(o[1] - .18,
     o[2] - .32,
     label = TeX('$\\rightarrow$'),
     pos = 1,
     col = colors[[6]],
     cex = .8)

# difference vector
arrows(o[1],
       o[2],
       p[1],
       p[2],
       col = colors[[7]],
       lwd = 2.5,
       length = 0.1)

# Label for difference vector
text(
  p[1] + ((o[1] - p[1]) / 2) + .097,
  p[2] + ((o[2] - p[2]) / 2) - 1.1,
  label = TeX('$\\textit{{s}_{po}}$'),
  pos = 3,
  col = colors[[7]],
  cex = .8
)

# Arrow for difference vector
text(
  p[1] + ((o[1] - p[1]) / 2) + .03,
  p[2] + ((o[2] - p[2]) / 2) - .93,
  label = TeX('$\\rightarrow$'),
  pos = 3,
  col = colors[[7]],
  cex = .8
)

# -------- Inset: Viewing Geometry Schematic (top-right) --------
par(fig = c(0.5, 1, 0.65, 1), new = TRUE)
par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(-2, 6), ylim = c(-2, 6), asp = 1, axes = FALSE, xlab = "", ylab = "")

# Define points
origin <- c(1,1.5)
tissueX <- c(4.5, 1.5)
tissueY <- c(1, 5)
tissueZ <- c(3.5, 2.8)
diagonal <- c(4, 4.5)

# Draw RAP plane (parallelogram for 3D effect)
polygon(
  x = c(-2, 1, 4.3, 1.3),
  y = c(2.8, 4, 0.2, -1),
  border = "gray50",
  col = "gray90"
)

# Origin
# points(origin[1], origin[2], pch = 16)

# # Axes arrows
arrows(origin[1], origin[2], tissueX[1], tissueX[2], length = 0.08, lwd = 1.2)
arrows(origin[1], origin[2], tissueY[1], tissueY[2], length = 0.08, lwd = 1.2)
arrows(origin[1], origin[2], tissueZ[1], tissueZ[2], length = 0.08, lwd = 1.2)
arrows(origin[1], origin[2], diagonal[1], diagonal[2], length = 0, lwd = 1.2, lty = 2, col=colors[[4]])

# # Labels
text(tissueX[1] + 0.1, tissueX[2], TeX("\\textit{tissue X}"), cex = 0.6, pos = 4)
text(tissueY[1], tissueY[2] - 0.5, TeX("\\textit{tissue Y}"), cex = 0.6, pos = 3)
text(tissueZ[1] + 0.1, tissueZ[2], TeX("\\textit{tissue Z}"), cex = 0.6, pos = 4)
text(diagonal[1] - 0.7, diagonal[2] + 0.1, "Diagonal", cex = 0.6, pos = 4, col=colors[[4]])
text(3.8, -1, "RAP", font = 1, cex = 0.8)

# # Angle arc from X to diagonal

origin <- c(1, 1.5)

# angle 
theta <- seq(2.5, pi/4, length.out = 100)
# angle radius
r <- 1.2

# angle points
x_arc <- origin[1] + r * cos(theta)
y_arc <- origin[2] + r * sin(theta)

# Draw angle
lines(x_arc, y_arc, col = "deeppink", lwd = 2)

dev.off()