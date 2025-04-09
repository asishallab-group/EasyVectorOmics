require(circlize)
require(RColorBrewer)
require(latex2exp)

addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

colors <- brewer.pal(6, "Dark4")
alpha.cols <- sapply(colors, addAlpha)

pdf('changeInTissuePreference_phi.pdf', width=4, height=4)


plot(
  c(1, .25, 1),
  c(1, 2, .25),
  col = c(colors[[4]], colors[[1]], colors[[2]]),
  bg = c(alpha.cols[[4]], alpha.cols[[1]], alpha.cols[[2]]),
  pch = c(20, 18, 18),
  lwd = 2.5,
  cex = 3,
  xlim = c(0, 2.5),
  ylim = c(0, 2.5),
  axes = FALSE,
  xlab = '',
  ylab = ''
)
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
arrows(1,
       1,
       1,
       .25,
       col = colors[[2]],
       lwd = 1.3,
       length = .15)
arrows(1,
       1,
       .25,
       2,
       col = colors[[1]],
       lwd = 1.3,
       length = .15)
lines(c(1, 1), c(1, (1 + sqrt(2))), col = "black", lwd = 1.5)
lines(c(1, 2), c(1, 0), col = "black", lwd = 1.5)
lines(c(1, 0), c(1, 0), col = "black", lwd = 1.5)
text(
  1,
  1,
  label = TeX('$\\bar{\\textit{D}}$'),
  pos = 4,
  offset = 1.25,
  col = colors[[4]],
  cex = .8
)
text(
  .25,
  2,
  label = TeX('$\\bar{\\textit{o}}$'),
  pos = 3,
  offset = 1.25,
  col = colors[[1]],
  cex = .8
)
text(
  1.2,
  .25,
  label = TeX('$\\bar{\\textit{p}}$'),
  pos = 1,
  offset = .755,
  col = colors[[2]],
  cex = .8
)
text(
  .6,
  .8,
  label = TeX('$\\varphi$'),
  pos = 2,
  offset = 1.25,
  col = colors[[3]],
  cex = 2
)
text(
  1,
  .55,
  label = TeX('$\\perp_{D}$'),
  pos = 4,
  offset = .5,
  col = colors[[2]],
  cex = .8
)
text(
  .5,
  1.45,
  label = TeX('$\\perp_{O}$'),
  pos = 4,
  offset = 1.0,
  col = colors[[1]],
  cex = .8
)
text(
  2,
  0,
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
  .35,
  -.1,
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
arrows(o[1] - .1,
       o[2],
       o[1] + unit_vec[1] - .1,
       o[2] + unit_vec[2],
       col = colors[[5]],
       lwd = 2,
       lty = 1,
       length = 0.15)

# Label for unit vector
text(o[1],
     o[2] - .4,
     label = TeX('$\\bar{u}_{po}$'),
     pos = 1,
     col = colors[[5]],
     cex = .8)

# difference vector
arrows(o[1],
       o[2],
       p[1],
       p[2],
       col = colors[[6]],
       lwd = 2,
       length = 0.15)

text(
  p[1] + ((o[1] - p[1]) / 2),
  p[2] + ((o[2] - p[2]) / 2) - 1,
  label = TeX('$\\bar{p} - \\bar{o}$'),
  pos = 3,
  col = colors[[6]],
  cex = .8
)

dev.off()