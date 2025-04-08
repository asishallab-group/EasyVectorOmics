require(circlize)
require(RColorBrewer)
require(latex2exp)

addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

colors <- brewer.pal(4, "Dark2")
alpha.cols <- sapply(colors, addAlpha)

pdf('changeInTissuePreference_phi.pdf', width=4, height=4)


plot(
  c(1, .5, 1),
  c(1, 1.5, .253),
  col = c(colors[[4]], colors[[1]], colors[[2]]),
  bg = c(alpha.cols[[4]], alpha.cols[[1]], alpha.cols[[2]]),
  pch = c(20, 18, 18),
  lwd = 4,
  cex = 4.25,
  xlim = c(0, 2.5),
  ylim = c(0, 2.5),
  axes = FALSE,
  xlab = '',
  ylab = ''
)
draw.sector(
  135,
  270,
  clock.wise = FALSE,
  center = c(1, 1),
  border = colors[[3]],
  col = alpha.cols[[3]],
  lwd = 3,
  rou1 = .55
)
arrows(1,
       .25,
       1,
       1,
       col = colors[[2]],
       lwd = 3,
       length = .25)
arrows(.5,
       1.5,
       1,
       1,
       col = colors[[1]],
       lwd = 3,
       length = .25)
lines(c(1, 1), c(1, (1 + sqrt(2))), col = "black", lwd = 1.5)
lines(c(1, 2), c(1, 0), col = "black", lwd = 1.5)
lines(c(1, 0), c(1, 0), col = "black", lwd = 1.5)
text(
  1,
  1,
  label = TeX('$\\bar{\\textit{\\nu}}$'),
  pos = 4,
  offset = 1.25,
  col = colors[[4]]
)
text(
  .6,
  1.7,
  label = TeX('$\\bar{\\textit{o}}$'),
  pos = 2,
  offset = 1.25,
  col = colors[[1]]
)
text(
  1.2,
  .25,
  label = TeX('$\\bar{\\textit{d}}$'),
  pos = 1,
  offset = .755,
  col = colors[[2]]
)
text(
  .6,
  .8,
  label = TeX('$\\varphi$'),
  pos = 2,
  offset = 1.25,
  col = colors[[3]],
  cex = 3
)
text(
  1,
  .55,
  label = TeX('$\\perp_{D}$'),
  pos = 4,
  offset = .5,
  col = colors[[2]]
)
text(
  .5,
  1.45,
  label = TeX('$\\perp_{O}$'),
  pos = 4,
  offset = 1.0,
  col = colors[[1]]
)
text(
  2,
  0,
  label = TeX('\\textit{tissue} Z'),
  pos = 3,
  offset = 1.25,
  cex = 1.0
)
text(
  1,
  2.2,
  label = TeX('\\textit{tissue} X'),
  pos = 4,
  offset = .2,
  cex = 1.0
)
text(
  .5,
  0,
  label = TeX('\\textit{tissue} Y'),
  pos = 3,
  offset = .1,
  cex = 1.0
)


dev.off()
