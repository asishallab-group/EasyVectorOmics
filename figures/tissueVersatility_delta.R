require(circlize)
require(RColorBrewer)
require(latex2exp)

addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

eucLen <- function(x) {
  sqrt(sum(sapply(x, function(i)
    i ^ 2)))
}

angleVecs <- function(x, y, conv.fac = (180 / pi)) {
  acos(sum(sapply(1:length(x), function(i)
    x[i] * y[i])) / (eucLen(x) * eucLen(y))) * conv.fac
}

colors <- brewer.pal(4, "Dark2")
alpha.cols <- sapply(colors, addAlpha)


dupl <- c(1.9, 1.2)
l.dupl <- eucLen(dupl)
orth <- c(.2, 1.75)
l.orth <- eucLen(orth)


pdf('tissueVersatility_delta.pdf',
    width = 4,
    height = 4)

old.par <- par(mgp = c(.4, 0.4, 0))

plot(
  rbind(orth, dupl),
  pch = 18,
  cex = 4,
  xlim = c(0, 2.2),
  ylim = c(0, 2.2),
  xlab = TeX('\\textit{tissue} X'),
  # '',
  ylab = TeX('\\textit{tissue} Y'),
  # '',
  cex.lab = 1.45,
  xaxt = 'n',
  yaxt = 'n',
  bty = 'l',
  # axes = FALSE,
  col = c(colors[[1]], colors[[2]])
)
# lines(c(0,0,0,2.15), c(0,2.15,0,0))
draw.sector(
  45,
  45 + angleVecs(c(1, 1), orth),
  clock.wise = FALSE,
  center = c(0, 0),
  border = colors[[1]],
  col = alpha.cols[[1]],
  lwd = 3,
  rou1 = l.orth
)
draw.sector(
  45,
  45 - angleVecs(c(1, 1), dupl),
  clock.wise = TRUE,
  center = c(0, 0),
  border = colors[[2]],
  col = alpha.cols[[2]],
  lwd = 3,
  rou1 = l.dupl
)
abline(0, 1, lwd = 4, col = colors[[4]])
# text(
#   .1,
#   2.15,
#   label = TeX('$\\textit{tissue} X$'),
#   pos = 3,
#   offset = .3,
#   cex = .75
# )
# text(
#   2.05,
#   0,
#   label = TeX('$\\textit{tissue} Y$'),
#   pos = 3,
#   offset = .5,
#   cex = .75
# )
text(
  orth[1],
  orth[2],
  label = TeX('$\\bar{\\textit{o}}$'),
  pos = 3,
  offset = 1,
  col = colors[[1]]
)
text(
  orth[1],
  orth[2],
  label = TeX('$\\delta(\\bar{\\textit{o}})$'),
  pos = 4,
  offset = 2.25,
  col = colors[[1]],
  cex = 1.5
)
text(
  dupl[1],
  dupl[2],
  label = TeX('$\\bar{\\textit{p}}$'),
  pos = 1,
  offset = 1.25,
  col = colors[[2]]
)
text(
  dupl[1] + .1,
  dupl[2],
  label = TeX('$\\delta(\\bar{\\textit{p}})$'),
  pos = 3,
  offset = 1.25,
  col = colors[[2]],
  cex = 1.5
)
text(
  1.75,
  1.75,
  label = TeX('$\\bar{\\textit{D}}$'),
  pos = 3,
  offset = 1.25,
  col = colors[[4]],
  cex = 1.5
)

# Calculate vectors
diff_vec <- dupl - orth
unit_vec <- diff_vec / eucLen(diff_vec)

# Difference vector
arrows(
  orth[1],
  orth[2],
  dupl[1],
  dupl[2],
  col = "darkblue",
  lwd = 2,
  length = 0.15
)

# Unit vector
arrows(
  orth[1],
  orth[2] - .1,
  orth[1] + unit_vec[1],
  orth[2] + unit_vec[2] - .1,
  col = "blue",
  lwd = 2,
  length = 0.1
)


# labels
text(
  (orth[1] + dupl[1]) / 2 + .2,
  (orth[2] + dupl[2]) / 2,
  label = TeX('$\\bar{p} - \\bar{o}$'),
  pos = 3,
  col = "darkblue",
  cex = 1.0
)

text(
  orth[1] + .15,
  orth[2] - .1,
  label = TeX('$\\bar{u}_{po}$'),
  pos = 1,
  col = "blue",
  cex = 1.0
)


par(old.par)

dev.off()