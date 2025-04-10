require(circlize)
require(RColorBrewer)
require(latex2exp)

addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

colors <- brewer.pal(7, "Dark2")
alpha.cols <- sapply(colors, addAlpha)

pdf('expressionChange.pdf', width=5, height=4)

plot(
  NA,
  pch = c(20, 18, 18),
  lwd = 2.5,
  cex = 3,
  xlim = c(0, 6),
  ylim = c(0, 4),
  axes = FALSE,
  xlab = '',
  ylab = ''
)

for (i in 0:10){
    segments(x0 = i/2 + .6, y0 = .6, x1 = i/2 + .6, y1 = 4, col = "grey", lwd = .8)
}

for (i in 0:6){
    segments(x0 = .6, y0 = i/2 + .6, x1 = 6, y1 = i/2 + .6, col = "grey", lwd = .8)
}

lines(c(.6, .6), c(.6, 4), col = "black", lwd = 1.5)              
lines(c(.6, 6), c(.6, .6), col = "black", lwd = 1.5)              

# Tissue labels
text(0, 3.9, label = TeX('\\textit{leaf}'), pos = 4, offset = .2, cex = 1.0)
text(5.9, .3, label = TeX('\\textit{root}'), pos = 3, offset = .1, cex = 1.0)


centroid_1 <- c(1.5, 3.5)
centroid_2 <- c(4.5, 1.5)
# expression_vectors <- 

points(centroid_1[1], centroid_1[2], pch = 21, bg = alpha.cols[[1]], col = colors[[1]], cex = 1)
points(centroid_2[1], centroid_2[2], pch = 21, bg = alpha.cols[[1]], col = colors[[1]], cex = 2)


centroid_1_pts <- list(
  p1 = centroid_1 + c(0.0, 0.3),
  p2 = centroid_1 + c(0.2, 0.1),
  p3 = centroid_1 + c(0.1, -0.15),
  p4 = centroid_1 + c(-0.15, -0.2),
  p5 = centroid_1 + c(-0.1, 0.2)
)

for (tp in names(centroid_1_pts)) {
  pt <- centroid_1_pts[[tp]]
  points(pt[1], pt[2], pch = 21, col = colors[sample(1:2, 1)], cex = 0.5)
}

centroid_2_pts <- list(
  p1 = centroid_2 + c(-0.25, -0.25),
  p2 = centroid_2 + c(0.2, 0.4),
  p3 = centroid_2 + c(0.3, -0.2),
  p4 = centroid_2 + c(-0.2, -0.5),
  p5 = centroid_2 + c(-0.3, 0.3)
)

for (tp in names(centroid_2_pts)) {
  pt <- centroid_2_pts[[tp]]
  points(pt[1], pt[2], pch = 21, col = colors[sample(1:2, 1)], cex = 1)
}

centroid_2_paras <- list(
    p1 = centroid_2 + c(-0.2, 1.5),
    p2 = centroid_2 + c(0.1, 1.7),
    p3 = centroid_2 + c(1.4, 0.3)
)

for (tp in names(centroid_2_paras)) {
  pt <- centroid_2_paras[[tp]]
  arrows(centroid_2[1], centroid_2[2], pt[1], pt[2],
         col = colors[2],
         lwd = 1.3,
         length = 0.1)
  points(pt[1], pt[2], pch = 21, col = colors[2], cex = 1)
}

text(
  centroid_2_paras[[1]][1] - .3,
  centroid_2_paras[[1]][2] - .2,
  label = TeX('$\\textit{pea}$'),
  pos = 3,
  offset = 1,
  col = colors[[2]],
  cex = .8
)

text(
  centroid_2_paras[[2]][1] + .3,
  centroid_2_paras[[2]][2] - .2,
  label = TeX('$\\textit{pea}$'),
  pos = 3,
  offset = 1,
  col = colors[[2]],
  cex = .8
)

text(
  centroid_2_paras[[3]][1] - .25,
  centroid_2_paras[[3]][2] - .15,
  label = TeX('$\\textit{bean}$'),
  pos = 3,
  offset = 1,
  col = colors[[2]],
  cex = .8
)

rect(0.63, 0.63, 4.1, 3.1, col = "white", border = NA)


usr <- par("usr") 

# Deine gewünschten Plot-Koordinaten (in "user space")
x_lo <- 0.9
x_hi <- 4.4
y_lo <- 0.7
y_hi <- 3.4

# Umrechnung in relative fig-Koordinaten (zwischen 0 und 1)
fig_x_lo <- (x_lo - usr[1]) / (usr[2] - usr[1])
fig_x_hi <- (x_hi - usr[1]) / (usr[2] - usr[1])
fig_y_lo <- (y_lo - usr[3]) / (usr[4] - usr[3])
fig_y_hi <- (y_hi - usr[3]) / (usr[4] - usr[3])

par(fig = c(fig_x_lo, fig_x_hi, fig_y_lo, fig_y_hi), new = TRUE)
par(mgp = c(1, 0.5, 0))

orth_1 <- c(1, 2, 1, 3, 2)
diff_vecs <- c(3, 4, 3, 4, 3)

plot(
  1:5,
  diff_vecs,
  type = "l",
  pch = 16,
  cex = 0.5,
  xaxt = "n",
  yaxt = "n",
  ylab = "distance",
  xlab = "time",
  main = "",
  col = colors[[2]],
  cex.axis = 0.5,
  cex.lab = 0.5,
  ylim = c(1,5)
)
axis(1, at = 1:5, labels = paste0("t", 1:5), cex.axis = 0.5)
axis(2, at = 1:5, labels = 1:5, cex.axis = 0.5)

orth_2 <- c(3, 2, 2, 1, 1)
mean <- c(mean(c(orth_1[1], orth_2[1])), 
    mean(c(orth_1[2], orth_2[2])), 
    mean(c(orth_1[3], orth_2[3])), 
    mean(c(orth_1[4], orth_2[4])),
    mean(c(orth_1[5], orth_2[5])))

lines(1:5, orth_2, type = "l", col = colors[[1]], pch = 16)
lines(1:5, mean, type = "l", col = "red", pch = 16)
lines(1:5, orth_1, type = "l", col = colors[[1]], pch = 16)

points(1:5, orth_1, pch = 16, cex = 0.5, col = colors[1])
points(1:5, mean, pch = 16, cex = 0.5, col = "red")
points(1:5, orth_2, pch = 16, cex = 0.5, col = colors[1])
points(1:5, diff_vecs, pch = 16, cex = 0.5, col = colors[2])

box()

dev.off()