require(circlize)
require(RColorBrewer)
require(latex2exp)

addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

colors <- brewer.pal(7, "Dark2")
alpha.cols <- sapply(colors, addAlpha)

pdf('figure_1a_bis_medical.pdf', width = 8, height = 6)

plot(
  NA,
  pch = c(20, 18, 18),
  lwd = 2.5,
  cex = 3,
  xlim = c(0, 6),
  ylim = c(0, 6),  # antes era 4
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
text(0, 3.7, label = TeX('\\textit{intestine}'), pos = 4, offset = -.2, cex = 1.0)
text(5.7, .15, label = TeX('\\textit{liver}'), pos = 3, offset = .1, cex = 1.0)


centroid_1 <- c(1.5, 3.5)
centroid_2 <- c(4.5, 1.5)
# expression_vectors <- 

points(centroid_1[1], centroid_1[2], pch = 21, bg = alpha.cols[[1]], col = colors[[1]], cex = 3)
points(centroid_2[1], centroid_2[2], pch = 21, bg = alpha.cols[[1]], col = colors[[1]], cex = 4)


centroid_1_pts <- list(
  p1 = centroid_1 + c(0.0, 0.4),
  p2 = centroid_1 + c(0.4, 0.1),
  p3 = centroid_1 + c(0.3, -0.15),
  p4 = centroid_1 + c(-0.2, -0.2),
  p5 = centroid_1 + c(-0.25, 0.2)
)

set.seed(1)

for (tp in names(centroid_1_pts)) {
  pt <- centroid_1_pts[[tp]]
  points(pt[1], pt[2], pch = 21, col = colors[sample(1:2, 1)], cex = 1.5)
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
  points(pt[1], pt[2], pch = 21, col = colors[sample(1:2, 1)], cex = 2)
}

centroid_2_paras <- list(
    p1 = centroid_2 + c(-0.2, 1.2),
    p2 = centroid_2 + c(0.1, 1.8),
    p3 = centroid_2 + c(1.1, 0.3)
)

for (tp in names(centroid_2_paras)) {
  pt <- centroid_2_paras[[tp]]
  arrows(centroid_2[1], centroid_2[2], pt[1], pt[2],
         col = colors[2],
         lwd = 2.2,
         length = 0.1)
  points(pt[1], pt[2], pch = 21, col = colors[2], cex = 2)
}

text(
  centroid_2_paras[[1]][1] - 0.8,
  centroid_2_paras[[1]][2] - 0.1,
  label = TeX('$\\textit{susceptible\\ patients}$'),
  pos = 3,
  offset = 1,
  col = colors[[2]],
  cex = 1.1
)

text(
  centroid_2_paras[[1]][1] + .59,
  centroid_2_paras[[1]][2] - .5,
  labels = TeX(paste0("${s}_{max}$")),
  pos = 3,
  offset = 1,
  col = colors[[2]],
  cex = 1.4
)

text(
    centroid_2_paras[[1]][1] + .47,
    centroid_2_paras[[1]][2] + .2,
    labels = TeX("$\\rightarrow$"),
    cex = 1.4,
    col = colors[[2]]
  )

text(
  centroid_2_paras[[2]][1] + .2,
  centroid_2_paras[[2]][2] - .1,
  label = TeX('$\\textit{susceptible\\ patients}$'),
  pos = 3,
  offset = 1,
  col = colors[[2]],
  cex = 1.1
)

text(
  centroid_2_paras[[3]][1] - .15,
  centroid_2_paras[[3]][2] - .11,
  label = TeX('$\\textit{resistant\\ patients}$'),
  pos = 3,
  offset = 1,
  col = colors[[2]],
  cex = 1.1
)

# Ortholog trajectories (2 examples) and average
orth_1 <- c(1, 2, 1, 3, 2)
orth_2 <- c(3, 2, 2, 1, 1)
mean <- (orth_1 + orth_2) / 2
diff_vecs <- c(5, 4.5, 4, 4.3, 5)

# Map time to space
map_time_to_coord <- function(x, y_base, step=0.4) {
  coords <- lapply(1:length(x), function(i) c(0.6 + i*step, y_base + x[i] * 0.3))
  return(do.call(rbind, coords))
}

orth_1_pts <- map_time_to_coord(orth_1, 1.2)
orth_2_pts <- map_time_to_coord(orth_2, 1.2)
mean_pts  <- map_time_to_coord(mean, 1.2)
diff_pts <- matrix(c(
  0.8, 2.7,
  1.3, 2.5,
  1.8, 1.4,
  2.5, 1.1,
  3.1, 0.7
), ncol = 2, byrow = TRUE)

# Smooth lines
lines(spline(orth_1_pts), col = colors[[1]], lwd = 3)
lines(spline(orth_2_pts), col = colors[[1]], lwd = 3)
lines(spline(mean_pts), col = alpha.cols[[1]], lwd = 6)
lines(spline(diff_pts), col = colors[[2]], lwd = 3)

# Points
points(orth_1_pts, pch = 16, col = colors[[1]])
points(orth_2_pts, pch = 16, col = colors[[1]])
points(mean_pts, pch = 16, col = alpha.cols[[1]])
points(diff_pts, pch = 16, col = colors[[2]])

# Shift vectors from mean to diff points
for (i in 1:nrow(mean_pts)) {
  arrows(mean_pts[i, 1], mean_pts[i, 2],
         diff_pts[i, 1], diff_pts[i, 2],
         length = 0.1, col = colors[[3]], lwd = 2.2)
}

s_vec <- matrix(c(
  0, 0,
  0.05,- 0.22 ,
  0.05,- 0.05 ,
  0.05, -0.14,
  0.4, 0.8
), ncol = 2, byrow = TRUE)

for (i in 1:nrow(mean_pts)) {
  
  x_coord <- (mean_pts[i, 1] + diff_pts[i, 1]) / 2 - 0.1
  y_coord <- (mean_pts[i, 2] + diff_pts[i, 2]) / 2 - 0.5

  if (i == 5 ) {
    # print(s_vec[i, 1])
    subscript <- paste0("t", i, ",max") 
    # x_coord <- x_coord + s_vec[i, 1]
    # y_coord <- y_coord + s_vec[i, 2]
    arrow_coord <- x_coord + 0.2
   
  }
  else {
    subscript <- paste0("t", i)
    arrow_coord <- x_coord

  }

#   if(i == 2 || i == 3 ) {
#     y_coord <- y_coord - 0.3  
#     x_coord <- x_coord + 0.2
#     arrow_coord <- x_coord 

# }

  text(
    x = x_coord + s_vec[i, 1],
    y = y_coord + s_vec[i, 2],
    labels = TeX(paste0("${s}_{", subscript, "}$")),
    cex = 1.5,
    col = colors[[3]]
  )

  text(
    x = arrow_coord,
    y = y_coord+ s_vec[i, 2]+0.22,
    labels = TeX("$\\rightarrow$"),
    cex = 1.4,
    col = colors[[3]]
  )
}



dev.off()