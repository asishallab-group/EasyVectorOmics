library(latex2exp)

pdf("figure_1f.pdf", width = 10, height = 4)
par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(0, 6), ylim = c(0, 5), axes = FALSE, xlab = "", ylab = "", asp = 1)

# Base positions
rows <- 3
cols <- 5
x0 <- -1
y0 <- 3
dx <- 0.9
dy <- 0.7

# Draw m_ij labels (no boxes)
for (i in 1:rows) {
  for (j in 1:cols) {
    text(x0 + (j - 0.5)*dx,
         y0 - (i - 0.5)*dy,
         labels = bquote(m[.(i)*.(j)]),
         cex = 1.5)
  }
}

# Draw vertical double bars (matrix notation)
bottom_y <- y0 - rows*dy + 0.2
top_y    <- y0 - 0.2

segments(x0 - 0.3, bottom_y, x0 - 0.3, top_y, lwd = 2)                 # Left bar
segments(x0 - 0.5, bottom_y, x0 - 0.5, top_y, lwd = 2)                 # Left bar
segments(x0 + cols*dx + 0.3, bottom_y, x0 + cols*dx + 0.3, top_y, lwd = 2)  # Right bar
segments(x0 + cols*dx + 0.5, bottom_y, x0 + cols*dx + 0.5, top_y, lwd = 2)  # Right bar

# Custom bracket drawing function
bracket <- function(x, width, y, height){
  data.frame(
    x = (c(0,1,4,5,6,9,10)/10 - 0.5) * width + x,
    y = c(0,1,1,2,1,1,0)/2 * height + y
  )
}

# Generate bracket below m31–m35
b <- bracket(x = 1.2, width = 4.1, y = 1, height = -0.4)
lines(b$x, b$y, lwd = 2)

# Label for rank spread
x_start <- x0 + 0.5 * dx
x_end <- x0 + 4.5 * dx
y_brace <- y0 - 3 * dy - 0.2
brace_height <- 0.8

text((x_start + x_end)/2, y_brace - brace_height * 0.5 -0.1,
     labels = "rank(M):\ndirectional spread", cex = 1.3)

# Top-right label
text(4, 4.5,
     labels = "Subfield", cex = 1.5)

# Frobenius norm label
text(6, 1,
     labels = TeX("Frobenius norm $|| M ||_F$:"), cex = 1.5)

text(6, 1 - 0.4,
     labels = "Signal strength",
     cex = 1.5)

# Subfield vector coordinates
subfield_x <- seq(x0 + 0.5 * dx, x0 + 4.5 * dx, length.out = 5)
subfield_y <- y0 + 0.8

# Draw subfield numbers
for (i in 1:5) {
  text(subfield_x[i]+0.3, subfield_y-0.1, labels = as.character(i), cex = 1)
}

# Vector directions (angled arrows in the cloud)
dx_list <- c(0.2, 0.3,  0.2,  0.27, 0.01)
dy_list <- c( 0.2, 0.1, -0.2, 0.1,  0.3)

for (i in 1:5) {
  arrows(
    x0 = subfield_x[i],
    y0 = subfield_y,
    x1 = subfield_x[i] + dx_list[i],
    y1 = subfield_y + dy_list[i],
    length = 0.06,
    lwd = 1.3
  )
}

# Dashed arrows pointing from cloud to matrix
for (i in 1:5) {
  arrows(
    x0 + (i - 0.5) * dx, subfield_y - 0.2,
    x0 + (i - 0.5) * dx, y0 + 0.1,
    length = 0.08, lwd = 1.2, lty = 3, col = "gray40"
  )
}

library(plotrix)

# Draw ellipse (the “cloud”)
x_centro <- mean(subfield_x)
y_centro <- subfield_y + 0.1

draw.ellipse(
  x = x_centro + 0.2,
  y = y_centro,
  a = 2.4,
  b = 0.6,
  border = "gray60",
  col = adjustcolor("gray90", alpha.f = 0.3),
  lwd = 1.5
)

# Principal flow vector
arrows(
  x0 = x_start + 5, y0 = y0,
  x1 = x_start + 7, y1 = y0 + 1.3,
  col = "blue", lwd = 7, length = 0.1
)
text(x_start + 7, y0 + 0.5 , labels = expression(f[m]), col = "blue", cex = 2)

# Rotation vector (horizontal with round head)
arrows(
  x0 = x_start + 5, y0 = y0 - 0.5,
  x1 = x_start + 7, y1 = y0 - 0.5,
  col = "red", lwd = 7, length = 0
)

points(
  x = x_start + 7,
  y = y0 - 0.5,
  pch = 16,
  col = "red",
  cex = 1.8
)

# Curved rotation arc
x0 <- x_start + 6
y0 <- y0 + 0.1
x1 <- x_start + 8
y1 <- y0 - 0.5

cx <- (x0 + x1) / 2
cy <- y0 - 0.6
r <- 0.6

# Angle for curved arc (250 degrees, curved toward the vector)
theta <- seq(-3*pi/4, -3*pi/4 + (250 * pi / 170), length.out = 200)

lines(
  x = cx + r * cos(theta),
  y = cy + r * sin(theta),
  lwd = 4,
  col = "red"
)

arrows(
  x0 = cx + r * cos(theta[2]),
  y0 = cy + r * sin(theta[2]),
  x1 = cx + r * cos(theta[1]),
  y1 = cy + r * sin(theta[1]),
  length = 0.07,
  col = "red",
  lwd = 4
)

# Label for rotation vector
text(x_start + 8, y0 - 0.5 , labels = expression(r[m]), col = "red", cex = 2)

dev.off()
