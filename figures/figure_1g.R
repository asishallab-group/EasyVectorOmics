require(RColorBrewer)
colors <- brewer.pal(12,"Paired")

pdf("figure_1g.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(1, 1, 1, 1))

# --- Panel Izquierdo: Campo vectorial con E1 y E2 ---
plot(NA, xlim = c(-2, 2), ylim = c(-2, 2), asp = 1, axes = FALSE, xlab = "", ylab = "")


# --- Parámetros configurables ---
a <- 1.8  # largo (eje E1)
b <- 1.5  # ancho (eje E2)
theta <-  -2*pi/3   # rotación en radianes (30°)

# --- Generar contorno tipo ameba ---
t <- seq(0, 2*pi, length.out = 200)
r <- 1 + 0.15 * sin(3 * t)  # forma ondulada
x <- a * r * cos(t)  # escalado horizontal (E1)
y <- b * r * sin(t)  # escalado vertical (E2)

# --- Rotar ---
rotation_matrix <- matrix(c(cos(theta), -sin(theta),
                            sin(theta),  cos(theta)), ncol = 2)
rotated <- rotation_matrix %*% rbind(x, y)

# --- Dibujar ---
lines(rotated[1, ] - 0.3, rotated[2, ] , col = "gray60", lwd = 2)

# Ejes con flechas en ambos extremos
# Ángulo de rotación en radianes
theta <- - pi / 6  # 30 grados

# Matriz de rotación
rotation_matrix <- matrix(c(cos(theta), -sin(theta),
                            sin(theta),  cos(theta)), ncol = 2)

# E1 original: (0, 2) → vertical
e1 <- rotation_matrix %*% c(0, 2) * 0.8
# E2 original: (2, 0) → horizontal
e2 <- rotation_matrix %*% c(2, 0) * 0.65

# E1 con flechas en ambos extremos
arrows(-e1[1]-0.6, -e1[2]-0.2, e1[1]-0.6, e1[2]-0.2, col = "black", lwd = 2, length = 0.2)
arrows(e1[1]-0.6, e1[2]-0.2,-e1[1]-0.6, -e1[2]-0.2, col = "black", lwd = 2, length = 0.2)

# E2 con flechas en ambos extremos
arrows(-e2[1], -e2[2], e2[1], e2[2], col = "black", lwd = 2, length = 0.2)
arrows(e2[1], e2[2], -e2[1], -e2[2], col = "black", lwd = 2, length = 0.2)

arrows(e1[1]*0.3-0.2, e1[2]*0.3-0.2,-e1[1]*0.3-0.2, -e1[2]*0.3-0.2, col = colors[[4]], lwd = 1.3, lty = 2,length = 0.2)


# Etiquetas
text(0.2, -1.2, labels = expression(E[1]), pos = 3, cex = 2, srt = 30)

# E2 label rotado (por ejemplo, -30 grados)
text(0.7, 0.1, labels = expression(E[2]), pos = 4, cex = 2, srt = 30)
# Curva gruesa que sigue E1
# 1. Datos de la curva (forma de distribución normal)
curve_y <- seq(-1, 1, length.out = 100)  # rango amplio sin escalar
curve_x <- dnorm(curve_y, mean = 0, sd = 0.4)  # forma de campana

# Escalamos para hacerlo más chico
curve_x <- curve_x * 0.4  # achica en eje X
curve_y <- curve_y * 0.5  # achica en eje Y si querés también

# Centrar la curva en el origen
curve_x <- curve_x - mean(curve_x)

# 2. Rotar la curva 30 grados
theta <- -pi / 6  # 30 grados en radianes
rotation_matrix <- matrix(c(cos(theta), -sin(theta),
                            sin(theta),  cos(theta)), ncol = 2)

rotated_coords <- rotation_matrix %*% rbind(curve_x, curve_y)

# 3. Dibujar la curva rotada
lines(rotated_coords[1, ]-0.7, rotated_coords[2, ]+0.6, lwd = 3, col = colors[[4]])

draw_vector <- function(x0, y0, length, angle_deg, col = "black", lwd = 3, arrow_size = 0.1) {
  # Convertir ángulo a radianes
  angle_rad <- angle_deg * pi / 180
  
  # Calcular punto final
  x1 <- x0 + length * cos(angle_rad)
  y1 <- y0 + length * sin(angle_rad)
  
  # Dibujar flecha
  arrows(x0, y0, x1, y1, length = arrow_size, col = colors[[6]], lwd = lwd)
}

draw_vector(-1.3, 0.7, 0.25, 270)    
draw_vector(-1.35, 0.1, 0.3, 250)    
draw_vector(-1.1, -.4, 0.35, 60)   
draw_vector(-.3, 1.1, 0.4, 30)   
draw_vector(0.5, 0.7, 0.3, 30)   
draw_vector(-.2, 0.4, 0.2, -30)   
draw_vector(0.3, 0.5, 0.25, -120)   
draw_vector(-.5, -0.8, 0.2, -60) 
draw_vector(-.5, -1.2, 0.3, -60) 
draw_vector(0.3, -1.3, 0.25, -60) 
draw_vector(0.5, -0.7, 0.25, -60) 
draw_vector(0.7, -0.2, 0.3, -60) 

# --- Panel Derecho: Proyección sobre E1 ---
plot(NA, xlim = c(0, 10), ylim = c(0, 1.2), type = "n",
     xlab = "", ylab = "", axes = FALSE, asp = 1)

# Ejes con flechas
arrows(0, -3, 8, -3, length = 0.1, lwd = 1.5)  # Eje E1
arrows(0, -3, 0, 4, length = 0.1, lwd = 1.5)   # Eje signal strength

# Etiquetas de los ejes
text(8, -3.5, expression(E[1]), cex = 1.8)
text(1.4, 4.5, "Signal strength", cex = 1.6)

# Secuencia de x
x <- seq(1, 8, length.out = 200)

# Curva 1: Distribución normal sesgada a la derecha
y1 <- dnorm(x, mean = 1.6, sd = 1.3)
y1 <- y1 / max(y1) * 3 # Escalar altura

# Curva 2: Forma de "S" acostada (tanh)
y2 <- 1 + 0.9 * tanh((x - 3))  # Ajustar ubicación y altura

# Dibujar ambas curvas
lines(x, y1, lwd = 4, col = colors[[10]])

# Valores de x
x <- seq(1, 8, length.out = 500)

# Curva compuesta: dos normales (montaña y valle)
y <- dnorm(x, mean = 1.8, sd = 1.2) - 
     0.5 * dnorm(x, mean = 4, sd = 1.5) + 
     0.4 * dnorm(x, mean = 6.5, sd = 0.7)

# Normalizar a altura máxima 1
y <- y / max(y) * 2.7

# Dibujar curva
lines(x, y, lwd = 4, col = colors[[12]])

dev.off()
