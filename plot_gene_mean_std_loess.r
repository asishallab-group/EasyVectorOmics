library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(tidyr)

# ---------
# 1) Cargar datos
# ---------
df <- read.csv("gene_mean_std_loess.csv")

# columnas de Fortran
f_cols <- grep("^y_smooth_span_", names(df), value = TRUE)

# Extraer spans numéricos desde el nombre.
# Ajusta este parser a tu convención exacta.
# Ejemplos esperados: y_smooth_span_p200 -> 0.2
span_from_name <- function(nm) {
  # toma el número después de "_p"
  m <- str_match(nm, "_p([0-9]+)$")[,2]
  as.numeric(m) / 1000
}

spans <- span_from_name(f_cols)
names(spans) <- f_cols

# ---------
# 2) Calcular LOESS en R sobre los mismos x/y
# ---------
# Recomendación: ordenar por x para estabilidad visual
df2 <- df %>% arrange(x_mean)

# Si hay duplicados de x_mean (muy común), R loess puede volverse inestable.
# Una solución sencilla para comparar curvas: agrupar por x_mean y promediar y_std.
# (Opcional, pero recomendado)
agg <- df2 %>%
  group_by(x_mean) %>%
  summarize(y_std = mean(y_std, na.rm = TRUE), .groups = "drop") %>%
  arrange(x_mean)

x <- agg$x_mean
y <- agg$y_std

# Puedes cambiar degree a 1 o 2 para que coincida con tu Fortran
degree_used <- 2

r_pred <- lapply(spans, function(sp) {
  fit <- loess(y_std ~ x_mean,
               data = agg,
               span = sp,
               degree = degree_used,
               family = "symmetric",)
  predict(fit, newdata = data.frame(x_mean = x))
})

# Convertir a dataframe (mismas x que agg)
r_df <- data.frame(x_mean = x, y_std = y)
for (i in seq_along(f_cols)) {
  colname <- f_cols[i]
  r_df[[paste0(colname, "_R")]] <- as.numeric(r_pred[[i]])
}

# ---------
# 3) Pasar a formato largo: Fortran vs R
# ---------
# Fortran (en tu df original) lo ploteamos sobre agg$x_mean también.
# Para eso, hacemos una tabla "fortran" re-muestreada/interpolada si hace falta.
# Si tu archivo ya tiene filas por punto (gene) y no por x único,
# lo mejor es construir una curva por x_mean único también:

fortran_long <- df2 %>%
  select(x_mean, y_std, all_of(f_cols)) %>%
  group_by(x_mean) %>%   # si hay duplicados, promediamos la salida Fortran
  summarize(
    y_std = mean(y_std, na.rm = TRUE),
    across(all_of(f_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(x_mean) %>%
  pivot_longer(
    cols = all_of(f_cols),
    names_to = "span_col",
    values_to = "y_smooth"
  ) %>%
  mutate(
    source = "Fortran",
    span = span_from_name(span_col)
  )

r_long <- r_df %>%
  select(x_mean, y_std, ends_with("_R")) %>%
  pivot_longer(
    cols = ends_with("_R"),
    names_to = "span_col",
    values_to = "y_smooth"
  ) %>%
  mutate(
    source = "R",
    span = span_from_name(str_replace(span_col, "_R$", ""))
  )

plot_long <- bind_rows(fortran_long, r_long) %>%
  mutate(
    span = factor(span, levels = sort(unique(span)))
  )

# ---------
# 4) Plot: Fortran sólido, R dashed
# ---------
p <- ggplot() +
  geom_point(data = agg, aes(x = x_mean, y = y_std),
             color = "grey70", alpha = 0.35, size = 0.7) +
  geom_line(
    data = plot_long,
    aes(x = x_mean, y = y_smooth, color = span, linetype = source, group = interaction(span, source)),
    linewidth = 1
  ) +
  scale_linetype_manual(values = c(Fortran = "solid", R = "dashed")) +
  labs(
    title = paste0("Mean vs Std: LOESS comparison (degree=", degree_used, ")"),
    x = "x_mean",
    y = "y_std / smoothed",
    color = "Span",
    linetype = "Source"
  ) +
  theme_minimal()

ggsave("gene_mean_std_loess_fortran_vs_r.png", plot = p, width = 11, height = 7, dpi = 200)
print(p)
