# Parcial-1--Ecolog-a-Acuatica
Evaluación de la Calidad de Agua y Macroinvertebrados Indicadores en un Humedal Impactado por Actividad Industrial
# Ecology Matrix Analysis

Este proyecto analiza una matriz ecológica almacenada en Excel utilizando Python y pandas.

# ANÁLISIS DE DATOS FISICOQUÍMICOS Y BIOLÓGICOS
# Autor: [Tu nombre]
# Fecha: [Fecha]

# Cargar librerías necesarias
library(tidyverse)
library(ggplot2)
library(corrplot)
library(patchwork)

# 1. Cargar los datos
fisicoquimicas <- read_csv("data/parametros_fisicoquimicos.csv")
biologicas <- read_csv("data/parametros_biologicos.csv")

# Unir ambos datasets por ID Muestra
datos_completos <- fisicoquimicas %>% 
  left_join(biologicas, by = c("ID Muestra", "Sitio", "Tipo"))

# 2. Análisis exploratorio
# Resumen estadístico por tipo de sitio
resumen_fisico <- datos_completos %>%
  group_by(Tipo) %>%
  summarise(
    across(c(Pb (µg/g), Cd (µg/g), Zn (µg/g), pH, Oxigeno disuelto (mg/L)),
           list(mean = mean, sd = sd, min = min, max = max), na.rm = TRUE)
  )

resumen_bio <- datos_completos %>%
  group_by(Tipo) %>%
  summarise(
    across(c(EPT (ind./m²), Oligochaeta (ind./m²), BMWP, Shannon (H')),
           list(mean = mean, sd = sd, min = min, max = max), na.rm = TRUE)
  )

# Guardar resúmenes
write_csv(resumen_fisico, "results/tablas/resumen_fisicoquimico.csv")
write_csv(resumen_bio, "results/tablas/resumen_biologico.csv")

# 3. Visualizaciones

# Boxplots para metales pesados
pb_plot <- ggplot(datos_completos, aes(x = Tipo, y = Pb (µg/g), fill = Tipo)) +
  geom_boxplot() +
  labs(title = "Concentración de Plomo por Tipo de Sitio",
       y = "Pb (µg/g)", x = "") +
  theme_minimal()

cd_plot <- ggplot(datos_completos, aes(x = Tipo, y = Cd (µg/g), fill = Tipo)) +
  geom_boxplot() +
  labs(title = "Concentración de Cadmio por Tipo de Sitio",
       y = "Cd (µg/g)", x = "") +
  theme_minimal()

zn_plot <- ggplot(datos_completos, aes(x = Tipo, y = Zn (µg/g), fill = Tipo)) +
  geom_boxplot() +
  labs(title = "Concentración de Zinc por Tipo de Sitio",
       y = "Zn (µg/g)", x = "") +
  theme_minimal()

# Combinar gráficos
metales_plot <- (pb_plot | cd_plot | zn_plot) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')

ggsave("results/graficos/boxplot_metales.png", metales_plot, width = 12, height = 6)

# Histogramas para parámetros fisicoquímicos
ph_hist <- ggplot(datos_completos, aes(x = pH, fill = Tipo)) +
  geom_histogram(binwidth = 0.1, alpha = 0.7, position = "identity") +
  labs(title = "Distribución de pH", x = "pH", y = "Frecuencia") +
  theme_minimal()

od_hist <- ggplot(datos_completos, aes(x = Oxigeno disuelto (mg/L), fill = Tipo)) +
  geom_histogram(binwidth = 0.2, alpha = 0.7, position = "identity") +
  labs(title = "Distribución de Oxígeno Disuelto", 
       x = "Oxígeno Disuelto (mg/L)", y = "Frecuencia") +
  theme_minimal()

ggsave("results/graficos/histogramas_fisicoquimicos.png", 
       (ph_hist | od_hist), width = 12, height = 6)

# Correlación entre variables
cor_vars <- datos_completos %>%
  select(Pb (µg/g), Cd (µg/g), Zn (µg/g), pH, Oxigeno disuelto (mg/L),
         EPT (ind./m²), BMWP, Shannon (H'), Dist. Industria (km))

cor_matrix <- cor(cor_vars, use = "complete.obs")

png("results/graficos/correlacion_variables.png", width = 800, height = 800)
corrplot(cor_matrix, method = "circle", type = "upper", 
         tl.col = "black", tl.srt = 45, diag = FALSE)
dev.off()

# 4. Análisis de índices biológicos vs contaminación
ggplot(datos_completos, aes(x = Pb (µg/g), y = EPT (ind./m²), color = Tipo)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Relación entre Plomo y Abundancia de EPT",
       x = "Pb (µg/g)", y = "EPT (ind./m²)") +
  theme_minimal()

ggsave("results/graficos/EPT_vs_Pb.png", width = 8, height = 6)

# 5. Análisis de distancia a industria
ggplot(datos_completos, aes(x = Dist. Industria (km), y = Pb (µg/g), color = Tipo)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Concentración de Plomo vs Distancia a Industria",
       x = "Distancia a industria (km)", y = "Pb (µg/g)") +
  theme_minimal()

ggsave("results/graficos/Pb_vs_distancia.png", width = 8, height = 6)

# 6. Análisis multivariado (PCA)
pca_data <- datos_completos %>%
  select(Pb (µg/g), Cd (µg/g), Zn (µg/g), pH, Oxigeno disuelto (mg/L),
         EPT (ind./m²), BMWP, Shannon (H')) %>%
  scale()

pca_result <- prcomp(pca_data, scale. = TRUE)

pca_plot <- as.data.frame(pca_result$x) %>%
  bind_cols(Tipo = datos_completos$Tipo) %>%
  ggplot(aes(x = PC1, y = PC2, color = Tipo)) +
  geom_point(size = 3) +
  labs(title = "Análisis de Componentes Principales",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  theme_minimal()

ggsave("results/graficos/PCA.png", pca_plot, width = 8, height = 6)
