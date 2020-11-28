# Gestion del conocimiento - Proyecto Univalle: cargar y explorar datos
# H. Achicanoy & C. Saavedra
# 2017

# R options
options(warn = -1); options(scipen = 999)

# Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(FactoMineR))
suppressMessages(library(factoextra))
suppressMessages(library(foreign))
suppressMessages(library(corrplot))
suppressMessages(library(polycor))
suppressMessages(library(psych))
suppressMessages(library(gplots))
suppressMessages(library(viridis))
suppressMessages(library(lsr))
suppressMessages(library(DescTools))

# ------------------------------------------------------- #
# Loading data
# ------------------------------------------------------- #
km_data <- read.spss(file = "../_data/Base GConocimiento PymeS  Valle_2017.sav", to.data.frame = T, use.value.labels = T) # F
names(km_data)

pilares <- km_data %>% select(P3_1:P3_8, P5_1:P9_6, P11_1:P13_16, P22_1:P25_5)
for(j in 1:ncol(pilares)){
  pilares[,j] <- factor(pilares[,j], levels = c("Totalmente en desacuerdo",
                                                "En desacuerdo",
                                                "Ni de acuerdo ni en desacuerdo",
                                                "De acuerdo",
                                                "Totalmente de acuerdo"), ordered = T)
}; rm(j)
pilares %>% glimpse
pilares %>% str

# ------------------------------------------------------- #
# Pilar aprendizaje
# ------------------------------------------------------- #
aprendizaje <- pilares %>% dplyr::select(P5_4, P5_5, P6_1, P6_2, P6_5, P7_2, P8_1, P13_13:P13_16)
aprendizaje %>% glimpse

cat.metrics <- function(var, var.name){
  mode.var <- DescTools::Mode(x = var)[1]
  entropy.var <- DescTools::Entropy(x = table(var), base = exp(1))/DescTools::Entropy(x = rep(24.8, 5), base = exp(1))
  results <- data.frame(Variable = var.name, Mode = mode.var, Entropy = entropy.var)
  return(results)
}

varList <- pilares %>% names %>% unique
varMetrics <- lapply(1:length(varList), function(i){
  var.name <- varList[i]
  eval(parse(text = paste0('var <- pilares$', varList[i])))
  return(cat.metrics(var = var, var.name = var.name))
})
varMetrics <- do.call(rbind, varMetrics)
varMetrics %>% ggplot(aes(x = reorder(Variable, -Entropy), y = Entropy)) + geom_bar(stat = "identity") +
  coord_flip() + theme_bw()

hist(varMetrics$Entropy)


# Entropy
DescTools::Entropy(x = table(aprendizaje$P5_4), base = exp(1))/DescTools::Entropy(x = rep(24.8, 5), base = exp(1))
DescTools::Entropy(x = table(aprendizaje$P13_16), base = exp(1))/DescTools::Entropy(x = rep(24.8, 5), base = exp(1))
DescTools::Entropy(x = table(aprendizaje$P7_2), base = exp(1))/DescTools::Entropy(x = rep(24.8, 5), base = exp(1))
DescTools::Entropy(x = table(aprendizaje$P6_1), base = exp(1))/DescTools::Entropy(x = rep(24.8, 5), base = exp(1))

# Analisis descriptivo
fqTable <- aprendizaje %>%
  gather(measure, value) %>%
  count(measure, value)
names(fqTable) <- c("Variable", "Categoria", "Frecuencia")
fqTable <- fqTable %>% dplyr::mutate(Porcentaje = Frecuencia/nrow(aprendizaje))
fqTable$Categoria <- factor(fqTable$Categoria, levels = c("Totalmente en desacuerdo",
                               "En desacuerdo",
                               "Ni de acuerdo ni en desacuerdo",
                               "De acuerdo",
                               "Totalmente de acuerdo"), ordered = T)
fqTable$Variable <- factor(fqTable$Variable, levels = c("P5_4", "P5_5", "P6_1",
                                                        "P6_2", "P6_5", "P7_2",
                                                        "P8_1", paste0("P13_", 13:16)),
                           ordered = T)

# Plot it (SAVE IT!!!!)
# x11()
gg <- fqTable %>% ggplot(aes(x = Categoria, y = Porcentaje*100)) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("Porcentaje (%)") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100)) +
  facet_wrap(~ Variable) + # , scales = "free"
  theme_bw() +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_text(size = 13, face = 'bold'),
        axis.text = element_text(size = 12))
ggsave("../_results/_descriptive_analysis/frecuencias_aprendizaje.png", plot = gg, width = 22, height = 10, units = "in"); rm(fqTable, gg)

# Chi-square test
options(warn=-1)
p.chisq = matrix(0, nrow=ncol(aprendizaje), ncol=ncol(aprendizaje), byrow=T)
for(i in 1:ncol(aprendizaje)){
  for(j in 1:ncol(aprendizaje)){
    p.chisq[i,j] = round(chisq.test(aprendizaje[,i],aprendizaje[,j])$p.value,3)
  }
}; rm(i); rm(j)

diag(p.chisq) = NA
colnames(p.chisq) = colnames(aprendizaje)
rownames(p.chisq) = colnames(aprendizaje)

color_scale = colorRampPalette(c("tomato3","lightyellow","lightseagreen"), space="rgb")(50)
png('../_results/_descriptive_analysis/chi_square_test_aprendizaje.png', height = 7, width = 7, units = "in", res = 300)
heatmap.2(p.chisq,
          main="Aprendizaje",
          key.title="Chi-square test",
          key.xlab="p-value",
          Rowv=F,
          Colv=NULL,
          col=color_scale,
          linecol=NULL,
          tracecol=NULL,
          density.info="density",
          denscol="blue",
          margins=c(11,11))
dev.off()#; rm(catVar, p.chisq, color_scale)

# Cramer's V test
options(warn=-1)
p.cramer = matrix(0, nrow=ncol(aprendizaje), ncol=ncol(aprendizaje), byrow=T)
for(i in 1:ncol(aprendizaje)){
  for(j in 1:ncol(aprendizaje)){
    p.cramer[i,j] = round(lsr::cramersV(aprendizaje[,i],aprendizaje[,j]),3)
  }
}; rm(i); rm(j)

# diag(p.cramer) = NA
colnames(p.cramer) = colnames(aprendizaje)
rownames(p.cramer) = colnames(aprendizaje)

x11()
png(file = '../_results/_descriptive_analysis/cramer_v_test_aprendizaje.png', height = 7, width = 7, units = "in", res = 300)
corrplot::corrplot(corr = p.cramer, method = "square")
dev.off()

# ------------------------------------------------------- #
# Pilar Tecnologia 
# ------------------------------------------------------- #
tecnologia <- pilares %>% dplyr::select(P3_4, P5_6, P7_4, P7_5, P11_1:P11_5, P13_5, P13_6, P13_7, P13_8)
tecnologia %>% glimpse

# Analisis descriptivo
fqTable2 <- tecnologia %>%
  gather(measure, value) %>%
  count(measure, value)
names(fqTable2) <- c("Variable", "Categoria", "Frecuencia")
fqTable2 <- fqTable2 %>% dplyr::mutate(Porcentaje = Frecuencia/nrow(aprendizaje))
fqTable2$Categoria <- factor(fqTable2$Categoria, levels = c("Totalmente en desacuerdo",
                                                          "En desacuerdo",
                                                          "Ni de acuerdo ni en desacuerdo",
                                                          "De acuerdo",
                                                          "Totalmente de acuerdo"), ordered = T)
fqTable2$Variable <- factor(fqTable2$Variable, levels = c("P3_4", "P5_6", "P7_4", "P7_5",
                                                          paste0("P11_", 1:5),
                                                          "P13_5", "P13_6", "P13_7", "P13_8"),
                           ordered = T)

# Plot it (SAVE IT!!!!)
# x11()
gg <- fqTable2 %>% ggplot(aes(x = Categoria, y = Porcentaje*100)) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("Porcentaje (%)") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100)) +
  facet_wrap(~ Variable) +
  theme_bw() +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_text(size = 13, face = 'bold'),
        axis.text = element_text(size = 12))
ggsave("../_results/_descriptive_analysis/frecuencias_tecnologia.png", plot = gg, width = 22, height = 10, units = "in"); rm(fqTable2, gg)


options(warn=-1)
p.chisq2 = matrix(0, nrow=ncol(tecnologia), ncol=ncol(tecnologia), byrow=T)
for(i in 1:ncol(tecnologia)){
  for(j in 1:ncol(tecnologia)){
    p.chisq2[i,j] = round(chisq.test(tecnologia[,i],tecnologia[,j])$p.value,3)
  }
}; rm(i); rm(j)

diag(p.chisq2) = NA
colnames(p.chisq2) = colnames(tecnologia)
rownames(p.chisq2) = colnames(tecnologia)

color_scale = colorRampPalette(c("tomato3","lightyellow","lightseagreen"), space="rgb")(50)
png('../_results/_descriptive_analysis/chi_test_tecnologia.png', height = 7, width = 7, units = "in", res = 300)
heatmap.2(p.chisq2,
          main="Tecnologia",
          key.title="Chi-square test",
          key.xlab="p-value",
          Rowv=F,
          Colv=NULL,
          col=color_scale,
          linecol=NULL,
          tracecol=NULL,
          density.info="density",
          denscol="blue",
          margins=c(11,11))
dev.off()

# Cramer's V test
options(warn=-1)
p.cramer2 = matrix(0, nrow=ncol(tecnologia), ncol=ncol(tecnologia), byrow=T)
for(i in 1:ncol(tecnologia)){
  for(j in 1:ncol(tecnologia)){
    p.cramer2[i,j] = round(lsr::cramersV(tecnologia[,i],tecnologia[,j]),3)
  }
}; rm(i); rm(j)

# diag(p.cramer) = NA
colnames(p.cramer2) = colnames(tecnologia)
rownames(p.cramer2) = colnames(tecnologia)

x11()
png(file = '../_results/_descriptive_analysis/cramer_v_test_tecnologia.png', height = 7, width = 7, units = "in", res = 300)
corrplot::corrplot(corr = p.cramer2, method = "square")
dev.off()

# ------------------------------------------------------- #
# Pilar Liderazgo 
# ------------------------------------------------------- #
liderazgo <- pilares %>% dplyr::select(P3_1:P3_3, P3_5:P3_8, P8_3, P8_4, P8_5, P13_1, P13_2, P13_3, P13_4)
liderazgo  %>% glimpse

# Analisis descriptivo
fqTable3 <- liderazgo  %>%
  gather(measure, value) %>%
  count(measure, value)
names(fqTable3) <- c("Variable", "Categoria", "Frecuencia")
fqTable3 <- fqTable3 %>% dplyr::mutate(Porcentaje = Frecuencia/nrow(liderazgo))
fqTable3$Categoria <- factor(fqTable3$Categoria, levels = c("Totalmente en desacuerdo",
                                                            "En desacuerdo",
                                                            "Ni de acuerdo ni en desacuerdo",
                                                            "De acuerdo",
                                                            "Totalmente de acuerdo"), ordered = T)
fqTable3$Variable <- factor(fqTable3$Variable, levels = c(paste0("P3_", 1:3), paste0("P3_", 5:8),
                                                          "P8_3", "P8_4", "P8_5", "P13_1", "P13_2",
                                                          "P13_3", "P13_4"),
                            ordered = T)

# Plot it (SAVE IT!!!!)
# x11()
gg <- fqTable3 %>% ggplot(aes(x = Categoria, y = Porcentaje*100)) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("Porcentaje (%)") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100)) +
  facet_wrap(~ Variable) +
  theme_bw() +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_text(size = 13, face = 'bold'),
        axis.text = element_text(size = 12))
ggsave("../_results/_descriptive_analysis/frecuencias_liderazgo.png", plot = gg, width = 22, height = 10, units = "in"); rm(fqTable3, gg)


options(warn=-1)
p.chisq3 = matrix(0, nrow=ncol(liderazgo), ncol=ncol(liderazgo), byrow=T)
for(i in 1:ncol(liderazgo)){
  for(j in 1:ncol(liderazgo)){
    p.chisq3[i,j] = round(chisq.test(liderazgo[,i],liderazgo[,j])$p.value,3)
  }
}; rm(i); rm(j)

diag(p.chisq3) = NA
colnames(p.chisq3) = colnames(liderazgo)
rownames(p.chisq3) = colnames(liderazgo)

color_scale = colorRampPalette(c("tomato3","lightyellow","lightseagreen"), space="rgb")(50)
png('../_results/_descriptive_analysis/chi_test_liderazgo.png', height = 7, width = 7, units = "in", res = 300)
heatmap.2(p.chisq3,
          main="Liderazgo",
          key.title="Chi-square test",
          key.xlab="p-value",
          Rowv=F,
          Colv=NULL,
          col=color_scale,
          linecol=NULL,
          tracecol=NULL,
          density.info="density",
          denscol="blue",
          margins=c(11,11))
dev.off()

# Cramer's V test
options(warn=-1)
p.cramer3 = matrix(0, nrow=ncol(liderazgo), ncol=ncol(liderazgo), byrow=T)
for(i in 1:ncol(liderazgo)){
  for(j in 1:ncol(liderazgo)){
    p.cramer3[i,j] = round(lsr::cramersV(liderazgo[,i],liderazgo[,j]),3)
  }
}; rm(i); rm(j)

# diag(p.cramer) = NA
colnames(p.cramer3) = colnames(liderazgo)
rownames(p.cramer3) = colnames(liderazgo)

x11()
png(file = '../_results/_descriptive_analysis/cramer_v_test_liderazgo.png', height = 7, width = 7, units = "in", res = 300)
corrplot::corrplot(corr = p.cramer3, method = "square")
dev.off()

# ------------------------------------------------------- #
# Pilar Organizacion 
# ------------------------------------------------------- #
organizacion <- pilares %>% dplyr::select(P5_1:P5_3, P6_3, P6_4, P7_1, P7_3, P8_2, P13_9, P13_11)
organizacion  %>% glimpse

# Analisis descriptivo
fqTable4 <- organizacion  %>%
  gather(measure, value) %>%
  count(measure, value)
names(fqTable4) <- c("Variable", "Categoria", "Frecuencia")
fqTable4 <- fqTable4 %>% dplyr::mutate(Porcentaje = Frecuencia/nrow(organizacion))
fqTable4$Categoria <- factor(fqTable4$Categoria, levels = c("Totalmente en desacuerdo",
                                                            "En desacuerdo",
                                                            "Ni de acuerdo ni en desacuerdo",
                                                            "De acuerdo",
                                                            "Totalmente de acuerdo"), ordered = T)
fqTable4$Variable <- factor(fqTable4$Variable, levels = c(paste0("P5_", 1:3), "P6_3", "P6_4",
                                                          "P7_1", "P7_3", "P8_2", "P13_9", "P13_11"),
                            ordered = T)

# Plot it (SAVE IT!!!!)
# x11()
gg <- fqTable4 %>% ggplot(aes(x = Categoria, y = Porcentaje*100)) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("Porcentaje (%)") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100)) +
  facet_wrap(~ Variable) +
  theme_bw() +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_text(size = 13, face = 'bold'),
        axis.text = element_text(size = 12))
ggsave("../_results/_descriptive_analysis/frecuencias_organizacion.png", plot = gg, width = 22, height = 10, units = "in"); rm(fqTable4, gg)


options(warn=-1)
p.chisq4 = matrix(0, nrow=ncol(organizacion), ncol=ncol(organizacion), byrow=T)
for(i in 1:ncol(organizacion)){
  for(j in 1:ncol(organizacion)){
    p.chisq4[i,j] = round(chisq.test(organizacion[,i],organizacion[,j])$p.value,3)
  }
}; rm(i); rm(j)

diag(p.chisq4) = NA
colnames(p.chisq4) = colnames(organizacion)
rownames(p.chisq4) = colnames(organizacion)

color_scale = colorRampPalette(c("tomato3","lightyellow","lightseagreen"), space="rgb")(50)
png('../_results/_descriptive_analysis/chi_test_organizacion.png', height = 7, width = 7, units = "in", res = 300)
heatmap.2(p.chisq4,
          main="Organizacion",
          key.title="Chi-square test",
          key.xlab="p-value",
          Rowv=F,
          Colv=NULL,
          col=color_scale,
          linecol=NULL,
          tracecol=NULL,
          density.info="density",
          denscol="blue",
          margins=c(11,11))
dev.off()

# Cramer's V test
options(warn=-1)
p.cramer4 = matrix(0, nrow=ncol(organizacion), ncol=ncol(organizacion), byrow=T)
for(i in 1:ncol(organizacion)){
  for(j in 1:ncol(organizacion)){
    p.cramer4[i,j] = round(lsr::cramersV(organizacion[,i],organizacion[,j]),3)
  }
}; rm(i); rm(j)

# diag(p.cramer) = NA
colnames(p.cramer4) = colnames(organizacion)
rownames(p.cramer4) = colnames(organizacion)

png(file = '../_results/_descriptive_analysis/cramer_v_test_organizacion.png', height = 7, width = 7, units = "in", res = 300)
corrplot::corrplot(corr = p.cramer4, method = "square")
dev.off()

# ------------------------------------------------------- #
# Pilar Gestion Del Conocimiento
# ------------------------------------------------------- #
gconocimiento <- pilares %>% dplyr::select(P9_1:P9_6)
gconocimiento  %>% glimpse

# Analisis descriptivo
fqTable5 <- gconocimiento  %>%
  gather(measure, value) %>%
  count(measure, value)
names(fqTable5) <- c("Variable", "Categoria", "Frecuencia")
fqTable5 <- fqTable5 %>% dplyr::mutate(Porcentaje = Frecuencia/nrow(gconocimiento))
fqTable5$Categoria <- factor(fqTable5$Categoria, levels = c("Totalmente en desacuerdo",
                                                            "En desacuerdo",
                                                            "Ni de acuerdo ni en desacuerdo",
                                                            "De acuerdo",
                                                            "Totalmente de acuerdo"), ordered = T)
fqTable5$Variable <- factor(fqTable5$Variable, levels = c(paste0("P9_", 1:6)),
                            ordered = T)

# Plot it (SAVE IT!!!!)
x11()
gg <- fqTable5 %>% ggplot(aes(x = Categoria, y = Porcentaje*100)) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("Porcentaje (%)") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100)) +
  facet_wrap(~ Variable) +
  theme_bw() +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_text(size = 13, face = 'bold'),
        axis.text = element_text(size = 12))
ggsave("../_results/_descriptive_analysis/frecuencias_gestion_conocimiento.png", plot = gg, width = 22, height = 10, units = "in"); rm(fqTable5, gg)


options(warn=-1)
p.chisq5 = matrix(0, nrow=ncol(gconocimiento), ncol=ncol(gconocimiento), byrow=T)
for(i in 1:ncol(gconocimiento)){
  for(j in 1:ncol(gconocimiento)){
    p.chisq5[i,j] = round(chisq.test(gconocimiento[,i],gconocimiento[,j])$p.value,3)
  }
}; rm(i); rm(j)

diag(p.chisq5) = NA
colnames(p.chisq5) = colnames(gconocimiento)
rownames(p.chisq5) = colnames(gconocimiento)

color_scale = colorRampPalette(c("tomato3","lightyellow","lightseagreen"), space="rgb")(50)
png('../_results/_descriptive_analysis/chi_test_gestion_conocimiento.png', height = 7, width = 7, units = "in", res = 300)
heatmap.2(p.chisq5,
          main="Gestion Conocimiento",
          key.title="Chi-square test",
          key.xlab="p-value",
          Rowv=F,
          Colv=NULL,
          col=color_scale,
          linecol=NULL,
          tracecol=NULL,
          density.info="density",
          denscol="blue",
          margins=c(11,11))
dev.off()

# Cramer's V test
options(warn=-1)
p.cramer5 = matrix(0, nrow=ncol(gconocimiento), ncol=ncol(gconocimiento), byrow=T)
for(i in 1:ncol(gconocimiento)){
  for(j in 1:ncol(gconocimiento)){
    p.cramer5[i,j] = round(lsr::cramersV(gconocimiento[,i],gconocimiento[,j]),3)
  }
}; rm(i); rm(j)

# diag(p.cramer) = NA
colnames(p.cramer5) = colnames(gconocimiento)
rownames(p.cramer5) = colnames(gconocimiento)

png(file = '../_results/_descriptive_analysis/cramer_v_test_gestion_conocimiento.png', height = 7, width = 7, units = "in", res = 300)
corrplot::corrplot(corr = p.cramer5, method = "square")
dev.off()

# ------------------------------------------------------- #
# Pilar Tecnologia (aparte)
# ------------------------------------------------------- #
tegnologia2 <- pilares %>% dplyr::select(P12_1:P12_10)
tegnologia2  %>% glimpse

# Analisis descriptivo
fqTable6 <- tegnologia2  %>%
  gather(measure, value) %>%
  count(measure, value)
names(fqTable6) <- c("Variable", "Categoria", "Frecuencia")
fqTable6 <- fqTable6 %>% dplyr::mutate(Porcentaje = Frecuencia/nrow(tegnologia2))
fqTable6$Categoria <- factor(fqTable6$Categoria, levels = c("Totalmente en desacuerdo",
                                                            "En desacuerdo",
                                                            "Ni de acuerdo ni en desacuerdo",
                                                            "De acuerdo",
                                                            "Totalmente de acuerdo"), ordered = T)
fqTable6$Variable <- factor(fqTable6$Variable, levels = c(paste0("P12_", 1:10)),
                            ordered = T)

# Plot it (SAVE IT!!!!)
x11()
gg <- fqTable6 %>% ggplot(aes(x = Categoria, y = Porcentaje*100)) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("Porcentaje (%)") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100)) +
  facet_wrap(~ Variable) +
  theme_bw() +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title.x = element_text(size = 13, face = 'bold'),
        axis.title.y = element_text(size = 13, face = 'bold'),
        axis.text = element_text(size = 12))
ggsave("../_results/_descriptive_analysis/frecuencias_tecnologia_aparte.png", plot = gg, width = 22, height = 10, units = "in"); rm(fqTable6, gg)


options(warn=-1)
p.chisq6 = matrix(0, nrow=ncol(tegnologia2), ncol=ncol(tegnologia2), byrow=T)
for(i in 1:ncol(tegnologia2)){
  for(j in 1:ncol(tegnologia2)){
    p.chisq6[i,j] = round(chisq.test(tegnologia2[,i],tegnologia2[,j])$p.value,3)
  }
}; rm(i); rm(j)

diag(p.chisq6) = NA
colnames(p.chisq6) = colnames(tegnologia2)
rownames(p.chisq6) = colnames(tegnologia2)

color_scale = colorRampPalette(c("tomato3","lightyellow","lightseagreen"), space="rgb")(50)
png('../_results/_descriptive_analysis/chi_test_tecnologia_aparte.png', height = 7, width = 7, units = "in", res = 300)
heatmap.2(p.chisq6,
          main="Gestion Conocimiento",
          key.title="Chi-square test",
          key.xlab="p-value",
          Rowv=F,
          Colv=NULL,
          col=color_scale,
          linecol=NULL,
          tracecol=NULL,
          density.info="density",
          denscol="blue",
          margins=c(11,11))
dev.off()

# Cramer's V test
options(warn=-1)
p.cramer6 = matrix(0, nrow=ncol(tegnologia2), ncol=ncol(tegnologia2), byrow=T)
for(i in 1:ncol(tegnologia2)){
  for(j in 1:ncol(tegnologia2)){
    p.cramer6[i,j] = round(lsr::cramersV(tegnologia2[,i],tegnologia2[,j]),3)
  }
}; rm(i); rm(j)

# diag(p.cramer) = NA
colnames(p.cramer6) = colnames(tegnologia2)
rownames(p.cramer6) = colnames(tegnologia2)

png(file = '../_results/_descriptive_analysis/cramer_v_test_tegnologia_aparte.png', height = 7, width = 7, units = "in", res = 300)
corrplot::corrplot(corr = p.cramer6, method = "square")
dev.off()

# ------------------------------------------------------- #
# Aprendizaje, Tecnologia & Gestion del conocimiento
# ------------------------------------------------------- #

AprTecGC <- pilares %>% dplyr::select(P5_4, P5_5, P6_1, P6_2, P6_5, P7_2, P8_1, P13_13:P13_16, # Aprendizaje
                                      P3_4, P5_6, P7_4, P7_5, P11_1:P11_5, P13_5, P13_6, P13_7, P13_8, # Tecnologia
                                      P9_1:P9_6) # Gestion del conocimiento

options(warn = -1)
p.chisq7 = matrix(0, nrow=ncol(AprTecGC), ncol=ncol(AprTecGC), byrow=T)
for(i in 1:ncol(AprTecGC)){
  for(j in 1:ncol(AprTecGC)){
    p.chisq7[i,j] = round(chisq.test(AprTecGC[,i],AprTecGC[,j])$p.value,3)
  }
}; rm(i); rm(j)

diag(p.chisq7) = NA
colnames(p.chisq7) = colnames(AprTecGC)
rownames(p.chisq7) = colnames(AprTecGC)

write.csv(x = p.chisq7, file = "../_results/_descriptive_analysis/p_chisq_Apr_Tec_GC.csv", row.names = T)

color_scale = colorRampPalette(c("tomato3","lightyellow","lightseagreen"), space="rgb")(50)
png('../_results/_descriptive_analysis/chi_test_apr_tec_GC.png', height = 7, width = 7, units = "in", res = 300)
heatmap.2(p.chisq7,
          main = "Aprendizaje, Tecnologia,\nGestion Conocimiento",
          key.title = "Chi-square test",
          key.xlab = "p-value",
          Rowv = F,
          Colv = NULL,
          col = color_scale,
          linecol = NULL,
          tracecol = NULL,
          density.info = "density",
          denscol = "blue",
          margins = c(11,11))
dev.off()

# Cramer's V test
options(warn = -1)
p.cramer7 = matrix(0, nrow = ncol(AprTecGC), ncol = ncol(AprTecGC), byrow = T)
for(i in 1:ncol(AprTecGC)){
  for(j in 1:ncol(AprTecGC)){
    p.cramer7[i,j] = round(lsr::cramersV(AprTecGC[,i],AprTecGC[,j]),3)
  }
}; rm(i); rm(j)

# diag(p.cramer) = NA
colnames(p.cramer7) = colnames(AprTecGC)
rownames(p.cramer7) = colnames(AprTecGC)

write.csv(x = p.cramer7, file = "../_results/_descriptive_analysis/p_cramer_Apr_Tec_GC.csv", row.names = T)

png(file = '../_results/_descriptive_analysis/cramer_v_test_apr_tec_GC.png', height = 7, width = 7, units = "in", res = 300)
corrplot::corrplot(corr = p.cramer7, method = "square")
dev.off()

# ------------------------------------------------------- #
# Organizacion, Liderazgo & Gestion del conocimiento
# ------------------------------------------------------- #

OrgLidGC <- pilares %>% dplyr::select(P5_1:P5_3, P6_3, P6_4, P7_1, P7_3, P8_2, P13_9, P13_11, # Organizacion
                                      P3_1:P3_3, P3_5:P3_8, P8_3, P8_4, P8_5, P13_1, P13_2, P13_3, P13_4, # Liderazgo
                                      P9_1:P9_6) # Gestion del conocimiento

options(warn = -1)
p.chisq8 = matrix(0, nrow=ncol(OrgLidGC), ncol=ncol(OrgLidGC), byrow=T)
for(i in 1:ncol(OrgLidGC)){
  for(j in 1:ncol(OrgLidGC)){
    p.chisq8[i,j] = round(chisq.test(OrgLidGC[,i],OrgLidGC[,j])$p.value,3)
  }
}; rm(i); rm(j)

diag(p.chisq8) = NA
colnames(p.chisq8) = colnames(OrgLidGC)
rownames(p.chisq8) = colnames(OrgLidGC)

write.csv(x = p.chisq8, file = "../_results/_descriptive_analysis/p_chisq_Org_Lid_GC.csv", row.names = T)

color_scale = colorRampPalette(c("tomato3","lightyellow","lightseagreen"), space="rgb")(50)
png('../_results/_descriptive_analysis/chi_test_org_lid_GC.png', height = 7, width = 7, units = "in", res = 300)
heatmap.2(p.chisq8,
          main = "Organizacion, Liderazgo,\nGestion Conocimiento",
          key.title = "Chi-square test",
          key.xlab = "p-value",
          Rowv = F,
          Colv = NULL,
          col = color_scale,
          linecol = NULL,
          tracecol = NULL,
          density.info = "density",
          denscol = "blue",
          margins = c(11,11))
dev.off()

# Cramer's V test
options(warn = -1)
p.cramer8 = matrix(0, nrow = ncol(OrgLidGC), ncol = ncol(OrgLidGC), byrow = T)
for(i in 1:ncol(OrgLidGC)){
  for(j in 1:ncol(OrgLidGC)){
    p.cramer8[i,j] = round(lsr::cramersV(OrgLidGC[,i],OrgLidGC[,j]),3)
  }
}; rm(i); rm(j)

# diag(p.cramer) = NA
colnames(p.cramer8) = colnames(OrgLidGC)
rownames(p.cramer8) = colnames(OrgLidGC)

write.csv(x = p.cramer8, file = "../_results/_descriptive_analysis/p_cramer_Org_Lid_GC.csv", row.names = T)

png(file = '../_results/_descriptive_analysis/cramer_v_test_org_lid_GC.png', height = 7, width = 7, units = "in", res = 300)
corrplot::corrplot(corr = p.cramer8, method = "square")
dev.off()

#----------------------------------------------------------#
# Analisis Canonico Tecnologia / Aprendizaje 
#----------------------------------------------------------#
suppressMessages(library(GGally))
suppressMessages(library(CCA))

km_data <- read.spss(file = "../_data/Base GConocimiento PymeS  Valle_2017.sav", to.data.frame = T, use.value.labels = F) # F
aprendizaje <- km_data %>% dplyr::select(P5_4, P5_5, P6_1, P6_2, P6_5, P7_2, P8_1, P13_13:P13_16)
tecnologia <- km_data %>% dplyr::select(P3_4, P5_6, P7_4, P7_5, P11_1:P11_5, P13_5, P13_6, P13_7, P13_8)

# x11()
corrplot(cor(cbind(aprendizaje, tecnologia)), method = "square")

# X: tecnologia
# Y: aprendizaje
cc1 <- cc(tecnologia, aprendizaje)
cc1$cor[1]
cc1$cor[1:3]

cc1[3:4]

cc2 <- comput(tecnologia, aprendizaje, cc1)
#cc2 <- comput(liderazgo, organizacion, cc1)
cc2[3:6]
cc2[3:6]$corr.X.xscores %>% View
cc2[3:6]$corr.Y.yscores %>% View

x11()
ggplot(data = data.frame(P7_2 = aprendizaje$P7_2, P7_5 = tecnologia$P7_5), aes(x = P7_2, y = P7_5, size = 28, alpha = .6)) + geom_point()

# Tests of canonical dimensions
ev <- (1 - cc1$cor^2)

n <- dim(aprendizaje)[1]
p <- length(aprendizaje)
q <- length(tecnologia)
k <- min(p, q)
m <- n - 3/2 - (p + q)/2

w <- rev(cumprod(rev(ev)))

# initialize
d1 <- d2 <- f <- vector("numeric", k)

for (i in 1:k) {
  s <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
  si <- 1/s
  d1[i] <- p * q
  d2[i] <- m * s - p * q/2 + 1
  r <- (1 - w[i]^si)/w[i]^si
  f[i] <- r * d2[i]/d1[i]
  p <- p - 1
  q <- q - 1
}

pv <- pf(f, d1, d2, lower.tail = FALSE)
(dmat <- cbind(WilksL = w, F = f, df1 = d1, df2 = d2, p = pv))

# standardized psych canonical coefficients diagonal matrix of apre sd's
s1 <- diag(sqrt(diag(cov(aprendizaje))))
s1 %*% cc1$xcoef
# standardized psych canonical coefficients diagonal matrix of tec sd's
s2 <- diag(sqrt(diag(cov(tecnologia))))
s2 %*% cc1$ycoef

#----------------------------------------------------------#
# Analisis Canonico Liderazgo / Organizacion 
#----------------------------------------------------------#

liderazgo <- km_data %>% dplyr::select(P3_1:P3_3, P3_5:P3_8, P8_3, P8_4, P8_5, P13_1, P13_2, P13_3, P13_4)
organizacion <- km_data %>% dplyr::select(P5_1:P5_3, P6_3, P6_4, P7_1, P7_3, P8_2, P13_9, P13_11)

# x11()
corrplot(cor(cbind(liderazgo, organizacion)), method = "square")

# X: liderazgo
# Y: organizacion
cc3 <- cc(liderazgo, organizacion)
cc3$cor[1]
cc3$cor[1:3]

cc3[3:4]


cc4 <- comput(liderazgo, organizacion, cc3)
cc4[3:6]
cc4[3:6]$corr.X.xscores %>% View
cc4[3:6]$corr.Y.yscores %>% View

#x11()
#ggplot(data = data.frame(P7_2 = aprendizaje$P7_2, P7_5 = tecnologia$P7_5), aes(x = P7_2, y = P7_5, size = 28, alpha = .6)) + geom_point()

# Tests of canonical dimensions
ev <- (1 - cc3$cor^2)

n <- dim(organizacion)[1]
p <- length(organizacion)
q <- length(liderazgo)
k <- min(p, q)
m <- n - 3/2 - (p + q)/2

w <- rev(cumprod(rev(ev)))

# initialize
d1 <- d2 <- f <- vector("numeric", k)

for (i in 1:k) {
  s <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
  si <- 1/s
  d1[i] <- p * q
  d2[i] <- m * s - p * q/2 + 1
  r <- (1 - w[i]^si)/w[i]^si
  f[i] <- r * d2[i]/d1[i]
  p <- p - 1
  q <- q - 1
}

pv <- pf(f, d1, d2, lower.tail = FALSE)
(dmat <- cbind(WilksL = w, F = f, df1 = d1, df2 = d2, p = pv))

# standardized psych canonical coefficients diagonal matrix of apre sd's
s1 <- diag(sqrt(diag(cov(organizacion))))
s1 %*% cc3$xcoef
# standardized psych canonical coefficients diagonal matrix of tec sd's
s2 <- diag(sqrt(diag(cov(liderazgo))))
s2 %*% cc3$ycoef


#----------------------------------------------------------#
# Alfa de Cronbach 
#----------------------------------------------------------#

suppressMessages(library(nlme))
suppressMessages(library(MASS))
suppressMessages(library(multilevel))

km_data <- read.spss(file = "../_data/Base GConocimiento PymeS  Valle_2017.sav", to.data.frame = T, use.value.labels = F)
pilares <- km_data %>% select(P3_1:P3_8, P5_1:P9_6, P11_1:P13_16, P22_1:P25_5)

# Alfa aprendizaje
aprendizaje <- pilares %>% dplyr::select(P5_4, P5_5, P6_1, P6_2, P6_5, P7_2, P8_1, P13_13:P13_16)
reliability(cov(aprendizaje[,c("P5_4","P5_5","P6_1","P6_2","P6_5","P7_2",
                               "P8_1","P13_13","P13_14","P13_15","P13_16")], use="complete.obs"))

# Alfa tecnologia
tecnologia <- pilares %>% dplyr::select(P3_4, P5_6, P7_4, P7_5, P11_1:P11_5, P13_5, P13_6, P13_7, P13_8)
reliability(cov(tecnologia[,c("P3_4", "P5_6", "P7_4", "P7_5", "P11_1", "P11_2", "P11_3", "P11_4", "P11_5", "P13_5", "P13_6", "P13_7", "P13_8")], 
                use="complete.obs"))

#Alfa Liderazgo
liderazgo <- pilares %>% dplyr::select(P3_1:P3_3, P3_5:P3_8, P8_3, P8_4, P8_5, P13_1, P13_2, P13_3, P13_4)
reliability(cov(liderazgo[,c("P3_1", "P3_2", "P3_3", "P3_5", "P3_6", "P3_7", "P3_8", "P8_3", "P8_4", "P8_5", "P13_1", "P13_2", "P13_3", "P13_4")], 
                use="complete.obs"))

#Alfa Orgamizacion 
organizacion <- pilares %>% dplyr::select(P5_1:P5_3, P6_3, P6_4, P7_1, P7_3, P8_2, P13_9, P13_11)
reliability(cov(organizacion[,c("P5_1", "P5_2", "P5_3", "P6_3", "P6_4", "P7_1", "P7_3", "P8_2", "P13_9", "P13_11")], 
                use="complete.obs"))

#----------------------------------------------------------#
# Cuantificacion optima de variables
#----------------------------------------------------------#

library(homals)

### Aprendizaje
faData <- aprendizaje
faHomals <- homals(faData, active=rep(TRUE, ncol(faData)), level='ordinal', sets=list(1:ncol(faData)))
faQuan   <- faHomals$low.rank # Quantified variables
faDataqF <- faData

for(i in 1:ncol(faDataqF)){
  
  faDataqF[,i] <- as.character(faDataqF[,i])
  dfQuan <- as.data.frame(faQuan[[i]])
  
  for(j in 1:length(rownames(dfQuan))){
    mtch <- which(faDataqF[,i]==rownames(dfQuan)[j])
    if(length(mtch)>0){
      faDataqF[,i][which(faDataqF[,i]==rownames(dfQuan)[j])] <- dfQuan[j,1]
    }
  }
  
  faDataqF[,i] <- as.numeric(faDataqF[,i])
  
}; rm(i,j,mtch)
rm(faData, faHomals, dfQuan, faQuan)

x11()
corrplot(cor(cbind(faDataqF)), method = "number")

