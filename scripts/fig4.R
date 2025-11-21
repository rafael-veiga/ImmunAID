#####################################################################################################################
# BASE

library(ape)
library(WCluster)
library(ggExtra)
library(tidyverse)
library("ggpubr")
library(ggtree)
library(patchwork)
library(gghalves)
library(randomForest)
library(ggrepel)
library(ggtext)
library(extrafont)
library("RColorBrewer")
set.seed(12345)
#font_import(paths = "c:/windows/Fonts/")
#loadfonts(device = "win")
#windowsFonts()
font ="serif"

#####################################################################################################################
color_eu = c("Healthy" ="#66CD00",
             "Autoinflammation of unknown origin" = "#009ACD",
             "Still's disease"="#FF4500",
             "FMF"="#9A32CD",
             "Behcet"="#F5E400")
#####################################################
dis = list(var = c("Healthy","Autoinflammation of unknown origin","Still's disease","FMF","Behcet"),
           size = c(NA,40,82,55,61),
           color = color_eu[1:5],
           ref1 = c(NA,"CD38+ CD8+ T Cells in CD8+ T Cells",
                    "CD38+ CD8+ T Cells in CD8+ T Cells",
                    "NK Cells in Live Cells",
                    "NK Cells in Live Cells"),
           ref2 = c(NA,"IgM+ IgD+ Memory B Cells in Memory B Cells",
                    "BAFF-Rhigh Naive B Cells in Naive B Cells",
                    "CD40+ CD21- B Cells in CD21- B Cells",
                    "CD94+ CD16+ NK Cells in CD16+ NK Cells"),
           pos_refx1 = c(NA,0.5,0.4,1,0.5),
           pos_refy1 =c(NA,0,1.7,1,1),
           pos_refx2 = c(NA,0.8,0.4,0.9,0.25),
           pos_refy2 =c(NA,1.5,0,-0.2,0))
new_labels = function(vet){
  aux = sapply(as.character(vet), function(x) {
    x = gsub("\\+", "<sup>+</sup>", x)
    x = gsub("\\-", "<sup>-</sup>", x)
    x = gsub("Non<sup>-</sup>Switched", "Non-Switched", x)
    x = gsub("low", "<sup>low</sup>", x)
    x = gsub("high", "<sup>high</sup>", x)
    x = gsub("dim", "<sup>dim</sup>", x)
    x = gsub("bright", "<sup>bright</sup>", x)
    return(x)
  })
  return(aux)
}
get_peach_color <- function(value, min_val, max_val) {
  num_colors <- 9
  colors <- brewer.pal(num_colors, "YlOrRd")
  normalized_value <- (value - min_val) / (max_val - min_val)
  color_index <- round(normalized_value * (num_colors - 1)) + 1
  return(colors[color_index])
}

centroid_dist <- function(x, y, groups) {
  stopifnot(length(groups) == 2)
  # Filtrar apenas as instâncias dos dois grupos
  sel <- y %in% groups
  Xsub <- x[sel, , drop = FALSE]
  ysub <- droplevels(y[sel])
  # Centróides de cada grupo
  mu1 <- colMeans(Xsub[ysub == groups[1], , drop = FALSE])
  mu2 <- colMeans(Xsub[ysub == groups[2], , drop = FALSE])
  # Distância euclidiana entre os centróides
  dist <- sqrt(sum((mu1 - mu2)^2))
  return(dist)
}

plot_triangle_dist <- function(labels, d12, d13, d23,
                               colors = c("red","green","blue"),
                               edge_label_digits = 2,
                               point_size = 5, edge_size = 0.8,
                               dashed = TRUE, normalize = TRUE) {
  stopifnot(length(labels) == 3, length(colors) == 3)
  
  # normalizar distâncias se pedido
  if (normalize) {
    s <- d12 + d13 + d23
    d12 <- d12 / s
    d13 <- d13 / s
    d23 <- d23 / s
  }
  
  # checagem de triângulo
  eps <- 1e-9
  if ((d12 + d13 <= d23 + eps) || (d12 + d23 <= d13 + eps) || (d13 + d23 <= d12 + eps)) {
    stop("As distâncias não formam um triângulo (violam desigualdade do triângulo).")
  }
  
  # A=labels[1], B=labels[2], C=labels[3]
  a <- d23; b <- d13; c <- d12
  xC <- (b^2 + c^2 - a^2) / (2*c)
  yC <- sqrt(max(b^2 - xC^2, 0))
  
  nodes <- data.frame(
    label = labels,
    x = c(0, c, xC),
    y = c(0, 0, yC),
    color = colors
  )
  
  edges <- data.frame(
    x    = c(0,     0,     c),
    y    = c(0,     0,     0),
    xend = c(c,     xC,    xC),
    yend = c(0,     yC,    yC),
    dist = c(d12,   d13,   d23)
  )
  edges$xm <- (edges$x + edges$xend)/2
  edges$ym <- (edges$y + edges$yend)/2
  
  # padding
  pad <- 0.15 * max(c(d12, d13, d23))
  xlim <- range(c(nodes$x, nodes$xend)) + c(-pad, pad)
  ylim <- range(c(nodes$y, nodes$yend)) + c(-pad, pad)
  
  ggplot() +
    geom_segment(data = edges,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = edge_size,
                 linetype = if (dashed) "dashed" else "solid",
                 color = "grey40") +
    geom_point(data = nodes, aes(x = x, y = y), size = point_size) +
    geom_label(data = nodes, aes(x = x, y = y, label = label, fill = color),
               nudge_y = 0.03*max(c(d12,d13,d23)), label.size = 0.2, alpha = 0.9,
               color = "black") +
    geom_text(data = edges,
              aes(x = xm, y = ym, label = round(dist, edge_label_digits)),
              vjust = -0.6) +
    scale_fill_identity() +
    coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, clip = "off") +
    theme_void(base_size = 12) +    # remove eixos e grids
    theme(legend.position = "none")
}

get_protein_names <- function(ids, prot_df) {
  idx <- match(ids, prot_df$`Majority protein IDs`)
  prot_df$`Protein names`[idx]
}

##################################################################################
# Figure 4
## violin B
df = read.csv("./RF_imp2_f.csv")
res = as.character(df$vari)
res = res[res!="age"]
res = res[res!="sex"]
res = res[1:16]
data = readRDS(file = "./data.rds")
df = data$df_f
df$treated[df$disease=="Healthy"] = 0
df = df %>% 
  mutate(disease = case_when(
    disease=="Inflammation of unknown origin" ~ "Autoinflammation of unknown origin",
    disease=="AOSD" ~ "Still's disease",
    TRUE ~ disease
  ))
df = df[df$disease %in% dis$var,]
df$disease = factor(df$disease,levels = dis$var)
df = df %>% select(c("disease","treated",res))
df = gather(df, marks,value,res)
df = df[!is.na(df$value),]


df$marks = new_labels(df$marks)
df$marks = factor(df$marks, levels=new_labels(res))
aux = df
fig_3 <- aux %>% 
  ggplot(aes(x = disease, y = value, fill = disease)) +
  facet_wrap(~ marks, scales = "free_y", ncol = 2) +
  geom_half_violin(
    data = aux %>% filter(treated == 0), side = "l",
    trim = FALSE, scale = "width", width = 1
  ) +
  geom_half_violin(
    data = aux %>% filter(treated == 1), side = "r",
    trim = FALSE, scale = "width", width = 1
  ) +
  geom_jitter(
    aes(shape = factor(treated), group = factor(treated)),
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.6, jitter.height = 0),
    size = 1
  ) +
  # barra de média dashed e fina
  stat_summary(
    aes(group = factor(treated)),
    fun      = mean, fun.min = mean, fun.max = mean,
    geom     = "crossbar", width = 0.3,
    position = position_dodge(width = 0.8),
    color    = "black",
    linetype = "dashed",
    size     = 0.3
  ) +
  scale_shape_manual(
    values = c("0" = 19, "1" = 4),
    labels = c("0" = "untreated", "1" = "treated"),
    name   = "Treated"
  ) +
  scale_fill_manual(values = dis$color, limits = dis$var) +
  scale_x_discrete(limits = dis$var) +
  #scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Quantification") +
  theme_bw(base_family = font) +
  theme(
    axis.title.x      = element_blank(),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(size = 8, family = font),
    axis.title.y      = element_text(size = 10, family = font),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 11, family = font),
    strip.text        = ggtext::element_markdown(size = 10, family = font),
    strip.background  = element_rect(fill = "bisque"),
    legend.position   = "top",
    legend.justification = "left",
    legend.direction  = "horizontal",
    plot.margin       = unit(c(0,2,2,2), "mm")
  )

pdf("./fig4.pdf",width = 9,height = 12)
print(fig_3)
dev.off()

##################################################################################
