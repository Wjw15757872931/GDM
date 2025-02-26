
rm(list = ls())

library(EpiDISH)
library(ggpubr)
library(viridis)
library(plotly)
library(tidyverse)
library(ggsci)

load("Rdata/mydata_beta_pd.Rdata")

BloodFrac.m <- epidish(beta.m = mydata_beta, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
boxplot(BloodFrac.m)

hub_cg = c("cg19037167", "cg22998811", "cg02249039", "cg15150348", "cg17537719", "cg26538349", "cg26792694")

mydata_beta = mydata_beta[hub_cg,]
mydata_beta = as.data.frame(t(mydata_beta))

mydata_immune = BloodFrac.m %>% 
  as.data.frame()

all_data = cbind(mydata_immune, mydata_beta)
all_data$group = mydata_pd$Sample_Group

y_positions <- all_data[c(1:7,15)] %>%
  pivot_longer(cols = -group, names_to = "Immune_Cells", values_to = "Proportion") %>%
  group_by(Immune_Cells) %>%
  summarise(max_value = max(Proportion) * 1.3) %>% # 增加偏移量到 1.3 倍
  pull(max_value)

color_labels <- function(x) {
  ifelse(x %in% c("CD4T"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

tiff(file = "./Figure/step4_3.tiff", 
     width = 10 * 300,   # 20 英寸 × 300 DPI = 6000 像素
     height = 8 * 300,  # 假设高度按比例为 12 英寸
     res = 300, 
     compression = "lzw")

plot.format=theme(plot.background=element_blank(),
                  panel.grid=element_blank(),
                  panel.background=element_blank(),
                  panel.border=element_rect(color="black",linewidth=0.5,fill=NA),
                  axis.line=element_blank(),
                  axis.ticks=element_line(color="black",linewidth=0.5),
                  axis.text=element_text(color="black",size=15),
                  axis.title=element_text(color="black",size=15),
                  plot.title=element_text(color="black",size=30, face = "bold"),
                  legend.background=element_blank(),
                  legend.key=element_blank(),
                  legend.text=element_text(color="black",size=15),
                  legend.title=element_text(color="black",size=15)
)

all_data[c(1:7,15)] %>% 
  pivot_longer(cols = -group, names_to = "Immune_Cells", values_to = "Proportion") %>% 
  ggplot(aes(x = Immune_Cells, y = Proportion)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4) +
  stat_compare_means(aes(group = group),
                     method = "t.test",
                     label = "p.signif",
                     label.y = y_positions) + # 使用更高的 y 位置
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  scale_fill_cosmic() +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  ggtitle("hnsfy_methydata") +
  labs(y = "Score") +
  plot.format
dev.off()

data_long <- all_data %>%
  select(c(hub_cg, 'Eosino')) %>%
  pivot_longer(cols = all_of(hub_cg), names_to = 'hub_cg', values_to = 'value')

# Create the plot
tiff(file = "./Figure/step4_4_Eosino.tiff", 
     width = 27 * 300,
     height = 4 * 300, 
     res = 300, 
     compression = "lzw")
ggplot(data_long, aes(x = value, y = Eosino)) +
  geom_point(size = 2, color = '#EC0101', alpha = 0.5) +
  theme_bw() +  
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()) +
  geom_smooth(method = 'lm', se = T, color = '#F9B208', size = 1.5, fill = '#FEA82F') +
  stat_cor(method = "pearson", digits = 3, size = 6) +
  facet_wrap(~hub_cg, scales = 'free', ncol = 7)+
  #ylab("Eosino") +
  xlab("Methylation value")
dev.off()

###
rm(list = ls())

load("Rdata/mydata_beta_pd.Rdata")

BloodFrac.m <- epidish(beta.m = mydata_beta, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
boxplot(BloodFrac.m)

mydata_immune = BloodFrac.m %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  arrange(rowname) %>% 
  mutate(Sample_Name = rowname) %>% 
  inner_join(mydata_pd,by="Sample_Name") %>% 
  dplyr::select(Sample_Name,B,NK,CD4T,CD8T,Mono,Neutro,Eosino,Sample_Group) %>% 
  pivot_longer(
    cols = B:Eosino,
    names_to = "class",
    values_to = "ratio"
  ) 

tiff(file = "./Figure/step4_5.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
ggplot(mydata_immune,aes(Sample_Name,ratio,fill = class)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.key.size = unit(15, "pt")
  ) + 
  scale_y_continuous(expand = c(0.01,0)) +
  guides(fill=guide_legend(nrow = 1,reverse = T,keywidth = 0))+
  scale_fill_aaas()
dev.off()

### hnsfy_methydata 
rm(list = ls())
load("Rdata/GSE153219_step3_combat_norm_data.Rdata")
BloodFrac.m_GSE153219 <- epidish(beta.m = myCombat, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
BloodFrac.m_GSE153219 = as.data.frame(BloodFrac.m_GSE153219)
BloodFrac.m_GSE153219$group = "GSE153219"

load("Rdata/mydata_beta_pd.Rdata")
BloodFrac.m_hnsfy_methydata <- epidish(beta.m = mydata_beta, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
BloodFrac.m_hnsfy_methydata = as.data.frame(BloodFrac.m_hnsfy_methydata)
BloodFrac.m_hnsfy_methydata$group = "hnsfy_methydata"

alldat = rbind(BloodFrac.m_hnsfy_methydata, BloodFrac.m_GSE153219)

y_positions <- alldat %>%
  pivot_longer(cols = -group, names_to = "Immune_Cells", values_to = "Proportion") %>%
  group_by(Immune_Cells) %>%
  summarise(max_value = max(Proportion) * 1.3) %>% # 增加偏移量到 1.3 倍
  pull(max_value)

color_labels <- function(x) {
  ifelse(x %in% c("CD4T", "B", "CD8T", "Eosino", "Mono", "Neutro", "NK"), 
         paste0('<span style="color:red">', x, '</span>'), 
         x)
}

plot.format=theme(plot.background=element_blank(),
                  panel.grid=element_blank(),
                  panel.background=element_blank(),
                  panel.border=element_rect(color="black",linewidth=0.5,fill=NA),
                  axis.line=element_blank(),
                  axis.ticks=element_line(color="black",linewidth=0.5),
                  axis.text=element_text(color="black",size=15),
                  axis.title=element_text(color="black",size=15),
                  plot.title=element_text(color="black",size=30, face = "bold"),
                  legend.background=element_blank(),
                  legend.key=element_blank(),
                  legend.text=element_text(color="black",size=15),
                  legend.title=element_text(color="black",size=15)
)

tiff(file = "./Figure/step4_1.tiff", 
     width = 10 * 300,   # 20 英寸 × 300 DPI = 6000 像素
     height = 8 * 300,  # 假设高度按比例为 12 英寸
     res = 300, 
     compression = "lzw")
alldat %>% 
  pivot_longer(cols = -group, names_to = "Immune_Cells", values_to = "Proportion") %>% 
  ggplot(aes(x = Immune_Cells, y = Proportion)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4) +
  stat_compare_means(aes(group = group),
                     method = "t.test",
                     label = "p.signif",
                     label.y = y_positions) + # 使用更高的 y 位置
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  scale_fill_cosmic() +
  scale_x_discrete(labels = color_labels) +
  theme(axis.text.x = element_markdown()) +
  #ggtitle("CIBERSORT") +
  labs(y = "Score") +
  plot.format
dev.off()

###
rm(list = ls())
load("Rdata/GSE153219_step3_combat_norm_data.Rdata")
BloodFrac.m_GSE153219 <- epidish(beta.m = myCombat, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
BloodFrac.m_GSE153219 = as.data.frame(BloodFrac.m_GSE153219)

BloodFrac.m_GSE153219$group = myLoad$pd$Sample_Group
BloodFrac.m_GSE153219$group = ifelse(BloodFrac.m_GSE153219$group == "NGT", "Control", "GDM")

alldat = BloodFrac.m_GSE153219

y_positions <- alldat %>%
  pivot_longer(cols = -group, names_to = "Immune_Cells", values_to = "Proportion") %>%
  group_by(Immune_Cells) %>%
  summarise(max_value = max(Proportion) * 1.3) %>% # 增加偏移量到 1.3 倍
  pull(max_value)

plot.format=theme(plot.background=element_blank(),
                  panel.grid=element_blank(),
                  panel.background=element_blank(),
                  panel.border=element_rect(color="black",linewidth=0.5,fill=NA),
                  axis.line=element_blank(),
                  axis.ticks=element_line(color="black",linewidth=0.5),
                  axis.text=element_text(color="black",size=15),
                  axis.title=element_text(color="black",size=15),
                  plot.title=element_text(color="black",size=30, face = "bold"),
                  legend.background=element_blank(),
                  legend.key=element_blank(),
                  legend.text=element_text(color="black",size=15),
                  legend.title=element_text(color="black",size=15)
)

tiff(file = "./Figure/step4_2.tiff", 
     width = 10 * 300,   # 20 英寸 × 300 DPI = 6000 像素
     height = 8 * 300,  # 假设高度按比例为 12 英寸
     res = 300, 
     compression = "lzw")
alldat %>% 
  pivot_longer(cols = -group, names_to = "Immune_Cells", values_to = "Proportion") %>% 
  ggplot(aes(x = Immune_Cells, y = Proportion)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.4) +
  stat_compare_means(aes(group = group),
                     method = "t.test",
                     label = "p.signif",
                     label.y = y_positions) + # 使用更高的 y 位置
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  scale_fill_cosmic() +
  #scale_x_discrete(labels = color_labels) +
  #theme(axis.text.x = element_markdown()) +
  ggtitle("GSE153219") +
  labs(y = "Score") +
  plot.format
dev.off()

####
rm(list = ls())

activation_scores = c('CD40LG','CD69','HLA-DPA1','HLA-DPB1','HLA-DQB1','HLA-DRB5','ICOS',
                      'TNFRSF4','NFAT5','NFATC3','NFATC2','LCP2','LAT','ZAP70','LCK')
cytokines_scores = c('IFNG','TNF','IL2','CSF2','IL10','IL4','IL5','IL13','IL17A','IL21',
                     'IL22','TGFB1','IL1B','IL6','EBI3','IFNB1','IFNA1')
resident_scores = c('CXCR3','CXCR6','ICAM1','ITGA1','ITGAE','RBPJ','ZEB2','ZNF683')
IFN_scores = c('IFIH1','IFIT1','IFIT2','IFIT3','FOS','JUN','JUND','REL','IRF1','IRF7',
               'ISG15','MX1','STAT1','NFKBIA','TRIM21','OAS1','OAS2')
inhibitory_scores = c('CTLA4','PDCD1','LAG3','HAVCR2','BTLA','TIGIT')

genes = c(activation_scores, cytokines_scores, resident_scores, IFN_scores, inhibitory_scores)

library(jjAnno)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(stringr)
library(tidyverse)
library("RColorBrewer")
library(dendextend)
library(reshape2)

load("./Rdata/mydata_DMP_DMR.Rdata")
myDMP = myDMP$Control_to_GDM
rm(myDMR)

dmp1_up = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% activation_scores) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  filter(logFC > 0 )

dmp1_down = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% activation_scores) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  filter(logFC < 0 )

rm(myDMP)

load("Rdata/mydata_beta_pd.Rdata")

mydata_beta = as.data.frame(mydata_beta)

data1 = mydata_beta %>% 
  rownames_to_column("cgname") %>% 
  as_tibble() %>% 
  inner_join(dmp1_down, by = "cgname") %>% 
  as.data.frame()

data2 = mydata_beta %>% 
  rownames_to_column("cgname") %>% 
  as_tibble() %>% 
  inner_join(dmp1_up, by = "cgname") %>% 
  as.data.frame()

data = rbind(data1, data2)
rownames(data) = data$gene
data = data[-c(1,32:51)]

GDM_beta_pd = mydata_pd %>% 
  as_tibble() %>% 
  dplyr::filter(Sample_Group == "GDM")

Control_beta_pd = mydata_pd %>% 
  as_tibble() %>% 
  dplyr::filter(Sample_Group == "Control")

data = cbind(data[GDM_beta_pd$Sample_Name], data[Control_beta_pd$Sample_Name])

#对数据进行归一化；
#由于scale函数默认对列进行归一化，因此这里做了两次转置；
beta_scaled <- t(scale(t(data)))

df_long <- melt(beta_scaled)
colnames(df_long) <- c("Gene", "Sample", "Expression")

mycols <- scale_fill_gradient2(low = "#2F509EFF", mid = "white" ,high = "#CF4E9CFF", midpoint = 0)

p = ggplot(df_long, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile() +
  mycols +  # 应用颜色渐变
  theme_minimal()+
  coord_cartesian(clip = 'off')+
  scale_y_discrete(position = "right") +
  theme(
    axis.text.x = element_blank(),  # 去除x轴的刻度标签
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.title.y = element_blank(),
    plot.margin = margin(4, 4, 4, 4, "cm"),
   # legend.position = "none"
  )

p1 = annoRect(object = p,
         annoPos = 'top',
         annoManual = T,
         xPosition = list(c(1,21),
                        c(20,30)),
         yPosition = c(10,10.5),
         lty = 'solid',
         lwd = 2,
         addText = T,
         textLabel = c("GDM", "Control"),
         textHVjust = 0.5,
         textCol = rep('black',2)
         )

tiff(file = "./Figure/step4_6.tiff", 
     width = 12 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
p2 = annoRect(object = p1,
              annoPos = 'left',
              annoManual = T,
              xPosition = c(-1.5,-1),
              yPosition = list(c(1),c(9)),
              addText = T,
              textLabel = "activation",
              pFill = '#F5F0BB',
              pCol = '#F5F0BB',
              textHVjust = -1,
              textCol = 'black',
              textRot = 90,
              textSize = 15
)
dev.off()

###

load("./Rdata/mydata_DMP_DMR.Rdata")
myDMP = myDMP$Control_to_GDM
rm(myDMR)

dmp1_up = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% cytokines_scores) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  filter(logFC > 0 )

dmp1_down = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% cytokines_scores) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  filter(logFC < 0 )

rm(myDMP)

load("Rdata/mydata_beta_pd.Rdata")

mydata_beta = as.data.frame(mydata_beta)

data1 = mydata_beta %>% 
  rownames_to_column("cgname") %>% 
  as_tibble() %>% 
  inner_join(dmp1_down, by = "cgname") %>% 
  as.data.frame()

data2 = mydata_beta %>% 
  rownames_to_column("cgname") %>% 
  as_tibble() %>% 
  inner_join(dmp1_up, by = "cgname") %>% 
  as.data.frame()

data = rbind(data1, data2)
rownames(data) = data$gene
data = data[-c(1,32:51)]

GDM_beta_pd = mydata_pd %>% 
  as_tibble() %>% 
  dplyr::filter(Sample_Group == "GDM")

Control_beta_pd = mydata_pd %>% 
  as_tibble() %>% 
  dplyr::filter(Sample_Group == "Control")

data = cbind(data[GDM_beta_pd$Sample_Name], data[Control_beta_pd$Sample_Name])

beta_scaled <- t(scale(t(data)))

df_long <- melt(beta_scaled)
colnames(df_long) <- c("Gene", "Sample", "Expression")

mycols <- scale_fill_gradient2(low = "#2F509EFF", mid = "white" ,high = "#CF4E9CFF", midpoint = 0)

p = ggplot(df_long, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile() +
  mycols +  # 应用颜色渐变
  theme_minimal()+
  coord_cartesian(clip = 'off')+
  scale_y_discrete(position = "right") +
  theme(
    axis.text.x = element_blank(),  # 去除x轴的刻度标签
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.title.y = element_blank(),
    plot.margin = margin(4, 4, 4, 4, "cm"),
    legend.position = "none"
  )

tiff(file = "./Figure/step4_7.tiff", 
     width = 10 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
p2 = annoRect(object = p,
              annoPos = 'left',
              annoManual = T,
              xPosition = c(-1.5,-1),
              yPosition = list(c(1),c(10)),
              addText = T,
              textLabel = "cytokines",
              pFill = '#F37B7D',
              pCol = '#F37B7D',
              textHVjust = -1,
              textCol = 'black',
              textRot = 90,
              textSize = 15
)
dev.off()

##

load("./Rdata/mydata_DMP_DMR.Rdata")
myDMP = myDMP$Control_to_GDM
rm(myDMR)

dmp1_up = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% resident_scores) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  filter(logFC > 0 )

dmp1_down = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% resident_scores) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  filter(logFC < 0 )

rm(myDMP)

load("Rdata/mydata_beta_pd.Rdata")

mydata_beta = as.data.frame(mydata_beta)

data1 = mydata_beta %>% 
  rownames_to_column("cgname") %>% 
  as_tibble() %>% 
  inner_join(dmp1_down, by = "cgname") %>% 
  as.data.frame()

data2 = mydata_beta %>% 
  rownames_to_column("cgname") %>% 
  as_tibble() %>% 
  inner_join(dmp1_up, by = "cgname") %>% 
  as.data.frame()

data = rbind(data1, data2)
rownames(data) = data$gene
data = data[-c(1,32:51)]

GDM_beta_pd = mydata_pd %>% 
  as_tibble() %>% 
  dplyr::filter(Sample_Group == "GDM")

Control_beta_pd = mydata_pd %>% 
  as_tibble() %>% 
  dplyr::filter(Sample_Group == "Control")

data = cbind(data[GDM_beta_pd$Sample_Name], data[Control_beta_pd$Sample_Name])

beta_scaled <- t(scale(t(data)))

df_long <- melt(beta_scaled)
colnames(df_long) <- c("Gene", "Sample", "Expression")

mycols <- scale_fill_gradient2(low = "#2F509EFF", mid = "white" ,high = "#CF4E9CFF", midpoint = 0)

p = ggplot(df_long, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile() +
  mycols +  # 应用颜色渐变
  theme_minimal()+
  coord_cartesian(clip = 'off')+
  scale_y_discrete(position = "right") +
  theme(
    axis.text.x = element_blank(),  # 去除x轴的刻度标签
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.title.y = element_blank(),
    plot.margin = margin(4, 4, 4, 4, "cm"),
    legend.position = "none"
  )

tiff(file = "./Figure/step4_8.tiff", 
     width = 10 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
p2 = annoRect(object = p,
              annoPos = 'left',
              annoManual = T,
              xPosition = c(-1.5,-1),
              yPosition = list(c(1),c(5)),
              addText = T,
              textLabel = "resident",
              pFill = '#89C75F',
              pCol = '#89C75F',
              textHVjust = -1,
              textCol = 'black',
              textRot = 90,
              textSize = 15
)
dev.off()

##
load("./Rdata/mydata_DMP_DMR.Rdata")
myDMP = myDMP$Control_to_GDM
rm(myDMR)

dmp1_up = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% IFN_scores) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  filter(logFC > 0 )

dmp1_down = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% IFN_scores) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  filter(logFC < 0 )

rm(myDMP)

load("Rdata/mydata_beta_pd.Rdata")

mydata_beta = as.data.frame(mydata_beta)

data1 = mydata_beta %>% 
  rownames_to_column("cgname") %>% 
  as_tibble() %>% 
  inner_join(dmp1_down, by = "cgname") %>% 
  as.data.frame()

data2 = mydata_beta %>% 
  rownames_to_column("cgname") %>% 
  as_tibble() %>% 
  inner_join(dmp1_up, by = "cgname") %>% 
  as.data.frame()

data = rbind(data1, data2)
rownames(data) = data$gene
data = data[-c(1,32:51)]

GDM_beta_pd = mydata_pd %>% 
  as_tibble() %>% 
  dplyr::filter(Sample_Group == "GDM")

Control_beta_pd = mydata_pd %>% 
  as_tibble() %>% 
  dplyr::filter(Sample_Group == "Control")

data = cbind(data[GDM_beta_pd$Sample_Name], data[Control_beta_pd$Sample_Name])

beta_scaled <- t(scale(t(data)))

df_long <- melt(beta_scaled)
colnames(df_long) <- c("Gene", "Sample", "Expression")

mycols <- scale_fill_gradient2(low = "#2F509EFF", mid = "white" ,high = "#CF4E9CFF", midpoint = 0)

p = ggplot(df_long, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile() +
  mycols +  # 应用颜色渐变
  theme_minimal()+
  coord_cartesian(clip = 'off')+
  scale_y_discrete(position = "right") +
  theme(
    axis.text.x = element_blank(),  # 去除x轴的刻度标签
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.title.y = element_blank(),
    plot.margin = margin(4, 4, 4, 4, "cm"),
    legend.position = "none"
  )

tiff(file = "./Figure/step4_9.tiff", 
     width = 10 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
p2 = annoRect(object = p,
              annoPos = 'left',
              annoManual = T,
              xPosition = c(-1.5,-1),
              yPosition = list(c(1),c(12)),
              addText = T,
              textLabel = "IFN",
              pFill = '#90D5E4',
              pCol = '#90D5E4',
              textHVjust = -1,
              textCol = 'black',
              textRot = 90,
              textSize = 15
)
dev.off()

###
load("./Rdata/mydata_DMP_DMR.Rdata")
myDMP = myDMP$Control_to_GDM
rm(myDMR)

dmp1_up = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% inhibitory_scores) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  filter(logFC > 0 )

dmp1_down = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(cgname = rowname) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(gene %in% inhibitory_scores) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  filter(logFC < 0 )

rm(myDMP)

load("Rdata/mydata_beta_pd.Rdata")

mydata_beta = as.data.frame(mydata_beta)

data1 = mydata_beta %>% 
  rownames_to_column("cgname") %>% 
  as_tibble() %>% 
  inner_join(dmp1_down, by = "cgname") %>% 
  as.data.frame()

data2 = mydata_beta %>% 
  rownames_to_column("cgname") %>% 
  as_tibble() %>% 
  inner_join(dmp1_up, by = "cgname") %>% 
  as.data.frame()

data = rbind(data1, data2)
rownames(data) = data$gene
data = data[-c(1,32:51)]

GDM_beta_pd = mydata_pd %>% 
  as_tibble() %>% 
  dplyr::filter(Sample_Group == "GDM")

Control_beta_pd = mydata_pd %>% 
  as_tibble() %>% 
  dplyr::filter(Sample_Group == "Control")

data = cbind(data[GDM_beta_pd$Sample_Name], data[Control_beta_pd$Sample_Name])

beta_scaled <- t(scale(t(data)))

df_long <- melt(beta_scaled)
colnames(df_long) <- c("Gene", "Sample", "Expression")

mycols <- scale_fill_gradient2(low = "#2F509EFF", mid = "white" ,high = "#CF4E9CFF", midpoint = 0)

p = ggplot(df_long, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile() +
  mycols +  # 应用颜色渐变
  theme_minimal()+
  coord_cartesian(clip = 'off')+
  scale_y_discrete(position = "right") +
  theme(
    axis.text.x = element_blank(),  # 去除x轴的刻度标签
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.title.y = element_blank(),
    plot.margin = margin(4, 4, 4, 4, "cm"),
    legend.position = "none"
  )

tiff(file = "./Figure/step4_10.tiff", 
     width = 10 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
p2 = annoRect(object = p,
              annoPos = 'left',
              annoManual = T,
              xPosition = c(-1.5,-1),
              yPosition = list(c(1),c(6)),
              addText = T,
              textLabel = "inhibitory",
              pFill = '#7E1416',
              pCol = '#7E1416',
              textHVjust = -1,
              textCol = 'black',
              textRot = 90,
              textSize = 15
)
dev.off()

###

hub_cg = c("cg19037167", "cg22998811", "cg02249039", "cg15150348", "cg17537719", "cg26538349", "cg26792694")

load("./Rdata/mydata_DMP_DMR.Rdata")
myDMP = myDMP$Control_to_GDM

hub_cg_info = myDMP %>% 
  rownames_to_column("cgname") %>% 
  filter(cgname %in% hub_cg)

write.xlsx(hub_cg_info, file = "./Figure/hub_cg_info.xlsx")
