## author:wjw
## description:annotation_DMR_gene_data. 
## time:20250127

library(tidyverse)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("org.Hs.eg.db")
library(xlsx)
library(colorspace)
library(ggvenn)
library(RColorBrewer)
library("ggsci")
library(tidyverse)
library("clusterProfiler")
library(ggsankey)
library(ggplot2)
install.packages("cols4all")
library(cols4all)
library(ggnewscale)
library(ggfun)
library(grid)
library(ggh4x)

rm(list = ls())

load("Rdata/mydata_DMP_DMR.Rdata")
myDMR = myDMR$BumphunterDMR
rm(myDMP)

myDMR$group = ifelse(myDMR$value > 0, "Hypermethylated", "Hypomethylated")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

write.table(myDMR,file = "myDMR.txt",row.names = F,sep = "\t",quote = F)

myDMR <- annotatePeak("./myDMR.txt", tssRegion = c(-3000, 3000), TxDb = txdb,annoDb = 'org.Hs.eg.db')

myDMR = as.data.frame(myDMR@anno@elementMetadata@listData)

myDMR %>% 
  select(SYMBOL) %>% 
  distinct()

myDMR %>%
  count(SYMBOL) %>%
  filter(n > 1) 

myDMR_anno = myDMR

##
geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=myDMR_anno$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 20
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:20)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:20) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_1.tiff", 
     width = 10 * 300,
     height = 15 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:40,
                                   labels = plot_df$Description)) +
  ggtitle(label = "GO and KEGG annotation") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
                    )
dev.off()

##CHR1

GENES = myDMR_anno %>% 
  filter(geneChr == 1) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:6) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_2.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:12,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr1") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

###chr2

GENES = myDMR_anno %>% 
  filter(geneChr == 2) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:6) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_3.tiff", 
     width = 9 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:12,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr2") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

##chr3

GENES = myDMR_anno %>% 
  filter(geneChr == 3) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:6) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_4.tiff", 
     width = 9 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:12,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr3") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(0.5,0),limits = c(1, 3), breaks = seq(1, 3, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr4

GENES = myDMR_anno %>% 
  filter(geneChr == 4) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:4) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_5.tiff", 
     width = 9 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:10,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr4") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(0.5,0),limits = c(1, 3), breaks = seq(1, 3, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr5

GENES = myDMR_anno %>% 
  filter(geneChr == 5) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:2) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_6.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:8,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr5") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(0.5,0),limits = c(1, 3), breaks = seq(1, 3, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr6
GENES = myDMR_anno %>% 
  filter(geneChr == 6) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:6) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_7.tiff", 
     width = 10 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:12,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr6") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(0.5,0),limits = c(2, 3), breaks = seq(2, 3, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr7
GENES = myDMR_anno %>% 
  filter(geneChr == 7) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:4) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_8.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:10,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr7") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr8
GENES = myDMR_anno %>% 
  filter(geneChr == 8) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:6) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_9.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:12,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr8") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr9
GENES = myDMR_anno %>% 
  filter(geneChr == 9) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:6) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- GO_top8 %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_10.tiff", 
     width = 8 * 300,
     height = 6 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:6,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr9") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = "#CC0000"),
    ggh4x.axis.nesttext.y = element_text(colour = "#CC0000"),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 1), breaks = seq(1, 1, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr10
GENES = myDMR_anno %>% 
  filter(geneChr == 10) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:3) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_11.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:9,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr10") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr11
GENES = myDMR_anno %>% 
  filter(geneChr == 11) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:6) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_12.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:12,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr11") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr12
GENES = myDMR_anno %>% 
  filter(geneChr == 12) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:1) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_13.tiff", 
     width = 8 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:7,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr12") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr13
GENES = myDMR_anno %>% 
  filter(geneChr == 13) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:4) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_14.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:10,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr13") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr14
GENES = myDMR_anno %>% 
  filter(geneChr == 14) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:3) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_15.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:9,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr14") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr15
GENES = myDMR_anno %>% 
  filter(geneChr == 15) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:3) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_16.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:9,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr15") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 1), breaks = seq(1, 1, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr16
GENES = myDMR_anno %>% 
  filter(geneChr == 16) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:1) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_17.tiff", 
     width = 9 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:7,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr16") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(0.5,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr17
GENES = myDMR_anno %>% 
  filter(geneChr == 17) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:5) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_18.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:11,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr17") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr18
GENES = myDMR_anno %>% 
  filter(geneChr == 18) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:1) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_19.tiff", 
     width = 8 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:7,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr18") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 1), breaks = seq(1, 1, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr19
GENES = myDMR_anno %>% 
  filter(geneChr == 19) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:2) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_20.tiff", 
     width = 8 * 300,
     height = 9 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:8,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr19") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr20
GENES = myDMR_anno %>% 
  filter(geneChr == 20) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:6) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_21.tiff", 
     width = 8 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:12,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr20") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 1), breaks = seq(1, 1, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()

### chr22
GENES = myDMR_anno %>% 
  filter(geneChr == 22) %>% 
  dplyr::select(SYMBOL)

geneIDselect <-AnnotationDbi::select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                                     keys=GENES$SYMBOL,
                                     columns=c("ENTREZID"), 
                                     keytype="SYMBOL" )

go_kk <- enrichGO(gene = geneIDselect$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =1, 
                  qvalueCutoff = 1,
                  ont="all",
                  readable =T)

deg_go_2 <- setReadable(go_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

deg_kegg_2 <- setReadable(kegg_kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame()

# top 6
GO_top8 <- deg_go_2 %>%
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::slice(1:6)

KEGG_top8 <- deg_kegg_2 %>%
  dplyr::slice(1:2) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

plot_df <- rbind(GO_top8, KEGG_top8)%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "KEGG")), ordered = T)) %>%
  dplyr::arrange(ONTOLOGY, desc(Count)) %>%
  #dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered =T))

tiff(file = "./Figure/step2_22.tiff", 
     width = 8 * 300,
     height = 9 * 300, 
     res = 300, 
     compression = "lzw")
plot_df %>%
  ggplot() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#CCCCFF", high = "#0000CC", name = "KEGG pvalue", guide = guide_colorbar(order = 2)) +
  ggnewscale::new_scale_fill() + 
  geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) +
  scale_fill_gradient(low = "#FFCCCC", high = "#CC0000", name = "BP pvalue", guide = guide_colorbar(order = 1)) + 
  guides(y = "axis_nested",
         y.sec = guide_axis_manual(breaks = 1:8,
                                   labels = plot_df$Description)) +
  ggtitle(label = "Chr22") +
  labs(x = "Count", y = "Description") +
  scale_size(range = c(3, 7),
             guide = guide_legend(override.aes = list(fill = "#000000"), order = 3)) +
  theme_bw() +
  theme(
    ggh4x.axis.nestline.y = element_line(size = 3, color = c("#0000CC", "#CC0000")),
    ggh4x.axis.nesttext.y = element_text(colour = c("#0000CC", "#CC0000")),
    legend.background = element_roundrect(color = "#969696"),
    panel.border = element_rect(size = 0.5),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    axis.text = element_text(color = "#000000", size = 11),
    #axis.text.y = element_text(color = rep(c("#4DB6AC", "#FDD835"), each = 8)),
    axis.text.y.left = element_blank(),
    axis.ticks.length.y.left = unit(8, "pt"),
    axis.ticks.y.left = element_line(color = NA),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(1,0),limits = c(1, 2), breaks = seq(1, 2, by = 1))+
  coord_cartesian(clip = "off") +
  annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                         gp = gpar(col = "#969696", lwd = 1.5)),
                    xmin = unit(1, "native"),
                    xmax = unit(12, "native"),
                    ymin = unit(40.85, "native"),
                    ymax = unit(41.75, "native")
  )
dev.off()
