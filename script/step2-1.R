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
library(cols4all)
library(ggnewscale)
library(ggfun)
library(grid)
library(ggh4x)
library(DOSE)
library(patchwork)

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

write.xlsx(GO_top8, file = "./Figure/GO_top8.xlsx")
write.xlsx(KEGG_top8, file = "./Figure/KEGG_top8.xlsx")

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
     height = 25 * 300, 
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
  coord_cartesian(clip = "off")
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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p1 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("GO Chr1")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k1 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("KEGG Chr1")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p2 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr2")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k2 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr2")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p3 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr3")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k3 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr3")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p4 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr4")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k4 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr4")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p5 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr5")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k5 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr5")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p6 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr6")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k6 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr6")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p7 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr7")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k7 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr7")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p8 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr8")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k8 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr8")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p9 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr9")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k9 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr9")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p10 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr10")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k10 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr10")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p11 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr11")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k11 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr11")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p12 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr12")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k12 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr12")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p13 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr13")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k13 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr13")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p14 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr14")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k14 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr14")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p15 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr15")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k15 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr15")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p16 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr16")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k16 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr16")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p17 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr17")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k17 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr17")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p18 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr18")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k18 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr18")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p19 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr19")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k19 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr19")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p20 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr20")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k20 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr20")

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

kegg_kk <- enrichKEGG(gene =geneIDselect$ENTREZID,
                      organism = 'hsa', #KEGG可以用organism = 'hsa'
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)

p22 = cnetplot(go_kk, showCategory = 6,foldChange=geneIDselect$ENTREZID, circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr22")

kegg_kk_symbol <- setReadable(
  kegg_kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"  # 原始ID类型
)

k22 = cnetplot(kegg_kk_symbol, showCategory = 6,
              foldChange=geneIDselect$SYMBOL, 
              circular = TRUE, colorEdge = TRUE) +
  ggtitle("Chr22")

tiff(file = "./Figure/step2_ch1_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p1)
dev.off()

tiff(file = "./Figure/step2_ch2_go.tiff", 
     width = 12 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p2)
dev.off()

tiff(file = "./Figure/step2_ch3_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p3)
dev.off()

tiff(file = "./Figure/step2_ch4_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p4)
dev.off()

tiff(file = "./Figure/step2_ch5_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p5)
dev.off()

tiff(file = "./Figure/step2_ch6_go.tiff", 
     width = 14 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p6)
dev.off()

tiff(file = "./Figure/step2_ch7_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p7)
dev.off()

tiff(file = "./Figure/step2_ch8_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p8)
dev.off()

tiff(file = "./Figure/step2_ch9_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p9)
dev.off()

tiff(file = "./Figure/step2_ch10_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p10)
dev.off()

tiff(file = "./Figure/step2_ch11_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p11)
dev.off()

tiff(file = "./Figure/step2_ch12_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p12)
dev.off()

tiff(file = "./Figure/step2_ch13_go.tiff", 
     width = 14 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p13)
dev.off()

tiff(file = "./Figure/step2_ch14_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p14)
dev.off()

tiff(file = "./Figure/step2_ch15_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p15)
dev.off()

tiff(file = "./Figure/step2_ch16_go.tiff", 
     width = 12 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p16)
dev.off()

tiff(file = "./Figure/step2_ch17_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p17)
dev.off()

tiff(file = "./Figure/step2_ch18_go.tiff", 
     width = 12 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p18)
dev.off()

tiff(file = "./Figure/step2_ch19_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p19)
dev.off()

tiff(file = "./Figure/step2_ch20_go.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p20)
dev.off()

tiff(file = "./Figure/step2_ch22_go.tiff", 
     width = 12 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(p22)
dev.off()

###
tiff(file = "./Figure/step2_ch1_kegg.tiff", 
     width = 12 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k1)
dev.off()

tiff(file = "./Figure/step2_ch2_kegg.tiff", 
     width = 12 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k2)
dev.off()

tiff(file = "./Figure/step2_ch3_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k3)
dev.off()

tiff(file = "./Figure/step2_ch4_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k4)
dev.off()

tiff(file = "./Figure/step2_ch5_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k5)
dev.off()

tiff(file = "./Figure/step2_ch6_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k6)
dev.off()

tiff(file = "./Figure/step2_ch7_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k7)
dev.off()

tiff(file = "./Figure/step2_ch8_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k8)
dev.off()

tiff(file = "./Figure/step2_ch10_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k10)
dev.off()

tiff(file = "./Figure/step2_ch11_kegg.tiff", 
     width = 12 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k11)
dev.off()

tiff(file = "./Figure/step2_ch12_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k12)
dev.off()

tiff(file = "./Figure/step2_ch13_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k13)
dev.off()

tiff(file = "./Figure/step2_ch14_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k14)
dev.off()

tiff(file = "./Figure/step2_ch15_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k15)
dev.off()

tiff(file = "./Figure/step2_ch16_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k16)
dev.off()

tiff(file = "./Figure/step2_ch17_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k17)
dev.off()

tiff(file = "./Figure/step2_ch18_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k18)
dev.off()

tiff(file = "./Figure/step2_ch19_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k19)
dev.off()

tiff(file = "./Figure/step2_ch20_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k20)
dev.off()

tiff(file = "./Figure/step2_ch22_kegg.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
print(k22)
dev.off()
