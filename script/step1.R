## author:wjw
## description:annition_DMP_gene_data. 
## time:20250125

rm(list = ls())

library(tidyverse)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("org.Hs.eg.db")
library(xlsx)

load("./Rdata/mydata_DMP_DMR.Rdata")
myDMP = myDMP$Control_to_GDM
myDMR = myDMR$BumphunterDMR

##G-N up9327 down10908 850K
myDMP_up = myDMP %>% 
  as_tibble() %>% 
  filter(logFC > 0) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(gene != '') %>% 
  dplyr::select(gene) %>% 
  distinct()

myDMP_down = myDMP %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  filter(logFC < 0) %>% 
  mutate(gene = as.character(gene)) %>% 
  dplyr::filter(gene != '') %>% 
  dplyr::select(gene) %>% 
  distinct()

myDMR$group = ifelse(myDMR$value > 0, "Hypermethylated", "Hypomethylated")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

write.table(myDMR,file = "myDMR.txt",row.names = F,sep = "\t",quote = F)

myDMR <- annotatePeak("./myDMR.txt", tssRegion = c(-3000, 3000), TxDb = txdb,annoDb = 'org.Hs.eg.db')

myDMR = as.data.frame(myDMR@anno@elementMetadata@listData)
myDMR$SYMBOL

library(ggrepel)
library(ggsci)
library(tidyverse)

load("Rdata/Manhattan_data.Rdata")

cg = c("cg11723077", "cg17588003", "cg11187204", "cg10139436", "cg20433858", "cg00812770", "cg22791932", "cg10986462", "cg15626350", "cg11540997", "cg26385222", "cg25551168", "cg26267561", "cg21842274", "cg26799474", "cg18548103", "cg04036196", "cg10507965", "cg20198384", "cg00922748", "cg05216211", "cg05376185", "cg06617468", "cg17097119", "cg22385669", "cg11388673", "cg08799779", "cg22790973", "cg19169154", "cg23376861", "cg08077807", "cg22865713", "cg16995742", "cg02823329", "cg17283620", "cg14060113", "cg09101062", "cg15626350", "cg11540997", "cg26385222", "cg25551168", "cg26267561", "cg21842274", "cg26799474", "cg07052015", "cg12890605", "cg24396358", "cg26643856", "cg00575744", "cg16236263", "cg27478707", "cg01889448", "cg09555323", "cg15820961", "cg07223713", "cg15074838", "cg14428733", "cg23508786", "cg24865779", "cg06969118", "cg11187204", "cg20433858", "cg00812770", "cg10986462", "cg22790973", "cg23643951", "cg06319822", "cg26952618", "cg04988367", "cg08407434", "cg01617280", "cg13911801", "cg09304293", "cg23109897", "cg06928695", "cg18255813", "cg10431340", "cg19565171", "cg09109411")

genes = c("Peg3", "H19", "LTF", "DUSP22", "OR2L13", "CYP2E1", "CCDC110", "KALRN", "PAG1", "GNRH1", "SLC2A9", "CSRP2BP", "HIVEP1", "RALGDS", "DHX37", "SCNN1D", "HOOK2", "LCE3C", "TMEM63B", "AVP", "OXT", "CRHBP", "ESR1", "DUOX2", "TMEM176A", "CASP8", "SCD", "ITM2C", "NT5C3A", "NPEPL1", "PIGCP1", "ADAM3A", "ZSCAN12P1", "IGF2", "H19", "MEST", "NR3CI", "S100A8", "MS4A3", "MMP8", "PPAGRC1A", "HIF3A", "CYP27B1", "CYP24A2", "VDR", "NRIP2", "MIR1185-1", "MOBKL2C", "NAP1L5", "DEFB133", "AQR", "ZNF696", "PLAC8", "TFCP2", "GNAS", "ATP5A1", "MFAP4", "PRKCH", "SLC17A4", "HIF3A", "S100A8", "MS4A3", "MMP8", "COPS8", "PIK3R5", "HAAO", "CCDC124", "C5orf34", "H3K27", "H3K4", "TFCP2", "IL-10", "TCF7L2", "ADIPOQ", "MTNR1B", "GCK", "GCKR", "FTO", "IRS1", "KCNQ1", "SLC30A8", "CDKAL1", "CAPN10", "KCNJ11", "RBP4", "GC", "STK11", "MIF", "CDKN2A/2B", "IGF2BP2", "AVP", "OXT", "CRHBP", "ESR1", "DUOX2", "TMEM176A", "CASP8", "IGF2", "H19", "MEST", "NR3CI", "HLA-DOA", "ABCA1", "DLGAP2", "KCNQ1", "ADCYAPI", "PRDM16", "BMP7", "PPARGC1α", "LDLR", "BLM", "PDE4B", "TNFRSF1B", "CYP27B1", "CYP24A2", "VDR", "FTO", "EIF3F", "ABLIM1", "GRHL1", "HLA-F", "NDRG1", "SASH1", "Oas1", "Ppie", "Polr2g", "BAMBI", "INSR", "IRSI", "SOD2", "MAP4K3", "DUSP5", "PAK2", "SERPINE1", "PPP1R15B", "IGF1R", "ATG7", "DICER1", "RANBP2", "LPL", "PGC-1a", "PDX1", "STAT1", "HLA-DOA", "HLA-DMB", "HLA-DQB1", "HLA-DRB1", "HLA-DMA", "HLA-DRA","HLA-DPB1", "IFNGR2", "EIF2AK2", "GLUT3", "RETN", "RBP4", "GLUT3", "PPARα", "CTNND2", "HNF4A", "MFHAS1", "LOC92973", "PACRG", "PRKDC", "HBM", "PITPNM3", "RREB1", "FAM18A", "KHDC1L", "MPZ", "CAMTA1", "TXNIP", "PDE6A", "IR", "TXNIP", "PI3KR1", "ADIPOQ", "RETN", "C10orf10", "FSTL1", "GSTT1", "HLA-DPB1", "HLA-DRB5", "HSPA6", "HSPA6", "C10orf10", "FSTL1", "GSTT1", "HLA-DPB1", "HLA-DRB5", "HIF3A", "ESR1", "SLC2A4", "ADIPOQ", "LEP", "ABCA1","SLC2A1","GLUT1", "SLC2A3", "GLUT3", "TXNIP", "Rarres2")

qdx_genes = c("LTF", "DUSP22", "OR2L13", "CYP2E1", "CCDC110", "KALRN", "PAG1", "GNRH1", "SLC2A9", "CSRP2BP", "HIVEP1", "RALGDS", "DHX37", "SCNN1D", "HOOK2", "LCE3C", "TMEM63B", "AVP", "OXT", "CRHBP", "ESR1", "DUOX2", "TMEM176A", "CASP8", "SCD", "ITM2C", "NT5C3A", "NPEPL1", "PIGCP1", "ADAM3A", "ZSCAN12P1", "IGF2", "H19", "MEST", "NR3CI", "S100A8", "MS4A3", "MMP8", "PPAGRC1A", "HIF3A", "CYP27B1", "CYP24A2", "VDR", "NRIP2", "MIR1185-1", "MOBKL2C", "NAP1L5", "DEFB133", "AQR", "ZNF696", "PLAC8", "TFCP2", "GNAS", "ATP5A1", "MFAP4", "PRKCH", "SLC17A4", "HIF3A","PRTN3","HIST1H1B", "HIST1H4A")

print_gene = read.csv("./data/PRINT_GENE.csv")


dmp_1 %>% 
  as_tibble() %>% 
  dplyr::filter(gene %in% print_gene$Gene) %>% 
  dplyr::select(gene) %>% 
  distinct()


library("CMplot") 
gwasResults <- don %>% dplyr::select(probeID,CHR,bp,adj.P.Val) %>% 
  dplyr::rename(
    SNP = probeID,
  )


SNPs = don %>% dplyr::select(probeID,CHR,bp,adj.P.Val,gene) %>% 
  dplyr::rename(
    SNP = probeID,
  ) %>% filter(SNP %in% cg) %>% 
  filter(adj.P.Val < 0.05) %>% 
  distinct(gene, .keep_all = T)

print_gene_DMP = don %>% dplyr::select(probeID,CHR,bp,adj.P.Val,gene) %>% 
  dplyr::rename(
    SNP = probeID,
  ) %>% filter(gene %in% myDMR$SYMBOL) %>% 
  filter(gene %in% print_gene$Gene) %>%
  filter(adj.P.Val < 0.05) %>% 
  distinct(gene, .keep_all = T)

all_gene_DMP = don %>% dplyr::select(probeID,CHR,bp,adj.P.Val,gene) %>% 
  dplyr::rename(
    SNP = probeID,
  ) %>% filter(gene %in% myDMR$SYMBOL) %>% 
  filter(gene %in% genes) %>%
  filter(adj.P.Val < 0.05) %>% 
  distinct(gene, .keep_all = T)

all_gene_DMP = rbind(print_gene_DMP, all_gene_DMP) %>% 
  distinct(gene, .keep_all = T)

tiff(file = "./Figure/step1_1.tiff", 
     width = 20 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
CMplot(gwasResults,
       col = c("#FFFFE0", "#FFFACD","#FAFAD2",
               "#FFEFD5", "#FFE4B5", "#EEE8AA","#F0E68C"),
       #col = c("#DCDCDC", "grey"),
       highlight = all_gene_DMP$SNP,
       highlight.col = NULL,
       highlight.cex = 1,
       highlight.text = all_gene_DMP$gene,      
       highlight.text.col = "black",
       type="p",
       plot.type="m",
       LOG10=TRUE,
       threshold=NULL,
       file="tiff",
       dpi=300,
       file.output=F,
       verbose=F,
       width=20,height=8,
       chr.labels.angle=0,
       axis.cex = 1,
       lab.cex = 1.5,
       lab.font= 4,
       ylab.pos=2,
       ylab = expression(-log[2](italic(adj.P.Val))))
dev.off()

##DMP_distribution_plot

dmp_1$group = ifelse(dmp_1$logFC > 0,"Hypermethylated","Hypomethylated")
dmp_1$group = factor(dmp_1$group,levels = c("Hypermethylated","Hypomethylated"))
dmp_1$CHR = factor(dmp_1$CHR)
dmp_1$feature = factor(dmp_1$feature,levels = c("ExonBnd","3'UTR","1stExon","TSS200",
                                                "5'UTR","TSS1500","IGR","Body")) 

len_data <- dmp_1 %>% 
  group_by(feature,group) %>% 
  summarise(n()) %>% 
  dplyr::rename(len = `n()`)

len_data


tiff(file = "./Figure/step1_2.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
len_data %>% 
  ggplot(aes(x=feature,y = len,fill = group)) +
  geom_bar(stat="identity", position=position_dodge(0.9),colour = "black") +
  theme_classic() + 
  scale_y_continuous(name = NULL ,expand = c(0,0),limits = c(0,10500)) +
  scale_x_discrete(name = NULL) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.key.size = unit(25, "pt"),
        legend.position = "none"
        #  panel.grid.major.y = element_blank()
  ) +
  guides(fill = guide_legend(title = "DMPs Group")) +
  #  scale_fill_brewer(palette="Accent") +
  geom_text(aes(label=len), vjust=-0.2, color="black",
            position = position_dodge(0.9), size=4,) +
  scale_fill_cosmic()
dev.off()

data = dmp_1 %>% as_tibble() %>% 
  group_by(cgi, group) %>% 
  summarise(n()) %>% 
  dplyr::rename(value = `n()`)

tiff(file = "./Figure/step1_3.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
data %>% 
  ggplot(aes(x=cgi,y = value,fill = group)) +
  geom_bar(stat="identity", position=position_dodge(0.9),colour = "black") +
  theme_classic() + 
  scale_y_continuous(name = NULL, expand = c(0,0),
                     limits = c(0,17000)
  ) +
  scale_x_discrete(name = NULL) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.key.size = unit(25, "pt"),
        legend.position = "right"
        #  panel.grid.major.y = element_blank()
  ) +
  guides(fill = guide_legend(title = "DMPs Group")) +
  #  scale_fill_brewer(palette="Accent") +
  geom_text(aes(label=value), vjust=-0.2, color="black",
            position = position_dodge(0.9), size=4,) +
  scale_fill_cosmic()
dev.off()


chr_data = dmp_1 %>% as_tibble() %>% 
  group_by(CHR, group) %>% 
  summarise(n()) %>% 
  dplyr::rename(Chr = `n()`)


tiff(file = "./Figure/step1_4.tiff", 
     width = 20 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
chr_data %>% 
  ggplot(aes(x=CHR,y = Chr,fill = group)) +
  geom_bar(stat="identity", position=position_dodge(0.9),colour = "black") +
  theme_classic() + 
  scale_y_continuous(name = "Counts of DMPs",expand = c(0,0),
                     limits = c(0,3100)
  ) +
  scale_x_discrete(name = NULL) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.key.size = unit(25, "pt"),
        legend.position = "none"
        #  panel.grid.major.y = element_blank()
  ) +
  guides(fill = guide_legend(title = "DMPs Group")) +
  #  scale_fill_brewer(palette="Accent") +
  geom_text(aes(label=Chr), vjust=-0.2, color="black",
            position = position_dodge(0.9), size=4,) +
  scale_fill_cosmic()
dev.off()


##DMR_distribution_plot
rm(list = ls())
load("Rdata/mydata_DMP_DMR.Rdata")
rm(myDMP)

DMR_bed = myDMR$BumphunterDMR %>% as_tibble() %>% dplyr::select(seqnames,start,end)
write.table(DMR_bed,file = "myDMR_bed.txt",row.names = F,sep = "\t",quote = F)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

myDMR_anni <- annotatePeak("./myDMR_bed.txt", tssRegion = c(-3000, 3000), TxDb = txdb,annoDb = 'org.Hs.eg.db')
myDMR$BumphunterDMR$anno = myDMR_anni@anno$annotation
myDMR$BumphunterDMR$group = ifelse(myDMR$BumphunterDMR$value < -0.2, "Hypomethylated", 
                                   ifelse(myDMR$BumphunterDMR$value > 0.2, "Hypermethylated", 'ns'))

table(myDMR_anni@anno$annotation)

DMR = myDMR$BumphunterDMR %>% as_tibble() %>% 
  arrange(seqnames)
DMR %>% group_by(group) %>% 
  summarise(n())
DMR_UP = DMR %>% dplyr::filter(group == "Hypermethylated") %>% 
  group_by(seqnames) %>% 
  summarise(n()) %>% 
  mutate(group = "Hypermethylated")
colnames(DMR_UP) = c("seqnames","number","group")

supp_dmr_up = data.frame(seqnames = c("chr13","chr16","chr20","chr21"),
                           number = NA,
                           group = "Hypermethylated") 
DMR_UP = rbind(DMR_UP,supp_dmr_up)

DMR_down = DMR %>% dplyr::filter(group == "Hypomethylated") %>% 
  group_by(seqnames) %>% 
  summarise(n()) %>% 
  mutate(group = "Hypomethylated")
colnames(DMR_down) = c("seqnames","number","group")

DMR = rbind(DMR_UP,DMR_down)

DMR$group = factor(DMR$group,levels = c("Hypermethylated","Hypomethylated"))
DMR$seqnames = substr(DMR$seqnames, start = 4, stop = 6)
DMR$seqnames = factor(DMR$seqnames,levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))



tiff(file = "./Figure/step1_5.tiff", 
     width = 20 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
DMR %>% 
  ggplot(aes(x=seqnames,y=number,fill = group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  scale_y_continuous(name = "Counts of DMRs",expand = c(0,0), limits = c(0,25)) +
  scale_x_discrete(name = NULL) +
  theme(axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.key.size = unit(25, "pt"),
        legend.position = "right") +
  guides(fill = guide_legend(title = "DMRs Group")) +
  geom_text(aes(label=number), vjust=-0.2, color="black",
            position = position_dodge(0.9), size=4,)+
  scale_fill_cosmic()
dev.off()


####
DMR = myDMR$BumphunterDMR %>% as_tibble() %>% 
  arrange(seqnames)
DMR = DMR %>% 
  mutate(sub_group = ifelse(anno %in% c("Promoter (1-2kb)", "Promoter (2-3kb)", "Promoter (<=1kb)"), "Promoter", "Others"))

DMR_data = DMR %>% group_by(sub_group) %>% 
  summarise(n()) %>% 
  dplyr::rename(value = `n()`)
DMR_data$sub_group = factor(DMR_data$sub_group, levels = c("Promoter", "Others"))

tiff(file = "./Figure/step1_6.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
ggplot(DMR_data, aes(x = "", y = value, fill = sub_group)) +
  geom_col(color = "grey") +
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5),
            size = 12, color = "white") +
  coord_polar(theta = "y") +
  #scale_fill_brewer(palette="Accent") +
  theme_void() +
  guides(fill = guide_legend(title = "DMRs Group")) +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.key.size = unit(25, "pt"),
        legend.position = "right") +
  scale_fill_cosmic()
dev.off()

rm(list = ls())

#### volcano
library(tidyverse)
library(ggrepel)
library(ggfun)
library(grid)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("org.Hs.eg.db")

####----load data----####
load("./Rdata/allDMR.Rdata")
allmyDMR = myDMR$BumphunterDMR
allmyDMR$value = allmyDMR$value * -1
allmyDMR$group = ifelse(allmyDMR$value > 0,"Hypermethylated","Hypomethylated")

allmyDMR %>% 
  as_tibble() %>% 
  filter(p.value < 0.05) %>% 
  filter(value < -0.2)

load("Rdata/mydata_DMP_DMR.Rdata")
rm(myDMP)

myDMR = myDMR$BumphunterDMR

subDMR = allmyDMR %>% 
  rownames_to_column() %>%
  as_tibble() %>%
  filter(start %in% myDMR$start) %>% 
  filter(end %in% myDMR$end)

allmyDMR = allmyDMR %>% 
  rownames_to_column() %>%
  as_tibble()

allmyDMR$change = ifelse(allmyDMR$rowname %in% subDMR$rowname,"Change","Normal")

DMR_bed = allmyDMR %>% dplyr::select(seqnames,start,end)
write.table(DMR_bed,file = "all_myDMR_bed.txt",row.names = F,sep = "\t",quote = F)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

myDMR_anni <- annotatePeak("./all_myDMR_bed.txt", tssRegion = c(-3000, 3000), TxDb = txdb,annoDb = 'org.Hs.eg.db')
allmyDMR$anno = myDMR_anni@anno$annotation
allmyDMR$gene = myDMR_anni@anno$SYMBOL

allmyDMR = allmyDMR %>% 
  as_tibble() %>% 
  mutate(Foldchange = value)

filt_gene = unique(c(genes, print_gene$Gene)) 

tiff(file = "./Figure/step1_7.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
ggplot(data = allmyDMR) + 
  geom_point(aes(x = Foldchange, y = -log2(p.value), 
                 color = Foldchange,
                 size = -log2(p.value))) + 
  geom_point(data =  allmyDMR %>%
               dplyr::filter(change != "Normal") %>% 
               dplyr::filter(gene %in% filt_gene),
             aes(x = Foldchange, y = -log2(p.value),
                 # fill = log2FoldChange,
                 size = -log2(p.value)),
             shape = 21, show.legend = F, color = "#000000") +
  geom_text_repel(data =  allmyDMR %>% 
                    dplyr::filter(change != "Normal") %>%
                    dplyr::filter(gene %in% filt_gene) %>% 
                    dplyr::filter(group == "Hypermethylated"),
                  aes(x = Foldchange, y = -log2(p.value), label = gene),
                  box.padding = 0.5,
                  nudge_x = 0.5,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  # segment.angle = 10,
                  direction = "y", 
                  hjust = "left"
  ) + 
  geom_text_repel(data =  allmyDMR %>% 
                    dplyr::filter(change != "Normal") %>%
                    dplyr::filter(gene %in% filt_gene) %>% 
                    dplyr::filter(group == "Hypomethylated"),
                  aes(x = Foldchange, y = -log2(p.value), label = gene),
                  box.padding = 0.5,
                  nudge_x = -0.2,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  direction = "y", 
                  hjust = "right"
  ) + 
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  scale_size(range = c(1,7)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.text = element_text(size = 15, color = "#000000"),
        axis.title = element_text(size = 20),
  ) + 
  scale_y_continuous(limits = c(-0.1,25)) +
  scale_x_continuous(limits = c(-1.5,1.5)) + 
  annotate(geom = "text", x = 1, y = 0, label = "p.value = 0.05", size = 5) + 
  coord_cartesian(clip = "off") + 
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#74add1")
    ), 
    xmin = -1.4, 
    xmax = -0.4,
    ymin = 22,
    ymax = 22
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Down",
      gp = grid::gpar(col = "#74add1")
    ),
    xmin = -0.7, 
    xmax = -1.2,
    ymin = 23,
    ymax = 23
  ) +
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
      gp = grid::gpar(lwd = 3, col = "#d73027")
    ), 
    xmin = 0.4, 
    xmax = 1.4,
    ymin = 22,
    ymax = 22
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Up",
      gp = grid::gpar(col = "#d73027")
    ),
    xmin = 0.7, 
    xmax = 1.2,
    ymin = 23,
    ymax = 23
  ) 

dev.off()

