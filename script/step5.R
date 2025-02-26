library(Seurat)
library(tidyverse)
library(harmony)
library(DoubletFinder)

source("Rfunction_DoubletFinder.R")

folders=list.files('./',pattern='filtered')
names(folders) = c('GDM1', 'GDM2', 'GDM3',
                   'GDM4', 'NC1', 'NC2',
                   'NC3')

scRNAlist <- list()
for(i in 1:length(folders)){
  counts <- Read10X(data.dir = folders[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells = 3, min.features = 200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]]$log10GenesPerUMI <- log10(scRNAlist[[i]]$nFeature_RNA)/log10(scRNAlist[[i]]$nCount_RNA)
}

sce_sc <- merge(scRNAlist[[1]], 
                y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],
                      scRNAlist[[5]],scRNAlist[[6]],
                      scRNAlist[[7]]
                ))

pdf("QC-1.pdf", width = 15)
VlnPlot(sce_sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI"), ncol = 4, group.by="orig.ident")
dev.off()

### del
scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 200 &
                log10GenesPerUMI > 0.7 &
                percent.mt < 10 & 
                nCount_RNA < quantile(nCount_RNA,0.97) & 
                nCount_RNA > 1000)})

for (i in 1:length(scRNAlist)){
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- ScaleData(scRNAlist[[i]], verbose = FALSE)
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], verbose = FALSE)
  scRNAlist[[i]] <- RunPCA(scRNAlist[[i]], npcs = 30, verbose = FALSE)
  scRNAlist[[i]] <- RunUMAP(scRNAlist[[i]], reduction = "pca", dims = 1:30)
  scRNAlist[[i]] <- FindNeighbors(scRNAlist[[i]], reduction = "pca", dims = 1:30)
  scRNAlist[[i]] <- FindClusters(scRNAlist[[i]], resolution = 0.8)
  
  Doubletrate = ncol(scRNAlist[[i]])*8*1e-6
  sweep.res.list <- paramSweep(scRNAlist[[i]], PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  sweep.stats[order(sweep.stats$BCreal),]
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  homotypic.prop <- modelHomotypic(scRNAlist[[i]]$seurat_clusters)
  nExp_poi <- round(Doubletrate*nrow(scRNAlist[[i]]@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  scRNAlist[[i]] <- doubletFinder(scRNAlist[[i]], PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  scRNAlist[[i]]$doubFind_res = scRNAlist[[i]]@meta.data %>% select(contains('DF.classifications'))
  scRNAlist[[i]]$doubFind_score = scRNAlist[[i]]@meta.data %>% select(contains('pANN'))
  scRNAlist[[i]] = subset(scRNAlist[[i]], doubFind_res == "Singlet")
}

sce <- merge(scRNAlist[[1]], 
             y = c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],
                   scRNAlist[[5]],scRNAlist[[6]],
                   scRNAlist[[7]]
             ), 
             project = "GDM_qdx")

sce <- NormalizeData(sce)

test.sec <- CellCycleScoring(sce, s.features = Seurat::cc.genes.updated.2019$s.genes,
                             g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)
sce$Phase = test.sec$Phase
sce$S.Score = test.sec$S.Score
sce$G2M.Score = test.sec$G2M.Score

sce = FindVariableFeatures(sce, selection.method = "vst",nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score"))
sce <- RunPCA(sce, npcs = 50, verbose = FALSE)

pdf("ElbowPlot1.pdf")
ElbowPlot(sce, ndims=50, reduction="pca")
dev.off()

pc.num = 1:12
sce_harmony <- RunHarmony(sce, group.by.vars = "orig.ident")
sce_harmony <- RunUMAP(sce_harmony, reduction = "harmony", dims = pc.num)
sce_harmony <- RunTSNE(sce_harmony, reduction = "harmony", dims = pc.num)
sce_harmony <- FindNeighbors(sce_harmony, reduction = "harmony", dims = pc.num)
sce_harmony <- FindClusters(sce_harmony, reduction = "harmony", resolution = 0.8)

saveRDS(sce_harmony, "GDM_step1data.rds")

## annotation celltype
rm(list = ls())

setwd("F:/obesity-GDM")
sce_harmony = readRDS("GDM_step1data.rds")

library(Seurat)
library(tidyverse)
library(SCP)
library(BiocParallel)
library(SingleR)
library(celldex)
library(ggsc)
library(Ragas)

register(MulticoreParam(workers = 8, progressbar = TRUE))
sce_harmony = readRDS("GDM_step1data.rds")

## T cells "CD3E", "CD3D" 0,11,12,15,16,19,2,21,4,7
## plasma "MZB1", "JCHAIN", "IGKC" 24
## B  "MS4A1", 'CD79B', 'CD79A' 3,23,25
## neutrophils 'FCGR3B', 10,13
## NK 'FGFBP2' 'KLRF1' 8
## γδ T cells  'TOP2A', 'ASPM', 'CENPF' 22
## Mac 'CD68', 'CD14' 1,14,20,5
## mono 'FCN1' 6,9 
## DC 'IL3RA', 'FCER1A' 18,24
## mast cel 'KIT', 'CPA3' 17

## T cells "CD3E", "CD3D" 11,12,15,16,21
## CD4 T cells 0,2,4
## CD8a 7
## CD4-CD8 19


VlnPlot(sce_harmony,
        features = c("CD3E", "CD3D",
                     "MZB1", "JCHAIN", "IGKC",
                     "MS4A1", 'CD79B', 'CD79A',
                     'FCGR3B', 
                     'FGFBP2', 'KLRF1',
                     'TOP2A', 'ASPM', 'CENPF',
                     'CD68', 'CD14',
                     'FCN1',
                     'IL3RA', 'FCER1A',
                     'KIT', 'CPA3'
                     ),
        split.by = 'seurat_clusters',stack = TRUE, 
        flip = T) + 
  theme(legend.position='none')+labs(x='',y='')

VlnPlot(sce_harmony,
        features = c("MZB1", "JCHAIN", "IGKC",
                     "MS4A1", 'CD79B', 'CD79A',
                     'IL3RA', 'FCER1A'
                     
        ),
        split.by = 'seurat_clusters',stack = TRUE, 
        flip = T) + 
  theme(legend.position='none')+labs(x='',y='')

### sub T cell

## T cells "CD3E", "CD3D" 0,11,12,15,16,19,2,21,4,7

## plasma "MZB1", "JCHAIN", "IGKC" 24
## B  "MS4A1", 'CD79B', 'CD79A' 3,23,25
## neutrophils 'FCGR3B', 10,13
## γδ T cells  'TOP2A', 'ASPM', 'CENPF' 22
## Mac 'CD68', 'CD14' 1,14,20,5
## mono 'FCN1' 6,9 
## DC 'IL3RA', 'FCER1A' 18
## mast cell 'KIT', 'CPA3' 17
## NK T cells 'KLRD1', "GZMH" 'FGFBP2' 'KLRF1' 8
## CD4+ native T cells  'TCF7', 'IL7R', 'CCR7', 'LEF1'  0, 2, 4
## CD8+ native T cells  'TCF7', 'IL7R', 'CCR7', 'LEF1'  7
## native T cells  'TCF7', 'IL7R', 'CCR7', 'LEF1' 11,12,15,16,21
## CD4+CD8+ native T cells  'TCF7', 'IL7R', 'CCR7', 'LEF1' 19
{
sce_harmony$seurat_clusters <- as.character(sce_harmony$seurat_clusters)
sce_harmony$celltype <- sce_harmony$seurat_clusters

sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('24')] <- "Plasma cells"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('3','23','25')] <- "B cells"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('10','13')] <- "Neutrophils"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('8')] <- "NK T cells"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('22')] <- "gamma-delta T cells"  ## γδ T cells
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('1','14','20','5')] <- "Macrophage"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('6','9')] <- "Monocyte"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('18')] <- "Dendritic Cells"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('17')] <- "Mast Cells"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('0', "2", "4")] <- "Native CD4+ T cells"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('7')] <- "Native CD8+ T cells"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('11', "12", "15", "16", "21")] <- "Native T cells"
sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('19')] <- "Native CD4+CD8+ T cells"
}


"0"  "1"  "10" "11" "12"
"13" "14" "15" "16" "17"
"18" "19" "2"  "20" "21"
"22" "23" "24" "25" "3"
"4"  "5"  "6"  "7"  "8"  "9"

## Assigning cell type identity to clusters
cluster.ids <- c("Native CD4+ T cells", "Macrophage", "Neutrophils", "Native T cells", "Native T cells",
                 "Neutrophils", "Macrophage", "Native T cells", "Native T cells", "Mast Cells",
                 "Dendritic Cells", "Native CD4+CD8+ T cells", "Native CD4+ T cells", "Macrophage", "Native T cells",
                 "γδ T cells", "B cells", "Dendritic Cells", "B cells", "B cells",
                 "Native CD4+ T cells", "Macrophage", "Monocyte", "Native CD8+ T cells", "NK T cells", "Monocyte")
names(cluster.ids) <- levels(sce_harmony)
sce_harmony <- RenameIdents(sce_harmony, cluster.ids)
{
  labels <- enc2utf8("γδ T cells")
  sce_harmony$celltype[sce_harmony$seurat_clusters %in% c('22')] <- labels
}

zzm60colors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
                 '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
                 '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
                 '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
                 '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
                 '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
                 '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
                 '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
                 '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
                 '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
                 '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
                 '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')


sce_harmony$Histopathology = sce_harmony$orig.ident
sce_harmony$CellType = sce_harmony$celltype
sce_harmony$CellType = factor(sce_harmony$CellType, levels = c("Plasma cells", "B cells", "Neutrophils",
                                                               "NK T cells", "gamma-delta T cells",
                                                               "Native T cells", "Native CD4+ T cells",
                                                               "Native CD8+ T cells", "Native CD4+CD8+ T cells",
                                                               "Macrophage", "Monocyte", "Dendritic Cells",
                                                               "Mast Cells"))
sce_harmony$Cluster = sce_harmony$seurat_clusters
sce_harmony$Cluster = factor(sce_harmony$Cluster, levels = c("0", "1", "2", "3", "4",
                                                             "5", "6", "7", "8", "9",
                                                             "10", "11", "12", "13", "14",
                                                             "15", "16", "17", "18", "19",
                                                             "20", "21", "22", "23", "24",
                                                             "25"))


col <- zzm60colors

pdf(file="./step5_1.pdf", width=15)
CellDimPlot(                                                                           
  srt = sce_harmony, group.by = c("Histopathology", "Phase", "Cluster", "CellType"),
  reduction = "tsne", palcolor = col                    
)
dev.off()


features = c("MZB1", "JCHAIN", "IGKC", 
             "MS4A1", 'CD79B', 'CD79A', 
             'FCGR3B',
             'KLRD1', "GZMH", 'FGFBP2', 'KLRF1',
             'TOP2A', 'ASPM', 'CENPF',
             "CD4", "CD8A", "CD8B",
             'TCF7', 'IL7R', 'CCR7', 'LEF1',
             'CD68', 'CD14',
             'FCN1',
             'IL3RA', 'FCER1A',
             'KIT', 'CPA3')

tiff(file = "step5_1.tiff", 
     width = 18 * 300,
     height = 10 * 300, 
     res = 300, 
     compression = "lzw")
sc_dot(sce_harmony, features=features, group.by = "CellType")
dev.off()


sce_harmony$Group = ifelse(sce_harmony$orig.ident %in% c("NC1", "NC2", "NC3"), "NC", "GDM")
RunProportionPlot(sce_harmony, ## input can be a Seurat object
                  ident = "CellType",
                  group.by = "Group",
                  method = "pooled",
                  pooled.prop.by = "cluster", ## default
                  axis.text.angle = 45,
                  axis.text.size = 8,
                  return.value = 'ggplot'
)

RunProportionPlot(sce_harmony, ## if input is a Pi object, an updated Pi object will be the returned by default
                  ident = "Group",
                  group.by = "CellType",
                  method = "pooled",
                  pooled.prop.by = "cluster", ## default
                  axis.text.size = 8)

RunProportionPlot(sce_harmony,
                  ident = "CellType",
                  group.by = "Group",
                  method = "pooled",
                  pooled.prop.by = "group",
                  axis.text.size = 8)

sce_harmony$Names = factor(sce_harmony$orig.ident, levels = c("GDM1", "GDM2", "GDM3", "GDM4",
                                                              "NC1", "NC2", "NC3"))
RunProportionPlot(sce_harmony,
                  ident = "CellType",
                  group.by = "Group",
                  method = "unpooled",
                  unpool.by = "Names",
                  unpool.ncol = 3,
                  title.text.size = 6)

