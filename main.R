library(DoubletFinder)
data_dir <- 'C://Bioinf/HUMAN_FETAL_RETINA/GW/'
list.files(data_dir)
D59fetal <- Read10X(data.dir = data_dir)
D59fetalS = CreateSeuratObject(counts = D59fetal)

meta <- read.csv('C://Bioinf/HUMAN_FETAL_RETINA/GSE138002_All_barcodes.csv')
rownames(meta) <- meta$X
D59fetalS <- AddMetaData(D59fetalS, metadata = meta)

table(D59fetalS$Age)

HGW <- subset(D59fetalS, subset = Age == 'Hgw27')

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat) #could be replaced with SCTransform
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:50)
  Seurat <- FindClusters(Seurat, resolution = 1)
  Seurat <- RunUMAP(Seurat, dims = 1:50)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:50 )
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

RDoublet <- function(tmp){
  sweep.res.list <- paramSweep_v3(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters) 
  nExp_poi <- round(0.1*length(colnames(tmp)))  ## Assuming 10% doublet formation rate 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp) 
}


HGW[["percent.rb"]] <- PercentageFeatureSet(HGW, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
HGW[["percent.mt"]] <- PercentageFeatureSet(HGW, pattern = "^MT-")
HGW <- CellCycleScoring(HGW, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE, nbin = 6)

VlnPlot(HGW, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)

HGW <- subset(HGW, subset = nCount_RNA > 300 & nCount_RNA < 15000 & nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 15 & percent.rb < 40)
HGW <- ScaleData(HGW, verbose = T, vars.to.regress = c('nCount_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score")) 

HGW <- ProcessSeu(HGW)

#Basic visualisation
FeaturePlot(HGW, features = c('RBPMS'))
DimPlot(HGW, reduction = 'umap', label = TRUE, repel = TRUE)

HGW <- RDoublet(HGW)
HGW <- subset(HGW, cells = colnames(HGW )[which(HGW [[]][17] == 'Singlet')])
HGW <- subset(HGW , cells = colnames(HGW )[which(HGW [[]][18] == 'Singlet')])

HGW <- ProcessSeu(HGW)

HGW$origin <- 'Lu'
HGW$timepoint <- 'Week 27'

SaveH5Seurat(HGW, 'C://Bioinf/HUMAN_FETAL_RETINA/HGW27.h5Seurat')
saveRDS(HGW, 'C://Bioinf/HUMAN_FETAL_RETINA/HGW27.rds')
