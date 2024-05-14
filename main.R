

RDoublet <- function(tmp){
  sweep.res.list <- paramSweep_v3(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters)
  nExp_poi <- round(0.05*length(colnames(tmp)))  ## Assuming 10% doublet formation rate
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp)
}

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 3000)
  Seurat <- ScaleData(Seurat) #could be replaced with SCTransform
  Seurat <- RunPCA(Seurat, npcs = 100)
  Seurat <- FindNeighbors(Seurat, dims = 1:100)
  Seurat <- FindClusters(Seurat, resolution = 1)
  Seurat <- RunUMAP(Seurat, dims = 1:100)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:100)
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

data_dir <- 'C://Bioinf/HUMAN_FETAL_RETINA/D59_fetal_filtered_gene_bc_matrices/GRCh38'
list.files(data_dir)
D59fetal <- Read10X(data.dir = data_dir)
D59fetalS = CreateSeuratObject(counts = D59fetal)

D59fetalS[["percent.rb"]] <- PercentageFeatureSet(D59fetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
D59fetalS[["percent.mt"]] <- PercentageFeatureSet(D59fetalS, pattern = "^MT-")
D59fetalS <- CellCycleScoring(D59fetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
VlnPlot(D59fetalS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
D59fetalS <- subset(D59fetalS, subset = nCount_RNA > 300 & nCount_RNA < 15000 & nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 15 & percent.rb < 40)
D59fetalS <- ScaleData(D59fetalS, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score"))

D59fetalS <- ProcessSeu(D59fetalS)
D59fetalS <- RDoublet(D59fetalS)
D59fetalS <- subset(D59fetalS, cells = colnames(D59fetalS )[which(D59fetalS [[]][12] == 'Singlet')])
D59fetalS <- subset(D59fetalS , cells = colnames(D59fetalS )[which(D59fetalS [[]][13] == 'Singlet')])
D59fetalS <- ProcessSeu(D59fetalS)
D59fetalS$origin <- 'Thomas Reh Lab'
D59fetalS$timepoint <- 'Day 59'
D59fetalS$tissue <-  'Total retina'
D59fetalS$donor <- '1'
SaveH5Seurat(D59fetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/v2/FD59.h5Seurat', overwrite = TRUE)
saveRDS(D59fetalS, 'C://Bioinf/HUMAN_FETAL_RETINA/v2/FD59.rds')

gc()
