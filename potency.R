library(CytoTRACE2)
fetal <- LoadH5Seurat('C://Bioinf/HUMAN_FETAL_RETINA/COMBINED_EKPB_v1_clean.h5Seurat')


cytotrace2_result <- cytotrace2(fetal,   
                                species = "human",
                                is_seurat = TRUE,
                                slot_type = "counts",
                                full_model = FALSE,
                                batch_size = 10000,
                                smooth_batch_size = 1000,
                                parallelize_models = TRUE,
                                parallelize_smoothing = TRUE,
                                ncores = NULL,
                                max_pcs = 200,
                                seed = 14)  


annotation <- data.frame(phenotype = fetal@meta.data$EK_PB_annov1) %>% set_rownames(., colnames(fetal))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = TRUE)
gc()

emb_1 <- fetal@reductions$umap@cell.embeddings[,1]  # First dimension of your "harmony" reduction embedding
emb_2 <- fetal@reductions$umap@cell.embeddings[,2]  # Second dimension of your "harmony" reduction embedding

# replace CytoTRACE2_UMAP embeddings
plots[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2

plots$CytoTRACE2_UMAP
