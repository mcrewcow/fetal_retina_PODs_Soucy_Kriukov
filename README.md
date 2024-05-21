# IGFBPL1 is a Master Driver of Microglia Homeostasis and Resolution of Inflammation that Limits Glaucomatous Neurodegeneration
<a href="https://zenodo.org/badge/latestdoi/600119531"><img src="https://zenodo.org/badge/600119531.svg" alt="DOI"></a> <br />
This repository has the information on the single-cell RNA-seq analysis described in the paper for IGFBPL1 WT/KO and glaucoma PBS/IGFBPL1 mouse microglia datasets.
![GA1_Github](https://user-images.githubusercontent.com/77118598/219808488-bb277362-bc15-4820-84cb-9eb30d7e5972.jpg)


The analysis is separated into the following parts:
1. Igfbpl1 WT/KO - glaucoma PBS/Igfbpl1 Seurat R processing - basic Seurat pipeline, availabe upon request (gahu@mit.edu)
2. ForceAtlas2 + scanpy processing (Jupyter Notebook)
3. RNA Velocity - bash preprocessing to .loom - file: [fastq_to_loom_processing.sh](https://github.com/mcrewcow/Pan_et_al_DFChen_IGFPL1_paper/blob/main/fastq_to_loom_processing.sh) <br />
                - scvelo utilization (Jupyter Notebook)
  Files: [Glaucoma datasets analysis](https://github.com/mcrewcow/Pan_et_al_DFChen_IGFPL1_paper/blob/main/glaucoma.ipynb), [IGFPL1 WT/KO datasets analysis](https://github.com/mcrewcow/Pan_et_al_DFChen_IGFPL1_paper/blob/main/IGFBPL1.ipynb) <br />
4. Python to R data transformation
Files: [Saving the object in Python](https://github.com/mcrewcow/Pan_et_al_DFChen_IGFPL1_paper/blob/main/h5ad_to_h5seurat_p1.py), [Building the Seurat object in R](https://github.com/mcrewcow/Pan_et_al_DFChen_IGFPL1_paper/blob/main/h5ad_to_h5seurat_p2.R)
                
The datasets are uploaded on GEO: <br />
[WT/IGFBPL1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176204) <br />


The processed datasets are available upon request. <br />

If you have any questions regarding the data, analysis or further improvements, here are the contacts: <br />
dongfeng_chen@meei.harvard.edu - PI\
gahu@mit.edu - Bioinformatics, part 1 + data availability\
ekriukov@meei.harvard.edu - Bioinformatics, part 2-4 + data availability<br />
<br />
