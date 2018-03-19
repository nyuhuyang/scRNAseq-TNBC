########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################

#======1.1 Setup the Seurat objects =========================
# Load the BH dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

samples <- c("272-HF-L","272-HF-R","272-NC-B","272-NC-R")
projects <- rep("EC-BH-4709",4)
conditions <- c("HF-L","HF-R","NC-B","NC-R")
BH_raw <- list()
for(i in 1:length(samples)){
  BH_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                              samples[i],"/outs/filtered_gene_bc_matrices/mm10/"))
  colnames(BH_raw[[i]]) <- paste0(conditions[i],
                                          "_",colnames(BH_raw[[i]]))
}
BH_list <- lapply(BH_raw, CreateSeuratObject,min.cells = 3,
            min.genes = 200, project = projects)
BH_list <- lapply(BH_list, FilterCells,subset.names = "nGene",low.thresholds = 500)
BH_list <- lapply(BH_list, NormalizeData)
BH_list <- lapply(BH_list, FindVariableGenes, do.plot = FALSE)
BH_list <- lapply(BH_list, ScaleData)
for(i in 1:length(samples)) BH_list[[i]]@meta.data$conditions <- conditions[i]

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
genes.use <- lapply(BH_list, function(x) head(rownames(x@hvg.info), 1000))
genes.use <- unique(unlist(genes.use))
for(i in 1:length(conditions)){
  genes.use <- intersect(genes.use, rownames(BH_list[[i]]@scale.data))
}
length(genes.use) # 1/10 of total

#======1.2 Perform a canonical correlation analysis (CCA) =========================
# Run multi-set CCA
BH <- RunMultiCCA(BH_list,genes.use = genes.use,num.cc = 30)
save(BH, file = "./data/BH_alignment.Rda")
remove(BH_raw)
remove(BH_list)

# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = BH, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = BH, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
plot_grid(p1, p2)

# CC Selection
p3 <- MetageneBicorPlot(BH, grouping.var = "conditions", dims.eval = 1:30, 
                        display.progress = FALSE)

#======1.3 QC ==================================
# Run rare non-overlapping filtering
BH <- CalcVarExpRatio(object = BH, reduction.type = "pca",
                      grouping.var = "conditions", dims.use = 1:10)
BH <- SubsetData(BH, subset.name = "var.ratio.pca",accept.low = 0.5)

#======1.4 align seurat objects =========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned

# Alignment
BH <- AlignSubspace(object = BH, reduction.type = "cca", grouping.var = "conditions", 
                            dims.align = 1:10)
# t-SNE and Clustering
BH <- FindClusters(object = BH, reduction.type = "cca.aligned", dims.use = 1:10, 
                   resolution = 1.2, force.recalc = T, save.SNN = TRUE)
BH <- RunTSNE(object = BH, reduction.use = "cca.aligned", dims.use = 1:10, do.fast = TRUE)

p1 <- TSNEPlot(BH, do.return = T, pt.size = 1, group.by = "conditions")
p2 <- TSNEPlot(BH, do.label = F, do.return = T, pt.size = 1)

plot_grid(p1, p2)


#Now, we annotate the clusters as before based on canonical markers.

TSNEPlot(object = BH,do.label = TRUE, group.by = "ident", 
         do.return = TRUE, no.legend = TRUE,
         pt.size = 1,label.size = 8 )+
  ggtitle("TSNE plot of all clusters")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle

# Compare clusters for each dataset
cell.all <- FetchData(BH,"conditions")
cell.subsets <- lapply(conditions, function(x) 
        rownames(cell.all)[cell.all$conditions == x])

BH.subsets <- list()
for(i in 1:length(conditions)){
        BH.subsets[[i]] <- SubsetData(BH, cells.use =cell.subsets[[i]])
}

table(BH.subsets[[1]]@ident)

p <- list()
for(i in 1:length(conditions)){
        p[[i]] <- TSNEPlot(object = BH.subsets[[i]],do.label = TRUE, group.by = "ident", 
                           do.return = TRUE, no.legend = TRUE,
                           pt.size = 1,label.size = 4 )+
                ggtitle(samples[i])+
                theme(text = element_text(size=20),     #larger text including legend title							
                      plot.title = element_text(hjust = 0.5)) #title in middle
}
do.call(plot_grid, p)
save(BH, file = "./data/BH_alignment.Rda")
remove(BH.subsets)
