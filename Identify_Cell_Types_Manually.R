library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/BH_alignment.Rda")
lnames

# all marker genes 
marker.genes <- c("SLC36A2","P2RX5","TRIP4","ASCC1","MYF5","UCP1","Cdh5","Pecam1",
                  "Flt1","Vwf","Plvap","Kdr","EMCN","Car4","ptprb","KRT19","Epcam",
                  "KRT5","MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1","Rpe65",
                  "Rlbp1","Msln","Upk3b","Lrrn4","Ccl25","Cxcl12","Ctsl","Psmb11","Aire",
                  "HLA-DMA","Krt5","Gas1","Plet1","Ly6d","Spink5","Reg3g","Bpifa1","FGF1",
                  "FGF9","SFRP1","PTPRC","LAPTM5","SRGN","PPBP","GNG11","HBA2","HBB",
                  "Cma1","Mcpt4","Tpsb2","Cpa3","LYZ","S100A9","CD14","CCL2","Emr1","CD68",
                  "MARCO","LYZ","FCGR3A","MS4A7","VMO1","Itgax","GPR183","CST3","HLA-DQA1",
                  "FCER1A","TSPAN13","IL3RA","IGJ","CD3G","CD3D","CD2","CREM","Cd4","CD62L",
                  "IL7R","IL2RG","GIMAP5","SELL","CREM","Cd8a","CCL5","CACYBP","Foxp3","Cd19",
                  "CD79A","MS4A1","CREM","MIR155HG","NME1","HLA-DQA1","CD27","SDC1","IL6R",
                  "SLAMF7","GNLY","Ncr1","CCL5","KLRD1","NKG7","Pmel","Mlana","Pdgfrb","Vim",
                  "Has2","Dcn","Eng","CD44","Nt5e","Alcam","Kit","Ly6a","Thy1","Itgb1","Vcam1",
                  "Icam2","CD72","CD2","MBP","MPZ","Rbfox3","Pdgfrb","Cspg4","Anpep","Rgs5",
                  "Myh11","Mylk","Sost","Des","Vtn","Ifitm1","Acta2","Myh11","POU5F1","FUT4",
                  "CD34","PROM1","ABCG2","ATXN1","NES","NCAM","NGFR","GPX3","DCN","COL6A1",
                  "TIMP3","PDGFRA")
marker.genes <- MouseGenes(BH, marker.genes)
for(i in 1:(length(marker.genes))) print(paste0(i," ",marker.genes[i]))
# Adipocytes
Adipocytes <- marker.genes[1:5]
FeaturePlot(object = BH, 
            features.plot = Adipocytes, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

# Endothelial Cells
Endothelial <- marker.genes[6:14]
FeaturePlot(object = BH, 
            features.plot = Endothelial, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

# Epithelium
Epithelium <- marker.genes[c(15:16,19:28)] 
FeaturePlot(object = BH, 
            features.plot = Epithelium, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

#fibroblasts
fibroblasts <- marker.genes[29:31] 
FeaturePlot(object = BH, 
            features.plot = fibroblasts, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
# Hematopoietic cells
Hematopoietic <- marker.genes[32:34]
myeloid <- marker.genes[c(36:47)]
Tcell <- marker.genes[c(48:59)]
Bcell <- marker.genes[c(62:70)]
NK <- marker.genes[c(71:73)]
FeaturePlot(object = BH, 
            features.plot = NK, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
#Pericytes and Mesenchymal Cells
Mesenchymal <- marker.genes[c(75:90)] 
FeaturePlot(object = BH, 
            features.plot = Mesenchymal, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
# Myelinating Schwann cells ,Pericytes
Myelinating <- marker.genes[c(91:93)] 
Pericytes <- marker.genes[c(94:102)]
muscle <- marker.genes[c(103:104)]
FeaturePlot(object = BH, 
            features.plot = muscle, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
# Stem cell
stem <- marker.genes[105:112]
FeaturePlot(object = BH, 
            features.plot = stem, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)

# Rename ident
table(BH@ident)
idents <- as.data.frame(table(BH@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("0.Unknown",
                     "1.Unknown",
                     "2.Unknown",
                     "3.Unknown",
                     "4.Myeloid cells",
                     "5.Unknown",
                     "6.Unknown",
                     "7.Unknown",
                     "8.Unknown",
                     "9.Epithelium",
                     "10.Mesenchymal Cells",
                     "11.Unknown",
                     "12.Myeloid cells",
                     "13.\Lymphoid cells",
                     "14.Endothelium",
                     "15.Epithelium",
                     "16.Endothelium")

BH@ident <- plyr::mapvalues(x = BH@ident,
                            from = old.ident.ids,
                            to = new.cluster.ids)

# The SplitDotPlotGG function can be useful for viewing conserved cell type markers
# across conditions, showing both the expression level and the percentage of cells
# in a cluster expressing any given gene. 
# Here we plot 1-3 strong marker genes for each of our 13 clusters.
markers.to.plot <- c("Cdh5","Pecam1","Flt1","Vwf","Plvap","Kdr","EMCN","ptprb","KRT19",
                     "Epcam","PTPRC","LAPTM5","SRGN","CD14","CD68","CD3D","CCL5","Pdgfrb",
                     "Dcn")
markers.to.plot <- MouseGenes(BH,markers.to.plot, unique =T)
sort(markers.to.plot)
sdp <- SplitDotPlotGG(BH, genes.plot = rev(markers.to.plot),
                      cols.use = c("grey","blue"), x.lab.rot = T, plot.legend = T,
                      dot.scale = 8, do.return = T, grouping.var = "conditions")


# BH <- RenameIdentBack(BH)
# How many cells are in each cluster
table(BH@ident)
TSNEPlot(object = BH, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
  ggtitle("Major cell types")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle

#=====2.2 - A table with the number of cells of each cluster and subcluster, for both B6 and 129_B6 strains.
# We can also compare proportional shifts in the data. As can be seen in the barplot, 
freq_table <- prop.table(x = table(BH@ident, BH@meta.data[, "conditions"]), 
                         margin = 2)
barplot(height = freq_table)

freq_table


#====== 2.3 Compare cell type changes across conditions  ==========================================
# the two patients profiled have very different composition
lnames = load(file = "./data/BH_alignment.Rda")
idents <- as.data.frame(table(BH@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Unknown",
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Myeloid cells",
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Unknown",
                     "Epithelium",
                     "Mesenchymal Cells",
                     "Unknown",
                     "Myeloid cells",
                     "Lymphoid cells",
                     "Endothelium",
                     "Epithelium",
                     "Endothelium")
BH@ident <- plyr::mapvalues(x = BH@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
TSNEPlot(object = BH, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6,
         colors.use = c("#DB72FB","#619CFF","#D39200","#00B9E3",
                        "#D39200","#F8766D"))+
        ggtitle("Major cell types")+
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
                           pt.size = 1,label.size = 4)+
                ggtitle(samples[i])+
                theme(text = element_text(size=20),     #larger text including legend title							
                      plot.title = element_text(hjust = 0.5)) #title in middle
}
do.call(plot_grid, p)
