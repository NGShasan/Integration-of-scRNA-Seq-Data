
###Project: Integrated Analysis of scRNA-Seq Data across species

library("SummarizedExperiment")
library("scran")
library("BiocGenerics")
library("SingleCellExperiment")
library("scRNAseq")
library("Seurat")
library("dplyr")
library("patchwork")
library("cellranger")
library("cowplot")
library("ggplot2")
library("pcaMethods")
library("PCAtools")

## set directory
getwd()
setwd("C:/Users/saifmehadi/Desktop/MS_Project")
rawPath <- "C:/Users/*****/Desktop/MS_Project/raw_data/"
saveIn <- "C:/Users/******/Desktop/MS_Project/graphics/"

##Load the count data of human
hum = read.table("human_count.csv",header=TRUE,sep=",")
rownames(hum)
colnames(hum)

#remove first column, i.e. cell indexes, from rawD
hum_outcome = data.frame(hum[,-1])
hum_outcome = data.frame(hum_outcome[,-1], row.names=hum_outcome[,1])

#reorder dataframes
order_hum = hum_outcome[ order(row.names(hum_outcome)), ]

#remove genes in less than 30 cells
rm_genes = which(rowSums(order_hum > 0) <30)
r1 = order_hum[-rm_genes,]

#create SingleCellExperiment object
sce = SingleCellExperiment(list(counts = data.matrix(r1)))
# To access the count data we just supplied
assay(sce, "counts")

#normalization
# using the scater package to compute a normalized and log-transformed representation of the initial primary data
sce <- scater::logNormCounts(sce)

#the log-transformed normalized data
logcounts(sce) 
assays(sce)

# We can access our column data with the colData() function
colData(sce) 

#the scater package contains the addPerCellQC() function that appends a lot of quality control data
sce <- scater::addPerCellQC(sce)
colData(sce)

#do clustering to reduce heterogeneity
clusters = quickCluster(sce, min.size=100)

# computeSumFactors function to scale normalization of sparse count data
sce = computeSumFactors (sce, cluster = clusters)

# calculate a PCA representation of our data using the runPCA() function from scater
sce <- scater::runPCA(sce)

#The reducedDims to store reduced dimensionality representations of the primary data obtained by PCA 
reducedDim(sce, "PCA")

#We can also calculate a tSNE representation using the scater package function runTSNE()
sce <- scater::runTSNE(sce, perplexity = 0.1)

# reducedDims to store reduced dimensionality representations of the primary data obtained by t-SNE
reducedDim(sce, "TSNE") 

#view the names of all our entries in the reducedDims
reducedDims(sce) 

#The sizeFactors() function allows us to get or set a numeric vector of per-cell scaling factors used for normalization 
sce <- scran::computeSumFactors(sce)
sizeFactors(sce)

#create Seurat object using scran object
s_obj= CreateSeuratObject(counts=log(counts(sce) +1))

#find variable features to reduce run time
s_obj=FindVariableFeatures(s_obj)

#regress out batch
s_obj_nobatch = ScaleData(s_obj)

##### STEP: Centering and scaling data matrix ######

#add metadata
s_obj_nobatch@meta.data[, "protocol"] <- "human"

#s_obj= ScaleData(s_obj, vars.to.regress = "orig.ident")
s_obj = RunPCA(s_obj_nobatch, features = VariableFeatures(s_obj))

#s_obj_nobatch = RunPCA (s_obj_nobatch, features=VariableFeatures(s_obj_nobatch))
s_obj = FindNeighbors(s_obj)

# Find cluster
s_obj = FindClusters(s_obj, resolution = 0.12)
table(Idents(s_obj))

# Run tSNE
s_obj = RunTSNE(s_obj) 

#Plot tSNE
DimPlot (s_obj, reduction = "tsne", label = TRUE)
ggsave('plot114.png')
dev.off()

###### data visualization by UMAP instead of tSNE ####

#reticulate::py_install(packages = 'umap-learn') ## Installation
s_obj <- RunUMAP(s_obj, dims = 1:10)

# visualize by UMAP plot
DimPlot(s_obj, reduction = "umap", label = TRUE)
ggsave('plot124.png')
dev.off()

#############################################################

#label clusters based on marker genes
current.cluster.ids = c(0, 1, 2, 3, 4, 5, 6)

new.cluster.ids = c("Glutamatergic neuron", "Glutamatergic neuron", "Oligodendrocyte", "GABAergic neuron",
                    "Glutamatergic neuron", "Astrocyte", "GABAergic neuron")

names(x = new.cluster.ids) <- levels(x = s_obj)
pbmc <- RenameIdents(object = s_obj, new.cluster.ids)

#Plot pbmc
DimPlot(object = pbmc, reduction = "tsne", label = TRUE, pt.size = 0.8) + NoLegend()
ggsave('plot140.png')
dev.off()

#Plot pbmc without legend
DimPlot(object = pbmc, reduction = "tsne", group.by = "protocol", pt.size = 0.5) + NoLegend()
ggsave('plot144.png')
dev.off()

pbmc_human = pbmc

######## mouse ###########

##import the count data
mus = read.table("mouse_count.csv",header=TRUE,sep=",")
#remove first column, i.e. cell indexes, from rawD
mus_outcome = data.frame(mus[,-1])
mus_outcome = data.frame(mus_outcome[,-1], row.names=mus_outcome[,1])

#reorder dataframes
order_mus = mus_outcome[ order(row.names(mus_outcome)), ]
#remove genes in less than 30 cells
rm_genes = which(rowSums(order_mus > 0) <30)
r2 = order_mus[-rm_genes,]

#create SingleCellExperiment object
sce = SingleCellExperiment(list(counts = data.matrix(r2)))
# To access the count data we just supplied, we can do any one of the following:
assay(sce, "counts")
counts(sce)

# using the scater package to compute a normalized and log-transformed representation of the initial primary data
sce <- scater::logNormCounts(sce)
#the log-transformed normalized data
logcounts(sce) 
assays(sce)

#the scater package contains the addPerCellQC() function that appends a lot of quality control data
sce <- scater::addPerCellQC(sce)
colData(sce)

#do clustering to reduce heterogeneity
clusters = quickCluster(sce, min.size=100)
# computeSumFactors function to scale normalization of sparse count data
sce = computeSumFactors(sce, cluster = clusters)

#normalization
sce <- scater::logNormCounts(sce)
#we can calculate a PCA representation of our data using the runPCA() function from scater.
sce <- scater::runPCA(sce)

#The reducedDims to store reduced dimensionality representations of the primary data obtained by PCA and t-SNE
reducedDim(sce, "PCA")
#We can also calculate a tSNE representation using the scater package function runTSNE()
sce <- scater::runTSNE(sce, perplexity = 0.1)

#The reducedDims to store reduced dimensionality representations of the primary data obtained by t-SNE
reducedDim(sce, "TSNE")
reducedDims(sce) #view the names of all our entries in the reducedDims

#The sizeFactors() function allows us to get or set a numeric vector of per-cell scaling factors used for normalization 
sce <- scran::computeSumFactors(sce)
sizeFactors(sce)

#create Seurat object using scran object
s_obj= CreateSeuratObject(counts = log(counts(sce) +1))
#find variable features to reduce run time
s_obj=FindVariableFeatures(s_obj)
#regress out batch
s_obj_nobatch = ScaleData(s_obj)


#### STEP: Centering and scaling data matrix ####
#add metadata
s_obj_nobatch@meta.data[, "protocol"] <- "mouse"
#s_obj= ScaleData(s_obj, vars.to.regress = "orig.ident")
s_obj = RunPCA(s_obj_nobatch, features = VariableFeatures(s_obj))
#s_obj_nobatch = RunPCA (s_obj_nobatch, features=VariableFeatures(s_obj_nobatch))
s_obj = FindNeighbors(s_obj)

# Finding cluster
s_obj = FindClusters(s_obj, resolution = 0.12)

## Find Indent numbers
table(Idents(s_obj))

# Run tSNE
s_obj = RunTSNE(s_obj, dims = 1:10) 
# Plot tSNE
DimPlot (s_obj, reduction = "tsne", label = TRUE)
ggsave('plot229.png')
dev.off()

######## Data visualization by UMAP #######

# Run UMAP
s_obj <- RunUMAP(s_obj, dims = 1:10)

# UMAP plot
DimPlot(s_obj, reduction = "umap", label = TRUE) 
ggsave('plot239.png')
dev.off()

############################################################################

# label clusters based on marker genes

current.cluster.ids = c(0, 1, 2, 3, 4, 5, 6)
new.cluster.ids = c("Glutamatergic neuron", "GABAergic neuron", "Glutamatergic neuron", "GABAergic neuron",
                    "Glutamatergic neuron","Oligodendrocyte","GABAergic neuron")

names(x = new.cluster.ids) <- levels(x = s_obj) 
pbmc <- RenameIdents(object = s_obj, new.cluster.ids) 
DimPlot(object = pbmc, reduction = "tsne", label = TRUE, pt.size = 0.7) + NoLegend()
ggsave('plot253.png')
dev.off()

pbmc_mouse = pbmc

#######Integration of human and mouse##########

# Combined human object and mouse object and separated by group protocol
pbmc.combined = merge(pbmc_mouse, y = pbmc_human, add.cell.ids = c("mouse", "human"), project = "protocol")
# Split object by protocol
species.list <- SplitObject(pbmc.combined, split.by = "protocol")

# Split between human and mouse
reference.list <- species.list[c("human", "mouse")]
# Integration to find anchoring point between two species and finding anchoring points using first 30 dimensions
species.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30) 

###### I did this to increase my storage capacity to 7GB
#To know the current storage capacity
memory.limit()
#To increase the storage capacity
memory.limit(size=56000)

# We then pass these anchors to the IntegrateData function, which returns a Seurat object.
species.integrated <- IntegrateData(anchorset = species.anchors, dims = 1:30)
# After running IntegrateData, the Seurat object now contain a new Assay with the integrated expression matrix.
DefaultAssay(species.integrated) <- "integrated"

##### Run the standard workflow for visualization and clustering####

# Scaling the integrated data
species.integrated <- ScaleData(species.integrated, verbose = FALSE)
# Run PCA
species.integrated <- RunPCA(species.integrated, npcs = 30, verbose = FALSE)
#Run t-SNE
species.integrated <- RunTSNE(species.integrated, dims = 1:30)

# Plot t-SNE and show the differences based on protocol labeled human and mouse
p1 <- DimPlot(species.integrated, reduction = "tsne", group.by = "protocol")
#Plot the integrated object
plot(p1)
ggsave('plot294.png')
dev.off()

# Plot the integrated object separated by clusters
p2 <- DimPlot(species.integrated, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +NoLegend()
plot(p2)
ggsave('plot300.png')
dev.off()

# Show two t-SNE plot side by side (one by species and another one by Seurat cluster)
plot_grid(p1, p2)
ggsave('plot305.png')
dev.off()

############## Run UMAP instead of t-SNE) ####################################
# species.integrated.umap <- RunUMAP(species.integrated, reduction = "pca", dims = 1:30)
# p3 <- DimPlot(species.integrated.umap, reduction = "umap", group.by = "protocol")
# p4 <- DimPlot(species.integrated.umap, reduction = "umap", group.by = "umap_clusters", label = TRUE, 
#               repel = TRUE) + NoLegend()
# plot_grid(p3, p4)
# ggsave('plot303.png')
# dev.off()

#################annotate based on markers###############

#Cluster cell Type labeling using some genetic markers to determine integrated object

#gabaergic
FeaturePlot(object = species.integrated, features= c("2572", "2571"), min.cutoff = "q9", pt.size = 0.5)
ggsave('plot323.png')
dev.off()

#oligodendrocyte
FeaturePlot(object = species.integrated, features = c("4155", "4340","6507"), min.cutoff = "q9", pt.size = 0.5)
ggsave('plot328.png')
dev.off()

#glutametergic
FeaturePlot(object = species.integrated, features = c("57084", "2904","2744"), min.cutoff = "q9", pt.size = 0.5)
ggsave('plot333.png')
dev.off()

#Astrocyte
FeaturePlot(object = species.integrated, features = c("361"), min.cutoff = "q9", pt.size = 0.5)
ggsave('plot338.png')
dev.off()

# FeaturePlot
FeaturePlot(object = species.integrated, features = c("5156"), min.cutoff = "q9", pt.size = 0.5)
ggsave('plot343.png')
dev.off()

# Assigned cell ids which I want from different clusters
new.ident <- c("Glutamatergic neuron", "GABAergic neuron", "Oligodendrocyte", "GABAergic neuron", "?", "Oligodendrocyte", "GABAergic neuron")
levels(x = species.integrated)

# used the labels to in t-SNE plot
p2 <- DimPlot(species.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + NoLegend()
# plot combined integrated object
plot(p2)
ggsave('plot354.png')
dev.off()

# Show two separated plot, one is separated by species type and another one is by cell ID type
plot_grid(p1, p2)
ggsave('plot359.png')
dev.off()

#### create a proportion table or frequency table to compare####

# differential gene expression across species
freq_table <- prop.table(x = table(species.integrated@active.ident, species.integrated@meta.data[, "protocol"]),margin = 2)
# Show the proportion table in barplot
barplot(height = freq_table)
ggsave('plot368.png')
dev.off()

# Show the actual proportion breakdown for human and mouse across all the cell IDs respectively
freq_table

