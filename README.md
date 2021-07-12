### Integration-of-scRNA-Seq-Data
Integrated Analysis of scRNA-Seq Data across species using R/Bioconductor libraries.

### Introduction to scRNA-seq integration
The joint analysis of two or more single-cell datasets poses unique challenges. In particular, identifying cell populations that are present across multiple datasets can be problematic under standard workflows. Seurat v4 includes a set of methods to match (or ‘align’) shared cell populations across datasets. These methods first identify cross-dataset pairs of cells that are in a matched biological state (‘anchors’), can be used both to correct for technical differences between datasets (i.e. batch effect correction), and to perform comparative scRNA-seq analysis of across experimental conditions.

### Seurat 
Seurat is a R package developed by Satija Lab, NYU, USA.

1. [satijalab](https://satijalab.org/seurat/)

### OSCA
“Orchestrating Single-Cell Analysis with Bioconductor”, a book that teaches users some common workflows for the analysis of single-cell RNA-seq data (scRNA-seq). 

2. [osca] https://bioconductor.org/books/release/OSCA/

### Integrated Analysis of scRNA-Seq Data across species

We integrate across two species as we have Single-cell Seq datasets for multiple species; also going to compare across different Seurat object by using Seurat library developed by Satijalab. We compared across two species human [human brain region-Human Lateral Geniculate Nucleus (LGN)] and mouse[Mouse Brain Region-Dorsolateral Geniculate Complex (LGd)]. 

In the first part, we do pre-processing and importing of data for the human single-cell RNA seq data. We do normal analysis as if it is one species and we created a t-SNE plot for the human and did cell type labelling. And I do the same thing for the mouse dataset. Then produce a combined t-SNE plot that integrates the two species to see the similarities and differences, also do a cluster labeling. And at the end, produce a frequency table to show frequencies of the different cell types across species

### For Integration section

We add metadata for each dataset and creat two Seurat object one called pbmc_human and another one is pbmc_mouse. So, we have two Seurat objects that are individually analyzed and then the next step we do is integration. So, for integration, we create a combined object and use the function called merge and then we merge two Seurat objects. We call the first object mouse and the second object human. For this next part, we do a split object and split it by protocol and the split is going to split into human and mouse. We want to split it by two species that we created from each individual analysis and then we do the anchoring step, so the anchoring step is the integration. 

### Result

![294](https://user-images.githubusercontent.com/65890522/123652351-82d3a080-d82c-11eb-90d5-b0cd6a1d261b.png)


### Reference
1. https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8
2. https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub
