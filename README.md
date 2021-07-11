### Integration-of-scRNA-Seq-Data
Integrated Analysis of scRNA-Seq Data across species using R/Bioconductor libraries.

### Seurat 
Seurat is a R package developed by Satija Lab, NYU, USA.

### Integrated Analysis of scRNA-Seq Data across species using Seurat

I am going to integrate across two species as I have Single-cell Seq datasets for multiple species; also going to compare across different Seurat object 
by using Seurat library developed by Satijalab. I am comparing across two species human [human brain region-Human Lateral Geniculate Nucleus (LGN)] 
and mouse[Mouse Brain Region-Dorsolateral Geniculate Complex (LGd)]. I have separate count matrix datasets from each of this species and I do a combined 
analysis. Prior to running the analysis, I have filtered the data on homologous genes of human and mouse.

In the first part I do pre-processing and importing of data for the human single-cell RNA seq data. I do normal analysis as if it is one species that I am 
comparing and I create a t-SNE plot for the human and I do cell type labelling. And I do the same thing for the mouse dataset. I am going to produce a combined
t-SNE plot that integrates the two species to see the similarities and differences, also do a cluster labeling. 
And at the end, produce a frequency table to show frequencies of the different cell types across species

### For Integration section

I add metadata for each dataset and created two Seurat object one called pbmc_human and another one is pbmc_mouse. 
So, I have two Seurat objects that are individually analyzed and then the next step I want to do is integration. So, for integration, 
I am going to create a combined object and use the function called merge and then I am going to merge my two Seurat objects. I call the first object mouse 
and the second object human. For this next part, I am doing a split object and splitting it by protocol and the split is going to split into human and mouse. 
I want to split it by two species that I created from each individual analysis and then I am going to do the anchoring step, so the anchoring step is the 
integration. 

### Result

![294](https://user-images.githubusercontent.com/65890522/123652351-82d3a080-d82c-11eb-90d5-b0cd6a1d261b.png)

### Contributor
- Mehadi Hasan <isaifmehadi@gmai.com>

### Licence and Copyright
© Hasan, Mehadi

### Reference
[satijalab](https://satijalab.org/seurat/)
