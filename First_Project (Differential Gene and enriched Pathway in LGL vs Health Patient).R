.libPaths("/standard/cphg-RLscratch/share/R/local/4.3.1/")
library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)

pbmc <- readRDS('/standard/cphg-RLscratch/ar7jq/LGL_single_cell/lawrence/seurat.rds')

colnames(pbmc)
rownames(pbmc)

dim(pbmc)

num_cell<-ncol(pbmc)
num_feature<-nrow(pbmc)

VariableFeatures(pbmc)
Idents(pbmc)
table(Idents(pbmc))

#Separate the metadata into diagnosis group
pbmc_Normal<-subset(x=pbmc, subset= diagnosis=='Normal')
pbmc_TLGL<-subset(x=pbmc, subset= diagnosis=='TLGL')

#Using Frequency Table to determine the abundance of each cell type in Normal Samples
Idents(object = pbmc_Normal) <- "celltype"
Idents(pbmc_Normal)
table(Idents(pbmc_Normal))

#Make a barplot of cell plot where 
ggplot(pbmc, aes=(x=diagnosis)) + geom_bar(aes(fill=celltype))

#UMAP clustering of Normal cell type 
DimPlot(pbmc_Normal)

#Same Procedure for TLGL Patients
Idents(object = pbmc_TLGL) <- "celltype"
Idents(pbmc_TLGL)
DimPlot(pbmc_TLGL)

#Create a new class that combine celltype and diagnosis status together
pbmc$celltype.stim <-paste(pbmc$celltype, pbmc$diagnosis, sep="_")
Idents(pbmc) <- "celltype.stim"

#Sample differential Gene Analaysis (Treg)
Treg.de <-FindMarkers(pbmc, ident.1='Treg_TLGL', ident.2="Treg_Normal", test.use="MAST", verbose=TRUE)
Treg.de$celltype <- rep('Treg', each=nrow(Treg.de))

unique(pbmc$celltype)
unique(pbmc_Normal$celltype)
unique(pbmc_TLGL$celltype)


#Create a meta dataset that combine RNA differential information of all celltype that present in both Normal and TLGL Patients.
#Then plot of a volcano plot that colored by celltype. Label the signficant genes.

Idents(object = pbmc) <- "celltype"

for(cell in unique(pbmc$celltype)){
  if(cell == 'TLGLL' | cell == 'NKLGLL'){
    next
  }
  Normal <-subset(x=pbmc_Normal, idents=cell)
  TLGL <- subset(x=pbmc_TLGL, idents=cell)
  FindMarkers(pbmc, ident.1='Normal', ident.2='TLGL', test.use='MAST', verbose=TRUE)
 
}

CD14_Mono.de <-FindMarkers(pbmc, ident.1 = 'CD14 Mono_Normal', ident.2="CD14 Mono_TLGL", test.use="MAST", verbose=TRUE)
CD14_Mono_Treg.de <-merge()

#ClusterProfiler

#Prepare Background Annotation
organism = "org.Hs.eg.db"
BiocManager:: install (organism, character.only = TRUE)
library(organism, character.only = TRUE)

#import the metadata list used in the differential gene program
original_gene_list <- 

