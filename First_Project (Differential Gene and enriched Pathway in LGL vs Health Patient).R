.libPaths("/standard/cphg-RLscratch/share/R/local/4.3.1/")
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
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
#ggplot(pbmc, aes=(x=diagnosis)) + geom_bar(aes(fill=celltype))

#UMAP clustering of Normal cell type 
#DimPlot(pbmc_Normal)

#Same Procedure for TLGL Patients
Idents(object = pbmc_TLGL) <- "celltype"
Idents(pbmc_TLGL)
table(Idents(pbmc_TLGL))
#DimPlot(pbmc_TLGL)

#Create a new class that combine celltype and diagnosis status together
pbmc$celltype.stim <-paste(pbmc$celltype, pbmc$diagnosis, sep="_")
Idents(pbmc) <- "celltype.stim"

#Sample differential Gene Analaysis (Treg)
Treg.de <-FindMarkers(pbmc, ident.1='Treg_TLGL', ident.2="Treg_Normal", latent.var=c("batch","donor","tissue"), test.use="MAST", verbose=TRUE)
Treg.de$celltype <- rep('Treg', each=nrow(Treg.de))


#Create a meta dataset that combine RNA differential information of all celltype that present in both Normal and TLGL Patients.
#Then plot of a volcano plot that colored by celltype. Label the signficant genes.


Meta.de<- data.frame()
latent_features<- 
for(cell in unique(pbmc$celltype)){
  if(cell == 'TLGLL' | cell == 'NKLGLL'){
    next
  }
  
  Normal_Iden<-paste(cell, "_Normal", sep="")
  TLGL_Iden<- paste(cell, "_TLGL", sep="")
  
  
  
  Temp.de<-FindMarkers(pbmc, ident.1=Normal_Iden, ident.2=TLGL_Iden, test.use='MAST', latent.var=c("batch","donor","tissue"), verbose=TRUE)
  Temp.de$celltype <- rep(cell, nrow(Temp.de))
  Meta.de<- rbind(Meta.de, Temp.de)
}
Meta.de$delabel <- NA
Meta.de$delabel[which(Meta.de$p_val_adj<10^(-10))] <- row.names(Meta.de)[which(Meta.de$p_val_adj<10^(-10))]

#Making Volcano Plot
ggplot(data=Meta.de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=celltype, label=delabel)) + 
  geom_point() + 
  geom_text_repel()
  theme_minimal() +
  geom_hline(yintercept=10, col="red")

#ClusterProfiler

#Prepare Background Annotation
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)


#import the metadata list used in the differential gene program
original_gene_list <- Meta.de$avg_log2FC
names(original_gene_list)<-Meta.de$X
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
#Create GSE
gse<-gseGO(geneList=gene_list, 
      ont ="ALL", 
      keyType = "GO", 
      nPerm = 3, 
      minGSSize = 3, 
      maxGSSize = 400, 
      pvalueCutoff = 10^(-10), 
      verbose = TRUE, 
      OrgDb = organism, 
      pAdjustMethod = "none")

