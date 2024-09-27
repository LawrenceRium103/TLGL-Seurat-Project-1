.libPaths("/standard/cphg-RLscratch/share/R/local/4.3.1/")
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(doMC)


pbmc <- readRDS('/standard/cphg-RLscratch/ar7jq/LGL_single_cell/lawrence/seurat.rds')
DefaultAssay(object = pbmc)<-'RNA'

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
Treg.de <-FindMarkers(pbmc, ident.1='Treg_TLGL', ident.2="Treg_Normal", latent.var="batch", test.use="MAST", verbose=TRUE)
Treg.de$celltype <- rep('Treg', each=nrow(Treg.de))
Treg.de$diffexpressed <- "NO"
Treg.de$diffexpressed[Treg.de$avg_log2FC > 0.6 & Treg.de$p_val < 0.05] <- "UP"
Treg.de$diffexpressed[Treg.de$avg_log2FC < -0.6 & Treg.de$p_val < 0.05] <- "DOWN"

Treg.de$delabel <- NA
Treg.de$delabel[which(Treg.de$p_val<0.05)] <- row.names(Treg.de)[which(Treg.de$p_val<0.05)]

ggplot(data=Treg.de, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  geom_text_repel()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=-log10(0.05), col="red")

#Sample GSE Analysis
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

original_gene_list <- Treg.de$avg_log2FC
names(original_gene_list)<-rownames(Treg.de)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
ID<-mapIds(org.Hs.eg.db, keys = names(gene_list),
           column = "ENTREZID", keytype = "ENSEMBL")
names(gene_list)<-ID
#Create GSE
gse<-gseGO(geneList=gene_list, 
           ont ="ALL",
           keyType = "ENTREZID",
           minGSSize = 5, 
           maxGSSize = 500, 
           pvalueCutoff = 0.05, 
           verbose = TRUE, 
           OrgDb = organism, 
           pAdjustMethod = "none")


require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)



#Do it for every Celltype
library(BiocParallel)
register(SerialParam())

Meta.de<-NA
for(cell in unique(pbmc$celltype)){
  if(cell == 'TLGLL' | cell == 'NKLGLL'){
    next
  }
  
  Normal_Iden<-paste(cell, "_Normal", sep="")
  TLGL_Iden<- paste(cell, "_TLGL", sep="")
  
  
  
  Temp.de<-FindMarkers(pbmc, ident.1=Normal_Iden, ident.2=TLGL_Iden, test.use='MAST', latent.var="batch", verbose=TRUE)
  Temp.de$celltype <- rep(cell, nrow(Temp.de))

  
  Temp.de$diffexpressed <- "NO"
  Temp.de$diffexpressed[Temp.de$avg_log2FC > 0.6 & Temp.de$p_val < 0.05] <- "UP"
  Temp.de$diffexpressed[Temp.de$avg_log2FC < -0.6 & Temp.de$p_val < 0.05] <- "DOWN"
  
  Temp.de$delabel <- NA
  Temp.de$delabel[which(Temp.de$p_val<0.05)] <- row.names(Temp.de)[which(Temp.de$p_val<0.05)]
  
  title_volcano<- paste(cell,"Volcano_plot.png", sep="_")
  ggplot(data=Temp.de, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    ggtitle(title)+
    geom_text_repel()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    geom_hline(yintercept=-log10(0.05), col="red")
  
  ggsave(title_dot,width=10, height=7, dpi=300,limitsize = FALSE)
  
  Meta.de<-rbind(Meta.de, Temp.de)
  
  original_gene_list <- Temp.de$avg_log2FC
  names(original_gene_list)<-rownames(Temp.de)
  gene_list<-na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  ID<-mapIds(org.Hs.eg.db, keys = names(gene_list),
             column = "ENTREZID", keytype = "SYMBOL")
  names(gene_list)<-ID
  #Create GSE
  gse<-gseGO(geneList=gene_list, 
             ont ="ALL",
             keyType = "ENTREZID",
             minGSSize = 5, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

  require(DOSE)
  title_dot<-paste(cell,"_dotplot.png",sep="")
  dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  ggsave(title_dot,width=7, height=10, dpi=300, limitsize = FALSE)
  
}

#Making Volcano Plot

#ClusterProfiler

#Prepare Background Annotation



