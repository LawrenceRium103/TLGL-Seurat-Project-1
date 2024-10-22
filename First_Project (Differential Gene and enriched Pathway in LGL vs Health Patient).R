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


#Do it for every Celltype
library(BiocParallel)
register(SerialParam())

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

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
  

  # title_volcano<- paste(cell,"Volcano_plot.png", sep="_")
  # ggplot(data=Temp.de, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) + 
  #   geom_point() + 
  #   geom_text_repel()+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #   geom_hline(yintercept=-log10(0.05), col="red")
  # 
  # ggsave(title_volcano,width=10, height=7, dpi=300,limitsize = FALSE)
  
  Meta.de<-rbind(Meta.de, Temp.de)
  
  original_gene_list <- Temp.de$avg_log2FC
  names(original_gene_list)<-rownames(Temp.de)
  gene_list<-na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  ID<-mapIds(org.Hs.eg.db, keys = names(gene_list),
             column = "ENTREZID", keytype = "SYMBOL")
  names(gene_list)<-ID
  
  # #Create GSE
  # gse<-gseGO(geneList=gene_list, 
  #            ont ="ALL",
  #            keyType = "ENTREZID",
  #            minGSSize = 5, 
  #            maxGSSize = 500, 
  #            pvalueCutoff = 0.05, 
  #            verbose = TRUE, 
  #            OrgDb = organism, 
  #            pAdjustMethod = "none")
  # 
  # require(DOSE)
  # title_dot<-paste(cell,"_dotplot.png",sep="")
  # dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  # ggsave(title_dot,width=7, height=10, dpi=300, limitsize = FALSE)
  
}
data<-read.csv("Meta.csv")
data = data[-1,]
rownames(data) <- data$X
data=data[,-1]

original_gene_list <- rownames(data)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
ID<-mapIds(org.Hs.eg.db, keys = gene_list,
           column = "ENTREZID", keytype = "SYMBOL")

rownames(data)<-ID
data$ENTREZID<-ID
data_by_celltype<-group_split(data,data$celltype)
cell_names<-NA

for (cell in data_by_celltype){
    cell_names<-c(cell_names,cell$celltype[1])  
}
cell_names<-cell_names[-1]
names(data_by_celltype)<-cell_names

nested_list<-vector("list", 21)  
for (i in 1:21) {
  # Create an inner list for multiples of i
  inner_list <- data_by_celltype[[i]][[9]]
  # Assign the inner list to the outer list
  nested_list[[i]]<-inner_list
}

# Name the outer list elements
names(nested_list) <- cell_names

ck <- compareCluster(geneCluster = nested_list, fun = enrichKEGG)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck)

ck_B_cell <- compareCluster(geneCluster = nested_list[1:3], fun = enrichKEGG)
ck_B_cell <- setReadable(ck_B_cell, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck_B_cell) 

ck_CD4_T_cell <- compareCluster(geneCluster = nested_list[6:9], fun = enrichKEGG)
ck_CD4_T_cell <- setReadable(ck_CD4_T_cell, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

ck_CD8_T_cell <- compareCluster(geneCluster = nested_list[10:13], fun = enrichKEGG)
ck_CD8_T_cell <- setReadable(ck_CD8_T_cell, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(ck)


