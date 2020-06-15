#This script is dedicated to the study of gene expression heterogeneity thorugh
#dimension reduction-based clustering, differential gene expression and
#pathway enrichment analysis. 

library(Seurat)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(htmlTable)
library(DescTools)
library(corto)


#Load the datsets.
load('kelly-rawcounts.rda')
k_rc <- rawcounts
load('be2c-rawcounts.rda')
b_rc <- rawcounts

#Remove genes that are unexpressed in both cells.
keep <- !(apply(k_rc,1,sum)==0 & apply(b_rc,1,sum)==0)
k_rc <- k_rc[keep,]
b_rc <- b_rc[keep,]


#Let us set up the Seurat object of the joint matrix for clustering.
data_matrix <- cbind(k_rc,b_rc)
dim(data_matrix)

#We will append a "k" to all colnames KELLY cells and a "b" to all colnames of BE(2)C cells to create two distinct groups in
#the Seurat object.
colnames(data_matrix) <- c(paste(colnames(data_matrix)[1:1105],'-k',sep=''),
                           paste(colnames(data_matrix)[1106:ncol(data_matrix)],'-b',sep=''))

seurat_expmat <-  CreateSeuratObject(data_matrix,names.field = 3,names.delim = '-')

#Clustering, variable feature selection and filtering.
seurat_expmat <- NormalizeData(seurat_expmat,normalization.method = 'LogNormalize',scale.factor = 10000)
seurat_expmat <- FindVariableFeatures(seurat_expmat, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(seurat_expmat)
seurat_expmat <- ScaleData(seurat_expmat, features = all.genes)

#Do Seurat-based clustering
PCA <- RunPCA(seurat_expmat, features = VariableFeatures(seurat_expmat))
DimPlot(PCA, reduction = "pca")+ggtitle('PCA clustering of the two cell lines')+
  theme(plot.title = element_text(hjust = 0.5))

#TSNE dimensional reduction.
set.seed(1)
TSNE <- RunTSNE(PCA)
png('tsne_seurat_clusterng.png',width=800,height=800, res=150)
DimPlot(TSNE, reduction = "tsne")+ggtitle('TSNE clustering KELLY vs SK-N-BE(2)-C')+  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#UMAP dimensional reduction.
set.seed(1)
UMAP <- RunUMAP(PCA, dims = 1:10)
png('umap_be2c_kelly.png',width=800,height=800, res=150)
DimPlot(UMAP, reduction = 'umap')+ggtitle('UMAP clustering KELLY vs SK-N-BE(2)-C')+  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Find differentially expressed features between 2 conditions through fold change calculations.
diff_markers <- FindMarkers(seurat_expmat, ident.1 = "b", ident.2 = "k",
                               logfc.threshold = 0,min.pct = 0,pseudocount.use = 0.25,slot='data')
p.adj <- diff_markers$p_val_adj #Check the adjusted p-values.
FC <- diff_markers$avg_logFC #Check the average log fold changes.
names(FC) <- rownames(diff_markers)

#Produce a MAplot of log average expression vs log average fold change.
data_matrix_mean <- apply(data_matrix,1,mean)[names(FC)]
png('mean_expr_fold_change.png',width=600,height = 600, res=120)
col=ifelse(p<0.01,'red','black')
plot(log10(data_matrix_mean),signature,pch=20,xlab='log10 mean expression',
     ylab='Average log fold change',main='Mean expression vs Fold change plot',col=col)
legend('topright',legend=c('significant','non-significant'),pch=20,col=c('red','black'))
dev.off()




#For differential gene expression analysis, we will treat the datasets separately.
k_seu <- CreateSeuratObject(k_rc)
b_seu <- CreateSeuratObject(b_rc) 
k_seu <- NormalizeData(k_seu)
b_seu <- NormalizeData(b_seu)
k_seu <- ScaleData(k_seu,features=rownames(k_seu)) 
b_seu <- ScaleData(b_seu,features=rownames(b_seu))
k_seu <- RunPCA(k_seu, features = rownames(k_seu))
b_seu <- RunPCA(b_seu, features = rownames(b_seu))


#Perform a t-test to draw gene signatures. 
manual_set <- cbind(GetAssayData(b_seu,slot='data'),GetAssayData(k_seu,slot='data'))
ps<-apply(manual_set,1,function(x){
  tt<-t.test(x[1:ncol(b_seu)],x[(ncol(b_seu)+1):ncol(manual_set)])
  return(cbind(p.adjust(tt$p.value),tt$statistic))
})
ps <- t(ps)
colnames(ps) <- c('p.value','t_statistic')

#Keep only the significant results.
signature_t <- ps[ps[,1]<0.05,2]


#Now draw out the list of genes of interest from the msigdb.
msigdbr_show_species()
mdf<-msigdbr(species="Homo sapiens") # Retrieve all human gene sets
mlist<-mdf %>% split(x=.$gene_symbol,f=.$gs_name) #%>% works as pipe operator.
plength<- sapply(mlist,length) #See how many genes for each path.

#Run fast GSEA (fgsea) on the drawn gene signatures.
results_t<- fgsea(pathways=mlist,stats=signature_t,nperm=1E5,minSize=5,maxSize=Inf,nproc=7)
significant_t <- results_t[results_t$padj<0.05]
View(significant_t)


#Plot the most enriched pathways with t-test
png('ES_before_regression.png',width=2000,height=1500, res=100)
par(mar=c(4,1,3,1))
plot_sign <- cbind(significant_t$NES,significant_t$padj)
rownames(plot_sign) <- significant_t$pathway
plot_sign <- plot_sign[order(plot_sign[,1],decreasing = F),]
to_plot <- rbind(plot_sign[1:15,],plot_sign[(nrow(plot_sign)-14):nrow(plot_sign),])[,1]
b <- barplot(to_plot,horiz = T,col=rep(c("cyan3","coral"),each=15),xlim=c(-4,2.5),
             main='Pathway enrichment scores before regression',xlab='NES',cex.main=2, cex.lab=1.6,font=2)
text(0,b[16:30],names(to_plot)[16:30],pos=2,font=2,cex=1.2)
text(0,b[1:15],names(to_plot)[1:15],pos=4,font=2,cex=1.2)
dev.off()


#Export html tables.
#Differential gene expression table
#Table differentially expressed genes
mean1 <- apply(manual_set[,1:ncol(b_seu)],1,mean)
mean2 <- apply(manual_set[,(ncol(b_seu)+1):ncol(manual_set)],1,mean)
logfc <- log((mean1+0.25)/(mean2+0.25))
logfc <- logfc[names(signature_t)]
diff_exp <- cbind(signature_t,logfc,ps[,1][names(signature_t)])
keep <- names(signature_t[!is.na(signature_t)])
colnames(diff_exp) <- c('t','logFC','p-value')
diff_exp <- diff_exp[order(diff_exp[,'t'],decreasing = T),]
diff_exp <- round(diff_exp[diff_exp[,'p-value']<0.05,],digits=2)
View(diff_exp)
table_diff <- htmlTable(diff_exp,cgroup='Differentially expressed genes')
write(table_diff,file='table_be2c_kelly_diff_exp.html')

#Enriched pathways table
sorted_significant <- as.data.frame(significant_t)
sorted_significant <- sorted_significant[order(sorted_significant$NES,decreasing = T),1:5]
sorted_significant_f <- signif(sorted_significant[,2:5],digits=2)
rownames(sorted_significant_f) <- sorted_significant[,1]

View(sorted_significant_f)

table_paths <- htmlTable(sorted_significant_f,ctable=T,rowlabel='Path names',
                   align.header='l',cgroup='Differentially expressed pathways')
write(table_paths,file='diff_pathways_kelly_be2c.html')

























