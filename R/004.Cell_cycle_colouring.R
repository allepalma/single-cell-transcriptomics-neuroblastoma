#In this script we will use the basic functionalities of Seurat to estimate and
#remove the effects of the cell cycle from the expression counts analysed in the
#previous scripts.

library(ggplot2)
library(Seurat)
library(msigdbr)
library(fgsea)
library(corto)
library(htmlTable)

#Load the dataset.
load('kelly-rawcounts.rda')
k_rc <- rawcounts
load('be2c-rawcounts.rda')
b_rc <- rawcounts

#The mitigation of the effects of the cell cycle on the expression counts is performed 
#treating the two datasets separately.

#First, import a list of cell cycle markers of S phase and G2M phase.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


#Do separate cell cycle regression analysis
k_seu <- CreateSeuratObject(k_rc)
k_seu <- NormalizeData(k_seu)
k_seu <- FindVariableFeatures(k_seu, selection.method = "vst",nfeatures = 5000)
k_seu <- ScaleData(k_seu, features = rownames(k_seu))
#Assign to each cell a cell cycle score for s and g2
k_seu <- CellCycleScoring(k_seu, s.features = s.genes, 
                           g2m.features = g2m.genes, set.ident = TRUE)

PCA_k <- RunPCA(k_seu, features = VariableFeatures(k_seu))

png('KELLY_before_regression.png',width=700,height=700, res=120)
plot <- DimPlot(PCA_k,pt.size = 2)
plot + ggtitle('Gene expression before regression in KELLY cells')+theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Regress out effects of cell cycle
k_seu_reg <- ScaleData(k_seu, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(k_seu))
k_seu_reg <- RunPCA(k_seu_reg, features = VariableFeatures(k_seu_reg), nfeatures.print = 10)
png('KELLY_after_regression.png',width=700,height=700, res=120)
DimPlot(k_seu_reg,pt.size = 2)+ggtitle('Gene expression after regression in KELLY cells')+theme(plot.title = element_text(hjust = 0.5))
dev.off()


#Repeat the pipeline for BE(2)C cells.
b_seu <- CreateSeuratObject(b_rc)
b_seu <- NormalizeData(b_seu)
b_seu <- FindVariableFeatures(b_seu, selection.method = "vst",nfeatures = 5000)
b_seu <- ScaleData(b_seu, features = rownames(b_seu))
b_seu <- CellCycleScoring(b_seu, s.features = s.genes, 
                          g2m.features = g2m.genes, set.ident = TRUE)
PCA_b <- RunPCA(b_seu, features = VariableFeatures(b_seu))

png('BE2C_before_regression.png',width=700,height=700, res=120)
plot <- DimPlot(PCA_b,pt.size = 2)
plot + ggtitle('Gene expression before regression in BE(2)C cells')+theme(plot.title = element_text(hjust = 0.5))
dev.off()

b_seu_reg <- ScaleData(b_seu, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(b_seu))
b_seu_reg <- RunPCA(b_seu_reg, features = VariableFeatures(b_seu_reg), nfeatures.print = 10)

png('BE2C_after_regression.png',width=700,height=700, res=120)
DimPlot(b_seu_reg,pt.size = 2)+ggtitle('Gene expression after regression in BE(2)C cells')+theme(plot.title = element_text(hjust = 0.5))
dev.off()



#Check the proportion of cells assigned to both phases of the cell cycle.
phases_k <- k_seu$Phase
phases_b <- b_seu$Phase
perc_k <- table(phases_k)/sum(table(phases_k))*100
perc_b <- table(phases_b)/sum(table(phases_b))*100


png(file='cell-cycle_phases.png',width = 700,height=700, res=90)
par(mfrow=c(1,2))
k <-plot(phases_k, col=c('coral','chartreuse3','cyan3'),main='KELLY cells phases',ylim=c(0,600),ylab='Number of cells'
         ,cex.lab=1.5,cex.main=1.5,font=2)
text(k,table(phases_k)+20,labels =paste0(round(perc_k),rep('%',3)),font=2)
b <- plot(phases_b, col=c('coral','chartreuse3','cyan3'),main='BE(2)-C cells phases',ylim=c(0,600),ylab='Number of cells'
          ,cex.lab=1.5,cex.main=1.5,font=2)
text(b,table(phases_b)+20,labels =paste0(round(perc_b),rep('%',3)),font=2)
dev.off()


#Check if the two SK-N-BE(2)-C clusters are preserved.

#First we load the two external objects containing the cluster labels of the two populations
#of BE(2)C cells.
load('cluster1.rda')
load('cluster2.rda')

#Colour the dimension reduction plot of the BE(2)C cells after the regression on the base of the cluster
#labels loaded before. 
PC1 <- Embeddings(b_seu_reg)[,1] 
PC2 <- Embeddings(b_seu_reg)[,2] 
col <- ifelse(names(PC1)%in%names(cluster1),'cornflowerblue','chartreuse3')

png('PCA_sep_be2c_after_reg.png',width=700,height=700)
plot(PC1,PC2,pch=20,col=col,xlab ='PC1',ylab='PC2',main='Clusters of SK-N-BE(2)-C after cell cycle regression',
     cex.lab=1.5,font.lab=2, font=2,cex.main=2.2,cex=2)
legend('topright',legend = c('group1','group2'),pch=20,col=c('cornflowerblue','chartreuse3'),cex=2)
dev.off()











