#This script is dedicated to the study of expression of heterogeneity of two
#sub-populations of SK-N-BE(2)-C cells. We first individuate the two groups
#through unsupervised k-means clustering. Then we submit the result to basic pathway enrichment analysis.


library(Seurat)
library(msigdbr)
library(htmlTable)
library(fgsea)
library(scatterplot3d)

#Import the BE(2)C dataset.
load('be2c-rawcounts.rda')
b_rc <- rawcounts
load('be2c-tpms.rda')
b_tpms <- tpms

b_rc <- b_rc[apply(b_rc,1,sum)!=0,]


#Create a Seurat object.
b_seu <- CreateSeuratObject(b_rc)
b_seu <- NormalizeData(b_seu)
b_seu <- ScaleData(b_seu)
b_seu <- FindVariableFeatures(b_seu,nfeatures=10000)
b_seu <- RunPCA(b_seu)

#Plot the PCA.
DimPlot(b_seu)


#Fetch the normalized data.
b_norm <- GetAssayData(b_seu,'data')
b_norm <- as.matrix(b_norm)
b_norm_scale <- scale(b_norm,T,T)

#Launch k-means clustering.
set.seed(12)
kmeans <- kmeans(t(b_norm_scale), centers=2, iter.max = 1000, nstart = 1)

#Isolate the two cluseters.
clusters <- kmeans$cluster
cluster1 <- clusters[clusters==1]
cluster2 <- clusters[clusters==2]


#Extract the PCAs and check yo what extent the two groups represent the
#dimensional reduction-based clusters.
PC1 <- Embeddings(b_seu)[,1]
PC2 <- Embeddings(b_seu)[,2]
PC3 <- Embeddings(b_seu)[,3]
col <- ifelse(names(PC1)%in%names(cluster1),'cornflowerblue','chartreuse3')

png('PCA_be2c.png',width=700,height=700)
plot(PC1,PC2,pch=20,col='coral',xlab ='PC1',ylab='PC2',main='Clusters of SK-N-BE(2)-C before k-means',cex.lab=1.5,font.lab=2,
     font=2,cex.main=2.2,cex=2)
dev.off()

png('PCA_sep_be2c.png',width=700,height=700)
plot(PC1,PC2,pch=20,col=col,xlab ='PC1',ylab='PC2',main='Clusters of SK-N-BE(2)-C after k-means',cex.lab=1.5,font.lab=2,
     font=2,cex.main=2.2,cex=2)
text(PC1[name],PC2[name],labels = name,font = 2)
legend('topright',legend = c('group1','group2'),pch=20,col=c('cornflowerblue','chartreuse3'),cex=2)
dev.off()

png('3d_pca.png',width=800,height=800)
scatterplot3d(PC1,PC2,PC3,color=col,pch=20,xlab = 'PC1',ylab= 'PC2' , zlab='PC3',main='3D principal component analysis',
              cex.lab=1.3,angle=130,tick.marks = T,cex.symbols = 1.5,cex.main=2.5,ylim=c(-25,40),
              xlim = c(-50,45),zlim = c(-40,25))
dev.off()

#Now we have got the two groups, call them b1 and b2
b1 <- b_norm[,names(cluster1)]
b2 <- b_norm[,names(cluster2)]

#Differential expression analysis
#Differential expression analysis with a t-test
manual_set <- cbind(b1,b2)
ps<-apply(manual_set,1,function(x){
  tt<-t.test(x[1:ncol(b1)],x[(ncol(b1)+1):ncol(manual_set)])
  return(cbind(p.adjust(tt$p.value),tt$statistic))
})
ps <- t(ps)
colnames(ps) <- c('p.value','t_statistic')
signature_b <- ps[ps[,1]<0.05,2]


#Perform pathway enrichment analysis.
mdf<-msigdbr(species="Homo sapiens") # Retrieve all human gene sets
mlist<-mdf %>% split(x=.$gene_symbol,f=.$gs_name) #%>% works as pipe operator.

#Run fgsea.
results_b<- fgsea(pathways=mlist,stats=signature_b,nperm=1E5,minSize=5,maxSize=Inf,nproc=7)
significant_b <- results_b[results_b$padj<0.05 ]
View(significant_b)

  
#Graph the significant pathways
png('b_score.png',width=3000,height=2000, res=90)
par(mar=c(4,1,3,1))
plot_sign <- significant_b$NES
names(plot_sign) <- significant_b$pathway
plot_sign <- plot_sign[order(plot_sign,decreasing = F)]
to_plot <- c(plot_sign[1:15],plot_sign[(length(plot_sign)-14):length(plot_sign)])
b <- barplot(to_plot,horiz = T,col=rep(c('cornflowerblue','chartreuse3'),each=15),xlim=c(-5,5),
             main='Pathway enrichment scores before regression',xlab='NES',cex.main=2, cex.lab=1.6,font=2)
text(0,b[16:30],names(to_plot)[16:30],pos=2,font=2,cex=1.2)
text(0,b[1:15],names(to_plot)[1:15],pos=4,font=2,cex=1.2)
legend(-4,b[29],legend=c('group 1','group 2'),pch=20,col=c('cornflowerblue','chartreuse3'),cex=2)
dev.off()
  

#Produce tables
#Table differentially expressed genes
mean1 <- apply(manual_set[,1:ncol(b1)],1,mean)
mean2 <- apply(manual_set[,(ncol(b1)+1):ncol(manual_set)],1,mean)
logfc <- log((mean1+0.25)/(mean2+0.25))
diff_exp <- cbind(signature_b,logfc,ps[,1])
keep <- names(signature_b[!is.na(signature_b)])
colnames(diff_exp) <- c('t','logFC','p-value')
diff_exp <- diff_exp[order(diff_exp[,'t'],decreasing = T),]
diff_exp <- round(diff_exp[diff_exp[,'p-value']<0.05,],digits=2)
View(diff_exp)
table_diff <- htmlTable(diff_exp, cgroup='Differentially expressed genes in BE(2)C groups')
write(table_diff,file='table_be2c_diff_exp.html')

#Table of differentially expressed pathways
sorted_significant <- as.data.frame(significant_b)
sorted_significant <- sorted_significant[order(sorted_significant$NES,decreasing = T),1:5]
sorted_significant_b <- signif(sorted_significant[,2:5],digits=2)
rownames(sorted_significant_b) <- sorted_significant[,1]
table_paths <- htmlTable(sorted_significant_b,ctable=T,rowlabel='Path names',
                   align.header='l', cgroup='Differentially expressed pathways in BE(2)C groups')
write(table_paths,file='b_table_paths.html')


#Cell cycle score assignments in SK-N-BE(2)-C using Seurat.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Create the Seurat objects.
group1_seu <- CreateSeuratObject(b_rc[,names(cluster1)])
group2_seu <- CreateSeuratObject(b_rc[,names(cluster2)])
#Normalize and scale the Seurat matrices.
group1_seu <- NormalizeData(group1_seu)
group2_seu <- NormalizeData(group2_seu)
group1_seu <- ScaleData(group1_seu)
group2_seu <- ScaleData(group2_seu)
#Cell-cycle score assignment.
group1_seu <- CellCycleScoring(group1_seu, s.features = s.genes, 
                              g2m.features = g2m.genes, set.ident = TRUE)
group2_seu <- CellCycleScoring(group2_seu, s.features = s.genes, 
                               g2m.features = g2m.genes, set.ident = TRUE)



#Plot the proportion of cells in the two distinct phases. 
phases_1 <- group1_seu$Phase
phases_2 <- group2_seu$Phase
perc_1 <- table(phases_1)/sum(table(phases_1))*100
perc_2 <- table(phases_2)/sum(table(phases_2))*100

png(file='cell-cycle_phases_be2c.png',width = 700,height=700, res = 120)
par(mfrow=c(1,2))
g1 <-plot(phases_1, col=c('coral','chartreuse3','cyan3'),main='Group 1 cells phases',ylim=c(0,300),ylab='Number of cells'
         ,cex.lab=1.5,cex.main=1.5,font=2)
text(g1,table(phases_1)+20,labels =paste0(round(perc_1),rep('%',3)),font=2)
g2 <- plot(phases_2, col=c('coral','chartreuse3','cyan3'),main='Group 2 cells phases',ylim=c(0,300),ylab='Number of cells'
          ,cex.lab=1.5,cex.main=1.5,font=2)
text(g2,table(phases_2)+20,labels =paste0(round(perc_2),rep('%',3)),font=2)
dev.off()



#Map the mean TPMs on chromosomal bands.
tpms1 <- tpms[,names(cluster1)]
tpms2 <- tpms[,names(cluster2)]

load('annotation.rda')
msigdb <- msigdb

#Since it contains a lot of information, reduce what's inside to only chromosomal band
#gene information
chrom_bands <- msigdb[grep("c1_all",names(msigdb))] 

#bands is an empty vector that will contain, for each gene, the band in which it appeared
bands <- rep(NA,20440)
names(bands) <- rownames(tpms1)


#assign genes to the respective band only if it meets a standard format
for (i in 1:length(names(chrom_bands))){
  common <- intersect(names(bands),unlist(msigdb[i]))
  split <- unlist(strsplit(names(msigdb)[i],split=';_;'))[2]
  #bands[common] <- split
  if(length(unlist(strsplit(split,'p|q')))==2){
    bands[common] <- split  
  }
  else{
    bands[common]<- NA
  }
}

#The missing values represent those genes assigned to chromosomal bands without a further
#numberical identifier after the p or the q. we only keep those elements that map
#to bands written as chrNUMBERp|qNUMBER
bands <- bands[!is.na(bands)] 

#Code for the cycle that will assign the mean TPM expression of each band
band_expr_1 <- data.frame(expr=apply(tpms1,1,mean)[names(bands)],bands=bands,stringsAsFactors = F) 
band_expr_2 <- data.frame(expr=apply(tpms2,1,mean)[names(bands)],bands=bands,stringsAsFactors = F)
unique_band_expr_1 <- data.frame(meanExpr=NA,bands=unique(bands),chr=NA)
unique_band_expr_2 <- data.frame(meanExpr=NA,bands=unique(bands),chr=NA)
for (i in 1:length(unique(bands))){
  unique_band_expr_1[i,1] <-mean(band_expr_1[band_expr_1$bands==unique(bands)[i],1])
  unique_band_expr_1[i,3] <- unlist(strsplit(unique(bands)[i],split='chr|p|q'))[2]
  unique_band_expr_2[i,1] <-mean(band_expr_2[band_expr_2$bands==unique(bands)[i],1])
  unique_band_expr_2[i,3] <- unlist(strsplit(unique(bands)[i],split='chr|p|q'))[2]
}


#Sort the data of kelly. Sort the entries for chr 1-22 independently from x and y to treat
#them as numeric 
index <- c('x','y')
no_xy_1 <- unique_band_expr_1[!unique_band_expr_1$chr %in% index,]
xy_1 <- unique_band_expr_1[unique_band_expr_1$chr %in% index,]
no_xy_1 <- no_xy_1[order(as.numeric(no_xy_1$chr),no_xy_1$bands),]
xy_1 <- xy_1[order(xy_1$chr,xy_1$bands),]
final_1 <- as.data.frame(rbind(no_xy_1,xy_1),stringAsFactor=F)

#Sort the data of be2c Sort the entries for chr 1-22 independently from x and y to treat
#them as numeric 
index <- c('x','y')
no_xy_2 <- unique_band_expr_2[!unique_band_expr_2$chr %in% index,]
xy_2 <- unique_band_expr_2[unique_band_expr_2$chr %in% index,]
no_xy_2 <- no_xy_2[order(as.numeric(no_xy_2$chr),no_xy_2$bands),]
xy_2 <- xy_2[order(xy_2$chr,xy_2$bands),]
final_2 <- as.data.frame(rbind(no_xy_2,xy_2),stringAsFactor=F)

#Plot

#Produce the chromosome labels for the plot and the respective coordinates
labels <- c(paste0('chr',seq(1,22)),'chrx','chry')
lab_coord <- c()
count <-0
for(i in 1:length(unique(final_1$chr))){
  l <- nrow(final_1[final_1$chr==unique(final_1$chr)[i],])
  lab_coord[i] <- count + l%/%2
  count <- count + l
}

#Produce a bar chart of the results
random_col <- c(seq(1,length(unique(final_1$chr))-2,by=2),'x')
colors = ifelse(final_1$chr%in%random_col,'cornflowerblue','blue')
png("chromosome_distribution_group1.png",w=8000,h=1000,res=300)
par(las=2)
b <- barplot(final_1$meanExpr,names.arg = final_1$bands,las=2,width=rep(100,296),
             cex.names = 0.5,col=colors,ylab='Mean expressivity',
             main='Group 1 gene expression distribution',ylim=c(0,1000),cex.main=1.5)
text(b[lab_coord],800,labels = labels,cex =0.7)
dev.off()

random_col <- c(seq(1,length(unique(final_2$chr))-2,by=2),'x')
colors = ifelse(final_2$chr%in%random_col,'chartreuse3','chartreuse4')
png("chromosome_distribution_group2.png",w=8000,h=1000,res=300)
par(las=2)
b <- barplot(final_2$meanExpr,names.arg = final_2$bands,las=2,width=rep(100,296),
             cex.names = 0.5,col=colors,ylab='Mean expressivity',
             main='Group 2 gene expression distribution',ylim=c(0,1000),cex.main=1.5)
text(b[lab_coord],800,labels = labels,cex =0.7)
dev.off()


#Save the two objects assigning the clusters to the single-cells as two separate rda objects.
save(cluster1,'cluster1.rda')
save(cluster2,'cluster2.rda')




  