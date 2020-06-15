#In this script we get acquainted with the datasets under study performing
#basic analysis and inspection of the most expressed and variable genes.

library(Seurat)
library(ggplot2)
library(wordcloud)
library(fgsea)
library(msigdbr)

#Load the msigdb containing the partition of genes across chromosomal bands.
load('annotation.rda')

#Load the datasets and assign them to variables.
load('kelly-rawcounts.rda')
load('kelly-tpms.rda')
k_rawcounts <- rawcounts
k_tpms <- tpms

load('be2c-rawcounts.rda')
load('be2c-tpms.rda')
b_rawcounts <- rawcounts
b_tpms <- tpms


#Remove the features that are unexpressed in both conditions.
keep <- !(apply(k_tpms,1,sum)==0 & apply(b_tpms,1,sum)==0)
k_tpms <- k_tpms[keep,]
b_tpms <- b_tpms[keep,]

#Check the 100 most expressed genes in terms of cumulative and average TPMs.

#Kelly
cum_tpms_k <- apply(k_tpms,1,sum)
mean_tpm_k <- apply(k_tpms,1,mean)
View(sort(cum_tpms_k,decreasing=T)[1:100])
View(sort(mean_tpm_k,decreasing=T)[1:100])

#be2c
cum_tpms_b <- apply(b_tpms,1,sum)
mean_tpm_b <- apply(b_tpms,1,mean)
View(sort(cum_tpms_b,decreasing=T)[1:100])
View(sort(mean_tpm_b,decreasing=T)[1:100])


#Check the 100 most expressed genes in terms of cumulative and average raw counts
#Kelly
cum_rc_k <- apply(k_rawcounts,1,sum)
mean_rc_k <- apply(k_rawcounts,1,mean)
View(sort(cum_rc_k,decreasing=T)[1:100])
View(sort(mean_rc_k,decreasing=T)[1:100])

#be2c
cum_rc_b <- apply(b_rawcounts,1,sum)
mean_rc_b <- apply(b_rawcounts,1,mean)
View(sort(cum_rc_b,decreasing=T)[1:100])
View(sort(mean_rc_b,decreasing=T)[1:100])

#Is ACE2 expressed?
cum_rc_k['ACE2']
cum_rc_b['ACE2']
#It is unexpressed.

#Assess anticorrelation of MYC and MYCN expression and expression diff in other genes
#through the production of error bars.
bar_names <- c('MYC','RPS2','GAPDH','ACTB','RPL37A','MYCN','TEAD4',"HIST1H4C","MIR92B",'NPY',
               'MYBL2','PRDM8','HMGB2')
#Code for the confidence interval of the transformed counts
bar_k_mean <- apply(log2(k_tpms[bar_names,]+1),1,mean)
bar_k_sd <- apply(log2(k_tpms[bar_names,]+1),1,sd)

bar_b_mean <- apply(log2(b_tpms[bar_names,]+1),1,mean)
bar_b_sd <- apply(log2(b_tpms[bar_names,]+1),1,sd)

#A function to add arrows on the chart (for confidence interval representation).
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#Plot the error bars.
png('error_bars_k.png',width=700,heigh=700,  res=120)
k <- barplot(bar_k_mean,col='cyan3',ylim=c(-3,18),las=2,main = 'KELLY expression bars',
             ylab = 'Average lognorm expression')
error.bar(k,bar_k_mean, bar_k_sd)
dev.off()

png('error_bars_b.png',width=700,heigh=700, res=120)
b <- barplot(bar_b_mean,col='coral',ylim=c(-3,18),las=2,main='SK-N-BE(2)-C expression bars',
             ylab = 'Average lognorm expression')
error.bar(b,bar_b_mean, bar_b_sd)
dev.off()


#Plot the number of cells in which a gene is expressed vs the log mean expression of
#the gene.

k_bool <- ifelse(k_tpms>0,1,0)
k_no_of_cells <- apply(k_bool,1,sum)
sort(k_no_of_cells,decreasing=T)[1:20]

b_bool <- ifelse(b_tpms>0,1,0)
b_no_of_cells <- apply(b_bool,1,sum)
sort(b_no_of_cells,decreasing=T)[1:20]

#The gene expressed in more cells is RPL22 for both cell lines 
genes_scatter<- c('MYCN','ALK','TEAD4','PHOX2B','PHOX2A',
                  'RPS2','GAPDH','ACTB',
                  'LGALS1','VIM','S100A6', 'DLK1', 'IGF2','CHGA')

png(filename ='no_cells_vs_log10_tpms.png',width = 800, height = 700, res=100)
par(mfrow=c(1,2))
plot(x=b_no_of_cells,y=log10(mean_tpm_b),xlab='Cells with Gene',ylab = 'log10 mean expression in TPM',
     main='be2c scRNA seq', pch=16,col='coral',xlim = c(0,1100),cex.lab=1.5,cex.main=1.6)
textplot(b_no_of_cells[genes_scatter],log10(mean_tpm_b)[genes_scatter],
         words=genes_scatter,cex=1,font=c(rep(2,5),rep(3,4),rep(4,6)),new=F,ylim=c(-3,4),xlim=c(0,1000))
legend('topleft',legend= c('NBL genes','Houskeeping','DE genes'), text.font = c(2,3,4))

plot(x=k_no_of_cells,y=log10(mean_tpm_k),xlab='Cells with Gene',ylab = 'log10 mean expression in TPM',
     main='Kelly scRNA seq', pch=16,col='cyan3',xlim = c(0,1200),cex.lab=1.5,cex.main=1.6)
textplot(k_no_of_cells[genes_scatter],log10(mean_tpm_k)[genes_scatter],
         words=genes_scatter,cex=1,font=c(rep(2,5),rep(3,4),rep(4,6)),new=F,ylim=c(-3,4),xlim=c(0,1000))
legend('topleft',legend= c('NBL genes','Houskeeping','DE genes'), text.font = c(2,3,4))
par(mfrow=c(1,1))
dev.off()


#Plot the correlation of expression between the 2 cells 
png('correlation_expressions_kelly_be2c.png',width = 1000,height=1000, res=100)
plot(log10(mean_tpm_k),log10(mean_tpm_b),pch=16,col='cadetblue',
     xlab = 'log10(Mean TPMs Kelly)',ylab='log10(Mean TPMs Be2c)',main='Correlation of mean TPM counts',
     cex.lab=1.6,cex.main=2)
grid()
textplot(log10(mean_tpm_k[genes_scatter]),log10(mean_tpm_b[genes_scatter]),
         words=genes_scatter,cex=1.5,new=F,font=2)
cor <- cor.test(log10(mean_tpm_k),log10(mean_tpm_b),method='spearman')
text(0,4,labels = paste('p.value <','e-108'),cex=2)
text(0,4.5,labels = paste('r =',as.character(round(cor$estimate,digits = 2))),cex=2)
dev.off()


#Let us analyze what are the most variant genes in the pool
marker_genes <- c('GAPDH','ACTB','NPY','MYCN','ALK','PHOX2B','HAND2','GATA3','TP53','RPS2','RPS27',
                'RPL37A','RPL41')

#Generate two separate Seurat objects for KELLY and SK-N-BE(2)-C
k_seurat <- CreateSeuratObject(k_tpms)
b_seurat <- CreateSeuratObject(b_tpms)
#Find the most variable features thorugh Seurat.
k_seurat <- FindVariableFeatures(k_seurat,nfeatures=1000)
b_seurat <- FindVariableFeatures(b_seurat,nfeatures=1000)
#Plot the most variabl features.
png('variable_genes_seurat_k.png',width=800,height=600, res=110)
plot1 <- VariableFeaturePlot(k_seurat,cols=c('red','coral'))
plot2 <- LabelPoints(plot = plot1, points = marker_genes,repel=T,fontface='bold',size=5)
plot2+ggtitle('Variable features KELLY')+theme(plot.title = element_text(hjust = 0.5,size=22))
dev.off()
png('variable_genes_seurat_b.png',width=800,height=600, res=110)
plot1 <- VariableFeaturePlot(b_seurat,cols=c('red','coral'))
plot2 <- LabelPoints(plot = plot1, points = marker_genes, repel = TRUE,fontface='bold',size=5)
plot2+ggtitle('Variable features BE(2)C')+theme(plot.title = element_text(hjust = 0.5,size=22))
dev.off()

#Now we will map the mean TPMs on the respective chromosomal bands.

#Fetch the msigdb object saved under annotation.rda.
load('annotation.rda')
msigdb <- msigdb

#Since it contains a lot of information, reduce what is inside to the gene pertitions
#across chromosomal bands.
chrom_bands <- msigdb[grep("c1_all",names(msigdb))] 

#bands is an empty vector that will contain, for each gene, the band on which it lays.
bands <- rep(NA,20440)
names(bands) <- rownames(k_tpms)


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
#numerical identifier after the p or the q. we only keep those elements that map
#to bands written as chrNUMBERp|qNUMBER
bands <- bands[!is.na(bands)] 

#Code for the cycle that will assign the mean TPM expression to each band.
band_expr_k <- data.frame(expr=apply(k_tpms,1,mean)[names(bands)],bands=bands,stringsAsFactors = F) 
band_expr_b <- data.frame(expr=apply(b_tpms,1,mean)[names(bands)],bands=bands,stringsAsFactors = F)
unique_band_expr_k <- data.frame(meanExpr=NA,bands=unique(bands),chr=NA)
unique_band_expr_b <- data.frame(meanExpr=NA,bands=unique(bands),chr=NA)
for (i in 1:length(unique(bands))){
  unique_band_expr_k[i,1] <-mean(band_expr_k[band_expr_k$bands==unique(bands)[i],1])
  unique_band_expr_k[i,3] <- unlist(strsplit(unique(bands)[i],split='chr|p|q'))[2]
  unique_band_expr_b[i,1] <-mean(band_expr_b[band_expr_b$bands==unique(bands)[i],1])
  unique_band_expr_b[i,3] <- unlist(strsplit(unique(bands)[i],split='chr|p|q'))[2]
  }


#Sort the data of kelly. Sort the entries for chr 1-22 independently from x and y to treat
#them as numeric 
index <- c('x','y')
no_xy_k <- unique_band_expr_k[!unique_band_expr_k$chr %in% index,]
xy_k <- unique_band_expr_k[unique_band_expr_k$chr %in% index,]
no_xy_k <- no_xy_k[order(as.numeric(no_xy_k$chr),no_xy_k$bands),]
xy_k <- xy_k[order(xy_k$chr,xy_k$bands),]
final_k <- as.data.frame(rbind(no_xy_k,xy_k),stringAsFactor=F)

#Sort the data of be2c Sort the entries for chr 1-22 independently from x and y to treat
#them as numeric 
index <- c('x','y')
no_xy_b <- unique_band_expr_b[!unique_band_expr_b$chr %in% index,]
xy_b <- unique_band_expr_b[unique_band_expr_b$chr %in% index,]
no_xy_b <- no_xy_b[order(as.numeric(no_xy_b$chr),no_xy_b$bands),]
xy_b <- xy_b[order(xy_b$chr,xy_b$bands),]
final_b <- as.data.frame(rbind(no_xy_b,xy_b),stringAsFactor=F)

#Plot

#Produce the chromosome labels for the plot and the respective coordinates
labels <- c(paste0('chr',seq(1,22)),'chrx','chry')
lab_coord <- c()
count <-0
for(i in 1:length(unique(final_k$chr))){
  l <- nrow(final_k[final_k$chr==unique(final_k$chr)[i],])
  lab_coord[i] <- count + l%/%2
  count <- count + l
}

#Produce a bar chart of the results
random_col <- c(seq(1,length(unique(final_k$chr))-2,by=2),'x')
colors = ifelse(final_k$chr%in%random_col,'cyan3','blue')
png("chromosome_distribution_kelly.png",w=8000,h=1000,res=300)
par(las=2)
b <- barplot(final_k$meanExpr,names.arg = final_k$bands,las=2,width=rep(100,296),
        cex.names = 0.5,col=colors,ylab='Mean expressivity',
        main='KELLY gene expression distribution',ylim=c(0,1000),cex.main=1.5)
text(b[lab_coord],800,labels = labels,cex =0.7)
dev.off()

random_col <- c(seq(1,length(unique(final_b$chr))-2,by=2),'x')
colors = ifelse(final_b$chr%in%random_col,'coral','red')
png("chromosome_distribution_be2c.png",w=8000,h=1000,res=300)
par(las=2)
b <- barplot(final_b$meanExpr,names.arg = final_b$bands,las=2,width=rep(100,296),
        cex.names = 0.5,col=colors,ylab='Mean expressivity',
        main='SK-N-BE(2)-C gene expression distribution',ylim=c(0,1000),cex.main=1.5)
text(b[lab_coord],800,labels = labels,cex =0.7)
dev.off()




