#In the present script, we will compare the average normalized counts of 
#our single-cell-sequenced data and the ones collected from the Harenza's
#bulk RNA-seq dataset of neuroblastoma cell lines.

library(Seurat)
library(wordcloud)
library(msigdbr)
library(fgsea)

#Load the data.
load('kelly-tpms.rda')
k_tpms <- tpms
load('be2c-tpms.rda')
b_tpms <- tpms
load('harenza-vstnorm.rda')
bulk_vst <- vstnorm
load('harenza-tpms.rda')
bulk_tpms <- tpms
load('harenza-fpkms.rda')
bulk_fpkm <- fpkms
load('harenza-rawcounts.rda')
bulk_rawcounts <- rawcounts

#Check the dimension of the TPM-normalized bulk counts.
dim(bulk_tpms)

#Setup the names to assign to the single-cell mean TPMs columns.
b <- 'sc-SK-N-BE-2--C'
k <- 'sc-KELLY'

#Reduce the single-cell TPMs to average counts.
sckelly <- apply(k_tpms,1,mean)
scbe2c <- apply(b_tpms, 1,mean)

#Generate a unique dataset containing as columns the bulk TPMs together with the average TPMs for KELLY and
#SK-N-BE(2)-C cells.
total_dataset <- bulk_tpms
nrow(bulk_tpms)

keep <- intersect(rownames(bulk_tpms), rownames(k_tpms))
total_dataset <- total_dataset[keep,]

k_tpms_new <- sckelly[keep]
b_tpms_new <- scbe2c[keep]
total_dataset <- cbind(total_dataset,k_tpms_new)
total_dataset <- cbind(total_dataset,b_tpms_new)
colnames(total_dataset)[c(41,42)] <- c(k,b)

dim(total_dataset)

#Now we generate a Spearman correlation table associating the cells sequenced with our
#method with the ones derived from the bulk-RNA-seq dataset.

#Calculate the correlation coefficients for both KELLY and SK-N-BE(2)-C cells.
cor_kelly <-apply(total_dataset,2,function(x){
  sc <- cbind(round(cor.test(x,total_dataset[,'sc-KELLY'],method='spearman')$estimate,digits=2),
              round(cor.test(x,total_dataset[,'sc-KELLY'])$p.value,digits=2))
  return(sc)
})
cor_kelly<- t(cor_kelly)[-c(41,42),]
colnames(cor_kelly) <- c('SCC','p.value')
cor_kelly <- cor_kelly[order(cor_kelly[,1],decreasing=T),]


cor_be2c <-apply(total_dataset,2,function(x){
  sc <- cbind(round(cor.test(x,total_dataset[,"sc-SK-N-BE-2--C"],method='spearman')$estimate,digits=2),
              round(cor.test(x,total_dataset[,"sc-SK-N-BE-2--C"])$p.value,digits=2))
  return(sc)
})
cor_be2c<- t(cor_be2c)[-c(41,42),]
colnames(cor_be2c) <- c('SCC','p.value')
cor_be2c <- cor_be2c[order(cor_be2c[,1],decreasing=T),]

#Export html table of the results.
library(htmlTable)
k_html <- htmlTable(cor_kelly,cgroup = c('CORRELATION scKELLY'))
b_html <- htmlTable(cor_be2c,cgroup = c('CORRELATION scBE2C'))
write(k_html,file='kelly_table.html')
write(b_html,file='be2c_table.html')


#Perform hierarchical clustering on the joint dataset of KELLY and SK-N-BE(2)-C.
png('HM_bulk_ss.png',width=800,height=850)
cormat <- cor(total_dataset,method='spearman')
distmat<-as.dist(1-cormat)
hcc<-hclust(distmat,method='average')
hcd<-as.dendrogram(hcc,horiz=T)
plot(hcd,main='Hierarchical clustering on the joint sc and bulk dataset',
     cex.main=2,font=2, ylab = 'Correlation distance', cex.lab=1.5)
dev.off()


#Notice that the dendrogram individuate a separate single cell cluster

#Perform the same dendrogram drafting considering single-cell counts separately from 
#each other. 
png('Sep_hierarhical_k.png',width=1000,height=1000)
cormat_kelly <- cor(total_dataset[,-42],method='spearman')
distmat_kelly<-as.dist(1-cormat_kelly)
hcc_kelly<-hclust(distmat_kelly,method='average')
hcd_kelly<-as.dendrogram(hcc_kelly,horiz=T)
plot(hcd_kelly,main='Hierarchical clustering of bulk RNA dataset and single-cell KELLY',
     cex.main=2,font=2)
dev.off()

png('Sep_hierarhical_b.png',width=1000,height=1000)
cormat_kelly <- cor(total_dataset[,-41],method='spearman')
distmat_kelly<-as.dist(1-cormat_kelly)
hcc_kelly<-hclust(distmat_kelly,method='average')
hcd_kelly<-as.dendrogram(hcc_kelly,horiz=T)
plot(hcd_kelly,main='Hierarchical clustering of bulk RNA dataset and single-cell SK-N-BE(2)-C',
       cex.main=2,font=2)
dev.off()
#Seems like the single-cell bias is not removable.


