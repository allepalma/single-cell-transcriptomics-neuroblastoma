library(corto)
library(Seurat)
library(biomaRt)

#Load the datsets.
load('kelly-rawcounts.rda')
k_rc <- rawcounts
load('be2c-rawcounts.rda')
b_rc <- rawcounts
#Load the datsets.
load('kelly-tpms.rda')
k_rc <- tpms
load('be2c-tpms.rda')
b_rc <- tpms

#Remove genes that are unexpressed in both cells.
keep <- !(apply(k_rc,1,sum)==0 & apply(b_rc,1,sum)==0)
k_rc <- k_rc[keep,]
b_rc <- b_rc[keep,]

#For differential gene expression analysis, we will treat the datasets separately.
k_seu <- CreateSeuratObject(k_rc)
b_seu <- CreateSeuratObject(b_rc) 
k_seu <- NormalizeData(k_seu)
b_seu <- NormalizeData(b_seu)

#Try MRA analysis on the pre cell cycle data
load('target_NBLrnaseq_regulon.rda')
TARGET <- regulon
load('kocak_NBL_regulon.rda')
kocak <- regulon

#Launch master regulators analysis with the Seurat-normalized data
b <- GetAssayData(b_seu,slot='data')
k <- GetAssayData(k_seu,slot='data')


#Crossroad

#To launch the mra run these commands. 
#mra1_before_cc <- mra(b,k,regulon=TARGET, nthreads = 2)
#mra2_before_cc <- mra(b,k,regulon=kocak, nthreads = 2)

#To load the precomputed mra objects run the following
#load(mra.rda)

png('MRA_target.png',width=2000,height=4000,res=200)
mraplot(mra1_before_cc,mrs=10)
dev.off()

png('MRA_kokack.png',width=2000,height=4000,res=200)
mraplot(mra2_before_cc,mrs=10)
dev.off()

#Calculate the correlation between the normalized enrichment scores obtained by running the
#master regulator analysis using the two different regulons.
keep <- intersect(names(mra1_before_cc$nes),names(mra2_before_cc$nes))
nes1_before_cc <- mra1_before_cc$nes[keep]
nes2_before_cc <- mra2_before_cc$nes[keep]
cor_before_cc <- cor.test(nes1_before_cc ,nes2_before_cc,method='spearman')
cor.stats_before_cc <- c(cor_before_cc$estimate,cor_before_cc$p.value)

png('corr.png',width=1000,height=1000, res=130)
plot(nes1_before_cc,nes2_before_cc,pch=20,cex=1.5,col='darkseagreen3',main='Correlation of
     enrichement scores',cex.lab=1.4,cex.main=2,xlab='NES TARGET',
     ylab='NES KOCAK',cex.font=2)
abline(lm(nes2_before_cc~nes1_before_cc),col='black')
text(-20,20,paste0('rho = ',round(cor.stats_before_cc[1],digits=2)),font=2,cex=1.4)
text(-20,15,paste0('p.value = 1.15e-42'),font=2,cex=1.4)
dev.off()

View(nes1_before_cc)


#We can use biomaRt for the annotation of the transription factors we retrieved.     
library(biomaRt)        
# look at top 10 databases      
head(biomaRt::listMarts(host = "www.ensembl.org"), 10)      


# 1) select a mart and data set        
mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",         
                            mart    = useMart("ENSEMBL_MART_ENSEMBL",       
                                              host = "www.ensembl.org"))       

# 2) run a biomart query using the getBM() function        
# and specify the attributes and filter arguments      
geneSet <- names(nes1_before_cc)        

resultTable <- biomaRt::getBM(attributes = c("external_gene_name","description","chromosome_name","band"),       
                              filters    = "hgnc_symbol",       
                              values     = geneSet,         
                              mart       = mart)        

View(resultTable)









