# START HERE to skip image processing ------------------------------------------
library(tidyverse)
library(Rfast)
library(regress)
library(data.table)
library(foreach)
library(doMC)

ComplexTraitToolsPath='/home/jbloom/Dropbox/code/ComplexTraitTools/'
base.path=c('/data/yeast/Leslie/02272023/')
cross.file='/data/yeast/Leslie/cross.RDS'
sample_sheet_to_seg_name.file='/data/yeast/Leslie/newsample_sheet_leslie_200_1.csv'

#change this 
qtl.fx.file=paste0(ComplexTraitToolsPath, 'plateQTL/plateQTL_fx.R')
source(qtl.fx.file)

results.df=readRDS(paste(base.path, 'ResultsTable.RDS', sep =''))

#filter out plate2 
lsegs=results.df %>% filter(file!='LS_0002_2023-02-17_14-16-26.JPG') %>% 
    filter(file!='LS_0001_2023-02-17_14-15-06.JPG') %>% filter(layout=='LY20.csv')
 lseg.split=split(lsegs, lsegs$condition)

#calculate broad-sense heritability
H2=sapply(lseg.split, function(x) { 
         #  P=model.matrix(~x$file-1)   
           Z=model.matrix(~x$StrainName-1)
           ZtZ=Z%*%t(Z)
           rr=regress(x$s.radius.mean~1, ~ZtZ,verbose=T)
           return(rr$sigma[1]/sum(rr$sigma))
                      })

barplot(H2, ylab='Broad-sense heritability', main='H^2')

#calculate strain means
avg.pheno=sapply(lseg.split, function(x) {
                sapply(split(x$s.radius.mean, x$StrainName), mean, na.rm=T) })

#read in genotype data 
cross=readRDS(cross.file)
cross$pheno$id= gsub('.bam.txt', '', cross$pheno$id ) 

# read in key to map sample names in sample sheet and cross object to sample names in layout files
seg.key=read.delim(sample_sheet_to_seg_name.file, header=F, sep=',')
#seg.key=read.delim('/data/yeast/Leslie/sample_sheet_leslie_200_1.csv', header=F, sep=',')

#get Leslie's strain names into the cross object
cross$pheno$id2=seg.key[,5][match(cross$pheno$id, seg.key[,1])]

#reorder the matrix of average strain values given their order in the genotype cross object
avg.pheno.m=avg.pheno[match(cross$pheno$id2, rownames(avg.pheno)),]

geno=extractGenotype.argmax(cross)

#get chr boundaries as marker indices
chr.bounds=cumsum(c(0,rle(tstrsplit(colnames(geno), '_')[[1]])$lengths))

#scale the genotypes
geno.s=scale(geno)

#scale the phenotypes 
avg.pheno.m.s=scale(avg.pheno.m)

#calculate pairwise correlation between the segregants
A=tcrossprod(geno.s)/(ncol(geno.s)-1)

#calculate narrow-sense heritability (optional but a good 
h2=apply(avg.pheno.m.s, 2, function(x) { 
         #  P=model.matrix(~x$file-1)   
           #Z=model.matrix(~x$StrainName-1)
           #ZtZ=Z%*%t(Z)
           rr=regress(x~1, ~A ,verbose=T)
           print(summary(rr))
           return(rr$sigma[1]/sum(rr$sigma))
                      })
barplot(h2, ylab='narrow-sense heritability', main='H^2')


#quick scan for QTL, input phenotype matrix and genotype need to be scaled with no NAs
L=fasterLOD(nrow(avg.pheno.m), avg.pheno.m.s, geno.s)

#plot results
x11()
par(mfrow=c(7,1))
for(i in 1:7){
    plot(L[i,], ylab='LOD', xlab='marker index', main=rownames(L)[i])
    abline(v=chr.bounds, lty=2, col='lightblue')
    #approx threshold for sig
    abline(h=3.5, col='red')
}

#normalize trait values by the control plates , this is just one of the ways to do that 
YNB.traits=scale(residuals(lm(avg.pheno.m[,2:3]~avg.pheno.m[,'YNB'])))
YPD.traits=scale(residuals(lm(avg.pheno.m[,c(4,5,7)]~avg.pheno.m[,'YPD'])))
YPD=scale(avg.pheno.m[,'YPD'])
YNB=scale(avg.pheno.m[,'YNB'])
#Ldelta_YNB=fasterLOD(nrow(avg.pheno.m), YNB.traits, geno.s)
#Ldelta_YPD=fasterLOD(nrow(avg.pheno.m), YPD.traits, geno.s)

#combine results together
modTraits=cbind(YNB, YPD, YNB.traits, YPD.traits)
colnames(modTraits)[c(1,2)]=c('YNB', 'YPD')

# retrun mapping for QTL on normalized data (1D-scan)
fL=fasterLOD(nrow(avg.pheno.m), scale(modTraits), geno.s)

#permutation based threshold
nperm=1000
ntraits=ncol(modTraits)
pL=matrix(NA,ntraits,nperm)
for(i in 1:nperm){
    print(i)
    pL[,i]=apply(fasterLOD(nrow(avg.pheno.m), scale(modTraits[sample(1:nrow(modTraits)),]), geno.s), 1, max)
}
    
apply(pL, 1, quantile, .95)
#[1] 3.652 3.650 3.612 3.341 2.553 3.282 3.614
 
#remake plots
par(mfrow=c(7,1))
for(i in 1:7){
    plot(fL[i,], ylab='LOD', xlab='marker index', main=rownames(fL)[i])
    abline(v=chr.bounds, lty=2, col='lightblue')
    abline(h=quantile(pL[i,],.95), col='red')
}

#slightly more sophisticated approach with an FDR based cutoff and forward-selection for QTL, 10,000 permutations per regression step
detectedQTL=list()
for(i in 1:ncol(modTraits)) {
    print(colnames(modTraits)[i])
    detectedQTL[[colnames(modTraits)[i]]]=doTraitFDR(scale(modTraits[,i]), geno.s, geno.s, FDR_thresh=.05, nperm=1e4, doLODdrop=T)
}
dQ=rbindlist(detectedQTL, idcol='trait')
write.table(dQ[order(dQ$trait,dQ$LOD, decreasing=T),], '/data/yeast/Leslie/sigQTL.tsv', quote=F, row.names=F, sep='\t')
write.csv(dQ[order(dQ$trait,dQ$LOD, decreasing=T),], '/data/yeast/Leslie/sigQTL.csv', quote=F, row.names=F)




