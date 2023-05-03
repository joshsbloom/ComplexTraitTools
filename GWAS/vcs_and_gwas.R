library(BEDMatrix)
library(rrBLUP)
library(mvnpermute)

#library(NAM)
#library(milorGWAS)

source('/home/jbloom/Dropbox/code/ComplexTraitTools/GWAS/R/preprocessing.R')
source('/home/jbloom/Dropbox/code/ComplexTraitTools/GWAS/R/dofastLMM_mvnpermute.R')

pheno=readRDS('/data/yeast/Colonies/2022-09-04/AreaPhenotypes_GJ-SDM.rds')
    
path.to.gwas.data='/data/yeast/1002genomes/1011GWASMatrix/'

#get gwas variants from bed/bim/bam
GWAS.variants=getGWASvariants(path.to.gwas.data)

#sync up geno and pheno matrices
synced=syncGenoPheno(pheno, GWAS.variants)

pheno=synced$pheno
GWAS.variants=synced$GWAS.variants

#qc on gwas variants
GWAS.variants=filterGWASvariants(GWAS.variants, af.cutoff=.05, na.frq=.05)

#calc additive relatedness matrices 
GWAS.variants=calcA(GWAS.variants)

A=GWAS.variants$A

#tidy this stuff above 
#x=pheno%>%
#   pivot_longer(cols=starts_with(rep.phenos[1]), names_to='rep', values_to="area") %>%
#   select("Strain_Name", "Plate", "rep", "area") %>% unite("plate_rep", "Plate", "rep", remove=F)
    #replace NAs with 0s
#
svdA=fastLMMsvd(GWAS.variants$A)

gwas.results=list()

for(p in colnames(pheno)[-c(1,2)]){
   #an example of some preprocessing ----------------
   up= unlist(pheno[,p])
   upnorm=sqrt(up-min(up))
   #go ahead and remove plate effect 
   ry=residuals(lm(upnorm~pheno$Plate))
   names(ry)=rownames(GWAS.variants$g)
   #-------------------------------------------------

   gwas.results[[p]]=dofastLMM_mvnpermute(y=ry,X=NULL,GWAS.variants, svdA=svdA, nperm=500, REML=T) 
}



#calculate narrow-sense heritabilities -------------------------------------------
library(data.table)
#library(ggplot2)
library(regress)

    nh2.raw=list()
    nh2.sqrt=list()
    for(p in colnames(pheno)[-c(1,2)]) {
        y=as.vector(unlist(pheno[,p]))
        sqrty=sqrt(y-min(y))
        
        nh2.raw[[p]]=regress(y~pheno$Plate, ~A, pos=c(T,T), verbose=T)
        nh2.sqrt[[p]]=regress(sqrty~pheno$Plate, ~A, pos=c(T,T), verbose=T)
    }  

#results --------------------------------------------------------------------------
t(sapply(nh2.raw,function(x) x$sigma))/(sapply(nh2.raw, function(x) sum(x$sigma)))
t(sapply(nh2.sqrt,function(x) x$sigma))/(sapply(nh2.sqrt, function(x) sum(x$sigma)))



# repeated measures mixed model, we need to better standardize input here ----------
rep.phenos=c('GJonly_s.area', 'GJ-lp_s.area', 'GJ-lpdelta_s.area')

rrh2.raw=list()
rrh2.sqrt=list()

for( r in rep.phenos) {

    to.stack=grep(r, colnames(pheno))

    fixed.cols=c('Strain_Name', 'Plate')

    stacked=pheno[,c(fixed.cols, colnames(pheno)[to.stack[1]])]
    stacked$Plate=paste0(stacked$Plate,':', 1)

    colnames(stacked)[ncol(stacked)]='area'
    for(i in 2:length(to.stack) ){
        temp=pheno[,c(fixed.cols, colnames(pheno)[to.stack[i]])]
        colnames(temp)[ncol(temp)]='area'
        temp$Plate=paste0(temp$Plate,':', i)
        stacked=rbind(stacked, temp)
    }


    stacked$Strain_Name=paste0(stacked$Strain_Name, '_', stacked$Strain_Name)

    #let's just make this manually for now to avoid accidental scrambling
    Z=matrix(0, nrow(stacked),ncol=nrow(A))
    colnames(Z)=colnames(A)
    for(i in 1:nrow(stacked)){
        Z[i,stacked$Strain_Name[i]]=1
    }

    ZZt=Z%*%t(Z)
    ZAZt=Z%*%A%*%t(Z)
    sqrty=sqrt(stacked$area-min(stacked$area))
    

   # rrh2.raw[[r]]=regress(area~Plate, ~ZAZt+ZZt, pos=c(T,T,T), verbose=T,data=stacked)
   rrh2.sqrt[[r]]=regress(sqrty~Plate, ~ZAZt+ZZt, pos=c(T,T,T), verbose=T,data=stacked)

}

t(sapply(rrh2.raw,function(x) x$sigma))/(sapply(rrh2.raw, function(x) sum(x$sigma)))
t(sapply(rrh2.sqrt,function(x) x$sigma))/(sapply(rrh2.sqrt, function(x) sum(x$sigma)))

#-------------------------------------------------------------------------------------------


