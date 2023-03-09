library(BEDMatrix)
library(rrBLUP)
library(NAM)
library(milorGWAS)
library(data.table)
library(ggplot2)
library(mvnpermute)
#library(regress)
 
pheno=readRDS('/data/yeast/Colonies/2022-04-29/results/sRadiusMean_SNmcrREPs.RDS') #SNR1p5Metrics.RDS')
pheno=pheno[,-c(6,7)]
colnames(pheno)[4:ncol(pheno)]=paste0('ph',colnames(pheno)[4:ncol(pheno)])

#remove FY4 for now given we don't have its genotype (stupid, could fix)
pheno=pheno[-which(pheno[,1]=='FY4'),]

pheno=pheno[!is.na(pheno[,1]),]

#rename annoying columns
colnames(pheno)[3:ncol(pheno)]=paste0('ph', colnames(pheno)[3:ncol(pheno)])

mapping.matrix=model.matrix(~pheno$StrainName.JS.short-1)
colnames(mapping.matrix)=gsub("pheno\\$StrainName.JS.short", '', colnames(mapping.matrix))



path.to.gwas.data='/data/yeast/1002genomes/1011GWASMatrix/'
setwd(path.to.gwas.data)

    #read marker data (BED/BIM/BAM plink formatted)
    x=BEDMatrix('1011GWAS_matrix.bed')
    g=as.matrix(x)
    mnames=read.delim('1011GWAS_matrix.bim', header=F, sep='\t', stringsAsFactors=F)


    #remove blank and missing strains -----------------
    g.subset=g

    #match to phenotypes
    mind=pmatch(colnames(mapping.matrix), rownames(g.subset) )

 #   pheno=pheno[-which(is.na(mind)),]
#    mind=pmatch(colnames(mapping.matrix), rownames(g.subset) )

    g.subset=g.subset[mind,]
    #morejunk=which(rownames(pheno)=='NA')
    #-------------------------------------------------

    #calculate relatedness matrix ----------------------------
    Amat=A.mat(g.subset-1, return.imputed=T)
    A=Amat$A
    g.subset.imp=Amat$imputed

    g.s=scale(g.subset.imp)

    mnames=mnames[match(data.table::tstrsplit(colnames(g.s), '_')[[1]], mnames[,2]),]
    bm=which(is.na(mnames[,1]))
    g.subset.imp=g.subset.imp[,-bm]
    g.s=g.s[,-bm]
    mnames=mnames[-bm,]


    g.s=g.subset[,colnames(g.subset.imp)]

    #totalNA=sum(is.na(g.s))
    af=colSums(g.s)/(2*nrow(g.s))

    af=colSums(g.s, na.rm=T)/(2*apply(g.s,2,function(x) sum(!is.na(x)))) 
    af.cutoff=.05


    notna.cnt=apply(g.s, 2, function(x) sum(!is.na(x)))
    notna.frq=notna.cnt/nrow(g.s)
    notna.frq.cut=.8

    #gvar=apply(g.s,2,var, na.rm=T)
    bm=which(af<af.cutoff | notna.frq<notna.frq.cut)

    #final marker set
    g.s=g.s[,-bm]
    #corresponding marker info
    mnames=mnames[-bm,]

    #replace NAs with 0s
    g.s0=g.s
    g.s0[is.na(g.s0)]=0

      #calculate pariwise correlation between strains
    GRM=tcrossprod(g.s0)
    GRM=GRM/(ncol(g.s0))
    #---------------------------------------------------------

Z=mapping.matrix

ZZt=Z%*%t(Z)
ZAZt=Z%*%GRM%*%t(Z)
ZAAZt=Z%*%(A*A)%*%t(Z)
library(regress)
library(lme4)


vcms=list()

for(phn in names(pheno)[-c(1:3)]) {
    p=unlist(pheno[,phn])
    m0=lmer(p~1+(1|pheno$StrainName.JS.short))
    print(summary(m0))

    m1=regress(p~1,~ZZt, verbose=T, pos=c(T,T,T))
    m2=regress(p~1,~ZAZt+ZZt, verbose=T, pos=c(T,T,T,T))
    m3=regress(p~1,~ZAZt+ZAAZt+ZZt, verbose=T, pos=c(T,T,T,T,T))

    vcms[[phn]]$m0=m0
    vcms[[phn]]$m1=m1
    vcms[[phn]]$m2=m2
    vcms[[phn]]$m3=m3
}

rawVC0=sapply(vcms, function(x) x$m1$sigma)
normVC0=t(rawVC0)/colSums(rawVC0)

rawVC=sapply(vcms, function(x) x$m2$sigma)
normVC=t(rawVC)/colSums(rawVC)


rawVC2=sapply(vcms, function(x) x$m3$sigma)
normVC=t(rawVC2)/colSums(rawVC2)

