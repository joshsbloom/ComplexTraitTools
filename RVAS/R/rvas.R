# some dependencies 
library(BEDMatrix)
library(rrBLUP)
library(NAM)
library(tidyverse)
library(regress)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(VariantAnnotation)
library(GenomicFeatures)
library(vcfR)
library(SKAT)
library(Rfast)
library(seqinr)


# make our own txdb structure (motivated by annotation discrepancy found for GRX3 ORF and issues with dubious ORFs)
# pulled from http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz and gunzipped it (beware of changes)
# gff file for genome (contains gene coordinate information)
reference.gff='/data/yeast/rvas/ref/saccharomyces_cerevisiae.gff'

# I split the previously curated vcf (variant call format files with genotype information) into sub vcfs for each chromosome here 
vcf.prefix='/data/yeast/chr_rename-014.'


# output directory
out.dir = '/data/yeast/rvas/ref/'

# list of unique yeast chr names 
unique.chrs=paste0('chr', as.roman(1:16))

# get yeast genome sequence 
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3

# make our own txdb structure (motivated by annotation discrepancy found for GRX3 ORF and issues with dubious ORFs)
# pulled from http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz and gunzipped it (beware of changes)
txdb=makeTxDbFromGFF(file=reference.gff)




addProvean=function(coding.effects, meltedProveanScores) {
    goi=coding.effects$GENEID
    refAA=sapply(as.character(coding.effects$REFAA), s2c)
    varAA=sapply(as.character(coding.effects$VARAA), s2c)
    pdiff=mapply(function(x,y)which(x!=y), refAA,varAA)
    loi = mapply(function(x,y) {
                     x+y-1}, 
                     x=sapply(coding.effects$PROTEINLOC, function(x)x[1]),
                    y=pdiff)
    varoi=mapply(function(x,y) {
                     x[y]
                    },x=varAA,y=pdiff)
    refoi=mapply(function(x,y) {
                     x[y]
                    },x=refAA,y=pdiff)

    lookupchange=mapply(function(g,l,v){
               paste0(g,':',l, ':', v)
            #   meltedProveanScores$value[meltedProveanScores$gene==g & meltedProveanScores$ALT==v & meltedProveanScores$POS==l]
                    }, goi, loi, varoi)


    provScore=relist(meltedProveanScores$value[match(unlist(lookupchange), meltedProveanScores$lookup)],lookupchange)
    coding.effects$proveanEffects=provScore
    coding.effects$POSAA.delta=loi
    coding.effects$REFAA.delta=refoi
    coding.effects$VARAA.delta=varoi
    #coding.effects$gEdit= coding.effects$gEditInfo[,ns]
    #coding.effects$gEditSeqFwd=coding.effects$editedSequenceFwdL[,ns]
    return(coding.effects) #[,c('guideIndex','gEdit', 'gEditSeqFwd', 'GENEID',  'POSAA.delta', 'REFAA.delta', 'VARAA.delta', 'proveanEffects', 'CONSEQUENCE')])
}





# Global Variables
working.dir='/data/yeast/yeast_oligos/'

provean.score.file=paste0(working.dir, 'ref/provean_scores.RData') 
load(provean.score.file)

meltedProveanScores=lapply(provean_scores, 
                         function(x) {
                             y=reshape::melt(x)
                             names(y)[c(1,2)]=c('REF', 'ALT', 'score')
                             y$POS=rep(1:nrow(x), ncol(x))
                             return(y)
                         })
meltedProveanScores=data.table::rbindlist(meltedProveanScores, idcol='gene')
meltedProveanScores$lookup=paste0(meltedProveanScores$gene, ':', meltedProveanScores$POS, ':', meltedProveanScores$ALT)

# for more details about ORF annotations see http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_dna.README
orf_coding=read.fasta(paste0(working.dir, 'ref/orf_coding.fasta.gz'))
e_orfs=read.delim(paste0(working.dir, 'ref/E_orfs_column.txt'), stringsAsFactors=F, header=F, sep='\t')[,1]




# Jiayu you can skip this bit from here to *!! below ------------------------------------------------------------------------------------------
# I've split the vcf into one for each chromosome, I've shared with you chrI as an example 

# for each chromosome make 3 RDS objects

# gr_chr*.RDS a GRanges object with the coordinates, alternate allele, and a unique index containing pos_ref/alt for each variant
# pC_chr*.RDS a GRanges object containing the predicted coding consequences of variants falling within genes 
# Z_chr*.RDS a numeric matrix of allelic dosage for each variant 0 (homozygous ref), 1 (heterozygous ref/alt), 2 (homozygous alt/alt), 
#            extended to contain multi-allelic variants

# granges object to store coordinates and indices
for(chr in unique.chrs ) { 
    print(chr)
    vcf.in=read.vcfR(paste0(vcf.prefix, chr, '.vcf'))

    # vcfR function to calculate minor allele frequency (why is this f(x) so slow, let's just diy this)
    #mafdf=maf(vcf.in, element=2)

    # a bunch of annoying code to deal with multiallelic sites
    alt.var.list=list()
    gr.obj.list=list()
    multi.allele.index=list()
     #up to 10 alleles, will handle fewer ggracefully
    ref=getREF(vcf.in)
    for(i in 1:10) {
        alt.var.list[[as.character(i)]]=as.character(
                 sapply(sapply(getALT(vcf.in), function(x) strsplit(x, ',')) ,
                        function(x) x[i]))
        
        #length of alternate variant allele minus 1
        alt.var.lenm1=as.vector(sapply(alt.var.list[[i]], function(x) nchar(x)))-1 
        if(sum(is.na( alt.var.lenm1))==length(alt.var.lenm1)) { break;} 
        #hack to construct GRanges object     
        alt.var.lenm1.hack=alt.var.lenm1
        alt.var.lenm1.hack[is.na( alt.var.lenm1)]=0 
        # keep track of sites that don't contain a variant for a given allele 
        multi.allele.index[[as.character(i)]]=!is.na(alt.var.lenm1)
            
        gr.obj.list[[as.character(i)]] = GRanges(
                   seqnames=chr, 
                   ranges=IRanges(start=getPOS(vcf.in), 
                                  end=  getPOS(vcf.in) + alt.var.lenm1.hack),
                       strand=c(rep('+', length(alt.var.list[[i]]))) #PAMp)),rep('-', length(PAMm)) ) 
                    )

       alt.var.list[[as.character(i)]]=  alt.var.list[[as.character(i)]][!is.na(alt.var.lenm1)]
       gr.obj.list[[as.character(i)]]= gr.obj.list[[as.character(i)]][!is.na(alt.var.lenm1)]
       gr.obj.list[[as.character(i)]]$alt.var=alt.var.list[[as.character(i)]] 
       gr.obj.list[[as.character(i)]]$idx=paste0( as.character(seqnames(gr.obj.list[[as.character(i)]])), 
                                                 '_', 
                                                 start(gr.obj.list[[as.character(i)]]),
                                                 ':', 
                                                 ref[!is.na(alt.var.lenm1)], '/',
                                                 alt.var.list[[as.character(i)]] )
    }

    #GRanges object merging
    gr.obj=unlist(as(gr.obj.list, "GRangesList"))

    # gr.obj$maf=mafdf[,4]
    # gr.obj$idx=paste0( as.character(seqnames(gr.obj)), '_', start(gr.obj))
    genome(gr.obj)='sacCer3'

    
    
    gr.obj$alt.var[gr.obj$alt.var=="*"]=""

    #map variants to genes and predicted coding consequences
    pC=predictCoding(gr.obj,
              txdb, sacCer3,
              DNAStringSet(gr.obj$alt.var), 
              ignore.strand=TRUE)

    # extracts genotypes from vcf
    Z.all= t( extract.gt(vcf.in, as.numeric=F))

    # parse genotypes, see comment above for structure of Z
    Z.subsets=list()
    for(i in names(multi.allele.index)){
        Ztemp=Z.all[,multi.allele.index[[i]]]
        Ztemp[Ztemp=='0/0']='0'
        Ztemp[Ztemp==paste0(i,'/',i)]='2'
        cm1=paste0(i)
        #cm2=paste0('*/',i)
        Ztemp[grepl(cm1, Ztemp)]='1'
        Ztemp=matrix(as.numeric(Ztemp), ncol=ncol(Ztemp), dimnames=dimnames(Ztemp))
        Ztemp[is.na(Ztemp)]=0
        colnames(Ztemp)=gr.obj.list[[as.character(i)]]$idx
        Z.subsets[[as.character(i)]]=Ztemp
    }

    Z=do.call('cbind', Z.subsets)

    saveRDS(gr.obj, file=paste0(out.dir, 'gr_', chr, '.RDS'))
    saveRDS(Z, file=paste0(out.dir, 'Z_', chr, '.RDS'))
    saveRDS(pC, file=paste0(out.dir , 'pC_', chr, '.RDS'))
}

#-------------------------------------------------------------------------------------------------------------*!!
for(chr in unique.chrs ) { 
   print(chr)
   pC=readRDS(paste0(out.dir , 'pC_', chr, '.RDS'))
   pC=addProvean(pC, meltedProveanScores)
   saveRDS(pC, file=paste0(out.dir , 'pC_', chr, '.RDS'))
}


#for(chr in unique.chrs ) { 
#   print(chr)
#   pC=readRDS(paste0(out.dir , 'pC_', chr, '.RDS'))
#}




# Extract submatrices of genotypes ---------------------------------------------------------------

# extract all the submatrices per gene for variants that don't have synonymous coding consequences 
nonsyn.matrices=list()
annot.matrices=list()
# extract all the variants per chromosome with minor allele > 5%
for(chr in unique.chrs ) { 
   print(chr)

   #gr.obj=readRDS(gr.obj, file=paste0(out.dir, 'gr_', chr, '.RDS'))
   pC=readRDS(paste0(out.dir , 'pC_', chr, '.RDS'))
   Z=readRDS(paste0(out.dir, 'Z_', chr, '.RDS'))
   
   #pC$af.cnt=colSums(Z)
   pC.subset=pC[pC$CONSEQUENCE!='synonymous',]

   pC.subset$maxProv=(sapply(pC.subset$proveanEffects, function(x) {
                i=which.max(abs(x))  
                if(length(i)==0L) { return(NA) }
              # if()  { return(NA) }
               else { return (x[i]) }
              }))

   Zsub=Z[,match(pC.subset$idx, colnames(Z))]
   pC.subset$af=colSums(Zsub, na.rm=T)/(nrow(Zsub)*2)
   pC.subset$flipped=pC.subset$af>.5
   maf=pC.subset$af
   maf[maf>.5]=(1-maf[maf>.5])
   pC.subset$maf = maf 

 #code to flip alleles
#   flipme=which(pC.subset$flipped)
#   for(f in flipme) {
#    tmp=Z[,f]
#    tmp[Z[,f]==0]=2
#    tmp[Z[,f]==2]=0
#    Z[,f]=tmp
#   }

# convert to minor allele frequency 
#   maf=af
#   maf[maf>.5]=(1-maf[maf>.5])


   #split out allt the nonsynonymous and nonsense variants 
   spC=split(pC.subset, pC.subset$GENEID)

   for(gene in names(spC)){
        m = colnames(Z) %in% spC[[gene]]$idx
        nonsyn.matrices[[gene]]=Z[,m]
        annot.matrices[[gene]]=spC[[gene]]
    }
}
#----------------------------------------------------------------------------------------------------
saveRDS(annot.matrices, file='/home/jbloom/Downloads/peter2018_vcf/nonsyn_matrices_annot.RDS')
saveRDS(nonsyn.matrices, file='/home/jbloom/Downloads/peter2018_vcf/nonsyn_matrices.RDS')
#all variants
#gvec.corA=gvec.cor

gvecs=list()
for(gene in names(nonsyn.matrices)){
       Z= nonsyn.matrices[[gene]]
       if(!is.null(dim(Z)[2])) {
         #maf=sum(Z)/2
       #}
       #else {
         maf= colSums(Z)/(nrow(Z)*2)
       #}
       sqrtw=dbeta(maf,1,25)
       w=sqrtw^2
      gvecs[[gene]]= (Z)%*%diag(w) 
      colnames(gvecs[[gene]])=colnames(Z)
      #%*%rep(1,length(w))
       }
}
gvecr=do.call('cbind', gvecs)
#gvecr=do.call('cbind', nonsyn.matrices)

#gvec.cor2=gvec.cor
sgvr=scale(gvecr)
bi=which(is.na(sgvr[1,]))
sgvr=sgvr[,-bi]
gvecr=gvecr[,-bi]

#scaled genotypes / n markers 
gvec.cor=tcrossprod(sgvr)/(ncol(sgvr)-1)

gvec.corL=tcrossprod(sgvr[,!grepl('chrI_', colnames(sgvr))])/(sum(!grepl('chrI_', colnames(sgvr)))-1)


cor.mat=readRDS('~/cor.RDS')
library(RSpectra)
r10pc=eigs(cor.mat, k=5)



causal_prop=0.2



gene=names(nonsyn.matrices)[10]
Z= nonsyn.matrices[[gene]]
maf= colSums(Z)/(nrow(Z)*2)


rar=(which(maf<.01))
ncaus=round(length(rar)*causal_prop)

B=rep(0, length(maf))
B[sort(sample(rar, ncaus))]=-1

#get predicted strain effect
XB=Z%*%B

#flooring operation for strains with multiple variants
XB[XB<0]=-1

# phenotypic variance explained for a single covariate
# PVE=(beta^2*VarX)/((beta^2*VarX)+pheno.sd^2)
error.sd=.5

#adjust error.sd such that pve ranges from 0.01 - 0.5
#phenotypic variance explained, generalized to multiple covariates
pve=var(XB)/(var(XB)+error.sd^2)
print(pve)

simy=XB+rnorm(nrow(gvec.cor),mean=0, sd=error.sd)
#double-check paramaterization
summary(lm(simy~Z))

#PVE=
#sqrtw=dbeta(maf,1,25)
#w=sqrtw^2

#y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))

simy= sqrt(h2)*scale(XB) + rnorm(nrow(gvec.cor), mean=0, sd=sqrt((1-h2)/(h2*var(sqrt(h2)*scale(XB)))))

X=matrix(1, length(simy))

nmll=regress(simy~X, ~gvec.cor, verbose=T,pos=c(T,T))$llik
fmll=regress(simy~X, ~gvec.cor+geneKernel,verbose=T, pos=c(T,T,T))
calc.pval.LRtest(nmll, fmll$llik)

#simy=XB + rnorm(1011,sd=.25) ##,1) 
#simy[simy<(-5)]=-5
#hist(simy)
#scale((Z)%*%diag(w)%*%rep(1,length(w)))+rnorm(1011,sd=1)

X=matrix(1, length(simy))
X2=cbind(X, r10pc$vectors)
#specify null model
obj=SKAT_NULL_emmaX(simy ~ X, K=gvec.cor) #cor.mat) #gvec.cor) #cor.mat) #vec.cor) 
obj2=SKAT_NULL_emmaX(simy ~ X, K=gvec.corL) #cor.mat) #gvec.cor) #cor.mat) #vec.cor) 
obj3=SKAT_NULL_emmaX(simy ~ X, K=cor.mat) #gvec.cor) #cor.mat) #vec.cor) 

#_emmaX
##gvec.cor) #cor.mat) #gvec.cor) #cor.mat)

#null model
calc.pval.LRtest=function(null,full) { pchisq(-2*(null-full),1, lower.tail=FALSE) }

nmll=regress(simy~X, ~gvec.cor, verbose=T,pos=c(T,T))$llik


#full model
fll=list()
sll=list()
rll=list()
cll=list()

for(gene in names(nonsyn.matrices)[1:100]) { # names(which(-log10(vcps)>3)) ) { # names(nonsyn.matrices)[1:100]) {
     print(gene)   
     Z=nonsyn.matrices[[gene]]
     if(!is.null(dim(Z))){
            if(ncol(Z)>0) {
            #Z=nonsyn.matrices[[1]]
            maf= colSums(Z)/(nrow(Z)*2)
            #using skat weights
            sqrtw=dbeta(maf,1,25)
            w=sqrtw^2

            geneKernel=Z%*%diag(w)%*%t(Z)    
            #important we constrain the vc estimates to be positive ... 
            #this will cause near zero estimates to be upwardly biased but otherwise we sometimes get highly significant nonsense
        #    fll[[gene]]= regress( simy~X, ~gvec.cor+geneKernel, verbose=T, pos=c(T,T,T))$llik              

            sll[[gene]]=SKAT(Z,obj, is_dosage=T,  is_check_genotype=F)$p.value  #,  #method="optimal.adj")$p.value
            
          #  rll[[gene]]=SKAT(Z,obj2, is_dosage=T,  is_check_genotype=F)$p.value
            cll[[gene]]=SKAT(Z,obj3, is_dosage=T,  is_check_genotype=F)$p.value

            }
        }

}
vcps=calc.pval.LRtest(nmll, unlist(fll))
unlist(sll)

x11()
par(mfrow=c(4,1))
plot(-log10(unlist(cll)), main='skat with common variants')
abline(h=-log10(.05/6500))
plot(-log10(unlist(sll)), main='skat with all variants')
abline(h=-log10(.05/6500))

plot(-log10(unlist(sll)), main='skat loco')
abline(h=-log10(.05/6500))

plot(-log10(vcps), main='2 component mm')
abline(h=-log10(.05/6500))



ps=list()
    for(gene in names(nonsyn.matrices)[1:100]){
    print(gene)   
        Z=nonsyn.matrices[[gene]]
        if(!is.null(dim(Z))){
            if(ncol(Z)>0) {
            #Z=nonsyn.matrices[[1]]
            ps[[gene]]=SKAT(Z,obj, is_dosage=T,  is_check_genotype=F)$p.value  #,  #method="optimal.adj")$p.value
            }
        }
    }
  #  skat.results[[p]]=unlist(ps)


ups=unlist(ps)
exp.ps=(rank(ups, ties.method='first')+.5)/(length(ups)+1)
plot(-log10(exp.ps), -log10(ups), ylab='observed -log10(p)', 
     xlab='expected -log10(p)',
     main='with random effect of maf(Beta,1,25) weighted sum of variants'
)
abline(0,1)




     sort(-log10(runif(length(unlist(ps))))),
          sort(-log10(unlist(ps))) ,ylab='observed -log10(p)', xlab='expected -log10(p)', main='with random effect of common variants')



plot(
     sort(-log10(runif(length(unlist(ps))))),
          sort(-log10(unlist(ps))) ,ylab='observed -log10(p)', xlab='expected -log10(p)', main='with random effect of maf(Beta,1,25) weighted sum of variants')
abline(0,1)


plot(
     sort(-log10(runif(length(unlist(ps))))),
          sort(-log10(unlist(ps))) ,ylab='observed -log10(p)', xlab='expected -log10(p)', main='with random effect of maf(Beta,1,25) weighted sum of variants')
abline(0,1)


#with common variant relatedness correction
psr=ps

#psg = gvec.cor
psg=ps











nstack=do.call('cbind', nonsyn.matrices)
pvec=lapply(annot.matrices, function(x) data.frame(x[,c('idx', 'maxProv')]))
pvecs=do.call('rbind',pvec)
pvecM=pvecs[match(colnames(nstack), pvecs$idx),]

provMat=apply(nstack, 1, function(x) pvecM$maxProv)
provMat=provMat/2
bigprovCnt=apply(provMat, 2, function(x) {
            x[is.na(x)]=0
            sum(abs(x)>5) })      
medprovCnt=apply(provMat, 2, function(x) {
            x[is.na(x)]=0
            sum(abs(x)>3) })      
                 
                  
                  sum(is.na(x)))
hprovCnt=apply(provMat, 2, function(x) { y= abs(x[!is.na(x)]); return(sum(y>5)) })


#CYR1
gene='YJL005W'
gene='YMR037C'
provdist=rep(NA,length(nonsyn.matrices))
names(provdist)=names(nonsyn.matrices)
provCnt=rep(NA,length(nonsyn.matrices))


for(gene in names(nonsyn.matrices)) {
    am=annot.matrices[[gene]]
    mprov=am$maxProv
    mprov[is.na(mprov)]=-20
    zm=nonsyn.matrices[[gene]]
    if(length(mprov)>1 & !is.null(dim(zm)) ) {
    #hist(apply(mprov * zm ,1, min))
    minprovs=(apply(apply(zm, 1, function(x) x * mprov), 2, min))
    provdist[gene]=sum(minprovs<(-5))
    }else{provdist[gene]=NA}
}

mProvAll=unlist(lapply(annot.matrices, function(x) x$maxProv))
hist(mProvAll)



nvar=sapply(zm, ncol)



all.variants.cor=matrix(0, 1011, 1011)
all.variants.cnt=0

rare.variants.cor=matrix(0, 1011, 1011)
rare.variants.cnt=0

common.variants.cor=matrix(0,1011,1011)
common.variants.cnt=0

all.variants.maf=list()

all.variants=list()
common.variants=list()
rare.variants=list()

for(chr in unique.chrs ) { 
   print(chr)

   Z=readRDS(paste0(out.dir, 'Z_', chr, '.RDS'))
   Z[is.na(Z)]=0
   #calculate variant allele frequency 
   af=colSums(Z)/(nrow(Z)*2)
   which(af==1)

   #convert to minor allele frequency 
   maf=af
   maf[maf>.5]=(1-maf[maf>.5])

   acnt=colSums(Z>0)


   sZ=Rfast::standardise(Z)
   print('Done Standardizing')
   
   badvar=which(colSums(is.na(sZ))>0)
   if(length(badvar)>0){
    Z=Z[,-badvar]
    maf=maf[-badvar]
   # sZ=sZ[,-badvar]
   }

   all.variants.maf[[chr]]=maf
   common.variants[[chr]]=Z[,maf >= .05]
   rare.variants[[chr]]=Z[,maf < .05]

#   common.variants.cor=common.variants.cor+tcrossprod(sZ[,maf >= .05] ) #scale(common.variants[[chr]]))
#   common.variants.cnt=common.variants.cnt+ sum(maf >= .05)
#   print(str(common.variants.cor))

#   rare.variants.cor=rare.variants.cor+tcrossprod(sZ[,maf < .05] ) #scale(rare.variants[[chr]]))
#   rare.variants.cnt=rare.variants.cnt+ sum(maf < .05)
#   print(str(rare.variants.cor))
   
#   all.variants.cor=all.variants.cor+tcrossprod(sZ)
#   all.variants.cnt=all.variants.cnt+ncol(sZ)
#   print(str(all.variants.cor))

   
}


avc=all.variants.cor/(all.variants.cnt-1)
cvc=common.variants.cor/(common.variants.cnt-1)
rvc=rare.variants.cor/(rare.variants.cnt-1)



# calculate a kinship matrix from common variants ------------------------------------------------------
Zcommon=do.call('cbind', common.variants)
#Amat=A.mat(Zcommon-1, return.imputed=T)
#A=Amat$A
#g.subset.imp=Amat$imputed

##scale marker effects
g.s=Rfast::standardise(Zcommon)


cor.mat=tcrossprod(g.s) #[,1:82869])
cor.mat=cor.mat/(ncol(g.s)-1) 
#------------------------------------------------------------------------------------------------------

k=20
library(RSpectra)
r20pc=svds(t(g.s), k=k)

M=ncol(g.s)

ChiSq=matrix(0, ncol(g.s), k) 

for(sig in 1:20) {
    print(sig)
    for(j in 1:M) {
    sigma=r20pc$d[sig]
    ChiSq[j,sig]=(M/(sigma^2)) * ((g.s[,j]) %*% (r20pc$v[,sig]))^2
    }
}




plot(-log10(pchisq(rowSums(ChiSq[,1:20]), df=20,lower.tail=F)))

Zrare=do.call('cbind', rare.variants)
g.s=Rfast::standardise(Zrare)
rcor.mat=tcrossprod(g.s)/(ncol(g.s)-1)

#Zall=cbind(Zcommon, Zrare) #do.call('cbind', rare.variants)
#g.s=Rfast::standardise(Zall)
#acor.mat=tcrossprod(g.s)/(ncol(g.s)-1)

###load Joseph's prefiltered GWAS genotype data and pre-processed phenotypes ---------------------------
###https://www.nature.com/articles/s41586-018-0030-5
###read marker data (BED/BIM/BAM plink formatted)
### load the subset of common variants across the genome (maf>5%) 
x=BEDMatrix('/data/yeast/1002genomes/1011GWASMatrix/1011GWAS_matrix.bed')
g=as.matrix(x)
###load marker information
mnames=read.delim('/data/yeast/1002genomes/1011GWASMatrix/1011GWAS_matrix.bim', 
                  header=F, sep='\t', stringsAsFactors=F)
##
nacnt=apply(g, 2, function(x) sum(is.na(x)))

mnames=mnames[nacnt<500,]
g2=g[,nacnt<500]


Amat=A.mat(g2-1, return.imputed=T)
#A=Amat$A
g.subset.imp=Amat$imputed
##
###scale marker effects
g.s=scale(g.subset.imp)
##
###calculate pariwise correlation between strains
cor.mat=tcrossprod(g.s)/ncol(g.s)  #[,1:82869])
#cor.mat2=cor.mat2/(ncol(g.s)-1) #(82869-1)
###-------------------------------------------------------------------------------------------------------------

rcor.mat=readRDS('~/rcor.RDS')
cor.mat=readRDS('~/cor.RDS')
library(RSpectra)
r10pc=eigs(rcor.mat, k=10)

pheno=read.delim('/data/yeast/rvas/data/pheno_35Conditions_NormalizedByYPD.txt', sep='\t')
names(pheno)[1]='Strain'


full.names=rownames(nonsyn.matrices[[1]])

#example calculating the SKAT statistic for each trait for each gene 
skat.results=list()
for(p in colnames(pheno)[-1]) {
    print(p)
    #initialize y with NAs
    y=rep(NA, nrow(nonsyn.matrices[[1]]))
    names(y)=full.names
    # propagate phenotype information
    y[match(pheno$Strain, full.names)]=pheno[,p]

  
    #additional covariate matrix (here we're just specifying an intercept)         
    X=matrix(1, length(y))
   # X2=cbind(X, r10pc$vectors)
    #specify null model
    obj=SKAT_NULL_emmaX(y ~ X2, K=cor.mat) #gvec.cor) #cor.mat)
 
    ps=list()
    for(gene in names(nonsyn.matrices)){
        #print(gene)   
        Z=nonsyn.matrices[[gene]]
        if(!is.null(dim(Z))){
            if(ncol(Z)>0) {
            #Z=nonsyn.matrices[[1]]
            ps[[gene]]=SKAT(Z,obj, is_dosage=T, is_check_genotype=F, r.corr=1)$p.value
            }
        }
    }
    skat.results[[p]]=unlist(ps)
}

# with gvec.cor
skat.results.r0=skat.results
skat.results.r1=skat.results

# simulation sketch -------------------------------------------------------------------------

#pick a gene 
gene = names(nonsyn.matrices)[1]

#get genotype matrix for non-synonymous variants for that gene across the yeast population
Z=nonsyn.matrices[[gene]]

# an example simulated architecture
# calculate counts for each minor allele
l=apply(Z,2, sum, na.rm=T)

# get logical vector indicating which are singletons
singletons=l==1

# pick n sites as casual
n.causals=20
# pick n.causals at random from singletons to have an effect 
pickcausals=names(sample(which(singletons), n.causals, replace=F)           )

#simulate those causals as having an effect of 1, the rest having an effect of 0 (no effect) and some additional random noise (residual error)

y=(names(singletons) %in% pickcausals)%*%t(Z) + rnorm(nrow(cor.mat), sd=resid.error)
y=t(y) 
obj=SKAT_NULL_emmaX(y ~ X, K=cor.mat)
test=SKAT(Z,obj, is_dosage=T, is_check_genotype=F)
print(test$p.value)

#---------------------------------------------------------------------------------------------
                 














        # completely random phenotype --------------------------------------------
        resid.error=.1
       
        y=rnorm(nrow(cor.mat), sd=resid.error)
        
        #  obj=SKAT_Null_Model(y ~ X, out_type="C")
        obj=SKAT_NULL_emmaX(y ~ X, K=cor.mat)
        test=SKAT(Z,obj, is_dosage=T, is_check_genotype=F)
        print(test$p.value)

YPGALACTOSE

        
        #











        
    


    
