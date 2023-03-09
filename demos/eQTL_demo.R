library(MASS)

# 5720 transcripts (rows) vs 1012 segregants (columns) 
count.matrix=readRDS('/data/single_cell_eQTL/yeast/Bulk/data/countMatrix.RDS')
#peaks.per.gene=readRDS('/data/single_cell_eQTL/yeast/Bulk/data/ppg.RDS')
# batch factors 
gbatch.fact=readRDS('/data/single_cell_eQTL/yeast/Bulk/data/batch_factor.RDS')
# OD covariate 
OD.cov=readRDS('/data/single_cell_eQTL/yeast/Bulk/data/OD_cov.RDS')
# genotype data 
#1,012 segregants (rows) vs 11,530 markers (columns) reduced from ~50k markers 
# -1 is BY and +1 is RM
gdata=readRDS('/data/single_cell_eQTL/yeast/Bulk/data/gdata.RDS')


strain.counts=colSums(count.matrix)

# do some normalization for covariates --------------------------------------------------
#MED14 or RGR1
#YLR071C
#negative binomial regression
ymodel=glm.nb(count.matrix['YLR071C',]~offset(log(strain.counts))+gbatch.fact+scale(OD.cov))
#poisson regression
ymodel1=glm(count.matrix['YLR071C',]~offset(log(strain.counts))+gbatch.fact+scale(OD.cov), family='poisson')
#ols regression 
ymodel2=lm(log2(count.matrix['YLR071C',]+1)~log(strain.counts)+gbatch.fact+scale(OD.cov))
#-------------------------------------------------------------------------------------------------
# all these models give about the same results
# plot(residuals(ymodel), residuals(ymodel2))

# calculate LOD score -----------------------------------------------
# normalized trait data 
ycorrected=residuals(ymodel)
lm.null=lm(ycorrected~1)

lm.full=rep(NA, ncol(gdata))
for(i in 1:ncol(gdata)) {
    print(i)
    m=lm(ycorrected~gdata[,i])
    #plot(jitter(gdata[,i]), ycorrected)
    lm.full[i]=as.numeric(logLik(m))
    # keeping likelihood here, could use 
    # z or t statistic and p-value (wald test)
}
LR=(lm.full-logLik(lm.null))
pchisqLR=pchisq(2*LR, 1, lower.tail=F)
#-log10 p-value
plot(-log10(pchisqLR))

#convert LRS to LOD
LOD=(LR/log(10))[1,]
#----------------------------------------------------------------------
LODhack=(-1012*log(1-(crossprod(scale(ycorrected),scale(gdata))/(1012-1))^2))/(2*log(10))


# calculate permutation threshold (see Churchill 1994) for FPR control (or GWER or FWER)
# input is scaled phenotypes and genotypes
#do permutations faster
yperm=replicate(1000, sample(scale(ycorrected)))
r=crossprod(yperm,scale(gdata))/(1012-1)
LODperm=(-1012*log(1-r^2))/(2*log(10))

maxStats=apply(LODperm,1, max) 
plot(LOD)
abline(h=quantile(maxStats, .95), col='red')

# max QTL position is CIS
colnames(gdata)[which.max(LOD)]
abline(v=which.max(LOD), col='blue')

# rescan for next biggest QTL, regressing out the effect of the biggest QTL
rnext=crossprod(scale(residuals(lm(ycorrected~gdata[,"chrXII:277604_A/T"]))), scale(gdata))/(1012-1)
abline(v=which.max(rnext[1,]^2), col='blue')

# plot phenotype by genotype at biggest QTL
plot(gdata[,"chrXII:277604_A/T"], ycorrected)
