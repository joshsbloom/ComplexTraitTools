library(BEDMatrix)
library(rrBLUP)

#load GWAS genotype data
setwd('/data/rrv2/1002genomes/1011GWASMatrix/')
x=BEDMatrix('1011GWAS_matrix.bed')
g=as.matrix(x)
#load marker information
mnames=read.delim('1011GWAS_matrix.bim', header=F, sep='\t', stringsAsFactors=F)

# function to calculate effective number of tests given LD matrix
#https://neurogenetics.qimrberghofer.edu.au/SNPSpD/Li2005.pdf
getMeff_Li_and_Ji=function(cor.mat) {
    evals = eigen(cor.mat,symmetric=T)$values
    M = length(evals)
    L = M-1
    # Equation 5 from Li 
    intevals=ifelse(evals>=1, 1, 0)
    # modification for negative eigenvalues JB
    # add up the contribution of fractional eigenvalues 
    nonintevals=c(evals-floor(evals))[evals>0]
    Meff.li=sum(intevals+nonintevals)
    print(Meff.li)
    return(Meff.li)
}

# restructure genotypes by chromosome
g.chr=list()
for(i in 1:17){
    print(i)
   g.chr[[as.character(i)]]=g[,which(mnames[,1]==i)]
}

# relatedness per chromosome, but called here for side effect of imputing missing data
A.chr=list()
for(i in 1:17){
    print(i)
    A.chr[[as.character(i)]]=A.mat(g.chr[[as.character(i)]]-1, return.imputed=T)
}

# LD per chromosome
LD.chr=list()
for(i in 1:17){
    print(i)
    chr=as.character(i)
    # faster than base R cor() but equivalent
    LD.chr[[chr]]=tcrossprod(t(scale(A.chr[[chr]]$imputed)))/(nrow(A.chr[[chr]]$imputed)-1)
     #A.mat(g.chr[[chr]]-1, return.imputed=T)
}

# calculate effective number of tests per chromosome 
Meff.li.chr=sapply(LD.chr, function(x) getMeff_Li_and_Ji(x))
Meff=sum(Meff.li.chr)

# controls FPR
Bonf=1-(0.95^(1/Meff))


# controls FDR
# for example, if there are n.markers (let's say 80,000)
i=1:n.markers
M=n.markers
#p=sorted(small to large p-values from association test)
# fdr threshold (let's say 5%)
q=.05
#find  the max i, where observed p-value is less than 
p < (q/Meff+(((i-1)/(M-1))*(q-(q/Meff))))
 
