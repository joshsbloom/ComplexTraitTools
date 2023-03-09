# Joshua Bloom 11/03/17

lodToPval =  function(x){  pchisq(x*(2*log(10)),df=1,lower.tail=FALSE)/2 }
PvalTolod =  function(x){  qchisq(x*2,df=1,lower.tail=FALSE)/(2*log(10)) } 

RtoLOD=function(n.pheno, r) {(-n.pheno*log(1-r^2))/(2*log(10)) }

# p-value threshold using the procedure of Li and Ji. 
# https://www.nature.com/hdy/journal/v95/n3/full/6800717a.html
# see equation 5

# needs marker correlation matrix (or back this out of genetic map)

# imagine we have 5 independent variables and a bunch of correlated variables (50 each))
corevecs=replicate(5, rnorm(1000))

# 1000 'individuals' by 250 'markers'
redunvecs=list()
for(i in 1:5) {
    redunvecs[[as.character(i)]] =replicate(50, { corevecs[,i]+rnorm(1000,sd=.1) } )
}
redunvecs=do.call('cbind', redunvecs)

# example with completely uncorrelated 
#redunvecs=replicate(250, rnorm(1000))

library(fields)
# visualize
cor.mat=(cor(redunvecs))
image.plot(cor.mat)
load('/home/jbloom/Desktop/snps.Rda')
g=snps[,-c(1:4)]
g=t(g)
g=apply(g, 2, scale)
g=t(apply(g,1,function(x) {x[is.na(x)]=0; return(x);}))


cor.mat=tcrossprod(t(g))
cor.mat=cor.mat/(nrow(g)-1)
geno=read.delim('/media/jbloom/d1/misc/Meff/FileS1.csv', header=T, sep=',')
g=geno[,-c(1:4)]
cor.mat=cor(t(g))
sg=split(g,geno$chr)
cg=lapply(sg, function(x) cor(t(x)))

getMeff_Li_and_Ji=function(cor.mat) {
    evals = eigen(cor.mat,symmetric=T)$values
    M = length(evals)
    L = M-1
    # Equation 5 from Li 
    intevals=ifelse(evals>=1, 1, 0)
    # modification for negative eigenvalues JB
    nonintevals=c(evals-floor(evals)) #[evals>0]
    Meff.li=sum(intevals+nonintevals)
    print(Meff.li)
    return(Meff.li)
}

Meff.li=getMeff_Li_and_Ji(cor.mat)

Meff.li.chr=sapply(cg, function(x) getMeff_Li_and_Ji(x))
sum(Meff.li.chr)

Bonf=1-(0.95^(1/Meff.li))
LOD.li=PvalTolod(Bonf)

rmvn=rmvnorm(1000, sigma=cor.mat)
z2r=tanh(rmvn)
sLOD=RtoLOD(359,z2r)


plot(-log10(2*pnorm(-abs(rmvn[9,]))))
z2p=2*pnorm(-abs(rmvn[9,]))
plot(PvalTolod(z2p))

# using method from
#http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000456
# and 
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3325408/

# we want to find the effective number of tests 

# see equation 3 from Eskin paper
library(mvtnorm)
cinv <- function (x) qnorm(x/2, lower = FALSE) 
toMin=function(n.eff, alpha, cor.mat) {
    alpha.prime=alpha/n.eff
    phi.inv=cinv(alpha.prime)
    l=rep(-phi.inv, nrow(cor.mat))
    u=rep(phi.inv, nrow(cor.mat))
    x=pmvnorm(lower=-phi.inv, upper=phi.inv,mean=rep(0, nrow(cor.mat)), corr=cor.mat)
    fwer=1-x
    print(paste(n.eff, fwer))
    return( (fwer-alpha) )
}

intv=c(6, 200)
mvn.meff=lapply(cg, function(cmm) { uniroot(toMin, intv, alpha=.05, cor.mat=cmm ) } )
Meff.mvn=sum(sapply(mvn.meff, function(x) x$root))
Bonf.mvn=1-(0.95^(1/Meff.mvn))

LOD.mvn=PvalTolod(Bonf.mvn)

x=t(g)
n=359

y=replicate(1000, rnorm(359))
p=cor(y,x)
tt=p/(sqrt((1-p^2)/(n-2)))
pvals=2*pt(-abs(tt), n-2)
apply(pvals,2, m




quantile(minp, .05)
Bonf.mvn

pt(tt[1], n-2, lower=FALSE)



neff=366

mvn_optim=function(neff, cor.mat){
        alpha=.05/neff
        phi=abs(qnorm(alpha/2))

        mvn.cd=pmvnorm(mean=rep(0, nrow(cor.mat)), sigma=cor.mat, 
                lower=rep(-phi, nrow(cor.mat)),
                 upper=rep(phi, nrow(cor.mat)))

        tail.area=(1-mvn.cd)
        print(neff)
        print(tail.area)
        return(tail.area-.05)
} 
mvn.neff=uniroot(mvn_optim, lower=5, upper=300, cor.mat=cor.mat)
