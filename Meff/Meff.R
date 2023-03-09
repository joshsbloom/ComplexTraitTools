# Joshua Bloom 11/03/17
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

cor.mat=c

evals = eigen(cor.mat,symmetric=T)$values
M = length(evals)
L = M-1

# Equation 5
intevals=ifelse(evals>=1, 1, 0)
# modification for negative eigenvalues JB
nonintevals=c(evals-floor(evals))[evals>0]
Meff.li=sum(intevals) + sum(nonintevals)

print(Meff.li)


Bonf=1-(0.95^(1/Meff.li))

library('mvtnorm')

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
