# general mixed model code for variance component analysis
# should link R to multithreaded blas library (like openblas or MKL) for best peformance

library(Matrix)
library(regress)

basedir='/home/jbloom/Dropbox/code/ComplexTraitTools/VarianceComponents/'
source(paste0(basedir,'functions.R'))
source(paste0(basedir,'calcMM.R'))

# pull 1000 BYxRM genotype data---------------------------------------------------------------------------------------
load(url("http://genomics-pubs.princeton.edu/YeastCross_BYxRM/data/cross.Rdata"))
#extracts and recodes genotype as -1 (BY) and +1 RM ... input genotypes could also be coded as 0,1,2 for F2 cross
#extractGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$data }))*2)-3 }

#extract genotype data  where G is (n x M) matrix (n = individuals, M = markers)
#!!!could use your own genotypes here
G  =  extractGenotype(cross)

# scale genotypes
G.s = scale(G)
#calculate strain relatedness as (G %*% t(G))/n 
A = tcrossprod(G.s)/ncol(G)

#for example simulate h2 = 60% and 20 markers with effects, constant effect size with random sign, could be whatever

#!!!could use your own phenotypes here
set.seed(5)
simY=simPhenotypes(G, h2=.6, nadditive=20, nsims=50) 


example_trait=simY[,7]
# 1) --- using random regression---------------------
    #take one of the simulations, for example the 7th one
    r=regress(example_trait~1, ~A)

    #additive variance
    r$sigma[1]
    # 0.5611
    #standard error
    sqrt(diag(r$sigma.cov)[1])
    #0.07371 


# 2) --- using tutorial code --------------------------
    #note vector needs to be named for this function
    r=calcMM(y=example_trait, B=list(A=A))
    #additive variance
    r$Var[1]
    #0.5611
    sqrt(diag(r$invI))[1]
    #0.07371
#----------------------------------------------------


# 3) using Eskin 2008 trick if genotype matrix is the same for different phenotypes ---
    eigA= doEigenA_forMM(A)

    #initialize matrix for output
    vA=matrix(0,ncol(simY),2)
    colnames(vA)=c('additive_var', 'error_var')

    #per phenotype
    for(i in 1:ncol(simY)){
        vA[i,]=m.S(simY[,i],theta=eigA$theta, Q=eigA$Q)
    }

    print(vA[7,1])
    #additive_var 
    #  0.5611 
# ------------------------------------------------------




# example F2 genotypes and phenotypes
# transpose to make it individuals (rows) X markers (cols)
G=t(read.delim(paste0(basedir, 'genos.txt'), header=F))
simY=read.delim(paste0(basedir, 'phenos.txt'), header=F)
G.s=scale(G)
A=tcrossprod(G.s)/ncol(G)
#1)
r=regress(simY[,1]~1, ~A)

    #additive variance
    r$sigma[1]
    #18.48 

    #standard error
    sqrt(diag(r$sigma.cov)[1])
    #1.093e-18 
    
namedY=simY[,1]
names(namedY)=seq_along(namedY)

#2)
#nr usually a bit more stable , but all of these algorithms fail with low N
r=calcMM(y=namedY, B=list(A=A), alg='nr')

#3)
    eigA= doEigenA_forMM(A)

    #initialize matrix for output
    vA=matrix(0,ncol(simY),2)
    colnames(vA)=c('additive_var', 'error_var')

    #per phenotype
    for(i in 1:ncol(simY)){
        vA[i,]=m.S(simY[,i],theta=eigA$theta, Q=eigA$Q)
    }
    print(vA[1,1])

