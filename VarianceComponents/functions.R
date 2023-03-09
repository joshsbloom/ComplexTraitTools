

#calculate strain blups (see calcMM for variable definitions)
calc.BLUPS= function(G,Z,Vinv,y,X,B ){    ((G%*%t(Z)) %*% Vinv) %*%( y - X%*%B)     }



#these two functions are for using Eskin 2008 SVD trick for speedup when genotype matrix is 
# super optimized for one VC and fixed effects ~1000X speedup by precomputing eigen decomp
m.S=function (y, K = NULL, bounds = c(1e-09, 1e+09), theta=NULL, Q=NULL, X=NULL ) 
{
    n <- length(y)
    y <- matrix(y, n, 1)
    if(is.null(X) ) {  p <- 1    } else { p = ncol(X) }
    Z <- diag(n)
    m <- ncol(Z)
       
    omega <- crossprod(Q, y)
    omega.sq <- omega^2
    
    f.REML <- function(lambda, n.p, theta, omega.sq) {
        n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + lambda))
    }
    soln <- optimize(f.REML, interval = bounds, n - p, theta,  omega.sq)
    lambda.opt <- soln$minimum
    
    df <- n - p
    Vu.opt <- sum(omega.sq/(theta + lambda.opt))/df
    Ve.opt <- lambda.opt * Vu.opt
    VCs=c(Vu.opt, Ve.opt)
    return(VCs)
}

doEigenA_forMM=function(A ,X=NULL ) {
        n=nrow(A)
        if(is.null(X) ) {  X = matrix(rep(1, n), n, 1); p=1 } else {p=ncol(X) }
        XtX = crossprod(X, X)
        XtXinv = solve(XtX)
        S = diag(n) - tcrossprod(X %*% XtXinv, X)
        SHbS = S %*% A %*% S
        SHbS.system = eigen(SHbS, symmetric = TRUE)
        theta = SHbS.system$values[1:(n - p)] 
        Q = SHbS.system$vectors[, 1:(n - p)]
        return(list(theta=theta, Q=Q))
        }
#---------------------------------------------------------------------------------------------------------------



#simple phenotype simulator
# G is (n x M) matrix (n = individuals, M = markers)
# h2 is additive variance 
# nadditive is number of markers with effects (here randomly + or -)
# nsims is number of simulated phenotypes
simPhenotypes=function(G, h2, nadditive, nsims) {
    nsample  =nrow(G)
    nmarker  =ncol(G)

    simY=replicate(nsims, {
        #total number of additive loci
        a.eff  = rep(0,nsample)
        #markers with effects
        add.qtl.ind  = sort(sample(nmarker, nadditive))
        add.qtl.sign = sample(ifelse(runif(nadditive)>.5,1,-1),replace=T)

        for(i in 1:nadditive){ a.eff=a.eff+add.qtl.sign[i]*G[,add.qtl.ind[i]] }
        a.eff=scale(a.eff)
        g=sqrt(h2)*a.eff 
        y=as.vector(g+rnorm(nrow(G),mean=0,sd=sqrt((1-h2)/h2*var(g))))
        return(y)
    })
    #for tutorial code name the strains
    rownames(simY)=seq_along(simY[,1])
    return(simY)
}


