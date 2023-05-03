fastLMMsvd=function(A){
    USUt=svd(A)
    U=USUt$u
    Ut=t(USUt$u)
    return(list(USUt=USUt, U=U, Ut=Ut))
}

#pre-compute once for speedup if there are multiple traits 



#y named vector
#X fixed effects
#GWAS.variants list of genotypes (g) , annotations (mnames), and additive relatedness matrix (A)

dofastLMM_mvnpermute=function(y,X=NULL,GWAS.variants, svdA=NULL, nperm=100, REML=T) {

    #print(p)
    
    #ytemp=pheno[,p]
    #names(ytemp)=pheno[,1]
    #ytemp=ups
    
    A=GWAS.variants$A
    mnames=GWAS.variants$mnames

    #y=scale(ytemp)
    
    # update and sanity check NA handling 

    g=GWAS.variants$g
    g[is.na(g)]=0

    if(is.null(X)){
      X = matrix(1,length(y))
    }

    if(REML) {method='REML'} else {method='ML'}


    nmm=rrBLUP::mixed.solve(y,X=X,K=A, method=method)
  

    delta=nmm$Ve/nmm$Vu

    #function(y,g, delta, K) 
    #delta=sigmaE/sigmaG
    #eig=eigen(cor.mat)
    
    if(is.null(svdA)){
        svdA=fastLMMsvd(A)
    } 
    USUt=svdA$USUt
    U=svdA$U
    Ut=svdA$Ut


    Sdelt=USUt$d+delta

    n=length(y)
    Uty=Ut%*%cbind(y, mvnpermute::mvnpermute(y,matrix(1,length(y)),nmm$Vu*A+nmm$Ve*diag(length(y)),nperm))
    sds=solve(diag(Sdelt))


    llconst0=log(2*pi*nmm$Vu)
    llconst1=log(det(diag(Sdelt)))
    if(!is.finite(llconst1)) { llconst1=0 }
    llconst=n*llconst0+llconst1

    llconst2=(1/nmm$Vu)

    XBase=matrix(1,n)
    d=ncol(XBase)+1
    XComb=cbind(1,g)

    UtXall=Ut%*%cbind(1,g)

    #null model
    UtX=UtXall[,1]
    XX=XComb[,c(1)]
    tUS=t(UtX/Sdelt)
    B=solve(tUS %*% UtX) %*% (tUS %*%Uty)
    uubSq=(Uty-(UtX)%*%B)^2
    redll=-.5*(llconst+llconst2*colSums(uubSq/Sdelt))
    remlredll=redll+.5*(d*llconst0+log(det(t(XX)%*%XX))-log(det(tUS%*%UtX)))

    Bmat=matrix(NA,ncol(g),(nperm+1))
    ll=matrix(NA,ncol(g),(nperm+1))
    rll=ll
    #i=47718
    #system.time(
    for(i in 2:(ncol(g)+1) ) {
        if(i%%1000==0) {print(i) }
        UtX=UtXall[,c(1,i)] #cbind(1,UtXall[,i])

        tUS=t(UtX/Sdelt)
        #solve for betas
        B=solve(tUS %*% UtX) %*% (tUS %*%Uty)
        uubSq=(Uty-(UtX)%*%B)^2
        #genet  var exp vG=colSums(uubSq/Sdelt)/nrow(K)
        #calculate log likelihood
        #colSums
        #ll[i-1]=-.5*(llconst+llconst2*sum(uubSq/Sdelt))

        Bmat[i-1,]=B[2,]
        ll[i-1,]=-.5*(llconst+llconst2*colSums(uubSq/Sdelt))
        if(REML==TRUE){
            XX=XComb[,c(1,i)]
            rll[i-1,]=ll[i-1,]+.5*(d*llconst0+log(det(t(XX)%*%XX))-log(det(tUS%*%UtX)))
                #genet  var exp vG=colSums(uubSq/Sdelt)/nrow(K)-

        }

    }
    dff=data.frame(marker=colnames(g),
               Beta=Bmat[,1],
               ML.nlp=-log10(pchisq(-2*(redll[1]-ll[,1]),df=1, lower.tail=F)),
               REML.nlp=-log10(pchisq(-2*(remlredll[1]-rll[,1]),df=1, lower.tail=F))
               )
    attr(dff, "ML.perm")=apply(-log10(pchisq(-2*(redll[-1]-ll[,-1]),df=1, lower.tail=F)),2,max, na.rm=T)
    attr(dff,"REML.perm")=apply(-log10(pchisq(-2*(remlredll[-1]-rll[,-1]),df=1, lower.tail=F)),2,max, na.rm=T)

    dff$chrom=mnames$chrom
    dff$pos=mnames$pos
    dff$MLsig=dff$ML.nlp>quantile(attr(dff,  "ML.perm"), .95)
    dff$REMLsig=dff$REML.nlp>quantile(attr(dff,  "REML.perm"), .95)

    return(dff)
}

#let's make sure to resurrect some of this lme4 stuff 
#library(tidyverse)
#library(lme4)
#x=read_csv('~/Downloads/CFUs_per_condition.csv')
#xx=split(x, x$medium)
#lapply(xx, function(y) {
#    z=y[!is.na(y$"Colonies Counted"),]
#    tb=as.data.frame(VarCorr(lmer(log(CFUs)~1+(1|BruniPlate_ID),data=z) ))
#    print( (tb$vcov/sum(tb$vcov))[1])
#})

