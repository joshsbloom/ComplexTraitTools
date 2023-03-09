library(BEDMatrix)
library(RcppML)
#mats=list()
#for(chr in c('I','II','III', 'IV', 'V', 'X')){
#    print(chr)
#    g=as.matrix(BEDMatrix(paste0('/data/elegans/CeNDR20210121_Plink/', chr, '.bed') ))
#    g=t(g)
#    g[is.na(g)]=0
#    #to simplify for now
#    g[g==2]=1
#    mats[[chr]]=g
#}

#g=rbind(mats[[1]], mats[[2]], mats[[3]], mats[[4]], mats[[5]], mats[[6]])
 

g=as.matrix(BEDMatrix(paste0('/data/yeast/noahbedbimbam/converted.bed') ))
g=t(g)
g[is.na(g)]=0

#g1=g[,1:1010]

subvec=c(1:ncol(g))[-1227]


GG=crossprod(g[,subvec])

GGp=crossprod(g)
#GGp2=crossprod(g[,1011:1276])
#
#test0=svd(GG)
#test1=svd(GGp)
#test2=svd(GGp2)



efitness=extraDistr::rinvchisq(length(subvec),12,1)
hist(efitness)
#efitness=rep(1, ncol(g))

eCount=rowSums(g[,subvec]%*%efitness)

#gCount=Rfast::rowsums(g)

eFreq=(.5*eCount)/sum(efitness)

efitness.norm=efitness/sum(efitness)


read.length=300
genome.size=12e6
depth.vec=c(5000,2500, 1000, 500,100,50) #,3,1)

yield=(depth.vec*genome.size)/read.length


t1=list()
for (depth in depth.vec) {
    t1[[as.character(depth)]]=rbinom(n=length(eFreq),size=depth, prob=eFreq)
}

t1=do.call('cbind',t1)

Gy=crossprod(g[,subvec], t1)
predictions=apply(Gy,2,function(x) {
       as.vector(RcppML::nnls(GG,matrix(x), fast_nnls=T)) #fast_nnls=F, cd_maxit=10000, cd_tol=1e-10))

})
predictions=apply(predictions, 2, function(x) x/sum(x))


par(mfrow=c(2,3))
for(i in 1:ncol(predictions)){
    plot(predictions[,i], efitness.norm,
         xlab='prediction',ylab='expected',
         main=paste0(depth.vec[i], 'X  ', yield[i]/1e6, ' million reads'),
        cex.lab=1.2, cex.axis=1.5, cex.main=1.5, cex.sub=1.5  )

    legend('topleft', legend= paste('r=', round(cor(predictions[,i], efitness.norm),2)), cex=1.5)
        #  abline(0,1)
}

plot(predictions





cor.test((unlist(predicted.list)), (unlist(expected.list))
lm(unlist(expected.list)~unlist(predicted.list))


predictions0=apply(predictions0, 2, function(x) x/sum(x))
plot(predictions0[,1], efitness0.norm)
for(i in 2:99) {
    points(predictions0[,i], efitness0.norm)
}


#------------------------------------------------------------------------------------------------------



#Elegans 
library(doMC)
registerDoMC(cores = 10)
#sub 100
set.seed(50)

expected.list=list()
predicted.list=list()
for( i in 1:100) {
    pick100=sort(sample.int(ncol(g),540))
    g.s=g[,pick100]
    GGs=crossprod(g.s)

    efitness0=runif(ncol(g.s)) # rep(1,ncol(g))
    eCount0=rowSums(g.s%*%efitness0)
    eFreq0=eCount0/sum(efitness0)
    efitness0.norm=efitness0/sum(efitness0)
    expected.list[[as.character(i)]]=efitness0.norm

    #depth.vec=c(500,100,50,30,10,5) #,3,1)
    yield=(depth.vec*genome.size)/read.length
    depth=50
    t0=list()
    print(i)
    t0[[as.character(i)]]=rbinom(n=length(eFreq),size=depth, prob=eFreq0)
    t0=do.call('cbind',t0)
    Gy0=crossprod(g.s, t0)
    predictions0=apply(Gy0,2,function(x) {
           as.vector(RcppML::nnls(GGs,matrix(x), fast_nnls=T)) #fast_nnls=F, cd_maxit=10000, cd_tol=1e-10))

    })
    predicted.list[[as.character(i)]]=predictions0/sum(predictions0)
}

plot((unlist(predicted.list)), (unlist(expected.list)), xlab='estimated', ylab='true', main='540 individuals')

cor.test((unlist(predicted.list)), (unlist(expected.list))
lm(unlist(expected.list)~unlist(predicted.list))


predictions0=apply(predictions0, 2, function(x) x/sum(x))
plot(predictions0[,1], efitness0.norm)
for(i in 2:99) {
    points(predictions0[,i], efitness0.norm)
}





#equalAF
#Acount=rowSums(g)
#Afreq=Acount/ncol(g)

#unequal AF
#efitness=rexp(ncol(g), rate=.6)




set.seed(100)
#simulate 10% missing 
efitness0=(runif(ncol(g))>.25)+0 # rep(1,ncol(g))

eCount0=rowSums(g%*%efitness0)
eFreq0=eCount0/sum(efitness0)
efitness0.norm=efitness0/sum(efitness0)



t0=list()
for (depth in depth.vec) {
    t0[[as.character(depth)]]=rbinom(n=length(eFreq),size=depth, prob=eFreq0)
}
t0=do.call('cbind',t0)
Gy0=crossprod(g, t0)
predictions0=apply(Gy0,2,function(x) {
       as.vector(RcppML::nnls(GGp,matrix(x), fast_nnls=T)) #fast_nnls=F, cd_maxit=10000, cd_tol=1e-10))

})
predictions0=apply(predictions0, 2, function(x) x/sum(x))


par(mfrow=c(2,3))
for(i in 1:ncol(predictions0)){
    plot(predictions0[,i], jitter(efitness0.norm),
         xlab='prediction',ylab='expected',
         main=paste(depth.vec[i], 'X  ', yield[i], 'reads'),
        cex.lab=1.2, cex.axis=1.5, cex.main=1.5, cex.sub=1.5  )
    #legend('topleft', legend= paste('r=', round(cor(predictions0[,i], efitness0.norm),2)), cex=1.5)
        #  abline(0,1)
}






hist(efitness.norm)
x11()
par(mfrow=c(2,3))

for(i in 1:ncol(predictions)){
    plot(predictions[,i], efitness.norm, xlab='prediction',ylab='expected',
         main=paste(depth.vec[i], 'X  ', yield[i], 'reads'),
        cex.lab=1.2, cex.axis=1.5, cex.main=1.5, cex.sub=1.5  )
    legend('topleft', legend= paste('r=', round(cor(predictions[,i], efitness.norm),2)), cex=2)
           abline(0,1)
}


test2=mdatools::mcrals.fcnnls(x, GGp)



test0=mdatools::mcrals.fcnnls(t1, g)

GGp=crossprod(g)
Gy=crossprod(g, t1)
test1=RcppML::nnls(GGp,Gy, fast_nnls=T)
plot(test1[,1], efitness)
cor(test1[,1], efitness)

test2=mdatools::mcrals.fcnnls(Gy, GGp)
plot(test2[1,], efitness)
cor(test2[1,], efitness)

test3=nnls::nnls(GGp,Gy)
plot(test3$x, efitness)
cor(test3$x, efitness)

test4=nnls::nnls(g,t1[,2])



simAlt=list()
  for (depth in depth.vec) {
        #depth=30
        #assuming PE150
        
        #yield=depth*genome.size/read.length #ead.length
        simAlt[[as.character(depth)]]=rbinom(n=length(eFreq),size=depth, prob=eFreq)
    }

    simCounts[[colnames(phenop)[p]]]=do.call('cbind', simAlt)
}

estFreqs=lapply(simCounts, function(s) {
    mcrals.fcnnls(s, g.s)
} )









pheno.short=list.files("/data/elegans/traits_with_validated_qtl/", full.names=F)
pheno.short=gsub('.tsv', '', pheno.short)
pheno.files=list.files("/data/elegans/traits_with_validated_qtl/", full.names=T)
phenos=lapply(pheno.files, readr::read_tsv)
names(phenos)=pheno.short
#blech
phenos=plyr::join_all(phenos, by='strain', type="full")

findme=match(paste0(phenos$strain,'_', phenos$strain), colnames(g))
phenos=phenos[which(!is.na(findme)),]

findme=match(paste0(phenos$strain,'_', phenos$strain), colnames(g))
g.s=g[,findme]

rownames(phenos)=phenos$strain
phenos=phenos[,-1]

phenos=(as.matrix(phenos))
phenop=apply(phenos,2,function(x) x - min(x,na.rm=T))
phenop[is.na(phenop)]=0

simCounts=list()
for(p in 1:ncol(phenop)){

    efitness=phenop[,p]
    eCount=rowSums(g.s%*%efitness)
    eFreq=eCount/sum(efitness)


    #1000x coverage
    depth.vec=c(500,100,50,30,10,5,3,1)

    simAlt=list()
    for (depth in depth.vec) {i

        #depth=30
        #assuming PE150
        read.length=300
        genome.size=120e6

        #yield=depth*genome.size/read.length #ead.length
        yield=(depth*genome.size)/read.length
        simAlt[[as.character(depth)]]=rbinom(n=length(eFreq),size=depth, prob=eFreq)
    }

    simCounts[[colnames(phenop)[p]]]=do.call('cbind', simAlt)
}

estFreqs=lapply(simCounts, function(s) {
    mcrals.fcnnls(s, g.s)
} )

saveRDS(estFreqs, file="/data/elegans/plots/nnls_estFreqs.RDS")

















for(i in 1:ncol(phenop)){
    print(cor(t(estFreqs[[i]]), phenop[,i])^2)
    pdf(paste0("/data/elegans/plots/", colnames(phenop)[i],'.pdf'), width=15,height=5)
    par(mfrow=c(1,8))
    for(j in 1:8) {
        plot(estFreqs[[i]][j,],phenop[,i],main=rownames(estFreqs[[i]])[j], 
             ylab='expected', xlab='observed', sub=paste('r^2=', round(cor(estFreqs[[i]][j,], phenop[,i])^2,2)))
    }
    dev.off()
    pairs(t(estFreqs[[i]]), phenop[,i])

}
cor(t(estFreqs[[2]]), phenop[,1])^2
cor(t(estFreqs[[1]]), phenop[,1])^2
cor(t(estFreqs[[1]]), phenop[,1])^2
cor(t(estFreqs[[1]]), phenop[,1])^2
cor(t(estFreqs[[1]]), phenop[,1])^2
cor(t(estFreqs[[1]]), phenop[,1])^2





















































library(doMC)
registerDoMC(cores = 10)
library(BEDMatrix)
g=as.matrix(BEDMatrix('/data/yeast/1002genomes/1011GWASMatrix/1011GWAS_matrix.bed'))
g=t(g)
g[is.na(g)]=0
#to simplify for now
g[g==2]=1

#equalAF
Acount=rowSums(g)
Afreq=Acount/ncol(g)

#unequal AF
efitness=rexp(1011, rate=.3)
eCount=rowSums(g%*%efitness)
eFreq=eCount/sum(efitness)


#1000x coverage
depth=2000
#assuming PE150
read.length=300
genome.size=1.2e7

#yield=depth*genome.size/read.length #ead.length
yield=(depth*genome.size)/read.length
    
t1=rbinom(n=length(eFreq),size=depth, prob=eFreq)

oFreq=t1/depth

#g=g[!duplicated(g),]
nnreg.result=cv.glmnet(g, oFreq, lower.limit=0, alpha=0, intercept=F, parallel=T, lambda=seq(.000001,.1, by=.001),nfolds=5) 






oFreq=t1/depth
#g[1,1]=NA

i=NMF::.fcnnls(g, as.matrix(t1)) #, pseudor) # oFreq)
cor.test(i$x[,1], efitness)

plot(h$x[,1], efitness)
plot(i$x[,1]/sum(i$x[,1]), efitness/sum(efitness),xlab='observed', ylab='expected', main='540 elegans wild isolates, 30X coverage pooled seq' )

cor.test(i$x[,1]/sum(i$x[,1]), efitness/sum(efitness),xlab='observed proportion', ylab='expected proportion', main='540 elegans wild isolates')

l=.fcnnls(as.matrix(g), as.matrix(t1), verbose=TRUE, pseudo=TRUE)

plot(l$coef[,1]/sum(l$coef[,1]), efitness/sum(efitness),xlab='observed', ylab='expected', main='540 elegans wild isolates, 30X coverage pooled seq' )



#g=g[!duplicated(g),]
nnreg.result=cv.glmnet(g3, oFreq, lower.limit=0, alpha=0, intercept=F, parallel=T, lambda=seq(.000001,.1, by=.001),nfolds=5) 

myRidgeRSS <- function(X,y,b, lambda){ 
                return( sum( (y - X%*%b)^2 ) + lambda * sum(b^2) ) 
              }
bfgsOptimConst2 = optim(myRidgeRSS, par = rep(1,540), X = g3, y = oFreq,
                        method = 'L-BFGS-B', lower = c(0,0), lambda = 1e-6)



XtX=crossprod(g3)
XtY=crossprod(g3,oFreq)

diagM=diag(1e-5,nrow=ncol(g3),ncol=ncol(g3))
q=qr.solve(XtX+diagM)%*%XtY
cor.test(q[,1], efitness) #q)

h=NMF::fcnnls(g, oFreq)
cor.test(h$x[,1], efitness)

plot(coef(nnreg.result)[-1], efitness)
cor.test(coef(nnreg.result)[-1], efitness)

#c(.01,.001,.0005,.0001))







test=cv.glmnet(g, Afreq, lower.limit=0, alpha=0, intercept=F, parallel=T)


library(vcfR)
x=read.vcfR('/data/yeast/1011Matrix-015.gvcf', nrows=5e4)
g=extract.gt(x, as.numeric=T)
library(Rfast)
rmg=rowMaxs(g, value=T)
g2=g[rmg<3,]
g2[is.na(g2)]=0

g2.dup=duplicated(g2)
g2=g2[which(!g2.dup),]

Acount=rowSums(g2)
Afreq=Acount/(2*ncol(g2))

#PE150 * 30e6 reads 
#(300*30e6)/1.2e7
efitness=rexp(1011, rate=.1)
eCount=rowSums(g2%*%efitness)
eFreq=eCount/(2*sum(efitness))

library(multiway)
XtX=crossprod(g2)
Xty=crossprod(g2,eFreq)
test0=fnnls(XtX, Xty)

test1=cv.glmnet(g2, eFreq, lower.limit=0, alpha=0, intercept=T, parallel=T, lambda=c(1,.01,.001,.0005,.0001)) #, lambda=c(0.0005,0.000559))

test2=cv.glmnet(g2, eFreq, lower.limit=0, alpha=0, intercept=F, parallel=T)
test3=cv.nnlasso(g2, eFreq, family="normal", tau=0)


test4=nnlasso(g2, eFreq, family="normal", tau=0, lambda=0.000559) #, path=F, lambda=0.000559, intercept=F)

#see algorithms-snmf.R from the NMF package

.fcnnls <- function(x, y, verbose=FALSE, pseudo=FALSE, eps=0){
	
	# check arguments
	if( any(dim(y) == 0L) ){
		stop("Empty target matrix 'y' [", paste(dim(y), collapse=' x '), "]")
	}
	if( any(dim(x) == 0L) ){
		stop("Empty regression variable matrix 'x' [", paste(dim(x), collapse=' x '), "]")
	}
	
	# map arguments
	C <- x
	A <- y
# NNLS using normal equations and the fast combinatorial strategy
	#
	# I/O: [K, Pset] = fcnnls(C, A);
	# K = fcnnls(C, A);
	#
	# C is the nObs x lVar coefficient matrix
	# A is the nObs x pRHS matrix of observations
	# K is the lVar x pRHS solution matrix
	# Pset is the lVar x pRHS passive set logical array
	#
	# M. H. Van Benthem and M. R. Keenan
	# Sandia National Laboratories
	#
	# Pset: set of passive sets, one for each column
	# Fset: set of column indices for solutions that have not yet converged
	# Hset: set of column indices for currently infeasible solutions
	# Jset: working set of column indices for currently optimal solutions
	#
	# Check the input arguments for consistency and initializeerror(nargchk(2,2,nargin))
	nObs = nrow(C); lVar = ncol(C);
	if ( nrow(A)!= nObs ) stop('C and A have imcompatible sizes')
	pRHS = ncol(A);
	W = matrix(0, lVar, pRHS);
	iter=0; maxiter=3*lVar;
	# Precompute parts of pseudoinverse
	#CtC = t(C)%*%C; CtA = t(C)%*%A;
	CtC = crossprod(C); CtA = crossprod(C,A);
	
	# Obtain the initial feasible solution and corresponding passive set
	K = .cssls(CtC, CtA, pseudo=pseudo);
	Pset = K > 0;
	K[!Pset] = 0;
	D = K;
	# which columns of Pset do not have all entries TRUE?
	Fset = which( colSums(Pset) != lVar );
	#V+# Active set algorithm for NNLS main loop
	oitr=0; # HKim
	while ( length(Fset)>0 ) {
		
		oitr=oitr+1; if ( verbose && oitr > 5 ) cat(sprintf("%d ",oitr));# HKim
		
		#Vc# Solve for the passive variables (uses subroutine below)				
		K[,Fset] = .cssls(CtC, CtA[,Fset, drop=FALSE], Pset[,Fset, drop=FALSE], pseudo=pseudo);
		
		# Find any infeasible solutions
		# subset Fset on the columns that have at least one negative entry
		Hset = Fset[ colSums(K[,Fset, drop=FALSE] < eps) > 0 ];
		#V+# Make infeasible solutions feasible (standard NNLS inner loop)
		if ( length(Hset)>0 ){
			nHset = length(Hset);
			alpha = matrix(0, lVar, nHset);
			while ( nHset>0  && (iter < maxiter) ){
				iter = iter + 1; 
				alpha[,1:nHset] = Inf;
				#Vc# Find indices of negative variables in passive set
				ij = which( Pset[,Hset, drop=FALSE] & (K[,Hset, drop=FALSE] < eps) , arr.ind=TRUE);			
				i = ij[,1]; j = ij[,2]
				if ( length(i)==0 ) break;			
				hIdx = (j - 1) * lVar + i; # convert array indices to indexes relative to a lVar x nHset matrix
				negIdx = (Hset[j] - 1) * lVar + i; # convert array indices to index relative to the matrix K (i.e. same row index but col index is stored in Hset)
				
				alpha[hIdx] = D[negIdx] / (D[negIdx] - K[negIdx]);				
				alpha.inf <- alpha[,1:nHset, drop=FALSE]
				minIdx = max.col(-t(alpha.inf)) # get the indce of the min of each row
				alphaMin = alpha.inf[minIdx + (0:(nHset-1) * lVar)]
				alpha[,1:nHset] = matrix(alphaMin, lVar, nHset, byrow=TRUE);
				D[,Hset] = D[,Hset, drop=FALSE] - alpha[,1:nHset, drop=FALSE] * (D[,Hset, drop=FALSE]-K[,Hset, drop=FALSE]);			
				idx2zero = (Hset - 1) * lVar + minIdx; # convert array indices to index relative to the matrix D
				D[idx2zero] = 0;
				Pset[idx2zero] = FALSE;
				K[, Hset] = .cssls(CtC, CtA[,Hset, drop=FALSE], Pset[,Hset, drop=FALSE], pseudo=pseudo);
				# which column of K have at least one negative entry?
				Hset = which( colSums(K < eps) > 0 );
				nHset = length(Hset);
			}
		}
		#V-#
		
		#Vc# Make sure the solution has converged
		#if iter == maxiter, error('Maximum number iterations exceeded'), end
		# Check solutions for optimality
		W[,Fset] = CtA[,Fset, drop=FALSE] - CtC %*% K[,Fset, drop=FALSE];
		# which columns have all entries non-positive
		Jset = which( colSums( (ifelse(!(Pset[,Fset, drop=FALSE]),1,0) * W[,Fset, drop=FALSE]) > eps ) == 0 );
		Fset = setdiff(Fset, Fset[Jset]);
		
		if ( length(Fset) > 0 ){				
			#Vc# For non-optimal solutions, add the appropriate variable to Pset						
			# get indice of the maximum in each column
			mxidx = max.col( t(ifelse(!Pset[,Fset, drop=FALSE],1,0) * W[,Fset, drop=FALSE]) )
			Pset[ (Fset - 1) * lVar + mxidx ] = TRUE;
			D[,Fset] = K[,Fset, drop=FALSE];
		}		
	}
	#V-#
	
	# return K and Pset
	list(coef=K, Pset=Pset)
}
# ****************************** Subroutine****************************
#library(corpcor)
.cssls <- function(CtC, CtA, Pset=NULL, pseudo=FALSE){
	
	# use provided function
	if( is.function(pseudo) ){
		pseudoinverse <- pseudo
		pseudo <- TRUE
	}
	
	# Solve the set of equations CtA = CtC*K for the variables in set Pset
	# using the fast combinatorial approach
	K = matrix(0, nrow(CtA), ncol(CtA));	
	if ( is.null(Pset) || length(Pset)==0 || all(Pset) ){		
		K <- (if( !pseudo ) solve(CtC) else pseudoinverse(CtC)) %*% CtA;
		# K = pseudoinverse(CtC) %*% CtA;
		#K=pinv(CtC)*CtA;
	}else{
		lVar = nrow(Pset); pRHS = ncol(Pset);
		codedPset = as.numeric(2.^(seq(lVar-1,0,-1)) %*% Pset);
		sortedPset = sort(codedPset)
		sortedEset = order(codedPset)
		breaks = diff(sortedPset);
		breakIdx = c(0, which(breaks > 0 ), pRHS);
		for( k in seq(1,length(breakIdx)-1) ){
			cols2solve = sortedEset[ seq(breakIdx[k]+1, breakIdx[k+1])];
			vars = Pset[,sortedEset[breakIdx[k]+1]];			
			K[vars,cols2solve] <- (if( !pseudo ) solve(CtC[vars,vars, drop=FALSE]) else pseudoinverse(CtC[vars,vars, drop=FALSE])) %*% CtA[vars,cols2solve, drop=FALSE];
			#K[vars,cols2solve] <-  pseudoinverse(CtC[vars,vars, drop=FALSE])) %*% CtA[vars,cols2solve, drop=FALSE];
			#TODO: check if this is the right way or needs to be reversed
			#K(vars,cols2solve) = pinv(CtC(vars,vars))*CtA(vars,cols2solve);
		}
	}
	
	# return K
	K
}


