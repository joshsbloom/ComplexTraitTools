#1000 BYxRM data
load(url("http://genomics-pubs.princeton.edu/YeastCross_BYxRM/data/pheno_raw.Rdata"))
#pheno_raw
load(url("http://genomics-pubs.princeton.edu/YeastCross_BYxRM/data/cross.Rdata"))
#4000 BYxRM data
#load('~/Dropbox/cross.RData')
library(regress)
library(rrBLUP)
library(qtl)
library(fields)
#library(gputools)
#source('~/Dropbox/calcMM.R')

# some useful functions
extractGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$data }))*2)-3 }
extractGenotype.argmax=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$argmax }))*2)-3 }
countStrainsPerTrait = function(pheno) {apply(pheno, 2, function(x){sum(!is.na(x))})}
get.LOD.by.COR = function(n.pheno, pheno, gdata) {
    # Lynch and Walsh p. 454
    return( (-n.pheno*log(1-cor(pheno, gdata, use='pairwise.complete.obs')$coefficients^2))/(2*log(10)) ) }
calc.pval.LRtest=function(null,full) { pchisq(-2*(null-full),1, lower.tail=FALSE) }
calc.BLUPS= function(G,Z,Vinv,y,X,B ){    G%*%t(Z)%*%Vinv%*%(y- X%*%B)     }
extractVarCompResults = function(r) {list(sigma=r$sigma, sigma.cov=r$sigma.cov, llik=r$llik) }

extractScaledPhenotype=function(impcross,scaleVar=FALSE){apply(impcross$pheno, 2, scale, scale=scaleVar)}


#extract genotype data
gdata     = extractGenotype(cross)

gcoord=as.numeric(do.call('rbind', strsplit(colnames(gdata), '_'))[,1])
# chr name
chr=as.character(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])

#for 1000 BYxRM-------------------------------------------------------------
chr=gsub('chr0', '', chr)
chr=gsub('chr', '', chr)
chr=as.numeric(chr)
# chr pos (ignore for now)
# chr name again, fixed for so/extrrting
chr.lab=do.call('rbind', strsplit(colnames(gdata), '_'))[,2]
#--------------------------------------------------------------------------

#for 4000 BYxRM-------------------------------------------------------------
#chr = gsub("chr", "", chr)
#chr.lab = match(chr, as.character(as.roman(1:16)))

pos=as.numeric(do.call('rbind', strsplit(colnames(gdata), '_'))[,3])

pdata  = extractScaledPhenotype(cross, TRUE)

A = A.mat(gdata, shrink=FALSE)/2

nsample=1008
nmarker=11623
nadditive=30
H2=.75
h2=.75
a.eff  = rep(0,nsample)
add.qtl.ind  = sort(sample(nmarker, nadditive))
add.qtl.sign = sample(c(1,-1),nadditive,replace=T)
for(i in 1:nadditive){ a.eff=a.eff+add.qtl.sign[i]*gdata[,add.qtl.ind[i]] } 
a.eff=scale(a.eff)
g=sqrt(h2)*a.eff 
y2=g+rnorm(nrow(gdata),mean=0,sd=sqrt((1-H2)/H2*var(g)))
y2=scale(y2)

h2=.75
H2=.75
nadditive=30 #
add.qtl.ind.2= c(add.qtl.ind[1:15], sort(sample(nmarker, nadditive-15)))
add.qtl.sign.2 =c(add.qtl.sign[1:15], sample(c(1,-1),nadditive-15,replace=T))
a.eff.2  = rep(0,nsample)
for(i in 1:nadditive){ a.eff.2=a.eff.2+add.qtl.sign.2[i]*gdata[,add.qtl.ind.2[i]] } 
a.eff.2=scale(a.eff.2)
g=sqrt(h2)*a.eff.2
y3=g+rnorm(nrow(gdata),mean=0,sd=sqrt((1-H2)/H2*var(g)))
y3=scale(y3)

summary(regress(y2~1, ~A, verbose=T))
summary(regress(y3~1, ~A, verbose=T))
Z0=matrix(0, 1008,1008)

Z1=rbind(diag(1008),Z0)
Z2=rbind(Z0, diag(1008))

g1=Z1%*%A%*%t(Z1)
g21=Z2%*%A%*%t(Z1)
g12=Z1%*%A%*%t(Z2)
g22=Z2%*%A%*%t(Z2)
E1=Z1%*%t(Z1)
E2=Z2%*%t(Z2)
E12=Z1%*%t(Z2)
E21=Z2%*%t(Z1)

r=regress(c(y2,y3)~1,~g1+g12+g21+g22+E1+E2, identity=F, verbose=T)

y4=c(y2,y3)
um=unique(sort(c(add.qtl.ind, add.qtl.ind.2)))

tfact=as.factor(c(rep('A', 1008), rep('B',1008)))
sfact=as.factor(c(seq(1,1008), seq(1,1008)))
g2=rbind(gdata[,um],gdata[,um])

g2=(g2+1)/2

Einc=cbind(c(rep(1,1008), rep(0,1008)), c( rep(0,1008), rep(1,1008)))

rlm=lm(y4~tfact*g2-1)
XX=model.matrix(rlm)





y=c(y1,y2)
#summary(regress(y1~1, ~A, verbose=T))
#summary(regress(y2~1, ~A, verbose=T))

y1=y

#y3=c(y,y2)
#y=c(y1,y2)

Z0=matrix(0, 1008,1008)
Z1=rbind(diag(1008),Z0)
Z2=rbind(Z0, diag(1008))
#summary(
s1=regress(y1~1,~A, verbose=T, pos=c(T,T))
s2=regress(y2~1,~A, verbose=T, pos=c(T,T))

g1=Z1%*%A%*%t(Z1)
g21=Z2%*%A%*%t(Z1)
g12=Z1%*%A%*%t(Z2)
g22=Z2%*%A%*%t(Z2)
E1=Z1%*%t(Z1)
E2=Z2%*%t(Z2)
E12=Z1%*%t(Z2)
E21=Z2%*%t(Z1)

yr=c(pdata[,2], pdata[,43])

r=regress(yr~1,~g1+g12+g21+g22+E1+E2, identity=F, verbose=T)


yr=c(pdata[,1], pdata[,2])

Gmat=matrix(NA,46,46)
sigmas=list()
for(i in 1:45) {
    for(j in (i+1):46){
    print((  paste(colnames(pdata)[i],colnames(pdata)[j],sep='_')))
    yr=c(pdata[,i], pdata[,j])
    yrr=regress(yr~1,~g1+g12+g21+g22+E1+E2, identity=F, verbose=T)
    simgas[[paste(i,j,sep='_')]]=yrr$sigma
    Gmat[i,j]=as.numeric(yrr$sigma[2])
    }
    #yrr$sigma[2]
}
save(Gmat, file='~/Dropbox/1000BYxRM_gmat.RData')
save(sigmas, file='~/Dropbox/1000BYxRM_bivar.RData')



for(i in 1:45) {
    for(j in (i+1):46){
 #   print((  paste(colnames(pdata)[i],colnames(pdata)[j],sep='_')))
 #   yr=c(pdata[,i], pdata[,j])
 #   yrr=regress(yr~1,~g1+g12+g21+g22+E1+E2, identity=F, verbose=T)
 #   simgas[[paste(i,j,sep='_')]]=yrr$sigma
    Gmat[j,i]=Gmat[i,j]
                     #as.numeric(yrr$sigma[2])

    }
    #yrr$sigma[2]
}


pcor=cor(pdata, pdata, use='pairwise.complete.obs')
hp=hclust(dist(pcor))
plot(hp)

image.plot(pcor[hp$order,hp$order])

Gmat2=Gmat
Gmat2[Gmat>1]=.9
#Gmat2[lower.tri(Gmat2)]=Gmat2[upper.tri(Gmat2)]
diag(Gmat2)=1
colnames(Gmat2)=colnames(pdata)
rownames(Gmat2)=colnames(pdata)

hp2=hclust(dist(Gmat2))
image.plot(Gmat2[hp$order,hp$order])

corrplot(Gmat2, method='ellipse', type='upper', diag=FALSE, 
         addCoef.col=NULL, order='FPC', tl.col='black')

library(corrplot)
scor = cor(pdata, use='pairwise.complete.obs', method='pearson')
corrplot(scor, method='ellipse', type='upper', 
         diag=FALSE, addCoef.col=NULL, order='FPC', tl.col='black')
#scor = cor(pdata, use='pairwise.complete.obs', method='spearman')
corrplot(Gmat2, method='shade', type='upper', diag=FALSE, addCoef.col=FALSE, order='hclust', tl.col='black')


load('/media/kserver/kruglyak/raid1/home/jbloom/1000BYxRM/QTL/peakindex.bin')

g1=Z1%*%A%*%t(Z1)
g21=Z2%*%A%*%t(Z1)
g12=Z1%*%A%*%t(Z2)
g22=Z2%*%A%*%t(Z2)


# investigate the negative pleiotropy between ypd37 and other traits
peak.index[[44]]$markerIndex

As=A.mat(gdata[,peak.index[[44]]$markerIndex], shrink=F)/2

yr=c(pdata[,44], pdata[,31])
g1a=Z1%*%As%*%t(Z1)
g21a=Z2%*%As%*%t(Z1)
g12a=Z1%*%As%*%t(Z2)
g22a=Z2%*%As%*%t(Z2)

g1=Z1%*%A%*%t(Z1)
g21=Z2%*%A%*%t(Z1)
g12=Z1%*%A%*%t(Z2)
g22=Z2%*%A%*%t(Z2)

yrra=regress(scale(yr)~1,~g1a+g12a+g21a+g22a+E1+E2, identity=F, verbose=T)
yrr=regress(scale(yr)~1,~g1+g12+g21+g22+E1+E2, identity=F, verbose=T)

lm(yr[1:1008]~gdata[,peak.index[[44]]$markerIndex])
lm(yr[1009:2016]~gdata[,peak.index[[44]]$markerIndex])
coef(lm(yr[1009:2016]~gdata[,peak.index[[44]]$markerIndex]))
plot(coef(lm(yr[1:1008]~gdata[,peak.index[[44]]$markerIndex])), coef(lm(yr[1009:2016]~gdata[,peak.index[[44]]$markerIndex])),
     xlab='effect size growth in YPD 37C', ylab='effect size growth in trehalose'  
     )

#with QTL
0.6043 

#withoutQTL

summary(lm(yr~c(rep(-1,1008),rep(1,1008))*rbind(gdata[,peak.index[[44]]$markerIndex], gdata[,peak.index[[44]]$markerIndex])))

load('/media/kserver/kruglyak/raid1/home/jbloom/1000BYxRM/QTL/fQTLs_FDR05r.bin')


sapply(fQTLs_FDR05r, function(x) x$fqtl$ests$ests[-1])

eff.size.r=sapply(fQTLs_FDR05r, function(x) x$fqtl$ests$ests[-1])

eff.size.cil=sapply(fQTLs_FDR05r, function(x) x$CIs$leftgpos)
eff.size.cir=sapply(fQTLs_FDR05r, function(x) x$CIs$rightgpos)


unlist(eff.size.r)

ii=Intervals(cbind(unlist(eff.size.cil), unlist(eff.size.cir)))
plot(ii, col=ifelse(unlist(eff.size.r)<0,'orange', 'purple'), names_cex=.5)






df2=data.frame(unlist(eff.size.r),unlist(eff.size.cil),unlist(eff.size.cir))
df2[order(df2[,2]),]




#scale sigmas


#sigmas[['6_14']]
"Congo_red"
"Hydrogen_Peroxide" 
R> length(m1)
[1] 979
R> length(m2)
[1] 769
0.6907 0.2244 0.2244 0.5946 0.3732 0.4665
m1=which(!is.na(pdata[,6]))
m2=which(!is.na(pdata[,14]))

Z0=matrix(0, 769,769 )
Z1=rbind(diag(979),Z0)
Z2=rbind(Z0, diag(769))
#summary(

g1=Z1%*%A%*%t(Z1)
g21=Z2%*%A%*%t(Z1)
g12=Z1%*%A%*%t(Z2)
g22=Z2%*%A%*%t(Z2)
E1=Z1%*%t(Z1)
E2=Z2%*%t(Z2)







summary(

X2=cbind(c(rep(1,1008), rep(0,1008)), c(rep(0,1008), rep(1,1008)))
        r=regress(y~X2,~g1+g12+g21+g22+E1+E2, identity=F, verbose=T)
str

calcMM = function(y, B=NULL,X=NULL, Z=NULL, Ze=NULL, reps=FALSE,
                     alg='ai', conv.val=1e-4, Var=NULL){
    strain.names=(names(y))
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)

 #   strain.ind  = seq_along(strain.names)
 #   strain.cnt  = length(unique.sn)
 #   #for constructing Strain Variance component
 #   Strain      = Matrix(diag(strain.cnt), sparse=T)

    # If B is NULL then add a covariance term that is the identity matix - will then calculate effect of strain (broad-sense heritability)
 #   if(is.null(B)) { B=list(Strain=Strain);  } else{
 #       if (reps) { B=c(B, list(Strain=Strain))  } }    
 #   # If Z is null and there are no replicates this will make Z a diagonal incidence matrix, otherwise this constructs an incidence matrix based on strain names
 #   if(is.null(Z)) {   Z=Matrix(0, length(y), strain.cnt,sparse=T);   Z[cbind(strain.ind, n.to.m)]=1 }
 #   # If X is null assume one fixed effect of population mean
 #   if(is.null(X) ) {  X=model.matrix(y~1)}
 #   # If Ze is null assume no error structure
 #   if(is.null(Ze)) {  Ze=Matrix((diag(length(y))),sparse=T) }

    #number of terms in the structured covariance
 #   VC.names=paste('sigma', c(names(B), 'E'), sep='')
        # second derivatives of V with respect to the variance components (Lynch and Walsh 27.15)

 #   VV = list()
 #   for(i in 1:N.s.c) {
 #            VV[[i]]=Z %*% tcrossprod(B[[i]],Z)     }
 #   VV[[ Vcmp.cnt ]]=Ze 
 
    
    X=model.matrix(y~1)
    which.pos=c(1,2,5,6)
    VV= list(g1=g1,g22=g22, g12=g12,g21=g21, E1=E1,E2=E2)

    Vcmp.cnt=length(VV)
    #N.s.c+1
    # starting values for VC estimates as 1 / (#of VCs including residual error term)
  # Var=rep(1/Vcmp.cnt, Vcmp.cnt) 
    Var=c(.5,.5,0,0,.5,.5)
    I = matrix(0, ncol= Vcmp.cnt, nrow= Vcmp.cnt)
	s = matrix(0, ncol=1, nrow= Vcmp.cnt)
    
    diffs=rep(10,  Vcmp.cnt)

#g1=Z1%*%A%*%t(Z1)
#g21=Z2%*%A%*%t(Z1)
#g12=Z1%*%A%*%t(Z2)
#g22=Z2%*%A%*%t(Z2)
#E1=Z1%*%t(Z1)
#E2=Z2%*%t(Z2)
VC.names=names(VV)

    i = 0
    # while the differences haven't converged 
    while ( sum(ifelse(diffs<conv.val, TRUE,FALSE)) <  Vcmp.cnt ) { 
		i = i + 1
        V=matrix(0,length(y), length(y))
	    for( vcs in 1:length(VV)) {  V=V+(VV[[vcs]]*Var[vcs]) }
             print('Inverting V')
             Vinv = solve(V)
             print('Done inverting V')
            tXVinvX=t(X) %*% Vinv %*% X
            inv.tXVinvX = solve(tXVinvX)
            itv = inv.tXVinvX %*% t(X)%*%Vinv
            P = Vinv - Vinv %*% X %*% itv 

        if(alg=='fs') {print("Fisher scoring algorithm: calculating expected VC Hessian") }
        if(alg=='nr') {print("Netwon rhapson algorithm: calculating observed VC Hessian") }
        if(alg=='ai') {print("Average information algorithm: calculating avg of expected and observed VC Hessians") }
        
        for(ii in 1:Vcmp.cnt) {
           for(jj in ii:Vcmp.cnt) {
                 if (alg=='fs') {    I[ii,jj]= 0.5*sum(diag( ((P%*%VV[[ii]]) %*%P )%*%VV[[jj]])) }
                 if (alg=='nr') {    I[ii,jj]=-0.5*sum(diag(P%*%VV[[ii]]%*%P%*%VV[[jj]])) + (t(y)%*%P%*%VV[[ii]]%*%P%*%VV[[jj]]%*%P%*%y)[1,1] }
                 if (alg=='ai') {   I[ii,jj]= 0.5*( t(y)%*%P%*%VV[[ii]]%*%P%*%VV[[jj]]%*%P%*%y)[1,1]                 }
                 print(paste(ii, jj))
                 I[jj,ii]=I[ii,jj]
           }
             s[ii,1]= -0.5*sum(diag(P%*%VV[[ii]])) + .5*(t(y)%*%P%*%VV[[ii]]%*%P%*%y )[1,1]
        }
        invI = solve(I)
        print(invI)
        print(s) 
        newVar = Var + invI%*%s
        # remove constraint that newvar has to be greater than 0 
       # for(nn in which.pos) {  if( newVar[nn]<0 ) {newVar[nn]=2.6e-9} }
        #newVar[newVar<0]=2.6e-9

        for(d in 1:length(diffs)) { diffs[d]=abs(Var[d] - newVar[d]) }
		Var = newVar
        
        cat('\n')
        cat("iteration ", i, '\n')
        cat(VC.names, '\n')
        cat(Var, '\n')
        cat("SE \n")
        cat(sqrt(diag(invI)), '\n')

        
        Bhat= itv %*% y
        cat("Fixed Effects, Bhat = ", as.matrix(Bhat), '\n')
        det.tXVinvX=determinant(tXVinvX, logarithm=TRUE)
        det.tXVinvX=det.tXVinvX$modulus*det.tXVinvX$sign
        det.V =determinant(V, logarithm=TRUE)
        det.V=det.V$modulus*det.V$sign

        LL = -.5 * (det.tXVinvX + det.V + t(y) %*% P %*% y )
        cat("Log Likelihood = " , as.matrix(LL), '\n')
        cat("VC convergence vals", '\n')
        cat(diffs, '\n')
	}
    	cat('\n')
	return(list(Var=Var, invI=invI, W=Vinv, Bhat=Bhat, llik=LL))
}


























y3=c(pdata[,c(2,3)])
y1=pdata[,2]
y2=pdata[,3]
Z0=matrix(0, 1008,1008)
Z1=rbind(diag(1008),Z0)
Z2=rbind(Z0, diag(1008))

summary(regress(y1~1,~A, verbose=T, pos=c(T,T)))
summary(regress(y2~1,~A, verbose=T, pos=c(T,T)))

g1=Z1%*%A%*%t(Z1)
g21=Z2%*%A%*%t(Z1)
g12=Z1%*%A%*%t(Z2)
g22=Z2%*%A%*%t(Z2)


summary(regress(y3~1,~g1+g12+g22, verbose=T, pos=c(T,T,T,T)))










# fix this
all.strain.names=names(pheno_raw$Copper)

#nsample  =nrow(gdata)
#nmarker  =ncol(gdata)

newMM.results=list()
newMM.cblups =list()
newMM.peaks =list()

for(phenotype in names(pheno_raw)[-c(1,39,42)] ) {
    s=(pheno_raw[[phenotype]])
    # fix this -----------------------------------
    srle=sapply(s, length)
    sr=list()
    sr$lengths=as.vector(srle)
    sr$values=names(srle)
    attr(sr, 'class')='rle'

    ny=as.vector(unlist(s))
    names(ny)=inverse.rle(sr)
    #-----------------------------------------------

    y=ny[!is.na(ny)]

    strain.names=(names(y))
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)
    strain.ind  = seq_along(strain.names)
    strain.cnt  = length(unique.sn)
    #for constructing Strain Variance component
    Strain      = diag(strain.cnt)
    Z=matrix(0, length(y), strain.cnt);   Z[cbind(strain.ind, n.to.m)]=1 

    strains.with.phenos=match(unique.sn, all.strain.names)

    n.strains=length(strains.with.phenos)
    # calculate relatedness matrix from genomewide marker data
    # also, in this case, this is equivalent to (gdata %*% t(gdata))
    A = A.mat(gdata[strains.with.phenos,], shrink=FALSE)/2
    AA=A*A

    #list to store relatedness matrices built per chromosome
    A.chr=list()
    # cov matrices for genome minus a chromosome
    A.chr.m = list()
    ZA=list()
    ZAm=list()
    for(i in unique(chr.lab)){
        print(i)
        A.chr[[i]]=A.mat(gdata[strains.with.phenos,chr.lab==i], shrink=FALSE)/2
        ZA[[i]]=(Z %*% A.chr[[i]] %*% t(Z))
        A.chr.m[[i]]=(A.mat(gdata[strains.with.phenos,chr.lab!=i], shrink=FALSE)/2)
        ZAm[[i]]=(Z %*% A.chr.m[[i]] %*% t(Z))
    }
    names(A.chr.m)=paste('m_', names(A.chr), sep='')
    names(ZA)=names(A.chr)
    names(ZAm)=names(A.chr.m)

    ZAtZ=Z%*%A%*%t(Z)
    ZAAtZ=Z%*%AA%*%t(Z)
    ZStraintZ=Z %*% Strain %*% t(Z)
    # or 
    #~factor(names(y)) to regress instead of ZStraintZ
    y.avg=as.vector(by(y, names(y), mean, na.rm=T))
    # Calculate LODS, avg data
    LODS.avg=get.LOD.by.COR(n.strains, 
                            y.avg, gdata[strains.with.phenos,])
    options(warn = -1)

   
    # Broad - sense heritability
    mm.broad = regress(y~1, ~ZStraintZ, pos=c(TRUE,TRUE),verbose=T)
    mm.broad = extractVarCompResults(mm.broad)
   
    mm.A = regress(y~1, ~ZAtZ+ZStraintZ, pos=c(TRUE,TRUE,TRUE),verbose=T)
    mm.A     = extractVarCompResults(mm.A)

    mm.AA    = regress(y~1, ~ZAtZ+ZAAtZ+ZStraintZ, pos=c(TRUE,TRUE,TRUE,TRUE), verbose=T)
    mm.AA     = extractVarCompResults(mm.AA)
    
    # equivalent, maybe even faster and more memory efficient
    #mm.AA = calcMM(y, B=list(A=A, AA=AA), reps=TRUE)
    print(phenotype)
    print(mm.AA)
    #newMM.results[[phenotype]]=mm.AA
    # for each chromosome calc model with epistasis and rest of genome with and without that chromosome
    
    m.chr.blups.wo = list()
    mm.null.lliks = rep(NA,16)

    for ( i in 1:16) {
     print(i)
       # for calcMM
       #mm.null = calcMM(y, B=list(A.m=A.chr.m[[i]], AA=AA), conv.val=1e-4)
       #b = calc.BLUPS(mm.null$Var[1]*A.chr.m[[i]]+ mm.null$Var[2]*AA, Z, mm.null$W,  y, rep(1,length(y)),mm.null$Bhat)         
       #y.all = calc.BLUPS(mm.null$Var[1]*A.chr.m[[i]]+ mm.null$Var[2]*AA + mm.null$Var[3]*Strain+mm.null$Var[4]*diag(n.strains), Z, mm.null$W,  y, rep(1,length(y)),mm.null$Bhat)         
        mm.null = regress(y~1, ~ZAm[[i]]+ZAAtZ+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)
        b    = calc.BLUPS(mm.null$sigma['ZAm[[i]]']*A.chr.m[[i]] + mm.null$sigma['ZAAtZ']*AA , Z, mm.null$W,  y, rep(1,length(y)),mm.null$beta)
        y.all= calc.BLUPS(mm.null$sigma['ZAm[[i]]']*A.chr.m[[i]] + mm.null$sigma['ZAAtZ']*AA + mm.null$sigma['ZStraintZ']*Strain+ mm.null$sigma['In']*diag(n.strains), Z, mm.null$W, y, rep(1,length(y)),mm.null$beta)
        m.chr.blups.wo[[i]]=y.all-b
        mm.null.lliks[i]=mm.null$llik
    }

    c.blups=do.call('cbind', m.chr.blups.wo)
    newMM.cblups[[phenotype]]=c.blups
    #mm.full.all.chr=regress(y~1, ~ZA[[1]]+ZA[[2]]+ZA[[3]]+ZA[[4]]+ZA[[5]]+ZA[[6]]+ZA[[7]]+ZA[[8]]+ZA[[9]]+ZA[[10]]+ZA[[1]]+ZA[[12]]+ZA[[13]]+ZA[[14]]+ZA[[15]]+ZA[[16]]+ZStraintZ, pos=rep(TRUE,18), verbose=9)

   
    LODS=get.LOD.by.COR(n.strains, c.blups, gdata[strains.with.phenos,])

    #plot(LODS[1,],ylim=c(0, max(LODS)) )
    #for(i in 2:16) { points(LODS[i,], col=i) }

    uc=unique(chr.lab)
    peaks=c()  

    #h2 nqtlh2 nqtlint epph2
    pdfout=paste('~/Dropbox/new_mm/', phenotype, '.pdf', sep='')
    pdf(file=pdfout, width=11, height=8)
    for ( i in 1:16){
        print(i)
        chr.ind=seq(1,11623)[chr.lab==uc[i] ]
        # for simulated QTL
        #sim.add.ind =  which(chr.ind %in% add.qtl.ind) 
        #sim.add.sign =add.qtl.sign[add.qtl.ind %in% chr.ind]
        
        g.s = gdata[strains.with.phenos,chr.ind]
        r.vec= c.blups[,i]
        plot.max=max(LODS[i,][which(chr.lab==uc[i])])
        plot(LODS[i,][which(chr.lab==uc[i])], ylim=c(0, plot.max+2), type='l',lwd=3, main=paste('Chr', i), ylab='LOD')
        #for simulated QTL
        #text(sim.add.ind, 0, '*', col=sim.add.sign+3, cex=4)
        mp  = which.max(LODS[i,][chr.ind])
        m.l = max(LODS[i,][chr.ind])
        perm.thresh=quantile(apply(get.LOD.by.COR(n.strains, replicate(1000,sample(r.vec)), g.s),1,max), .975)
        peaks.picked=0
        while(m.l>perm.thresh) {
            abline(v=mp, col=peaks.picked+1)
            peaks.picked=peaks.picked+1
            peaks=c(peaks, chr.ind[mp])
            r.vec=residuals(lm(r.vec~g.s[,mp]-1))
            lod.resid = as.vector(get.LOD.by.COR(n.strains, r.vec, g.s) )
            points(lod.resid, type='l', col=peaks.picked+1)
            mp  = which.max(lod.resid)
            m.l = max(lod.resid)
            perm.thresh=quantile(apply(get.LOD.by.COR(n.strains, replicate(1000,sample(r.vec)), g.s),1,max), .975)
        }
       # naive
       points(LODS.avg[1,][which(chr.lab==uc[i])], type='l',lwd=2, lty=3)
    }
    dev.off()
    newMM.peaks[[phenotype]]=peaks
    
    peak.A=A.mat(gdata[strains.with.phenos,peaks], shrink=FALSE)/2
    Zp = Z %*% peak.A %*% t(Z)
    pAA=peak.A*peak.A
    ZpAA = Z %*% pAA %*% t(Z)
    pAG = peak.A *A 
    ZpAG = Z %*% pAG %*% t(Z)
    
    mm.peaks = regress(y~1, ~Zp+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)
    mm.peaks= extractVarCompResults(mm.peaks)

    mm.peaks.AA = regress(y~1, ~Zp+ZpAA+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)
    mm.peaks.AA = extractVarCompResults(mm.peaks.AA)
    mm.peaks.AAg =  regress(y~1, ~Zp+ZpAG+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)
    mm.peaks.AAg = extractVarCompResults(mm.peaks.AAg)
 
    newMM.results[[phenotype]]=list(mm.broad=mm.broad, mm.A=mm.A, mm.AA=mm.AA, 
                                    mm.peaks=mm.peaks, mm.peaks.AA=mm.peaks.AA, mm.peaks.AAg =mm.peaks.AAg)

    results=list(newMM.peaks, newMM.results, newMM.cblups)

    save(file='~/Dropbox/new_mm/Results.RData',results )
}



    error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y, angle=90, code=3, length=length, ...)
    }




sf = t(sapply(lapply(newMM.results, function(x) {x$mm.AA$sigma}), function(x) {sum(x) }))

ses=sapply(newMM.results, function(x) {sqrt(diag(x$mm.AA$sigma.cov))})
for(i in 1:ncol(ses) ) {ses[,i]=ses[,i]/sf[i] } 

otb = order(txx.b)


xx=t(sapply(lapply(newMM.results, function(x) {x$mm.AA$sigma}), function(x) {x/sum(x)}))
txx= t(xx)
txx.b=apply(txx, 2, function(x) {sum(x[1:3])})
tb=txx[,otb]
ses=sapply(newMM.results, function(x) {sqrt(diag(x$mm.AA$sigma.cov))})
sf = t(sapply(lapply(newMM.results, function(x) {x$mm.AA$sigma}), function(x) {sum(x) }))
for(i in 1:ncol(ses) ) {ses[,i]=ses[,i]/sf[i] } 
ses.o = ses[,otb]
ys=apply(tb,2, cumsum)

pdf(file='~/Dropbox/barplot_whole_genome_A_AA_1000.pdf', width=11, height=8)
par(mar=c(12,4,4,2))
b=barplot(tb,las=2, col=c('blue', 'red', 'pink', 'grey'), ylab='fraction of phenotypic variance')
error.bar(b, ys[1,], ses.o[1,], col='blue')
error.bar(b, ys[2,], ses.o[2,], col='red')
error.bar(b, ys[3,], ses.o[3,], col='pink')
dev.off()

mmAA.table=tb
mmAA.table.se=ses.o

xx=t(sapply(lapply(newMM.results, function(x) {x$mm.peaks.AA$sigma}), function(x) {x/sum(x)}))
txx= t(xx)
txx.b=apply(txx, 2, function(x) {sum(x[1:3])})
tb=txx[,otb]
ses=sapply(newMM.results, function(x) {sqrt(diag(x$mm.peaks.AA$sigma.cov))})
sf = t(sapply(lapply(newMM.results, function(x) {x$mm.peaks.AA$sigma}), function(x) {sum(x) }))
for(i in 1:ncol(ses) ) {ses[,i]=ses[,i]/sf[i] } 
ses.o = ses[,otb]
ys=apply(tb,2, cumsum)

pdf(file='~/Dropbox/barplot_peaks_A_AA__1000.pdf', width=11, height=8)
par(mar=c(12,4,4,2))
b=barplot(tb,las=2, col=c('blue', 'red', 'pink', 'grey'), ylab='fraction of phenotypic variance')
error.bar(b, ys[1,], ses.o[1,], col='blue')
error.bar(b, ys[2,], ses.o[2,], col='red')
error.bar(b, ys[3,], ses.o[3,], col='pink')
dev.off()

mmAApeak.table=tb
mmAApeak.table.se=ses.o



xx=t(sapply(lapply(newMM.results, function(x) {x$mm.peaks.AAg$sigma}), function(x) {x/sum(x)}))
txx= t(xx)
txx.b=apply(txx, 2, function(x) {sum(x[1:3])})
tb=txx[,otb]
ses=sapply(newMM.results, function(x) {sqrt(diag(x$mm.peaks.AAg$sigma.cov))})
sf = t(sapply(lapply(newMM.results, function(x) {x$mm.peaks.AAg$sigma}), function(x) {sum(x) }))
for(i in 1:ncol(ses) ) {ses[,i]=ses[,i]/sf[i] } 
ses.o = ses[,otb]
ys=apply(tb,2, cumsum)

pdf(file='~/Dropbox/barplot_peaks_A_peaks_genomeAA_1000.pdf', width=11, height=8)
par(mar=c(12,4,4,2))
b=barplot(tb,las=2, col=c('blue', 'red', 'pink', 'grey'), ylab='fraction of phenotypic variance')
error.bar(b, ys[1,], ses.o[1,], col='blue')
error.bar(b, ys[2,], ses.o[2,], col='red')
error.bar(b, ys[3,], ses.o[3,], col='pink')
dev.off()


mmAAGpeak.table=tb
mmAAGpeak.table.se=ses.o
library(gplots)


pdf(file='~/Dropbox/1000BYxRMvarcomp.pdf', width=10, height=10)
plotCI( mmAA.table[1,], mmAApeak.table[1,], 
        uiw=mmAA.table.se[1,] ,err='x',
       xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',
     xlab='A genome', ylab='A QTL', barcol='grey')
plotCI( mmAA.table[1,], mmAApeak.table[1,], 
        uiw=mmAApeak.table.se[1,] ,err='y',
       xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',
     xlab='A genome', ylab='A QTL', add=T, barcol='grey')
abline(0,1)


plotCI( mmAA.table[2,], mmAApeak.table[2,], 
        uiw=mmAA.table.se[2,] ,err='x',
       xlim=c(0,.4), ylim=c(0,.4), xaxs='i', yaxs='i',
     xlab='AA genome', ylab='AA QTL', barcol='grey')
plotCI( mmAA.table[2,], mmAApeak.table[2,], 
        uiw=mmAApeak.table.se[2,] ,err='y',
       xlim=c(0,.4), ylim=c(0,.4), xaxs='i', yaxs='i',
     xlab='AA genome', ylab='AA QTL', add=T, barcol='grey')
abline(0,1)

plotCI( mmAA.table[2,], mmAAGpeak.table[2,], 
        uiw=mmAA.table.se[2,] ,err='x',
       xlim=c(0,.4), ylim=c(0,.4), xaxs='i', yaxs='i',
     xlab='AA genome', ylab='AA (QTL*genome)', barcol='grey')
plotCI( mmAA.table[2,], mmAAGpeak.table[2,], 
        uiw=mmAAGpeak.table.se[2,] ,err='y',
       xlim=c(0,.4), ylim=c(0,.4), xaxs='i', yaxs='i',
     xlab='AA genome', ylab='AA (QTL*genome)', add=T, barcol='grey')
abline(0,1)
dev.off()





p1=cross$pheno[,1]
p1[is.na(p1)]=0

1	25	9535	-0.1974	0.1284	-0.0801	0.0493	0.9999

1	1	208	-0.1996	0.1548	0.0159	0.1758	0.2134
cross2= cross
cross2$pheno[,1]=p1

test=scantwo(cross2, pheno.col=1, chr=1, method='mr')



plot( mmAA.table[2,], mmAApeak.table[2,], xlim=c(0,.4), ylim=c(0,.4), xaxs='i', yaxs='i',
     ylab='AA genome', xlab='AA QTL')
abline(0,1)

plot( mmAA.table[2,], mmAAGpeak.table[2,], xlim=c(0,.4), ylim=c(0,.4), xaxs='i', yaxs='i',
     ylab='AA genome', xlab='A QTL * A genome')
abline(0,1)


xxg=t(sapply(lapply(newMM.results, function(x) {x$mm.peaks.AAg$sigma}), function(x) {x/sum(x)}))






extractScaledPhenotype=function(impcross,scaleVar=FALSE){apply(impcross$pheno, 2, scale, scale=scaleVar)}

gdata     = extractGenotype(cross)
gdata.sub=cbind(gdata, gdata)
#downsample gdata by 4
pdata.01  = extractScaledPhenotype(cross, TRUE)
pdata.01  = rbind(pdata.01, pdata.01,pdata.01,pdata.01)
n.pheno   = countStrainsPerTrait(cross$pheno)*4
#gdata.sub = gdata[,seq(1, ncol(gdata), 4)]
#rm(cross)
#rm(gdata)

gsub.name=do.call('rbind', strsplit(colnames(gdata) ,'_'))
gsub.pos=as.numeric(gsub.name[,1])
tick.spot=c(0, cumsum(rle(gsub.name[,2])$lengths)[-16])+1
tick.name=rle(gsub.name[,2])$values

gsub.ind.chr.split= split(1:nrow(gsub.name), match(gsub('chr', '', gsub.name[,2]), as.character(as.roman(c(1:16)))))
#mcor=rep(NA,ncol(gdata)-1)
#for(i in 2:ncol(gdata) ){ print(i); mcor[i]=cor(gdata[,i-1], gdata[,i], use='pairwise.complete.obs') } 

# calculate 2-locus model with residuals from additive polygenic model
max.marker=ncol(gdata.sub)


#foreach(i=1:(max.marker-50)) %dopar% {
for(i in 1:(max.marker-50)) {
system.time({
    for(i in 1:10){
    print(i)
    s2=get.LOD.by.COR(n.pheno, pdata.01, gdata.sub[,i]*gdata.sub[,(i+50):max.marker])
   # save(s2, file=paste('/data/4000BYxRM/s2_addrm/', i, sep=''))
} })




pdata.addrm.perm=cbind(pdata.addrm[sample(1:nrow(pdata.addrm)),],
                       pdata.addrm[sample(1:nrow(pdata.addrm)),],
                       pdata.addrm[sample(1:nrow(pdata.addrm)),],
                       pdata.addrm[sample(1:nrow(pdata.addrm)),],
                       pdata.addrm[sample(1:nrow(pdata.addrm)),])
n.pheno.perm=rep(n.pheno,5)
# calculate 2-locus model with residuals from additive polygenic model for permuted data 
foreach(i=1:(max.marker-50)) %dopar% {
    print(i)
    s2=get.LOD.by.COR(n.pheno.perm, pdata.addrm.perm, gdata.sub[,i]*gdata.sub[,(i+50):max.marker])
    save(s2, file=paste('/data/4000BYxRM/s2_addrm_perm5/', i, sep=''))
}

s1=get.LOD.by.COR(n.pheno, pdata.01, gdata.sub)
s2.array= array(NA,c(21,7055,7055))
for(i in 1:7005) {
    print(i)
    load(paste('/data/4000BYxRM/s2_addrm/', i, sep=''))
    s2.array[,i,(50+i):max.marker]=s2
    rm(s2)
}
s2.peaks = getIntPeaks(s2.array, gdata.sub, gsub.ind.chr.split)

#load permutation data
perm.ind=rep(1:5,each=21)
s2.array.perm = array(NA,c(21,7055,7055))
s2.perm.peaks=list()
for(p in unique(perm.ind)) {
    for(i in 1:7005) {
      print(i)
      load(paste('/data/4000BYxRM/s2_addrm_perm5/', i, sep=''))
      s2.array.perm[,i,(50+i):max.marker]=s2[perm.ind==p,]
      rm(s2)
    }
    s2.perm.peaks[[p]]=getIntPeaks(s2.array.perm, gdata.sub, gsub.ind.chr.split)
}
rm(s2.array.perm)
#----------------------------
# calc avg expected
perm.lods.s2=lapply(s2.perm.peaks, function(x) {sapply(x, function(y) {y[,3] } ) })
iecnt=sapply(perm.lods.s2, function(x) { sapply(seq(2,7,.05), function(tt) {sum(x>tt)}) })
rownames(iecnt)=seq(2,7,.05)
iecnt=apply(iecnt,1,mean)

obscnt=sapply(seq(2,7,.05), function(tt) {sum(sapply(s2.peaks, function(x){x[,3]})>tt) } )
iecnt/obscnt















#X11()
    #pdf(file='~/Dropbox/testLODS.pdf', width=11, height=8)
    #plot(LODS.avg[1,], xlab='marker index', ylab='LOD',  col='black', type='l')
    #dev.off()
    #abline(v=match(unique(chr.lab), chr.lab),lty=2)
    #text(add.qtl.ind, 0, '*', col=add.qtl.sign+3, cex=2)
    #abline(v=add.qtl.ind, col=add.qtl.sign+3)




#testing 2-locus scan strategies

#y.avg
n.strains=length(strains.with.phenos)
y.fe=lm(y.avg~gdata[strains.with.phenos,peaks]-1)
y.resid.fixed=residuals(y.fe)
peak.A=A.mat(gdata[strains.with.phenos,peaks], shrink=FALSE)/2
Zp = Z %*% peak.A %*% t(Z)
pAA=peak.A*peak.A
ZpAA = Z %*% pAA %*% t(Z)
pAG = peak.A *A 
ZpAG = Z %*% pAG %*% t(Z)
g.w=gdata[strains.with.phenos,peaks]

#for(i in 1:ncol(g.w)){g.w[,i]=coefficients(y.fe)[i]*g.w[,i] }
#peak.Aw=A.mat(g.w, shrink=FALSE)/2
#Zpw = Z %*% peak.Aw %*% t(Z)
#pAAw=peak.Aw*peak.Aw
#ZpAAw = Z %*% pAAw %*% t(Z)
#pAGw = peak.Aw *A 
#ZpAGw = Z %*% pAGw %*% t(Z)



mm.AA    = regress(y~1, ~ZAtZ+ZAAtZ+ZStraintZ, pos=c(TRUE,TRUE,TRUE,TRUE), verbose=T)
mm.peaks = regress(y~1, ~Zp+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)
mm.peaks.AA = regress(y~1, ~Zp+ZpAA+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)
mm.peaks.AAg =  regress(y~1, ~Zp+ZpAG+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)


#mm.peaksw = regress(y~1, ~Zpw+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)
#mm.peaks.AAw = regress(y~1, ~Zpw+ZpAAw+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)
#mm.peaks.AAgw =  regress(y~1, ~Zpw+ZpAGw+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)


s1=get.LOD.by.COR(n.strains, y.avg, gdata[strains.with.phenos,])
s1.i197 =get.LOD.by.COR(n.strains, y.avg, gdata[strains.with.phenos,]*gdata[strains.with.phenos,197])
y2=y.avg
y2[1:900]=NA
s1.i197 =get.LOD.by.COR(n.strains, cbind(y.avg, y2), gdata[strains.with.phenos,]*gdata[strains.with.phenos,197])

s2.test= scantwo(cross, pheno.col=7, method='mr', chr=c(1,1))


[1] "194861_chr01_194861_A_T"
#197
R> find.marker(cross,1,114)
[1] "196192_chr01_196192_A_G"
#201

gg=gdata[strains.with.phenos, c(197, 201)]
aov(lm(y.avg~gg), lm(y.avg~gg[,1]*gg[,2]))

R> summary(s2.test, what='int')
      pos1 pos2 lod.full lod.fv1 lod.int lod.add lod.av1
c1:c1  113  114     7.62    5.07    2.61       5    2.45


i=197
j=201
#anova(lm(cbind(y.avg)~g.sub[,i]+g.sub[,j]), lm(cbind(y.avg)~g.sub[,i]*g.sub[,j]))
#Analysis of Variance Table

#Model 1: cbind(y.avg) ~ g.sub[, i] + g.sub[, j]
#Model 2: cbind(y.avg) ~ g.sub[, i] * g.sub[, j]
#  Res.Df  RSS Df Sum of Sq    F Pr(>F)
#1    969 3138                         
#2    968 3138  1     0.157 0.05   0.83

lm(y.avg~g.sub[,i]*g.sub[,j])
200*200
r1=rnorm(972)


ph=cbind(cross$pheno[,1], cross$pheno[,2],cross$pheno[,3])

for(i in 1:208){
    print(i)
    for(j in i:208){
        #aa=lm.fit(g.sub[,c(i,j)], y.avg)
        #ii=lm.fit(cbind(g.sub[,c(i,j)], gdata[,i]*gdata[,j]) ,y.avg)  
        #aa=lm(ph~ gdata[,c(i,j)])
        #ii=lm(y.avg ~ cbind( g.sub[,c(i,j)], g.sub[,i]*g.sub[,j]))    
        # .avg
         anova(lm(y.avg~g.sub[,i]+g.sub[,j]), lm(y.avg~g.sub[,i]*g.sub[,j]))[[6]][2]        
    
    }
}


anova(lm(y.avg~g.sub[,i]*g.sub[,j]))


mm.null = regress(y~1, ~ZAm[[i]]+ZAAtZ+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)
            b    = calc.BLUPS(mm.null$sigma['ZAm[[i]]']*A.chr.m[[i]] + mm.null$sigma['ZAAtZ']*AA + mm.null$sigma['ZStraintZ']*Strain, Z, mm.null$W,  y, rep(1,length(y)),mm.null$beta)
            y.all= calc.BLUPS(mm.null$sigma['ZAm[[i]]']*A.chr.m[[i]] + mm.null$sigma['ZAAtZ']*AA + mm.null$sigma['ZStraintZ']*Strain+ mm.null$sigma['In']*diag(length(strains.with.phenos)), Z, mm.null$W, y, rep(1,length(y)),mm.null$beta)
        m.chr.blups.wo[[i]]=y.all-b





