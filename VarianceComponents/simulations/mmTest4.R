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
library(Rfast)
#library(gputools)
#source('~/Dropbox/calcMM.R')
# some useful functions
doTraitFDR=function(trait, genos,  FDR_thresh=.05, nperm=100) {
    f.found=c()
    p.found=c()
    q.found=c()
    m.found=c()

    n=length(trait)
    L= (crossprod(trait,genos)/(n-1))^2 
    mLi=which.max(L)
    mL=max(L)
    
    yperm=replicate(nperm, sample(trait))
    nullD=(crossprod(yperm,genos)/(n-1))^2
    
    permMax=rowMaxs(nullD,value=T)
    pNull=1-ecdf(permMax)(mL)
    if(pNull==0) {pNull=1/nperm}
    
    step=1
    
    repeat{
       p.temp=c(p.found, pNull)
       q=-mean(log(1-p.temp))
       if(q>FDR_thresh) {break;}
       p.found=c(p.found, pNull)
       q.found=c(q.found, q)
       m.found=c(m.found, colnames(genos)[mLi])
       f.found=c(f.found, mLi)
       #print(paste('step=', step, 'max index=', colnames(genos)[mLi], 'max r^2=', mL, 'pnull=', pNull, 'fdr=', q))
       yr=scale(residuals(lm(trait~genos[,f.found]) ))
       L=(crossprod(yr,genos)/(n-1))^2 
       mLi=which.max(L)
       mL=max(L)
       yperm=replicate(nperm, sample(yr))
       nullD=(crossprod(yperm,genos)/(n-1))^2
       permMax=rowMaxs(nullD, value=T) 
       pNull=1-ecdf(permMax)(mL)
       if(pNull==0) {pNull=1/nperm}
       step=step+1
   }
   results=data.frame(fscan.markers=m.found, index=f.found, p=p.found, q=q.found, stringsAsFactors=F) 
   return(results)
}
extractGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$data }))*2)-3 }
extractGenotype.argmax=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$argmax }))*2)-3 }
countStrainsPerTrait = function(pheno) {apply(pheno, 2, function(x){sum(!is.na(x))})}
get.LOD.by.COR = function(n.pheno, pheno, gdata) {
    # Lynch and Walsh p. 454
    return( (-n.pheno*log(1-cor(pheno, gdata, use='pairwise.complete.obs')^2))/(2*log(10)) ) }
get.LOD.by.COR.gpu = function(n.pheno, pheno, gdata) {
    # Lynch and Walsh p. 454
    return( (-n.pheno*log(1-gpuCor(pheno, gdata, use='pairwise.complete.obs')$coefficients^2))/(2*log(10)) ) }
calc.pval.LRtest=function(null,full) { pchisq(-2*(null-full),1, lower.tail=FALSE) }
calc.BLUPS= function(G,Z,Vinv,y,X,B ){    ((G%*%t(Z)) %*% Vinv) %*%( y - X%*%B)     }
#extract genotype data
gdata     = extractGenotype(cross)
gdata.s=scale(gdata)
A=tcrossprod(gdata.s)/ncol(gdata)
AA=A*A

library(MASS)
MA=mvrnorm(6000,mu=rep(0,1008),Sigma=.4*A+.6*diag(1008))

fm=crossprod(scale(t(MA)), gdata.s)/1007
LOD=-1008*(log(1-fm^2))/(2*log(10))
flod=which(LOD>6, arr.ind=T)
 plot(flod[,2], flod[,1])

gcoord=as.numeric(do.call('rbind', strsplit(colnames(gdata), '_'))[,1])
# chr name
chr=as.character(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])

#for 1000 BYxRM-------------------------------------------------------------
chr=gsub('chr0', '', chr)
chr=gsub('chr', '', chr)
chr=as.numeric(chr)
# chr pos (ignore for now)
# chr name again, fixed for sorting
chr.lab=do.call('rbind', strsplit(colnames(gdata), '_'))[,2]
#--------------------------------------------------------------------------

#for 4000 BYxRM-------------------------------------------------------------
#chr = gsub("chr", "", chr)
#chr.lab = match(chr, as.character(as.roman(1:16)))

pos=as.numeric(do.call('rbind', strsplit(colnames(gdata), '_'))[,3])

nsample  =nrow(gdata)
nmarker  =ncol(gdata)


# for simulations #################################################################
all.marker.combo.ind=t(which(upper.tri(matrix(0,nmarker,nmarker)), arr.ind=T))
#mchosen=seq(1,ncol(all.marker.combo.ind), 10)
#bigI=matrix(0, 1008,length(mchosen))
#for(i in 1:length(mchosen)){
#    bigI[,i]=gdata[all.marker.combo.ind[1, mchosen[i]]]*gdata[all.marker.combo.ind[2, mchosen[i]]]
#    if(i%%1000==0) {print(i)}
#}
#nn=apply(all.marker.combo.ind[,mchosen], 2, paste, collapse=':')
#colnames(bigI)=nn

H2=.6
h2=.6

nadditive = 20
nint      = 20
only.incl = FALSE
#total number of additive loci
#total number of interacting loci
a.eff  = rep(0,nsample)
aa.eff = rep(0,nsample)
#markers with effects
add.qtl.ind  = sort(sample(nmarker, nadditive))
#add.qtl.sign = sample(rnorm(nadditive),replace=T)
add.qtl.sign = sample(ifelse(runif(nadditive)>.5,1,-1),replace=T)

#for assigning QTL
#int.qtl.combos=combn(add.qtl.ind,2) 
#if(ncol(int.qtl.combos)<nint)  {
#    int.qtl.combos=cbind(int.qtl.combos, all.marker.combo.ind[,sort(sample(ncol(all.marker.combo.ind), nint-ncol(int.qtl.combos)))])
#}

#int.qtl.ind=sort(sample(ncol(int.qtl.combos), nint))
#int.qtl.sign=sample(c(1,-1),nint,replace=T)


for(i in 1:nadditive){ a.eff=a.eff+add.qtl.sign[i]*gdata[,add.qtl.ind[i]] }
#for(i in 1:nint) {aa.eff=aa.eff+int.qtl.sign[i]*gdata[,int.qtl.combos[,int.qtl.ind[i]][1]]*
#                                                gdata[,int.qtl.combos[,int.qtl.ind[i]][2]]  }

a.eff=scale(a.eff)
#aa.eff=scale(aa.eff)

g=sqrt(h2)*a.eff #+ sqrt(H2-h2)*aa.eff
y=g+rnorm(nrow(gdata),mean=0,sd=sqrt((1-H2)/H2*var(g)))
y=scale(as.vector(y))

library(knockoff)
Xf=create.gaussian(gdata, rep(0,11623), diag(11623)  )

library(spaMM)
id=as.factor(rownames(A))
Ai=solve(A)
obs=data.frame(y=y, id=id)
prec=as_precision(A)
AA=A*A
fitme(y~1+corrMatrix(1|id), covStruct=list(precision=Ai), data=obs, method='REML')
test=fitme(y~1+corrMatrix(1|id), covStruct=list(corrMatrix=A), data=obs, method='PQL/L')

library(lme4qtl)
dat=data.frame(y=y[,1], id=as.character(1:1008))
colnames(A)=as.character(1:1008)
rownames(A)=as.character(1:1008)

test=relmatLmer(y~1+(1|id), relmat=list(id=A), dat) #, family='gaussian') #'poisson') #gaussian') #poisson')


Betas=cor(y,gdata.s)[1,]
csplit=split(data.frame(t(gdata.s)),chr)
lcor=lapply(csplit, function(x) crossprod(as.matrix(t(x)))/(ncol(x)-1))

#LD=cor((gdata.s))
V=crossprod(gdata.s)/(nrow(gdata.s)-1)
svd.V=svd(V)


library(L0Learn)
L0fit=L0Learn.cvfit(gdata.s, y, penalty="L0", algorithm='CDPSI') 
lambdaIndex = which.min(L0fit$cvMeans[[1]]) 
L0coef = as.numeric(coef(L0fit$fit, lambda = L0fit$fit$lambda[[1]][lambdaIndex]))
effect.beta = L0coef[which(L0coef!=0)][-1]
effect.index = (which(L0coef!=0)-1)[-1] 
length(effect.beta)
set.seed(1)
s.init = susie_init_coef(effect.index, effect.beta, 11623)
# could attempt using FDR procedure for s.init
s.init2 = susie_init_coef(test$index, as.vector(coef(lm(y~gdata.s[,test$index]-1))), 11623)
susieL0.fit = susie(gdata.s,y,s_init=s.init,coverage=.99,verbose=T,compute_univariate_zscore=T)
susieL0.fit2 = susie(gdata.s,y,s_init=s.init2,coverage=.95,verbose=T,compute_univariate_zscore=T)

susie_plot(susieL0.fit2,'z')
abline(v=add.qtl.ind, col='red')
 sum(sapply(sapply(susieL0.fit2$sets$cs, function(x) x%in% add.qtl.ind), sum))


sapply(susieL0.fit2$sets$cs, function(x) gcoord[x])
sapply(sapply(susieL0.fit2$sets$cs, function(x) gcoord[x]), function(y) max(y)-min(y))
library(susieR)





calib <- susie_auto(gdata.s,y,verbose=T, compute_univariate_zscore=T)
res = susie(gdata.s,y,L=32, verbose=T, coverage=.99,estimate_prior_variance=F, track_fit=T, compute_univariate_zscore=T)
#susie_plot(calib,'z')
susie_plot(res,'z')
abline(v=add.qtl.ind, lwd=2, col='grey', lty=2)
points(test$index, -log10(pnorm(-abs(res$z[test$index]))), col='red', pch='*', cex=2)
points(abs((L0fit$fit$beta[[1]][,27])*30), cex=3, col='blue')   

       3*abs(res$z[test$index])

       
       -log10(res$z[test$index]+.5) , col='red', pch="*")


Z=cor(y,gdata.s)[1,]
library(bigstep)
test2=selectModel(gdata.s, y, fastST=T, crit='mbic2', p=ncol(gdata.s), maxf=100)


library(selectiveInference)
test=fs(gdata.s,y, maxsteps=100,intercept=F, normalize=F, verbose=T)
t2=fsInf(test)
forwardStop(t2, alpha=0.05)

library(PhenotypeSimulator)

test=regress(y~1, ~A, verbose=T)
Z=cor(y,gdata.s)[1,]

# test grouped lasso
library(gglasso)
test=cv.gglasso(gdata.s,y, chr, pred.loss='L1')
#BAGS
t2=BAGS(X=gdata.s, y=y, pi=1)


g=gdata+1
library(PhenotypeSimulator)
writeStandardOutput('/home/jbloom/Downloads/sra/', genotypes=g, phenotypes=matrix(y), kinship=A,
                    id_samples=seq(1,1008), id_snps=colnames(g), id_phenos='y', outstring='t1000', 
                    format='bimbam')
writeStandardOutput('/home/jbloom/Downloads/sra/', genotypes=g, phenotypes=cbind(y,y), kinship=A,
                    id_samples=as.character(seq(1,1008)), id_snps=colnames(g), id_phenos=c('y1', 'y2'), outstring='p1000', 
                    format='plink')
# then convert
# plink --noweb --bfile genotypesp1000 --recode --make-bed

~/Local/DPR -g genotypest1000.bimbam -p Ysimt1000_bimbam.txt -dpr 2 -nk 2 -o test

~/Local/DPR -g genotypest1000.bimbam -p Ysimt1000_bimbam.txt -k relatedness.sXX.txt -dpr 1 -nk 4 -o test -m 0

~/Local/DPR -bfile genotypesp1000 -p Ysimt1000_bimbam.txt -k relatedness.sXX.txt -dpr 1 -nk 4 -o test -m 0
~/Local/DPR -bfile plink -p Ysimt1000_bimbam.txt -k relatedness.sXX.txt -dpr 1 -nk 4 -o test -m 0


~/Local/DPR -g genotypest1000.bimbam -p Ysimt1000_bimbam.txt -k relatedness.sXX.txt -dpr 3 -nk 4 -o test -m 0 -o mcmc
#~/Local/bayesR/bin/bayesR -bfile plink -n 1 -numit 10000 -burnin 5000

-out br -model br1 -n 1 -nthreads 12

x=read.delim('~/Downloads/sra/output/test.param.txt', header=T)
abline(v=add.qtl.ind)


test=r.REML(y, matrix(1,1008), A, method='h2', intercept=FALSE, ubound=3)

gt=t(gdata)

E=diag(1008)
#Ainv=solve(A)
result=rep(NA,50)
j=1
for(i in seq(.3,.95, length.out=50)) {
    print(i)
    V=i*A+(1-i)*E
    Vinv=solve(V)
    blup=(i*A)%*%(Vinv)%*%(y-mean(y))
    #tr=rmvnorm(10000, sigma=V)
    result[j]=cor(y,blup)^2
        #sum((y-blup)^2)
    #print(sum((y-blup)^2))
    j=j+1
}
plot(seq(.3,.95, length.out=50), log(result))
library(brms)
library(rstan)
rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run on multiple cores

data=data.frame(y=y,G=as.factor(seq(1,1008)))
rownames(A)=data$G
colnames(A)=data$G 
AA=A*A
brm(y~1 + (1|G), family=gaussian(), cov_ranef=list(G=A), data=data)

Gnew=matrix(0, 1008, 1008)
for(i in 1:1008){
    print(i)
    gvec=gdata[i,]
    Gnew[i,]=colSums(gvec==gt)
}

P=diag(1008)-((matrix(1,1008)%*%t(matrix(1,1008)))/11623)
denom=sum(diag(P%*%Gnew%*%P))

Ak=2*((1007*Gnew)/denom)


GnewN=(Gnew/ncol(gdata)) #*2
GnewN=2*GnewN/mean(diag(GnewN))

regress(y~1, ~A, verbose=T)
regress(y~1, ~Ak, verbose=T)

regress(y~1, ~GnewN, verbose=T)

GN=A/(sum(diag(A))/1008)

GN=Gnew/(sum(diag(Gnew))/1008)


Sn=(1008*Gnew)/denom
regress(y~1, ~Sn, verbose=T)
,

lm(as.vector(A)~as.vector(GnewN))



yds=y[seq(1,1008,2)]
GnewNds=GnewN[seq(1,1008,2),seq(1,1008,2)]
regress(y~1, ~GnewN, verbose=T)
regress(yds~1, ~GnewNds, verbose=T)

regress(y~1, ~A, verbose=T)

    gt=t(gdata)






##https://github.com/cheuerde/BayesianAlphabet/blob/master/JagsBayesR.r
#Z=gdata
#X <- array(1, dim = c(nrow(Z), 1))
#library('VIGoR')

#library(regress)
#regress(y~1, ~A+AA, verbose=T)
#summary(lm(y~gdata[,add.qtl.ind]-1))
#logodds <- seq(-3,-1,0.1)
#XX=cbind(gdata,bigI)

yy=cross$pheno[,7]
yy[is.na(yy)]=mean(yy,na.rm=T)
y=yy
#fitR=varbvsmix(X=gdata, Z=NULL,sa=c(0,0.1,.2236,.5)^2,  y=yy,verbose=T, tol=1e-5) #, sigma=.2, sa=0.8)#, logodds=logodds)
#plot(rowSums(fitR$alpha*fitR$mu))
#fitR$w*11623
#max(fitR$logZ)


#library(vbsr)
library('varbvs')
par(mfrow=c(3,1))
gLOD=get.LOD.by.COR(1008,y,gdata)[1,]
plot(gLOD)
#segments(add.qtl.ind, 0, add.qtl.ind, 5, col=ifelse(add.qtl.sign>0, 3, 2)) #add.qtl.sign+3)
abline(h=3, col='grey')
#fit1=varbvsmix(X=gdata, Z=NULL,sa=c(0,0.2,0.5)^2,  y=y,verbose=T, tol=1e-5) #, sigma=.2, sa=0.8)#, logodds=logodds)
fit1=varbvsmix(X=gdata, Z=NULL, sa=c(0,0.1,0.2,0.5)^2,  y=y,verbose=T, tol=1e-6, drop.threshold=1e-8 )#, sigma=.2, sa=0.8)#, logodds=logodds)
c(0,1e-3, 1e-2,.05)
fit2=varbvsmix(X=gdata, Z=NULL, sa=c(0,.03163,.1,.2,.5)^2,  y=y,verbose=T, tol=1e-12, drop.threshold=1e-9 )#, sigma=.2, sa=0.8)#, logodds=logodds)
sum(c(fit1$w*11623)[-1])
max(fit1$logZ)
max(fit2$logZ)

#fit1=varbvsmix(X=gdata, Z=NULL,sa=c(0,0.5)^2,  y=y,verbose=T, tol=1e-5) #, sigma=.2, sa=0.8)#, logodds=logodds)
fit1$w*11623
plot(-log10(fit1$lfsr))
plot(rowSums(fit1$alpha*fit1$mu), col=ifelse(fit1$lfsr< .05, 'red', 'black'), cex=2)
max(fit1$logZ)
#fit=varbvs(X=gdata, Z=NULL, y=y,verbose=T) #, sigma=.2, sa=0.8)#, logodds=logodds)
#w <- normalizelogweights(fit$logw)
#pip  <- fit$alpha %*% c(w)
#beta <- fit$mu %*% c(w)
#sa <- sum(fit$sa * w)
#y.fit <- predict(fit,gdata)
#print(cor(y,y.fit))
#varbvs
#vars=which(pip>.9)










max(fit1$logZ)
fsignset=which(fit1$lfsr< .50)
closestA=sapply(fsignset, function(x) add.qtl.ind[which.min(abs(x-add.qtl.ind))] )
#cbind(closestA, fsignset, closestA-fsignset)

#points(add.qtl.ind, coef(lm(y~gdata[,add.qtl.ind]-1)), col='blue')
segments(add.qtl.ind, 0, add.qtl.ind, ifelse(add.qtl.sign>0, 1, -1)) #, col=ifelse(add.qtl.sign+3)

points( fsignset, rep(0, length(fsignset)), col='black', cex=2, pch=20)

plot(gdata %*% (rowSums(fit1$alpha*fit1$mu)), y)
cor.test(gdata %*% (rowSums(fit1$alpha*fit1$mu)), y)
#
fit1$w*11623
sum(fit1$w[-1]*11623)


library(BGLR)
nIter=20000; burnIn=2000; thin=20
g02=gdata+1

#fmBA=BGLR(y=y,ETA=list( list(X=g02,model='BayesA')),             nIter=nIter,burnIn=burnIn, R2=.8)
fmBB=BGLR(y=y,ETA=list( list(X=g02,model='BayesB')),             nIter=nIter,burnIn=burnIn, thin=thin, R2=.8)
fmBC=BGLR(y=y,ETA=list(list(X=g02,model='BayesC')), nIter=nIter,burnIn=burnIn, thin=thin)

fmBC$ETA[[1]]$probIn*11623 
#+/-
1.96*fmBC$ETA[[1]]$SD.probIn*11623 


0.00328

#fmBL=BGLR(y=y,ETA=list( list(X=g02,model='BL')),             nIter=nIter,burnIn=burnIn, R2=.8)
#library(NAM)
#test=emBC(y, g02, R2=.8, Pi=.01)

plot(abs(fmBC$ETA[[1]]$b),col=4,cex=.5, type='o',main='BayesC');

library(gdmp)
bc=BayesCpi(g02, y=y, numiter=10000)
plot(bc$meanb[-1])













fit2=vbsr(y, gdata #, sigma=.2, sa=0.8)#, logodds=logodds)
segments(add.qtl.ind, 0, add.qtl.ind, 5, col=add.qtl.sign+3)
plot(-log10(fit2$pval))
segments(add.qtl.ind, 0, add.qtl.ind, 5, col=add.qtl.sign+3)
#points(add.qtl.ind, coef(lm(y~gdata[,add.qtl.ind]-1)), col='blue', pch=20)
replicate(100, vbsr(sample(y), gdata)
plot(-log10(fit3$pval))

library(EBglmnet)
fitEB=cv.EBglmnet(gdata, scale(y), verbose=5)

library(varbvs)
fit=varbvs(X=gdata, Z=NULL, y=y,sigma=0.2, sa=0.8, verbose=T) #, sigma=.2, sa=0.8)#, logodds=logodds)
fit2=varbvs(X=gdata, Z=NULL, y=y,sigma=0.2, sa=0.8, logodds=gLOD,verbose=T) #, sigma=.2, sa=0.8)#, logodds=logodds)








plot(rowSums(fit1$alpha*fit1$mu))

w <- normalizelogweights(fit$logw)
pip  <- fit$alpha %*% c(w)
beta <- fit$mu %*% c(w)
sa <- sum(fit$sa * w)
y.fit <- predict(fit,gdata)
print(cor(y,y.fit))
#varbvs
vars=which(pip>.9)

varbvscoefcred(fit, vars)

cis=varbvscoefcred(fit, vars, cred.int = 0.95, nr = 10000)
ignoreCor=varbvsindep(fit, gdata, Z=NULL, y)
par(mfrow=c(2,1))
plot(gLOD)
segments(add.qtl.ind, 0, add.qtl.ind, 1, col=add.qtl.sign+3)
abline(v=vars)

plot(beta)
points(add.qtl.ind, coef(lm(y~gdata[,add.qtl.ind]-1)), col='blue', pch=20)

plot(-log10(fit1$lfsr))
abline(v=add.qtl.ind)


library(gdmp)
bc=BayesCpi(g02, y=y)
bc2=BayesCpi(g02, y=y, Pi=0.99)

plot(bc2$meanb[-1])

bca=BayesAB(g02, y=y)

library(BGLR)
nIter=6000; burnIn=1000
fmBA=BGLR(y=y,ETA=list( list(X=g02,model='BayesA')),             nIter=nIter,burnIn=burnIn, R2=.8)
fmBB=BGLR(y=y,ETA=list( list(X=g02,model='BayesB')),             nIter=nIter,burnIn=burnIn, R2=.8)
fmBC=BGLR(y=y,ETA=list( list(X=g02,model='BayesC')),             nIter=nIter,burnIn=burnIn, R2=.8)
fmBL=BGLR(y=y,ETA=list( list(X=g02,model='BL')),             nIter=nIter,burnIn=burnIn, R2=.8)



 plot(abs(fmBC$ETA[[1]]$b),col=4,cex=.5, type='o',main='BayesC');
abline(v=QTL,col=2,lty=2)

AAA=A*A*A
AAAA=A*A*A*A
AAAAA=A*A*A*A*A

test.pheno=sapply(pheno_raw[[2]], function(x) mean(x, na.rm=T))
test.pheno=sapply(pheno_raw[[46]], function(x) mean(x, na.rm=T))
test.pheno[is.na(test.pheno)]=mean(test.pheno,na.rm=T)

y=scale(test.pheno)
r=regress(y~1, ~A,verbose=T)

lmb=r$sigma[2]/(r$sigma[1]/ncol(gdata))
lr=lm.ridge(y~gdata.s-1, lambda=lmb)

mm=solve(gdata.s%*%t(gdata.s)+lmb*diag(1008))
mb=t(gdata.s)%*%mm%*%y
library(MASS)
library(pls)
pc2=pcr(y~A,1007)


pacc=rep(0,1007)
for(i in 1:1007) {
pacc[i]=cor(y,predict(pc2)[,,i])
}

lr2=lm.ridge(y~A, lambda=lmb)
plot(lr2$coef, BLUP(r)$Mean)

r2=regress(y~1, ~AA,verbose=T)
r3=regress(y~1, ~AAA,verbose=T)
r4=regress(y~1, ~AAAA,verbose=T)
r5=regress(y~1, ~AAAAA,verbose=T)

r2.1=regress(y~1,~A+AA,verbose=T)
r2.2=regress(y~1, ~A+AA+AAA,verbose=T)
D1=diag(1008)
r2.3=regress(y~1, ~D1,verbose=T)

plot(r$predicted, y)
plot(r2$predicted, y)
plot(r3$predicted, y)
plot(r4$predicted, y)
plot(r5$predicted, y)


gcoord=as.numeric(do.call('rbind', strsplit(colnames(gdata), '_'))[,1])
# chr name
chr=as.character(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])

#for 1000 BYxRM-------------------------------------------------------------
chr=gsub('chr0', '', chr)
chr=gsub('chr', '', chr)
chr=as.numeric(chr)
# chr pos (ignore for now)
# chr name again, fixed for sorting
chr.lab=do.call('rbind', strsplit(colnames(gdata), '_'))[,2]
#--------------------------------------------------------------------------

#for 4000 BYxRM-------------------------------------------------------------
#chr = gsub("chr", "", chr)
#chr.lab = match(chr, as.character(as.roman(1:16)))

pos=as.numeric(do.call('rbind', strsplit(colnames(gdata), '_'))[,3])

nsample  =nrow(gdata)
nmarker  =ncol(gdata)


# for simulations #################################################################
all.marker.combo.ind=t(which(upper.tri(matrix(0,nmarker,nmarker)), arr.ind=T))

n.replicates=1e6

H2=.02
h2=.01

nadditive = 40
nint      = 80
only.incl = FALSE
#total number of additive loci
#total number of interacting loci
a.eff  = rep(0,nsample)
aa.eff = rep(0,nsample)
#markers with effects
add.qtl.ind  = sort(sample(nmarker, nadditive))
add.qtl.sign = sample(c(1,-1),nadditive,replace=T)

#for assigning QTL
int.qtl.combos=combn(add.qtl.ind,2) 
if(ncol(int.qtl.combos)<nint)  {
    int.qtl.combos=cbind(int.qtl.combos, all.marker.combo.ind[,sort(sample(ncol(all.marker.combo.ind), nint-ncol(int.qtl.combos)))])
}

int.qtl.ind=sort(sample(ncol(int.qtl.combos), nint))
int.qtl.sign=sample(c(1,-1),nint,replace=T)


for(i in 1:nadditive){ a.eff=a.eff+add.qtl.sign[i]*gdata[,add.qtl.ind[i]] }
for(i in 1:nint) {aa.eff=aa.eff+int.qtl.sign[i]*gdata[,int.qtl.combos[,int.qtl.ind[i]][1]]*
                                                gdata[,int.qtl.combos[,int.qtl.ind[i]][2]]  }

a.eff=scale(a.eff)
aa.eff=scale(aa.eff)

g=sqrt(h2)*a.eff + sqrt(H2-h2)*aa.eff
s=replicate(n.replicates, {g+rnorm(nrow(gdata),mean=0,sd=sqrt((1-H2)/H2*var(g)))})[,1,]
rownames(s)=paste(letters[1:26], 1:1008, sep='_')

#y=as.vector(do.call('c', data.frame(s)))
#names(y)=rep(rownames(s),n.replicates)

#subs=c(1,2,3)
#y=as.vector(do.call('c', data.frame(s[,subs])))
#names(y)=rep(rownames(s),length(subs))

#######################################################################################

strain.names=(names(y))
unique.sn=unique(strain.names)
n.to.m=match(strain.names, unique.sn)
strain.ind  = seq_along(strain.names)
strain.cnt  = length(unique.sn)
#for constructing Strain Variance component
Strain      = diag(strain.cnt)
Z=matrix(0, length(y), strain.cnt);   Z[cbind(strain.ind, n.to.m)]=1 

# calculate relatedness matrix from genomewide marker data
# also, in this case, this is equivalent to (gdata %*% t(gdata))
A = A.mat(gdata, shrink=FALSE)/2
AA=A*A


#list to store relatedness matrices built per chromosome
#A.chr=list()
# cov matrices for genome minus a chromosome
#A.chr.m = list()
#for(i in unique(chr.lab)){
#    print(i)
#    A.chr[[i]]=A.mat(gdata[,chr.lab==i], shrink=FALSE)/2
#    A.chr.m[[i]]=(A.mat(gdata[,chr.lab!=i], shrink=FALSE)/2)
#}
#names(A.chr.m)=paste('m_', names(A.chr), sep='')

ZAtZ=Z%*%A%*%t(Z)
ZAAtZ=Z%*%AA%*%t(Z)
ZStraintZ=Z %*% Strain %*% t(Z)
# or 
#~factor(names(y)
# cov matrices for each chromosome separately 
#ZA=list()
#ZAm=list()
#for(i in 1:length(A.chr)){
#    print(i)
#    ZA[[i]]=(Z %*% A.chr[[i]] %*% t(Z))
#    ZAm[[i]]=(Z %*% A.chr.m[[i]] %*% t(Z))
#}
#names(ZA)=names(A.chr)
#names(ZAm)=names(A.chr.m)

infl.sim=list()
#y.avg=y
for(i in c(1e1,1e2,1e3,1e4,1e5,1e6)) {
    y.avg= apply(s[,1:i],1,mean)
    #y.avg=s[,1]
    #y.avg=apply(s[,1],1,mean)
#run mixed model with y.avg
    mm.ymean=regress(y.avg~1, ~A+AA, pos=c(T,T,T), verbose=T)
    infl.sim[[i]]=(mm.ymean$sigma/sum(mm.ymean$sigma))
}
cis=do.call('cbind', infl.sim)
pdf(file='~/Dropbox/rep_infl_3.pdf', width=10, height=8)
b=barplot(cis,las=1, col=c('blue', 'red', 'grey'), 
          ylab='fraction of phenotypic variance', names.arg=c(1,1e1,1e2,1e3,1e4,1e5,1e6),main='effect of  averaging replicates on VC estimates'
          )
abline(h=.01)
abline(h=.02)
dev.off()

# Calculate LODS, avg data
LODS.avg=get.LOD.by.COR(nsample, y.avg, gdata)
#X11()
plot(LODS.avg[1,], xlab='marker index', ylab='LOD',  col='black', type='l')
abline(v=match(unique(chr.lab), chr.lab),lty=2)
text(add.qtl.ind, 0, '*', col=add.qtl.sign+3, cex=2)
#abline(v=add.qtl.ind, col=add.qtl.sign+3)

# Calculate LODS, one replicate
#LODS=get.LOD.by.COR(nsample, s[,1], gdata)
#plot(LODS[1,], xlab='marker index', ylab='LOD', main='One Replicate' ,col='grey',type='l')
#abline(v=match(unique(chr.lab), chr.lab),lty=2)
#abline(v=add.qtl.ind, col=add.qtl.sign+3)



#save.image(file='~/Desktop/mm_map/mm.sim.RData')
options(warn = -1)

# nullest null
mm=regress(y~1,verbose=T)
#-1593

# Broad - sense heritability
mm.broad = regress(y~1, ~ZStraintZ, pos=c(TRUE,TRUE),verbose=T)
mm.A = regress(y~1, ~ZAtZ+ZStraintZ, pos=c(TRUE,TRUE,TRUE),verbose=T)
mm.AA = regress(y~1, ~ZAtZ+ZAAtZ+ZStraintZ, pos=c(TRUE,TRUE,TRUE,TRUE),verbose=T)

# comparisons of predictions from blup models and QTL models 

B.blups = calc.BLUPS(mm.broad$sigma['ZStraintZ']*Strain, Z, mm.broad$W, y, rep(1,length(y)), mm.broad$beta)

A.blups=calc.BLUPS(mm.A$sigma['ZAtZ']*A, Z, mm.A$W, y, rep(1,length(y)),mm.A$beta)
AA.blups=calc.BLUPS(mm.AA$sigma['ZAtZ']*A+mm.AA$sigma['ZAAtZ']*AA, Z, mm.AA$W, y, rep(1,length(y)),mm.AA$beta)

AQ = A.mat(gdata[,add.qtl.ind], shrink=FALSE)/2
ZAQtZ=(Z%*%AQ) %*%t(Z)
mm.Q = regress(y~1, ~ZAQtZ+ZStraintZ, pos=c(TRUE,TRUE,TRUE), verbose=T)
QTL.blups = calc.BLUPS(mm.Q$sigma['ZAQtZ']*AQ, Z, mm.Q$W, y, rep(1,length(y)),mm.Q$beta)

mm.QAA = regress(y~1, ~ZAQtZ+ZAAtZ+ZStraintZ, pos=c(TRUE,TRUE,TRUE), verbose=T)
QTLAA.blups = calc.BLUPS(mm.QAA$sigma['ZAQtZ']*AQ +mm.QAA$sigma['ZAAtZ']*AA , Z, mm.QAA$W, y, rep(1,length(y)),mm.QAA$beta)


cor(A.blups,    B.blups)
cor(AA.blups,   B.blups)
cor(QTL.blups,  B.blups)
cor(QTLAA.blups,  B.blups)


# for each chromosome calc model with epistasis and rest of genome with and without that chromosome
#m.chr.blups.w  = list()
m.chr.blups.wo = list()
mm.null.lliks = rep(NA,16)
mm.full.lliks = rep(NA,16)

for ( i in 1:16) {
    print(i)
    mm.null = regress(y~1, ~ZAm[[i]]+ZAAtZ+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE),verbose=2)
    b    = calc.BLUPS(mm.null$sigma['ZAm[[i]]']*A.chr.m[[i]] + mm.null$sigma['ZAAtZ']*AA , Z, mm.null$W,  y, rep(1,length(y)),mm.null$beta)
    y.all= calc.BLUPS(mm.null$sigma['ZAm[[i]]']*A.chr.m[[i]] + mm.null$sigma['ZAAtZ']*AA + mm.null$sigma['ZStraintZ']*Strain+ mm.null$sigma['In']*diag(length(s[,1])), Z, mm.null$W, y, rep(1,length(y)),mm.null$beta)
    m.chr.blups.wo[[i]]=y.all-b
    mm.null.lliks[i]=mm.null$llik
    rm(mm.null)
    gc()
    #save(mm.null, file=paste('~/Desktop/mm_map/chr_null_',i, sep=''))
    #mm.full = regress(y~1, ~ZA[[i]]+ZAm[[i]]+ZAAtZ+ZStraintZ, pos=c(TRUE, TRUE, TRUE, TRUE,TRUE),verbose=2)
    #    b    = calc.BLUPS(mm.full$sigma['ZAm[[i]]']*A.chr.m[[i]] + mm.full$sigma['ZStraintZ']*Strain+mm.full$sigma['ZAAtZ']*AA, Z, mm.full$W, y, rep(1,length(y)),mm.full$beta)
    #    y.all= calc.BLUPS(mm.full$sigma['ZA[[i]]']*A.chr[[i]] + mm.full$sigma['ZAm[[i]]']*A.chr.m[[i]] + mm.full$sigma['ZAAtZ']*AA + mm.full$sigma['ZStraintZ']*Strain+mm.full$sigma['In']*diag(length(s[,1])), Z, mm.full$W, y, rep(1,length(y)),mm.full$beta)
    #    print(mm.full)
    #m.chr.blups.w[[i]]=y.all-b
    #mm.full.lliks[i]=mm.full$llik
    #save(mm.full, file=paste('~/Desktop/mm_map/chr_full_',i, sep=''))
}
#chr.ps=calc.pval.LRtest(mm.null.lliks,mm.full.lliks)

c.blups=do.call('cbind', m.chr.blups.wo)


#mm.full.all.chr=regress(y~1, ~ZA[[1]]+ZA[[2]]+ZA[[3]]+ZA[[4]]+ZA[[5]]+ZA[[6]]+ZA[[7]]+ZA[[8]]+ZA[[9]]+ZA[[10]]+ZA[[1]]+ZA[[12]]+ZA[[13]]+ZA[[14]]+ZA[[15]]+ZA[[16]]+ZStraintZ, pos=rep(TRUE,18), verbose=9)

LODS=get.LOD.by.COR(nsample,c.blups, gdata)

plot(LODS[1,],ylim=c(0, max(LODS))) 
for(i in 2:16) { points(LODS[i,], col=i) }

uc=unique(chr.lab)
peaks=c()  

#h2 nqtlh2 nqtlint epph2
pdf(file='~/Desktop/mm_map/sim_oneqtl_blups_2.pdf', width=11, height=8)
for ( i in 1:16){
    print(i)
    chr.ind=seq(1,11623)[chr.lab==uc[i] ]
    # for simulated QTL
    sim.add.ind =  which(chr.ind %in% add.qtl.ind) 
    sim.add.sign =add.qtl.sign[add.qtl.ind %in% chr.ind]
    
    g.s = gdata[,chr.ind]
    r.vec= c.blups[,i]
    plot.max=max(LODS[i,][which(chr.lab==uc[i])])
    plot(LODS[i,][which(chr.lab==uc[i])], ylim=c(0, plot.max+2), type='l',lwd=3, main=paste('Chr', i), ylab='LOD')
    text(sim.add.ind, 0, '*', col=sim.add.sign+3, cex=4)
    mp  = which.max(LODS[i,][chr.ind])
    m.l = max(LODS[i,][chr.ind])
    perm.thresh=quantile(apply(get.LOD.by.COR(nsample, replicate(1000,sample(r.vec)), g.s),1,max), .975)
    peaks.picked=0
    while(m.l>perm.thresh) {
        abline(v=mp, col=peaks.picked+1)
        peaks.picked=peaks.picked+1
        peaks=c(peaks, chr.ind[mp])
        r.vec=residuals(lm(r.vec~g.s[,mp]-1))
        lod.resid = as.vector(get.LOD.by.COR(nsample, r.vec, g.s) )
        points(lod.resid, type='l', col=peaks.picked+1)
        mp  = which.max(lod.resid)
        m.l = max(lod.resid)
        perm.thresh=quantile(apply(get.LOD.by.COR(nsample, replicate(1000,sample(r.vec)), g.s),1,max), .975)
    }
   # naive
   points(LODS.avg[1,][which(chr.lab==uc[i])], type='l',lwd=2,lty=3)
}
dev.off()





##understanding Eskin method
#y=s[,1]
#ml=mixed.solve(y,K=A,method='ML')
##ml$Vu/(ml$Vu+ml$Ve)
##[1] 0.4428
#G=.4638
#delta=ml$Vu/ml$Ve
#
#eig.A=eigen(A)
#
#ei
#
#calc.mmodel.nllik=function(par,data){
#   #B=0
#   G=par[1]
#   delta=par[2]
#   with(data, {
#   #H=Z %*% A %*% t(Z)+delta*II
#   H = ((eig.A$vectors %*% diag(eig.A$values+delta)) %*% t(eig.A$vectors))
#   Hinv = (eig.A$vectors %*% diag(1/(eig.A$values+delta)) %*% t(eig.A$vectors))
#   V=(1/G) * H
#   Vinv= G * Hinv
#   ldetV= determinant(V)$modulus * determinant(V)$sign
#   # maximize log likelihood (minimize the -1*(log likelihood))
#   nll=-(-(n/2)*log(2*pi)-(1/2)*ldetV-(1/(2))*t(y-(X*B)) %*% (G*Hinv) %*% (y-(X*B)))
#   print(paste(
#            'h2= ', round(G/(G+(G*delta)),3),
#            'B = ', round(B, 2),
#            'sigma.g =', round(G,2),
#            'delta =', round(delta,2),
#            'llik =', round(-nll,6)))
#   return(nll)
#   })
#  }
#
#calc.mmodel.nllik.l=function(par,data){
#   #B=0
#   #G=par[1]
#   delta=par[1]
#   with(data, {
#   #H=Z %*% A %*% t(Z)+delta*II
#   H = ((eig.A$vectors %*% diag(eig.A$values+delta)) %*% t(eig.A$vectors))
#   Hinv = (eig.A$vectors %*% diag(1/(eig.A$values+delta)) %*% t(eig.A$vectors))
#   Px = Hinv- Hinv %*%X %*% (solve(t(X)%*%Hinv%*%X)) %*% t(X) %*% Hinv
#   ldetH = determinant(H)$modulus * determinant(H)$sign
#   nll = (n/2) * log(n/(2*pi)) - (n/2) -.5 *ldetH - (n/2) *  log( (t(y) %*% Px) %*% y)
#   tau=n/((t(y)%*%Px%*%y))
#   # maximize log likelihood (minimize the -1*(log likelihood))
#   print(paste(
#            'delta =' , round(delta,4), 
#            'tau =' ,   round(tau,4),
#            'llik =', round(-nll,6)))
#   return(nll=nll) 
#   })
#  }
#
#
#mmjb= optim(par=c(.0001), method='Brent',
#            calc.mmodel.nllik.l, 
#            lower=c(1e-3), 
#            upper=c(1),
#            data=list(X=matrix(rep(1,length(y))),n=length(y), eig.A=eig.A))
#
#
#
#II=diag(length(y)),Z=diag(length(y)), A=A))
















