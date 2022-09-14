library(Matrix)
library(data.table)
library(mvnfast)
library(snpStats)
library(Rfast)

#https://github.com/biotangweiqi/PQHE/blob/master/estimate_QTL_heritability.v0.4.10.R
## for genetic map
crossRdata.file='/data/single_cell_eQTL/yeast/reference/cross.list.RData'
#pearson R to LOD given sample size 
RtoLOD=function(n.pheno, r) {(-n.pheno*log(1-r^2))/(2*log(10)) }

LODToPval =  function(x){  pchisq(x*(2*log(10)),df=1,lower.tail=FALSE)/2 }
PvalToLOD =  function(x){ l= qchisq(x*2,df=1,lower.tail=FALSE)/(2*log(10)) 
                          l[is.na(l)]=0
return(l)
} 

#some functions for processing genetic maps from 2019 Bloom Elife paper 
extractGeneticMapDf=function(cross){
    preliminary.map=sapply(cross$geno, function(x) x$map)
    # we need to add some jitter here 
    gmap=data.frame(marker=as.vector(unlist(sapply(preliminary.map, names))),
                    gpos=as.vector(unlist(sapply(preliminary.map, function(x)x))),stringsAsFactors=F)
    gmap2=tstrsplit(gmap$marker,'_')
    gmap$chrom=gmap2[[1]]
    gmap$ppos=as.numeric(gmap2[[2]])
    return(gmap) 
} 
jitterGmap=function(themap, amount=1e-6) {
    for (i in 1:length(themap)) {
         n <- nrow(themap[[i]])
         themap[[i]]$map <- themap[[i]]$map + c(0, cumsum(rep(amount, n - 1)))
    }
    return(themap)
}
extractGmaps=function(crossRdata.file) {
    load(crossRdata.file)
    gmaps=list()
    for(cname in names(cross.list)){
        print(cname)
        A=cross.list[[cname]]
        #rm(cross.list)
        a.gmap=extractGeneticMapDf(A)
        str(a.gmap)
        chroms=paste0('chr', as.roman(1:16)) 
        gmap.s=a.gmap; colnames(gmap.s)[2]='map'
        gmap.s=split(gmap.s, gmap.s$chrom)
        gmap.s=gmap.s[chroms]
        gmap.s=jitterGmap(gmap.s)
        gmaps[[cname]]=gmap.s
    }
    return(gmaps)
}


# calculate the effective number of test given a pairwise relatedness matrix or LD matrix  
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


def_block_pos <- function(chr, size){
	# block number
	n <- as.integer(2*chr/size)+2;
	if( n%%2 != 0 ) n <- n+1;
	n <- as.integer(n/2);
	# block index and the middle position of each block
	i <- c(1:n);
	pos <- (i-1)*size+floor(size/2); # middle position
	if(pos[n]>chr) pos[n]<-chr;
	return(pos);
}

# meta value of block
#https://academic.oup.com/bioinformatics/article/36/7/2150/5631910
cal_block_meta <- function(loc, val, chr, size, depth, MIN){
	# input: location vector, value vector
	# input: chr length, block size, location depth vector
	pos <- def_block_pos(chr, size);
	idx <- as.integer(0.5+loc/size)+1;
	#
	avg <- c();
	blockDepth <- c();
	for(i in 1:length(pos)){
		k  <- which(idx==i);
		no <- length(k);
		a <- NA;
		n <- 0;
		if (no > 0) {n <- sum(depth[k])};
		if( n >= MIN ) a <- sum(val[k])/n;
		avg <- c(avg, a);
	}
	return( list(pos=pos, avg=avg) );
}

#https://academic.oup.com/bioinformatics/article/36/7/2150/5631910
calcBlockLoess=function(p1,p2,pos,bin.width=1000,min.depth=10,DEG=2,doSE=T) {
    opt.span <- function(model, criterion = c("aicc", "gcv"), span.range = c(0.05, 0.95)) {
	    as.crit <- function(x) {
	        span   <- x$pars$span;
	        traceL <- x$trace.hat;
	        sigma2 <- sum(x$residuals^2)/(x$n - 1);
	        aicc   <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL - 2);
	        gcv    <- x$n * sigma2/(x$n - traceL)^2;
	        result <- list(span = span, aicc = aicc, gcv = gcv);
	        return(result);
	    }
	    criterion <- match.arg(criterion);
	    fn <- function(span) {
	        mod <- update(model, span = span);
	        as.crit(mod)[[criterion]];
	    }
	    result <- optimize(fn, span.range);
	    return(list(span = result$minimum, criterion = result$objective));
	}
	# return

    #pmap binning 
    val=p1
    chr=max(pos,na.rm=T)+bin.width
    depth=p1+p2
    MIN=min.depth

    
    block <- cal_block_meta(pos, val,chr, bin.width, depth,MIN) 
    xx0 <- as.numeric(block$pos);
 #   print(xx0)
    yy <- as.numeric(block$avg);
    #print(x0)
    jdx <- which(!is.na(yy));
#    print(jdx)
 #   print(xx0[jdx])
    #print(jdx)
    fit0  <- loess(yy[jdx]~xx0[jdx], degree=DEG);

    span1 <- opt.span(fit0, criterion="aicc")$span;
    print(span1)
 #   print(span1)
   fit1  <- loess(yy[jdx]~xx0[jdx], degree=DEG, span=span1);
    plo <- predict(fit1, xx0, se=doSE);
    if(doSE){
        value <- plo$fit;
        return(list(x0=xx0, value=value, se=plo$se.fit))
    } else {
        value=plo
        return(list(x0=xx0, value=value))
        
    }
}




#generate progeny given a reference LD (pairwise pearson R between markers) matrix 
generateSegs=function(block.size, sample.size, LD.chol, threads, cutoff) { 
    # draw from multivariate normal
    holdme=matrix(0,block.size,dim(LD.chol)[1])
    system.time({
    rmvn(block.size, mu=rep(0,dim(LD.chol)[1]),sigma=LD.chol, isChol=T,ncores=threads,A=holdme)
    })
    # could cut af .5 (qnorm(.5)=0) but this generalizes to cases where allele frequencies 
    # could deviate from that 
    sim.x=apply(holdme, 1, function(y) ifelse(y<cutoff,0,1))
    sx=new("XSnpMatrix", sim.x, diploid=FALSE)
    colnames(sx)
     
    while(ncol(sx)<sample.size){
        rmvn(block.size, mu=rep(0,dim(LD.chol)[1]),sigma=LD.chol, isChol=T,ncores=threads,A=holdme)
        # could cut af .5 (qnorm(.5)=0) but this generalizes to cases where allele frequencies 
        # could deviate from that 
        sim.x=apply(holdme, 1, function(y) ifelse(y<cutoff,0,1))
        sx2=new("XSnpMatrix", sim.x, diploid=FALSE)
        current.indv=colnames(sx)[ncol(sx)]
        print(current.indv)
        l.indv=as.numeric(current.indv)+1
        colnames(sx2)=as.character(seq(from=l.indv, length.out=block.size))
        sx=cbind(sx,sx2)
    }

    return(sx)
}


#Calculate effect size given ratio of alleles
#K is selection threshold
#OR is observed P2/P1
# calculate x given OR and K
calcX.R = function (K, OR ) {
    T=qnorm(K, lower.tail = FALSE)
    fx=function(x,T,OR) { exp((pnorm(T, -x/2,sd=1, lower.tail=FALSE,log.p=T)
                              -pnorm(T,  x/2,sd=1, lower.tail=FALSE,log.p=T)))-OR }
    abs(uniroot ( 
        fx,
        c(-50,50), 
        T=T,
        OR=OR)$root)
}

# Leonid's clever approximation
calcX.LK = function(K,OR) {
    T=qnorm(K, lower.tail = FALSE)
    i=dnorm(T)/K
    (log(OR)/i)
}

#VE=x^2/4
#kvals = c(.1, .05, .01, .001, .0001, .00001, .000001)


#https://academic.oup.com/bioinformatics/article/36/7/2150/5631910#201693880
#read in gmaps
gmaps=readRDS('/data/xQTL/yeast_gmaps.RDS')

#Lg=sum(sapply(gmaps[['A']], function(x) max(x$map)))
calcEffectiveNumberOfTests=function(Lg=4900,Lp=1.2e4, el=8.4, nchr=16){
    #Lg=4900
    #Lp=sum(sapply(gmaps[['A']], function(x) max(x$ppos)))/1e3
    #Lp=1.2e4
    rM=Lp/Lg # ~3.5kb per cM
    #el=8.4
    eff.tests=16+Lp/(rM*el)
    return(eff.tests)
}

eff.tests=calcEffectiveNumberOfTests()



#reference genotypes for BYxRM bulk eQTL
g42.RDS='/data/single_cell_eQTL/yeast/Bulk/data/gdata42.RDS'

gdata42=readRDS(g42.RDS)

tg=data.frame(t(gdata42))
chr=tstrsplit(rownames(tg), ':')[[1]]
g.chr=split(tg, chr)
ppos=split(colnames(gdata42), chr)


# assuming no noteable long-range LD
cc='chrVII'
print(cc)

# lookup genetic map pos 
#more bits for handling genetic map
pname7=ppos[['chrVII']]
pname7=gsub(':|/', '_', pname7) 
mvec=pmatch(pname7, gmaps[['A']][['chrVII']]$marker)
#jfc
agmap= gmaps[['A']][['chrVII']][round(approxfun(mvec)(1:length(mvec))),]
#agmap$marker=ppos[['chrVII']]


#for sims 

x=data.matrix(g.chr[[cc]])
   # calculate LD as r
LD=crossprod(scale(t(x)))/ncol(x)
    #eff.tests[[cc]]=getMeff_Li_and_Ji(LD)

     #for(cc in unique(chr)) {
       # force it to be positive definite
LD=nearPD(LD)
    # cholesky decomposition
LD.chol=chol(LD$mat)

af=apply(x,1, function(y) sum(y==-1)/nrow(gdata42))
cutoff=as.vector(qnorm(af))

# read in data for 42k snps genotype data 
sample.size=5e5
for(i in 1:19) {
    sx= generateSegs(block.size=1e4, sample.size=sample.size, LD.chol, 70, cutoff)
    saveRDS(sx, file = paste0('/data/xQTL/chrVII_500k', i,'.RDS'))
}


sx=readRDS(paste0('/data/xQTL/chrVII_500k.RDS'))
#rand.segs[[cc]]=apply(rs, 1, function(y) ifelse(y<cutoff,0,1))
    #image.plot(as.matrix(LD$mat))
    #saveRDS(LD, file='/data/single_cell_eQTL/yeast/Bulk/data/LD.RDS')
    #LD.chol=cholesky(LD, parallel=T)
    #rand.segs.discrete
#}
##rs=lapply(rand.segs, function(x) Matrix(x,sparse=T))
#saveRDS(rs, file='/data/single_cell_eQTL/yeast/Bulk/data/simSegs.RDS')

#simulate marker effects from multivariate normal
#sim.effects=MASS::mvrnorm(1, mu=rep(0,3629), Sigma=(.05/3629)*as.matrix(LD$mat)+(.95/3629)*diag(3629))

simulatePhenotypesForRefPop=function(sx, h2=.1, nQTL=1) {

    nmarker=nrow(sx)
    nsample=ncol(sx)
    #h2=.1
    #for mvnorm
    #nadditive = nsample
    #a.eff  = rnorm(nsample) # rep(0,nsample)
    a.eff=rep(0,nsample)
    nadditive=nQTL
    add.qtl.ind  = sort(sample(nmarker, nadditive))
    add.qtl.sign = sample(ifelse(runif(nadditive)>.5,1,-1),replace=T)
    for(i in 1:nadditive){ a.eff=a.eff+add.qtl.sign[i]*as.numeric(sx[add.qtl.ind[i],]) }
    a.eff=scale(a.eff)

    g=sqrt(h2)*a.eff #+ sqrt(H2-h2)*aa.eff
    y=g+rnorm(nsample,mean=0,sd=sqrt((1-h2)/h2*var(g)))
    y=scale(as.vector(y))

    attr(y, "add.qtl")=add.qtl.ind*add.qtl.sign

    return(y)
}

#hugex=matrix(as.numeric(t(sx@.Data)), dim(sx)[2], dim(sx)[1])
#lstat=Rfast::correls(y, hugex)

#summary(lm(y~as.numeric(sx[add.qtl.ind[1],])))
#summary(lm(y~as.numeric(sx[add.qtl.ind[2],])))

#plot(lstat[,4])
#get.LOD.by.COR = function(n.pheno, pheno, gdata) {
    # Lynch and Walsh p. 454
#    return( (-n.pheno*log(1-cor(pheno, gdata, use='pairwise.complete.obs')^2))/(2*log(10)) ) }

#LOD=get.LOD.by.COR(nsample, y,t(x))
#pearsonr=as.numeric(cor(y,t(x))[1,])
#plot(pearsonr)

eff.length=unlist(eff.tests*(sapply(gmaps[['A']], function(x) max(x$ppos))/12e6)[cc])




generateCounts=function(y, sx, sel.frac, depth, lower.tail=F) {
    #$sel.frac=.10
    #get positive tail
    if(lower.tail==F) {
        sel.indv=which(y> quantile(y,1-sel.frac))
    } else{
        sel.indv=which(y< quantile(y,sel.frac))
    }
    sel.indv.x=sx[,sel.indv]

    #assuming 0/1 coding 
    sel.indv.af=Rfast::rowsums(sel.indv.x)/ncol(sel.indv.x)
    #plot(sel.indv.af)

    r=rbinom(n=length(sel.indv.af),size=depth, prob=sel.indv.af)
    a=rbinom(n=length(sel.indv.af),size=depth, prob=1-sel.indv.af)
    return(   list(expected=sel.indv.af, p1=r, p2=a))
}



doBlockLoess=function(p1, p2, pos, bin.width=500) {
    loess.results=calcBlockLoess(p1,p2,pos, bin.width=bin.width, DEG=2)
    imputed.afd.delta=approxfun(loess.results$x0,loess.results$value)(pos)
    return(imputed.afd.delta) 
}



physical.position = as.numeric((tstrsplit(ppos[[cc]], ':|_'))[[2]])
 genetic.position = agmap

y = simulatePhenotypesForRefPop(sx, h2=.025, nQTL=2)

#get simulated high tail
sel.high=0.1
sel.low=0.1

high.tail = generateCounts(y,sx, sel.high, depth=25)
low.tail  = generateCounts(y,sx, sel.low, 25, lower.tail=T)

n1=sample.size*sel.high
n2=sample.size*sel.low


results.data=data.frame(
                      physical.position=physical.position,
                      genetic.postion=genetic.position,
                      expected.af.high=high.tail$expected,
                      expected.af.low=low.tail$expected,
                      p1.high=high.tail$p1,
                      p2.high=high.tail$p2,
                      p1.low=low.tail$p1,
                      p2.low=low.tail$p2)

results.data$afd.high=doBlockLoess(p1=results.data$p1.high, 
                                   p2=results.data$p2.high,
                                   pos=results.data$physical.position)
results.data$afd.low=doBlockLoess(p1=results.data$p1.low, 
                                  p2=results.data$p2.low,
                                  pos=results.data$physical.position)


effective.sample.n.two.tailed=1/((n1+n2)/(n1*n2))
effective.sample.n.low=n2
effective.sample.n.high=n1

#sample.n.tail=sample.size*sel.frac
#print(paste('actual sample size= ', sample.n.tail, ' |  estimated sample size', sum(r+a)/eff.length ))
depth.sample.size.low=sum(results.data$p1.low+results.data$p2.low)/eff.length 
depth.sample.size.high=sum(results.data$p1.high+results.data$p2.high)/eff.length 

n.low=min(effective.sample.n.low,depth.sample.size.low)
n.high=min(effective.sample.n.high, depth.sample.size.high)

n.both=1/((n.low+n.high)/(n.low*n.high))

results.data$afd.high.se=sqrt((results.data$afd.high*(1-results.data$afd.high))/n.high) 
results.data$afd.low.se=sqrt((results.data$afd.low*(1-results.data$afd.low))/n.low) 

#plots for the two tails 
par(mfrow=c(2,1))

with(results.data, plot(physical.position, p1.high/(p1.high+p2.high), ylim=c(0,1),main='high tail'))
with(results.data, points(physical.position,expected.af.high, col='red',lwd=2, type='b'))
with(results.data, abline(v=physical.position[abs(attr(y,"add.qtl"))]))
with(results.data, points(physical.position, afd.high, col='purple', type='l', lwd=3))

with(results.data, points(physical.position, afd.high+1.96* afd.high.se, col='purple', type='l', lwd=2, lty=2))
with(results.data, points(physical.position, afd.high-1.96* afd.high.se, col='purple', type='l', lwd=2, lty=2))


with(results.data, plot(physical.position, p1.low/(p1.low+p2.low), ylim=c(0,1),main='low tail'))
with(results.data, points(physical.position,expected.af.low, col='red',lwd=2, type='b'))
with(results.data, abline(v=physical.position[abs(attr(y,"add.qtl"))]))
with(results.data, points(physical.position, afd.low, col='purple', type='l', lwd=3))

with(results.data, points(physical.position, afd.low+1.96* afd.low.se, col='purple', type='l', lwd=2, lty=2))
with(results.data, points(physical.position, afd.low-1.96* afd.low.se, col='purple', type='l', lwd=2, lty=2))

# plots for the contrast 
#x11()
with(results.data, plot(physical.position, (p1.high/(p1.high+p2.high))-(p1.low/(p1.low+p2.low)) )  )
with(results.data, points(physical.position,expected.af.high-expected.af.low, col='red',lwd=2, type='b'))
with(results.data, abline(v=physical.position[abs(attr(y,"add.qtl"))]))
with(results.data, points(physical.position, afd.high-afd.low, col='purple', type='l', lwd=3))
with(results.data, points(physical.position, (afd.high-afd.low)+1.96*sqrt(afd.high.se^2+afd.low.se^2), col='purple', type='l', lwd=2, lty=2))
with(results.data, points(physical.position, (afd.high-afd.low)-1.96*sqrt(afd.high.se^2+afd.low.se^2), col='purple', type='l', lwd=2, lty=2))


results.data$z=with(results.data, (afd.high-afd.low)/sqrt(afd.high.se^2+afd.low.se^2))
results.data$p=with(results.data, 2*pnorm(abs(z), lower.tail=F) ) 
results.data$LOD=PvalToLOD(results.data$p)

pracma::findpeaks(-log10(results.data$p), nups=3, threshold=-log10(.05/600), sortstr=T, npeaks=2)

with(results.data, plot(physical.position, z))
with(results.data, plot(physical.position, -log10(p)) ) #LOD)) #-log10(p)))

with(results.data, plot(physical.position, LOD)) #-log10(p)))
with(results.data, abline(v=physical.position[abs(attr(y,"add.qtl"))]))

with(results.data, plot(physical.position, pchisq(-log10(results.data$p)))

with(results.data, abline(v=physical.position[abs(attr(y,"add.qtl"))]))
abline(h=-log10(.05/600))


peak.index=with(results.data, which.max(afd.high-afd.low))
peak.afd=with(results.data,(afd.high-afd.low)[peak.index] ) 
peak.afd.delta=with(results.data, (afd.high-afd.low)[peak.index]-1.96*sqrt(afd.high.se[peak.index]^2+afd.low.se[peak.index]^2))




























which.max()


RtoLOD(n.both, tanh(results.data$z))

sample.n=min(sample.n.tail, sum(r+a)/eff.length )
se=sqrt((afd*(1-afd))/sample.n) 



yy=(r/(r+a))
xx=1:length(r)

loc=as.numeric((tstrsplit(ppos[[cc]], ':|_'))[[2]]) #tstrsplit(pppos[[chr]], ':|_')[,2]
xdf=data.frame(ppos=loc, gpos=agmap, p1=r, p2=a, d=r+a)

#agmap$ppos
#gloc=agmap$map

plot(xdf$ppos, yy, ylim=c(0,1))
points(xdf$ppos, sel.indv.af, col='red',type='l', lwd=4)
abline(v=xdf$ppos[add.qtl.ind])

phys.loess=calcBlockLoess(xdf$p1, xdf$p2, xdf$ppos, bin.width=500, DEG=2)
points(phys.loess$x0, phys.loess$value, col='purple', type='l', lwd=3)
#SE from loess 
points(phys.loess$x0, phys.loess$value+1.96* phys.loess$se, col='purple', type='l', lwd=2, lty=2)
points(phys.loess$x0, phys.loess$value-1.96* phys.loess$se, col='purple', type='l', lwd=2,lty=2)

afd=(phys.loess$value)
sample.n.tail=sample.size*sel.frac
print(paste('actual sample size= ', sample.n.tail, ' |  estimated sample size', sum(r+a)/eff.length ))

sample.n=min(sample.n.tail, sum(r+a)/eff.length )
se=sqrt((afd*(1-afd))/sample.n) #sample.n) #median(r+a)) #median(r+a))  #(sample.n/10))
points(phys.loess$x0, phys.loess$value+1.96* se, col='green', type='l', lwd=2, lty=2)
points(phys.loess$x0, phys.loess$value-1.96* se, col='green', type='l', lwd=2,lty=2)

imp.raw=approxfun(xdf$ppos,sel.indv.af)(phys.loess$x0)
#points(phys.loess$x0, imp.raw, col='pink', type='l', lwd=2, lty=2)

checkSE=(imp.raw>(phys.loess$value-1.96* se) & imp.raw<(phys.loess$value+1.96* se))
sum(checkSE, na.rm=T)/sum(!is.na(checkSE))

OR=(sel.indv.af[add.qtl.ind])/(1-sel.indv.af[add.qtl.ind])
#OR=ifelse(OR<1,1/OR,OR)
#calcX.R(sel.frac, OR)

B_012=sapply(OR,function(x) calcX.R(sel.frac,x))/2
B_012
VE_012=(B_012^2)
VE_012

B_012=calcX.LK(sel.frac,OR)
B_012
VE_012=(B_012^2)/4
VE_012


# individs in pools
N=sample.n.tail
# resolution

#bin size for discrete model 
res=100

#length of a cm in base pairs 
cM=2800

pr= res/100.0/cM

dfs=cbind(xdf$ppos,xdf$p1,xdf$p2) #agmap$map, r, a)
sfactor=seq(.5,100, length.out=20)
ll=rep(0,20)
intm=0

for(ss in sfactor) {

    intm=intm+1
nout   = round(max(dfs[,1])/res)+1
bmatch = dfs[,1]/res
bmatch.int =findInterval(bmatch, 1:nout, all.inside=T)

means  = rep(0, nout)
means[bmatch.int]=dfs[,2]                # np  (i.e. 1000 * .5)
#means[means==0]=NA
   
counts = rep(0, nout)
counts[bmatch.int]= dfs[,2] + dfs[,3]      # n   (i.e. 1000)
   
p =rep(0,nout)
p[bmatch.int] = dfs[,2]/(dfs[,2]+dfs[,3])   # p   (i.e. .5)
 #  0.6013
variances = counts/ss*(p*(1-p)) #counts*(p *( 1-p))  # np(1-p)  (i.e  1000*.5(1-.5)
variances[variances==0]=NA

#means=p


y=means
y_var=variances
d=counts
TT=length(y)
k1=kalman(y,y_var,d,TT,N,pr)
#points(xdf$ppos, (k1$mu_pstr)[bmatch.int]/5000, col='blue', lwd=2)
    ll[intm]=k1$logLik

#L1=calcLOD( k1$mu_pstr, k1$V_pstr, TT, N)
plot(loc, yy, ylim=c(0,1))
points(loc, sel.indv.af, col='red',type='l', lwd=4)
points(xdf$ppos, (k1$mu_pstr)[bmatch.int]/N , col='blue', lwd=2)
points(xdf$ppos, (k1$mu_pstr)[bmatch.int]/N -1.96*sqrt(k1$V_pstr)[bmatch.int]/N,  col='lightblue', type='l', lwd=2, lty=2)
points(xdf$ppos, (k1$mu_pstr)[bmatch.int]/N +1.96*sqrt(k1$V_pstr)[bmatch.int]/N  , col='lightblue', type='l', lwd=2,lty=2)

}

# for gmap
#res=.001
#pr=max(xdf$gpos.map)/nrow(xdf)

#sfactor=seq(.0001,4, length.out=100) 


#for genetic map

res=.05
#ipr=6.364e-05
sfactor=seq(1e-3,1e-1, length.out=100) 

#sfactor=seq(.01,.8, length.out=100) 
ll=rep(0,100)
intm=0
sfactor=0.0006367
#sfactor=1
for(pr in sfactor) {
    intm=intm+1
    #push this into gmap space 
    #dfs=cbind(xdf$ppos,xdf$p1,xdf$p2) #agmap$map, r, a)
    dfs=cbind(xdf$gpos.map,xdf$p1,xdf$p2) #agmap$map, r, a)
    
    nout   = round(max(dfs[,1])/res)+1
    bmatch = dfs[,1]/res
    bmatch.int =findInterval(bmatch, 1:nout, all.inside=T)

    means  = rep(0, nout)
    means[bmatch.int]=dfs[,2]                # np  (i.e. 1000 * .5)
    #means[means==0]=NA
       
    counts = rep(0, nout)
    counts[bmatch.int]= dfs[,2] + dfs[,3]      # n   (i.e. 1000)
       
    p =rep(0,nout)
    p[bmatch.int] = dfs[,2]/(dfs[,2]+dfs[,3])   # p   (i.e. .5)
     #  0.6013
    variances = counts*(p*(1-p)) #counts*(p *( 1-p))  # np(1-p)  (i.e  1000*.5(1-.5)
    variances[variances==0]=NA

    #means=p
   

    y=means
    y_var=variances
    d=counts
    TT=length(y)
    k1=kalman(y,y_var,d,TT,N,pr)
    #points(xdf$ppos, (k1$mu_pstr)[bmatch.int]/5000, col='blue', lwd=2)

    #L1=calcLOD( k1$mu_pstr, k1$V_pstr, TT, N)
    plot(loc, yy, ylim=c(0,1))
    points(loc, sel.indv.af, col='red',type='l', lwd=4)
    points(xdf$ppos, (k1$mu_pstr)[bmatch.int]/N , col='blue', lwd=2)
    points(xdf$ppos, (k1$mu_pstr)[bmatch.int]/N -1.96*sqrt(k1$V_pstr)[bmatch.int]/N,  col='lightblue', type='l', lwd=2, lty=2)
    points(xdf$ppos, (k1$mu_pstr)[bmatch.int]/N +1.96*sqrt(k1$V_pstr)[bmatch.int]/N  , col='lightblue', type='l', lwd=2,lty=2)

    ll[intm]=k1$logLik
    
}

#KFAS



library(dlm)




plot(k1$mu_pstr)






#calcX.R(sel.frac,OR)

z=(afd-.5)/se
yimp=approxfun(phys.loess$x0,z)(loc)

approxp=2*pnorm(abs(yimp), lower.tail=F)

plot(-log10(approxp))
abline(h=-log10(.05/600))
abline(v=add.qtl.ind)

#impute stat at typed snps 
#(loc,yimp)

baf=replicate(100, {
    print('yo')
    rs=rbinom(n=length(sel.indv.af),size=r+a, prob=r/(r+a)) 
    as=rbinom(n=length(sel.indv.af),size=r+a, prob=(1-r/(r+a))) 
    ls=calcBlockLoess(rs,as,loc,bin.width=1000,doSE=F)
    return(ls$value)
})
baf.se=apply(baf,1,sd, na.rm=T)
median(baf.se, na.rm=T)



R=as.matrix(LD$mat)

test=susie_rss(yimp,R)




















DEG=2

m0=loess(yy~loc, degree=2);
m0span1 <- opt.span(m0, criterion="aicc")$span;
m1  <- loess(yy~loc, degree=DEG, span=m0span1);

m0g=loess(yy~gloc, degree=2)
m0span1g <- opt.span(m0g, criterion="aicc")$span;
m1g  <- loess(yy~gloc, degree=DEG, span=m0span1g);

points(loc, predict(m1), col='blue', type='l', lwd=2)
points(loc, predict(m1g), col='green', type='l', lwd=2)





phys.loess=calcBlockLoess(r,a,loc,bin.width=1000)












#gmap binning 
val=r/(r+a)
chr=max(gloc) #1070000 #68842
depth=r+a
MIN=10

block <- cal_block_meta(gloc, val,chr, .5, depth,10) 
x0 <- as.numeric(block$pos);
y <- as.numeric(block$avg);
jdx <- which(!is.na(y));
fits <- list();
DEG=2
fit0  <- loess(y[jdx]~x0[jdx], degree=DEG);
span1 <- opt.span(fit0, criterion="aicc")$span;
fit1  <- loess(y[jdx]~x0[jdx], degree=DEG, span=span1);
plo <- predict(fit1, x0, se=TRUE);
value <- plo$fit;

#points(x0, value,col='purple', type='l', lwd=2)

points(approxfun(gloc, loc)(x0), approxfun(x0,value)(x0),col='purple', type='l', lwd=2 )





points(x0, value,col='pink', type='l', lwd=2)

































































xx=agmap$map
b=gam(yy~s(xx,k=5), drop.intercept=T)

G=gam(yy~s(xx), fit=FALSE, drop.intercept=T)

Cv=LD.chol
w=solve(t(Cv))
w=matrix(w,dim(w)[1], dim(w)[2]) #matrix(w)
mgfit <- magic(G$y, G$X, G$sp, G$S, G$off, rank=G$rank, C=G$C, w=w)
mg.stuff <- magic.post.proc(G$X, mgfit, w)

b2 <- b ## copy b
b2$edf <- mg.stuff$edf
b2$Vp <- mg.stuff$Vb
b2$coefficients <- mgfit$b 






w=matrix(w,dim(w)[1], dim(w)[2]) #matrix(w)
mgfit <- magic(G$y, G$X, G$sp, G$S, G$off, rank=G$rank, C=G$C, w=w)


se.ds.af=sqrt((ds.af*(1-ds.af)))





depth=100
rsim=rbinom(n=length(sel.indv.af),size=100, prob=sel.indv.af)
asim=rbinom(n=length(sel.indv.af),size=100, prob=1-sel.indv.af)

tsim=rsim+asim

pR=dbinom(rsim, tsim, .01)
pA=dbinom(asim, tsim, .01)

map2rf = function(map, tol=1e-12) {
    rf = mf.h(diff(map))
    rf[rf < tol] = tol # don't let values be too small
    rf
}

g.rf=map2rf(agmap$map)
mTmatH=function(r) {   matrix(c((1-r),r,r,(1-r)),2,2) }

#convert to transmission matrices 
#generalize
#tMats=lapply(g.rf, mTmat)
tMats=lapply(g.rf, mTmatH)

gg
code.dir='/data/single_cell_eQTL/yeast/code/'
source(paste0(code.dir, 'HMM_fxs.R'))
source(paste0(code.dir, 'runHMM_custom.R'))

eMats=rbind(pR, pA)
ess=colSums(eMats)
eMats[,ess<1e-12]=1

f=calcForward(eMats,tMats)
b=calcBackward(eMats,tMats)
posteriorProb=calcPosterior(f,b)















    g=data.matrix(g.chr[[cc]])
    # calculate LD as r
    LD=crossprod(scale(t(g)))/ncol(g)
    LDaug=LD+diag(.1, nrow(g))
    
    # force it to be positive definite
    #LD=nearPD(LD)









#Xinv=chol2inv(LD.chol)
Xinv=chol2inv(LD.chol) #+diag(.2, nrow(x)))
SigIt=matrix(LD$mat@x,LD$mat@Dim[1], LD$mat@Dim[1]) #+diag(.2, nrow(x))

n=rep(NA, nrow(x))
for(i in 1:nrow(x)) {
    n[i]=unlist(SigIt[i,]%*%Xinv%*%ds.af)
}


Xinv=solve(as.matrix(LD.chol)+diag(1e-2, nrow(x)))


Xinv=chol2inv(LD.chol+diag(1e2, nrow(x)))

n=Xinv%*%(ds.af)

r=cor(as.numeric(n),t(x))
plot(r[1,])
