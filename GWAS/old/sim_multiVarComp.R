library(BEDMatrix)
library(regress)
library(Rfast)

path.to.gwas.data='/data/yeast/1002genomes/1011GWASMatrix/'
setwd(path.to.gwas.data)

x=BEDMatrix('1011GWAS_matrix.bed')
g=as.matrix(x)
mnames=read.delim('1011GWAS_matrix.bim', header=F, sep='\t', stringsAsFactors=F)

#calculate allele frequency
af=colSums(g, na.rm=T)/(2*apply(g,2,function(x) sum(!is.na(x)))) 

af.cutoff=.05

notna.cnt=apply(g, 2, function(x) sum(!is.na(x)))
notna.frq=notna.cnt/nrow(g)
notna.frq.cut=.8

#gvar=apply(g.s,2,var, na.rm=T)
bm=which(af<af.cutoff | notna.frq<notna.frq.cut)

#final marker set
g=g[,-bm]

#corresponding marker info
mnames=mnames[-bm,]

#replace NAs with 0s
g.0=g
g.0[is.na(g.0)]=0

GRM=tcrossprod(scale(g.0))/ncol(g.0)

aneu=readRDS('~/Desktop/aneu.rds')
#check that strain names match up
GRM=GRM[pmatch(rownames(aneu), rownames(GRM)), pmatch(rownames(aneu), rownames(GRM))]
g=g[pmatch(rownames(aneu), rownames(g)),]


ANEU.kernel=tcrossprod(scale(aneu))/ncol(aneu)

#simulate variance from each term
Aeff=.05
Geff=.6
Eeff=1-(Aeff+Geff)
#Total variance = Aeff+Geff+Eeff and should sum to 1

#simulate effects from the multivariate normal distribution
simy=rmvnorm(1,mu=rep(0,nrow(GRM)), sigma=Aeff*ANEU.kernel+Geff*GRM+Eeff*diag(nrow(GRM)))

#fit model, check if we can recover simulated effects
regress(simy[1,]~1, ~ANEU.kernel+GRM,pos=c(T,T,T), verbose=T)


colnames(GRM)=rownames(GRM)=1:length(simnb)

simy=rmvnorm(1,mu=rep(0,nrow(GRM)), sigma=.2*GRM+.8*diag(nrow(GRM)))
simnb=rnegbin(length(simy), exp(simy)*.5, theta=1/(.2^2))

pheno=data.frame(p=simnb,id=1:length(simnb))

regress(log(simnb+.5)~1,~GRM, verbose=T)

#test=glmmkin(p ~ 1, data=pheno, kins=GRM, id="id", family=poisson(link='log'),verbose=T)
test=glmmkin(p ~ 1, data=pheno, kins=GRM, id="id", family=quasipoisson(link='log'),verbose=T)

omega=test$theta[1]
#lambda=mean(test$fitted.values)
lambda=exp(test$coefficients[1]+0.5*test$theta[2])
test$theta[2]/(test$theta[2]+log(1+omega/lambda))



 test=glmmkin(p ~ lUMI + experiment, data=pheno, kins=GRM, id="id", family=quasipoisson(link='log'), verbose=T)

