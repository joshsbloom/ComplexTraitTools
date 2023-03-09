library(lme4)

#some observed data
p=readRDS('/data/yeast/Colonies/2022-04-29/ResultsTable.RDS')

ps=p[p$s.perimeter<600,]
cond.split=split(ps, ps$condition1)

test=cond.split[['9']]

#bm=lm(test$s.radius.mean~test$file + test$StrainName)

lbm=lmer(test$s.radius.mean~test$file +(1|test$StrainName))
sd(residuals(lbm))
#min is about 20
# max is about 80
# mean is about 40


#now the GWAS stuff ----------------------------------------------------
library(BEDMatrix)
library(rrBLUP)
library(NAM)
library(regress)

#strongly recommend linking your R to openblas or MKL

#load Joseph's prefiltered GWAS genotype data and pre-processed phenotypes 
#https://www.nature.com/articles/s41586-018-0030-5
path.to.gwas.data='/data/yeast/1002genomes/1011GWASMatrix/'

setwd(path.to.gwas.data)

#read marker data (BED/BIM/BAM plink formatted)
x=BEDMatrix('1011GWAS_matrix.bed')
g=as.matrix(x)
mnames=read.delim('1011GWAS_matrix.bim', header=F, sep='\t', stringsAsFactors=F)


#remove blank and missing strains -----------------
g.subset=g

#morejunk=which(rownames(pheno)=='NA')
#-------------------------------------------------

Amat=A.mat(g.subset-1, return.imputed=T)
A=Amat$A
g.subset.imp=Amat$imputed

g.s=scale(g.subset.imp)

coi=5

#calculate pariwise correlation between strains
cor.mat=tcrossprod(g.s)
cor.mat=cor.mat/(ncol(g.s)-1)

target.chr=g.subset.imp[,mnames[,1]==coi]
target.chr.s=g.s[,mnames[,1]==coi]

af=colSums(target.chr+1)/(2*nrow(target.chr))

lococor.mat=tcrossprod(g.subset.imp[,mnames[,1]!=coi]) #[,1:82869])
lococor.mat=lococor.mat/(ncol(g.subset.imp[,mnames[,1]!=coi])-1) #(82869-1)

#calculate eigen value decomposition of relatedness matrix 
eigout.wg=eigen(cor.mat, symmetric=T)
eigout.loco=eigen(lococor.mat, symmetric=T)


#H2=.05
pheno.sd=1.4

nsim=250

set.seed(150)
#pick 100 markers
causals=sort(sample(1:ncol(target.chr), 150))

gwas.sim=list()

betas=c(.1,.25,.5,.75,1,2,4)

for(beta in betas) {
    print(beta)
    #sim.ys=40+target.chr.s[,causals]*beta

    #sim.ys=sim.ys+error #10) #sqrt((1-H2)/H2*var(sim.y)))
    VarX=apply(target.chr[,causals],2,var)
    PVE=(beta^2*VarX)/((beta^2*VarX)+pheno.sd^2)

    #sqrt(PVE)
    #plot(af[causals], PVE)
    plot(af[causals], sqrt(PVE))

    #i=10
    #summary(lm(scale(sim.y[,i])~scale(target.chr[,causals[i]])))

    pmat=matrix(NA, nsim, length(causals))
    for( n in 1:nsim) {
       # print(n)
        error=rnorm(nrow(target.chr), mean=0, sd=pheno.sd)
        sim.y=40+target.chr[,causals]*beta
        sim.y=sim.y+error

        pvec=rep(NA, length(causals))
        #zvec=rep(NA, length(causals))
        for(i in 1:length(causals)) {
            summ=summary(lm(sim.y[,i]~target.chr[,causals[i]]+eigout.loco$vectors[,1:20]))
            pvec[i] = summ$coef[2,4]
            #zvec[i] = summ$coef[2,3]
        }
        pmat[n,]=pvec
    }
    powerR=colSums(pmat<1e-5)/nsim

    gwas.sim[[as.character(beta)]]=
        data.frame(causal=causals, af=af[causals], PVE=PVE, sqrtPVE=sqrt(PVE), beta=beta, sd=pheno.sd, power=powerR)
}

gwas.power=data.table::rbindlist(gwas.sim)

library(ggplot2)
ggplot(gwas.power, aes(x=sqrtPVE, y=power, color=factor(beta)))+geom_point()+theme_bw()+
    xlab('standarized QTL effect')+ xlim(0,1)


ggplot(gwas.power, aes(x=af, y=sqrtPVE))+facet_wrap(~beta)+geom_point()+theme_bw()+xlab('maf')+ ylab('standardized QTL effect')

saveRDS(gwas.power, file = '/data/yeast/rvas/data/gwaspower1.4.RDS')


gwas.power=readRDS(gwas.power, file = '/data/yeast/rvas/data/gwaspower2.75.RDS')

gwas.results=list()

for(n in colnames(pheno) ) { # colnames(pheno)[-1]){
    print(n)
    #implements gemma algorithm
    gwas.results[[n]]=
        test=gwas2(sim.y, target.chr, EIG=eigout.loco, fixed=F)
}

lm(sim.y~


rr=regress(sim.y~1, ~lococor.mat)
rresid=sim.y-rr$predicted

ypc=lm(sim.y~eigout.loco$vectors[,1:20])
gpc=residuals(lm(target.chr~eigout.loco$vectors[,1:20]))
t2=Rfast::regression(gpc, residuals(ypc)) #rresid) #residuals(ypc))

#plot(test$PolyTest$pval)

t2=Rfast::regression(target.chr, rresid) #residuals(ypc))

plot(-log10(t2[,2]))


y=g+rnorm(nrow(gdata),mean=0,sd=sqrt((1-H2)/H2*var(g)))
y=g+rnorm(nrow(gdata),mean=0,sd=sqrt((1-H2)/H2*var(g)))
y=scale(as.vector(y))











#you see what I'm doing here , just some ideas 
gw_GJ_norm=gwas2(pheno[,'GJ+LP']-pheno[,'GJ'], g.subset, EIG=eigout, fixed=F)
gwas.results[['GJ_norm']]=gw_GJ_norm


gw_SDM_norm=gwas2(pheno[,'SDM+LP']-pheno[,'SDM'], g.subset, EIG=eigout, fixed=F)
gwas.results[['SDM_norm']]=gw_GJ_norm

#ph4m6=gwas2(pheno[,'ph4']-pheno[,'ph6'], g.subset, EIG=eigout, fixed=F)
#gwas.results[['ph4m6']]=ph4m6

#ph8m6=gwas2(pheno[,'ph8']-pheno[,'ph6'], g.subset, EIG=eigout, fixed=F)
#gwas.results[['ph8m6']]=ph8m6

# fit a linear model of colony radius vs pH
#slopes=apply(pheno, 1, function(x) coef(lm(x~c(4,5,6,7,8))))
#phSlope=gwas2(slopes[2,], g.subset, EIG=eigout, fixed=F)

#gwas.results[['phSlope']]=phSlope

saveRDS(gwas.results, file =paste(base.path, 'GWASResults.RDS', sep =''))


 
