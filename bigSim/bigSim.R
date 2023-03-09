#devtools::install()

library(xQTLStats)
jitterGmapVector=function(themap, amount=1e-6) {
    for (i in 1:length(themap)) {
         n <- length(themap[[i]])
         themap[[i]] <- themap[[i]] + c(0, cumsum(rep(amount, n - 1)))
    }
    return(themap)
}
getGmapPositions=function(vcf.cross, gmap, uchr) {
    #get physical position, split by chromosome
    p.by.chr=split(vcfR::getPOS(vcf.cross),vcfR::getCHROM(vcf.cross))
    #keep things sorted (yay yeast chr names with roman numerals)
    p.by.chr=p.by.chr[uchr]

   #where to put the variant sites, impute onto gmap
    imputed.positions=mapply( 
           function(x, y){
                approxfun(y$ppos, y$map, rule=2)(x)
            },
            x=p.by.chr, y=gmap,
            SIMPLIFY=F)

}


#some more complex simulations here 

#optional for doing simulations 
#remotes::install_url("https://cran.microsoft.com/snapshot/2017-09-17/src/contrib/Meiosis_1.0.2.tar.gz")

data(crosses.to.parents)

#genetic maps 
gmaps=readRDS(system.file('reference', 'yeast_gmaps.RDS', package='xQTLStats'))

#reference vcf (assume ploidy=1)
ref.vcf=system.file('reference', 'parents_w_svar_sorted.vcf.gz', package='xQTLStats')

#assuming input vcf with haploid parents and haploid specific variant calls 
# [1] "273614xa" "BYa"      "CBS2888a" "CLIB219x" "CLIB413a" "I14a"    
# [7] "M22"      "PW5a"     "RMx"      "Y10x"     "YJM145x"  "YJM454a" 
#[13] "YJM978x"  "YJM981x"  "YPS1009x" "YPS163a"

#specify sample size
sample.size=1e5

#specify cross (see crosses.to.parents)
cross='A'

#get the two parent names 
p1.name=crosses.to.parents[[cross]][1]
p2.name=crosses.to.parents[[cross]][2]

#extract parental genotypes at variant sites between those parents from a reference vcf
vcf.cross=getCrossVCF(ref.vcf,p1.name, p2.name)

gmap=gmaps[[cross]]

uchr=paste0('chr', as.roman(1:16))
eg.GT=vcfR::extract.gt(vcf.cross)
imputed.positions=jitterGmapVector(getGmapPositions(vcf.cross, gmap, uchr=paste0('chr', as.roman(1:16))))


#how many chromosomes
n_chr=length(imputed.positions[5])
    #total gmap size
L=round(sapply(gmap, function(x) max(x$map)))[5]
n_chr=c(n_chr,1)
L=c(L,51)
#setup for Meiosis package
xoparam = Meiosis::create_xoparam(L,obligate_chiasma=T) 
    
    #parental haplotypes
    p1.ind=split(as.numeric(eg.GT[,1]), vcfR::getCHROM(vcf.cross))
    p1.ind=p1.ind[uchr]
    p2.ind=split(as.numeric(eg.GT[,2]), vcfR::getCHROM(vcf.cross))
    p2.ind=p2.ind[uchr]

    ind=list(list('chrVI'=p1.ind[[5]], 'chrL'=rep(0,51)),
             list('chrVI'=p2.ind[[5]],'chrL'=rep(1,51)))
    nsegs=1e5 #sample.size
    #total number of meiosis to simulate
    psize=round(nsegs/2)
    #do sim
    imp.aug=list('chrVI'=imputed.positions[[5]], 'chrL'=seq(0,51, length.out=51))
    p_geno = replicate(psize, Meiosis::cross_geno(father = ind, mother = ind,
                                                  positions = imp.aug, xoparam=xoparam), 
                       simplify=F)

    ngenerations=60
    if(ngenerations>2){
        for(ngen in 3:ngenerations) {
            print(ngen)

            p_geno=replicate(psize, {
                    n1=sample.int(psize,1)
                    n2=sample.int(psize,1)
                    #paternal and maternal are interchangeable here but might as well just choose one for each sampling
                    nind=list(paternal=p_geno[[n1]]$paternal,
                              maternal=p_geno[[n2]]$maternal)
                    return(Meiosis::cross_geno(father = nind, mother = nind, 
                                               positions = imp.aug,  xoparam=xoparam))
                           },simplify=F)

        }
    } 


    pp=sapply(p_geno,function(x) x$paternal[[1]])
    pm=sapply(p_geno,function(x) x$maternal[[1]])
    pg=cbind(pp,pm)      
   # pgx2=2*tcrossprod(pg)/ncol(pg)

    pgx=crossprod(scale(t(pg)))/ncol(pg)
    rm(pp)
    rm(pm)

    y=simBayesC(pg, target.size=.01, h2=.8/16)
    


    XtX=tcrossprod(p_geno)
    XtXinv=MASS::ginv(XtX)
   
    Xty=(p_geno %*%y)
    seB=sqrt((1-.8/16)*diag(XtXinv))


    r=(t(y) %*% scale(t(pg)))/(nrow(pg)-1)
    plot(r[1,])
    abline(v=attr(y,'causal.indices'))



    B=XtXinv %*% Xty
    simB=rep(0, nrow(p_geno))
    simB[attr(y, 'causal.indices')]=attr(y, "causal.betas")



#segments(B,B+1.96*seB,B,B-1.96*seB)
plot(simB*.05,B)
plot(B)
segments(1:length(B),B+1.96*seB,1:length(B),B-1.96*seB)
points(simB*.05, col='red')

test=cv.glmnet(x=t(p_geno), y=y,alpha=0, trace.it=T)


    p_geno=do.call('rbind', p_geno)

geno.matrix=simHaploidSegsFN(vcf.cross[which(vcf.cross@fix[,1]=='chrI'|vcf.cross@fix[,1]=='chrVI'),] ,
                 gmaps[[cross]][c('chrI', 'chrVI')],
                 sample.size, ngenerations=2,
                 uchr = c('chrI', 'chVI'))



g=matrix(as.numeric(geno.matrix@.Data), nrow(geno.matrix), ncol(geno.matrix))
tg=t(g)
dtg=duplicated(g)
tg=tg[,-which(dtg)]
g.s=(g)

GGs=crossprod(g.s)

efitness0=runif(ncol(g.s)) # rep(1,ncol(g))
eCount0=rowSums(g.s%*%efitness0)
eFreq0=eCount0/sum(efitness0)
efitness0.norm=efitness0/sum(efitness0)
depth=1e6
ytest=rbinom(n=length(eFreq0),size=depth, prob=eFreq0)
test=RcppML::nnls(GGs,matrix(crossprod(g.s,ytest)), fast_nnls=T)
 plot(test[,1],efitness0)

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

