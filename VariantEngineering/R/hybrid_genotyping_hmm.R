library(Matrix)
library(data.table)
library(magrittr)
library(multidplyr)
library(dplyr)

nthreads=35
input_file='/data/CRISPR_variant_engineering/guide_diploids/strain_list.RDS'
output_file='/data/CRISPR_variant_engineering/guide_diploids/ase_list_withHMM.RDS'

data=readRDS(input_file)
data$ase_list$contig=factor(data$ase_list$contig, levels=paste0('chr', as.roman(1:16)))

devtools::source_url("https://raw.githubusercontent.com/joshsbloom/single_cell_eQTL/master/elegans/code/HMM_fxs.R")

#input here isn't a sparse array across samples, so build a new f(x) to accomodate 
get_eMats=function(rcount, acount, error.rate=.005) {

    n=rcount+acount
    k=rcount
    # per UMI genotype-specific count error rate
    ee=error.rate

    # simpler to keep this as a vector for these functions
    hascounts=which(n>0)
    pRRt  = dbinom(n[hascounts]-k[hascounts],n[hascounts],ee) 
    pHett = choose(n[hascounts],k[hascounts])*(1/(2^n[hascounts]))
    pAAt  = dbinom(k[hascounts],n[hascounts],ee) 

    pRR=pHet=pAA=rep(0,length(n))
    pRR[hascounts]=pRRt
    pHet[hascounts]=pHett
    pAA[hascounts]=pAAt
    eMats=rbind(pRR,pHet,pAA)
    return(eMats)
}

doViterbi=function(r,a, error.rate=.002, tprob=1e-3, sprob=c(.25,.5,.25)) {
   eMats=get_eMats(r,a, error.rate=error.rate)
   eMats[is.na(eMats)]=1e-256
   eMats[eMats==0]=1e-256
   tMats=lapply(1:length(r), function(x) mTmat(tprob))
   return(viterbi(eMats,tMats,startProbs=sprob))
}


##for testing
##newdata[data$ase_list$sample=='one_A01.txt'  | data$ase_list$sample=='one_A02.txt'  ,]
##single threaded version here 
#newdata=data$ase_list %>%
#    group_by(sample, contig) %>% 
#    mutate(hmm=doViterbi(refCount, altCount)) %>%
#    ungroup()

#multi-threaded
cluster=new_cluster(nthreads)
cluster_copy(cluster, c('get_eMats', 'doViterbi', 'mTmat', 'viterbi'))
newdata=data$ase_list %>%
    group_by(sample, contig) %>% 
    partition(cluster)  %>%
    mutate(hmm=doViterbi(refCount, altCount, error.rate=.01, tprob=1e-3)) %>% collect() %>%  ungroup()
saveRDS(newdata,file=output_file)


# visualization
#ns=split(newdata, newdata$sample)

#for(i in names(ns)){
#x=ns[[i]]
#x=x[order(x$contig),]
##ewdata[newdata$sample=='one_A01.txt',]
#plot(-x$refCount, type='h', ylim=c(-250,250), main=i)
#points(-x$refCount, type='h', ylim=c(-250,250))
#points(x$altCount, type='h', ylim=c(-250,250))
#points(50*(x$hmm-2), type='l',col='red', lwd=3)
#abline(v= cumsum(rle(as.character(x$contig))$lengths), lty=2, col='lightblue')
#
#readline()
#}
