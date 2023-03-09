library(ggplot2)

# logic here is we keep the mean difference between the two homozygotes fixed at 1 and then vary the error (noise)
Beta=1
# expectation for total number of effective tests (number of genotype pcs) 
effective.ntests=100
# total number of simulations per paramter
nsims=1000
#set alpha for significance taking into account total effective number of tests
sig.level=.05/effective.ntests

# parameters to vary
n.indivs=c(25,50, 75,100,150,200,300,500,1000)
H2s=c(.01,.025,.05,.075,.1,.2,.3,.4,.5)

# data frame containing parameters that vary and power estimation that we calculate
params=expand.grid(n.indivs, H2s)
colnames(params)=c('n', 'h2')
params=data.frame(params)
params$power=0

#for each combo of parameters, simulate and estimate power 
for(j in 1:nrow(params) ) {
    n.indiv=params[j,1]
    H2=params[j,2]

    genN=rmultinom(1, n.indiv,c(.25,.5,.25))[,1]
    #allelic dosages
    genotype=c(rep(0,genN[1]), rep(1,genN[2]), rep(2,genN[3]))

    simP=replicate(nsims, {
        #allelic dosage
        p=genotype*Beta
        #add noise
        p=p+rnorm(n.indiv,mean=0,sd=sqrt((1-H2)/H2*var(genotype)))
        #stripchart(p~genotype, vertical=T, method='jitter', main=H2)
        return(p)
    })

    pvals=rep(1,nsims)
    # calculate a likelihood ratio statistic comparing a model with effect of genotype to model without
    for(i in 1:nsims) {
        nullM=lm(simP[,i]~1)
        fullM=lm(simP[,i]~genotype)
        pvals[i]=pchisq(-2*(logLik(nullM)-logLik(fullM)),1, lower.tail=F)
    }
    # significance is p<sig.level defined above (takes into account multiple testing)
    power.estimate = sum(pvals<sig.level)/nsims
    params[j,3]=power.estimate
    print(params[j,])
}

# rethink how to visualize 
ggplot(params, aes(x=h2, y=power))+facet_wrap(~n)+geom_point()+geom_line()
