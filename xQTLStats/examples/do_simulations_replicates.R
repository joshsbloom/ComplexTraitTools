#setwd('~/Dropbox/code/xQTLStats/')
#devtools::install()

library(xQTLStats)

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
sample.size=2.5e5

#specify cross (see crosses.to.parents)
cross='B'

#get the two parent names 
p1.name=crosses.to.parents[[cross]][1]
p2.name=crosses.to.parents[[cross]][2]

#extract parental genotypes at variant sites between those parents from a reference vcf
vcf.cross=getCrossVCF(ref.vcf,p1.name, p2.name)

doSimulateGenos=F
if(doSimulateGenos){
    #simulate (smaple.size) haploid segregants after n generations, be careful about memory requirements
    geno.matrix=simHaploidSegsFN(vcf.cross, gmaps[[cross]], sample.size, ngenerations=12)
    #saveRDS(geno.matrix, file = paste0('/data/xQTL/', cross, '_250K_12.RDS')) 
}



geno.matrix=readRDS( file = paste0('/data/xQTL/', cross, '_250K_24.RDS')) 
simout.dir='/data/xQTL/simulation/sims24/'
dir.create(simout.dir)
#precompute this 
splitit=cut(1:ncol(geno.matrix),5, labels=F)
saf=lapply(1:5, function(h){
                   idx=which(splitit==h)
                   return(Rfast::rowsums(geno.matrix[,idx])) #/ncol(sel.indv.x)
                               })
sel.indv.af=Rfast::rowsums(do.call('cbind', saf))/ncol(geno.matrix)
attr(geno.matrix, 'sel.indv.af')=sel.indv.af


doSimulatePhenos=T
#sel.high=0.2
#sel.low=0.2
#depth.high=20
#depth.low=20
#if(doSimulatePhenos){

#sims2
h2.set=0.2
nQTL.set=c(40,20,10,5)
sel.high.set=c(.5,.25,.1,.01,.001,.0001)
sel.low.set=1
depth.set=c(5,10,20,50,100,200,500)
nrep.set=c(1,2,3,4,5)

sim.space=expand.grid(h2=h2.set, nQTL=nQTL.set, sel.high=sel.high.set, sel.low=sel.low.set, depth=depth.set, nrep=nrep.set)
sim.space$ID=apply(sim.space, 1, paste, collapse=':')
saveRDS(sim.space, file=paste0(simout.dir, 'simSpace.RDS'  ))

for(n in 1:nrow(sim.space)){
    h2=sim.space$h2[n]
    nQTL=sim.space$nQTL[n]
    sel.high=sim.space$sel.high[n]
    sel.low=sim.space$sel.low[n]
    depth.low=sim.space$depth[n]
    depth.high=sim.space$depth[n]
    nrep=sim.space$nrep[n]
   # depth.low=20
   # depth.high=20

    y=simPhenotypesForRefPop(geno.matrix, h2=h2, nQTL=nQTL)

    xqtl_reps=list()
    for(i in 1:nrep){
        tf=tempfile()
        expected=simXQTLExperiment(y, geno.matrix, vcf.cross, sel.low=sel.low, sel.high=sel.high,depth.low=depth.low, depth.high=depth.high, vcf.out=tf)
        low.tail.expected=expected$low.tail.expected; high.tail.expected=expected$high.tail.expected

        #------------------------------------------------------------------------
        experiment.vcf=vcfR::read.vcfR(tf) #'/data/xQTL/sim.vcf.gz')

        #experiment names (should match col in vcf)
        low.tail.name='low.tail.sim'
        high.tail.name='high.tail.sim'

        low.tail=getBiallelicCounts(experiment.vcf, low.tail.name)
        high.tail=getBiallelicCounts(experiment.vcf, high.tail.name)

        low.tail$expected=low.tail.expected[low.tail$ID]
        high.tail$expected=high.tail.expected[high.tail$ID]


        #phase everything 
        low.tail  = phaseCounts(vcf.cross, p1.name, low.tail)
        high.tail = phaseCounts(vcf.cross, p1.name, high.tail)

        #----------------------------------------------------------------------------
        low.tail  = calcAFD(low.tail, experiment.name=paste0('low',i),sample.size=sample.size, sel.strength=sel.low)
        high.tail = calcAFD(high.tail, experiment.name=paste0('high',i),sample.size=sample.size, sel.strength=sel.high)

        xqtl_reps[[as.character(i)]]$low.tail=low.tail
        xqtl_reps[[as.character(i)]]$high.tail=high.tail
    }

    results=lapply(xqtl_reps, function(x) x$low.tail)
    results.low=calcMetaAFD(results, meta_name='meta_low')
    results=lapply(xqtl_reps, function(x) x$high.tail)
    results.high=calcMetaAFD(results, meta_name='meta_high')
    results=calcContrastStats(list(results.high, results.low), L='meta_high', R='meta_low')

    simulatedQTL=data.frame(chrom=results$chrom[abs(attr(y,"add.qtl"))], 
                    physical.position=results$physical.position[abs(attr(y,"add.qtl"))])
    attr(results, 'simulatedQTL')=simulatedQTL
    saveRDS(results, file=paste0(simout.dir, n))

}


#embed this in plotting function
    h1=plotIndividualExperiment(results, 'high1',attr(results, 'simulatedQTL'))
    l1=plotIndividualExperiment(results, 'low1',attr(results, 'simulatedQTL'))
    c1=plotContrast(results, 'high1', 'low1',attr(results, 'simulatedQTL'))

    ggpubr::ggarrange(h1, l1, c1, nrow=3) 

    x11()
    plotSummary(results, attr(results,'simulatedQTL'))


simout.dir='/data/xQTL/simulation/sims24/'
sim.space=readRDS(file=paste0(simout.dir, 'simSpace.RDS'  ))
results=list()
for(i in 1:nrow(sim.space)){
    results[[as.character(i)]]=readRDS(paste0(simout.dir, i))

}
sim.space$detected=0
sim.space$power=0
for(i in 1:nrow(sim.space)){
    h=results[[as.character(i)]]
    pp=attr(h,'simulatedQTL')$physical.position
    sim.space$power[i]=sum(h$p[h$physical.position%in%pp]<.05/600, na.rm=T)/sum(!is.na(h$p[h$physical.position%in%pp])) #[h$physical.position%in%pp]<.05/600))
}

#sim.space$power=sim.space$detected/sim.space$nQTL
ggplot(sim.space, aes(x=nrep, y=power, color=as.factor(sel.high)))+facet_grid(depth~nQTL)+geom_point()+geom_line()+ylim(0,1)

h=results[[as.character(i)]]


    ss=sim.space[sim.space$nQTL==40 & sim.space$nrep==1 & sim.space$depth==100,]
    i=rownames(ss)[3]
    print(sim.space[i,])
    h=results[[as.character(i)]]

    h1=plotIndividualExperiment(results[[as.character(i)]], 'high1',attr(results[[as.character(i)]], 'simulatedQTL'))
    l1=plotIndividualExperiment(results[[as.character(i)]], 'low1',attr(results[[as.character(i)]], 'simulatedQTL'))
    c1=plotContrast(results[[as.character(i)]], 'high1', 'low1',attr(results[[as.character(i)]], 'simulatedQTL'))
    ggpubr::ggarrange(h1, l1, c1, nrow=3) 
    plotSummary(results[[as.character(i)]], attr(results[[as.character(i)]],'simulatedQTL'))
z
