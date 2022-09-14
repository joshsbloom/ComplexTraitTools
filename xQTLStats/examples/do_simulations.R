#setwd('~/Dropbox/code/xQTLStats/')
#devtools::install()

library(xQTLStats)


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
sample.size=10e6

#specify cross (see crosses.to.parents)
cross='A'

#get the two parent names 
p1.name=crosses.to.parents[[cross]][1]
p2.name=crosses.to.parents[[cross]][2]

#extract parental genotypes at variant sites between those parents from a reference vcf
vcf.cross=getCrossVCF(ref.vcf,p1.name, p2.name)

doSimulateGenos=F
if(doSimulateGenos){
    #simulate (smaple.size) haploid segregants after n generations, be careful about memory requirements
    geno.matrix=simHaploidSegsFN(vcf.cross, gmaps[[cross]], sample.size, ngenerations=2)
    saveRDS(geno.matrix, file = paste0('/data/xQTL/', cross, '_250K_2.RDS')) 
}
doSimulatePhenos=T
if(doSimulatePhenos){
    geno.matrix=readRDS( file = paste0('/data/xQTL/', cross, '_250K_2.RDS')) 
    #y=simPhenotypesForRefPop(geno.matrix, h2=.2, nQTL=40)
}
library(hibayes)

g=matrix(as.numeric(geno.matrix@.Data), nrow(geno.matrix), ncol(geno.matrix))
tg=t(g)
dtg=duplicated(g)
tg=tg[,-which(dtg)]

y2=simBayesC(t(tg), target.size=.001, h2=0.4)

svF=sparsevb::svb.fit(tg, y2)

fCR=hibayes::bayes(y=y2[,1], X=t(g), model='BayesCpi',niter=20000, nburn=10000, outfreq=10, verbose=TRUE)

fRR=hibayes::bayes(y=y2[,1], X=t(g), model='BayesR',niter=20000, nburn=10000, outfreq=10, verbose=TRUE)
SNPeffect <- fRR$g
gebv <- t(g) %*% SNPeffect    # calculate the estimated genomic breeding value
pve <- apply(t(g),2,var) * (fRR$g^2) / var(y2[,1]) #pheno[,1])    # the phenotypic variance explained for each SNPs
pip <- 1-fRR$nzrate   # the rate of stepping into non-zero effects in MCMC iteration for each SNPs
abline(v=attr(y2,'causal.indices'))

y3=simBayesC(geno.matrix, target.size=.006, h2=0.4)
#g=matrix(as.numeric(geno.matrix@.Data), nrow(geno.matrix), ncol(geno.matrix))
fRR2=bayes(y=y3[,1], X=t(g), model='BayesR',niter=20000, nburn=10000, outfreq=10, verbose=TRUE)
SNPeffect2 <- fRR2$g
gebv2 <- t(g) %*% SNPeffect2    # calculate the estimated genomic breeding value
pve2 <- apply(t(g),2,var) * (fRR2$g^2) / var(y2[,1]) #pheno[,1])    # the phenotypic variance explained for each SNPs
pip2 <- 1-fRR2$nzrate   # the rate of stepping into non-zero effects in MCMC iteration for each SNPs
abline(v=attr(y2,'causal.indices'))






y1=simBayesC(geno.matrix)
y2=simBayesC(geno.matrix, target.size=.0002, h2=0.6)
y3=scale(as.vector(y1)+as.vector(y2))

sel.low=.4
sel.high=.4
tf=tempfile()
expected=simXQTLExperiment(as.vector(y1), geno.matrix, vcf.cross, sel.low=sel.low, sel.high=sel.high, depth.low=50, depth.high=50, vcf.out=tf)
low.tail.expected=expected$low.tail.expected
high.tail.expected=expected$high.tail.expected
#------------------------------------------------------------------------
experiment.vcf=vcfR::read.vcfR(tf) #'/data/xQTL/sim.vcf.gz')

#experiment names (should match col in vcf)
low.tail.name='low.tail.sim'
high.tail.name='high.tail.sim'

    low.tail =getBiallelicCounts(experiment.vcf, low.tail.name)
    high.tail=getBiallelicCounts(experiment.vcf, high.tail.name)

    low.tail$expected=low.tail.expected[low.tail$ID]
    high.tail$expected=high.tail.expected[high.tail$ID]


    #phase everything 
    low.tail  = phaseCounts(vcf.cross, p1.name, low.tail)
    high.tail = phaseCounts(vcf.cross, p1.name, high.tail)

low.tail1  = calcAFD(low.tail, experiment.name='low1',sample.size=sample.size, sel.strength=sel.low)
high.tail1 = calcAFD(high.tail, experiment.name='high1',sample.size=sample.size, sel.strength=sel.high)


results=calcContrastStats(results=list(high.tail1, low.tail1), L='_high1', R='_low1')


h1=plotIndividualExperiment(results, 'high')#,simulatedQTL)
l1=plotIndividualExperiment(results, 'low') #,simulatedQTL)
c1=plotContrast(results, 'high', 'low') #,simulatedQTL)

ggpubr::ggarrange(h1, l1, c1, nrow=3) 

x11()
simulatedQTL=data.frame(chrom=results3$chrom[(attr(y2,"causal.indices"))], 
                    physical.position=results3$physical.position[(attr(y2,"causal.indices")) ])


plotSummary(results, simulatedQTL)




tf2=tempfile()
sel.low=.05
sel.high=.05
expected=simXQTLExperiment(as.vector(y3), geno.matrix, vcf.cross, sel.low=sel.low, sel.high=sel.high, depth.low=50, depth.high=50, vcf.out=tf2)
low.tail.expected=expected$low.tail.expected
high.tail.expected=expected$high.tail.expected
#------------------------------------------------------------------------
experiment.vcf=vcfR::read.vcfR(tf2) #'/data/xQTL/sim.vcf.gz')

#experiment names (should match col in vcf)
low.tail.name='low.tail.sim'
high.tail.name='high.tail.sim'

    low.tail =getBiallelicCounts(experiment.vcf, low.tail.name)
    high.tail=getBiallelicCounts(experiment.vcf, high.tail.name)

    low.tail$expected=low.tail.expected[low.tail$ID]
    high.tail$expected=high.tail.expected[high.tail$ID]


    #phase everything 
    low.tail  = phaseCounts(vcf.cross, p1.name, low.tail)
    high.tail = phaseCounts(vcf.cross, p1.name, high.tail)

low.tail2  = calcAFD(low.tail, experiment.name='low2',sample.size=sample.size, sel.strength=sel.low)
high.tail2 = calcAFD(high.tail, experiment.name='high2',sample.size=sample.size, sel.strength=sel.high)


results2=calcContrastStats(results=list(high.tail2, low.tail2), L='_high2', R='_low2')

results3=calcContrastStats(results=list(high.tail2, high.tail1), L='_high2', R='_high1')


h1=plotIndividualExperiment(results3, 'high2')#,simulatedQTL)
l1=plotIndividualExperiment(results3, 'high1') #,simulatedQTL)
c1=plotContrast(results3, 'high2', 'high1') #,simulatedQTL)

ggpubr::ggarrange(h1, l1, c1, nrow=3) 

x11()
plotSummary(results3,simulatedQTL ) #,











































experiment.vcf.file=system.file('data','B_250K_h2_40_40_QTL_250Kpop_10percenttails.vcf.gz', package='xQTLStats')
experiment.vcf=vcfR::read.vcfR(experiment.vcf.file) #'/data/xQTL/sim.vcf.gz')

low.tail.name='low.tail.sim'
high.tail.name='high.tail.sim'

low.tail=getBiallelicCounts(experiment.vcf, low.tail.name)
high.tail=getBiallelicCounts(experiment.vcf, high.tail.name)

#low.tail$expected=low.tail.expected[low.tail$ID]
#high.tail$expected=high.tail.expected[high.tail$ID]


#get simulated high tail
sel.low  = 0.1
sel.high = 0.1
sample.size=250000
#phase everything 
low.tail  = phaseCounts(vcf.cross, p1.name, low.tail)
high.tail = phaseCounts(vcf.cross, p1.name, high.tail)


#----------------------------------------------------------------------------
low.tail  = calcAFD(low.tail, experiment.name='low',sample.size=sample.size, sel.strength=sel.low)
high.tail = calcAFD(high.tail, experiment.name='high',sample.size=sample.size, sel.strength=sel.high)


results=calcContrastStats(results=list(high.tail, low.tail), L='_high', R='_low')


h1=plotIndividualExperiment(results, 'high')#,simulatedQTL)
l1=plotIndividualExperiment(results, 'low') #,simulatedQTL)
c1=plotContrast(results, 'high', 'low') #,simulatedQTL)

ggpubr::ggarrange(h1, l1, c1, nrow=3) 

x11()
plotSummary(results) #, simulatedQTL)







#for(cross in names(crosses.to.parents)){

#geno.matrix=simHaploidSegsFN(vcf.cross, gmaps[[cross]], sample.size, ngenerations=24)
#saveRDS(geno.matrix, file = paste0('/data/xQTL/', cross, '_250K_24.RDS')) 

    #saveRDS(vcf.cross, file = paste0('/data/xQTL/', cross, '_vcf.RDS')) 
#}















#  saveRDS(vcf.cross, file = paste0('/data/xQTL/', cross, '_vcf.RDS')) 

    #simulation bit ----------------------------------------------------------
geno.matrix=simHaploidSegsFN(vcf.cross, gmaps[[cross]], sample.size, ngenerations=2)

    #saveRDS(geno.matrix, file = paste0('/data/xQTL/', cross, '_250K_24.RDS')) 
#}
    #'/data/xQTL/B_48.RDS')




#library(ggplot2)
#high.plot= 
#    ggplot(results.data, aes(x=physical.position,y=p1.high/(p1.high+p2.high)))+
#    geom_point(size=0.3,alpha=0.6, color='gray21')+
#    geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')+
#    facet_grid(~chrom, scales='free_x', space='free_x')+
#    geom_ribbon(aes(ymin=afd.high-1.96*afd.high.se, ymax=afd.high+1.96*afd.high.se, fill='grey'), 
#                alpha=0.7,linetype='dashed', color='grey')+
#    geom_line(aes(x=physical.position, y=afd.high),color='red', size=2, alpha=1)+ 
#    geom_line(aes(x=physical.position, y=expected.af.high),color='black', size=.5)+
#    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none')+ggtitle('high tail')
#
#low.plot=
#    ggplot(results.data, aes(x=physical.position,y=p1.low/(p1.low+p2.low)))+
#    geom_point(size=0.3,alpha=0.6, color='gray21')+
#    geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')+
#    facet_grid(~chrom, scales='free_x', space='free_x')+
#    geom_ribbon(aes(ymin=afd.low-1.96*afd.low.se, ymax=afd.low+1.96*afd.low.se, fill='grey'), 
#                alpha=0.7,linetype='dashed', color='grey')+
#    geom_line(aes(x=physical.position, y=afd.low),color='red', size=2, alpha=1)+ 
#    geom_line(aes(x=physical.position, y=expected.af.low),color='black', size=.5)+
#    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none')+ggtitle('low tail')
#
##ggpubr::ggarrange(high.plot, low.plot, nrow=2) 
#
#contrast.plot=ggplot(results.data, aes(x=physical.position,y=(p1.high/(p1.high+p2.high))-(p1.low/(p1.low+p2.low))))+
#    geom_point(size=0.3,alpha=0.6, color='gray21')+
#    geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')+
#    facet_grid(~chrom, scales='free_x', space='free_x')+
#    geom_ribbon(aes(ymin=afd.contrast-1.96*afd.contrast.se, ymax=afd.contrast+1.96*afd.contrast.se, fill='grey'), 
#                     alpha=0.7,linetype='dashed', color='grey')+
#    geom_line(aes(x=physical.position, y=afd.contrast),color='red', size=2, alpha=1)+ 
#    geom_line(aes(x=physical.position, y=expected.af.high-expected.af.low),color='black', size=.5)+
#    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none')+
#    ggtitle('allele frequence contrast')
#
##ggplot(results.data, aes(x=physical.position, y=afd.contrast.LOD))+geom_line(color='black', size=1.5)+ 
##        facet_grid(~chrom, scales='free_x', space='free_x')+
##        theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle('LOD')
#neglogp.plot=ggplot(results.data, aes(x=physical.position, y=-log10(afd.contrast.p)))+geom_line(color='black', size=1.5)+ 
#         geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')+
#        facet_grid(~chrom, scales='free_x', space='free_x')+
#        geom_hline(aes(yintercept=-log10(.05/600)), color='red')+
#        theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle('-log10(p)')
#
# ggpubr::ggarrange(high.plot, low.plot, contrast.plot, neglogp.plot, nrow=4) 














