library(xQTLStats)
library(vcfR)
data(crosses.to.parents)

gmaps=readRDS(system.file('reference', 'yeast_gmaps.RDS', package='xQTLStats'))
#reference vcf
ref.vcf=system.file('reference', 'parents_w_svar_sorted.vcf.gz', package='xQTLStats')
cross='3049'

p1.name=crosses.to.parents[[cross]][1]
p2.name=crosses.to.parents[[cross]][2]
   
vcf.cross=getCrossVCF(ref.vcf,p1.name, p2.name)

experiment.vcf=read.vcfR('/data/xQTL/Laila/xqtl_49_filtered.vcf.gz')

cross.chr=getFIX(vcf.cross)[,'CHROM']
cross.pos=getFIX(vcf.cross)[,'POS']
cross.ref=getFIX(vcf.cross)[,'REF']
cross.alt=getFIX(vcf.cross)[,'ALT']
cross.id=paste(cross.chr, '_' , cross.pos, '_', cross.ref, '_', cross.alt, sep='')

experiment.chr=getFIX(experiment.vcf)[,'CHROM']
experiment.pos=getFIX(experiment.vcf)[,'POS']
experiment.ref=getFIX(experiment.vcf)[,'REF']
experiment.alt=getFIX(experiment.vcf)[,'ALT']
experiment.id=paste(experiment.chr, '_' , experiment.pos, '_', experiment.ref, '_', experiment.alt, sep='')

#assuming experiment.id has fewer entries than cross.id 
experiment.id.to.keep=experiment.id[which(experiment.id %in% cross.id)]
cross.id.to.keep=experiment.id[which(cross.id %in% experiment.id)]

#first reduce 
vcf.cross=vcf.cross[which(cross.id %in% experiment.id.to.keep),]
#experiment.vcf.big=experiment.vcf[which(experiment.id %in% experiment.id.to.keep),]

#print experiment
colnames(experiment.vcf@gt)[-1]

experiment.vcf=addID(experiment.vcf)

cntrl1=getBiallelicCounts(experiment.vcf, '49_parental_1.bam')
cntrl2=getBiallelicCounts(experiment.vcf, '49_parental_2.bam')
eth15 =getBiallelicCounts(experiment.vcf, '49_15_etoh.bam')

cntrl1=phaseCounts(vcf.cross, p1.name, cntrl1)
cntrl2=phaseCounts(vcf.cross, p1.name, cntrl2)
eth15 =phaseCounts(vcf.cross, p1.name, eth15)

cntrl1=calcAFD(cntrl1, experiment.name='cntrl1', sample.size=1e6, sel.strength=1)
cntrl2=calcAFD(cntrl2, experiment.name='cntrl2', sample.size=1e6, sel.strength=1)
eth15 =calcAFD(eth15, experiment.name='eth15',   sample.size=1e6, sel.strength=.0004)

controls=list(cntrl1, cntrl2)
results.cntrl = calcMetaAFD(controls, meta_name='meta_cntrl')
results.exp   = calcMetaAFD(list(eth15), meta_name='meta_eth15')
results=calcContrastStats(list(results.exp, results.cntrl), L='meta_eth15', R='meta_cntrl')

p1=plotIndividualExperiment(results, 'cntrl1')#,simulatedQTL)
p2=plotIndividualExperiment(results, 'cntrl2')#,simulatedQTL)
e1=plotIndividualExperiment(results, 'eth15')
c1=plotContrast(results, 'eth15', 'cntrl1') #,simulatedQTL)

ggpubr::ggarrange(p1, p2, e1, c1, nrow=4) 

x11()
plotSummary(results)


