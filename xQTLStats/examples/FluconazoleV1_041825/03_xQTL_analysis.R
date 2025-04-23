library(xQTLStats)
library(qs)
library(stringr)
library(ggpubr)

#some global variables ----------------------------------------
# locations of things -----------------------------------------------------------------------
#location of vcf file from 1011 isolates, preprocessed as a qs object
vcf.qs.file='/u/project/kruglyak/jsbloom/PUBLIC_SHARED/yeast/1011/mega_filtered.qs'
vcf.qs.file='/media/hoffman2/PUBLIC_SHARED/yeast/1011/mega_filtered.qs'

#location of GATK ASEReadCount tables, output from 02_align.sh
count_dir='/data2/xQTL/FluconazoleV1_041825/count_files/'

#location of key file, grab key from github
key.file='/home/jbloom/Dropbox/code/ComplexTraitTools/xQTLStats/examples/FluconazoleV1_041825/flatkey.csv'
#---------------------------------------------------------------------------------------------

# file that contains existing genetic maps
#need to pass in a genetic map, extract one from xQTLStats
gmaps=readRDS(system.file('reference', 'yeast_gmaps.RDS', package='xQTLStats'))
gmap=gmaps[['A']]
#unique chromosomes 
uchr=paste0('chr', as.character(as.roman(1:16)))
#----------------------------------------------------------------

mega_filtered=qs::qread(vcf.qs.file) #'/data1/yeast/reference/pop_vcfs/1011/mega_filtered.qs')
vcf=mega_filtered$vcf
# pull out the genotype entries (num vector w/ marker + genotype info added as charac vectors)
gt=mega_filtered$gt
# remove the og file from memory, retain vcf object and gt object only
rm(mega_filtered)

samples=readr::read_csv(key.file)

countdfs=list()
for(n in 1:nrow(samples)) {
   print(n) 
    s=samples$'Sample ID'[n]
    sname = paste0(s, '.txt') ##- "Sample012.txt"
    p.names <- c(samples$parent_a[n],samples$parent_alpha[n]) #samples$parent_alpha[n]) #samc("S288C", "ACD")
    countdf <- makeCountTablesGATK(count_dir, sname, p.names, vcf, gt, gmap)
    #snn = str_c(paste(p.names, collapse="_"), "a")
    countdfs[[s]] = countdf
}

#calculate coverage per sample at the expected segregating sites 
tdepth.v=sapply(countdfs, function(x) sum(x$ref+x$alt)/nrow(x))

#calculate allele frequence differences smoothed per chromosome given expected sample size and selection strength 
afds=lapply(names(countdfs[tdepth.v>1]), function(snn) {
       calcAFD(countdfs[[snn]], experiment.name=snn,sample.size=1e4, sel.strength=.9) #, bin.width=3000, eff.length=2000, uchr=uchr)
      })
names(afds)=names(countdfs[tdepth.v>1])

# pre-generate allele frequency plots 
plots=lapply(names(afds), function(snn) {
    plotIndividualExperiment(afds[[snn]], snn) 
      })

#list samples
samples$oldLabel_parent_a[tdepth.v>1]
#visualize one of them 
plots[[7]]

# now calculate the contrast statistic between a drug condition and its dmso control
#87 vs 82
results=calcContrastStats(results=list(afds[[82]], afds[[87]]), #S2_N2_XZ1516_7_RNAi, LWM1_S56_N2_XZ1516_3_RNAi),
                          L=paste0('_', names(afds)[82]), R=paste0('_', names(afds)[87]) ) #'_high1', R='_unsel1')
#visualize contrast
plotContrast(results,names(afds)[82],names(afds)[87])
#visualize summary statistics from contrast
plotSummary(results) #ontrast(results,names(afds)[82],names(afds)[87])

# from xQTLStats/R/xQTL_stats.R
# use ggpubr functions to stack more plots together 
#a=ggarrange(plots[[1]],plots[[2]],  nrow=2)

