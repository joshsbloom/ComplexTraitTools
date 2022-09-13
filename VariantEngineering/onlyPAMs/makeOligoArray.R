library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(RFLPtools)
library(seqinr)
library(intervals)
library(VariantAnnotation)
library(GenomicFeatures)
## note, the bioconductor TxDb object is old and broken, we'll go ahead and build our own
##library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
##txdb=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene


# Global Variables
working.dir='/media/jbloom/d1/CRISPR_base_editor/'

# import some functions 
source(paste0(working.dir, 'code/designCRISPR_gRNAs_fx.R'))

sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
# make our own txdb structure (motivated by annotation discrepancy found for GRX3 ORF and issues with dubious ORFs)
# pulled from http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff (beware of changes)
txdb=makeTxDbFromGFF(file=paste0(working.dir, 'reference/saccharomyces_cerevisiae.gff'))

guides=readRDS(paste0(working.dir, 'results/guideTable.RDS'))
guidesNoCs=readRDS(paste0(working.dir, 'results/guideTable_noEdits.RDS'))
coding.effects=readRDS(paste0(working.dir, 'results/codingEffects.RDS'))

# write resulting tables as tsv files
## write.table(coding.effects, file=paste0(working.dir, 'results/codingEffectsTable.tsv'), sep = "\t", row.names = F, quote=F)
## write.table(guides, file=paste0(working.dir, 'results/guideTable.tsv'), sep = "\t", row.names = F, quote=F)

table(coding.effects$CONSEQUENCE, coding.effects$essential, coding.effects$dubious)
with(coding.effects[coding.effects$CONSEQUENCE=='nonsynonymous',], boxplot(maxAbsProvean~essential))
with(coding.effects[coding.effects$CONSEQUENCE=='nonsynonymous',], t.test(maxAbsProvean~essential))

## additional info (http://www.yeastgenome.org/phenotype/increased_competitive_fitness/overview)
# makeGRangesFromDataFrame(df)

# some guide localization stats
##gPAM=guides$X
##gPAM$guideIndex=guides$guideIndex
##codingGuides = locateVariants(gPAM, txdb, CodingVariants(), ignore.strand=TRUE)
##promoterGuides   = locateVariants(gPAM, txdb, PromoterVariants(455, 0), ignore.strand=TRUE)
##intergenicGuides = locateVariants(gPAM, txdb, IntergenicVariants(1000000, 1000000), ignore.strand=TRUE)


# Defining oligo subsets

# Essential stops 
    eStops = coding.effects[coding.effects$CONSEQUENCE=='stop_gained' & 
                            coding.effects$essential & !coding.effects$dubious & !coding.effects$U6terminator & !coding.effects$overlappingORFs &
                            coding.effects$guidePAM12PMcount<2, ] 
    # diagnostics
    par(mfrow=c(2,1))
    hist(table(eStops$GENEID), xlim=c(0,60), breaks=1000, xlab='gRNAs per gene')

# Non-essential stops
    neStops = coding.effects[coding.effects$CONSEQUENCE=='stop_gained' &
                             !coding.effects$essential & !coding.effects$dubious & !coding.effects$U6terminator & !coding.effects$overlappingORFs & 
                             coding.effects$guidePAM12PMcount<2, ] 
    hist(table(neStops$GENEID), xlim=c(0,60), breaks=1000, xlab='gRNAs per gene')

    # one idea, get one gRNA for each unique gene
    # length(unique(neStops$GENEID)) #4321
    #set.seed(1234)
    #neStopsRandomized=neStops[sample(1:nrow(neStops)),]
    #neStopsSubset=neStopsRandomized[match(unique(neStops$GENEID), neStopsRandomized$GENEID),]

    # another idea (just get n random oligos)
    set.seed(50)
    n.neStopOligos=5500
    neStopsSubset=neStops[sample(1:nrow(neStops), n.neStopOligos),]
    hist(table(neStopsSubset$GENEID))

# Provean Scores
    # Choice of 6 here is arbitrary. Check literature 
    eProv = coding.effects[coding.effects$CONSEQUENCE=='nonsynonymous' & 
                           coding.effects$essential & !coding.effects$dubious & !coding.effects$overlappingORFs & 
                           !coding.effects$U6terminator & 
                           coding.effects$guidePAM12PMcount<2 &
                           coding.effects$maxAbsProvean>5 & !is.na(coding.effects$maxAbsProvean),  ] 
    set.seed(4321)
    n.eProvOligos=5500
    eProvSubset=eProv[sample(1:nrow(eProv),n.eProvOligos),]

# fourth set from guideNoCs as controls either merged in individually or independently


# Primers used for divided oligo experiment ... F2/R2 for essential stops
skpF1='GGGTCACGCGTAGGA' #3’  sharon2012_F
skpF2='CGCGTCGAGTAGGGT' #3’  firstAUG_F
skpF3='CGATCGCCCTTGGTG' #3’  firstUTR_F
skpF4='GGTCGAGCCGGAACT' #3’  upstream_F
skpF5='TCCCGGCGTTGTCCT'
skpF6='CGCAGGGTCCAGAGT'

skpR1r='GTTCCGCAGCCACAC' #3’ sharon2012_R
skpR2r='GCCGTGTGAAGCTGG' #3’ firstAUG_R
skpR3r='GGTTTAGCCGGCGTG' #3’ firstUTR_R
skpR4r='GGATGCGCACCCAGA' #3’ upstream_R
skpR5r='GCTCCGTCACTGCCC'
skpR6r='GTTCGCGCGAAGGAA'

skpR1=reverseComplement(DNAString(skpR1r))
skpR2=reverseComplement(DNAString(skpR2r))
skpR3=reverseComplement(DNAString(skpR3r))
skpR4=reverseComplement(DNAString(skpR4r))
skpR5=reverseComplement(DNAString(skpR5r))
skpR6=reverseComplement(DNAString(skpR6r))
#-----------------------------------------------------------------------------

# recognition sites for restriction enzymes  
#1 or more instance of 
eagI='CGGCCG'
# search for 
#more than 1 instance of 
MluI='ACGCGT'
BstEII='GGTNACC'
SphI='GCATGC' 
# ----------------------------------

#  WARNING prototype, sequences here were specific to the essential PTC array, modify!!
#  Essential PTC experiment used skpF2/skpR2r
# for example 

chosenOligos=do.call('rbind', list(eStops, neStopsSubset, eProvSubset))

# sanity check not picking duplicate guides 
x=chosenOligos[which(duplicated(chosenOligos$guide, fromLast=FALSE)),]
y=chosenOligos[which(duplicated(chosenOligos$guide, fromLast=TRUE)),]
x[order(x$guide),]
y[order(y$guide),]

final.oligo=paste0(skpF2, 'GGTGACC', chosenOligos$guide,'GTTTTAGAGCATGC', 'CGATCGAT','ACGCGT',  skpR2)
final.oligoD= DNAStringSet(final.oligo)

# Reverse complement if more A's then T's ... Olga please sanity check this code 
Acnt=letterFrequency(final.oligoD, 'A')
Tcnt=letterFrequency(final.oligoD, 'T')

revCompOligos=ifelse(Acnt>Tcnt,TRUE,FALSE)
selectedRevComp=reverseComplement(final.oligoD[revCompOligos]) 
final.oligoDrc=final.oligoD
final.oligoDrc[revCompOligos]=selectedRevComp

#pb=txtProgressBar(min=1, max=length(final.oligoD), style=3)
#for(i in 1:length(final.oligoD) ){
#    setTxtProgressBar(pb,i)
#    if(Acnt[i]>Tcnt[i]) final.oligoDrc[i]=reverseComplement(final.oligoDrc[i]) 
#}
#close(pb)

# Change this as necessary given restriction enzymes used
eagI.cnt=vcountPattern(eagI, final.oligoDrc)
MluI.cnt=vcountPattern(MluI, final.oligoDrc)
SphI.cnt=vcountPattern(SphI, final.oligoDrc)
BstEII.cnt=vcountPattern(BstEII, final.oligoDrc, fixed=F)
bad.oligos=(eagI.cnt>0 | MluI.cnt>1 | SphI.cnt >1 | BstEII.cnt >1)

fos=data.frame(guideIndex=chosenOligos$guideIndex, final.oligo.toSynth=as.character(final.oligoDrc), stringsAsFactors=F)
fos=fos[-which(bad.oligos),]
write.table(fos, file=paste0(working.dir, 'results/final_oligos_to_synthesize.txt'), sep='\t', row.names=F, quote=F)
oligo.column=fos[,2]
oligo.column[!duplicated(oligo.column)]
write.table(fos[,2], file=paste0(working.dir, 'results/final_oligos_to_synthesize_for_oligo_printer.txt'),row.names=F, col.names=F, quote=F)
