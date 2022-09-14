# Code for the design of gRNA libraries for bulk editing of S.Cerevisiae
# targeting all variants to PAM sites 
# yeast genome sequence build sacCer3
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
    
#where is code
base.dir='~/Dropbox/code/ComplexTraitTools/'
source(paste0(base.dir, 'VariantAnnotation/functions/proveanFx.R'))


# import some functions 
source(paste0(base.dir, 'VariantEngineering/onlyPAMs/designCRISPR_gRNAs_fx.R'))

sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
    # make our own txdb structure (motivated by annotation discrepancy found for GRX3 ORF and issues with dubious ORFs)
    # pulled from http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz and gunzipped it (beware of changes)

ref.data.dir=paste0(base.dir, 'VariantAnnotation/yeast_ref/') 

txdb=GenomicFeatures::makeTxDbFromGFF(paste0(base.dir,'VariantAnnotation/yeast_ref/saccharomyces_cerevisiae.gff.gz'))

unique.chrs=paste0('chr', as.roman(1:16))
sc3.set=buildSacCer3_genome_dictionaryFR(sacCer3, unique.chrs)
provean.score.file=paste0(ref.data.dir,'provean_scores.RDS')
provean_scores=readRDS(provean.score.file)
meltedProveanScores=meltProvean(provean_scores)

# for more details about ORF annotations see http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_dna.README



# for each chromosome build GRanges objects of gRNA sequences and predict the effect of the base editor on coding sequence
# lists to hold the gRNA tables 
guide.tables=list()
#guide.tables.No.Cs=list()
coding.effect.tables=list()
for(chr.name in unique.chrs) {
    print(chr.name)
    guide.tables[[chr.name]]        =get_gRNAs_and_predict_edit(match(chr.name,unique.chrs),sacCer3, txdb, sc3.set)
#    guide.tables.No.Cs[[chr.name]]  =get_gRNAs_and_predict_edit(match(chr.name,unique.chrs),sacCer3, txdb, sc3.set, returnIdentical=TRUE)
    for(ns in names(guide.tables[[chr.name]]$editedSequenceFwdL)) {
        print(ns)
            
        coding.effects=predictCoding(guide.tables[[chr.name]],txdb, sacCer3, 
                                                         DNAStringSet(guide.tables[[chr.name]]$editedSequenceFwdL[,ns]), 
                                                         ignore.strand=TRUE)
        pC=mungeCodingEffects(coding.effects, ns, meltedProveanScores)
        #pC$gEditInfo=guide.tables[[chr.name]]$gEditInfo[,ns]
        #pC$editedSequenceFwdL=DNAStringSet(guide.tables[[chr.name]]$editedSequenceFwdL[,ns])
        coding.effect.tables[[paste0(chr.name,':',ns) ]]=pC
    }
}

for(ns in names(coding.effect.tables) ){   coding.effect.tables[[ns]]$PAMbreak=ns           }



# Build DataFrames concatenating per-chromosome tables (see package S4vectors, a fast flexible extension to data.frame)
guides=do.call('rbind', lapply(guide.tables, DataFrame)) 
#guidesNoCs=do.call('rbind', lapply(guide.tables.No.Cs, DataFrame))
coding.effects=do.call('rbind', lapply(coding.effect.tables, DataFrame)) 

# Add additional features to guides
guides = augment_guides(guides, addScores=F)
#guidesNoCs = augment_guides(guidesNoCs, addScores=F)

# Add additional features to coding.effects  
coding.effects$guidePAM12PMcount=guides$guidePAM12PMcount[match(coding.effects$guideIndex, guides$guideIndex)]
coding.effects$guideGCcontent=guides$guideGCcontent[match(coding.effects$guideIndex, guides$guideIndex)]
coding.effects$U6terminator=guides$U6terminator[match(coding.effects$guideIndex, guides$guideIndex)]
coding.effects$guideScore=guides$guideScore[match(coding.effects$guideIndex, guides$guideIndex)]

# Annotate 'Dubious' ORFs 
# for more details about ORF annotations see http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_dna.README

# For simplicity, here we are going to redefine dubious as (dubious or pseudogene). There are 784 dubious and 12 psuedogenes.
# get_pseudogenes(working.dir, orf_coding) #"YAR061W" "YCL074W" "YCL075W" "YFL056C" "YIL170W" "YIL171W" "YIL174W" "YIL175W" "YLL016W" "YLL017W" "YPL275W" "YPL276W"
##working.dir='/home/jbloom/Dropbox/code/ComplexTraitTools/VariantEngineering/onlyPAMs/'

orf_coding=seqinr::read.fasta(paste0(ref.data.dir,'orf_coding.fasta.gz'))
coding.effects$dubious=!(coding.effects$GENEID %in% names(orf_coding))

# Essential ORFs as annotated here (http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt)
e_orfs=read.delim(paste0(ref.data.dir, 'E_orfs_column.txt'), stringsAsFactors=F, header=F, sep='\t')[,1]
coding.effects$essential=coding.effects$GENEID %in% e_orfs
# NOTE!!! BE CAREFUL, A GENE CAN BE ANNOTATED AS ESSENTIAL AND DUBIOUS, USUALLY DUE TO OVERLAP OF A BONA FIDE ESSENTIAL GENE AND A DUBIOUS ORF

# flag the gRNAs that could edit overlapping ORFs 
coding.effects$overlappingORFs=as.vector(ifelse(sapply(coding.effects$CDSID, length)>1,TRUE,FALSE))

# reannotate nonsynonymous PTCs as PTCs
coding.effects=annotate_PTCs(coding.effects)

g=coding.effects$GENEID
r=coding.effects$REFAA.delta
p=coding.effects$POSAA.delta
v=coding.effects$VARAA.delta
grlv=paste0(g,':',r, ':' ,p, ':',v)

coding.effects$AAeffect=grlv
saveRDS(coding.effects, file='/data/yeast/yeast_oligos/onlyPams/tables/codingEffects.RDS')
saveRDS(guides, file='/data/yeast/yeast_oligos/onlyPams/tables/guides.RDS')








#nvar.coding=readRDS('/data/yeast/yeast_oligos/1002genomes/1002codingEffects.RDS')
#
##pedit=coding.effects$AAeffect %in% nvar.coding$AAeffect
#pedit = nvar.coding$AAeffect %in% coding.effects$AAeffect & nvar.coding$CONSEQUENCE!='synonymous'
#
#length(unique(coding.effects[coding.effects$CONSEQUENCE!='synonymous'& !grepl(':c\\(',coding.effects$AAeffect),]$AAeffect))
## 8371329
#
#ces=coding.effects[coding.effects$CONSEQUENCE!='synonymous'& !grepl(':c\\(',coding.effects$AAeffect),]
#
#
#length(unique(nvar.coding[pedit,]$AAeffect))
#[1] 133606
#length(unique(nvar.coding[nvar.coding$CONSEQUENCE!='synonymous',]$AAeffect))
#[1] 556683
#
#
#par(mfrow=c(3,1))
#hist(unlist(nvar.coding$proveanEffects), xlab='provean score', main='all natural coding variants', xlim=c(-15,15))
#hist(unlist(nvar.coding[pedit,]$proveanEffects), xlab='provean score', main='all NGG editable natural coding variants', xlim=c(-15,15))
#hist(unlist(ces$proveanEffects), xlab='provean score', main='all NGG editable coding variants', xlim=c(-15,15))
#
#
#rv=paste(r,v)
#trv=table(rv)
#
#eq=mapply(function(x,y) {x==y}, r, v)
#
#singleAA=trv[!grepl('\\"', names(trv))]
#
#singleAAdf=data.frame(refAA=sapply(strsplit(names(singleAA),' '), function(x) x[1]),
#           varAA=sapply(strsplit(names(singleAA),' '), function(x) x[2]),
#           Instances=singleAA)
#singleAAdf=singleAAdf[-13,]
#ggplot(singleAAdf, aes(x=refAA, y=varAA, fill=Instances.Freq))+geom_tile()+scale_fill_viridis_c()
#
#r=nvar.coding$REFAA.delta
#v=nvar.coding$VARAA.delta
#rvn=paste(r,v)
#trvn=table(rvn)
#singleAAn=trvn[!grepl('\\"', names(trvn))]
#singleAAdfn=data.frame(refAA=sapply(strsplit(names(singleAAn),' '), function(x) x[1]),
#           varAA=sapply(strsplit(names(singleAAn),' '), function(x) x[2]),
#           Instances=singleAAn)
#singleAAdfn=singleAAdfn[-44,]
#ggplot(singleAAdfn, aes(x=refAA, y=varAA, fill=Instances.Freq))+geom_tile()+scale_fill_viridis_c()
#
#
#
##loi= as.vector(sapply(coding.effects$PROTEINLOC, function(x)x[1])) + 
##
##
##lapply(coding.effects$REFAA, as.matrix)
##varAA=lapply(cooding.effects$VARAA, as.matrix)
##pdiff=which(refAA!=varAA)
##loi=sapply(coding.effects$PROTEINLOC, function(x) x[1])+pdiff-1
##             varoi=varAA[pdiff] 
##
##saveRDS(coding.effects[,-c(4,5)] , file = '/data/yeast/oligos/onlyPams/chrI.RDS')
##
#
#
#
#
#
## lookup provean scores (needs code optimization)
#proveanEffects =get_provean_scores(coding.effects, paste0(working.dir, 'ref/provean_scores.RData'))
### save(proveanEffects, file=paste0(working.dir, 'results/editor_provean_effects.RData') )
##load(paste0(working.dir, 'reference/editor_provean_effects.RData'))
#coding.effects$proveanEffects=proveanEffects
#
#
#
##add in column for absolute max value of provean effects 
#coding.effects=get_max_absolute_provean(coding.effects)
#
#
#
## save as RDS objects
#saveRDS(coding.effects, file=paste0(working.dir, 'results/codingEffects_chrI.RDS')) 
#saveRDS(guides, file=paste0(working.dir, 'results/guideTable.RDS'))
#saveRDS(guidesNoCs, file=paste0(working.dir, 'results/guideTable_noEdits.RDS'))
