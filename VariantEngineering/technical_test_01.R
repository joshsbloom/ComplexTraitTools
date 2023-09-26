library(Biostrings)
library(GenomicFeatures)
library(seqinr)
library('TxDb.Scerevisiae.UCSC.sacCer3.sgdGene')
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(SimRAD)

source('/home/jbloom/Dropbox/code/ComplexTraitTools/VariantEngineering/R/designCRISPR_gRNAs_fx.R')

#pre-annotated guides
g=readRDS('/data/yeast/yeast_oligos/onlyPams/tables/guides.RDS')

#pre-annotated coding variants 
#ce=readRDS('/data/yeast/yeast_oligos/onlyPams/tables/codingEffects.RDS')
#ce$nedit=nchar(data.table::tstrsplit(ce$gEdit, ':')[[3]])
#dupce=duplicated(ce$gEdit)
#ce=ce[!dupce,]
#saveRDS(ce, file = '/data/yeast/yeast_oligos/onlyPams/tables/codingEffects_dedup.RDS')
ce=readRDS('/data/yeast/yeast_oligos/onlyPams/tables/codingEffects_dedup.RDS')


# ---extract introns -------------------------------------------------------------------------------
ibt=intronsByTranscript(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, use.names=T)

gwi=sapply(ibt, function(x) length(x))

ibt=unlist(ibt[gwi>0])
ibt=ibt[seqnames(ibt)!='chrM',]
ibt=ibt[width(ibt)>150,]

ibt_p=narrow(ibt[strand(ibt)=='+'], start=10, end=-45)
ibt_m=narrow(ibt[strand(ibt)=='-'], start=45, end=-10)

guides_in_introns.index=findOverlaps(g$X, c(ibt_p,ibt_m), type='within')
guides_in_introns=g[queryHits(guides_in_introns.index),]
guides_in_introns=guides_in_introns[!guides_in_introns$U6terminator,]
guides_in_introns=guides_in_introns[!(guides_in_introns$guidePAM12PMcount>2),]
#------------------------------------------------------------------------------------introns------



# ---construnct randomer sequence, guides and repair templates -----------------------------------
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3

set.seed(1)
rndmer=sim.DNAseq(2000, GCfreq=0.38)
#write.fasta(rndmer, names='rndmer', file='~/Downloads/randomer_bsgenome/rndmer.fa')
rndmer=DNAString(rndmer)
write.fasta(rndmer, names='rndmer', file='/data/yeast/yeast_oligos/Shen/rndmer_2kb.fa')


unique.chrs=paste0('chr', as.roman(1:16))
sc3.set=buildSacCer3_genome_dictionaryFR(sacCer3, unique.chrs)

guides_in_rndmer=get_gRNAs_and_predict_string(rndmer, sc3.set, 'rndmer')
guides_in_rndmer=augment_guides(guides_in_rndmer)
guides_in_rndmer=guides_in_rndmer[!guides_in_rndmer$U6terminator,]
guides_in_rndmer=guides_in_rndmer[!(guides_in_rndmer$guidePAM12PMcount>2),]

#nchar(data.table::tstrsplit(unlist(as.vector(guides_in_rndmer$gEditInfo[1,])), ":")[[3]])
#5 1
#15 2
#36 3 
#--------------------------------------------------------------------------------------------------




# construct GRanges objects for UTRs 

# For project 2 , start here ---------------------------------------------------
makeGR_from_fasta_header=function(fasta.in) {
    f=read.fasta(fasta.in)

    f.annot=sapply(f, function(x) attr(x, 'Annot'))
    f.annot.s=data.table::tstrsplit(f.annot, '_')
    #5th element is yeast ORF
    #9th element has all the other goodies

    f.annot.s2=data.table::tstrsplit(f.annot.s[[9]], ' ')


    f.annot.s3=data.table::tstrsplit(f.annot.s2[[2]], '=|:|-')
    #2 is chr, #3 is start, #4 is end

    f.annot.strand=gsub('strand=', '',f.annot.s2[[5]])

    chr=f.annot.s3[[2]]
    #reorder
    levels(chr)=paste0('chr', as.roman(1:16))
    f.GR=makeGRangesFromDataFrame(data.frame(chr=chr, 
                  start=as.numeric(f.annot.s3[[3]]),
                  end= as.numeric(f.annot.s3[[4]]), 
                  strand= f.annot.strand))
    f.GR$gene_id=f.annot.s[[5]]
    #annoying roman numerals 
    f.GR=f.GR[order(seqnames(f.GR), start(f.GR)),]
    return(f.GR)
}


#parse coordinates in SGD UTR annotation files
fasta.in.5utr='/home/jbloom/Dropbox/code/ComplexTraitTools/VariantAnnotation/yeast_ref/SGD_all_ORFs_5prime_UTRs.fsa'
fasta.in.3utr='/home/jbloom/Dropbox/code/ComplexTraitTools/VariantAnnotation/yeast_ref/SGD_all_ORFs_3prime_UTRs.fsa'

utr5.GR=makeGR_from_fasta_header(fasta.in.5utr)
utr3.GR=makeGR_from_fasta_header(fasta.in.3utr)


#deal with overlapping 5' utrs
utr5.GRs=split(utr5.GR, utr5.GR$gene_id)
#oy sometimes there are more than one annotation per gene, extend to longest possible range, collapsing overlap
utr5.GRsr=sapply(utr5.GRs, reduce)
#annoying reformatting
utr5.GRsr=do.call('rbind', lapply(utr5.GRsr, as.data.frame))
#and back into a Granges object
utr5.GR=makeGRangesFromDataFrame(utr5.GRsr)
utr5.GR$gene_id=names(utr5.GR)
utr5.GR=utr5.GR[order(seqnames(utr5.GR), start(utr5.GR)),]


#deal with overlapping 3' utrs
utr3.GRs=split(utr3.GR, utr3.GR$gene_id)
#oy sometimes there are more than one annotation per gene, extend to longest possible range, collapsing overlap
utr3.GRsr=sapply(utr3.GRs, reduce)
#annoying reformatting
utr3.GRsr=do.call('rbind', lapply(utr3.GRsr, as.data.frame))
#and back into a Granges object
utr3.GR=makeGRangesFromDataFrame(utr3.GRsr)
utr3.GR$gene_id=names(utr3.GR)
utr3.GR=utr3.GR[order(seqnames(utr3.GR), start(utr3.GR)),]

#-------------------------------------------------------------------------------------------------



# Now get the intergenic regions between 3' UTRs --------------------------------------------------------------------------


# extend genes by their 3' UTRs and find genes with 3'UTRs pointing toward each other 
base.dir='~/Dropbox/code/ComplexTraitTools/'
#source(paste0(base.dir, 'VariantAnnotation/R/proveanFx.R'))

sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
    # make our own txdb structure (motivated by annotation discrepancy found for GRX3 ORF and issues with dubious ORFs)
    # pulled from http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz and gunzipped it (beware of changes)


ref.data.dir=paste0(base.dir, 'VariantAnnotation/yeast_ref/') 

# see here:
#https://github.com/joshsbloom/ComplexTraitTools/tree/main/VariantAnnotation/yeast_ref
txdb=GenomicFeatures::makeTxDbFromGFF(paste0(base.dir,'VariantAnnotation/yeast_ref/saccharomyces_cerevisiae.gff.gz'))

#download vcf from here 

#td=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
#gtd=genes(td)
gtd=genes(txdb)

#follow( gtd[strand(gtd)=='+'], gtd)


gtds=split(gtd, gtd$gene_id)
gtd=gtd[seqnames(gtd)!='chrM']

utr3.GRs=split(utr3.GR, utr3.GR$gene_id)
nn=names(gtds)[names(gtds) %in% names(utr3.GRs)]
nno=names(gtds)[!(names(gtds) %in% names(utr3.GRs))]
#n=nn[1]

#add 3' utr to gene 
gList=list()
for(n in nn) {
    print(n)
    w3utr=reduce(unlist(c(gtds[n], utr3.GRs[n])))
    w3utr$gene_id=n
    gList[[n]]=w3utr
}
for(n in nno) {
    gList[[n]]=gtds[[n]]
}

gList2=unlist(GRangesList(unlist(gList)))
gtd=gList2


gtd_pos=gtd[strand(gtd)=='+']
gtd_min=gtd[strand(gtd)=='-']

f=precede(gtd_pos, gtd_min, ignore.strand=T)
ff=precede(gtd_pos, gtd_pos)

uf3u=gtd_pos[which(!is.na(f)&!is.na(ff))]
uf3dp=gtd_pos[ff[!is.na(f)&!is.na(ff)]]
uf3d=gtd_min[f[!is.na(f)&!is.na(ff)]]

point5=uf3u[!(start(uf3dp)<start(uf3d))]
point3=uf3d[!(start(uf3dp)<start(uf3d))]
#median distance between genes that 

gdf=data.frame(chr=seqnames(point5), start=end(point5)+1, end=start(point3)-1,
           gene_id1=point5$gene_id, gene_id2=point3$gene_id)

gdfGR=makeGRangesFromDataFrame(gdf)
gdfGR$gene_id1=point5$gene_id
gdfGR$gene_id2=point3$gene_id

intergenic3=gdfGR[which(!(1:length(gdfGR)) %in% unique(queryHits(findOverlaps(gdfGR, gtd))) )]


intergenic3.300=intergenic3[width(intergenic3)>300]
guides_in_intergenic.index=findOverlaps(g$X, intergenic3.300, type='within')
guides_in_intergenic=g[queryHits(guides_in_intergenic.index),]
guides_in_intergenic=guides_in_intergenic[!guides_in_intergenic$U6terminator,]
guides_in_intergenic=guides_in_intergenic[!(guides_in_intergenic$guidePAM12PMcount>2),]
guides_in_intergenic$guide=as.character(subseq(DNAStringSet(guides_in_intergenic$guide), 1,20))

#---------------------------------------------------------------------------------------------------------------




# Shen gene sett -----------------------------------------------------------------------------------------------

#https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-11-260
#shen set
geneset=c('YDR448W',
          'YMR116C',
          'YOR198C',
          'YCR047C',
          'YLR110C',
          'YNL080C',
          'YGL020C',
          'YML094W',
          'YEL044W',
          'YJL124C',
          'YHL011C',
          'YGL058W',
          'YFR032C-A',
          'YJL189W',
          'YOR096W',
          'YHL025W',
          'YLR435W',
          'YGR020C',
          'YGR105W')



# synonymous PAM  breakers

#egene.set=na.omit(prev[prev[,3]<.5,])$GENEID
y=ce[ce$GENEID %in% geneset,]
#y=x[!grepl('\\(', x$AAeffect),]
y=y[y$CONSEQUENCE == 'synonymous',]
dds=y[!duplicated(y$repairTemplate),]
dds=dds[!dds$U6terminator,]
dds=dds[!(dds$guidePAM12PMcount>2),]

h=g[match(dds$guideIndex, g$guideIndex),]
dds$guide=as.character(subseq(DNAStringSet(h$guide), 1,20))

#non-synonymous PAM breakers 
y=ce[ce$GENEID %in% geneset,]
#y=x[!grepl('\\(', x$AAeffect),]
y=y[y$CONSEQUENCE == 'nonsynonymous',]
ddn=y[!duplicated(y$repairTemplate),]
ddn=ddn[ddn$guideIndex %in% dds$guideIndex,]
ddn=ddn[!grepl('\\(', ddn$AAeffect),]
ddn=ddn[!is.na(ddn$proveanEffects),]
ddn=ddn[!ddn$U6terminator,]
ddn=ddn[!(ddn$guidePAM12PMcount>2),]
h=g[match(ddn$guideIndex, g$guideIndex),]
ddn$guide=as.character(subseq(DNAStringSet(h$guide), 1,20))
#-----------------------------------------------------------------------------------------------------------------












# Pseudogenes ----------------------------------------------------------------------------------------------------
#https://wiki.yeastgenome.org/index.php/Category:pseudogene
#pseudogenes 
pseudo=rbind(c('YAR061W', 'chrI',   218140, 219145, '+'),
             c('YAR062W', 'chrI',   214888, 224046, '+'),
             c('YIL167W', 'chrIX',  29032, 30048, '+'),
             c('YIL174W', 'chrIX',  9469, 9696, '+'),
             c('YIL175W', 'chrIX',  9183, 9500, '+'),
             c('YIR043C', 'chrIX',  437043, 438179, '-'),
             c('YOL153C', 'chrXV',  36822, 38567, '-'),
             c('YLL016W', 'chrXII', 112847, 115993, '+'),
             c('YLL017W', 'chrXII', 112234, 112545, '+'))
         
             
#c('uURA3',   'chrV',   115944, 116165, '+'),
#c('uLEU2',   'chrIII', 91000, 91323, '+'),
# 'YCL074W',
# 'YCL075W',
# 'YDR134C',
# 'YIL080W',
# 'YIL168W',
# 'YIL170W',
# 'YIL171W',
# 'YIR044C',
# 'YNL054W-B',
# 'YPL060C-A',
# 'YPL275W',
# 'YPL276W')

pseudo=data.frame(pseudo)
pseudo[,3]=as.numeric(pseudo[,3])
pseudo[,4]=as.numeric(pseudo[,4])
names(pseudo)=c('gene', 'chr', 'start', 'end', 'strand')
pseudoDf=makeGRangesFromDataFrame(pseudo)

pseudoGuides=g[unique(queryHits(findOverlaps(g$X, pseudoDf, type='within'))),]
pseudoGuides=pseudoGuides[!pseudoGuides$U6terminator,]
pseudoGuides=pseudoGuides[!(pseudoGuides$guidePAM12PMcount>2),]
pseudoGuides$guide=as.character(subseq(DNAStringSet(pseudoGuides$guide), 1,20))

#see gEdit to eval number of variants
#prev=readxl::read_excel('/data/yeast/yeast_oligos/estopsV2/NIHMS945384-supplement-Sup_Table_12.xls')
#854 genes previously annotated as essential plus observed as essential in our PTC experiment 




# basic idea for breakdown 

#rndmer 
#1bp and 2bp
#530+1590
#2120


#intergenic
#sample 500 

#introns
#sample 500 

# pseudogenes 
# sample 500

# synonymous 
# full Shen set
# 1011

# non-synonymous Shen
# sample 2000

# non-synon high provean (<-10)
# 500
# non-synon med provean  (-10 to -5)
# 500
# non-synon essential stops
# 500

# + random synonymous 
# 1000

# + randon non-synonymous


stack.it=function(gtable) {
    #fix guide sequence b

    gtable$guide=as.character(subseq(DNAStringSet(gtable$guide), 1,20))
    #was this 
    #rstack=rep(gtable, each=ncol(gtable$repairTemplate[1,]))
    rstack=rep(gtable, ncol(gtable$repairTemplate[1,]))
    rstack$repairTemplate=stack(gtable$repairTemplate)$values
    rstack$editedSequenceFwdL=stack(gtable$editedSequenceFwdL)
    rstack$gEditInfo=stack(gtable$gEditInfo)
    rstack$nedit=nchar(data.table::tstrsplit(rstack$gEditInfo[,1], ':')[[3]])
    return(rstack)
}


o.array=list()

#ncol blah blah is 56 for all possible pam breakers 
rstack=stack.it(guides_in_rndmer)
o.array[['rndmer']]=DataFrame(rstack[rstack$nedit<3,])
   

istack=stack.it(guides_in_intergenic)
istack=istack[istack$nedit==1,]
set.seed(50)
istack=istack[sort(sample(1:nrow(istack), 500)),]

o.array[['intergenic']]=istack


istack2=stack.it(guides_in_introns)
istack2=istack2[istack2$nedit==1,]
set.seed(51)
istack2=istack2[sort(sample(1:nrow(istack2), 500)),]

o.array[['intron']]=istack2


pstack=stack.it(pseudoGuides)
pstack=pstack[pstack$nedit==1,]
set.seed(499)
pstack=pstack[sort(sample(1:nrow(pstack), 500)),]
o.array[['pseudo']]=pstack


c.array=list()
c.array[['shen_syn']]=dds
set.seed(53)
c.array[['shen_nonsyn']]=ddn[sort(sample(1:nrow(ddn), 2000)),]



#subset provean
nv.change=sapply(ce$proveanEffects, function(x) length(x))
ce.sub=ce[nv.change==1,]
ce.sub=ce.sub[!ce.sub$dubious & !ce.sub$U6terminator & !(ce.sub$guidePAM12PMcount>2) & !(ce.sub$GENEID %in% geneset),  ]
uce=unlist(ce.sub$proveanEffects)
#strong provean
#ce.sub[na.omit(which(uce<(-10))),]

set.seed(54)
strong.provean=ce.sub[na.omit(which(uce<(-10))),]
c.array[['strong_provean']]=strong.provean[sort(sample(1:nrow(strong.provean), 500)),]
h=g[match(c.array[['strong_provean']]$guideIndex, g$guideIndex),]
c.array[['strong_provean']]$guide=as.character(subseq(DNAStringSet(h$guide), 1,20))



set.seed(55)
med.provean=ce.sub[na.omit(which(uce>(-10)&uce<(-5))),]
c.array[['med_provean']]=med.provean[sort(sample(1:nrow(med.provean), 500)),]
h=g[match(c.array[['med_provean']]$guideIndex, g$guideIndex),]
c.array[['med_provean']]$guide=as.character(subseq(DNAStringSet(h$guide), 1,20))


set.seed(56)
#first half of essential stops 
y=ce[is.na(ce$proveanEffects) & ce$CONSEQUENCE=='nonsynonymous' & ce$essential==T &ce$U6terminator==F & ce$dubious==F & ce$guidePAM12PMcount<3,]
psplit=split(y$POSAA.delta, y$GENEID)
pmax=sapply(psplit, function(x) max(unlist(x)))

#y=y[y$VARAA.delta=='*',]
firsthalf=round(pmax/2)
y$firsthalf=firsthalf[y$GENEID]
y=y[y$POSAA.delta<y$firsthalf,]
y=y[,-ncol(y)]
y=y[!duplicated(y),]

c.array[['estops']]=y[sort(sample(1:nrow(y), 500)),]
h=g[match(c.array[['estops']]$guideIndex, g$guideIndex),]
c.array[['estops']]$guide=as.character(subseq(DNAStringSet(h$guide), 1,20))


ce1=ce[ce$nedit==1,]
set.seed(577)
y=ce1[ ce1$CONSEQUENCE=='synonymous' & ce1$U6terminator==F & ce1$dubious==F & ce1$guidePAM12PMcount<3 & !(ce1$GENEID %in% geneset) ,]
c.array[['syn_random']]=y[sort(sample(1:nrow(y), 1000)),]
h=g[match(c.array[['syn_random']]$guideIndex, g$guideIndex),]
c.array[['syn_random']]$guide=as.character(subseq(DNAStringSet(h$guide), 1,20))


set.seed(58)
y=ce1[ ce1$CONSEQUENCE=='nonsynonymous' &ce1$U6terminator==F & ce1$dubious==F & ce1$guidePAM12PMcount<3 & !(ce1$GENEID %in% geneset),]
c.array[['nonsyn_random']]=y[sort(sample(1:nrow(y), 1000)),]
h=g[match(c.array[['nonsyn_random']]$guideIndex, g$guideIndex),]
c.array[['nonsyn_random']]$guide=as.character(subseq(DNAStringSet(h$guide), 1,20))

saveRDS(o.array, file = '/data/yeast/yeast_oligos/Shen/o_array.RDS')
saveRDS(c.array, file = '/data/yeast/yeast_oligos/Shen/c_array.RDS')
#--------------------------------------------------------------------------------------------------------------------------------------

wd='/data/yeast/yeast_oligos/Shen/'
o.array=readRDS(paste0(wd, 'o_array.RDS'))
c.array=readRDS(paste0(wd, 'c_array.RDS'))

guides=c(do.call('c', sapply(o.array, function(x) x$guide)), do.call('c', sapply(c.array, function(x) x$guide)))
repTemps=c(do.call('c', sapply(o.array, function(x) x$repairTemplate)), do.call('c', sapply(c.array, function(x) x$repairTemplate)))



PBS_gRNA1=toupper('gcagtgaaagataaatgatc')
struct_region_terminator=toupper('gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtggtgcTTTTTTTGTTTTTTATGTCT')
rightside_can1_1=toupper('TTCCCGACGAGAGTAAATGGCGAGGATACGTTCTCTATGG')
rightside_can1_2='TTCCCGACGAGAGTAAATGG'

#pick 1
#final.oligo=paste0(PBS_gRNA1, y$guide, struct_region_terminator, y$repairTemplate, rightside_can1_1)
final.oligo=paste0(PBS_gRNA1, guides, struct_region_terminator, repTemps, rightside_can1_2)

final.oligoD= DNAStringSet(final.oligo)
Acnt=letterFrequency(final.oligoD, 'A')
Tcnt=letterFrequency(final.oligoD, 'T')
final.oligoDrc=final.oligoD
for(i in 1:length(final.oligoD) ){
    print(i)
    if(Acnt[i]>Tcnt[i]) final.oligoDrc[i]=reverseComplement(final.oligoDrc[i]) 
}

write.table(data.frame(oligo=as.character(final.oligoDrc)),
                       file='/data/yeast/yeast_oligos/Shen/syn_nonsyn_10k.txt', sep='\t', quote=F, row.names=F)




#
#length(unique(ddy$guideIndex))
#
#ddy[order(ddy$guideIndex),]
#edist.check=rep(NA, nrow(ddy))
#for(i in 1:nrow(ddy)){
#    print(i)
#    edist.check[i]=sum(stringdist(ddy$repairTemplate[i], ddy$repairTemplate[-i], method='lv')<3)
#}
#
#
#
#
#
#y=x[x$GENEID=='YAR015W' | x$GENEID=='YOR128C',]
#y=y[is.na(y$proveanEffects) & y$CONSEQUENCE=='nonsynonymous',]
#y=y[!duplicated(paste(y$guideIndex,y$repairTemplate)),]
#
##manually get first half of each gene   
#filtered=y[y$guideGCcontent>.4 & y$guideGCcontent<.6 & y$U6terminator==F & y$guidePAM12PMcount<3,]
#set.seed(100)
#isample=sort(sample(1:nrow(filtered), 500))
#
#
#estops500=filtered[isample,]
#h=g[match(estops500$guideIndex, g$guideIndex),]
#estops500$guide=as.character(subseq(DNAStringSet(h$guide), 1,20))
#saveRDS(estops500, '/data/yeast/yeast_oligos/estopsV2/estops500.RDS')
#

y=estops384
PBS_gRNA1=toupper('gcagtgaaagataaatgatc')
struct_region_terminator=toupper('gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtggtgcTTTTTTTGTTTTTTATGTCT')
rightside_can1_1=toupper('TTCCCGACGAGAGTAAATGGCGAGGATACGTTCTCTATGG')
rightside_can1_2='TTCCCGACGAGAGTAAATGG'

#pick 1
#final.oligo=paste0(PBS_gRNA1, y$guide, struct_region_terminator, y$repairTemplate, rightside_can1_1)
final.oligo=paste0(PBS_gRNA1, y$guide, struct_region_terminator, y$repairTemplate, rightside_can1_2)

final.oligoD= DNAStringSet(final.oligo)
Acnt=letterFrequency(final.oligoD, 'A')
Tcnt=letterFrequency(final.oligoD, 'T')
final.oligoDrc=final.oligoD
for(i in 1:length(final.oligoD) ){
    print(i)
    if(Acnt[i]>Tcnt[i]) final.oligoDrc[i]=reverseComplement(final.oligoDrc[i]) 
}

write.table(data.frame(oligo=as.character(final.oligoDrc)),
                       file='/data/yeast/yeast_oligos/onlyPams/tables/ade_opool_pilot_can1_2.txt', sep='\t', quote=F, row.names=F)








#automated, get the max POSAA.delta for each gene and then div/2 and logic it 

y=y[(y$GENEID=='YAR015W' &  (as.vector(unlist(y$POSAA.delta))<152)) |
    (y$GENEID=='YOR128C'&   (as.vector(unlist(y$POSAA.delta))<271)),    ]

#sync up with the guide sequences 

#
h=g[match(y$guideIndex, g$guideIndex),]
y$guide=as.character(subseq(DNAStringSet(h$guide), 1,20))

saveRDS(y, '/data/yeast/yeast_oligos/onlyPams/tables/first_half_of_ade1_2_stops.RDS')

y=readRDS('/data/yeast/yeast_oligos/onlyPams/tables/first_half_of_ade1_2_stops.RDS')

PBS_gRNA1=toupper('gcagtgaaagataaatgatc')
struct_region_terminator=toupper('gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtggtgcTTTTTTTGTTTTTTATGTCT')
rightside_can1_1=toupper('TTCCCGACGAGAGTAAATGGCGAGGATACGTTCTCTATGG')
rightside_can1_2='TTCCCGACGAGAGTAAATGG'

#pick 1
#final.oligo=paste0(PBS_gRNA1, y$guide, struct_region_terminator, y$repairTemplate, rightside_can1_1)
final.oligo=paste0(PBS_gRNA1, y$guide, struct_region_terminator, y$repairTemplate, rightside_can1_2)

final.oligoD= DNAStringSet(final.oligo)
Acnt=letterFrequency(final.oligoD, 'A')
Tcnt=letterFrequency(final.oligoD, 'T')
final.oligoDrc=final.oligoD
for(i in 1:length(final.oligoD) ){
    print(i)
    if(Acnt[i]>Tcnt[i]) final.oligoDrc[i]=reverseComplement(final.oligoDrc[i]) 
}

write.table(data.frame(oligo=as.character(final.oligoDrc)),
                       file='/data/yeast/yeast_oligos/onlyPams/tables/ade_opool_pilot_can1_2.txt', sep='\t', quote=F, row.names=F)



#zz=split(y, y$GENEID)
#par(mfrow=c(2,1))
#plot(start(zz[[1]]$X), zz[[1]]$POSAA.delta)
#plot(start(zz[[2]]$X), zz[[2]]$POSAA.delta)

# quick take on recurring guides

gsubseq=subseq(g$guide, start=1, end=20)
h=rle(sort(gsubseq))

gmax40=g[gsubseq%in% h$values[h$lengths>40], ] #="GCAAGGATTGATAATGTAAT",]



