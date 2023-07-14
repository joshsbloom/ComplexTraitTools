library(MutationalPatterns)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(seqinr)
library(GenomicFeatures)
#library(RFLPtools)
#library(intervals)
#library(VariantAnnotation)
    
#where is code
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
#https://cdn.elifesciences.org/articles/62669/elife-62669-supp6-v2.vcf.gz

vcf.file='~/Downloads/elife-62669-supp6-v2.vcf'

#convert VCF into GenomicRanges object
rmg=read_vcfs_as_granges(vcf.file, 'RMx', sacCer3, type='all')$RMx

chunksize=112

#extract coordinates for all promoters 
proms=promoters(genes(txdb), upstream=chunksize, downstream=0)

#find promoters with variants
pvar=findOverlaps(proms, rmg)

#remake promoters object for only the ones with variants
proms=proms[unique(queryHits(pvar))]

#remake object intersecting promoters with variants so indexing is less annoying
pvar=findOverlaps(proms, rmg)

#split variant indices by promoters 
vpg= split(subjectHits(pvar), queryHits(pvar))


#vectorized retrieval of promoter sequences
refSeq=getSeq(sacCer3, as.character(seqnames(proms)), 
              start=start(proms), 
              end=end(proms),strand='+') #as.character(strand(proms)))

# iterate through each promoter
#refs=list()
refSeqs=list()
altSeqs=list()
for(i in 1:length(vpg)){
    print(i) #i=2595

    #get the variants 
    test=rmg[vpg[[i]]]
  
    #when ranges in the vcf overlap we're hosed, ignore genes with such cases for now
    if(isDisjoint(test)) {
    
    #put into chunksize coordinates
    nstart=start(ranges(test))-start(proms[i])+1
    nend=end(ranges(test))-start(proms[i])+1
    
    #deal with annoying edge cases (literally edge here)
    hangsover=nend>chunksize
    hangsunder=min(nstart)<1
    
    #drop these
    if(hangsunder) {next;}

    #pad ref for these 
    if(sum(hangsover)>0 ) {
        nendpad=max(nend)-chunksize
        expref= getSeq(sacCer3, as.character(seqnames(proms[i])), 
              start=start(proms[i]), 
              end=end(proms[i])+nendpad,strand='+')   
    } else {
        expref=refSeq[i][[1]]
    }

    refSeqs[[proms[i]$gene_id]]=expref
    #the unlist business is to deal with cases where more than one alternate allele is listed in the vcf
    altSeqs[[proms[i]$gene_id]]=replaceAt(expref, at=IRanges(start=nstart, end=nend),
                                      value= unlist(DNAStringSetList(sapply(test$ALT, function(x) x[1]))))
    }
                                                                                  #value=unlist(test$ALT))
}

#we lose about a dozen due to annoying edge cases 
final.BY.W=DNAStringSet(refSeqs)
final.RM.W=DNAStringSet(altSeqs)

#go back and reverse complement when promoter is on minus strand 
final.BY=final.BY.W
final.RM=final.RM.W
mstrandcases=as.character(strand(proms[names(final.BY.W)]))=='-'
final.BY[mstrandcases]=reverseComplement(final.BY[mstrandcases])
final.RM[mstrandcases]=reverseComplement(final.RM[mstrandcases])
writeXStringSet(final.BY, file='~/Desktop/BYprom.fasta')
writeXStringSet(final.RM, file='~/Desktop/RMprom.fasta')

