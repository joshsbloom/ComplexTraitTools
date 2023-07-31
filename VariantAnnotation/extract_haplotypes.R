library(MutationalPatterns)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(seqinr)
library(GenomicFeatures)
library(data.table)
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



#given a Granges object of sequence ranges of interest and GRanges object constructed from a VCF
makeVarSeq=function(proms, rmg) {

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
       
        chunksize=width(proms[i]) 
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
    #return a list of the two parents sequences, push variants into second parent sequence 
    return(list(P1=final.BY, P2=final.RM))
}


# Specifically for Proj 3==========================================
#this is just specfic to extracting the first 112 bp for proj 3
chunksize=112

#extract coordinates for all promoters 
proms=promoters(genes(txdb), upstream=chunksize, downstream=0)

VSout=makeVarSeq(proms, rmg)

writeXStringSet(VSout$P1, file='~/Desktop/BYprom.fasta')
writeXStringSet(VSout$P2, file='~/Desktop/RMprom.fasta')
#=============================================================================



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






#construct Prom sequence 
#separate logic based on strand--------------------------------------------
    #+ strand
    utr5.GRf=utr5.GR[strand(utr5.GR)=='+']
    utr5.GRf$gene_id=names(utr5.GRf)
    #this seems to be the function we want, confusing documentation
    pg=follow(utr5.GRf, unstrand(genes(txdb)), ignore.strand=T)
    #yuck hack fix this
    pg[is.na(pg)]=1
    ug=genes(txdb)[pg]
    dist.up.check=(start(utr5.GRf)-end(ug))
    dist.up.check[dist.up.check>1000]=1000
    dist.up.check[dist.up.check<1]=1000


    igF=makeGRangesFromDataFrame(
            data.frame(chr=as.character(seqnames(utr5.GRf)),
               start=start(utr5.GRf)-dist.up.check,
               end=start(utr5.GRf)-1, strand='+'))
    names(igF)=utr5.GRf$gene_id
    igF$gene_id=utr5.GRf$gene_id
    #-----------------------------------------------------------------------------

    #-strand
    utr5.GRr=utr5.GR[strand(utr5.GR)=='-']
    utr5.GRr$gene_id=names(utr5.GRr)

    pg=precede(utr5.GRr, unstrand(genes(txdb)), ignore.strand=T)
    #yuck hack fix this 
    pg[is.na(pg)]=1
    ug=genes(txdb)[pg]
    dist.up.check=start(ug) - (end(utr5.GRr))
    dist.up.check[dist.up.check>1000]=1000
    dist.up.check[dist.up.check<1]=1000
    igR=makeGRangesFromDataFrame(
            data.frame(chr=as.character(seqnames(utr5.GRr)),
               start=end(utr5.GRr)+1,
               end=end(utr5.GRr)+dist.up.check, strand='-'))
    names(igR)=utr5.GRr$gene_id
    igR$gene_id=utr5.GRr$gene_id
Prom=c(igF,igR)
Prom=Prom[order(seqnames(Prom), start(Prom)),]
#---------------------------------------------------------------------------------------


#construct Term sequence 
#separate logic based on strand--------------------------------------------
#+ strand
utr3.GRf=utr3.GR[strand(utr3.GR)=='+']
utr3.GRf$gene_id=names(utr3.GRf)
#this seems to be the function we want, confusing documentation
pg=precede(utr3.GRf, unstrand(genes(txdb)), ignore.strand=T)
pg[is.na(pg)]=1
ug=genes(txdb)[pg]
dist.up.check=start(ug) - (end(utr3.GRf))
dist.up.check[dist.up.check>500]=500
dist.up.check[dist.up.check<1]=500
termF=makeGRangesFromDataFrame(
        data.frame(chr=as.character(seqnames(utr3.GRf)),
           start=end(utr3.GRf)+1,
           end=end(utr3.GRf)+dist.up.check, strand='+'))
names(termF)=utr3.GRf$gene_id
termF$gene_id=utr3.GRf$gene_id

#- strand
utr3.GRr=utr3.GR[strand(utr3.GR)=='-']
utr3.GRr$gene_id=names(utr3.GRr)

pg=follow(utr3.GRr, unstrand(genes(txdb)), ignore.strand=T)
#yuck hack fix this
pg[is.na(pg)]=1
ug=genes(txdb)[pg]
dist.up.check=(start(utr3.GRr)-end(ug))
dist.up.check[dist.up.check>1000]=500
dist.up.check[dist.up.check<1]=500
termR=makeGRangesFromDataFrame(
        data.frame(chr=as.character(seqnames(utr3.GRr)),
           start=start(utr3.GRr)-dist.up.check,
           end=start(utr3.GRr)-1, strand='-'))
names(termR)=utr3.GRr$gene_id
termR$gene_id=utr3.GRr$gene_id

Term=c(termF,termR)
Term=Term[order(seqnames(Term), start(Term)),]



#all sequences are 5' - 3' relative to CDS

#These are based on the experimentally defined 5' UTR start site amd nearest upstream CDS end
PromSeq=getSeq(sacCer3, as.character(seqnames(Prom)), 
                  start=start(Prom), 
                  end=end(Prom),strand=strand(Prom) )
names(PromSeq)=names(Prom)

#These are based on the experimentally defined 5' UTR lengths
# if you want to clip to 300 bp then take the 300bp at the end of the sequence 
# if you clip then make sure you append the bit clipped off into the end of Prom and
UTR5Seq=getSeq(sacCer3, as.character(seqnames(utr5.GR)), 
                  start=start(utr5.GR), 
                  end=end(utr5.GR),strand=strand(utr5.GR) )
names(UTR5Seq)=names(UTR5Seq)


#These are based on the experimentally defined 3' UTR lengths
# if you want to clip to 350 bp then take the 350bp at the start of the sequence  
UTR3Seq=getSeq(sacCer3, as.character(seqnames(utr3.GR)), 
                  start=start(utr3.GR), 
                  end=end(utr3.GR),strand=strand(utr3.GR) )
names(UTR3Seq)=names(UTR3Seq)
# if you clip then make sure you append the bit clipped off into the start of Term

#These are based on the experimentally defined 3' UTR end site and nearest  downstream CDS start
TermSeq=getSeq(sacCer3, as.character(seqnames(Term)), 
                  start=start(Term), 
                  end=end(Term),strand=strand(Term) )
names(TermSeq)=names(TermSeq)

writeXStringSet(PromSeq, file='~/BY_Prom.fasta')
writeXStringSet(UTR5Seq, file='~/BY_5UTR.fasta')
writeXStringSet(UTR3Seq, file='~/BY_3UTR.fasta')
writeXStringSet(TermSeq, file='~/BY_Term.fasta')



PromVar=makeVarSeq(Prom, rmg)
UTR5Var=makeVarSeq(utr5.GR, rmg)
UTR3Var=makeVarSeq(utr3.GR, rmg)
TermVar=makeVarSeq(Term, rmg)
writeXStringSet(PromVar$P2, file='~/RM_Prom.fasta')
writeXStringSet(UTR5Var$P2, file='~/RM_5UTR.fasta')
writeXStringSet(UTR3Var$P2, file='~/RM_3UTR.fasta')
writeXStringSet(TermVar$P2, file='~/RM_Term.fasta')
