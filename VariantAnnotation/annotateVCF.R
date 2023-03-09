library("BSgenome.Scerevisiae.UCSC.sacCer3")
#library(RFLPtools)
#library(intervals)

library(seqinr)
library(VariantAnnotation)
library(GenomicFeatures)
library(vcfR)

# Global Variables
#working.dir='/data/yeast/yeast_oligos/'


# I split the previously curated vcf (variant call format files with genotype information) into sub vcfs for each chromosome here 
vcf.prefix='/data/yeast/chr_rename-014.'

#where are vcf.gz files split by chr
#vcf.prefix='/home/jbloom/Downloads/'

#where to save output
out.dir='/home/jbloom/Downloads/peter2018_vcf/'

#where is code
base.dir='~/Dropbox/code/ComplexTraitTools/'

source(paste0(base.dir, 'VariantAnnotation/R/proveanFx.R'))

## note, the bioconductor TxDb object is old and broken, we'll go ahead and build our own
##library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
##txdb=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
# make our own txdb structure (motivated by annotation discrepancy found for GRX3 ORF and issues with dubious ORFs)
# pulled from http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz and gunzipped it (beware of changes)
#see VariantAnnotation/yeast_ref
ref.data.dir=paste0(base.dir, 'VariantAnnotation/yeast_ref/') 


#'/data/yeast/yeast_oligos/ref/'
txdb=GenomicFeatures::makeTxDbFromGFF(paste0(base.dir,'VariantAnnotation/yeast_ref/saccharomyces_cerevisiae.gff.gz'))

provean.score.file=paste0(ref.data.dir,'provean_scores.RDS')
provean_scores=readRDS(provean.score.file)
meltedProveanScores=meltProvean(provean_scores)

# for more details about ORF annotations see http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_dna.README
orf_coding=seqinr::read.fasta(paste0(ref.data.dir,'orf_coding.fasta.gz'))
e_orfs=read.delim(paste0(ref.data.dir, 'E_orfs_column.txt'), stringsAsFactors=F, header=F, sep='\t')[,1]

sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3

unique.chrs=paste0('chr', as.roman(1:16))

gr.list=list()
coding.list=list()
intergenic.list=list()

#extract incidence matrix for each variant in each strain
#need this for rare variant association mapping code
extractZ=T

# granges object to store coordinates and indices
for(chr in unique.chrs ) { 
    print(chr)
    #from James
    vcf.in=read.vcfR(paste0(vcf.prefix, chr, '.vcf'))
#    vcf.in=read.vcfR(paste0(vcf.prefix, chr, '.vcf.gz'))
    
    # vcfR function to calculate minor allele frequency (why is this f(x) so slow, let's just diy this)
    #mafdf=maf(vcf.in, element=2)

    # a bunch of annoying code to deal with multiallelic sites
    alt.var.list=list()
    gr.obj.list=list()
    multi.allele.index=list()
     #up to 10 alleles, will handle fewer ggracefully
    ref=getREF(vcf.in)
    for(i in 1:10) {
        alt.var.list[[as.character(i)]]=as.character(
                 sapply(sapply(getALT(vcf.in), function(x) strsplit(x, ',')) ,
                        function(x) x[i]))
        
        #length of alternate variant allele minus 1
        alt.var.lenm1=as.vector(sapply(alt.var.list[[i]], function(x) nchar(x)))-1 
        if(sum(is.na( alt.var.lenm1))==length(alt.var.lenm1)) { break;} 
        #hack to construct GRanges object     
        alt.var.lenm1.hack=alt.var.lenm1
        alt.var.lenm1.hack[is.na( alt.var.lenm1)]=0 
        # keep track of sites that don't contain a variant for a given allele 
        multi.allele.index[[as.character(i)]]=!is.na(alt.var.lenm1)
            
        gr.obj.list[[as.character(i)]] = GRanges(
                   seqnames=chr, 
                   ranges=IRanges(start=getPOS(vcf.in), 
                                  end=  getPOS(vcf.in) + alt.var.lenm1.hack),
                       strand=c(rep('+', length(alt.var.list[[i]]))) #PAMp)),rep('-', length(PAMm)) ) 
                    )

       alt.var.list[[as.character(i)]]=  alt.var.list[[as.character(i)]][!is.na(alt.var.lenm1)]
       gr.obj.list[[as.character(i)]]= gr.obj.list[[as.character(i)]][!is.na(alt.var.lenm1)]
       gr.obj.list[[as.character(i)]]$alt.var=alt.var.list[[as.character(i)]] 
       gr.obj.list[[as.character(i)]]$idx=paste0( as.character(seqnames(gr.obj.list[[as.character(i)]])), 
                                                 '_', 
                                                 start(gr.obj.list[[as.character(i)]]),
                                                 ':', 
                                                 ref[!is.na(alt.var.lenm1)], '/',
                                                 alt.var.list[[as.character(i)]] )
    }

    #GRanges object merging
    gr.obj=unlist(as(gr.obj.list, "GRangesList"))

    # gr.obj$maf=mafdf[,4]
    # gr.obj$idx=paste0( as.character(seqnames(gr.obj)), '_', start(gr.obj))
    genome(gr.obj)='sacCer3'

    
    ogo=order(gr.obj) 
    gr.obj$alt.var[gr.obj$alt.var=="*"]=""

    gr.obj=gr.obj[ogo,]
    #map variants to genes and predicted coding consequences
    pC=predictCoding(gr.obj,
              txdb, sacCer3,
              DNAStringSet(gr.obj$alt.var), 
              ignore.strand=TRUE)
    pC=addProvean(pC, meltedProveanScores)


    # For simplicity, here we are going to redefine dubious as (dubious or pseudogene). There are 784 dubious and 12 psuedogenes.
    # get_pseudogenes(working.dir, orf_coding) #"YAR061W" "YCL074W" "YCL075W" "YFL056C" "YIL170W" "YIL171W" "YIL174W" "YIL175W" "YLL016W" "YLL017W" "YPL275W" "YPL276W"
    pC$dubious=!(pC$GENEID %in% names(orf_coding))

    # Essential ORFs as annotated here (http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt)
    pC$essential=pC$GENEID %in% e_orfs


    lVP=locateVariants(gr.obj, txdb, IntergenicVariants(), strand=F)
    lVP$idx=paste0(as.character(seqnames(lVP)), ':' , start(lVP))

    gr.obj$Coding=gr.obj$idx %in% pC$idx
    gr.obj$Intergenic=gr.obj$idx %in% lVP$idx
    
    gr.list[[chr]]=gr.obj
    coding.list[[chr]]=pC
    intergenic.list[[chr]]=lVP


    if(extractZ) {
    # extracts genotypes from vcf as allelic dosage
        Z.all= t( extract.gt(vcf.in, as.numeric=F))

        # parse genotypes, see comment above for structure of Z
        Z.subsets=list()
        for(i in names(multi.allele.index)){
            Ztemp=Z.all[,multi.allele.index[[i]]]
            Ztemp[Ztemp=='0/0']='0'
            Ztemp[Ztemp==paste0(i,'/',i)]='2'
            cm1=paste0(i)
            #cm2=paste0('*/',i)
            Ztemp[grepl(cm1, Ztemp)]='1'
            Ztemp=matrix(as.numeric(Ztemp), ncol=ncol(Ztemp), dimnames=dimnames(Ztemp))
            #Ztemp[is.na(Ztemp)]=0
            colnames(Ztemp)=gr.obj.list[[as.character(i)]]$idx
            Z.subsets[[as.character(i)]]=Ztemp
        }

        Z=do.call('cbind', Z.subsets)
        Z=Z[,ogo]
 

        acnt=colSums(Z,na.rm=T) #/(nrow(Zsub)*2)
        non_missing=apply(Z,2,function(x) sum(!is.na(x)))
        missing=nrow(Z)-non_missing
        af=acnt/(2*(non_missing))
        flipped=af>.5
        maf=af
        maf[maf>.5]=(1-maf[maf>.5])

        gr.obj$acnt=acnt
        gr.obj$missing=missing
        gr.obj$af=af
        gr.obj$maf=maf

        # all.equal(colnames(Z), gr.obj$idx)
        saveRDS(Z, file=paste0(out.dir, 'Z_', chr, '.RDS'))
    }

    saveRDS(gr.obj, file=paste0(out.dir, 'gr_', chr, '.RDS'))
    saveRDS(pC, file=paste0(out.dir , 'pC_', chr, '.RDS'))
}



cL=do.call('rbind', lapply(coding.list, DataFrame)) 
gL=do.call('rbind', lapply(gr.list, DataFrame))
iL=do.call('rbind', lapply(intergenic.list, DataFrame))

g=cL$GENEID
r=cL$REFAA.delta
p=cL$POSAA.delta
v=cL$VARAA.delta
grlv=paste0(g,':',r, ':' ,p, ':',v)

#can use AAeffect or idx for lookup
cL$AAeffect=grlv
cL$gEdit=cL$idx
cL$gEdit=gsub('_', ':',cL$gEditInfo)
cL$gEdit=gsub('\\/', ':',cL$gEditInfo)

saveRDS(cL, file = paste0(out.dir, 'vcf_codingEffects.RDS'))
saveRDS(gL, file = paste0(out.dir, 'vcf_granges.RDS')) #'/data/yeast/yeast_oligos/1002genomes/1002granges.RDS')
saveRDS(iL, file = paste0(out.dir, 'vcf_intergenic.RDS')) #'/data/yeast/yeast_oligos/1002genomes/1002intergenic.RDS')



# some example lookups of variants in pre-computed gRNA and -------------------------------------------------------------
# download the RDS files from VariantEngineering/onlyPAMs/tables/tables_gdrive.txt

# coding.effects=readRDS('/data/yeast/yeast_oligos/onlyPams/tables/codingEffects.RDS')
# guides=readRDS('/data/yeast/yeast_oligos/onlyPams/tables/guides.RDS')

#load pre-computed tables
coding.effects=readRDS(paste0(base.dir,'VariantEngineering/onlyPAMs/tables/codingEffects.RDS'))
guides=readRDS(paste0(base.dir,'VariantEngineering/onlyPAMs/tables/guides.RDS'))


# lookup by specific pam breaking variants
matchsub=cL[cL$gEdit %in% coding.effects$gEdit,]
nvc=coding.effects[match(matchsub$gEdit, coding.effects$gEdit),]
all.equal(matchsub$gEdit, gvc$gEdit)

gvc=guides[match(nvc$guideIndex, guides$guideIndex),]

#make sure you don't bring in repairTemplate from the guides data structure here 
gsubset=cbind(matchsub,gvc[,c('guide','guideIndex', 'guidePAM12PMcount', 'guideGCcontent', 'U6terminator')],nvc)

filter.vec=gsubset$guidePAM12PMcount==1  & gsubset$U6terminator==F
  
gsubset=gsubset[filter.vec,]

#extract guide sequence
gsubset$guide.seq=substr(gsubset$guide,1,20)
#extract repair sequence 
gsubset$repairTemplate
