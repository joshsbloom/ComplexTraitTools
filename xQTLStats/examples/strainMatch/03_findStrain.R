library(xQTLStats)
library(qs)
library(stringr)
library(ggpubr)
library(dplyr)
library(readr)
library(data.table)

#some global variables ----------------------------------------
# locations of things -----------------------------------------------------------------------
#location of vcf file from 1011 isolates, preprocessed as a qs object
#vcf.qs.file='/u/project/kruglyak/jsbloom/PUBLIC_SHARED/yeast/1011/mega_filtered.qs'
vcf.qs.file='/media/hoffman2/PUBLIC_SHARED/yeast/1011/mega_filtered.qs'

#location of GATK ASEReadCount tables, output from 02_align.sh
count_dir='/data2/strainMatch/count_files/' #xQTL/FluconazoleV1_041825/count_files/'

#location of key file, grab key from github
key.file='/data2/strainMatch/SampleSheet.csv'
#---------------------------------------------------------------------------------------------




mega_filtered=qs::qread(vcf.qs.file) #'/data1/yeast/reference/pop_vcfs/1011/mega_filtered.qs')
vcf=mega_filtered$vcf
# pull out the genotype entries (num vector w/ marker + genotype info added as charac vectors)
gt=mega_filtered$gt
# remove the og file from memory, retain vcf object and gt object only
rm(mega_filtered)

#read in expected samples, remove first 417 lines
samples=readr::read_csv(key.file, skip=417)

#get expected biallelic variant sites and store as a data.frame
gt_id=data.frame(id=paste0(data.table::tstrsplit(rownames(gt), '_')[[1]], '_', 
                           data.table::tstrsplit(rownames(gt), '_')[[2]]))

#build list to contain informations about reads at expected variant sites 
countdfs=list()
for(n in 1:nrow(samples) ) { #nrow(samples)) {
    print(n) 
    s=samples$'Sample_ID'[n]
    sname = paste0(s, '.txt') ##- "Sample012.txt"
    scounts <- readr::read_tsv(stringr::str_c(count_dir, "/",  sname))
    #if file isn't empty
    if(nrow(scounts)>0) {
    #making sure sites we see from ASEreadcount match expected in vcf
    scounts.sub = scounts[paste0(scounts$contig, "_", scounts$position) %in%       gt_id, ]
    #pull out only the stuff we want and overwrite scounts
    scounts = data.frame(id = paste0(scounts.sub$contig, "_", 
        scounts.sub$position), ref = scounts.sub$refCount, alt = scounts.sub$altCount)
    scounts = dplyr::left_join(gt_id, scounts, by = "id")
    countdfs[[s]]=scounts
    }

}
    
#calculate coverage per sample at the expected segregating sites 
tdepth.v=sapply(countdfs, function(x) { sum(x$ref+x$alt, na.rm=T)/sum(!is.na(x$ref)) })

#hard call variant sites (if alt read counts  > ref ref
ref.vec=sapply(countdfs, function(x) ifelse(x$alt>x$ref, 1,0) ) #sum(x$ref+x$alt)/nrow(x))

#crap.sites=apply(gt, 1, function(x) sum(is.na(x)))
hdist.matrix=matrix(NA,ncol(gt),ncol(ref.vec)) #ncol(gt))
rownames(hdist.matrix)=colnames(gt)
colnames(hdist.matrix)=colnames(ref.vec)

#calculate hamming distance between observed variants and expected for each strain
for(i in 1:ncol(ref.vec) ){
    print(i)
    hdist.matrix[,i]=colSums(ref.vec[,i]!=gt, na.rm=T)
    #hist(hdist, breaks=100)
    #which.min(hdist)
}


