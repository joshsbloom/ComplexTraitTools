library(Biostrings)
library(VariantAnnotation)
library(GenomicFeatures)
#library(ensemblVEP)
#library('org.Sc.sgd.db')
library(S4Vectors)
library(seqinr)
#library(rbamtools)
library(ShortRead)
library(rBLAST)
library(stringdist)

#wd='/data/yeast/yeast_oligos/ADE2/230502_M01373_0326_000000000-KVMK9/out/'#/home/jbloom/Desktop/CC_Hiseq_053018/fastq/'
wd='/data0/yeast/estopsV2/bcls6/out/'
setwd(wd)
fq.files=list.files(wd, recursive=F)
#tweak this put output somewhere else 
samples=unique(gsub('_R\\d_001.fastq.gz', '', fq.files[grep('.gz', fq.files)]))[1:8]
for(sm in samples) {
    print(sm)
  #  F.in=fq.files[grep(t, fq.files)]
   # F.in=F.in[grep('fq.gz', F.in)]
    R1.in=paste0(sm,'_R1_001.fastq.gz') #F.in[grep('1.fq', F.in)] #paste0(t, '_R1.fastq.gz')
    R2.in=paste0(sm,'_R2_001.fastq.gz') #F.in[grep('2.fq', F.in)] #paste0(t, '_R2.fastq.gz')
    R1.out=sm
    # replace PEAR with FLASH2
    system(paste('/home/jbloom/Local/FLASH2/flash2 -m 20 -t 36 --compress-prog=pigz --suffix=gz -o', R1.out, R1.in, R2.in))
    # PEAR version is considerably slower
    # also caution when specifying memory parameter ... bug here that caps the total number of reads 
    #system(paste('pear -v 70 -j 48 -n 190 -m 230 -f', R1.in, '-r', R2.in, '-o', R1.out))
    #R1.out2=paste0(t, '.assembled.fastq')
    #system(paste('pigz', R1.out2))
}

for(sm in samples[1:8]) {
    print(sm)
    system(paste0("umi_tools extract --extract-method=regex --bc-pattern='^(?P<umi_1>.{6}).+(?P<cell_1>.{30})(?P<discard_1>AGTACAGAACG){s<=1}(?P<discard_2>.+)' -I ", sm, '.extendedFrags.fastq.gz -S ', sm, '_umi.fastq.gz'))
}       
#for(sm in samples[4:6]) {
#    print(sm)
#    system(paste0("umi_tools extract --extract-method=regex --bc-pattern='^(?P<umi_1>.{6})' -I ", sm, '.extendedFrags.fastq.gz -S ', sm, '_umi.fastq.gz'))
#}  

#first extract umi and barcode
#umi_tools extract --extract-method=regex --bc-pattern='^(?P<umi_1>.{6}).+(?P<cell_1>.{30})(?P<discard_1>AGTACAGAACG){s<=1}(?P<discard_2>.+)' -I A_A02_S1_L001_R1_001.extendedFrags.fastq > test.fastq

#could add this back as filter 
#cutadapt --minimum-length 100 --output out.fq  /data/yeast/yeast_oligos/ADE2/flash/test.fastq




#table of guides and repair templates 
#y=readRDS('/data/yeast/yeast_oligos/ADE2/design/first_half_of_ade1_2_stops.RDS')
y=readRDS('/data0/yeast/estopsV2/estops500.RDS')

uguides=unique(y$guide)
uguides=DNAStringSet(uguides)
names(uguides)=as.character(seq(1:length(uguides)))
#writeXStringSet(uguides,'/data/yeast/yeast_oligos/ADE2/ref_amplicons/ade_oligos_expected_guides.fasta') # mode='a', full=FALSE, compress=F) #TRUE)

utemps=unique(y$repairTemplate)
utemps=DNAStringSet(utemps)
names(utemps)=as.character(seq(1:length(utemps)))
#writeXStringSet(utemps,'/data/yeast/yeast_oligos/ADE2/ref_amplicons/ade_oligos_expected_repairs.fasta') # mode='a', full=FALSE, compress=F) #TRUE)

#some expectations 
#m=as.matrix(stringdistmatrix(eoligos))
#m[m>5]=NA
#image.plot(m)

#m=as.matrix(stringdistmatrix(utemps))
#m[m>5]=NA
#image.plot(m)
ecombos=cbind(match(y$guide, uguides), match(y$repairTemplate, utemps))
ecombos.str=paste0(ecombos[,1], ':', ecombos[,2])



#define expected amplicon structure 
umi.length=6

#sequence upstream of guide
#short version

#additional primer handle (only used in long version)
a0='CTTCTCC'

#expected oligo
#array handle left 
#a1='GCAGTGAAAGATAAATGATC' #ATACGACTCACTATAGGGCGAATTGGGTACC'
a1='GCAGTGAAAGATAGGTGACC'
nchar.a1=nchar(a1)

#long version
gRNA.length=20

struct_region=toupper('gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtggtgc')
struct_region.length=nchar(struct_region)

struct_region_terminator='TTTTTTTGTTTTTTATGTCT'
srt=paste0(struct_region, struct_region_terminator) 
struct_region_terminator.length=nchar(srt)

repair_template.length=90

#sequence upstream of CAN1 barcode 
a2='TTCCCGACGAGAGTAAATGGCGAGGATACGTTCTCTATGG'
nchar.a2=nchar(a2)

barcode.length=30

#handle downstream of barcode
a3='AGTACAGAACGCTGAAGTGAAGAGAGAGCTCAAGCAAAGAGGTACCCAATTCGCC'
nchar.a3=nchar(a3)
#additional primer handle (only used in long version)
a4='CTATAGTG'
nchar.a4=nchar(a4)

eoligos=DNAStringSet(paste0(a1, y$guide, srt, y$repairTemplate, a2))
names(eoligos)=as.character(seq(1:length(eoligos)))

#writeXStringSet(eoligos,'/data/yeast/yeast_oligos/ADE2/ref_amplicons/ade_oligos_expected.fasta')

#system(paste0('makeblastdb -in /data/yeast/yeast_oligos/ADE2/ref_amplicons/ade_oligos_expected.fasta -dbtype nucl'))
#system(paste0('makeblastdb -in /data/yeast/yeast_oligos/ADE2/ref_amplicons/ade_oligos_expected_guides.fasta -dbtype nucl'))
#system(paste0('makeblastdb -in /data/yeast/yeast_oligos/ADE2/ref_amplicons/ade_oligos_expected_repairs.fasta -dbtype nucl'))

in.file1='/data/yeast/yeast_oligos/ADE2/230502_M01373_0326_000000000-KVMK9/out/D_A01_S1_L001_umi.fastq.gz' #/data/yeast/yeast_oligos/ADE2/flash/out.fq'

in.file1='/data/yeast/yeast_oligos/ADE2/230502_M01373_0326_000000000-KVMK9/out/D_B01_S2_L001_umi.fastq.gz' #/data/yeast/yeast_oligos/ADE2/flash/out.fq'

in.file1='/data/yeast/yeast_oligos/ADE2/230502_M01373_0326_000000000-KVMK9/out/D_C01_S3_L001_umi.fastq.gz' 

in.file1='/data/yeast/yeast_oligos/ADE2/230502_M01373_0326_000000000-KVMK9/out/D_D01_S4_L001_umi.fastq.gz' #/data/yeast/yeast_oligos/ADE2/flash/out.fq'

in.file1='/data/yeast/yeast_oligos/ADE2/230502_M01373_0326_000000000-KVMK9/out/D_E01_S5_L001_umi.fastq.gz' #/data/yeast/yeast_oligos/ADE2/flash/out.fq'

in.file1='/data/yeast/yeast_oligos/ADE2/230502_M01373_0326_000000000-KVMK9/out/D_F01_S6_L001_umi.fastq.gz' #/data/yeast/yeast_oligos/ADE2/flash/out.fq'


#read in trimmed data
nbuffer=1e8
blocksize=1e9

eamplicon.size=umi.length+width(eoligos)[1]
#+barcode.length+nchar(a3)

# could be more aggressive here
min.width=eamplicon.size*.8
error.rate=.10


for( s in samples ) {
in.file1=paste0(wd, s,'_umi.fastq.gz') #)'/data/yeast/yeast_oligos/ADE2/230502_M01373_0326_000000000-KVMK9/out/D_A01_S1_L001_umi.fastq.gz' #/data/yeast/yeast_oligos/ADE2/flash/out.fq'

fi1 =FastqStreamer(in.file1, nbuffer, readerBlockSize=blocksize,verbose=T)
#fi2 =FastqStreamer(in.file2, nbuffer, readerBlockSize=blocksize,verbose=T)

  #  repeat {
        #ShortReadQ class 
        rfq1=yield(fi1) 
        #rfq2=yield(fi2) 

        #246867
        
        #badreads 
        if(length(rfq1)==0 ) {break }
        br=width(rfq1)<min.width #) | (width(rfq2)< 20)
        if(sum(br)>0) { rfq1=rfq1[-which(br)]  }
        #rfq2=rfq2[-which(br)]
        #min width 220
        #231150
        
        # read the whole thing into memory ... this can obviously go very badly but should provide a 
        #template for a more memory efficient version     #rfq=readFastq(rr.file, withIDs=TRUE)
        # ---------------------------------------------------------------------------------------------
        cread=sread(rfq1)
        v1.start=which.isMatchingEndingAt(a1, cread, ending.at=(nchar.a1-2):(nchar.a1+2), max.mismatch=1, follow.index=T, with.indels=F)+1
      
        r.loc=nchar.a1+gRNA.length+struct_region_terminator.length

        repair.start=which.isMatchingEndingAt(struct_region_terminator, cread, ending.at=(r.loc-5:r.loc+5), max.mismatch=2, follow.index=T, with.indels=F)+1

        r.loc2=nchar.a1+gRNA.length+struct_region_terminator.length+repair_template.length
        
        #was a2
        #fairly significant loss at this step, debug
        repair.end=which.isMatchingStartingAt('TTCCCGACGAGAGTAAATGG', cread,
              starting.at=(r.loc2-5):(r.loc2+5), max.mismatch=2, follow.index=T, with.indels=F)-1

        
        bad.reads=(is.na(v1.start)|is.na(repair.start)|is.na(repair.end))
        rfq1=rfq1[!bad.reads]
        v1.start=v1.start[!bad.reads]
        repair.start=repair.start[!bad.reads]
        repair.end=repair.end[!bad.reads]


        cread=sread(rfq1)

        br=width(cread)
        rid=id(rfq1)
        ridc=as.character(rid)
        ridc=gsub('\\s.*', '', ridc)
        ridcs=data.table::tstrsplit(ridc, '_')
       
        #if extracting barcode and umi
        if(length(ridcs)==3) {
            umi=ridcs[[3]]
            barcode=ridcs[[2]]
        } 
        #if no barcode
        if(length(ridcs)==2) {
            umi=ridcs[[2]]
        }
        
        
        v1=subseq(cread, start=v1.start, width=gRNA.length) #5)
        repTemp=subseq(cread, start=repair.start, end=repair.end)
           
        rmatch=amatch(repTemp,utemps,method='lv', maxDist=5)

      #  system.time({ rmatch=amatch(repTemp,utemps,method='lv', maxDist=5)})
    
        gmatch=  amatch(v1, uguides, method='lv', maxDist=2)
    #b1=blast(db="/data/yeast/yeast_oligos/ADE2/ref_amplicons/ade_oligos_expected_repairs.fasta")
    
   # b1=blast(db="/data/yeast/yeast_oligos/ADE2/ref_amplicons/ade_oligos_expected.fasta")
   # c1=predict(b1, cread, BLAST_args="-max_target_seqs 1 -num_threads 32")
#    c1$QueryID=as.numeric(gsub('Query_', '', c1$QueryID))
    #c1$QueryID=as.numeric(c1$QueryID)
    #c1$GuideIndexMatch=gmatch[c1$QueryID]

    c1=data.frame(QueryID=1:length(gmatch), SubjectID=rmatch, GuideIndexMatch=gmatch)

    c1$GuideExpectedMatch=as.vector(as.character(uguides)[gmatch[c1$QueryID]]) #uguides[gmatch[c1$QueryID]]
    c1$repTempExpectedMatch=as.vector(as.character(utemps)[c1$SubjectID])


    close(fi1)

    c2=DataFrame(c1)
    c2$UMI=umi[c1$QueryID]
    c2$gRNA=v1[c1$QueryID]
    if(length(ridcs)==3) {
        c2$barcode=barcode[c1$QueryID]
    }
    if(length(ridcs)==2) {
        c2$barcode=NA
    }
 #   c2$repTemp=repTemp[c1$QueryID]

    ocombos=cbind(c1$GuideIndexMatch, c1$SubjectID)
    ocombos.str=paste0(ocombos[,1], ':', ocombos[,2])
    c2$ocombos=ocombos.str
    c2$ecombos=ecombos.str[match(ocombos.str, ecombos.str)]

 
 print(s)
 print(sum(nrow(c2)))
 print(sum(is.na(c2$ecombos))/nrow(c2)) #sum(!is.na(c2$GuideIndexMatch)) ) #sum(!is.na(c2$ecombos))/nrow(c2))
}






c2$trimmed=cread[c1$QueryID]
    csplit=split(c2, c2$SubjectID)

 sum(csplit[[46]]$ocombos==na.omit(unique(csplit[[46]]$ecombos)))/nrow(csplit[[46]])
 sum(csplit[[58]]$ocombos==na.omit(unique(csplit[[58]]$ecombos)))/nrow(csplit[[58]])


    algn=list()
    for(o in names(csplit)[45:length(csplit)]){
        print(o)
       # o='46'
    #46 and 58
    ex=csplit[[o]]
    o='46'
#     which.max(apply(m, 1, function(x) sum(x>2))/224)
    #stringdist(uguides[19], uguides, method='lv')
    stringdist(ex$repTempExpectedMatch[1], utemps, method='lv')













#g=readRDS('/data/yeast/yeast_oligos/onlyPams/tables/guides.RDS')

#x=readRDS('/data/yeast/yeast_oligos/onlyPams/tables/codingEffects.RDS')
#y=x[x$GENEID=='YAR015W' | x$GENEID=='YOR128C',]
#y=y[is.na(y$proveanEffects) & y$CONSEQUENCE=='nonsynonymous',]
#y=y[!duplicated(paste(y$guideIndex,y$repairTemplate)),]

#manually get first half of each gene   
#y=y[(y$GENEID=='YAR015W' &  (as.vector(unlist(y$POSAA.delta))<152)) |
#    (y$GENEID=='YOR128C'&   (as.vector(unlist(y$POSAA.delta))<271)),    ]

#sync up with the guide sequences 

#
#h=g[match(y$guideIndex, g$guideIndex),]
#y$guide=as.character(subseq(DNAStringSet(h$guide), 1,20))
#define some anchors
#PBS_gRNA1=toupper('gcagtgaaagataaatgatc')
#rightside_can1_1=toupper('TTCCCGACGAGAGTAAATGGCGAGGATACGTTCTCTATGG')
#rightside_can1_2='TTCCCGACGAGAGTAAATGG'


in.file1='/data/yeast/yeast_oligos/ADE2/flash/A_A02_S1_L001_R1_001.extendedFrags.fastq'
out.file='/data/yeast/yeast_oligos/ADE2/flash/A_A02_S1_L001_R1_001.extendedFrags_tbs.fastq.gz'


in.file1='/data/yeast/yeast_oligos/ADE2/flash/A_B02_S2_L001_R1_001.extendedFrags.fastq'

#barcodetrimFastqs=function(in.file1, out.file, umi.length=6, barcode.length1 = 20, barcode.length2 = 12, nbuffer=1e6, blocksize=1e9) 
    #nbuffer=5e6
    #in.file1='/media/jbloom/d1/rr_variant_oligo_results/Annotation_Run/hiseq_annotation_run/fastq_hiseq/S2_R1_trimmed.fastq.gz'
    #in.file2='/media/jbloom/d1/rr_variant_oligo_results/Annotation_Run/hiseq_annotation_run/fastq_hiseq/S2_R2_trimmed.fastq.gz'
    #rfout= '/media/jbloom/d1/rr_variant_oligo_results/Annotation_Run/hiseq_annotation_run/fastq_hiseq/S2_trimmed_barcode_stripped.fastq.gz'
      
eamplicon.size=umi.length+width(eoligos)[1]+barcode.length+nchar(a3)

# could be more aggressive here
min.width=eamplicon.size*.8
error.rate=.10

#define handle mismatch tolerances 
m.mm1=round(nchar.a1*error.rate)
m.mm2=round(nchar.a2*error.rate)

nbuffer=1e8
blocksize=1e9

    fi1 =FastqStreamer(in.file1, nbuffer, readerBlockSize=blocksize,verbose=T)
    #fi2 =FastqStreamer(in.file2, nbuffer, readerBlockSize=blocksize,verbose=T)

    repeat {
        #ShortReadQ class 
        rfq1=yield(fi1) 
        #rfq2=yield(fi2) 

        #badreads 
        if(length(rfq1)==0 ) {break }
        br=width(rfq1)<min.width #) | (width(rfq2)< 20)
        if(sum(br)>0) { rfq1=rfq1[-which(br)]  }
        #rfq2=rfq2[-which(br)]
        
        # read the whole thing into memory ... this can obviously go very badly but should provide a 
        #template for a more memory efficient version     #rfq=readFastq(rr.file, withIDs=TRUE)
        # ---------------------------------------------------------------------------------------------
        
        cread=sread(rfq1)
        br=width(cread)



        # allow for +/-1 given expected stagger length
        #findStagger=which.isMatchingStartingAt(a1, cread, starting.at=9:20, max.mismatch=m.mm, follow.index=T, with.indels=F)
     
        #v1 is gRNA
        v1.start=which.isMatchingEndingAt(a1, cread, ending.at=(umi.length+nchar.a1-3):(umi.length+nchar.a1+3), max.mismatch=m.mm1, follow.index=T, with.indels=F)+1
      
        r.loc=umi.length+nchar.a1+gRNA.length+struct_region_terminator.length

        #should replace this section with matchLRpatterns
        repair.start=which.isMatchingEndingAt(struct_region_terminator, cread, ending.at=(r.loc-10:r.loc+10), max.mismatch=3, follow.index=T, with.indels=F)+1

        r.loc2=umi.length+nchar.a1+gRNA.length+struct_region_terminator.length+repair_template.length
        
        repair.end=which.isMatchingStartingAt(a2, cread,
              starting.at=(r.loc2-10):(r.loc2+10),
                                              max.mismatch=4, follow.index=T, with.indels=F)-1

        #should replace this section with matchLRpatterns and use expected right handle
        v2.loc=umi.length+nchar.a1+gRNA.length+struct_region_terminator.length+repair_template.length+nchar.a2
        
        v2.start=which.isMatchingEndingAt(a2, cread,
              ending.at=(v2.loc-25):(v2.loc+25),
                                              max.mismatch=4, follow.index=T, with.indels=F)+1
       
        toolong=v2.start+barcode.length>br
        
        bad.reads=(is.na(v1.start) | is.na(v2.start) | is.na(repair.start) | is.na(repair.end) | toolong)

        cread=cread[!bad.reads]
        #UMI
        umi=subseq(cread, start=1, width=umi.length)
       
        #here barcode 1 is gRNA
        v1=subseq(cread, start=v1.start[!bad.reads], width=gRNA.length) #5)
        #barcode2 is plasmid barcode
        v2=subseq(cread, start=v2.start[!bad.reads], width=barcode.length) #5)

        repTemp=subseq(cread,start=repair.start[!bad.reads],end=repair.end[!bad.reads])

        trimmed.read=subseq(cread, start=umi.length+1,   end=v2.start[!bad.reads]-1)  
        qtrimmed=subseq(quality(quality(rfq1[!bad.reads])), start=umi.length+1,   end=v2.start[!bad.reads]-1)

    #b1=blast(db="/data/yeast/yeast_oligos/ADE2/ade_oligos_expected_guides.fasta")
    #c2=predict(b1,barcodes1[1:1000], BLAST_args='-task blastn-short')

    #image.plot(as.matrix(stringdistmatrix(uguides, method='lv')))

    gmatch=  amatch(v1, uguides, method='lv', maxDist=2)
    b1=blast(db="/data/yeast/yeast_oligos/ADE2/ref_amplicons/ade_oligos_expected_repairs.fasta")
    c1=predict(b1,repTemp, BLAST_args="-max_target_seqs 1 -num_threads 32")
    c1$QueryID=as.numeric(gsub('Query_', '', c1$QueryID))
    c1$QueryID=as.numeric(c1$QueryID)
    c1$GuideIndexMatch=gmatch[c1$QueryID]
    c1$GuideExpectedMatch=as.vector(as.character(uguides)[gmatch[c1$QueryID]]) #uguides[gmatch[c1$QueryID]]
    c1$repTempExpectedMatch=as.vector(as.character(utemps)[c1$SubjectID])
    
    c2=DataFrame(c1)
    c2$UMI=umi[c1$QueryID]
    c2$gRNA=v1[c1$QueryID]
    c2$barcode=v2[c1$QueryID]
    c2$repTemp=repTemp[c1$QueryID]

    ocombos=cbind(c1$GuideIndexMatch, c1$SubjectID)
    ocombos.str=paste0(ocombos[,1], ':', ocombos[,2])
    c2$ocombos=ocombos.str
    c2$ecombos=ecombos.str[match(ocombos.str, ecombos.str)]

    c2$trimmed=trimmed.read[c1$QueryID]
    csplit=split(c2, c2$SubjectID)

    algn=list()
    for(o in names(csplit)[45:length(csplit)]){
        print(o)
       # o='46'
    #46 and 58
    ex=csplit[[o]]
    o='46'
    #stringdist(uguides[19], uguides, method='lv')
    stringdist(ex$repTempExpectedMatch[1], utemps, method='lv')



    #annotate variants 
    subject.seq=eoligos[as.numeric(o)]

    #    mat = nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = FALSE, type='DNA')
     mat = nucleotideSubstitutionMatrix(match = 5, mismatch = -4, baseOnly = FALSE, type='DNA')

    threads=36

    observed.seqs=as.character(ex$trimmed)
    fobserved.seqs=as.factor(observed.seqs)
    to.align.seqs=as.character(levels(fobserved.seqs))
    
  # water -gapopen 5 -gapextend 2 subject.fasta patterns.fasta -outfile out.aln

    writeXStringSet(subject.seq,'~/Desktop/subject.fasta')
    tas=DNAStringSet(to.align.seqs)
    names(tas)=as.character(seq(1:length(tas)))
    writeXStringSet(tas,'~/Desktop/patterns.fasta')


    mpas=lapply(to.align.seqs, function(x) {
        pairwiseAlignment(pattern=x , #to.align.seqs[1], 
    #                  subject=subject.seq,type='local', substitutionMatrix=mat, gapOpening=25, gapExtension=0) })

                      subject=subject.seq,type='local', substitutionMatrix=mat, gapOpening=5, gapExtension=2) })

    ampas.mat=do.call('rbind', mclapply(mpas, function(x) as.matrix(aligned(x)), mc.cores=threads))
   
    tot.seqs=as.vector(table(fobserved.seqs)) #observed.seqs))
    tot.seqs.vec=rep(1:length(tot.seqs),tot.seqs)
    tot.seqs.split=split(tot.seqs, test)


       mismatches=mcmapply(function(x,y) { 
                                m=mismatchTable(x) 
                                m$instances=y[m$PatternId]
                                return(m)
                     }
                     ,x=mpas, y=tot.seqs, mc.cores=threads, SIMPLIFY=FALSE)

        mismatches=rbindlist(mismatches, idcol='Id')
        mismatches$PatternId=mismatches$Id
    
#        mismatches=do.call('rbind', mismatches)
        mismatches=mismatches[order(mismatches$instances, decreasing=T),]
        mismatches=mismatches[,-'Id']

       mismatch.table=sapply( split(mismatches, mismatches$SubjectStart), function(x) table(rep(x$PatternSubstring, x$instances)) )
        #crude.plot=sapply( split(mismatches, mismatches$SubjectStart), function(x) sum(x$instances)) #/length(tot.seqs.vec)
        mismatch.subject.pos=as.numeric(names(split(mismatches, mismatches$SubjectStart)))
        names(mismatch.table)=s2c(as.character(subject.seq))[ mismatch.subject.pos]

        attr(mismatch.table, 'chr.pos')= mismatch.subject.pos  #+ samples[idx, 'start']
        attr(mismatch.table, 'total.reads')=length( tot.seqs.vec)

        subject.start.list=mclapply(mpas, function(x) {start(subject(x))}, mc.cores=threads)

        deletionsL=list()
        for(i in 1:length(tot.seqs)){
              pas=mpas[[i]]         
              subject.start=subject.start.list[[i]]
              tot.seq=tot.seqs[i]

             widths=width(deletion(pas)[[1]])
             deletion.intervals=data.frame((deletion(pas))[[1]])
             deletion.intervals$subject.start=rep(subject.start,nrow(deletion.intervals)) #as.numeric(deletion.intervals$space)]
             deletion.intervals$start=deletion.intervals$start+deletion.intervals$subject.start-1
             deletion.intervals$end= deletion.intervals$end+deletion.intervals$subject.start-1
             deletions=data.frame(start=deletion.intervals$start, end=deletion.intervals$end, width=widths)
             if(nrow(deletions)==0) {next;}

             deletions$instances=tot.seq #s[deletions$PatternID]
             deletions=deletions[order(deletions$instances, decreasing=T),]
             deletionsL[[as.character(i)]]=deletions
        }
        deletions=deletionsL
        rm(deletionsL)

      deletions=rbindlist(deletions, idcol='PatternId')
      deletions=deletions[order(deletions$instances, decreasing=T),]

      insertionsL=list()
      for(i in 1:length(tot.seqs)){
          pas=mpas[[i]]         
          subject.start=subject.start.list[[i]]
          tot.seq=tot.seqs[i]
          widths=width(insertion(pas)[[1]])
          insertion.intervals=data.frame((insertion(pas))[[1]])
          insertion.intervals$subject.start=rep(subject.start,nrow(insertion.intervals)) #as.numeric(insertion.intervals$space)]
          insertion.intervals$start=insertion.intervals$start+insertion.intervals$subject.start-1
          insertion.intervals$end= insertion.intervals$end+insertion.intervals$subject.start-1
          insertions=data.frame(start=insertion.intervals$start, end=insertion.intervals$end, width=widths)
          if(nrow(insertions)==0) {next;}
          insertions$instances=tot.seq #s[insertions$PatternID]
          insertions=insertions[order(insertions$instances, decreasing=T),]
          insertionsL[[as.character(i)]]=insertions
      }
      insertions=insertionsL
      rm(insertionsL)
      insertions=rbindlist(insertions, idcol='PatternId')
      insertions=insertions[order(insertions$instances, decreasing=T),]
        

      dexpand=rep(IRanges(start=deletions$start, end=deletions$end),deletions$instances)
      iexpand=rep(IRanges(start=insertions$start, end=insertions$end),insertions$instances)

      dCount=countOverlaps(1:nchar(subject.seq), dexpand)
      iCount=countOverlaps(1:nchar(subject.seq), iexpand)

      
      #quantify
      deletions[deletions$start>45,]
      indels=rbind(deletions,insertions)
      indelsf=indels[indels$start>45,]
      print('indels')
      print(sum(indelsf[!duplicated(indelsf$PatternId),]$instances)/sum(tot.seqs))

      mismatchesf=mismatches[mismatches$SubjectStart>45,]
      print('mismatches')

      print(sum(mismatchesf[!duplicated(mismatchesf$PatternId),]$instances)/sum(tot.seqs))
     
      mrf=mismatches[,c(1,5,6,8,8)]
      names(mrf)=names(indels)

      jnt=rbind(deletions,insertions,mrf)
      jntf=jnt[jnt$start>45,]
      print('errors')
      print(sum(jntf[!duplicated(jntf$PatternId),]$instances)/sum(tot.seqs))


      algn[[o]]=list(
      subject.seq=subject.seq,   
      mismatches=mismatches,
      mismatch.table=mismatch.table,
      deletions=deletions,
      insertions=insertions,
      dCount=dCount,
      iCount=iCount,
      tot.mm=sum(mismatchesf[!duplicated(mismatchesf$PatternId),]$instances),
      tot.in=sum(indelsf[!duplicated(indelsf$PatternId),]$instances),
      tot.err=sum(jntf[!duplicated(jntf$PatternId),]$instances),
      tot.seq=sum(tot.seqs)
      )
    }
sum(sapply(algn, function(x) x$tot.err))/ts
sum(sapply(algn, function(x) x$tot.mm))/ts
sum(sapply(algn, function(x) x$tot.in))/ts

(sapply(algn, function(x) x$dCount))

mjoint=sapply(algn[c(46,58)], function(x) {
      h=rep(0,269)
      y=x$mismatch.table
      h[attr(y, 'chr.pos')]=sapply(y,sum)
      return(h) }
)
par(mfrow=c(3,1))
#plot(1:269,rep(0,269), type='n')
plot(rowSums(mjoint), main='mismatches', ylim=c(0,40))
text(1:269, rep(40,269), s2c(as.character(subject.seq)))

plot(rowSums((sapply(algn[c(46,58)], function(x) x$dCount))),main='deletions', ylim=c(0,40))
plot(rowSums((sapply(algn[c(46,58)], function(x) x$iCount))), main='insertions', ylim=c(0,40))

te=sum(sapply(algn[c(46,58)], function(x) x$tot.err))
ts=sum(sapply(algn[c(46,58)], function(x) x$tot.seq))

      plot(mismatch.subject.pos, sapply(mismatch.table, sum))
      plot(1:nchar(subject.seq), countOverlaps(1:nchar(subject.seq), dexpand))
      plot(1:nchar(subject.seq), countOverlaps(1:nchar(subject.seq), iexpand))

te=sum(sapply(algn, function(x) x$tot.err))
ts=sum(sapply(algn, function(x) x$tot.seq))

    
     





      dsub=deletions[deletions$start>45,]







    expected.combos=match(ocombos.str, ecombos.str)


         
        
    cname=id(rfq1[!bad.reads])
   # note how this is constructed useful for stripping barcode back out
        cname=subseq(cname, end=width(cname)-28)
        cname=xscat(cname, ':', BStringSet(umi),':', BStringSet(v1), ':', BStringSet(v2))
        writeFastq( ShortReadQ(sread=trimmed.read, quality=qtrimmed, id= cname), file=out.file, mode='a', full=FALSE, compress=TRUE)

        writeFastq( ShortReadQ(sread=trimmed.read, quality=qtrimmed, id= cname), file=out.file, mode='a', full=FALSE, compress=TRUE)

    }
    close(fi1)















#expected fasta
input.fasta=readDNAStringSet('/data/yeast/yeast_oligos/ADE2/ade_pilot.fasta', 'fasta')
output.fasta=subseq(input.fasta, start=umi.length+1,end=b2.loc)
writeXStringSet(output.fasta,'/data/yeast/yeast_oligos/ADE2/ade_oligos_expected.fasta') # mode='a', full=FALSE, compress=F) #TRUE)


#        cread2=sread(rfq2)
#
#        barcodes1=subseq(cread, start=1, width=barcode.length1) #5)
#        trimmed.read1=subseq(cread, start=barcode.length1+1) 
#        
#        barcodes2=subseq(cread2, start=width(cread2)-(barcode.length2-1), width=barcode.length2)
#        qtrimmed=narrow(quality(rfq1), start=barcode.length1+1 )#, end=width(quality(rfq))-20)
#        cname=as.character(id(rfq1))
#        cname=gsub('\\s', ':', cname)
#        cname=paste(cname, barcodes1, barcodes2, sep=':')
#        newname=BStringSet(cname)
#        writeFastq( ShortReadQ(sread=trimmed.read1, quality=qtrimmed, id= newname), file=out.file, mode='a', full=FALSE, compress=TRUE)
#    }
#close(fi1)
#close(fi2)
#}
#ref.fasta='/data/CRISPR_variant_engineering/rr_variant_oligos/results/coding_variants_LOD30.fasta'
ref.fasta='/data/yeast/yeast_oligos/ADE2/ade_oligos_expected.fasta' #/media/jbloom/d1/rr_variant_oligo_results/reference/round2_expected_oligos_LOD30.fasta'
ref.fasta='/data/yeast/yeast_oligos/ADE2/ade_oligos_expected_repairs.fasta' 
for(t in traits) {
 stripped.fastq='/data/yeast/yeast_oligos/ADE2/flash/A_A02_S1_L001_R1_001.extendedFrags_tbs.fastq.gz'

     
   #  paste0(wd, t, '_trimmed_barcode_stripped.fastq.gz')
 bam.out= '/data/yeast/yeast_oligos/ADE2/flash/A_A02_S1_L001_R1_001.extendedFrags_tbs.bam' # paste0(wd, t, '_matched.bam')
 system( paste("bwa mem -t 48", ref.fasta, stripped.fastq,  "| samtools sort -@ 12 -O BAM -m 8G >", bam.out))
 system( paste("minimap2 -t 48 -a", ref.fasta, stripped.fastq,  "| samtools sort -@ 12 -O BAM -m 8G >", bam.out))

  # or perhaps this 
 #system( paste('bwa mem -t 64',  ref.fasta, out.file1, '| samtools sort -O BAM -m 4G >', out.bam))
  system(paste("samtools index", bam.out))
}

#samtools view -F0x900 A_A02_S1_L001_R1_001.extendedFrags_tbs.bam > tes
#samtools view A_A02_S1_L001_R1_001.extendedFrags_tbs.bam > test.sam
library(tidyverse)
sam=read_delim('/data/yeast/yeast_oligos/ADE2/flash/test_assign.sam', col_names=F, delim='\t')



nbuffer=1e9
fi1 =FastqStreamer(stripped.fastq,nbuffer, readerBlockSize=blocksize,verbose=T)
    #fi2 =FastqStreamer(in.file2, nbuffer, readerBlockSize=blocksize,verbose=T)
        rfq1=yield(fi1) 
close(fi1)






    test= cut(1:length(to.align.seqs), threads)
    tas=split(to.align.seqs,test)
    mpas=mclapply(tas, function(x) {
    pairwiseAlignment(pattern=x, 
                      subject=subject.seq,type='local', substitutionMatrix=mat, gapOpening=5, gapExtension=2) }, mc.cores=threads)
    
    ampas.mat=do.call('rbind', mclapply(mpas, function(x) as.matrix(aligned(x)), mc.cores=threads))
    tot.seqs=as.vector(table(fobserved.seqs)) #observed.seqs))
    tot.seqs.vec=rep(1:length(tot.seqs),tot.seqs)
    tot.seqs.split=split(tot.seqs, test)




    mismatches=mcmapply(function(x,y) { 
                                m=mismatchTable(x) 
                                m$instances=y[m$PatternId]
                                return(m)
                     }
                     ,x=mpas, y=tot.seqs.split, mc.cores=threads, SIMPLIFY=FALSE)

    
    
       mismatch.table=sapply( split(mismatches, mismatches$SubjectStart), function(x) table(rep(x$PatternSubstring, x$instances)) )
        #crude.plot=sapply( split(mismatches, mismatches$SubjectStart), function(x) sum(x$instances)) #/length(tot.seqs.vec)
        mismatch.subject.pos=as.numeric(names(split(mismatches, mismatches$SubjectStart)))
        names(mismatch.table)=s2c(as.character(subject.seq))[ mismatch.subject.pos]

       # attr(mismatch.table, 'chr.pos')= mismatch.subject.pos + samples[idx, 'start']
        attr(mismatch.table, 'total.reads')=length( tot.seqs.vec)

        subject.start.list=mclapply(mpas, function(x) {start(subject(x))}, mc.cores=threads)
        
        deletions=mcmapply(function(pas, subject.start, tot.seqs) {
                            deletion.intervals=RangedData(deletion(pas))
                            deletion.intervals$subject.start=subject.start[as.numeric(deletion.intervals$space)]
                            deletion.intervals$start=start(deletion.intervals)+deletion.intervals$subject.start-1
                            deletion.intervals$end= end(deletion.intervals)+deletion.intervals$subject.start-1
                            deletions=data.frame(PatternID=as.numeric(deletion.intervals$space), start=deletion.intervals$start, end=deletion.intervals$end, width=width(deletion.intervals))
                            deletions$instances=tot.seq[deletions$PatternID]
                            deletions=deletions[order(deletions$instances, decreasing=T),]
                            return(deletions) }, pas=mpas, subject.start=subject.start.list, tot.seq=as.list(tot.seqs),mc.cores=threads, SIMPLIFY=FALSE)
        deletions=do.call('rbind', deletions)
        deletions=deletions[order(deletions$instances, decreasing=T),]

        insertions=mcmapply(function(pas, subject.start, tot.seqs) {
                            insertion.intervals=RangedData(insertion(pas))
                            insertion.intervals$subject.start=subject.start[as.numeric(insertion.intervals$space)]
                            insertion.intervals$start=start(insertion.intervals)+insertion.intervals$subject.start-1
                            insertion.intervals$end= end(insertion.intervals)+insertion.intervals$subject.start-1
                            insertions=data.frame(PatternID=as.numeric(insertion.intervals$space), start=insertion.intervals$start, end=insertion.intervals$end, width=width(insertion.intervals))
                            insertions$instances=tot.seqs[insertions$PatternID]
                            insertions=insertions[order(insertions$instances, decreasing=T),]
                            return(insertions) }, pas=mpas, subject.start=subject.start.list, tot.seq=as.list(tot.seqs),mc.cores=threads, SIMPLIFY=FALSE)
        insertions=do.call('rbind', insertions)
        insertions=insertions[order(insertions$instances, decreasing=T),]
    
    
    bdist=stringdistmatrix(ex$barcode,method='lv') 












    talign=pairwiseAlignment(pattern=to.align.seqs[213], subject=subject.seq, type='local', substitutionMatrix=mat, gapOpening=5, gapExtension=2)

   x= smith_waterman(to.align.seqs[213], as.character(subject.seq), match=5, mismatch=-4, gap=-5,lower=F)


































# pilot processing ... ignoring barcodes for now 
#Step 1 after read pairing
#                                           [BARCODE_1 20bp]   [oligo primer 15bp]  [REPAIR TEMPLATE 101 bp]                                                                     [mluI + cut site + structural seq 28bp]  [GUIDE 20nt]     [guide promoter (low quality) ] 
stagger_seq ATACGACTCACTATAGGGCGAATTGGGTACC TTAACGATAATCTTACTGAA GCCGTGTGAAGCTGG AGTTCCTCTAATTACGAACGAGCAAGCAAATTAGTATTGTGTGGGAGACGGAATGAAATTATTTAACAAGGAAGAAGCGAGCTTCGAGACTCTGTTAAAAC ACGCGTATCGATCGGCATCCTCTAAAAC GTCTCCCACACAATACTAAT GGTCACCTATCTTTCACTGCGGAGAAGTTTCGAACG
stagger_seq ATACGACTCACTATAGGGCGAATTGGGTACC ACTCGTCCCCCCAACTGTTG GCCGTGTGAAGCTGG TATGATTCCATTTCATGCCCAAAATATAATAAAATACATGTATTTCAGGCTGTCACCTTTAATCCATCACTGGCAGAACAGCAAATTTCAACTTTTGATGA ACGCGTATCGATCGGCATCCTCTAAAAC CCTGAAATACATGTATTTTA GGTCACCTATCTTTCACTGCGGAGAAGTTTCGAACG


#Step 2 (also what gets sequenced for selection experiments 

#R1                                         [BARCODE_1 -20bp]    [fixed sequence]  [REPAIR TEMPLATE PT1  ~70bp]                                          R2_rc [REPAIR TEMPLATE PT2 25bp] [mlu1] [BARCODE_2 12bp]  [additional plasmid sequence]
stagger_seq ATACGACTCACTATAGGGCGAATTGGGTACC TTTGATAGATAATGGGTTAA GCCGTGTGAAGCTGGC TATTTAACTCGGTTTGGGTTTTCCAACTGCACTAGTTCCAAGATGATCATCCGAAGGAACTAGCTGTAA ------- ATCGTTAATGGCTGTCCCACTTAGA ACGCGT GATCGTGTTAAG  GGATAGTTGCCCTCCGAGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGAAATTGTAAGCGTTAATAATTCAGAAGAACTCGTCAAGAAGG
stagger_seq ATACGACTCACTATAGGGCGAATTGGGTACC TCAAATACATTGACGGTATG GCCGTGTGAAGCTGGC CATATCATTTCTTGATTGCAAATTCGCTTGGTTCTACGGTCATAAAATTGACAATATGTGTTGATAGT  ------- ATATAGTAGATTTGAAGAAAGATTT ACGCGT GATCACCCTCGC  GGATAGTTGCCCTCCGAGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGAAATTGTAAGCGTTAATAATTCAGAAGAACTCGTCAAGAAGG    




#S2_R2 [tail of repair template]  [mlu1]  [barcode 2]
ATCGTTAATGGCTGTCCCACTTAGA ACGCGT GATCGTGTTAAG
ATATAGTAGATTTGAAGAAAGATTT ACGCGT GATCACCCTCGC
GGGGATAACCGCTTAGTGAAGAGTA ACGCGT CTTGGAATCTTA



pear -f yn2_S3_L001_R1_001.fastq.gz -r yn2_S3_L001_R2_001.fastq.gz  -v 80 -j 48 -n 190 -o test.fq


library(VariantAnnotation)
library(GenomicFeatures)
#library(ensemblVEP)
library('org.Sc.sgd.db')
library(S4Vectors)
library(seqinr)
library(rbamtools)
library(ShortRead)
library(Biostrings)


txdb=makeTxDbFromGFF(file=paste0('/data/CRISPR_variant_engineering/rr_variant_oligos/',
                                 'reference/saccharomyces_cerevisiae.gff'))
tscl=as.list(txdb)
tscm=merge(tscl$transcripts, tscl$genes, by='tx_id')
library(org.Sc.sgd.db)
xx <- as.list(org.Sc.sgdGENENAME) #org.Sc.sgdCOMMON2ORF)
flatlist=sapply(xx, function(x) x)
flatlist[is.na(flatlist)]=''
gene.GR=GRanges(seqnames=tscm$tx_chrom, ranges=IRanges(start=tscm$tx_start, end=tscm$tx_end), 
             strand=tscm$tx_strand, ORF=tscm$gene_id)
gene.GR$NAME=as.character(flatlist[gene.GR$ORF])
load('/data/CRISPR_variant_engineering/rr_variant_oligos/rr.vcf.vep.RData')

# what if ref is also 
csq=parseCSQToGRanges(RR.vcf)
ccoding=csq[!is.na(csq$CDS_position),]

library(tidyr)
library(data.table)
#traits=c('yw1', 'yw2', 'yn1', 'yn2', 'sdswt', 'sdsnit', 'mnwt', 'mnnit')
traits=c('PRE_1', 'PRE_2', 
         'zero_hr_1_1',
         'zero_hr_1_2',
         'zero_hr_1_3',
         'zero_hr_2_1',
         'zero_hr_2_2',
         'zero_hr_2_3',
         'YPD_1_1',
         'YPD_1_2',
         'YPD_1_3',
         'YPD_2_1',
         'YPD_2_2',
         'YPD_2_3',
         'Cd_1', 'Cd_2', 'Co_1', 'Co_2', 'Paraquat_1', 'Paraquat_2',  'SDS_1', 'SDS_2', 'YNB_1', 'YNB_2')

#'Lactate', '15C', 'Cobalt', 'SDS', 'MN')
# adapter trimming 

#wd='/media/jbloom/d1/rr_variant_oligo_results/CC_Hiseq_053018/fastq/'
# move to SSD for faster IO
# modify with custom 

wd='/data/rr_variant_oligos/CC_Hiseq081318/fastq/'#/home/jbloom/Desktop/CC_Hiseq_053018/fastq/'
setwd(wd)
fq.files=list.files(wd, recursive=T)
for(t in traits[-c(1)]) {
    print(t)
    F.in=fq.files[grep(t, fq.files)]
    F.in=F.in[grep('fq.gz', F.in)]
    R1.in=F.in[grep('1.fq', F.in)] #paste0(t, '_R1.fastq.gz')
    R2.in=F.in[grep('2.fq', F.in)] #paste0(t, '_R2.fastq.gz')
    R1.out=t
    # replace PEAR with FLASH2
    system(paste('/home/jbloom/Local/FLASH2/flash2 -m 70 -M 100 -t 60 --compress-prog=pigz --suffix=gz -o', R1.out, R1.in, R2.in))
    # PEAR version is considerably slower
    # also caution when specifying memory parameter ... bug here that caps the total number of reads 
    #system(paste('pear -v 70 -j 48 -n 190 -m 230 -f', R1.in, '-r', R2.in, '-o', R1.out))
    #R1.out2=paste0(t, '.assembled.fastq')
    #system(paste('pigz', R1.out2))
}
   
#fq.files=list.files(wd)
#for(t in traits) {
#    F.in=fq.files[grep(t, fq.files)]
#    R1.in=F.in[grep('1.fq', F.in)] #paste0(t, '_R1.fastq.gz')
#    R2.in=F.in[grep('2.fq', F.in)] #paste0(t, '_R2.fastq.gz')
#
#    R1.out=paste0(t, '.fastq.gz')
#    R1t.out=paste0(t,'_R1_trimmed.fastq.gz')
#    system(paste('/usr/local/bin/cutadapt -j 64 -n 3 -g ATACGACTCACTATAGGGCGAATTGGGTACC -e 0.20  -o ', R1t.out, R1.in)) 
#    
#    R2t.out=paste0(t,'_R2_trimmed.fastq.gz')
#    RRout=paste0(t, '_R2_revc.fastq.gz')
#    system(paste("seqtk seq -r",R2.in, " | gzip > ", RRout))
#
#    system(paste("/usr/local/bin/cutadapt -j 64 -n 3",
#             "-a GGATAGTTGCCCTCCGAGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGAAATTGTAAGCGTTAATAATTCAGAAGAACTCGTCAAGAAGG",
#             "-e 0.20 -o ",
#             R2t.out, RRout))
#    system(paste('rm -rf', RRout))
#}


setwd(wd)
for(t in traits) {
     R1t.in=paste0(t, '.extendedFrags.fastq.gz')
     #paste0(t,'_R1_trimmed.fastq.gz')
     #R2t.in=paste0(t,'_R2_trimmed.fastq.gz')
     ot.out=paste0(t, '_trimmed_barcode_stripped.fastq.gz')
     barcodetrimFastqs(R1t.in, ot.out)
}

#in.file1=paste0(t, '.extendedFrags.fastq.gz')
#in.file1=paste0(t,'.assembled.fastq.gz')
# could ouput the 25 bp of rt in read 2 here too



















#bamtools filter -tag NM:i:0 -in yw1_matched.bam  -out test.bam
#-tag MD:Z:122 
#for(t in traits) {
#    bam.out=   paste0(wd, t, '_matched.bam')
#    txt.out=paste0(wd, t, '.txt')
#    system(paste("samtools idxstats", bam.out, ">", txt.out))
#}
for(t in traits)  {
    bam.in = paste0(wd, t, '_matched.bam')
    bam.out=paste0(wd, t, '.bam')
    txt.out=paste0(wd, t, '.txt')
    system(paste('bamtools filter -tag NM:i:0 -tag MD:122 -in', bam.in, '-out', bam.out))
    system(paste("samtools index", bam.out))
    system(paste("samtools idxstats", bam.out, ">", txt.out))
}

#---------------------------------------------------------------------------------------------------------
#bwa mem -t 35 /data/CRISPR_variant_engineering/rr_variant_oligos/results/coding_variants_LOD30.fasta \
#/media/jbloom/d1/rr_variant_oligo_results/Annotation_Run/hiseq_annotation_run/fastq_hiseq/S2_trimmed_barcode_stripped.fastq.gz \
#| samtools sort -@ 12 -O BAM -m 24G > S2_matched.bam
#samtools index  S2_matched.bam

# from annotation run
load('/media/jbloom/d1/rr_variant_oligo_results/Annotation_Run/hiseq_annotation_run/fastq_hiseq/bdf_ec.RData')
load('/media/jbloom/d1/rr_variant_oligo_results/Annotation_Run/hiseq_annotation_run/fastq_hiseq/bdf2_ec.RData')

expected.seqs=read.fasta(ref.fasta,as.string=T, forceDNAtolower=F)



#all.exp=list()

library(stringdist)
library(rbamtools)
library(tidyr)
library(dplyr)
barcode.counts=list()
#bam.file = paste0(wd, t, '_matched.bam')
#/media/jbloom/d1/rr_variant_oligo_results/Annotation_Run/hiseq_annotation_run/fastq_hiseq/S2_matched.bam'

bam.files=paste0(wd, traits, '_matched.bam')
readers=sapply(bam.files, bamReader, idx=TRUE)
#reader = bamReader(bam.file,idx=TRUE)
contigs = (refSeqDict(getHeaderText(readers[[1]])))@SN
contigs.length=(refSeqDict(getHeaderText(readers[[1]])))@LN
# construct a table to parse the information about each oligo construct
#ac=bamCountAll(reader)
#bamClose(reader)
#reader=bamReader(bam.file,idx=TRUE)
#pb =txtProgressBar(min = 1, max =length(contigs), style = 3)
   for(n in 710:length(contigs)) {
       print(paste(n, contigs[n]))
       #setTxtProgressBar(pb, n)
       # argh, 0 indexing
       # branges[[contigs[n]]]=
       brl=lapply(readers, bamRange, coords=c(n-1, 1, contigs.length[n]))  
       #br=bamRange(reader, coords=c(n-1, 1, contigs.length[n]))
       al=lapply(brl, function(x) DataFrame(data.frame(x)) ) 
       #12s
       al=lapply(al, function(x) {
           test= BStringSet(x$name)
           x$barcode1=subseq(test, end=width(test)-13, width=20)
           x$barcode2=subseq(test, end=width(test), width=12)
           x$megabarcode=DNAStringSet(xscat(x$barcode1, x$barcode2))
           return(x)
        })
       #3s
       u1=lapply(al, function(x) as.character(x$megabarcode))
       #1s 
       u1=DNAStringSet(do.call(c, u1),use.names=F)
       #2.5s
       u1.order=order(u1) 
       #2.5s
       if(length(u1.order)==0) {            next;       }
       ubc=rle(as.vector(u1)[u1.order])
       #1s
       ubl=sort(ubc$lengths,decreasing=T)
       if(length(ubl)<2 ) {next;}
       mlength=ubl[max(which((cumsum(ubl)/sum(ubl))<.98))]
       if(mlength<10) {mlength=10}
       if(length(ubl)<50) {mlength=1}
       # tweak this 
       core.barcodes=ubc$values[ubc$lengths>mlength]
       core.barcodes.counts=ubc$lengths[ubc$lengths>mlength]

       if(length(core.barcodes)<2) {next;}
       kd2=stringdistmatrix(core.barcodes, method='lv')
       
       hc=stats::hclust(kd2)
       b.groups=cutree(hc, h=5)
       ad=data.frame(barcode=core.barcodes, barcode.count=core.barcodes.counts, group=b.groups,stringsAsFactors=F)

       # true barcodes 
       core.barcodes=as.character(sapply(split(ad,ad$group), function(x) x$barcode[which.max(x$barcode.count)]))
       
       kd=stringdistmatrix(a=ubc$values, b=core.barcodes, method='lv') 
       #3.3s
       mind=apply(kd,1,min) 
       #1s
       if(length(mind)==0) {next;}

       which.mind=apply(kd,1,which.min)
       which.mind[mind>5]=NA

       #print(sum(ubc$lengths[mind>5])/sum(ubc$lengths))
       ad=data.frame(barcode=ubc[[2]], barcode.count=ubc[[1]], core.lookup=core.barcodes[which.mind],stringsAsFactors=F)
      
       al2=lapply(al, function(x) {
              x$barcode.group=ad$core.lookup[match(x$megabarcode, ad$barcode)]
              y=x[!is.na(x$barcode.group),]
              return(data.frame(y));  })
       names(al2)=traits
       df=rbindlist(al2, idcol='condition')
       df$condition=factor(df$condition, levels=traits)

       ss=split(df$seq, df$barcode.group)
       max.seq=sapply(ss, function(x){
                       if(length(x)==1) { (as.vector(x)) }
                       else {y=rle(x); return(as.vector(y$values[which.max(y$lengths)]))}}            )
       ss=split(df$megabarcode, df$barcode.group)
       max.barcode=sapply(ss, function(x){
                       if(length(x)==1) { (as.vector(x)) }
                       else {y=rle(x); return(as.vector(y$values[which.max(y$lengths)]))}}            )
       df$seq[match(names(max.seq), df$barcode.group)]=max.seq
       df$megabarcode[match(names(max.barcode), df$barcode.group)]=max.barcode
       tallydf=df %>% group_by(condition, barcode.group) %>% tally() %>% spread(., condition, n, 0, drop=F)
       print(tallydf)
       df2=data.frame(sequence=df$seq[match(tallydf$barcode.group, df$barcode.group)],
                      perfectMatch=df$seq[match(tallydf$barcode.group, df$barcode.group)]==expected.seqs[[contigs[n]]],
                      barcode=df$megabarcode[match(tallydf$barcode.group, df$barcode.group)],
                      tallydf,   stringsAsFactors=F)
       barcode.counts[[(contigs)[n]]]=df2
   }
#close(pb)
lapply(readers, bamClose)
        
#save(barcode.counts, file='/data/rr_variant_oligos/CC_Hiseq081318/barcodeCounts.RData')
load('/data/rr_variant_oligos/CC_Hiseq081318/barcodeCounts.RData')


# SKIP .... add more annotations to oligo list 
    # add provean score to oligo 
    load('/media/jbloom/d1/CRISPR_base_editor/reference/provean_scores.RData')

    #load array annotation information 
    load('/data/CRISPR_variant_engineering/rr_variant_oligos/results/SKP2_coding_variants_LOD30.RData')
    #(ompos3)
    # provean_scores
    ompos3$Provean=NA
    for(i in 1:nrow(ompos3)) {
        print(i)
        goi = ompos3$Gene[i] 
        ppos=as.numeric(ompos3$Protein_position[i])
        aa.change=tstrsplit(ompos3$Amino_acids[i], '/')
        ompos3$Provean[i]=tryCatch( { provean_scores[[goi]][ppos,aa.change[[2]]] }, error=function(e) {NA})
    }

    x1=sapply(ompos3$pam_seq, s2c)
    x2=sapply(ompos3$pam_alt_seq, s2c)
    ed=mapply(function(x,y) sum(x!=y), x=x1, y=x2)
    ompos3$pam_edit_distance=ed
    ompos3$syn_diff=grepl(';', ompos3$pam_alt_seq)
    save(ompos3, file='/data/rr_variant_oligos/CC_Hiseq081318/oligo_annotations.RData')
#------------------------------------------------------------------------------------------
    
load('/data/rr_variant_oligos/CC_Hiseq081318/oligo_annotations.RData')
  
load('/data/rrv2/genotyping/RData/QTGsorted.RData')
#QTGsortedYPD=QTGsorted[QTGsorted$trait.1=='YPD;;1',]

crosses.to.parents=list(
     '375'=c("M22", "BYa"),
     'A'  =c("BYa", "RMx"),
     '376'=c("RMx", "YPS163a"),
     'B'  =c("YPS163a", "YJM145x"),
     '377'=c("YJM145x", "CLIB413a"),
     '393'=c("CLIB413a", "YJM978x"),
     '381'=c("YJM978x", "YJM454a"),
    '3008'=c("YJM454a", "YPS1009x"),
    '2999'=c("YPS1009x", "I14a"),
    '3000'=c("I14a", "Y10x"),
    '3001'=c("Y10x", "PW5a"),
    '3049'=c("PW5a", "273614xa"),
    '3003'=c("273614xa", "YJM981x"),
    '3004'=c("YJM981x", "CBS2888a"),
    '3043'=c("CBS2888a", "CLIB219x"),
    '3028'=c("CLIB219x", "M22")
    )

load('/data/rrv2/genotyping/RData/parents.list.RData')
#source('/data/rrv2/genotyping/code/mapping_fx.R')
# add marker.name.n column
parents.list=lapply(parents.list, function(x) {
                  z=x;
                  z$marker.name.n=paste0(z$marker.name, '_', seq(1:nrow(z)))
                  return(z) })


per.variant.crosses=list()
per.variant.QTL.info=list()
for(test.marker in as.character(ompos3$var.name) ) {
    t.marker=as.character(test.marker)
    print(t.marker)
    # test.marker='chrXIV:467219_A/G' #chrI:203922_A/G'
    tmarker=gsub(':|/', '_', t.marker)
    per.variant.crosses[[t.marker]]=sapply(parents.list, function(x) sum(grepl(tmarker, x$marker.name)))
    # match against one of the two crosses
    otest0=as.logical(sapply(QTGsorted$pCausal, function(x) sum(grepl(tmarker, names(x)))))
    otest1=as.logical(sapply(QTGsorted$pCausal.1, function(x) sum(grepl(tmarker, names(x)))))
    otest=otest0 | otest1
    print(QTGsorted[otest,])
    per.variant.QTL.info[[t.marker]]=QTGsorted[otest,]
}
# 
sapply(per.variant.crosses, function(x) sum(x))
# private only variants 
which(sapply(per.variant.crosses, function(x) sum(x)==2)
#1159
#save(per.variant.crosses, file='/data/rr_variant_oligos/CC_Hiseq081318/perVariantCrosses.RData')
#save(per.variant.QTL.info, file='/data/rr_variant_oligos/CC_Hiseq081318/perVariantQTLinfo.RData')
sum(sapply(per.variant.QTL.info, nrow)==0)
# 1184 variants  ... from LOD>30 QTL 
# about 30% no longer overlap large effect QTL (QTL got tighter) or were not a high confidence variant

causalOligos=list()
for(t in unique(QTGsorted$trait.1)) {
    print(t)
    yp=lapply(per.variant.QTL.info, function(x) x[grep(t, x$trait),])
    ypm=names(which(sapply(yp, nrow)>0))
    #[1] 130
    ypdname=tstrsplit(ypm, ':|_', type.convert=T)
    # 2500 isn't quite enough to collapse QTL
    nQTL=reduce(GRanges(seqnames=ypdname[[1]], IRanges(start=(ypdname[[2]]), end=(ypdname[[2]])))+2500)
    print(ypm)
    print(nQTL)
    causalOligos[[t]]=list(set=ypm, nQTL=nQTL)
}
#save(causalOligos, file='/data/rr_variant_oligos/CC_Hiseq081318/causalOligos.RData')

load('/data/rr_variant_oligos/CC_Hiseq081318/causalOligos.RData')

sapply(causalOligos, function(x) length(x$set))



# restructure 
library(data.table)
library(tidyr)
library(lme4)
library(MASS)
library(foreach)
library(doMC)
library(glmmTMB)
library(edgeR)
registerDoMC(cores=70)

bc=rbindlist(barcode.counts, idcol='oligo')
names(bc)[grep('YPD_1_3', names(bc))]='YPD_2_4'
names(bc)[grep('YPD_2_1', names(bc))]='YPD_1_4'
names(bc)[grep('zero_hr_1_3', names(bc))]='zero_hr_2_4'
names(bc)[grep('zero_hr_2_1', names(bc))]='zero_hr_1_4'
names(bc)=gsub('o_hr', 'ohr', names(bc))

# analysis of all barcodes regardless of perfect match 
bca=data.frame(bc)
bcpl=gather(bca, key=experiment, value=count, PRE_1:YNB_2)
bcpl=separate(bcpl, experiment, c('condition', 'transformation'), sep='_', extra="drop", remove=F)
bcpl=separate(bcpl, oligo, c('chr', 'position', 'variant'), sep=':', extra="drop", remove=F)
bcpl$condition=factor(bcpl$condition, levels=unique(bcpl$condition))
# per barcode
dodge=position_dodge(width = .4)
ggplot(bcpl, aes(x=condition, y=log2(count+1), color=variant))+facet_wrap(~perfectMatch)+
    geom_violin(position=dodge, width=1)+
    geom_boxplot(width=0.05,position=dodge, alpha=.5)+
   scale_colour_manual(name = 'variant',values=c('red','blue'),  labels = c('synonymous','synonymous+variant'))

#only oligos with perfect match 


# filter on barcodes that decrease between t=0 and pre-gal
# see filter_PREvst0.R
#bcp2=bcp
#bcp2$good.barcode=good.barcodes[match(bcp2$barcode, names(good.barcodes))]
#bcp2=bcp2[!is.na(bcp2$good.barcode),]
#bcps=split(bcp2, bcp2$oligo)

bcp=bca[bca$perfectMatch,]
bcps=split(bcp, bcp$oligo)
nbcps=names(bcps)
onames=tstrsplit(nbcps, '_|:', type.convert=T)
oncore=unique(paste0(onames[[1]], ':', onames[[2]], '_', onames[[3]]))
oncount=sapply(oncore, function(x) length(grep(x, nbcps)))
oncore=names(oncount)[oncount==2]

bcl=gather(bcp, key=experiment, value=count, PRE_1:YNB_2)
bcl=separate(bcl, experiment, c('condition', 'transformation'), sep='_', extra="drop", remove=F)
bcl=separate(bcl, oligo, c('chr', 'position', 'variant'), sep=':', extra="drop", remove=F)
bcl$condition=factor(bcl$condition, levels=unique(bcl$condition))

# per barcode
ggplot(bcl, aes(x=condition, y=log2(count+1), color=variant))+
    geom_violin(position=position_dodge(width = .4), width=1)+
    geom_boxplot(width=0.05,position=position_dodge(width = 0.4), alpha=.5)+
   scale_colour_manual(name = 'variant',values=c('red','blue'),  labels = c('synonymous','synonymous+variant'))
ggsave(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/syn_dropout.png', width=11, height=5)

# per oligo
ocl.counts=data.frame(t(sapply(bcps, function(x) colSums(x[,c(6:ncol(x))]) )))
ocl=ocl.counts
ocl$oligo=rownames(ocl)
ocl=gather(ocl, key=experiment, value=count, PRE_1:YNB_2)
ocl=separate(ocl, experiment, c('condition', 'transformation'), sep='_', extra="drop", remove=F)
ocl=separate(ocl, oligo, c('chr', 'position', 'variant'), sep=':', extra="drop", remove=F)
ocl$condition=factor(ocl$condition, levels=unique(ocl$condition))
#per oligo
ggplot(ocl, aes(x=condition, y=log2(count+1), color=variant))+
    geom_violin(position=position_dodge(width = .5), width=1)+
    geom_boxplot(width=0.05, position=position_dodge(width = 0.5), alpha=.5)+
    scale_colour_manual(name = 'variant',values=c('red','blue'),  labels = c('synonymous','synonymous+variant'))

#ggplot(bcl, aes(x=condition, y=log2(count+1), color=variant))+facet_wrap(~experiment)+geom_boxplot()
#ggplot(bcl, aes(x=log2(count+1)))+facet_grid(~experiment+transformation)+geom_histogram(binwidth=.5)
#geom_histogram(binwidth=.5)

# 
#glmer.nb
#d.only=bcp2[,c(6:(ncol(bcp2)-1))]
# normalization factors based on SA only
d.only=bcp[,c(6:(ncol(bcp)))]
ncountsAll1=colSums(d.only)
ncountsAll2=apply(d.only, 2, function(x) median(x[x>quantile(x,.8) & x<quantile(x,.99)]) )
ncountsAll3=apply(d.only,2, function(x) median(x[x>30]))
#ncounts4=calcNormFactors(bcp[,c(6:ncol(bcp))])
ncountsAll5=calcNormFactors(d.only, refColumn='YPD_1_1')
names(ncountsAll5)=names(ncountsAll3)

d.only=bcp[grep('SA$', bcp$oligo),c(6:(ncol(bcp)))]
d.only.s=split(d.only, bcp$oligo[grep('SA$', bcp$oligo)])
ncountsSA3=apply(sapply(d.only.s, function(y) sapply(y, function(z) sum(z))),1 , median)
#ncountsSA3=apply(sapply(d.only.s, function(y) sapply(y, function(z) sum(z))),1 , quantile, .75)

d.only=bcp[grep('S$', bcp$oligo),c(6:(ncol(bcp)))]
d.only.sa=split(d.only, bcp$oligo[grep('S$', bcp$oligo)])
ncountsS3=apply(sapply(d.only.sa, function(y) sapply(y, function(z) sum(z))),1 , median)
#ncountsS3=apply(sapply(d.only.sa, function(y) sapply(y, function(z) sum(z))),1 , quantile, .75)


z1=d.only$zerohr_1_1+d.only$zerohr_1_2+d.only$zerohr_1_4
y1=d.only$YPD_1_1+d.only$YPD_1_2+d.only$YPD_1_4
t1=z1+y1


z2=d.only$zerohr_2_4+d.only$zerohr_2_2+d.only$zerohr_2_3
y2=d.only$YPD_2_4+d.only$YPD_2_2+d.only$YPD_2_3
t2=z2+y2
keep=t1>500 & t2>500
cor(c(log2(y1+1)-log2(t1+1))[keep], c(log2(y2+1)-log2(t2+1))[keep], method='spearman')
plot(c(log2(y1+1)-log2(t1+1))[keep], c(log2(y2+1)-log2(t2+1))[keep])




#ncountsSA1=colSums(d.only)
#ncountsSA2=apply(d.only, 2, function(x) median(x[x>quantile(x,.8) & x<quantile(x,.99)]) )
#ncountsSA3=apply(d.only,2, function(x) median(x[x>30]))
#ncounts4=calcNormFactors(bcp[,c(6:ncol(bcp))])
#ncountsSA5=calcNormFactors(d.only, refColumn='YPD_1_1')
#names(ncountsSA5)=names(ncountsSA3)

#dnorm.counts=(log2(d.only+1)-log2(ncounts2))
#cc=cor(dnorm.counts)
#dnorm.counts=(log2(d.only+1)-(ncounts5))
#cc=cor(dnorm.counts)
#dnorm.counts=(log2(d.only+1)-log2(ncounts3))
#cc=cor(dnorm.counts)
#col1 <- colorRampPalette(c("white", "white","red"),bias=.2)
#corrplot(cc, addCoef.col=T, method='shade', col=col1(100), order='hclust')


doNB=function(mm) {
        if(length(unique(mm$variant))<2 | length(unique(mm$condition))<2 | length(unique(mm$barcode))<2 ) {
        return(list(mu=NA, disp=NA, S_Cont=NA, SA_Cont=NA, S_Cond=NA, SA_Cond=NA, l2fc=NA, p=NA)) } else {
        gtb=glmmTMB(count~-1+variant+condition+variant*condition+offset(log(off))+(1|barcode),data=mm, family=list(family="nbinom2", link="log")) # , ziformula=~1)
        gtbs= summary(gtb)$coefficients$cond
        #l2fc=((gtbs[4,1]-gtbs[2,1])-(gtbs[3,1]-gtbs[1,1]))/log(2)
        
        py=as.vector(predict(gtb, type='response'))
        #gtb=glmmTMB(count~-1+variant+condition+variant*condition+offset(log(off))+(1|barcode)+(1|experiment),data=mm, family=list(family="nbinom2", link="log")) # , ziformula=~1)
        #plot(log2(getME(gmer, 'mu')), log(py))
        #cor(log2(py+1),log2(mm$y+1))^2
        #overdispersion
        #variance =μ(1 +μ/k):  Hardin and Hilbe(2007))
        
        disp=sigma(gtb)
        mu=sum(log(py))
        ps=drop1(gtb, test='Chisq')[[4]][2]
           
        conds=paste(mm$condition, mm$variant, sep='_')
        gmeans=sapply(split(log2(py), conds), mean)
        SA_Cond=gmeans[paste0(tt, '_synonymous+variant')]
        S_Cond=gmeans[paste0(tt, '_synonymous')]
        SA_Cont=gmeans[paste0(vcont, '_synonymous+variant')]
        S_Cont=gmeans[paste0(vcont, '_synonymous')]
        l2fc=(as.numeric((SA_Cond-S_Cond)-(SA_Cont-S_Cont)))
        return(list(mu=mu, disp=disp, S_Cont=S_Cont, SA_Cont=SA_Cont, S_Cond=S_Cond, SA_Cond=SA_Cond, l2fc=l2fc, p=ps))
        }
}







doNB.SA=function(mm2) {
        if(length(unique(mm2$condition))<2 | length(unique(mm2$barcode))<2 ) {
            return(list(l2fc=NA, p=NA)) } else {
        gtb=glmmTMB(count~condition+offset(log(off))+(1|barcode),data=mm2, family=list(family="nbinom2", link="log")) # , ziformula=~1)
        ps=drop1(gtb, test='Chisq')[[4]][2]
        gtbs= summary(gtb)$coefficients$cond
        l2fc=(gtbs[2,1])/log(2)
        return(list(l2fc=l2fc, p=ps))
        }
}


# all contrasts vs PRE
gtraits=c('zerohr') #, 'Cd_', 'Co_', 'Paraquat_', 'SDS_', 'YNB_', 'YPD_')
vcont='PRE'
vsPre3=list()

# all contrasts vs zero_hr
gtraits=c('Cd', 'Co', 'Paraquat', 'SDS', 'YNB', 'YPD')
vcont='zerohr'
vsZero3=list()

# all contrasts vs YPD
gtraits=c('Cd', 'Co', 'Paraquat', 'SDS', 'YNB')
vcont='YPD'
vsYPD2=list()

print(vcont)
for( tt in gtraits ) {
    print(tt)
    #nvec=list(         joint=list(mu=NA, disp=NA,S_Cont=NA, SA_Cont=NA, S_Cond=NA, SA_Cond=NA, l2fc=NA, p=NA),
    #          indep=list('1'=list(mu=NA, disp=NA, S_Cont=NA, SA_Cont=NA, S_Cond=NA, SA_Cond=NA, l2fc=NA, p=NA),
    #                     '2'=list(mu=NA, disp=NA, S_Cont=NA, SA_Cont=NA, S_Cond=NA, SA_Cond=NA,  l2fc=NA, p=NA)))
    nvec=list(         joint=list(l2fc=NA, p=NA),
              indep=list('1'=list(l2fc=NA, p=NA),
                         '2'=list(l2fc=NA, p=NA)))

    o.out=foreach(o = oncore) %dopar% {
        #for(oo in 1:length(oncore)) {
        #o=oncore[oo]
        print(o)
        SA=paste0(o,':SA')
        S=paste0(o,':S')
        xdf=rbindlist(bcps[c(SA,S)])
        xdf$oligo=tstrsplit(xdf$oligo, ':')[[3]]
        names(xdf)[grep('oligo', names(xdf))]='variant'
        xdf=data.frame(xdf,stringsAsFactors=T)
        
        cond=grep(tt, names(xdf))
        cont=grep(vcont, names(xdf))
        pre1=grep('PRE_1', names(xdf))
        pre2=grep('PRE_2', names(xdf))

        if(vcont=='PRE') { 
        Y=xdf[,c(1:5,cont,cond)]
        } else { 
        Y=xdf[,c(1:5,pre1,pre2,cont,cond)]
        }

        cont1=grep(paste0(vcont, '_1'), names(Y))
        cont2=grep(paste0(vcont, '_2'), names(Y))
        
        cond1=grep(paste0(tt, '_1'), names(Y))
        cond2=grep(paste0(tt, '_2'), names(Y))
        #rcount=rowSums(Y[,c(cont1,cont2, cond1,cond2)])
        #Y[rcount<5, ]=NA
        # if(length(cont1)>1) {
            rcount1=rowSums(Y[,c(cont1,cond1)])
        #} else {
        #    rcount1=Y[,c(cont1,cond1)]
        #}
        #if(length(cont2)>1) {
            rcount2=rowSums(Y[,c(cont2,cond2)])
        #} else {
        #    rcount2=Y[,c(cont2,cond2)]
        #}
        Y[rcount1<10, c(cont1,cond1)]=NA
        Y[rcount2<10, c(cont2,cond2)]=NA

        if(nrow(Y)==0) {return(nvec); } #ps[oo]=1; next;}
        mm=gather(Y, experiment, count, 6:ncol(Y), convert=T)
        mm$variant=ifelse(mm$variant=='SA', 'synonymous+variant', 'synonymous')
        mm$variant=factor(mm$variant)
        mm$barcode=factor(mm$barcode)
        if(vcont=='PRE') { 
                 mm$condition=factor(tstrsplit(mm$experiment, '_')[[1]],levels=c(vcont, tt))
        } else { 
                 mm$condition=factor(tstrsplit(mm$experiment, '_')[[1]],levels=c('PRE',vcont, tt))
        }
        mm$transformation=factor(tstrsplit(mm$experiment, '_')[[2]])
        mm$experiment=as.factor(mm$experiment)
        mm$off=NA
        mm$off[mm$variant=='synonymous']=ncountsS3[mm$experiment[mm$variant=='synonymous']]
        mm$off[mm$variant=='synonymous+variant']=ncountsSA3[mm$experiment[mm$variant=='synonymous+variant']]
        
        #mm$off=ncountsAll3[mm$experiment]
        #mm$offSA=ncountsSA3[mm$experiment]

        mm=mm[!is.na(mm$count),]
        #add pre-gal 
        #mm=mm[mm$variant=='synonymous',] 
        #mm=mm[mm$transformation=='1',]

     #   ggplot(mm, aes(x=condition, y=log2((count+.5)/off) , group=interaction(barcode,variant,transformation)))+
     #       theme(legend.position = "none")+
     #       scale_y_continuous(name="log2(normalized counts)")+
     #       facet_wrap(~transformation+variant)+
     #       stat_summary(geom="line", fun.y="median", aes(color=barcode))+
     #       geom_jitter(size=.5)+
     #       geom_boxplot(aes(group=experiment),alpha=.1)+
     #       theme(legend.position = "none")+ggtitle(o)
       # ggsave(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/01_7_pca1_result_barcodes.png', width=10, height=10)

          
       #geom_violin(draw_quantiles=c(.25,.5,.75), alpha=.2)+
        if(vcont!='PRE') { 
                 mm=mm[mm$condition!='PRE',]
        }

        if(length(levels(mm$variant))<2 | length(levels(mm$condition))<2 | length(levels(mm$barcode))<2 ) { return(nvec) } else { 
        #joint=doNB(mm)
        mm2=mm[mm$variant=='synonymous',]
        jointSA=doNB.SA(mm2)
        #print(paste(o, joint$mu, round(joint$disp,2), round(joint$l2fc,2), joint$p, ifelse(joint$p<.0001, "***", "")))
        #print(paste(o, round(jointSA$l2fc,2), jointSA$p, ifelse(jointSA$p<.0001, "***", "")))
        #mm.s=split(mm, mm$transformation)    
        mm2.s=split(mm2, mm2$transformation)    
        #indep=lapply(mm.s, doNB)
        indepSA=lapply(mm2.s, doNB.SA)
        #return(list(joint=joint, indep=indep))
        return(list(joint=jointSA, indep=indepSA))
        }
    }
    names(o.out)=oncore
    results.joint=rbindlist(lapply(o.out, function(x) x$joint),idcol='variant')
    results.joint=data.frame(results.joint, fdr=p.adjust(results.joint$p, method='fdr'))
    
    results.t1=rbindlist(lapply(o.out, function(x) x$indep[[1]]),idcol='variant')
    results.t1=data.frame(results.t1, fdr=p.adjust(results.t1$p, method='fdr'))

    results.t2=rbindlist(lapply(o.out, function(x) x$indep[[2]]),idcol='variant')
    results.t2=data.frame(results.t2, fdr=p.adjust(results.t2$p, method='fdr'))
    
    results=list(joint=results.joint, t1=results.t1, t2=results.t2)
   
    if(vcont =='PRE') { vsPre3[[tt]]=results  } # list(ps=ps, l2fc=l2fc)}
    if(vcont =='zerohr'){ vsZero3[[tt]]=results }#list(ps=ps, l2fc=l2fc)}
    if(vcont =='YPD') { vsYPD3[[tt]]=results } #list(ps=ps, l2fc=l2fc)}
}


vsPRE_S=list()
vsPRE_S[[tt]]=results
save(vsPRE_S, file='/data/rr_variant_oligos/CC_Hiseq081318/vPre_S.RData')
#save(vsPre3, file='/data/rr_variant_oligos/CC_Hiseq081318/vPre3_NB.RData')


#library(metap)
#save(vsPre, file='/data/rr_variant_oligos/CC_Hiseq081318/vPre_NB.RData')
#save(vsZero, file='/data/rr_variant_oligos/CC_Hiseq081318/vZero_NB.RData')
load('/data/rr_variant_oligos/CC_Hiseq081318/vZero_NB.RData')



cypd=unique(unlist(sapply(causalOligos[grep('YPD;;', names(causalOligos))], function(x) x$set)))

egenes=read.delim('/data/CRISPR_nonsense/E_orfs_column.txt',header=F)
ompos3$varInPG=ifelse((ompos3$var_pam_D<=-1 | ompos3$var_pam_D==2), FALSE,TRUE)

testPre=merge(vsPre3[['zerohr']]$joint, ompos3, by.x='variant', by.y='var.name')
testPre=merge(vsPRE_S[['zerohr']]$joint, ompos3, by.x='variant', by.y='var.name')
testPre[which(testPre$l2fc<(-2000)),]$l2fc=NA
testPre$essential= testPre$Gene %in% as.character(egenes[[1]])
testPre$gl=round(width(gene.GR[match(testPre$Gene, gene.GR$ORF)])/3)
testPre$STRAND[testPre$STRAND=='-1']="-"
testPre$STRAND[testPre$STRAND=='1']="+"
with(testPre, t.test(c(l2fc)~varInPG))
with(testPre, t.test(c(l2fc)~STRAND))
with(testPre, t.test(c(l2fc)~pam_strand))
with(testPre, t.test(c(l2fc)~pam_strand==STRAND))
with(testPre, t.test(c(l2fc)~essential))
with(testPre, t.test(c(l2fc)~syn_diff))
with(testPre, t.test(c(l2fc)~ifelse(pam_edit_distance>1,T,F)))

multiguides=unique(rpg$values[rpg$lengths>1])
tPmg=testPre[testPre$guide %in% multiguides,]

tb=ifelse(testPre$l2fc>(-1), T, F)
fisher.test(with(testPre, sapply(split(tb, varInPG), table)))
fisher.test(with(testPre, sapply(split(tb, STRAND), table)))
fisher.test(with(testPre, sapply(split(tb, pam_strand), table)))
fisher.test(with(testPre, sapply(split(tb, pam_strand==STRAND), table)))
fisher.test(with(testPre, sapply(split(tb, essential), table)))
fisher.test(with(testPre, sapply(split(tb, syn_diff), table)))
fisher.test(with(testPre, sapply(split(tb, ifelse(pam_edit_distance>1,T,F)), table)))


table(rbind(tb, testPre$varInPG))

with(testPre, .test(c(l2fc)~varInPG))
with(testPre, t.test(c(l2fc)~STRAND))
with(testPre, t.test(c(l2fc)~pam_strand))
with(testPre, t.test(c(l2fc)~pam_strand==STRAND))
with(testPre, t.test(c(l2fc)~essential))
with(testPre, t.test(c(l2fc)~syn_diff))
with(testPre, t.test(c(l2fc)~ifelse(pam_edit_distance>1,T,F)))


with(testPre, t.test(c(S_Cond-S_Cont)~varInPG))
with(testPre, t.test(c(S_Cond-S_Cont)~STRAND))
with(testPre, t.test(c(S_Cond-S_Cont)~pam_strand))
with(testPre, t.test(c(S_Cond-S_Cont)~pam_strand==STRAND))
with(testPre, t.test(c(S_Cond-S_Cont)~essential))
with(testPre, t.test(c(S_Cond-S_Cont)~syn_diff))
with(testPre, t.test(c(S_Cond-S_Cont)~ifelse(pam_edit_distance>1,T,F)))
#data:  c(S_Cond - S_Cont) by ifelse(pam_edit_distance > 1, T, F)
#t = -2.9, df = 660, p-value = 0.004
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -0.5627 -0.1042
#sample estimates:
#mean in group FALSE  mean in group TRUE 
#             -2.716              -2.383 
with(testPre, t.test(c(S_Cond-S_Cont)~U6terminator))
#data:  c(S_Cond - S_Cont) by U6terminator
#t = -8.1, df = 230, p-value = 2e-14
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -1.4992 -0.9152
#sample estimates:
#mean in group FALSE  mean in group TRUE 
#             -2.797              -1.590 
with(testPre, cor.test(c(S_Cond-S_Cont),Provean))
with(testPre, stripchart((S_Cond-S_Cont)~U6terminator, method='jitter', vertical=T, ylim=c(-10,5)))



ggplot(data.frame(testPre),aes(x =start,y=-log10(p), color=(fdr<.05))) + geom_point(size=2) + scale_color_brewer(palette='Set2')+
      theme_bw() + facet_grid(~ chr,scales = "free_x", space='free_x')+xlab('variant position')+ ggtitle('variant oligo results for t=0 vs Pre')

ggplot(data.frame(testPre),aes(y=S_Cont,x=var_pam_D, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylim(c(1,10))+ylab('log2 (conditional mean for synonymous(PRE))')
ggplot(data.frame(testPre),aes(y=SA_Cont,x=var_pam_D, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylab('log2 (conditional mean for synonymous + nonysynonymous (PRE))')
ggplot(data.frame(testPre),aes(y=S_Cond,x=var_pam_D, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylim(c(1,10))+ylab('log2 (conditional mean for synonymous(t=0))')
ggplot(data.frame(testPre),aes(y=SA_Cond,x=var_pam_D, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylab('log2 (conditional mean for synonymous + nonysynonymous (t=0))')
ggplot(data.frame(testPre),aes(y=SA_Cond-SA_Cont,x=var_pam_D, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylab('log2 (conditional mean for synonymous + nonysynonymous (t=0))')+ylim(-8,5)
ggplot(data.frame(testPre),aes(y=S_Cond-S_Cont,x=var_pam_D, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylab('log2 (conditional mean for synonymous + nonysynonymous (t=0))')+ylim(-8,5)
ggplot(data.frame(testPre),aes(y=S_Cond-S_Cont,x=STRAND, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylab('log2 (conditional mean for synonymous + nonysynonymous (t=0))')+ylim(-8,5)

ggplot(data.frame(testPre),aes(y=l2fc,x=var_pam_D, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylab('log2 (conditional mean for synonymous + nonysynonymous (t=0))')

yy=testPre$SA_Cond-testPre$SA_Cont


ggsave(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/synonymous_control_dropout_by_position.png', width=11, height=8)





allr=vsPre3[['zerohr']]$joint
ompos3$varInPG=ifelse((ompos3$var_pam_D<=-1 | ompos3$var_pam_D==2), FALSE,TRUE)
test=merge(allr, ompos3, by.x='variant', by.y='var.name')
testS=test[which(test$fdr<.1),]
plot(-log10(test$p), col=ifelse(test$variant %in% cypd, 'red', 'black'))
# enrichment for called significant for oligos where variant is in pam or guide 
fisher.test(rbind(table(testS$varInPG), table(test$varInPG)))

# genomic trends 
egenes=read.delim('/data/CRISPR_nonsense/E_orfs_column.txt',header=F)
test$Essential= test$Gene %in% as.character(egenes[[1]])
test2=data.frame(test)

save.image(file='/data/rr_variant_oligos/CC_Hiseq081318/working_090618.RData')


# plot by position
#CD
testCD=merge(vsZero$Cd$joint, ompos3, by.x='variant', by.y='var.name')
ggplot(data.frame(testCD),aes(x =start,y=-log10(p), color=(fdr<.05))) + geom_point(size=2) + scale_color_brewer(palette='Set2')+
      theme_bw() + facet_grid(~ chr,scales = "free_x", space='free_x')+xlab('variant position')+ ggtitle('variant oligo results for Cd vs t=0')
ggsave(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/Cd_results.png', width=20, height=8)
#YPD
ggplot(test2,aes(x =start,y=-log10(p), color=(fdr<.05))) + geom_point(size=2) + scale_color_brewer(palette='Set2')+
      theme_bw() + facet_grid(~ chr,scales = "free_x", space='free_x')+xlab('variant position')+ ggtitle('variant oligo results for YPD vs t=0')
ggsave(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/02_YPD_results.png', width=20, height=8)

test2$causalVariants=test2$variant %in% cypd
ggplot(test2,aes(x =start,y=-log10(p), color=causalVariants)) + geom_point(size=2) + scale_color_brewer(palette='Set2')+
      theme_bw() + facet_grid(~ chr,scales = "free_x", space='free_x')+xlab('variant position')+ ggtitle('variant oligo results for YPD vs t=0')
ggsave(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/YPD_possibly_causal_results.png', width=20, height=8)



#non-synonymous

ggplot(test2,aes(y=SA_Cont,x=var_pam_D, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylab('log2 (conditional mean for synonymous + nonysynonymous (t=0))')
ggsave(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/variant_control_dropout_by_position.png', width=11, height=8)

# synonymous 
ggplot(test2,aes(y=S_Cont,x=var_pam_D, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylim(c(1,10))+ylab('log2 (conditional mean for synonymous(t=0))')
ggsave(file='/home/jbloom/Dropbox/Lab Meeting - Presentations/090518/synonymous_control_dropout_by_position.png', width=11, height=8)



# u6 terminator explains enrichment for S
ggplot(test2,aes(y=SA_Cont,x=var_pam_D, color=varInPG)) + geom_jitter( width=.2)+geom_boxplot(alpha=.1)+scale_color_brewer(palette='Set1')+
    xlab("variant position relative to PAM  + NG[G] -")+ylab('log2 (conditional mean for synonymous + nonysynonymous (t=0))')



ggplot(test2,aes(x =pam_start,y=SA_Cont)) + geom_point(size=1) + ylim(c(0,25)) + geom_smooth(span=1)+
      theme_bw() +facet_grid(~ chr,scales = "free_x", space='free_x') 

t3=test2[which(test2$fdr<.1),]


ggplot(test2,aes(x =pam_start,y=-log10(p), color=Essential)) + geom_point(size=1) + ylim(c(0,25)) +
      theme_bw() + facet_grid(~ chr,scales = "free_x", space='free_x') 
ggplot(test2,aes(x =pam_start,y=S_Cont, color=Essential)) + geom_point(size=1) + ylim(c(0,12)) +
      theme_bw() + facet_grid(~ chr,scales = "free_x", space='free_x') 
ggplot(test2,aes(x =pam_start,y=S_Cond, color=Essential)) + geom_point(size=1) + ylim(c(0,12)) +
      theme_bw() + facet_grid(~ chr,scales = "free_x", space='free_x') 
ggplot(test2,aes(x =pam_start,y=SA_Cond, color=Essential)) + geom_point(size=1) + ylim(c(0,12)) +
      theme_bw() + facet_grid(~ chr,scales = "free_x", space='free_x') 



+ ylim(c(0,12)) +
      theme_bw() + facet_grid(~ varInPG,scales = "free_x", space='free_x') 

ggplot(test2,aes(x =S_Cont, color=Essential)) + geom_histogram(position='identity', fill='white')


+ ylim(c(0,25)) +
      theme_bw() + facet_grid(~ chr,scales = "free_x", space='free_x') 

  + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(size = 5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
      geom_vline(aes(xintercept=pos_fineMapping, colour=QTN), data=sct, linetype=1, size=.6, alpha=.6 )+
      geom_point(aes(x=pos_fineMapping, y=0, colour=dup), size=2, shape=17, data=sct )+scale_color_manual(values=c('FALSE'="blue", 'Identical'='green', 'TRUE'='red'))+
      ggtitle(paste(trait, ' - ' ,nrow(sct),  'total QTL (', sum(sct$peakWidth==1), 'QTN) identified by She and Jarosz', sep=' '))
      


#t0
t.test(test$SA_Cont~test$varInPG)
stripchart(test$SA_Cont~test$varInPG, vertical=T, method='jitter')

t.test(test$S_Cont~test$varInPG)
stripchart(test$S_Cont~test$varInPG, vertical=T, method='jitter')

#YPD
t.test(test$SA_Cond~test$varInPG)
stripchart(test$SA_Cond~test$varInPG, vertical=T, method='jitter')

t.test(test$S_Cond~test$varInPG)
stripchart(test$S_Cond~test$varInPG, vertical=T, method='jitter')

stripchart(c(test$SA_Cond-test$S_Cond)~test$varInPG, vertical=T, method='jitter')
t.test(c(test$S_Cont-test$S_Cond)~test$varInPG)
t.test(c(test$SA_Cont-test$SA_Cond)~test$varInPG)

stripchart(scale(test$S_Cond)~test$var_pam_D, vertical=T, method='jitter')
X=data.frame(y=c(test$S_Cond),vpd=factor(test$var_pam_D))
XX=model.matrix(lm(y~vpd-1,data=X))
y2=X$y[!is.na(X$y)]
summary(lm(scale(y2)~XX-1) )

x11()
stripchart(scale(test$SA_Cond)~test$var_pam_D, vertical=T, method='jitter')
X=data.frame(y=c(test$SA_Cond),vpd=factor(test$var_pam_D))
XX=model.matrix(lm(y~vpd-1,data=X))
y2=X$y[!is.na(X$y)]
summary(lm(scale(y2)~XX-1) )


x11()
test2=test[test$var_pam_D==0,]
x11()
stripchart(scale(test2$SA_Cond)~test2$syn_diff, vertical=T, method='jitter')
X=data.frame(y=c(test$SA_Cond),vpd=factor(test$var_pam_D))
XX=model.matrix(lm(y~vpd-1,data=X))
y2=X$y[!is.na(X$y)]
summary(lm(scale(y2)~XX-1) )

# annotation run use it
# variants getting in ????
# are more things happening than desired (clue from synonymous depletion)
# uncontrolled sources of noise 
# 

sapply(split(c(test$SA_Cond-test$S_Cond),test$var_pam_D), median, na.rm=T)

X=data.frame(y=c(test$SA_Cond-test$S_Cond),vpd=factor(test$var_pam_D))
XX=model.matrix(lm(y~vpd-1,data=X))
y2=X$y[!is.na(X$y)]
lm.fit(XX, y2) 
thsd=TukeyHSD(aov(y~vpd-1, data=X))
do.call('rbind', strsplit(rownames(thsd$vpd), '-'))



t.test(test$SA_Cont~(test$pam_edit_distance==1)) 
t.test(test$S_Cont~(test$pam_edit_distance==1)) 
t.test(test$SA_Cond~(test$pam_edit_distance==1)) 
t.test(test$S_Cond~(test$pam_edit_distance==1)) 

#wilcox.test(test$SA_Cont~(test$syn_diff)) 
#wilcox.test(test$S_Cont~(test$syn_diff)) 
#wilcox.test(test$SA_Cond~(test$syn_diff)) 
#wilcox.test(test$S_Cond~(test$syn_diff)) 





fisher.test(rbind(table(testS$variant %in% cypd), table(test$variant %in% cypd)))


allr=results.joint
r=allr[which(allr$fdr<.05),]
wilcox.test(ompos3$pam_edit_distance[ompos3$var.name %in% r$variant], ompos3$pam_edit_distance)
fisher.test(rbind(table(ompos3$syn_diff[ompos3$var.name %in% r$variant]), table(ompos3$syn_diff)))

, ompos3$pam_edit_distance

# variants that overlap a QTL enriched in significant oligos  
fisher.test(rbind(table(r$variant %in% ypm),table( oncore %in% ypm)))
r.overlap=r[r$variant %in% ypm,]

yp[r.overlap$variant]

overlap.yp.table=do.call('rbind', lapply(yp[r.overlap$variant], function(x) x[1,]))
ompos3[ompos3$var.name %in% r.overlap$variant,]

#calculate pam edit distance 

# pre ...is variant in targeting sequence. ... are differential oligos  concentrated in essential genes etc 
# pre that don't deplete could be problematic 
for(n in names(vsZero2)){
   print(n)
    print(sum(vsZero2[[n]]$joint$fdr<.05, na.rm=T))
   print( vsZero2[[n]]$joint[which(vsZero2[[n]]$joint$fdr<.05),])
}
h=vsZero2[[n]]$joint[which(vsZero2[[n]]$joint$fdr<.05),]


for(n in names(vsZero)){
   print(n)
    print(sum(vsZero[[n]]$joint$fdr<.05, na.rm=T))

   print( vsZero[[n]]$joint[which(vsZero[[n]]$joint$fdr<.05),])
}
#if increase at t0 vs pre toss 

plot(vsZero[[1]]$joint$disp, vsZero[[6]]$joint$disp)


# non-disruptive to pam or guide are depleting more 



msig=match(r$variant, rownames(ompos3))
t.test(ompos3$var_cut_D[msig], ompos3$var_cut_D)
# depletion 
x=sum(c(47,98,26,56,64))
fisher.test(rbind(c(1,36),c(x,1208-x)))
notPAMGuide=ifelse(ompos3$var_pam_D<0, FALSE,TRUE)

notpg=rownames(ompos3)[!notPAMGuide]

npg=allr[allr$variant %in% notpg,]
pg=allr[!(allr$variant %in% notpg),]


rp=cbind(r$t1$p, r$t2$p)
fmpj=r$joint$p
fmp=apply(rp,1, function(x) { if(sum(is.na(x))==0) {sumlog(x)$p} else {NA}  })
fmp2=apply(rp,1, function(x) { if(sum(is.na(x))==0) {x[1]*x[2]} else {NA}  })

r=results

rc=do.call('cbind', r)
rcsig=rc[which(rc$joint.fdr<.05),]



        SAF=rep(xdf$oligo, ncol(Y))
        bcode=rep(xdf$barcode,ncol(Y))
        #not log
        mm=data.frame(y=Ys$values,variant=SAF, bcode=bcode, condition=Ys$ind2, experiment=Ys$ind, off=ncounts2[as.character(Ys$ind)], transformation=Ys$ind3)
        #mm$transformation=as.factor(tstrsplit(mm$experiment, '_')[[2]])
        
        #log
        #mm=data.frame(y=log(Ys$values+.5),variant=SAF, bcode=bcode, condition=Ys$ind2, experiment=Ys$ind, off=ncounts2[as.character(Ys$ind)])
        # weight by counts 
        #w=(split(Ys$value,mm$bcode))
        #w=1/(sapply(w, function(x) sum(1/(x+.5))))
        #mm$weight=w[mm$bcode]
        
        #
        ggplot(mm, aes(x=condition, y=log((y+.5)/off)))+facet_wrap(~transformation+variant)+geom_boxplot()+geom_jitter(aes(color=bcode))

        if(length(levels(mm$variant))<2 | length(levels(mm$condition))<2 | length(levels(mm$bcode))<2 ) { return(nvec) } else { 
            #ps[oo]=1} else {
            mmcontrols=mm[mm$condition==gsub('_' ,'', vcont),]
            #thetaN=glm.nb(y~variant+offset(log(off)), data=mmcontrols)$theta
            # dispersion 
            # weighted
            #mer=lmer(y~condition+variant+condition*variant+offset(log(off))+(1|bcode),data=mm, weights=(weight))
            #gmer=glmer(y~condition+variant+condition*variant+offset(log(off))+(1|bcode),data=mm, family=negative.binomial(theta=thetaN))
            #gmer=glmer.nb(y~condition+variant+condition*variant+offset(log(off))+(1|bcode),data=mm) #, family=negative.binomial(theta=thetaN))
            #drop1(gmer, test='Chisq')
            #gmer@theta
            mm2=mm[mm$transformation=='2',]
            gtb=glmmTMB(y~-1+variant+condition+variant*condition+offset(log(off))+(1|bcode)+(1|experiment),
                        data=mm, family=list(family="nbinom2", link="log"), dispformula=~transformation) # , ziformula=~1)


            gtb=glmmTMB(y~-1+variant+condition+variant*condition+offset(log(off))+(1|bcode)+(1|experiment),data=mm2, family=list(family="nbinom2", link="log")) # , ziformula=~1)

            gtb=glmmTMB(y~-1+variant+condition+variant*condition+offset(log(off))+(1|bcode)+(1|experiment),data=mm, family=list(family="nbinom2", link="log")) # , ziformula=~1)
            py=as.vector(predict(gtb, type='response'))
            #plot(log2(getME(gmer, 'mu')), log(py))
            #cor(log2(py+1),log2(mm$y+1))^2
            #overdispersion
           
            disp=sigma(gtb)
            mu=sum(log(py))
            ps=drop1(gtb, test='Chisq')[[4]][2]
            
            #variance =μ(1 +μ/k):  Hardin and Hilbe(2007))
            #drop1(gtb, test='Chisq')

            # family=poisson)
            #drop1(gmer, test='Chisq')[[4]][2]
            # not weighted
            #mer=lmer(y~condition+variant+condition*variant+offset(log(off))+(1|bcode),data=mm)
            conds=paste(mm$condition, mm$variant, sep='_')
            gmeans=sapply(split(log2(py), conds), median)
            SA_Cond=gmeans[paste0(tt, 'SA')]
            S_Cond=gmeans[paste0(tt, 'S')]
            SA_Cont=gmeans[paste0(vcont, 'SA')]
            S_Cont=gmeans[paste0(vcont, 'S')]
            l2fc=(as.numeric((SA_Cond-S_Cond)-(SA_Cont-S_Cont)))
           
            print(paste(o, mu, round(disp,2), round(l2fc,2), ps, ifelse(ps<.0001, "***", "")))
            return(list(mu=mu, disp=disp, l2fc=l2fc, p=ps))
        }
    }
    names(o.out)=oncore
    results=rbindlist(o.out,idcol='variant')
    results=data.frame(results, fdr=p.adjust(results$p, method='fdr'))
    
    if(vcont =='PRE_') { vsPre[[tt]]=results  } # list(ps=ps, l2fc=l2fc)}
    if(vcont =='zero_'){ vsZero[[tt]]=results }#list(ps=ps, l2fc=l2fc)}
    if(vcont =='YPD_') { vsYPD[[tt]]=results } #list(ps=ps, l2fc=l2fc)}
}



par(mfrow=c(6,1))
for(i in 1:length(vsZero)){
    hist(vsZero[[i]]$p, breaks=10, main=names(vsZero)[i])
}

x11()
par(mfrow=c(6,1))
for(i in 1:length(vsZero)){
    plot(vsZero[[i]]$l2fc, col=ifelse(vsZero[[i]]$fdr<.2, 'red', 'black'),
         cex=ifelse(vsZero[[i]]$fdr<.2, 3, 1), main=names(vsZero)[i])
}
i=3
vsZero[[i]][which(vsZero[[i]]$fdr<.25),]


#save(vsPre, file='/data/rr_variant_oligos/CC_Hiseq081318/vPre_W.RData')
#save(vsZero, file='/data/rr_variant_oligos/CC_Hiseq081318/vZero_W.RData')
#save(vsYPD, file='/data/rr_variant_oligos/CC_Hiseq081318/vYPD_W.RData')

#save(vsPre, file='/data/rr_variant_oligos/CC_Hiseq081318/vPre.RData')
#save(vsZero, file='/data/rr_variant_oligos/CC_Hiseq081318/vZero.RData')
#save(vsYPD, file='/data/rr_variant_oligos/CC_Hiseq081318/vYPD.RData')

save(vsPre, file='/data/rr_variant_oligos/CC_Hiseq081318/vPre_NB.RData')
save(vsZero, file='/data/rr_variant_oligos/CC_Hiseq081318/vZero_NB.RData')
save(vsYPD, file='/data/rr_variant_oligos/CC_Hiseq081318/vYPD_NB.RData')

#load('/data/rr_variant_oligos/CC_Hiseq081318/vPre.RData')
#load('/data/rr_variant_oligos/CC_Hiseq081318/vZero.RData')
#load('/data/rr_variant_oligos/CC_Hiseq081318/vYPD.RData')

load('/data/rr_variant_oligos/CC_Hiseq081318/vPre_NB.RData')
load('/data/rr_variant_oligos/CC_Hiseq081318/vZero_NB.RData')
load('/data/rr_variant_oligos/CC_Hiseq081318/vYPD_NB.RData')

#load('/data/rr_variant_oligos/CC_Hiseq081318/vPre_W.RData')
#load('/data/rr_variant_oligos/CC_Hiseq081318/vZero_W.RData')
#load('/data/rr_variant_oligos/CC_Hiseq081318/vYPD_W.RData')

vsPreM=do.call('cbind', lapply(vsPre, function(x) x$ps))
vsZeroM=do.call('cbind',  lapply(vsZero, function(x) x$ps))
vsYPDM=do.call('cbind',  lapply(vsYPD, function(x) x$ps))

vsPreFC=do.call('cbind', lapply(vsPre, function(x) x$l2fc))
vsZeroFC=do.call('cbind',  lapply(vsZero, function(x) x$l2fc))
vsYPDFC=do.call('cbind',  lapply(vsYPD, function(x) x$l2fc))


#glmer.nb(y~condition+variant+offset(log(off))+1|bcode)
#my.contrast=makeContrasts(condVypd=(SA_Cond-S_Cond)-(SA_YPD-S_YPD), levels=design.mat)

#boxplot(log2((mm$y+1)/(mm$off))~mm$condition+mm$variant)

par(mfrow=c(ncol(vsYPDFC),1))
for(i in 1:ncol(vsYPDFC)) {
    plot(vsYPDFC[,i], col=ifelse(vsYPDM[,i]<5e-5, 'red', 'black'), main=colnames(vsYPDFC)[i])
}
par(mfrow=c(ncol(vsPreM),1))
for(i in 1:ncol(vsPreM)) {
    plot(vsPreFC[,i], col=ifelse(vsPreM[,i]<5e-5, 'red', 'black'), main=colnames(vsPreM)[i])
}

par(mfrow=c(ncol(vsZeroM),1))
for(i in 1:ncol(vsZeroM)) {
    plot(vsZeroFC[,i], col=ifelse(vsZeroM[,i]<5e-5, 'red', 'black'), main=colnames(vsZeroM)[i])
}
agg=list(vsPreM=vsPreM, vsPreFC=vsPreFC, vsZeroM=vsZeroM, 
              vsZeroFC=vsZeroFC, vsYPDM=vsYPDM, vsYPDFC=vsYPDFC)
agg=lapply(agg,data.frame)
WriteXLS(agg, "/data/rr_variant_oligos/CC_Hiseq081318/NegativeBinomial.xls",      row.names=T)





par(mfrow=c(ncol(vsYPDFC),1))
for(i in 1:ncol(vsYPDFC)) {
    plot(vsYPDFC[,i], col=ifelse(vsYPDM[,i]<.01, 'red', 'black'), main=colnames(vsYPDFC)[i])
}




library(stringdist)
library(fastcluster)
require(venneuler)

mo=c()
for(i in 11:20 ) {
    e1=bdf[[i]]
    e2=exp_r[[i]]
    u1=c(e1$barcode1, e2$barcode1)
    ubc=rle(sort(as.vector(u1)))
    zz=ubc$values
    #[ubc$lengths>2]
    kd=stringdistmatrix(a=zz, method='lv')
    hc=fastcluster::hclust(kd)
    b.groups=cutree(hc, h=4)
    length(unique(b.groups))
    ad=data.frame(barcode=ubc[[2]], barcode.count=ubc[[1]], bgroup=b.groups, stringsAsFactors=F)
    ad=ad[order(ad$bgroup, ad$barcode.count),]
  
    e1.vec= ad$bgroup[match(e1$barcode1, ad$barcode)]
    e2.vec= ad$bgroup[match(e2$barcode1, ad$barcode)]

    sum(e2.vec %in% e1.vec)/length(e2.vec)
    mo=c(mo, sum(unique(e2.vec) %in% unique(e1.vec))/length(unique(e2.vec)))

}


#bam.out=paste0(t,'.bam')
#    txt.out=paste0(t,'.txt')
#
#    system(paste('cutadapt -u 20 -o', t2.out, t1.out))
#    system(paste('bwa mem -t 70 /data/CRISPR_variant_engineering/rr_variant_oligos/results/coding_variants_LOD30.fasta', t2.out, '| samtools view -hSu - | samtools sort - -o',  bam.out))
#    system(paste('samtools index', bam.out))
#    system(paste('samtools idxstats', bam.out, '>',txt.out))
#}

library(data.table)
resultsL=list()
for(t in traits) {
    txt.out=paste0(wd,t,'.txt')
    print(t)
    x=read.delim(txt.out, header=F, stringsAsFactors=F)
    x=x[-nrow(x),]
    xn=do.call('rbind', strsplit(x[,1], ':'))
    nn=paste(xn[,1], xn[,2], sep=':')
    unn=split(x, (nn))
    unn=unn[sapply(unn, nrow)==2]
    sa=(sapply(unn, function(x) (x[2,3]+1)))
    s=(sapply(unn, function(x) (x[1,3]+1)))
    df=data.frame(trait=t, do.call('rbind', strsplit(names(unn), ':|_')), S=s, SA=sa, stringsAsFactors=F)
    df[,2]=as.character(df[,2])
    df[,3]=as.numeric(df[,3])
    names(df)[2:4]=c('chr', 'index', 'variant')
    df$pos=tstrsplit(rownames(df), '_|:', type.convert=T)[[2]]
    df$log2_E=log2(df$SA/df$S)
    df$chr=factor(df$chr, levels=paste0('chr', as.roman(1:16)))
    tot.S=sum(df$S)
    print(tot.S)
    print(tot.SA)
    tot.SA=sum(df$SA)
#    df$ps=rep(NA, nrow(df))
#    for(i in 1:nrow(df)){
#        mm=t(matrix(as.vector(unlist(c(df[i,cS], tot.S-df[i,cS],c(df[i,cSa],tot.SA-df[i,cSa])))),2,2))
#        df$ps[i]=fisher.test(mm)$p.value
#    }
#    df$ps.adjusted=p.adjust(df$ps, method='bonferroni')
    resultsL[[t]]=df
}
library(dplyr)

results=bind_rows(resultsL)
#results$log2_E=log2(results$SA/results$S)
#results$chr=factor(results$chr, levels=paste0('chr', as.roman(1:16)))
#results=results[,c(1,2,7,4,5,6,8,9,10),]
#results$nlog10p=-log10(results$ps.adjusted)
testD=do.call('cbind', lapply(resultsL, function(x) x[,c(5,6,8)]))
testD=cbind(resultsL[[1]][,c(2,3,4,7)], testD)
tD=(testD[,grep('log2', colnames(testD))])
tD=tD[,c(1,2,5,7,3,4,6,8)]


tD2=(testD[,grep('.S$', colnames(testD))])
tDS=tD2[,c(1,2,5,7,3,4,6,8)]

tD3=(testD[,grep('.SA$', colnames(testD))])
tDSA=tD3[,c(1,2,5,7,3,4,6,8)]


panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y,method='spearman'), digits=2)
    txt <- paste0("rho = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(tD, lower.panel=panel.cor)

to=test[order(abs(test$mnnit.log2_E), decreasing=T),]


test[abs(test$mnnit.log2_E-test$yn1.log2_E)>2 & test$yn1.S>20 & abs(test$mnwt.log2_E-test$yw1.log2_E)>2 & test$yw1.S>20 ,]
test[abs(test$mnnit.log2_E-test$yn1.log2_E)>6 & test$yn1.S>30,]

ggplot(results, aes(x=pos, y=log2_E))+  geom_point(size=.3) +  theme_bw() + facet_grid(trait ~ chr,scales = "free_x", space='free_x') + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(size = 5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())


library(edgeR)
r0=resultsL[['yw1']]
r1=resultsL[['yw2']]
#r0=resultsL[['yn1']]
#r1=resultsL[['yn2']]
#r2=resultsL[['yn1']]
#r3=resultsL[['yn2']]
#r2=resultsL[['15C']]
disp.set=DGEList(counts=cbind(r0[,c(5:6)],r1[,c(5:6)]),
                 group=c('S', 'SA', 'S', 'SA')) #, 'S', 'SA','S', 'SA') )
                 #,r2[,c(5:6)], r3[,c(5:6)]), 

hist(aveLogCPM(disp.set$counts, offset = getOffset(disp.set)))

disp.set=calcNormFactors(disp.set)
hist(rowSums(log(cpm(disp.set))/6))

design.mat=model.matrix(~0+disp.set$samples$group)
disp.set=estimateGLMCommonDisp(disp.set, verbose=T)
disp.set=estimateGLMTrendedDisp(disp.set)
disp.set=estimateGLMTagwiseDisp(disp.set)
plotBCV(disp.set, cex=1)
#tin=.1392
tin=1.709
tin=1.3
rmm=rowMeans(cpm(disp.set))
plot(log2(rmm), log2(apply(cpm(disp.set), 1,var)), xlab='log2(Mean)', ylab='log2(Var)')
abline(0,1, col='blue')
points(log2(rmm), log2(rmm+tin*rmm^2), col='red' , pch=20, cex=.5)
legend('topleft',c(expression(sigma^2==mu), 
                   expression(sigma^2==mu+phi*mu^2)),text.col=c('blue', 'red'),cex=1.25)

#test=estimateTagwiseDisp(test)
#d2 <- estimateGLMCommonDisp(test,design.mat)
#glmFit(d2, design.mat)
#$dispersion
#[1] 0.06879
wt.bck=names(resultsL)[c(5,7)]
#nt.bck=names(resultsL)[c(6,8)]
hacked.stats.wt=list()
#hacked.stats.nt=list()
for(t in wt.bck) { 
#for(t in nt.bck) { 
    print(t)
    #-c(1,2)] ) {
    r3=resultsL[[t]]
    r2=data.frame(r3,r0[,c(5:6)], r1[,c(5:6)])
    names(r2)[c((ncol(r2)-3):ncol(r2))]=c('S_YPD', 'SA_YPD','S_YPD', 'SA_YPD')
    test=DGEList(counts=cbind(r3[,c(5:6)],r0[,c(5:6)],r1[,c(5:6)]), group=c( 'S_Cond', 'SA_Cond', 'S_YPD', 'SA_YPD','S_YPD', 'SA_YPD'))
    test=calcNormFactors(test)
    design.mat=model.matrix(~0+test$samples$group)
    colnames(design.mat)=levels(test$samples$group)
    tfit=glmFit(test, design.mat, dispersion=disp.set$common.dispersion)
    #tfit=glmFit(test, design.mat, dispersion=disp.set$tagwise.dispersion)
    #tfit=glmFit(test, design.mat, dispersion=disp.set$trended.dispersion)
    my.contrast=makeContrasts(condVypd=(SA_Cond-S_Cond)-(SA_YPD-S_YPD), levels=design.mat)
    #hacked=glmLRT(tfit, contrast=c(1,1,-1,-1))
    #hacked1=glmLRT(tfit, contrast=c(1,-1,1,-1))
    #hacked=glmQLFit(test, design.mat,dispersion=disp.set$common.dispersion)
    hacked=glmLRT(tfit, contrast=my.contrast)
    r3=data.frame(r2, hacked$table)
    r3$PValue[r3$S<10 | r3$S_YPD<10]=1
    r4=DataFrame(r3, ccoding[match(rownames(r3), names(ccoding)),])
    r4=r4[order(r4$PValue, decreasing=F),]
    r4$p.adjusted=p.adjust(r4$PValue, method='BH')
#    hacked.stats.nt[[t]]=r4
    hacked.stats.wt[[t]]=r4
}
hs=lapply(hacked.stats.nt, data.frame)
WriteXLS(hs, '~/Desktop/nit1_hiseq.xls')

hs2=lapply(hacked.stats.wt, data.frame)
WriteXLS(hs2, '~/Desktop/wt_hiseq.xls')

lapply(


    
    r4=head(r3[order(r3$PValue, decreasing=F),], 20)

    head(hacked$table[order(hacked$table$PValue, decreasing=F),])




glmFit(test, 
hacked=(exactTest(test, pair=c(dispersion=.1))




#results$subject=apply(results,1, function(x) paste( c(x[2], x[3], x[4]), collapse='_', sep=''))
#results$subject=gsub(' ', '', results$subject)
#tg=gather(results, "type", "count", S, SA)
#stg=tg[tg$trait=='YPD',]

#library(DESeq2)
#DESeqDataSetFromMatrix(stg,  











par(mfrow=c(6,1))
for(t in traits){
    dff=results[[t]]

    with(dff, {
        plot(log10(SA/S),type='p',yaxt='n',ylab='Enrichment of non-synonymous variant', xlab='snp #')
        axis(2, at=c(-2,-1,0,1,2,3,4), c(-100,-10,0,10,100,1000,10000),las=2 )
        abline(h=0)
        abline(h=c(-2,2), lty=2)
    #points(92, 4.09, col='red', pch=20, cex=3)
    })
}


ggplot(dff,aes(x = position,y=LOD)) + geom_point(size=.3) + ylim(c(0,15)) +
      theme_bw() + facet_grid(trait ~ chr,scales = "free_x", space='free_x') + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(size = 5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
      geom_vline(aes(xintercept=pos_fineMapping, colour=QTN), data=sct, linetype=1, size=.6, alpha=.6 )+
      geom_point(aes(x=pos_fineMapping, y=0, colour=dup), size=2, shape=17, data=sct )+scale_color_manual(values=c('FALSE'="blue", 'Identical'='green', 'TRUE'='red'))+
      ggtitle(paste(trait, ' - ' ,nrow(sct),  'total QTL (', sum(sct$peakWidth==1), 'QTN) identified by She and Jarosz', sep=' '))
      

