# build whole genome dna string set with reverse complement
buildSacCer3_genome_dictionaryFR=function(sacCer3, unique.chrs) {
    sc3.dna=list()
    for(chr in  unique.chrs ) { 
        sc3.dna[[chr]]=DNAString(sacCer3[[chr]])
        rc.name=paste0(chr, '_rc')
        sc3.dna[[rc.name]]=reverseComplement(DNAString(sacCer3[[chr]]))
    }

    sc3.set=DNAStringSet(sc3.dna)
}

# build pdict for matching guide rna seed sequences to the genome
# includes GG of PAM sequence
# 1-12 left of the PAM (not including PAM)
build_4base_pdict=function(guide.seq) { 
    guide.seeds.N13.pdict=list()
    for(base in c('A', 'C', 'G', 'T') ) {
        gtable.seq=sapply(guide.seq, s2c)[9:23,]
        gtable.seq[13,]=base
        all.guide.seeds=as.vector(apply(gtable.seq,2, c2s))
        seed.set=DNAStringSet(unlist(all.guide.seeds))
        guide.seeds.N13.pdict[[base]]=PDict(seed.set)
    }
    return(guide.seeds.N13.pdict)
}


#Locate PAM sequences and predict the effect of the base  editor
#every targetable PAM site
# if returnIdentical then only output gRNAs that target a region without Cs
get_gRNAs_and_predict_edit=function(nchr, sacCer3, txdb, sc3.set ) {

    PAMp_string=DNAString('GG')
    PAMm_string=reverseComplement(PAMp_string)
    chr.name=names(sacCer3)[nchr]

    PAMp=matchPattern(PAMp_string, sacCer3[[nchr]], fixed=T)
    PAMm=matchPattern(PAMm_string, sacCer3[[nchr]], fixed=T)
    
   # construct granges object 
    gPAM = GRanges(
                   seqnames=names(sacCer3)[nchr], 
                   ranges=IRanges(start=c(start(ranges(PAMp))-1, start(ranges(PAMm))), 
                                           end=c(end(ranges(PAMp)), end(ranges(PAMm))+1)),
                       strand=c(rep('+', length(PAMp)),rep('-', length(PAMm)) ) 
                    )
    genome(gPAM)='sacCer3'
    gPAM$guide=rep('', length(gPAM))

    too.close.to.ends=which( ((start(gPAM)-50)<0 ) | ((end(gPAM)+50)> length(sacCer3[[nchr]])) )
    gPAM=gPAM[-too.close.to.ends]

    # extract PAM sequence and targetting sequence (convert CCN to NGG)
    # extract targetting sequence only
    soi=as.vector(strand(gPAM)=='+')
    gPAM$guide[soi]=getSeq(sacCer3,chr.name, 
                           start=start(gPAM)[soi]-20,
                           end=start(gPAM)[soi]+2,
                           as.character=T) 
   
    gPAM$repair[soi]=getSeq(sacCer3,chr.name, 
                           start=start(gPAM)[soi]-43,
                           end=start(gPAM)[soi]+46,
                           as.character=T) 


    soi=as.vector(strand(gPAM)=='-')
    gPAM$guide[soi]=as.character(reverseComplement(DNAStringSet(getSeq(sacCer3,chr.name, 
                           start=end(gPAM)[soi]-2,
                           end=end(gPAM)[soi]+20,
                           as.character=T) )))

    #this bit is confusing
    gPAM$repair[soi]=as.character(reverseComplement(DNAStringSet(getSeq(sacCer3,chr.name, 
                           start=start(gPAM)[soi]-44,
                           end=start(gPAM)[soi]+45,
                           as.character=T)))) 


    # modify gPAM such that ranges match where the sequences start and stop
    start(ranges(gPAM[strand(gPAM)=='+']))=start(ranges(gPAM[strand(gPAM)=='+']))-20
    #this isn't fixed yet
    end(ranges(gPAM[strand(gPAM)=='-']))=end(ranges(gPAM[strand(gPAM)=='-']))+20
    
    # edit all Cs to Ts for positions 4-8 
    gPAMseq=t(sapply(gPAM$guide, s2c))
    #for(p in 4:8) { gPAMseq[gPAMseq[,p]=='C' ,p]='T' }

    gPAMrepair=t(sapply(gPAM$repair, s2c))

    gPAMrepL=gPAMseqL=list()
    #targeted changes 

    Nchanges=c('A', 'C', 'T', 'G')
    GGchanges=c('CG', 'TG', 'CA', 'TA', 'AA', 'GA', 'AC', 'CC', 'TC', 'GC', 'CT', 'TT', 'AT', 'GT')
    
    for(g in Nchanges) {
            gPAMseqtemp=gPAMseq
            gPAMrepairtemp=gPAMrepair
            for(i in 1:length(GGchanges)){
                gPAMrepairtemp[,44]=gPAMseqtemp[,21]= g
                gPAMrepairtemp[,45]=gPAMseqtemp[,22]= s2c(GGchanges[i])[1]
                gPAMrepairtemp[,46]=gPAMseqtemp[,23]= s2c(GGchanges[i])[2]
                gPAMseqL[[ paste0(g,GGchanges[i]) ]]=gPAMseqtemp
                gPAMrepL[[ paste0(g,GGchanges[i]) ]]=gPAMrepairtemp
            }
    }

    editedSequencesL=lapply(gPAMseqL, function(x) apply(x, 1, c2s))


    editedSequences.stringsetL=lapply(editedSequencesL, function(x) DNAStringSet(x))
    editedSequences.stringset.fwdL=editedSequences.stringsetL
    editedSequences.stringset.fwdL=lapply(editedSequences.stringset.fwdL, 
           function(x) {
               x[as.vector(strand(gPAM)=='-')]=reverseComplement(x[as.vector(strand(gPAM)=='-')])
               return(x)
           })

    
    repSequences=lapply(gPAMrepL, function(x) apply(x, 1, c2s))
    repSequences=lapply(repSequences, function(x) DNAStringSet(x))

    repSequences=lapply(repSequences, 
           function(x) {
               x[as.vector(strand(gPAM)=='-')]=reverseComplement(x[as.vector(strand(gPAM)=='-')])
               return(x)
           })
   


    reffull=getSeq(sacCer3,chr.name, 
                           start=start(gPAM),
                           end=end(gPAM),
                           as.character=T)
    reffull=t(sapply(reffull, s2c))

    
    gEditInfo=lapply(editedSequences.stringset.fwdL, function(x) {
        altfull=t(sapply(as.character(x), s2c))
        gPOS=rep(0,length(gPAM))
        gREF=rep('',length(gPAM))
        gALT=rep('', length(gPAM))         
        for(i in 1:nrow(reffull)){
            nref.range=range(which(reffull[i,]!=altfull[i,] ))
            gPOS[i]=min(nref.range)
            gREF[i]=(c2s(reffull[i, seq(nref.range[1], nref.range[2])   ])) 
            gALT[i]=(c2s(altfull[i, seq(nref.range[1], nref.range[2])   ])) 
        }
        return(paste0(chr.name, ':', start(gPAM)+gPOS-1,':', gREF,':', gALT))
    })



# turn back into strings
#   editedSequences=apply(gPAMseq, 1, c2s)

    # reverse complement 
#    editedSequences.stringset=DNAStringSet(editedSequences)
#    editedSequences.stringset.fwd=editedSequences.stringset
#    editedSequences.stringset.fwd[as.vector(strand(gPAM)=='-')]=reverseComplement(editedSequences.stringset.fwd[as.vector(strand(gPAM)=='-')])
    
    # ranges and editedSequences.stringset correspond to sequence on forward strand
    # so change strand to forward
    strand(gPAM)[strand(gPAM)=='-']='+'

    # add in expected edited sequence as a DNA string set (all sequence on forward strand)
    #gPAM$forward.edited.sequence=editedSequences.stringset
    editedSequences.stringsetL=as.data.frame(editedSequences.stringsetL)
    editedSequences.stringset.fwdL=as.data.frame(editedSequences.stringset.fwdL)

    repairTemplate=as.data.frame(repSequences)
    #we never use this, so commenting it out
    #gPAM$editedSequenceAsIsL=editedSequences.stringsetL
    gPAM$repairTemplate=repairTemplate
   
    gPAM$editedSequenceFwdL=editedSequences.stringset.fwdL
    gPAM$gEditInfo=as.data.frame(gEditInfo)
    # remove guides that do not introduce any changes
   # identical=which(gPAM$guide == gPAM$editedSequenceAsIs)
   # if(returnIdentical) {
   #     gPAM=gPAM[identical]
   #     # create an index for each guide sequence ... note that the guideIndex is chromosome specific
   #     gPAM$guideIndexNoChange=paste0(chr.name, ':', seq(1, length(gPAM)))
   # } else {
   #     gPAM=gPAM[-identical]
   #     # create an index for each guide sequence ... note that the guideIndex is chromosome specific
        gPAM$guideIndex=paste0(chr.name, ':', seq(1, length(gPAM)))
   # }

    pdict_for_guide_seqs=build_4base_pdict(gPAM$guide)
    guide.12PM.genome.count.list=lapply(pdict_for_guide_seqs, function(seed.dict) {
                    vcountPDict(seed.dict,  sc3.set, max.mismatch=0, with.indels=F)  
                  })
    # describe this : guidePAM12PMcount
    gPAM$guidePAM12PMcount=rowSums(sapply(guide.12PM.genome.count.list, rowSums))
    gPAM$repair=NULL
    return(gPAM)
}


mungeCodingEffects=function(coding.effects, ns, meltedProveanScores) {
    goi=coding.effects$GENEID
    refAA=sapply(as.character(coding.effects$REFAA), s2c)
    varAA=sapply(as.character(coding.effects$VARAA), s2c)
    pdiff=mapply(function(x,y)which(x!=y), refAA,varAA)
    loi = mapply(function(x,y) {
                     x+y-1}, 
                     x=sapply(coding.effects$PROTEINLOC, function(x)x[1]),
                    y=pdiff)
    varoi=mapply(function(x,y) {
                     x[y]
                    },x=varAA,y=pdiff)
    refoi=mapply(function(x,y) {
                     x[y]
                    },x=refAA,y=pdiff)

    lookupchange=mapply(function(g,l,v){
               paste0(g,':',l, ':', v)
            #   meltedProveanScores$value[meltedProveanScores$gene==g & meltedProveanScores$ALT==v & meltedProveanScores$POS==l]
                    }, goi, loi, varoi)


    provScore=relist(meltedProveanScores$value[match(unlist(lookupchange), meltedProveanScores$lookup)],lookupchange)
    coding.effects$proveanEffects=provScore
    coding.effects$POSAA.delta=loi
    coding.effects$REFAA.delta=refoi
    coding.effects$VARAA.delta=varoi
    coding.effects$gEdit= coding.effects$gEditInfo[,ns]
    coding.effects$gEditSeqFwd=coding.effects$editedSequenceFwdL[,ns]
    coding.effects$repairTemplate=coding.effects$repairTemplate[,ns]
    return(coding.effects[,c('guideIndex','gEdit', 'gEditSeqFwd', 'repairTemplate', 'GENEID',  'POSAA.delta', 'REFAA.delta', 'VARAA.delta', 'proveanEffects', 'CONSEQUENCE')])
}





#Add addtional information to guides tables 
augment_guides=function(guides, addScores=FALSE) {
     # GC content 
    lf=letterFrequency(DNAStringSet(guides$guide), letters='ACGT', OR=0)
    guides$guideGCcontent=(lf[,2]+lf[,3])/rowSums(lf)
       
    # Contains putative pol3 terminator
    guides$U6terminator=grepl('TTTT', guides$guide)

    ## write.table(guides$guide, file='/media/jbloom/d1/CRISPR_base_editor/results/guides_to_score.txt', quote=F, row.names=F, col.names=F)
    # guide scores ... for example 'http://www.flyrnai.org/evaluateCrispr/' 
    if(addScores) {
        guides$guideScore = read.csv(gzfile(paste0(working.dir, 'results/guide_scores.csv.gz')), as.is = TRUE)$Score
    }
    # at this point `guides` is finished, could append additional metrics from BLAST see the designCRISPR_gRNAs_fx.R for more information
    return(guides)
}



# Run Local Blast
do_local_blast=function(guide.tables.df, working.dir) {
    ##! run once locally on the command line in the reference/ folder
    # makeblastdb -in sacCer3.fasta -dbtype 'nucl'
    bl = makeblast(db=paste0(working.dir, "reference/sacCer3.fasta"))
    # adjust number of threads if not running on Josh's computer (4)
    blast_all_guides=predictBLAST(bl, DNAStringSet(guide.tables.df$guide), 
                             BLAST_args='-evalue=100 -word_size=8 -gapopen=5 -gapextend=2 -reward=2 -penalty=-3 -num_threads=70',
                              custom_format='sseq qseq'
                              )
    save(blast_all_guides, file =paste0(working.dir, "results/BLAST_gRNAs.RData"))
    return(blast_all_guides)
}

# Filter Blast Junk
filter_blast_junk=function(blast_all_guides, guide.tables.df) {
    gi=paste0('Query_', 1:length(DNAStringSet(guide.tables.df$guide)) ) #length(guideRC))
    blast_all_guides$ID=match(as.character(blast_all_guides$QueryID), gi)
        ## EEK!
        #R> nrow(blast_all_guides)
        ##[1] 332,235,725
    blast_all_guides=blast_all_guides[blast_all_guides$Q.end==23 & blast_all_guides$Alignment.Length>14 & grepl('GG$', blast_all_guides$sseq) ,]
        #R> nrow(blast_all_guides)
        ##[1] 7,563,337
    save(blast_all_guides, file =paste0(working.dir, "results/BLAST_gRNAs_reduced.RData"))
    return(blast_all_guides)
}

# Catalog potential off-target activitity using BLAST --------------------------------------------------------------------------------------
    ##source(paste0(working.dir, 'code/BLAST/BLAST.R'))
    ##blast_all_guides=do_local_blast(guides, working.dir)
    ##load(paste0(working.dir, "results/BLAST_gRNAs.RData"))
    ##blast_all_guides=filter_blast_junk(blast_all_guides, guides)
    
    # Recommend starting here with blast output analysis
    #load(paste0(working.dir, "results/BLAST_gRNAs_reduced.RData"))

    #qseq=DNAStringSet(blast_all_guides$qseq)
    #sseq=DNAStringSet(blast_all_guides$sseq)

    # perfect match of 12 bp seed region 
    #qseq12=subseq(qseq, start=width(qseq)-14, end=width(qseq)-3)
    #sseq12=subseq(sseq, start=width(sseq)-14, end=width(sseq)-3)
    #check12=qseq12==sseq12
    #blast_all_guides=blast_all_guides[check12,]

    #filter only hits where the seed sequence is identical
    # lost 563 due to blast fail
    #bcount=rle(blast_all_guides$ID)
    #which(!(1:nrow(guides) %in% bcount$values))
#--------------------------------------------------------------------------------------------------------------------------------------------

# identify yeast pseudogenes
get_pseudogenes=function(working.dir, orf_coding) {
    orf_coding_all=read.fasta(paste0(working.dir, 'reference/orf_coding_all.fasta.gz'))
    orf_coding_dubious=read.fasta(paste0(working.dir, 'reference/orf_coding_dubious.fasta.gz'))

    n1=names(orf_coding_all)[!(names(orf_coding_all) %in% names(orf_coding))]
    pseudogenes=n1[!(n1 %in% names(orf_coding_dubious))]
    return(pseudogenes)
}

# reannotate nonsynonymous PTCs, note that in theory we could get VARAA ending in PTC and REFAA not, but this never occurs
annotate_PTCs=function(coding.effects){
        ## for example:   sum(grepl('[A-Z]$', as.character(coding.effects$REFAA)) & grepl('\\*$', as.character(coding.effects$VARAA)))
        ptcs = grepl('\\*[A-Z]', coding.effects$VARAA)
        coding.effects$CONSEQUENCE=as.character(coding.effects$CONSEQUENCE)
        coding.effects$CONSEQUENCE[ptcs]='stop_gained'
        coding.effects$CONSEQUENCE=as.factor(coding.effects$CONSEQUENCE)
        return(coding.effects)
}


#deprecated
# return list of provean scores per gRNA
get_provean_scores=function(coding.effects, provean.score.file) {
  
    #loads provean_effects (list of matrices)
    load(provean.score.file)
    provean_effects=list()
    
    # this is terribly slow and should be vectorized
    pb=txtProgressBar(min=1, max=nrow(coding.effects), style=3)

    for(i in 1:nrow(coding.effects))  {
        #    if(i%%1000==0) { print(i) }
        setTxtProgressBar(pb,i)
        if(coding.effects$CONSEQUENCE[i]!='synonymous') {
             goi=coding.effects$GENEID[i]
             refAA=as.matrix(coding.effects$REFAA[i])
             varAA=as.matrix(coding.effects$VARAA[i])
             pdiff=which(refAA!=varAA)
             loi=coding.effects$PROTEINLOC[[i]][1]+pdiff-1
             varoi=varAA[pdiff] 
             lmat=cbind(loi,match(varoi, colnames(provean_scores[[goi]])))
             # sanity checking
             # refoi=s2c(as.character(coding.effects$REFAA[i]))[pdiff] 
              prot.seq.rownames=rownames(provean_scores[[goi]])
             # prot.seq.rownames[loi]
             # to deal with an exception for HKR1
             if(!is.null(prot.seq.rownames)){
              provean_effects[[as.character(i)]]=c(provean_scores[[goi]][lmat] )
             } else {provean_effects[[as.character(i)]]=numeric(0) }
        } else {
            provean_effects[[as.character(i)]]=numeric(0)
        }
    }
    close(pb)
    return(provean_effects)
}








#also deprecated
# late night attempt at provean lookup speedup
#inputs are list of AAstringsets
#x=coding.effects$REFAA[1:10000]
#y=coding.effects$VARAA[1:10000]
get_AA_mismatches= function(x, y, ret.y=TRUE)
{
    x_width = width(x)
    y_width = width(y)
    ux=unlist(x)
    uy=unlist(y)
    mx=as.matrix(ux)
    my=as.matrix(uy)
    unlisted_ans = which(as.raw(ux) != as.raw(uy))
    breakpoints = cumsum(x_width)
    ans_eltlens = tabulate(findInterval(unlisted_ans - 1L, breakpoints) + 1L,  nbins=length(x))
    skeleton = PartitioningByEnd(cumsum(ans_eltlens))
    if(ret.y) {
        mm   = relist(my[unlisted_ans,1], skeleton)
    } else {
        mm   = relist(mx[unlisted_ans,1], skeleton)
    }
    return(mm)
    #loc  = relist(unlisted_ans, skeleton)
    #offsets = c(0L, breakpoints[-length(breakpoints)])
    #loc = loc - offsets
}

get_max_absolute_provean= function(coding.effects) {
        hasProveanScores=which(sapply(coding.effects$proveanEffects, length)>0)
        maxAbsProvean=sapply(coding.effects[hasProveanScores,]$proveanEffects, function(x) max(abs(x), na.rm=T))
        maxAbsProvean[!is.finite(maxAbsProvean)]=NA
        coding.effects$maxAbsProvean=NA
        coding.effects$maxAbsProvean[hasProveanScores]=maxAbsProvean
        return(coding.effects)
}
