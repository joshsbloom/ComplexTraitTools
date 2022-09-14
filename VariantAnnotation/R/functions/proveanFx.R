library(data.table)

addProvean=function(coding.effects, meltedProveanScores) {
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
    #coding.effects$gEdit= coding.effects$gEditInfo[,ns]
    #coding.effects$gEditSeqFwd=coding.effects$editedSequenceFwdL[,ns]
    return(coding.effects) #[,c('guideIndex','gEdit', 'gEditSeqFwd', 'GENEID',  'POSAA.delta', 'REFAA.delta', 'VARAA.delta', 'proveanEffects', 'CONSEQUENCE')])
}



meltProvean=function(provean_scores){ 
    meltedProveanScores=lapply(provean_scores, 
                         function(x) {
                             y=reshape::melt(x)
                             names(y)[c(1,2)]=c('REF', 'ALT', 'score')
                             y$POS=rep(1:nrow(x), ncol(x))
                             return(y)
                         })
    meltedProveanScores=data.table::rbindlist(meltedProveanScores, idcol='gene')
    meltedProveanScores$lookup=paste0(meltedProveanScores$gene, ':', meltedProveanScores$POS, ':', meltedProveanScores$ALT)
    return(meltedProveanScores)
}
