# take the larger vcf,  genotype calls, and a subset of parents 
# extracts segregating sites 
# return an alphaSimR founder population
createFounderPop=function(vcf, gt, p.names, X.only=F, X.drop=T, gmap) { 
    gt.sub=gt[,colnames(gt) %in% p.names]

    #monomorphic=apply(gt.sub, 1, function(x) all.equal(x))
    #monomorphic sites 
    #faster to do this with math
    rSg=rowSums(gt.sub)
    #sites with hets 
    sum(is.na(rSg))
    #sites all ref
    sum(rSg==0, na.rm=T)
    #sites all alt
    sum(rSg==length(p.names), na.rm=T)
    #sites mito
    sum(grepl('MtDNA', rownames(gt.sub)))

    bad.sites= is.na(rSg) | rSg==0 | rSg==length(p.names)  | grepl('MtDNA', rownames(gt.sub))
    if(X.only) {  bad.sites = bad.sites | !(grepl('X_', rownames(gt.sub))) }
    if(X.drop) {  bad.sites = bad.sites | (grepl('X_', rownames(gt.sub))) }
    gt.sub=gt.sub[-which(bad.sites),]
    vcf.cross=vcf[match(rownames(gt.sub), rownames(gt)), samples=colnames(gt.sub)]
    #generate sample ID
    vcf.cross=vcfR::addID(vcf.cross)


    uchrU=unique(vcfR::getCHROM(vcf.cross))

    imputed.positions=jitterGmapVector(getGmapPositions(vcf.cross, gmap[uchrU], uchrU)) 
    

    #genetic map positions must be in Morgans
    genMap=data.frame(markerName=paste0(vcfR::getCHROM(vcf.cross),'_',vcfR::getPOS(vcf.cross)), 
                      chromosome=vcfR::getCHROM(vcf.cross), position=unlist(imputed.positions)/100)

    teg.GT=t(gt.sub)
    #recode
    teg.GT[teg.GT==0]=-1
    colnames(teg.GT)=paste0(vcfR::getCHROM(vcf.cross),'_',vcfR::getPOS(vcf.cross))
    ped=data.frame(id=rownames(teg.GT), mother=rep(0, nrow(teg.GT)), father=rep(0,nrow(teg.GT)) ) #c(0,0), father=c(0,0))
    return(AlphaSimR::importInbredGeno(geno=teg.GT, genMap=genMap, ped=ped))
}


#' Phase reference and alt counts given the ref/alt calls for one of the parents, using AlphaSimR objects
#'
#' @param df data.frame output from getBiallelicCounts() should contain the columns: ID,ref,alt
#' @param p1.name name of parent to be designated parent 1 in vcf.cross (name must exist in vcf.cross) 
#' @param founderPop #### 
#' @return data.frame of variant ID, ref count, and alt count phased  
#' @export
phaseBiparental=function(df, p1.name, founderPop, genMap){

        gID=genMap$id
        p1.ref=AlphaSimR::pullMarkerGeno(founderPop, genMap$id)[p1.name,]==0
        p1.ref=p1.ref[gID]
        vname=names(p1.ref)

        p1=c(df$ref[p1.ref], df$alt[!p1.ref])
        vscramb=c(vname[p1.ref], vname[!p1.ref])
        names(p1)=vscramb
        p1=p1[vname]

        p2=c(df$ref[!p1.ref], df$alt[p1.ref])
        vscramb=c(vname[!p1.ref], vname[p1.ref])
        names(p2)=vscramb
        p2=p2[vname]
        if(!is.null(df$expected)) {
            expected.phased=ifelse(p1.ref, df$expected, 1-df$expected)
            df$expected.phased=expected.phased
         }
         df$p1=p1
         df$p2=p2
    return(df)
}

#' Subset GATK table for known segregating variants and phase given parental genotypes
#' 
#' 
makeCountTablesGATK=function(sample_dir, allele_counts, p.names, vcf, gt, gmap) {
	founderPop = createFounderPop(vcf,gt, p.names,X.only=F, X.drop=F, gmap)
	genMap=AlphaSimR::getGenMap(founderPop)
	scounts <- readr::read_tsv(stringr::str_c(sample_dir, "/", allele_counts))

	scounts.sub=scounts[paste0(scounts$contig, '_', scounts$position) %in% genMap$id,]
        scounts=data.frame(id=paste0(scounts.sub$contig, '_', scounts.sub$position),ref=scounts.sub$refCount, alt=scounts.sub$altCount)
        scounts=dplyr::left_join(genMap, scounts, by='id')

        scounts$ref[is.na(scounts$ref)]=0
        scounts$alt[is.na(scounts$alt)]=0
        names(scounts)[1]='ID'

	countdf=phaseBiparental(scounts, p.names[1], founderPop, genMap)
        
        #note, we need a better structure for keeping track of which parent is which 
        attr(countdf, 'p1')=p.names[1]
        attr(countdf, 'p2')=p.names[2]

	return(countdf)

}
