#add jitter to a genetic map 
jitterGmapVector=function(themap, amount=1e-6) {
    for (i in 1:length(themap)) {
         n <- length(themap[[i]])
         themap[[i]] <- themap[[i]] + c(0, cumsum(rep(amount, n - 1)))
    }
    return(themap)
}

#convert physical position to genetic map position 
getGmapPositions=function(vcf.cross, gmap, uchr) {
    #get physical position, split by chromosome
    p.by.chr=split(vcfR::getPOS(vcf.cross),vcfR::getCHROM(vcf.cross))
    #keep things sorted (yay yeast chr names with roman numerals)
    p.by.chr=p.by.chr[uchr]

   #where to put the variant sites, impute onto gmap
    imputed.positions=mapply( 
           function(x, y){
                approxfun(y$ppos, y$map, rule=2)(x)
            },
            x=p.by.chr, y=gmap,
            SIMPLIFY=F)

}
createFounderPop=function(vcf, gt, p.names, X.only=F, X.drop=T) { 
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




