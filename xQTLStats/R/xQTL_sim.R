#crosses.to.parents=list(
#     '375'=c("M22", "BYa"),           #1
#     'A'  =c("BYa", "RMx"),           #2
#     '376'=c("RMx", "YPS163a"),       #3
#     'B'  =c("YPS163a", "YJM145x"),   #4
#     '377'=c("YJM145x", "CLIB413a"),  #5
#     '393'=c("CLIB413a", "YJM978x"),  #6
#     '381'=c("YJM978x", "YJM454a"),   #7
#    '3008'=c("YJM454a", "YPS1009x"),  #8
#    '2999'=c("YPS1009x", "I14a"),     #9
#    '3000'=c("I14a", "Y10x"),         #10
#    '3001'=c("Y10x", "PW5a"),         #11
#    '3049'=c("PW5a", "273614xa"),     #12
#    '3003'=c("273614xa", "YJM981x"),  #13
#    '3004'=c("YJM981x", "CBS2888a"),  #14
#    '3043'=c("CBS2888a", "CLIB219x"), #15
#    '3028'=c("CLIB219x", "M22")       #16
#    )
#usethis::use_data(crosses.to.parents, crosses.to.parents)


#source('/home/jbloom/Dropbox/code/xQTLStats/code/xQTLStats.R')

#' Get a subset of a VCF for two input parents from a biparental cross
#
#' Read in a reference VCF file, assuming variant calls are assuming ploidy = 1, 
#' only extract information at biallelic variants, and ensure that only sites that vary
#' between the two parents are output.
#' 
#' @param ref.vcf Path to reference vcf file
#' @param p1.name string corresponding to parent designated as parent 1 (column name needs to exist in vcf )
#' @param p2.name string corresponding to parent designated as parent 2 (column name needs to exist in vcf )
#' @return a vcfR object for the biallelic variants for that cross 
#' @export
getCrossVCF=function(ref.vcf, p1.name, p2.name) {
    vcf=vcfR::read.vcfR(ref.vcf)
    #just get biallelic variants for now 
    vcf=vcf[vcfR::is.biallelic(vcf),]
    #get genotypes
    gt=vcfR::extract.gt(vcf)
    #extract genos for parents 
    gt.sub=gt[,c(p1.name,p2.name)]
    #find variant sites
    gt.sub=gt.sub[gt.sub[,1]!=gt.sub[,2],]
    #dump sites with missing info for now
    gt.sub=gt.sub[!( is.na(gt.sub[,1]) | is.na(gt.sub[,2])),]
    #remove more complex structural stuff or mitochondrial variants for now
    gt.sub=gt.sub[!grepl('simple|deletion|complex|duplication|chrM', rownames(gt.sub)),]
    #subset vcf ile
    vcf.cross=vcf[match(rownames(gt.sub), rownames(gt)), samples=c(p1.name,p2.name)]
    #generate sno ID
    vcf.cross=vcfR::addID(vcf.cross)
    return(vcf.cross)
}

#' Jitter a genetic map so map positions don't completely overlap
#' 
#' @param themap a list of vectors of marker positions in cM
#' @param amount amount of jitter
#' @return a list of vectors of jittered marker positions
jitterGmapVector=function(themap, amount=1e-6) {
    for (i in 1:length(themap)) {
         n <- length(themap[[i]])
         themap[[i]] <- themap[[i]] + c(0, cumsum(rep(amount, n - 1)))
    }
    return(themap)
}

#' Given a reference genetic map and physical positions of markers, return imputed genetic positions
#'
#' @param vcf.cross vcfR object for the biallelic variants segregating for a given cross
#' @param gmap a genetic map object
#' @param uchr a vector of chromosome names
#' @return a list of cM positions for given physical positions for markers on each chromosome
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


#' Given a vcfR object, get marker IDs sorted in intended chromosome order
#'
#' @param vcf.cross vcfR object for the biallelic variants segregating for a given cross
#' @param uchr vector of chromosome names in intended order
#' @return vector of marker IDs
getIDSorted=function(vcf.cross, uchr){
    gID=split(vcfR::getID(vcf.cross),vcfR::getCHROM(vcf.cross))
    gID=gID[uchr]
    return(stack(gID)$values)
}

#' Simulate genotypes for haploid progeny for a given cross given known biallelic variants and a genetic map
#'
#' @param vcf.cross vcfR object for the biallelic variants segregating for a given cross
#' @param gmap a genetic map object
#' @param nsegs number of desired haploid progeny
#' @param uchr a vector of chromosome names
#' @param ngenerations (default=2) number of generations of intercrossing, assuming only progeny from most recent generation mate etc
#' @return XSnpMatrix snpStats matrix object of 0/1 (relative to reference) coded genotypes of progeny 
#' @export
#' @import snpStats
simHaploidSegsFN=function(vcf.cross, gmap, nsegs, uchr=paste0('chr', as.roman(1:16)), ngenerations=2 ) {
 
    eg.GT=vcfR::extract.gt(vcf.cross)
    imputed.positions=jitterGmapVector(getGmapPositions(vcf.cross, gmap, uchr))

    #how many chromosomes
    n_chr=length(imputed.positions)
    #total gmap size
    L=round(sapply(gmap, function(x) max(x$map)))
    #setup for Meiosis package
    xoparam = Meiosis::create_xoparam(L,obligate_chiasma=T) 
    
    #parental haplotypes
    p1.ind=split(as.numeric(eg.GT[,1]), vcfR::getCHROM(vcf.cross))
    p1.ind=p1.ind[uchr]
    p2.ind=split(as.numeric(eg.GT[,2]), vcfR::getCHROM(vcf.cross))
    p2.ind=p2.ind[uchr]

    ind=list(p1.ind, p2.ind)
    #total number of meiosis to simulate
    psize=round(nsegs/2)
    #do sim
    p_geno = replicate(psize, Meiosis::cross_geno(father = ind, mother = ind, positions = imputed.positions, xoparam=xoparam), 
                       simplify=F)
 
    if(ngenerations>2){
        for(ngen in 3:ngenerations) {
            print(ngen)

            p_geno=replicate(psize, {
                    n1=sample.int(psize,1)
                    n2=sample.int(psize,1)
                    #paternal and maternal are interchangeable here but might as well just choose one for each sampling
                    nind=list(paternal=p_geno[[n1]]$paternal,
                              maternal=p_geno[[n2]]$maternal)
                    return(Meiosis::cross_geno(father = nind, mother = nind, positions = imputed.positions,  xoparam=xoparam))
                           },simplify=F)

        }
    } 
    
    #restructure and store in snpStats XSnpMatrix structure to shink it a bit 
    p_geno=lapply(1:n_chr, function(cc) {
                    new("XSnpMatrix",
                    cbind(sapply(p_geno,function(x) x$paternal[[cc]]),
                          sapply(p_geno,function(x) x$maternal[[cc]])), diploid=F)
                           })
    p_geno=do.call('rbind', p_geno)

    #agh, chr names 
    gID=getIDSorted(vcf.cross,uchr)
    rownames(p_geno)=gID
    
    return(p_geno)
}


#' Simulate phenotypes for haploid progeny for an XSnpMatrix with reference-based biallelic coding
#'
#' Simulate a number of QTL with equal effects (and random sign) with known summed total h^2 
#' @param sx XSnpMatrix
#' @param h2 total additive heritability
#' @param nQTL total number of QTL
#' @return y vector of simulated phenotype, scaled so total phenotypic variance is 1
#' @export
simPhenotypesForRefPop=function(sx, h2=.1, nQTL=1) {

    nmarker=nrow(sx)
    nsample=ncol(sx)
    #h2=.1
    #for mvnorm
    #nadditive = nsample
    #a.eff  = rnorm(nsample) # rep(0,nsample)
    a.eff=rep(0,nsample)
    nadditive=nQTL
    add.qtl.ind  = sort(sample(nmarker, nadditive))
    add.qtl.sign = sample(ifelse(runif(nadditive)>.5,1,-1),replace=T)
    for(i in 1:nadditive){ a.eff=a.eff+add.qtl.sign[i]*as.numeric(sx[add.qtl.ind[i],]) }
    a.eff=scale(a.eff)

    g=sqrt(h2)*a.eff #+ sqrt(H2-h2)*aa.eff
    y=g+rnorm(nsample,mean=0,sd=sqrt((1-h2)/h2*var(g)))
    y=scale(as.vector(y))

    attr(y, "add.qtl")=add.qtl.ind*add.qtl.sign

    return(y)
}

#' Simulate reference-based 0/1 counts for ref and alt alleles
#'
#' @param y vector of simulated phenotype, scaled so total phenotypic variance is 1
#' @param sx XSnpMatrix
#' @param sel.frac fraction of population selected
#' @param depth average depth per variant site (summed over both alleles)
#' @param lower.tail selecting lower tail (T/F, default=F)
#' @return data.frame of variant ID, expected allele frequence, ref counts and alt counts
generateCounts=function(y, sx, sel.frac, depth, lower.tail=F) {
    #$sel.frac=.10
    #get positive tail
    
    if(sel.frac==1) { 
    
    sel.indv.af=attr(geno.matrix, 'sel.indv.af')
    } else {
    
        if(lower.tail==F) {
            sel.indv=which(y> quantile(y,1-sel.frac))
        } else{
            sel.indv=which(y< quantile(y,sel.frac))
        }
        sel.indv.x=sx[,sel.indv]

        #lol R can't handle it, break it up
        splitit=cut(1:ncol(sel.indv.x),5, labels=F)
        saf=lapply(1:5, function(h){
                   idx=which(splitit==h)
                   return(Rfast::rowsums(sel.indv.x[,idx])) #/ncol(sel.indv.x)
                               })
        sel.indv.af=Rfast::rowsums(do.call('cbind', saf))/ncol(sel.indv.x)
    }

    #assuming 0/1 coding 
    #sel.indv.af=Rfast::rowsums(sel.indv.x)/ncol(sel.indv.x)
    #plot(sel.indv.af)

    r=rbinom(n=length(sel.indv.af),size=depth, prob=sel.indv.af)
    a=rbinom(n=length(sel.indv.af),size=depth, prob=1-sel.indv.af)

    countdf=data.frame(ID=rownames(sx), expected=sel.indv.af, ref=r, alt=a)
    return(countdf)
}

#' Simulate an xQTL experiment with reference-based 0/1 counts for ref and alt alleles
#'
#' @param y vector of simulated phenotype, scaled so total phenotypic variance is 1
#' @param geno.matrix XSnpMatrix
#' @param vcf.cross vcfR object for the biallelic variants segregating for a given cross
#' @param sel.low fraction of population selected for low tail
#' @param sel.high fraction of population selected for high tail
#' @param depth.low average depth per variant site for low tail (summed over both alleles)
#' @param depth.high average depth per variant site for low tail (summed over both alleles)
#' @param vcf.out filename of simulated vcf file output
#' @return expected allele frequencies for low and high tail 
#' @export
simXQTLExperiment=function(y, geno.matrix, vcf.cross, sel.low=0.1, sel.high=0.1,depth.low=50, depth.high=50, vcf.out=NULL){
    #get simulated high tail
    #sel.low=0.1
    #sel.high=0.1

    #simulated depth
    #depth.high=50
    #depth.low=50

    low.tail  = generateCounts(y,geno.matrix, sel.low, depth.low, lower.tail=T)
    high.tail = generateCounts(y,geno.matrix, sel.high, depth=depth.high)

    simulated.vcf=vcf.cross

    low.tail.GT=rep("0/1", nrow(vcf.cross))
    low.tail.AD=paste0(low.tail$ref,',', low.tail$alt)
    low.tail.DP=low.tail$ref+low.tail$alt
    low.tail.GQ=rep(99, nrow(vcf.cross))
    low.tail.PL=rep(paste(100,0,100, sep=','),nrow(vcf.cross))
    low.tail.sim=paste0(low.tail.GT, ':', low.tail.AD, ':', low.tail.DP, ':', low.tail.GQ, ':', low.tail.PL)
    names(low.tail.sim)=rownames(geno.matrix)

    high.tail.GT=rep("0/1", nrow(vcf.cross))
    high.tail.AD=paste0(high.tail$ref,',', high.tail$alt)
    high.tail.DP=high.tail$ref+high.tail$alt
    high.tail.GQ=rep(99, nrow(vcf.cross))
    high.tail.PL=rep(paste(100,0,100, sep=','),nrow(vcf.cross))
    high.tail.sim=paste0(high.tail.GT, ':', high.tail.AD, ':', high.tail.DP, ':', high.tail.GQ, ':', high.tail.PL)
    names(high.tail.sim)=rownames(geno.matrix)

    agh=vcfR::getID(simulated.vcf)

    simulated.vcf@gt[,2]=low.tail.sim[agh]
    simulated.vcf@gt[,3]=high.tail.sim[agh]
    colnames(simulated.vcf@gt)[2]='low.tail.sim'
    colnames(simulated.vcf@gt)[3]='high.tail.sim'

    if(!is.null(vcf.out)){
        vcfR::write.vcf(simulated.vcf, file=vcf.out)
        
        low.tail.expected=low.tail$expected
        names(low.tail.expected)=low.tail$ID
        
        high.tail.expected=high.tail$expected
        names(high.tail.expected)=high.tail$ID
        expected.af=list(low.tail.expected=low.tail.expected, high.tail.expected=high.tail.expected)
        return(expected.af)
    } else{
        return(list(low.tail=low.tail, high.tail=high.tail))
    }
}


#' Extract reference and alt counts for a given experiment from a VCF
#'
#' @param experiment.vcf vcfR object for the biallelic variants segregating for a given cross and allelic depths for reference and alt alleles
#' @param exp.name column name in vcf for experiment of interest
#' @param uchr vector of chromosome names in sorted desired order (optional)
#' @return data.frame of variant ID, ref count, and alt count 
#' @export
getBiallelicCounts=function(experiment.vcf, exp.name, uchr=paste0('chr', as.roman(1:16))) {
    uneffedvariants=getIDSorted(experiment.vcf,uchr)
    dfg=vcfR::extract.gt(experiment.vcf, 'AD')[,exp.name]
    dfg=dfg[uneffedvariants]
    df=data.frame(ID=names(dfg), 
                         ref=data.table::tstrsplit(dfg, ',',type.convert=T)[[1]],
                         alt=data.table::tstrsplit(dfg, ',',type.convert=T)[[2]])
    return(df)
}


#' Phase reference and alt counts given the ref/alt calls for one of the parents
#'
#' @param vcf.cross vcfR object for the biallelic variants segregating for a given cross 
#' @param p1.name name of parent to be designated parent 1 in vcf.cross (name must exist in vcf.cross) 
#' @param df data.frame output from getBiallelicCounts() should contain the columns: ID,ref,alt
#' @param uchr vector of chromosome names in sorted desired order (optional)
#' @return data.frame of variant ID, ref count, and alt count phased  
#' @export
phaseCounts=function(vcf.cross, p1.name, df, uchr=paste0('chr', as.roman(1:16))){
    gID=getIDSorted(vcf.cross, uchr)
    #indictor for whether parent 1 is reference 
    p1.ref=vcfR::extract.gt(vcf.cross)[,p1.name]=='0'
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

#' Simulate phenotypes for haploid prograny for an XSnpMatrix with reference-based biallelic coding coming from a BayesC genetic architecture
#' @param geno.matrix XSnpMatrix
#' @param h2 total additive heritability
#' @param target.size fraction of markers with effects sampled from N(0,h^2)
#' @return y vector of simulated phenotype, scaled so total phenotypic variance is 1
#' @export
simBayesC=function(geno.matrix,target.size=0.1, h2=0.4 ){
    #geno.matrix=matrix(as.numeric(geno.matrix@.Data), nrow(geno.matrix), ncol(geno.matrix))
    #target.size=.1
    nmarker=nrow(geno.matrix)
    nsample=ncol(geno.matrix)

    tindx=round(target.size*nrow(geno.matrix))
    causal.indices=sort.int(sample(1:nrow(geno.matrix), tindx))
    causal.betas=rnorm(tindx)
    #snorm=matrix(rnorm(nrow(geno.matrix)),1)
    a.eff=rep(0,ncol(geno.matrix))
    for(m in 1:length(causal.indices)) {
        if(m%%1000==0) {print(m)}
        a.eff=a.eff+as.numeric(geno.matrix[causal.indices[m],])*causal.betas[m] #snorm[m]
        #print(m)
       # g=as.matrix(geno.matrix@.Data)
    }
    a.eff=scale(a.eff)
    g=sqrt(h2)*a.eff #+ sqrt(H2-h2)*aa.eff
    y=g+rnorm(nsample,mean=0,sd=sqrt((1-h2)/h2*var(g)))
    y=scale(as.vector(y)) 

    attr(y, "causal.indices")=causal.indices
    attr(y, "causal.betas")=causal.betas
   
    return(y)

}


