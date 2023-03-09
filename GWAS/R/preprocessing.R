 
getGWASvariants=function(path.to.gwas.data){
    setwd(path.to.gwas.data)
    x=BEDMatrix('1011GWAS_matrix.bed')
    g=as.matrix(x)
    mnames=read.delim('1011GWAS_matrix.bim', header=F, sep='\t', stringsAsFactors=F)
    mnames=mnames[,c(1,2,4)]
    colnames(mnames)=c('chrom', 'mname', 'pos')
    return(list(g=g,mnames=mnames))
}


syncGenoPheno=function(pheno, GWAS.variants, strain.col='Strain_Name') {
    #remove blank and missing strains -----------------
    g.subset=GWAS.variants$g

    #match to phenotypes
    mind=pmatch(pheno$Strain_Name, rownames(g.subset) )
    pheno=pheno[-which(is.na(mind)),]
    mind=pmatch(pheno$Strain_Name, rownames(g.subset) )

    g.subset=g.subset[mind,]

    GWAS.variants$g=g.subset
    #morejunk=which(rownames(pheno)=='NA')
    #-------------------------------------------------
    return(list(pheno=pheno, GWAS.variants=GWAS.variants))

}

filterGWASvariants=function(GWAS.variants, af.cutoff=.05, na.frq=.1) {
   
    # calc af 
    g.s=GWAS.variants$g
    af=colSums(g.s, na.rm=T)/(2*apply(g.s,2,function(x) sum(!is.na(x)))) 
    af.cutoff=.05

    # calc missingness 
    notna.cnt=apply(g.s, 2, function(x) sum(!is.na(x)))
    notna.frq=notna.cnt/nrow(g.s)

    notna.frq.cut=1-na.frq

    #gvar=apply(g.s,2,var, na.rm=T)
    bm=which(af<af.cutoff | notna.frq<notna.frq.cut)

    #final marker set
    g.s=g.s[,-bm]
    #corresponding marker info
    mnames=GWAS.variants$mnames[-bm,]

    #filter variants that don't vary
    vx=apply(g.s,2,function(x) var(x[!is.na(x)]))
    bm=which(is.na(vx))

    if(length(bm)>0) {
        g.s=g.s[,-bm]
        mnames=mnames[-bm,]
    }
    
    return(list(g=g.s, mnames=mnames))    
}


calcA=function(GWAS.variants) {
    g.subset=GWAS.variants$g
    Amat=A.mat(g.subset-1, return.imputed=T)
    A=Amat$A
    g.subset.imp=Amat$imputed

    g.s=scale(g.subset.imp)

#    mnames=mnames[match(data.table::tstrsplit(colnames(g.s), '_')[[1]], mnames[,2]),]
#    bm=which(is.na(mnames[,1]))
#    g.subset.imp=g.subset.imp[,-bm]
#    g.s=g.s[,-bm]
#    mnames=mnames[-bm,]


    #calculate pariwise correlation between strains
    cor.mat=tcrossprod(g.s)
    cor.mat=cor.mat/(ncol(g.s))
    #---------------------------------------------------------

    GWAS.variants$A=cor.mat

    return(GWAS.variants)
}

