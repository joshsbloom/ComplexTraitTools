library(BEDMatrix)
library(rrBLUP)
library(NAM)
library(milorGWAS)
library(data.table)
library(ggplot2)
library(mvnpermute)
#library(regress)


#get the effective number of markers given the Li and Ji procedure
#https://neurogenetics.qimrberghofer.edu.au/SNPSpD/Li2005.pdf
getMeff_Li_and_Ji=function(cor.mat) {
    evals = eigen(cor.mat,symmetric=T)$values
    M = length(evals)
    L = M-1
    # Equation 5 from Li 
    intevals=ifelse(evals>=1, 1, 0)
    # modification for negative eigenvalues JB
    nonintevals=c(evals-floor(evals)) #[evals>0]
    Meff.li=sum(intevals+nonintevals)
    print(Meff.li)
    return(Meff.li)
}
#using genomewide LD matrix, pain to calc, don't recommend rerunning
Meff=42255
Meff.thresh=log10(.05/42255)
#5.927

# pre-processing genotypes and phenotypes -------------------------------------

    #get phenotypes 
    pheno=readRDS('/data/yeast/Colonies/SNR1p5Metrics.RDS')
    #remove NAs from phenotype for now 
    pheno=pheno[-c(47,485),]

    #load Joseph's prefiltered GWAS genotype data 
    #https://www.nature.com/articles/s41586-018-0030-5
    path.to.gwas.data='/data/yeast/1002genomes/1011GWASMatrix/'

    setwd(path.to.gwas.data)

    #read marker data (BED/BIM/BAM plink formatted)
    x=BEDMatrix('1011GWAS_matrix.bed')
    g=as.matrix(x)
    mnames=read.delim('1011GWAS_matrix.bim', header=F, sep='\t', stringsAsFactors=F)


    #remove blank and missing strains -----------------
    g.subset=g

    #match to phenotypes
    mind=pmatch(pheno$Row.names, rownames(g.subset) )
    pheno=pheno[-which(is.na(mind)),]
    mind=pmatch(pheno$Row.names, rownames(g.subset) )

    g.subset=g.subset[mind,]
    #morejunk=which(rownames(pheno)=='NA')
    #-------------------------------------------------

    #calculate relatedness matrix ----------------------------
    Amat=A.mat(g.subset-1, return.imputed=T)
    A=Amat$A
    g.subset.imp=Amat$imputed

    g.s=scale(g.subset.imp)

    mnames=mnames[match(data.table::tstrsplit(colnames(g.s), '_')[[1]], mnames[,2]),]
    bm=which(is.na(mnames[,1]))
    g.subset.imp=g.subset.imp[,-bm]
    g.s=g.s[,-bm]
    mnames=mnames[-bm,]


    #calculate pariwise correlation between strains
    cor.mat=tcrossprod(g.s)
    cor.mat=cor.mat/(ncol(g.s))
    #---------------------------------------------------------


    g.s=g.subset[,colnames(g.subset.imp)]

    #totalNA=sum(is.na(g.s))
    af=colSums(g.s)/(2*nrow(g.s))

    af=colSums(g.s, na.rm=T)/(2*apply(g.s,2,function(x) sum(!is.na(x)))) 
    af.cutoff=.05


    notna.cnt=apply(g.s, 2, function(x) sum(!is.na(x)))
    notna.frq=notna.cnt/nrow(g.s)
    notna.frq.cut=.8

    #gvar=apply(g.s,2,var, na.rm=T)
    bm=which(af<af.cutoff | notna.frq<notna.frq.cut)

    #final marker set
    g.s=g.s[,-bm]
    #corresponding marker info
    mnames=mnames[-bm,]

    #replace NAs with 0s
    g.s0=g.s
    g.s0[is.na(g.s0)]=0

#----------------------------------------------------------------------------------------

# do the logistic regression gwas
gwas.binary=list()
for(p in colnames(pheno)[which(grepl('b.p', colnames(pheno)))] ) {
        gwas.binary[[p]]= milorGWAS::association.test.logistic(x=as.bed.matrix(g.s0),
                                                  Y=pheno[,p], 
                                                  X=matrix(1,nrow(pheno),1),K=cor.mat, algorithm='offset')

}

#process results
gb=lapply(gwas.binary, function(x) {
           x$chr=mnames[,1]
           x$pos=mnames[,4]
           names(x)[1]='chrom'
           names(x)[3]='marker'
           x$nlogp=-log10(x$p)
           x$chrom=factor(x$chrom, levels=1:17)
           return(x)
                                 }) 
gbc=data.table::rbindlist(gb, idcol='trait')
gbc$sig=gbc$nlogp>Meff.thresh
gbc=gbc[gbc$trait!='b.ph7' | gbc$trait!='b.ph6',]
gbc=gbc[!is.na(gbc$nlogp),]
ggplot(gbc, aes(x=pos,y=nlogp,colour=sig))+scale_colour_manual(values=c('black', 'red'))+
       geom_point()+facet_grid(trait~chrom,scales="free", space='free')+
        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
        theme_classic()+
        theme(axis.text.x.bottom = element_text(angle = 45,hjust=1))+
        theme(panel.spacing=unit(0, "lines"))


saveRDS(gbc, file= '/data/yeast/Colonies/GB_pH_080922_logistic_binary_traits_gwas.RDS')



#see this for math derivations for factored matrices speedup 
#https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.1681/MediaObjects/41592_2011_BFnmeth1681_MOESM290_ESM.pdf
#novelty here is for the multitrait stuff
#and combination with mvnpermute, for thresholds robust to model violations
#https://github.com/markabney/MVNpermute
#to do, have delta be a vector 

otherphenos=read_csv('/data/yeast/Colonies/P2018_ST1.csv')
up=unlist(otherphenos[,13])
names(up)=paste0(otherphenos$"Standardized name", '_', otherphenos$"Standardized name")
ups=up[match(rownames(cor.mat), names(up))]

REML=TRUE

nperm=100   

K=cor.mat

USUt=svd(K)
U=USUt$u
Ut=t(USUt$u)

flmm=list()
for(p in colnames(pheno)[2:15]) {
    print(p)
    
    ytemp=pheno[,p]
    names(ytemp)=pheno[,1]
    #ytemp=ups
    

    y=scale(ytemp)
    g=g.s0

    nmm=mixed.solve(y,X=matrix(1,length(y)),K=K, method='ML')
  

    # nmmREML=mixed.solve(y,X=matrix(1,length(y)),K=K, method='REML')

    delta=nmm$Ve/nmm$Vu

    #function(y,g, delta, K) 
    #delta=sigmaE/sigmaG
    #eig=eigen(cor.mat)

    Sdelt=USUt$d+delta

    n=length(y)
    Uty=Ut%*%cbind(y, mvnpermute(y,matrix(1,length(y)),nmm$Vu*K+nmm$Ve*diag(length(y)),nperm))
    sds=solve(diag(Sdelt))


    llconst0=log(2*pi*nmm$Vu)
    llconst1=log(det(diag(Sdelt)))
    if(!is.finite(llconst1)) { llconst1=0 }
    llconst=n*llconst0+llconst1

    llconst2=(1/nmm$Vu)

    XBase=matrix(1,n)
    d=ncol(XBase)+1
    XComb=cbind(1,g)

    UtXall=Ut%*%cbind(1,g)

    #null model
    UtX=UtXall[,1]
    XX=XComb[,c(1)]
    tUS=t(UtX/Sdelt)
    B=solve(tUS %*% UtX) %*% (tUS %*%Uty)
    uubSq=(Uty-(UtX)%*%B)^2
    redll=-.5*(llconst+llconst2*colSums(uubSq/Sdelt))
    remlredll=redll+.5*(d*llconst0+log(det(t(XX)%*%XX))-log(det(tUS%*%UtX)))

    Bmat=matrix(NA,ncol(g),(nperm+1))
    ll=matrix(NA,ncol(g),(nperm+1))
    rll=ll
    #i=47718
    #system.time(
    for(i in 2:(ncol(g)+1) ) {
        if(i%%1000==0) {print(i) }
        UtX=UtXall[,c(1,i)] #cbind(1,UtXall[,i])

        tUS=t(UtX/Sdelt)
        #solve for betas
        B=solve(tUS %*% UtX) %*% (tUS %*%Uty)
        uubSq=(Uty-(UtX)%*%B)^2
        #genet  var exp vG=colSums(uubSq/Sdelt)/nrow(K)
        #calculate log likelihood
        #colSums
        #ll[i-1]=-.5*(llconst+llconst2*sum(uubSq/Sdelt))

        Bmat[i-1,]=B[2,]
        ll[i-1,]=-.5*(llconst+llconst2*colSums(uubSq/Sdelt))
        if(REML==TRUE){
            XX=XComb[,c(1,i)]
            rll[i-1,]=ll[i-1,]+.5*(d*llconst0+log(det(t(XX)%*%XX))-log(det(tUS%*%UtX)))
                #genet  var exp vG=colSums(uubSq/Sdelt)/nrow(K)-

        }

    }
    dff=data.frame(marker=colnames(g),
               Beta=Bmat[,1],
               ML.nlp=-log10(pchisq(-2*(redll[1]-ll[,1]),df=1, lower.tail=F)),
               REML.nlp=-log10(pchisq(-2*(remlredll[1]-rll[,1]),df=1, lower.tail=F))
               )
    attr(dff, "ML.perm")=apply(-log10(pchisq(-2*(redll[-1]-ll[,-1]),df=1, lower.tail=F)),2,max, na.rm=T)
    attr(dff,"REML.perm")=apply(-log10(pchisq(-2*(remlredll[-1]-rll[,-1]),df=1, lower.tail=F)),2,max, na.rm=T)
    flmm[[p]]=dff
}


flmm2=lapply(flmm, function(x) {
       x$chrom=mnames[,1]
       x$pos=mnames[,4]       
       x$MLsig=x$ML.nlp>quantile(attr(x,  "ML.perm"), .95)
       x$REMLsig=x$REML.nlp>quantile(attr(x,  "REML.perm"), .95)
       return(x)
               })


gr.combined=data.table::rbindlist(flmm2, idcol='trait')
gr.combined$chrom=factor(gr.combined$chrom,levels=1:17) 


#gr.combined$trait=factor(gr.combined$trait, 
#                                 levels=names(table(gr.combined$trait))[c(8,11:14,9,10,15:18,1:6,7,19,20)]
#                                 )
ggplot(gr.combined, aes(x=pos,y=ML.nlp,colour=MLsig))+scale_colour_manual(values=c('black', 'red'))+
       geom_point()+facet_grid(trait~chrom,scales="free", space='free')+
        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
        theme_classic()+
        theme(axis.text.x.bottom = element_text(angle = 45,hjust=1))+
        theme(panel.spacing=unit(0, "lines"))
ggplot(gr.combined, aes(x=pos,y=REML.nlp,colour=REMLsig))+scale_colour_manual(values=c('black', 'red'))+
       geom_point()+facet_grid(trait~chrom,scales="free", space='free')+
        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
        theme_classic()+
        theme(axis.text.x.bottom = element_text(angle = 45,hjust=1))+
        theme(panel.spacing=unit(0, "lines"))

saveRDS(gr.combined, file =  '/data/yeast/Colonies/GB_pH_080922_continuous_traits_perm_gwas.RDS')
gr.combined=readRDS( '/data/yeast/Colonies/GB_pH_080922_continuous_traits_perm_gwas.RDS')

f0=regress(pheno[,'ph6']~1, ~cor.mat, verbose=T)
f1=regress(pheno[,'ph6']~g.s0[,'rs27626_1'], ~cor.mat, verbose=T)
pchisq(-2*(f0$llik-f1$llik),df=1, lower.tail=F)

plot( pheno['ph6']-(BLUP(f0))




#  nmm1=mixed.solve(y,X=cbind(matrix(1,length(y)),g.s[,32413]),K=K, method='ML')
# 
#   
#    mnames[32413,]
#    g2=g.s[,32413]
#    g2[g2==1]=NA
#   # nmm20=mixed.solve(ytemp[!is.na(g2)],X=cbind(matrix(1,length(y)))[!is.na(g2),],K=K[which(!is.na(g2)),which(!is.na(g2))], method='ML')
#   # nmm2=mixed.solve(ytemp[!is.na(g2)],X=cbind(matrix(1,length(y)),g2)[!is.na(g2),],K=K[which(!is.na(g2)),which(!is.na(g2))], method='ML')
#
#    r0=regress(ytemp[!is.na(g2)]~1, ~K[which(!is.na(g2)),which(!is.na(g2))])
#    r1=regress(ytemp[!is.na(g2)]~g2[!is.na(g2)], ~K[which(!is.na(g2)),which(!is.na(g2))])
#
#    pchisq(-2*(r0$llik-r1$llik),df=1, lower.tail=F)
#
#
#
#    pchisq(-2*(nmm20$LL-nmm2$LL),df=1, lower.tail=F)

