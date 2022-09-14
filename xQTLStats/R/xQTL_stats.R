usethis::use_pipe(export =TRUE)

#' Convert pearson R to LOD score
#' 
#' @param n sample size 
#' @param r pearson r
#' @return LOD score
#' @export
RtoLOD=function(n, r) {(-n*log(1-r^2))/(2*log(10)) }

#' Convert LOD to p-value
#' 
#' Convert LOD to LRS, and then assume LRS is χ2 distributed with 1 d.f. 
#' @param x LOD score
#' @return p-value
#' @export
LODToPval =  function(x){  pchisq(x*(2*log(10)),df=1,lower.tail=FALSE)/2 }

#' Convert p-value to LOD score
#' 
#' Convert p-value to LOD score, inverse of LODToPval
#' @param x p-value
#' @return LOD score
#' @export
PvalToLOD =  function(x){ l= qchisq(x*2,df=1,lower.tail=FALSE)/(2*log(10)) 
                          l[is.na(l)]=0
                          return(l) }

#' Do loess regression of p1/(p1+p2) per block, in user defined blocks
#' 
#' @param p1 counts for parent 1 at each position
#' @param p2 counts for parent 2 at each position
#' @param pos positions for each marker 
#' @param bin.width collapse counts to bins every 500bp (default) or user specificied physical marker positions (optional)
#' @return imputed fitted allele frequency delta at each pos
#' @export
doBlockLoess=function(p1, p2, pos, bin.width=500) {
    loess.results=calcBlockLoess(p1,p2,pos, bin.width=bin.width, DEG=2)
    imputed.afd.delta=approxfun(loess.results$x0,loess.results$value)(pos)
    return(imputed.afd.delta) 
}

#SEE BRM https://github.com/huanglikun/BRM and
#https://academic.oup.com/bioinformatics/article/36/7/2150/5631910
# meta value of block
#https://academic.oup.com/bioinformatics/article/36/7/2150/5631910
def_block_pos <- function(chr, size){
	# block number
	n <- as.integer(2*chr/size)+2;
	if( n%%2 != 0 ) n <- n+1;
	n <- as.integer(n/2);
	# block index and the middle position of each block
	i <- c(1:n);
	pos <- (i-1)*size+floor(size/2); # middle position
	if(pos[n]>chr) pos[n]<-chr;
	return(pos);
}

cal_block_meta <- function(loc, val, chr, size, depth, MIN){
	# input: location vector, value vector
	# input: chr length, block size, location depth vector
	pos <- def_block_pos(chr, size);
	idx <- as.integer(0.5+loc/size)+1;
	#
	avg <- c();
	blockDepth <- c();
	for(i in 1:length(pos)){
		k  <- which(idx==i);
		no <- length(k);
		a <- NA;
		n <- 0;
		if (no > 0) {n <- sum(depth[k])};
		if( n >= MIN ) a <- sum(val[k])/n;
		avg <- c(avg, a);
	}
	return( list(pos=pos, avg=avg) );
}

calcBlockLoess=function(p1,p2,pos,bin.width=1000,min.depth=10,DEG=2,doSE=T) {
    opt.span <- function(model, criterion = c("aicc", "gcv"), span.range = c(0.05, 0.95)) {
	    as.crit <- function(x) {
	        span   <- x$pars$span;
	        traceL <- x$trace.hat;
	        sigma2 <- sum(x$residuals^2)/(x$n - 1);
	        aicc   <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL - 2);
	        gcv    <- x$n * sigma2/(x$n - traceL)^2;
	        result <- list(span = span, aicc = aicc, gcv = gcv);
	        return(result);
	    }
	    criterion <- match.arg(criterion);
	    fn <- function(span) {
	        mod <- update(model, span = span);
	        as.crit(mod)[[criterion]];
	    }
	    result <- optimize(fn, span.range);
	    return(list(span = result$minimum, criterion = result$objective));
	}
	# return

    #pmap binning 
    val=p1
    chr=max(pos,na.rm=T)+bin.width
    depth=p1+p2
    MIN=min.depth

    
    block <- cal_block_meta(pos, val,chr, bin.width, depth,MIN) 
    xx0 <- as.numeric(block$pos);
 #   print(xx0)
    yy <- as.numeric(block$avg);
    #print(x0)
    jdx <- which(!is.na(yy));
#    print(jdx)
 #   print(xx0[jdx])
    #print(jdx)
    fit0  <- loess(yy[jdx]~xx0[jdx], degree=DEG);

    span1 <- opt.span(fit0, criterion="aicc")$span;
    print(span1)
 #   print(span1)
   fit1  <- loess(yy[jdx]~xx0[jdx], degree=DEG, span=span1);
    plo <- predict(fit1, xx0, se=doSE);
    if(doSE){
        value <- plo$fit;
        return(list(x0=xx0, value=value, se=plo$se.fit))
    } else {
        value=plo
        return(list(x0=xx0, value=value))
        
    }
}

calcEffectiveNumberOfTests=function(Lg=4900,Lp=1.2e4, el=8.4, nchr=16){
    #Lg=4900
    #Lp=sum(sapply(gmaps[['A']], function(x) max(x$ppos)))/1e3
    #Lp=1.2e4
    rM=Lp/Lg # ~3.5kb per cM
    #el=8.4
    eff.tests=16+Lp/(rM*el)
    return(eff.tests)
}


#Calculate effect size given ratio of alleles
#K is selection threshold
#OR is observed P2/P1
# calculate x given OR and K
calcX.R = function (K, OR ) {
    T=qnorm(K, lower.tail = FALSE)
    fx=function(x,T,OR) { exp((pnorm(T, -x/2,sd=1, lower.tail=FALSE,log.p=T)
                              -pnorm(T,  x/2,sd=1, lower.tail=FALSE,log.p=T)))-OR }
    abs(uniroot ( 
        fx,
        c(-50,50), 
        T=T,
        OR=OR)$root)
}

# Leonid's clever approximation
calcX.LK = function(K,OR) {
    T=qnorm(K, lower.tail = FALSE)
    i=dnorm(T)/K
    (log(OR)/i)
}

#VE=x^2/4
#kvals = c(.1, .05, .01, .001, .0001, .00001, .000001)

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



#Lg=sum(sapply(gmaps[['A']], function(x) max(x$map)))
#eff.tests=calcEffectiveNumberOfTests()

#' Calculate allele frequency differences by block loess regression 
#' 
#' @param phasedCounts data frame with columns ID, p1, p2, which are marker ID, counts for parent 1 and counts for parent 2
#' @param experiment.name column name in vcf for experiment of interest
#' @param sample.size sample size of unselected population
#' @param sel.strength fraction of population selected
#' @param bin.width collapse counts to bins every 500bp (default) or user specificied physical marker positions (optional)
#' @param eff.length effective number of indepedent tests across the genome  
#' @param gmap genetic map object (optional) 
#' @param vcf.cross vcfR object containing the two parent variant calls (optional) 
#' @param uchr a vector of chromosome names (optional)
#' @return data.frame of smoothed results with SE
#' @export
calcAFD=function(phasedCounts, experiment.name='', 
                           sample.size=250000, 
                           sel.strength=.5,
                           bin.width=500,
                           eff.length=600,
                           gmap=NULL,
                           vcf.cross=NULL,
                           uchr = paste0('chr', as.roman(1:16))
                          ){
  
    chrom=data.table::tstrsplit(phasedCounts$ID, '_')[[1]]
    #uchr=paste0('chr', as.roman(1:16))
    chrom=factor(chrom, levels=uchr)
    physical.position=data.table::tstrsplit(phasedCounts$ID, '_', type.convert=T)[[2]]
    coreCols=data.frame(chrom=chrom, 
                        physical.position=physical.position)
    
    if(!is.null(gmap) & !is.null(vcf.cross)){
        genetic.position=stack(jitterGmapVector(getGmapPositions(vcf.cross, gmap, uchr)))$values
        coreCols$genetic.position= genetic.position
    }

    pC = dplyr::bind_cols(coreCols, phasedCounts) %>% 
            dplyr::mutate(nindv=round(sample.size*sel.strength)) %>%
            dplyr::mutate(ndepth=round(sum(p1+p2)/eff.length) ) %>%
            dplyr::mutate(n=min(nindv,ndepth)) %>%
            dplyr::group_by(chrom) %>%
            #dplyr::group_map(~dplyr::mutate(., 
            #   afd=doBlockLoess(.x$p1, .x$p2,.x$physical.position)), .keep=T) %>%
            #dplyr::bind_rows() %>% 
            dplyr::mutate(afd=doBlockLoess(p1, p2,physical.position, bin.width=bin.width)) %>% 
            dplyr::bind_rows() %>%
            dplyr::mutate(afd.se=sqrt((afd*(1-afd))/n)) 
       if(!is.null(gmap) & !is.null(vcf.cross)){
            pC=pC %>% dplyr::rename_with(.,~gsub('$', paste0('_', experiment.name), .x), .cols=c(-ID,-chrom, -physical.position,-genetic.position))
       } else{
            pC=pC %>% dplyr::rename_with(.,~gsub('$', paste0('_', experiment.name), .x), .cols=c(-ID,-chrom, -physical.position))

       }
    attr(pC, 'experiment.name')=experiment.name
    return(pC)
   
}


##' Use inverse variance weighting to combine information across replicates
##'
##' @param results data frame of results, should contain smoothed allele-frequency differences (afd) and SEs (afd.se)
##' @param meta_name specify name of meta analysis beta and se
##' @return combined.results meta analysis Betas and SEs added to results merged across the replicates 
##' @export 
calcMetaAFD=function(results, meta_name="meta"){
    combined.results = results %>% purrr::reduce(dplyr::left_join) %>% dplyr::ungroup()#, by=c('ID', 'chrom', 'physical.position') #%>% dplyr::ungroup()

    w=1/(combined.results %>% dplyr::select(starts_with('afd.se')))^2
    m=combined.results %>% dplyr::select(starts_with('afd_'))
    sw=apply(w, 1, sum)
    wm=apply(m*w, 1, sum)/sw
    wm.se=sqrt(1/sw)
    combined.results$wm=wm
    combined.results$wm.se=wm.se
    combined.results=combined.results %>% dplyr::rename(!!meta_name :=wm,!!paste0(meta_name,'.se') :=wm.se )
    return(combined.results)
}


##' Calculate allele frequency between tails, or tail and control. Convert to Z, p-value, and LOD.
##'
##' @param results list of two data frames of results to compare, should contain column names of betas and se with prefix L and R 
##' @param L prefix for column name containing beta (L) and SE L.se 
##' @param R prefix for column name containing beta (R) and SE R.se 
##' @return combined.results meta analysis Betas and SEs added to results merged across the replicates 
##' @export 
calcContrastStats=function(results, L='meta_high', R='meta_low'){
    
    combined.results = results %>% purrr::reduce(dplyr::left_join) %>% dplyr::ungroup() #, by=c('ID', 'chrom', 'physical.position') #%>% dplyr::ungroup()
    #contrast.label=paste0(L,'-',R)
    if(!is.null(R)){
        #calc contrasts without meta analysis 
        if(grepl('^_', L)){
                combined.results$contrast.beta=as.numeric(unlist(combined.results[,paste0('afd',L)]-combined.results[,paste0('afd',R)]))
                combined.results$contrast.beta.se = as.numeric(unlist(sqrt(combined.results[,paste0('afd.se',L)]^2+combined.results[,paste0('afd.se',R)]^2)))
        }else{
            combined.results$contrast.beta=as.numeric(unlist(combined.results[,L]-combined.results[,R]))
            combined.results$contrast.beta.se = as.numeric(unlist(sqrt(combined.results[,paste0(L, '.se')]^2+combined.results[,paste0(R,'.se')]^2)))
        }
       combined.results = combined.results %>% dplyr::mutate(z=contrast.beta/contrast.beta.se)
    } else {
        combined.results$z = as.numeric(unlist(combined.results[,L]/combined.results[,paste0(L, '.se')]))
    }
    combined.results=combined.results %>% 
        dplyr::mutate(p= 2*pnorm(abs(z), lower.tail=F) ) %>%
        dplyr::mutate(LOD=PvalToLOD(p)) %>%suppressWarnings()

    return(combined.results)
}

##' Plot individual experiment results
##'
##' @param results data.frame of experiment results 
##' @param suffix of experiment name 
##' @param simulatedQTL data frame of chromosome, and position of simulated QTL (optional)
##' @return ggplot object of experiment results  
##' @export 
##' @import ggplot2
plotIndividualExperiment=function(results, suffix,simulatedQTL=NULL){   
#suffix='high1'
    p1=rlang::sym(paste0('p1_',suffix))
    p2=rlang::sym(paste0('p2_',suffix))
    afd=rlang::sym(paste0('afd_', suffix))
    afd.se=rlang::sym(paste0('afd.se_', suffix))
    expected=rlang::sym(paste0('expected.phased_', suffix))
    if(paste0('expected.phased_', suffix)  %in% colnames(results) ) {Switch=T} else {Switch=F}


  ggplot(results, aes(x=physical.position,y={{p1}}/({{p1}}+{{p2}})))+    
    geom_point(size=0.3,alpha=0.6, color='gray21')+
    {if(!is.null(simulatedQTL))
        geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')}+
    facet_grid(~chrom, scales='free_x', space='free_x') +
    geom_ribbon(aes(ymin={{afd}}-1.96*{{afd.se}}, ymax={{afd}}+1.96*{{afd.se}}, fill='grey'), 
                alpha=0.7,linetype='dashed', color='grey')+
    geom_line(aes(x=physical.position, y={{afd}}),color='red', size=2, alpha=1)+
    {if(Switch) 
    geom_line(aes(x=physical.position, y={{expected}}),color='black', size=.5)}+
    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none')+ggtitle(suffix)

}

##' Plot contrast of experiment 1 - experiment 2
##'
##' @param results data.frame of experiment results 
##' @param suffix1 of experiment 1 name  
##' @param suffix2 of experiment 2 name 
##' @param simulatedQTL data frame of chromosome, and position of simulated QTL (optional)
##' @return ggplot object of experiment results  
##' @export 
##' @import ggplot2
plotContrast=function(results, suffix1, suffix2, simulatedQTL=NULL) {
    Lp1=rlang::sym(paste0('p1_',suffix1))
    Lp2=rlang::sym(paste0('p2_',suffix1))
    Lafd=rlang::sym(paste0('afd_', suffix1))
    Lafd.se=rlang::sym(paste0('afd.se_', suffix1))
    Lexpected=rlang::sym(paste0('expected.phased_', suffix1))
    if(paste0('expected.phased_', suffix1)  %in% colnames(results) ) {SwitchL=T} else {SwitchL=F}

    Rp1=rlang::sym(paste0('p1_',suffix2))
    Rp2=rlang::sym(paste0('p2_',suffix2))
    Rafd=rlang::sym(paste0('afd_', suffix2))
    Rafd.se=rlang::sym(paste0('afd.se_', suffix2))
    Rexpected=rlang::sym(paste0('expected.phased_', suffix2))
    if(paste0('expected.phased_', suffix2)  %in% colnames(results) ) {SwitchL=T} else {SwitchL=F}

    ggplot(results, aes(x=physical.position,y=({{Lp1}}/({{Lp1}}+{{Lp2}}))-({{Rp1}}/({{Rp1}}+{{Rp2}}))))+
    geom_point(size=0.3,alpha=0.6, color='gray21')+
    {if(!is.null(simulatedQTL))
        geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')}+
    facet_grid(~chrom, scales='free_x', space='free_x')+
    geom_ribbon(aes(ymin=({{Lafd}}-{{Rafd}})-1.96*sqrt( {{Lafd.se}}^2+{{Rafd.se}}^2),
                    ymax=({{Lafd}}-{{Rafd}})+1.96*sqrt( {{Lafd.se}}^2+{{Rafd.se}}^2), fill='grey'), 
                     alpha=0.7,linetype='dashed', color='grey')+
    geom_line(aes(x=physical.position, y={{Lafd}}-{{Rafd}}),color='red', size=2, alpha=1)+ 
    {if(SwitchL)
    geom_line(aes(x=physical.position, y={{Lexpected}}-{{Rexpected}}),color='black', size=.5)}+
    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none')+
    ggtitle(paste(suffix1, '-', suffix2))

}

##' Plot summary contrast and standard error, as well as -log10(p)
##'
##' @param results data.frame of experiment results 
##' @param simulatedQTL data frame of chromosome, and position of simulated QTL (optional)
##' @param effective.n.tests effective number of tests genomewide (integer)
##' @return ggplot object of experiment results  
##' @export 
##' @import ggplot2
plotSummary=function(results, simulatedQTL=NULL, effective.n.tests=600){   
    #suffix='high1'

    Z=ggplot(results, aes(x=physical.position,y=contrast.beta))+    
        #geom_line(size=0.3,alpha=0.6, color='gray21')+
        geom_line(color='black', size=.5)+ 
        {if(!is.null(simulatedQTL))
        geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')}+
        facet_grid(~chrom, scales='free_x', space='free_x') +
        geom_ribbon(aes(ymin=contrast.beta-1.96*contrast.beta.se, ymax=contrast.beta+1.96*contrast.beta.se, fill='grey'), 
                    alpha=0.7,linetype='dashed', color='grey')+
        #geom_line(aes(x=physical.position, y={{afd}}),color='red', size=2, alpha=1)+
        #{if(Switch) 
        #geom_line(aes(x=physical.position, y={{expected}}),color='black', size=.5)}+
        theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none')+ggtitle('contrast.beta')
     nlp=ggplot(results, aes(x=physical.position, y=-log10(p)))+geom_line(color='black', size=1.5)+ 
          {if(!is.null(simulatedQTL))
            geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')}+
            facet_grid(~chrom, scales='free_x', space='free_x')+
            geom_hline(aes(yintercept=-log10(.05/effective.n.tests)), color='red')+
            theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle('-log10(p)')

     ggpubr::ggarrange(Z,nlp, nrow=2) 
}


