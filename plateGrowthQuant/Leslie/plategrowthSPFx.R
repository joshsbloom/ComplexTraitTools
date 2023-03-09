

#for debugging 
#options(browser='chromium-browser')

# a couple of global variables -------------------------------------------
plate.types  <- list('96'=c(8,12), '384'=c(16,24), '1536'=c(32,48) )
dimensions = c(5184,3456)

#NUNC
corners=list()
#corners[['nunc']] = data.frame(X=c(525,4618,525,4624),
#                               Y=c(388,400,3060,3064))

corners[['nunc']] = data.frame(X=c(562,4674,554,4666),
                               Y=c(388,392,3063,3064))

# new 05/01/22
#VWR plates
corners[['vwr']]= data.frame( X=c(576, 4660, 572, 4660),
                     Y=c(392, 396,  3060, 3068))
#-----------------------------------------------------------------------



# Functions! 

# Part 1 functions start ##################################################################################
# make grid for expected colony locations
makeGrid = function(plate.type, corners) {
    x.c=sort(corners$X)
    y.c=sort(corners$Y)
    av=list( xmin=round(mean(x.c[c(1,2)])), xmax=round(mean(x.c[c(3,4)])),
                     ymin=round(mean(y.c[c(1,2)])), ymax=round(mean(y.c[c(3,4)])))
    gridme = expand.grid( 
                    round(seq(av$xmin, av$xmax, length.out=max(plate.type)) ), 
                    round(seq(av$ymin, av$ymax, length.out=min(plate.type))))
    colnames(gridme)=c('X','Y') 
    return(gridme)
}     

# segment plate into features and non-features using kmeans clustering
# hard-coded a lot of the dimensions here given the consistent behavior of the S&P pinning robot
segmentPlate=function(img, spot.seeds){
    print('Segmenting plate')
    img2=img
    #top
    img2[,1:307]=0
    #bottom
    img2[,3149:3456]=0
    #left
    img2[1:331,]=0
    #right
    img2[4834:5184,]=0 

    lrect=img2[331:528,2950:3149]
    lrect[upper.tri(lrect)]=0
    img2[331:528,2950:3149]=lrect
    rrect=img2[4650:4834,2950:3149]
    tr=flip(rrect)
    tr[lower.tri(tr)]=0
    rrect=flip(tr)
    img2[4650:4834,2950:3149]=rrect
    
    img=img2
    tplate=img
    
    intensity.vector=as.vector(img2)
    informative.area=which(intensity.vector!=0)
    #km = kmeans(as.vector(img), centers=c(.3906,.7812))
    km = kmeans( intensity.vector[informative.area], centers=c(.39,.7812))
    #km = kmeans(as.vector(img), centers=c(0,.7812))

    kmin = which.min(km$centers)
    kmax = which.max(km$centers)
    if(kmin==2) {km$cluster = abs(km$cluster-3)}
    matrix.vector=rep(0, length(intensity.vector))
    matrix.vector[informative.area]=km$cluster-1
    tplate@.Data=matrix(matrix.vector,nrow(tplate),ncol(tplate))

   # EBImage::display(matrix(km$cluster-1, 5184,3456))

    seeds = matrix(0, nrow(img), ncol(img))
    seeds[is.numeric(seeds)]=0
    ycirc    <- makeBrush(49, 'disc', step=TRUE)

    for(i in 1:nrow(spot.seeds) ){
            x = spot.seeds[i,'X']
            y = spot.seeds[i, 'Y']
            seeds[(x-24):(x+24), (y-24):(y+24) ] =ycirc*i
    }
    
    # find voroni region for each seed
    plate.mask = propagate(img, seeds, mask=tplate)
    return(plate.mask)
}

# calculate colony features and attach names
getColonyFeatures=function(plate.mask, img, layout.file){
    print('Getting features')
    strain.layout=read.delim(layout.file, sep=',', header=F, stringsAsFactors=F)
    mfeatures=computeFeatures.moment(plate.mask, img )
    sfeatures=computeFeatures.shape(plate.mask)
    nfeatures=data.frame(mfeatures, sfeatures)
    rownames(nfeatures)=paste(seq(1,384), as.vector(t(strain.layout)), sep=':')
    #as.vector(t(strain.layout))
    return(nfeatures)
}

# spline fit for position effects
# calculate colony features and attach names
localPfit = function(vecin, x,y) {
   m = locfit(vecin~lp(x,y, scale=FALSE),lfproc=locfit.robust)
   predx = predict(m, data.frame(vecin, x,y))
   return(vecin-predx) }

# do everything
# assuming a data frame key with rows i
# for each row i
processImages = function(i, key, image.path, layout.path, outdir, plate.type, corners, flip=F, save.labeled=F){
    print(i)
    image.file   = paste(image.path, key$Filename[i], sep='')
    layout.file  = paste(layout.path, key$StrainLayout[i], sep='')
    output.file =  paste(outdir, i, '.jpeg', sep='')
    if(flip) {
        raw.img   = rotate( readImage(image.file), 180)
    } else  {
        raw.img = readImage(image.file)
    }        
    img          = channel(raw.img, 'grey')
    spot.seeds   = makeGrid(plate.type, corners)
    plate.mask   = segmentPlate(img,spot.seeds)
    features     = getColonyFeatures(plate.mask, img, layout.file)
    #display(paintObjects(plate.mask, raw.img ), title=key$Filename[i] )
    pO=paintObjects(plate.mask, raw.img )
    
    if(save.labeled) {
        jpeg(output.file, width=dim(pO)[1],height=dim(pO)[2])
        display(pO, method='raster')
        text(x=features$m.cx, y=features$m.cy, 
             label=rownames(features), #data.table::tstrsplit(rownames(features),';')[[3]],
            , adj=c(0,1), col='blue', cex = 2 )
        dev.off()
    } else {

        writeImage(pO, file=output.file)
    }
    return(features)
}
#  text files ... 
writeTXTfiles = function(outD, Results.norm) {
    dir.create(outDir)
    for(n in names(Results.norm)){
        n2 =n
        n2=gsub('/', '__', n2)
        n2=gsub(':', '_', n2)
        write.table(Results.norm[[n]], file=paste(outDir,n2,'.txt',sep=''), 
                sep='\t', row.names=T,col.names=NA,quote=F) }
}
################### Part 1 functions end #####################################################################


# Part 2 functions start ####################################################################################
# reconstruct key from 
reconstructKey=function(Results) {
    ndf=c("StrainLayout" ,    "PermutationGroup", "Filename"   ,      "Condition"     ,   "Concentration"  ,  "Control")
    df=data.frame(do.call('rbind', strsplit(names(Results), '::'))) 
    names(df)=ndf
    return(df)
}

getResults.df=function(Results, key, ind, plate.type=384) {
    R_list=Results.processed[ind]
    key_sub=key[ind,]
    #key_sub$Concentration[key_sub$Concentration=='0']=''
    plate.factor=as.factor(rep(key_sub$Filename, each=plate.type))
    cond_conc.factor=rep(paste(key_sub$Condition, key_sub$Concentration, sep=';'), each=plate.type)
    R_list.df=do.call('rbind', R_list)
    strain.names = as.vector(unlist(lapply(R_list, function(x) rownames(x) )))
    strain.names=do.call('rbind', strsplit(strain.names, ':') )[,2] 
    R_list.df=data.frame(strain=strain.names, R_list.df, plate=plate.factor, condition=cond_conc.factor, stringsAsFactors=F)
    return(R_list.df)
}

filterEdges=function(Results, plate.type=c(16,24)){
    rows.n=plate.type[1]
    cols.n=plate.type[2]
    col.row = expand.grid(seq(1,cols.n), seq(1,rows.n))
    edges = list(    
        top.row = col.row[,2]==1,
       bottom.row = col.row[,2]==rows.n,
       left.col = col.row[,1]==1,
       right.col =col.row[,1]==cols.n)                                                      

    for( i in 1:length(Results) ){
        pa = Results[[i]]
        pagn = pa[,'s.radius.mean']
        e.ps= sapply(edges, function(x) {
                     tryCatch( {t.test(pagn[x], pagn)$p.value}, error=function(e) {return(1) } )  
         } )
     for(j in 1:4){   if(e.ps[j]<.05) {Results[[i]][edges[[j]], -c(1,2)]=NA;   } } }
    return(Results)
}

# within plate normalization --------------------------------------------------------------------------------------------
# this would be a good spot to adjust for plate effects, if desired
normalizePlate=function(Results){
   Results.processed=list()
   for( i in 1:length(Results) ){
        print(i)
        pa = Results[[i]]
        normalized =  apply(pa[,c('s.radius.mean', 's.area')], 2, function(x){ 
                    tryCatch( { return(localPfit(x, pa[,'m.cx'], pa[,'m.cy'])+mean(x,na.rm=T))}, 
                       error = function(e){return(x)} ) })
        colnames(normalized)=paste(colnames(normalized),'norm', sep='.')
        r = cbind(pa, normalized)
        Results.processed[[i]]=r
    }
    names(Results.processed)=names(Results)
    return(Results.processed)
}   
#-------------------------------------------------------------------------------------------------------------------------











#################### Part 2 functions end ######################################################################


# Part 3 Mapping functions ######################################################################################
processPhenos_for_MM=function(s, all.strain.names){
    #s=(pheno_raw[[phenotype]])
    # fix this -----------------------------------
    s=s[match(all.strain.names, names(s))]
    srle=sapply(s, length)
    sr=list()
    sr$lengths=as.vector(srle)
    sr$values=names(srle)
    attr(sr, 'class')='rle'
    ny=as.vector(unlist(s))
    names(ny)=inverse.rle(sr)
    #-----------------------------------------------
    y=ny[!is.na(ny)]
    strain.names=(names(y))
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)
    strain.ind  = seq_along(strain.names)
    strain.cnt  = length(unique.sn)
   
    #for constructing Strain Variance component
    Strain      = diag(strain.cnt)
    Z=Matrix(0, length(y), strain.cnt,sparse=T);   Z[cbind(strain.ind, n.to.m)]=1 
    strains.with.phenos=match(unique.sn, all.strain.names)
    n.strains=length(strains.with.phenos)
    return(list(y=y, Z=Z, Strain=Strain, strains.with.phenos=strains.with.phenos, n.strains=n.strains))
}  


# MAPPING
extractScaledPhenotype=function(impcross,scaleVar=FALSE){apply(impcross$pheno, 2, scale, scale=scaleVar)}
extractGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$data }))*2)-3 }
extractGenotype.argmax=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$argmax }))*2)-3 }
countStrainsPerTrait = function(pheno) {apply(pheno, 2, function(x){sum(!is.na(x))})}

get.LOD.by.COR = function(n.pheno, pheno, gdata,betas=FALSE, sdx=1) {
   # Lynch and Walsh p. 454
   r=cor(pheno, gdata, use='pairwise.complete.obs')
   LOD=(-n.pheno*log(1-r^2))/(2*log(10))
   if(betas==TRUE) {
       # convert pearson R to B from y=u+Bx+e model
       beta=r * apply(cbind(pheno),2, sd,na.rm=T)/sdx
       return(list(r=r, LOD=LOD, beta=beta))
   }
   else {   return( list(r=r, LOD=LOD) ) } 
}

calc.pval.LRtest=function(null,full) { pchisq(-2*(null-full),1, lower.tail=FALSE) }
calc.BLUPS= function(G,Z,Vinv,y,X,B ){    G%*%t(Z)%*%Vinv%*%(y- X%*%B)     }

extractVarCompResults = function(r) {list(sigma=r$Var, sigma.cov=r$invI, W=r$W,
                                          Bhat=as.vector(r$Bhat), llik=as.vector(r$llik)) }

# can generalize for more than two groups if desired
structured.perm=function(y, bck.split){
    yperm=y
    yperm[bck.split]=sample(y[bck.split])
    yperm[!bck.split]=sample(y[!bck.split])
    return(yperm)
}

# calculate averages for a named vector with replicates
avg.nvec=function(y, un=NULL){
    if(is.null(un)) {un=names(y) }
    y.avg.l=(by(y, un, mean, na.rm=T))
    nvec=names(y.avg.l)
    y.avg.l=as.vector(y.avg.l)
    names(y.avg.l)=nvec
    return(y.avg.l)
}
