
library(seqinr)
library(foreach)
library(doMC)
registerDoMC(cores=70)
#rf1=read.fasta('~/Desktop/isla.fasta')
rf1=read.fasta('/data/Databases/Yeast/orf_trans_all.fasta')
a.acids=a()[-1]

out.dir='/data/Databases/Yeast/orf_all_AA/'
dir.create(out.dir)
 
#Call provean
foreach (gene = names(rf1)) %dopar% {
    print(gene)
    aout=rf1[[gene]]
    protein.sequence = toupper(as.character(aout[-length(aout)]))
    query.file=paste0(out.dir, gene, '.fasta')
    write.fasta(protein.sequence, names=gsub('>', '', attr(aout, 'Annot')), file.out=query.file)

    dps=data.frame(AA.ref=protein.sequence, AA.pos=seq(1,length(protein.sequence)))
    dps.split=split(dps, dps$AA.pos)
    dps.df=do.call('rbind',lapply(dps.split, function(x)
    data.frame(x, AA.alt=a.acids) ) )

    dps.df$AA.ref=as.character(dps.df$AA.ref)
    dps.df$AA.alt=as.character(dps.df$AA.alt)
    dps.df=dps.df[dps.df$AA.ref!=dps.df$AA.alt,]

    dps.out=apply(dps.df, 1, function(x) paste(x[1], x[2], x[3], collapse='') )
    dps.out=gsub(' ', '', dps.out)

    last.aa=dps.split[[length(dps.split)]]
    
    dps.truncated.end=as.character(paste0(last.aa$AA.ref, last.aa$AA.pos))
    dps.truncated.end=paste(paste0(dps[,1], dps[,2]), dps.truncated.end, sep='_')
    dps.truncated.end=paste0(dps.truncated.end, 'del')

    dps.out=c(dps.out, dps.truncated.end)
    var.file = paste0(out.dir, gene, '.var')
    write.table(dps.out, file=var.file, col.names=F, row.names=F, sep='\t', quote=F)

    supporting_set=paste0(out.dir, gene, '_', 'set')
    output_file=paste0(out.dir, gene, '_', 'var.provean.txt')    

    system(
           paste(
                 '/home/jbloom/Local/bin/provean.sh -q' ,
                 query.file,
                 '-v',
                 var.file,
                '--num_threads 1 --save_supporting_set',
                 supporting_set, '>', output_file
                ) 
           )
}

# now read in provean output
in.dir='/data/Databases/Yeast/orf_all_AA/'
#gene='YAR042W'
in.files=list.files(in.dir, pattern='_var.provean.txt') #paste0(in.dir, gene, '_var.provean.txt')

provean_scores=list()
for(in.file in in.files[1:length(in.files)]) {
    pin = try( {read.delim(paste0(in.dir,in.file), skip=13, header=F, stringsAsFactors=F, sep='\t') })
    if(class(pin)=='try-error') {print(pin); next;} 
    # check YDR420W
    # YLR105C
    # YMR218C
    truncations=grep('del', pin[,1])
    AA.changes=pin[-truncations,]
    AA.endtruncations=pin[truncations,]
    
    AA.change.table=do.call('rbind', strsplit(AA.changes[,1], '\\d+'))
    AA.change.table=data.frame(pos=as.numeric(gsub('[A-Z]', '', AA.changes[,1])), 
                               ref=AA.change.table[,1], 
                               alt=AA.change.table[,2],
                               provean=AA.changes[,2], 
                               stringsAsFactors=F)
    if(sum(is.na(AA.change.table$pos))>0 ) {next;}

    AA.change.matrix=matrix(NA,max(AA.change.table$pos),20)

    A1.ind=match(as.character(AA.change.table[,2]), a.acids)
    A2.ind=match(as.character(AA.change.table[,3]), a.acids)
    inds=cbind(AA.change.table[,1], A2.ind)
    AA.change.matrix[inds]=AA.change.table$provean
    colnames(AA.change.matrix)=a.acids
      
    gene=strsplit(in.file, '_')[[1]][1]
    print(gene)
    fasta.in=read.fasta(paste0(in.dir, gene, '.fasta'))

    AA.seq=toupper(as.vector(fasta.in[[1]]))
    #    read.fasta( a.acids[A1.ind[match(unique(AA.change.table$pos), AA.change.table[,1])]]
    rownames(AA.change.matrix)=AA.seq
   
    provean_scores[[gene]]=AA.change.matrix 
}
saveRDS(provean_scores, file='/data/Databases/provean_scores.RDS')


# for visualizing provean scores per gene like 
#https://www.r-bloggers.com/adding-a-scale-to-an-image-plot/
image.scale <- function(z, zlim, col = heat.colors(12), ca=2, breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
     if(!missing(breaks)){
      if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
     }
     if(missing(breaks) & !missing(zlim)){
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
     }
     if(missing(breaks) & missing(zlim)){
      zlim <- range(z, na.rm=TRUE)
      zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
      zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
     }
     poly <- vector(mode="list", length(col))
     for(i in seq(poly)){
      poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
     }
     xaxt <- ifelse(horiz, "s", "n")
     yaxt <- ifelse(horiz, "n", "s")
     if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
     if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
     if(missing(xlim)) xlim=XLIM
     if(missing(ylim)) ylim=YLIM
     plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xlab='', ylab='', xaxs="i", cex.axis=ca, yaxs="i", ...)  
     for(i in seq(poly)){
      if(horiz){
       polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
      }
      if(!horiz){
       polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
      }
     }
}


# extra stuff with domain annotations ---------------------------------------------------------------------------------



domains2=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/pfam.tab2', header=F, sep='\t', stringsAsFactors=F)
#domains=domains[domains[,4]=='Pfam',]
dgsplit=split(domains2, domains2[,1]) 





#m=cbind(c(2:(total.oligos+1)),c((total.oligos+2):(total.oligos+total.oligos+1)))
library(fields)
m=matrix(1, 8, 30)
m[8,]=3
m[,30]=2
m[8,30]=4
layout(m)

#image.plot(AA.change.matrix[1:30,], axes=F)
gene=names(provean_scores)[700]

protein.length=nrow(provean_scores[[gene]])
par(mar=c(1,1,1,2), oma=c(1,1,4,1))
image(provean_scores[[gene]], col=(tim.colors(64)), axes=F)# , main=gene)
axis(2, at=seq(0,1,length.out=20), labels=colnames(provean_scores[[gene]]), cex=10, las=2)
axis(3, at=seq(0,1,length.out=protein.length), labels=rownames(provean_scores[[gene]]) , las=1 )
pretty.axis=c(pretty(1:nrow(provean_scores[[gene]]), 10, min.n=1), protein.length)
axis(1, at=pretty.axis/protein.length, labels=pretty.axis , las=1, cex.axis=1.5 )
box() 
#image.plot(AA.change.matrix,legend.only=TRUE, add=T)
image.scale(provean_scores[[gene]], col=tim.colors(64), horiz=FALSE, ca=1.5)
if(!is.null(dgsplit[[gene]])) {
    par(xaxs='i')
    plot(0,0, xlim=c(1, protein.length), ylim=c(0,1), type='n' ,yaxt='n')
    rect(dgsplit[[gene]][,6], 0, dgsplit[[gene]][,7], 1, col='#ff000033', lwd=2, density=20)
    text(mean(c(dgsplit[[gene]][,6], dgsplit[[gene]][,7])),.5, dgsplit[[gene]][,5], cex=2)
}
title(gene, outer=T)


library(Rpdb)
#http://foldxsuite.crg.eu/command/PositionScan

