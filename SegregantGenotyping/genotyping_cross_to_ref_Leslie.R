library(parallel)
devtools::source_url("https://raw.githubusercontent.com/joshsbloom/yeast-16-parents/master/genotyping/code/segregants_hmm_fx.R")

# filtering parameters -------------------------
parent.min.read.support=3
parent.gq.cutoff=98

segregant.quality.cutoff=100
segregant.gq.cutoff=30

tot.segs=96*2

segregant.missing.by.marker.cutoff=tot.segs*.9
segregant.missing.marker.fraction = .98
segregant.afH= .98
segregant.afL= .02
possible_het_ratio= 0.1
structual_variant_padding= 200
error.prob1=0.1
error.prob2=0.01 #075
max_XO_after_preimp=15
# -----------------------------------------------

cross.list=list()
parents.list=list()

# unusual input structure here
# in this case input here is a vcf containing pre-filtered biallelic heterozygote sites in a BY/XXX diploid 
# in the future we will want two columns in the vcf for the two haploid parents , could include diploid column as a reference
parent.vcf.file='/data/yeast/Leslie/het_sites_ly11.vcf.bg.gz'

#segregant data, long data frame with output of gatk ASEReadCounter, we'll want a module upstream here for running that and restructuring the output as a maatrix
seg.data='/data/yeast/Leslie/leslie_segregants.RDS'


# Read in parental data and reformat ----------------------------------------------------------------------------
#parents=buildParentsDF(parent.vcf.file, possible_het_ratio, parent.min.read.support, parent.gq.cutoff, blocked.regions)
pvcf=read.vcfR(parent.vcf.file)

#p.gt=extract.gt(pvcf, element='GT', as.numeric=T)
#p.dp=extract.gt(pvcf, element='DP', as.numeric=T)
#p.gq=extract.gt(pvcf, element='GQ', as.numeric=T)
#p.ad=extract.gt(pvcf, element='AD')
#p.ad1=masplit(p.ad, record=1, sort=0)
  # alt counts 
#p.ad2=masplit(p.ad, record=2, sort=0)
#get name of Leslie's strain

parents=DataFrame(chr=getCHROM(pvcf), 
                  pos=getPOS(pvcf), 
                  qual=getQUAL(pvcf),
                  id=getID(pvcf), 
                  ref=getREF(pvcf), 
                  alt=getALT(pvcf),
                  GT.BY=0, GT.AMH=1) #), GT=p.gt,DP=p.dp, GQ=p.gq,  AD1=p.ad1, AD2=p.ad2)
parents$marker.name=paste0(parents$chr, '_', parents$pos, '_', parents$ref, '_', parents$alt)
#----------------------------------------------------------------------------------------------------------

# Read in segregant data and reformat --------------------------------------------------------------------
segs=readRDS(seg.data)

segs$marker.name=paste0( segs$contig, '_',segs$position, '_', segs$refAllele, '_', segs$altAllele)
seg.s=split(segs, segs$sample)

# intitialize matrix to contain reformatted segregant data
seg.matrix=matrix(NA, nrow(parents), length(seg.s))

rownames(seg.matrix)=parents$marker.name
colnames(seg.matrix)=names(seg.s)

ref.matrix=seg.matrix
alt.matrix=seg.matrix

#chr11 - 570,000 #5621
# iterate through each segregant and fill matrix with read counts for the parental variants 
# also add a simple 
for(i in names(seg.s)){
    print(i)
    
    mind=match(seg.s[[i]]$marker.name, parents$marker.name)
    seg.matrix[mind[seg.s[[i]]$refCount>seg.s[[i]]$altCount],i]=0
    seg.matrix[mind[seg.s[[i]]$altCount>seg.s[[i]]$refCount],i]=1
    ref.matrix[mind,i]=seg.s[[i]]$refCount
    alt.matrix[mind,i]=seg.s[[i]]$altCount
}


# af.check=rowSums(ref.matrix, na.rm=T)/(rowSums(ref.matrix,na.rm=T)+rowSums(alt.matrix, na.rm=T))
# build a GRanges object for cross-referencing with parents genotype data 
segregants=makeGRangesFromDataFrame(DataFrame(chr=parents$chr, start=parents$pos, end=parents$pos))

# various QC filters here ---------------------------------------------------------------------------------
#NA out problematic regions 
seg.matrix[subjectHits(findOverlaps(blocked, segregants)),]=NA

#here we deal with converting reference and alt calls to parental 
g.recode=recode.as.parental(seg.matrix, parents)

s.cntNA=apply(g.recode,2, function(x) sum(is.na(x)))

s.fracNA=s.cntNA/nrow(g.recode)

# some aggregate stats across the panel, per marker
   m.cnt1=apply(g.recode,1, function(x) sum(x==1,na.rm=T))
   m.cnt2=apply(g.recode,1, function(x) sum(x==2,na.rm=T))
   m.cntNA=apply(g.recode,1, function(x) sum(is.na(x)))
   m.af = m.cnt1/(m.cnt1+m.cnt2)

   g.recode[which(m.cntNA>segregant.missing.by.marker.cutoff),]=NA
   g.recode[which(m.af>segregant.afH),]=NA
   g.recode[which(m.af<segregant.afL),]=NA
#   g.recode=g.recode[,-which(s.fracNA>segregant.missing.marker.fraction)]

#---------------------------------------------------------------------------------------------------------------

   
# use some helper functions from R/qtl for the HMM 

#manually remove region with the translocation, may need to redefine 
bad.int=101025:103418

#use the sum of the genotype calls in the translocation region to call the translocation event 
t.call=apply(g.recode[bad.int,],2, function(x) sum(x, na.rm=T)/sum(!is.na(x)))
t.discrete=as.numeric(cut(t.call, c(0,1.1,1.9,2)))-1
names(t.discrete)=names(t.call)


cross=hackCrossObject(g.recode[-bad.int,],parents[-bad.int,],unique.chrs)

cross=fake_gmap(cross)
cross =filterPoorQualitySegs(cross, error.prob1, max_XO_after_preimp=50)
# build a crude genetic map
cross =quickEstJB.mc(cross, ep=error.prob2)
cross =  argmax.geno(cross, step=0, map.function='kosambi', error.prob=error.prob2)
 #cross =  calc.genoprob(cross, step=0, map.function='kosambi', error.prob=.0075)
cross$geno = lapply(cross$geno, function(x) { rownames(x$argmax)=rownames(x$data); return(x) } )

cross$pheno$translocation=


#sanity check output 
genotable = hackScanOneObject(cross)
genotable$gpos=(gcoord.key[parents$chr]+parents$pos)[-bad.int]
    
pnames=names(parents)[grep('GT', names(parents))]
pnames=gsub('GT.', '', pnames)
plot(genotable$gpos, m.af[-bad.int], ylim=c(0,1), 
     main=cross.name, xlab='genome pos', ylab=paste(pnames[1], '/', pnames[1], '+', pnames[2]) )
points(genotable$gpos, genotable$AF, col='blue', ylim=c(0,1))
abline(v=gcoord.key, lty=2)



plotMap(cross, main='AMH x BY', sub=paste(pnames[1], 'X', pnames[2]) )
hmm.calls=do.call('cbind', lapply(cross$geno, function(x) x$argmax))
parents$gpos=(gcoord.key[parents$chr]+parents$pos)

i= cross$pheno$id[1]

#pdf(file='/data/yeast/Leslie/segplots.pdf', width=11, height=5)
for(l in 1:length(cross$pheno$id)){
    print(l)
    i= cross$pheno$id[l]
    png(file=paste0('/data/yeast/Leslie/segplots/', gsub('.bam.txt', '', i), '.png'), width=512*3, height=512)
    plot(parents$gpos, parents$GT.BY , ylim=c(-15,15), type='n', main=i, sub=
    paste('translocation genotype =', t.discrete[i]), 
         xlab='genomic coordinate', ylab='read depth (-BY)')
    text(gcoord.key+100000, rep(-12,16), names(gcoord.key))
    points(parents$gpos,-ref.matrix[,i], ylim=c(-15,15), cex=.35)
    points(parents$gpos, alt.matrix[,i], cex=.35) #, ylim=c(-20,20))
    points(genotable$gpos, 10*(hmm.calls[i,]-1.5), type='p', col=ifelse(hmm.calls[i,]==1,'orange', 'blue'))
    abline(v=gcoord.key, lty=2, col='lightblue')
dev.off()

}

cross$pheno$translocation=as.vector(t.discrete[cross$pheno$id])
cross$geno = lapply(cross$geno, function(x) { x$data=x$argmax; return(x) } )

write.cross(cross, format='csvsr', filestem='/data/yeast/Leslie/cross')
