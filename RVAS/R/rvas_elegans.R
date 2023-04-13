library(tidyverse)
library(data.table)
library(SKAT)

#whole genome variant annotation file
var.annot.file='/data/elegans/RVAS/WI.20220216.strain-annotation.tsv.gz'
#GRM info
GRM.file='/data/elegans/RVAS/gemmaGRM.sXX.txt'
GRM.rowcolnames.file='/data/elegans/RVAS/traits.sample'

#Phenotype(s)
pheno.file='/data/elegans/RVAS/ppat.1007226.s004.tsv'

#read in big annotation file
e.annot=read_delim(var.annot.file, delim='\t')
#pre-filter for speed 
e.annot.f=e.annot[!is.na(e.annot$WORMBASE_ID) & !is.na(e.annot$CONSEQUENCE) & e.annot$CONSEQUENCE!='synonymous', ]

#read in GRM
grm=as.matrix(read_delim(GRM.file, delim='\t', col_names=F))

grm_names=read_delim(GRM.rowcolnames.file, delim=' ')[-1,]

pheno=read_delim(pheno.file, delim='\t')

#add rownames to grm
colnames(grm)=rownames(grm)=grm_names$ID_1

#fewer genes than wormbase IDs
e.annot.fg=split(e.annot.f, e.annot.f$GENE)
nvar=sapply(e.annot.fg, function(x) nrow(x))

#filter genes with not enough variants, could be more strict here 
e.annot.fg=e.annot.fg[nvar>2]


#i see N2 has almost no variants relative to ref so we need to go deep into the genome before we see an N2 variant
# made this before I got vector of strain names from GRM
ustrains=sort(unique(unlist(sapply(e.annot.f$Strains[1:2e5], function(x) strsplit(x,',')))))

#cleanup
rm(e.annot)
rm(e.annot.f)


#reorder grm by ustrains
grm=(grm[ustrains,ustrains])


# build SNP incident matrices for each gene, store the matrices as sparse matrices
# i added a column name with a new snp id, might be useful for crossref at some point downstream
nonsyn.matrices=list()
for( g in names(e.annot.fg)) {
    print(g)
    gd=e.annot.fg[[g]]
    Z=matrix(0, length(ustrains), nrow(gd))
    rownames(Z)=ustrains
    colnames(Z)=paste(gd$CHROM, gd$POS, gd$REF, gd$ALT, sep="_")
    umatch=lapply(strsplit(gd$Strains, ','), function(x) match(x, ustrains))
    for(i in 1:length(umatch)){      Z[umatch[[i]],i]=1 }
    nonsyn.matrices[[g]]=Matrix(Z,sparse=T)
}


#remove phenotypes for genos that don't exist in the grm
pheno=pheno[(pheno$strain %in% rownames(grm)),]

#reorder grm again, yay
grm.s=grm[pheno$strain,pheno$strain]

#intercept effect in SKAT model, could have other experimental or PC based covariates here 
X=matrix(1, nrow(grm.s))

#get a phenotype

y=pheno$Albendazole_q90.TOF


#specify null model
obj=SKAT_NULL_emmaX(y ~ X, K=grm.s) #

#put skat p-vals here
skat_ps=list()
for( g in names(e.annot.fg)) {
    print(g)
    Z=nonsyn.matrices[[g]]
    Zs=Z[pheno$strain,]
    Zs=Matrix(Zs, sparse=T)

    skat_ps[[g]]=SKAT(Zs, obj, is_dosage=T, is_check_genotype=F)$p.value
}
