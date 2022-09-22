#sample_well is well ID formatted like A01, B02, B11 etc 
get_96_Quadrant_in_384=function(Sample_Well_384) {
    qA=apply(expand.grid(toupper(letters[seq(1,16,2)]), sprintf("%02d",seq(1,24,2))),1, paste0, collapse='')
    qB=apply(expand.grid(toupper(letters[seq(1,16,2)]), sprintf("%02d",seq(2,24,2))),1, paste0, collapse='')
    qC=apply(expand.grid(toupper(letters[seq(2,16,2)]), sprintf("%02d",seq(1,24,2))),1, paste0, collapse='')
    qD=apply(expand.grid(toupper(letters[seq(2,16,2)]), sprintf("%02d",seq(2,24,2))),1, paste0, collapse='')
    

    quadrant_96=rep('', length(Sample_Well_384))
    quadrant_96[Sample_Well_384 %in% qA]='A'
    quadrant_96[Sample_Well_384 %in% qB]='B'
    quadrant_96[Sample_Well_384 %in% qC]='C'
    quadrant_96[Sample_Well_384 %in% qD]='D'
    
    return(quadrant_96)
}

#given sample_well and quadrant return 96 well plate position
get_96_pos_from_384_pos=function(Sample_Well, Quadrant_96) {

    Col=as.factor(gsub('^.', '', Sample_Well))
    #assumes 384-well plate layout
    Row=factor(gsub('..$', '', Sample_Well), levels=rev(toupper(letters[1:16])))

    Row96=Row
    for(l in c('A','B','C','D')){
          y=droplevels(Row[Quadrant_96==l])
           levels(y)=toupper(rev(letters[1:8]))
           Row96[Quadrant_96==l]=y
    }

    Col96=Col
        for(l in c('A','B','C','D')){
            y=droplevels(Col[Quadrant_96==l])
            levels(y)=sprintf('%02d', 1:12)
            Col96[Quadrant_96==l]=y
    }
    return(paste0(Row96, Col96))
}

# adding quadrant info and 96 well plate position info ---------------------------
path='/data/yeast/strains/'

strain.file = paste0(path, '1011_Sace_strains_matrix_positions_384.csv')
index.key=read.delim(strain.file, sep=';', stringsAsFactors=F)

colnames(index.key)[1]='Plate_384'
colnames(index.key)[5]='YJS_name'

index.key$row_384_alpha=toupper(letters[1:16])[index.key$row_384]
index.key$Sample_Well_384=paste0(index.key$row_384_alpha, sprintf('%02d',index.key$col_384))
index.key$Quadrant_96=get_96_Quadrant_in_384(index.key$Sample_Well_384)

index.key$Plate_96=paste(index.key$Plate_384, index.key$Quadrant_96, sep='_')
index.key$Sample_Well_96=get_96_pos_from_384_pos(index.key$Sample_Well_384, index.key$Quadrant_96) 
index.key$row_96=match(gsub('..$', '', index.key$Sample_Well_96),  toupper(letters[1:8]))
index.key$col_96=as.numeric(gsub('^.', '', index.key$Sample_Well_96))
    
# output augmented table 
write.table(index.key, file = paste0(path,'1011_Sace_strains_matrix_positions_384_augmented.csv'), 
            row.names=F, quote=F, sep=',')

# ----------------------------------------------------------------------------------


# plater doesn't have any write functions so here's one way to interconvert

mat96=matrix('',8,12)
rownames(mat96)=toupper(letters[1:8])
colnames(mat96)=1:12

mat384=matrix('',16,24)
rownames(mat384)=toupper(letters[1:16])
colnames(mat384)=1:24

s=split(index.key, index.key$Plate_96)
split96=lapply(s, function(x) {
       mega.name=paste(x$Strain_Name_4SRA, x$YJS_name, x$Standardized_name, sep=';')
       n=mat96
       n[cbind(x$row_96, x$col_96)]=mega.name 
       return(n)
       })


s=split(index.key, index.key$Plate_384)
split384=lapply(s, function(x) {
       mega.name=paste(x$Strain_Name_4SRA, x$YJS_name, x$Standardized_name, sep=';')
       n=mat384
       n[cbind(x$row_384, x$col_384)]=mega.name 
       return(n)
       })
options(width=1000)
sink(file=paste0(path,'plates96.csv')
split96
sink()

sink(file=paste0(path, 'plates384.csv')
split384
sink()
# manually reformatted the rest in vim (\s to , put plate name on same line as column numbers, remove trailing ',')


# then for reading in plates in a standardized format 
plates96=plater::read_plates('/data/yeast/strains/plates96.csv', well_ids_column='Sample_Well')

plates384=plater::read_plates('/data/yeast/strains/plates384.csv', well_ids_column='Sample_Well')






