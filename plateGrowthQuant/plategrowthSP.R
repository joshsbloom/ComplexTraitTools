######## <PART 1>  ###################################################################################################################
# Image processing code V3 ... greatly simplified from V2 using new EBImage functions and kmeans for image segmentation
library(EBImage)
library(gdata)
library(locfit)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=6)


#location of file with functions 
fx.file='/data/yeast/Colonies/code/plategrowthSPFx.R'

#base directory
base.path=c('/data/yeast/Colonies/2022-07-08/')
#../base.path/
#  /base.path/Images/*.JPG
#  /base.path/keys/Key.csv
#  /base.path/keys/PlateLayoutM1.csv .. etc 
#  /base.path/out/

source(fx.file)


image.path=paste(base.path, 'Images/', sep='')
dir.create(image.path)
layout.path=paste(base.path, 'keys/', sep='')
dir.create(layout.path)
key.file = paste(layout.path, 'Key.csv', sep='')
Results.file =  paste(base.path, 'Results.RDS', sep='')
out.dir=paste(base.path, 'out/', sep='')
dir.create(out.dir)

key=read.delim(key.file, stringsAsFactors=F, header=T, sep=',')

# Run image processing with parallelization
# can replace %dopar% with %do% to run single-threaded
Results=foreach(i=1:nrow(key)) %dopar% { 
        processImages(i, key, image.path, layout.path, out.dir,
                      plate.types[['384']], corners[['vwr']], 
                      flip=F, save.labeled=T) 
}
plate.names=apply(key, 1, paste, collapse='::' )
names(Results)=plate.names
#------------------------------------------------------------------------------------------------------------------------

#save output as a list per plate 
saveRDS(Results, file=Results.file)

#Results=readRDS(Results.file)
results.df=data.table::rbindlist(lapply(Results, tibble::rownames_to_column), idcol='Plate')
# new 05/03/22
#dumb munging stuff for WI collection ---------------------------------------------------------------
    strain.id=data.frame(data.table::tstrsplit(results.df$rowname, ';', type.convert=T,names=T))
    colony.position=as.numeric(gsub('(\\d{1-3}):(.*)', '\\1' , strain.id[,1]))
    #ref.num=t(matrix(as.character(1:384), 24,16))

    ref.row=rep(toupper(letters[1:16]),each=24)
    ref.col=rep(c(1:24),times=16)

    #add row and col
    colony.row=ref.row[colony.position]
    colony.col=ref.col[colony.position]

    #parse yeast isolate strain names 
    StrainName=gsub('(\\d{1-3}):(.*)', '\\2' , strain.id[,1])
    strain.id=data.frame(colony.position, colony.row, colony.col, StrainName, strain.id[,c(2,3)])
    names(strain.id)[c(5,6)]=c('StrainName.JS.long', 'StrainName.JS.short')

    #parse plate names 
    parsed.plate.name=data.frame(data.table::tstrsplit(results.df$Plate, '::', type.convert=T, names=T))
    names(parsed.plate.name)=c('layout', 'permutation', 'file', 'condition', 'condition1', 'control')
    #specific to the isolate phenotyping 
    results.df=cbind(parsed.plate.name, strain.id, results.df[,-c(1,2)])
#-------------------------------------------------------------------------------------------

#save long formatted data as an RDS
saveRDS(results.df, file=paste(base.path, 'ResultsTable.RDS', sep =''))

# START HERE to skip image processing ------------------------------------------


library(tidyverse)
results.df=readRDS(paste(base.path, 'ResultsTable.RDS', sep =''))

