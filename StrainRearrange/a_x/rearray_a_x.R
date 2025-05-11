library(tidyverse)
library(plater)
library(data.table)
#libreoffice --headless --convert-to csv --outdir /home/jbloom/Dropbox/code/ComplexTraitTools/StrainRearrange/a_x/ a_x\ 500\ collection_0509.xlsx
plates=read_plates('/home/jbloom/Dropbox/code/ComplexTraitTools/StrainRearrange/a_x/a_x_500_collection_050925.csv')

allwells=plates$Wells
blanks=c(16,58)
allwells[blanks]
nblanks=allwells[-blanks]

plates_longer=data.frame(pivot_longer(plates[,-1], cols=-Wells, names_to="Source Plate", values_to="strain"))
plates_longer=plates_longer[!is.na(plates_longer[,3]),]
plates_longer$Source.Plate = gsub('Plate ', '', plates_longer$Source.Plate)
plates_longer$Source.Plate = gsub('plate ', '', plates_longer$Source.Plate)

unfuck=strsplit(plates_longer$Source.Plate , '\\.')
unfuck=sapply(unfuck, function(x)  paste0(sprintf('%02d',as.numeric(x[1])), '_', sprintf('%02d',as.numeric(x[2])) ))
plates_longer$Source.Plate=unfuck
plates
plates_longer=plates_longer[order(plates_longer$Source.Plate),]
plates_longer=plates_longer[!grepl('^ $', plates_longer$strain),]
plates_longer$mating=NA
plates_longer$mating[grepl(' x$', plates_longer$strain)]='x'
plates_longer$mating[grepl(' a$', plates_longer$strain)]='a'
plates_longer$drug=NA
plates_longer$drug[plates_longer$Source.Plate %in% unique(plates_longer$Source.Plate)[1:27]]='KAN'
plates_longer$drug[plates_longer$Source.Plate %in% unique(plates_longer$Source.Plate)[38:40]]='KAN'
plates_longer$drug[plates_longer$Source.Plate %in% unique(plates_longer$Source.Plate)[30:37]]='HYG'
plates_longer$drug[plates_longer$Source.Plate %in% unique(plates_longer$Source.Plate)[41:43]]='HYG'
plates_longer$drug[plates_longer$Source.Plate %in% unique(plates_longer$Source.Plate)[28:29] & grepl('A|B|C|D', plates_longer$Wells)]='HYG'
plates_longer$drug[plates_longer$Source.Plate %in% unique(plates_longer$Source.Plate)[28:29] & grepl('E|F|G|H', plates_longer$Wells)]='KAN'

plates_longer=plates_longer[!is.na(plates_longer$mating) & !is.na(plates_longer$drug),]
names(plates_longer)[1]='Source.Wells'

spl=split(plates_longer, paste0(plates_longer$mating, '_', plates_longer$drug))

for (n in names(spl)){
           x=spl[[n]]
           nplates=nrow(x)/length(nblanks)
           Target.Wells=rep(nblanks, ceiling(nplates))
           Target.Plate=paste0(n, '_', rep(1:ceiling(nplates), each=length(nblanks)))
           x$Target.Wells=Target.Wells[1:nrow(x)]
           x$Target.Plate=Target.Plate[1:nrow(x)]
           x$Volume=2
           spl[[n]]=x
}

spl_fat=data.frame(data.table::rbindlist(spl))
spl_fat=spl_fat[order(spl_fat$Source.Plate, spl_fat$Source.Wells),]
outfile=paste0('/home/jbloom/Dropbox/code/ComplexTraitTools/StrainRearrange/a_x/master_a_x_rearray', '.csv')
write.table(spl_fat, outfile, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE, eol = "\r\n")
dir.create('/home/jbloom/Dropbox/code/ComplexTraitTools/StrainRearrange/a_x/Biomek/')

spl_sf=split(spl_fat, spl_fat$Target.Plate)
for(p in names(spl_sf)){ 
    outfile=paste0('/home/jbloom/Dropbox/code/ComplexTraitTools/StrainRearrange/a_x/Biomek/', p, '.csv')
    write.table(spl_sf[[p]], outfile, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE, eol = "\r\n")
}

