library(reshape)
library(gdata)

# Save the excel file as a text file, tab delimited
# First tab should be a 96-well plate with plate reader intensities
# Second tab should be the columns of wells corresponding to the standards 

# Should have the same layout as the plate with concentrations in the "cells"

# Change this each time to the directory containing xlsx files
dir.in = '/data/rr/DNAquant/13/'
fs = list.files(dir.in)

for (filein in fs) { 

DNAfile = paste(dir.in, filein, sep='')
print(DNAfile)
outfile = gsub('xlsx', 'txt', DNAfile)

# How much will DNA stock be diluted?
# What is the absolute amount of DNA (in ng) and in what volume (rxn.vol)?
parameters <- list(total.DNA = 200, rxn.vol = 120)
attach(parameters)

intensity = read.xls(DNAfile, header=F,sheet=1)
standard  = read.xls(DNAfile, header=F,sheet=2)
# standard in ng/ul

#standardNGin = c(0,5,10,20,40,60,80,100)
standardNGin = c(0,5,10,20,40,60,80,100)

g1=apply(standard[,c(1:6)], 1, mean)
g1model = lm(standardNGin[-c(7,8)]~g1[-c(7,8)]-1)

#plot(g1[-c(7,8)], standardNGin[-c(7,8)], col='red')
#for(i in 1:3) { points(standard[-c(7,8),i], standardNGin[-c(7,8)]) }
#abline(g1model)

DNA  = cbind( intensity[,1:12]*as.numeric(g1model$coeff[1]))
DNA[DNA==0]=0.000001
# Use only plate 1. Determine the volume needed per individual well.
#DNAnew <- DNA[DNA$Plate == 1,]
DNAneeded <- total.DNA/DNA
DNAneeded[DNAneeded>120]=120
# Reconfigure info into lists compatible with the robot.
positions <- as.character(melt(
                               sapply(1:12, function(num) {
                                      sapply( toupper(letters[1:8]),
                                        function(letter) {paste(letter, num, sep = "")}, USE.NAMES = FALSE)}, 
                                                            USE.NAMES = FALSE))$value)
DNAneeded <- melt(DNAneeded)$value
#DNAneeded <- sapply(DNAneeded, function(well) {if (well > rxn.vol) {rxn.vol} else {well}})
H2Oneeded <- (rxn.vol - DNAneeded)
robot.df <- data.frame(sourceplate = "DNAsource", watersource = "watersource", receivingplate = "DNAreceiving", 
        startpos = positions, DNAvol = DNAneeded, H2Ovol = H2Oneeded, endpos = positions, stringsAsFactors = FALSE)
write.table(robot.df, outfile, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)

}
# Using terminal, drag the folder from Finder so it changes the working directory (cd .....) 
# use command "cat -t 'filename'" to view format (if in unix, doesn't show a line break character)
# use following command "unix2dos 'filename'" to change it to DOS format
# use "cat -t 'filename'" to check it's been changed (should have "^M" at the end of every line)
