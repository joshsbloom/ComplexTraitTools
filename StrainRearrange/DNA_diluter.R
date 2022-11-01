library(plater) 
library(tidyr)

qubit.data=read_plates('/data/yeast/oxford/Oxford_Plate_1234.csv')
outfile=paste0('/data/yeast/oxford/Oxford_Plate_1234_dil_4.csv')

#specify an output file

standardNGin = c(0,10,20,40,80,100)
std.wells=paste0('G',sprintf('%02d',1:6))

qubit.std.vals=qubit.data$Plate_4[qubit.data$Wells %in% std.wells]


g1model = lm(standardNGin~qubit.std.vals-1)
plot(qubit.std.vals, standardNGin)
abline(g1model,col='red')

DNA=qubit.data
DNA[,-c(1,2)]=DNA[,-c(1,2)]*as.numeric(g1model$coeff[1])
DNA[,-c(1,2)][DNA[,-c(1,2)]==0]=0.000001

total.DNA=200 # 200ng total DNA
rxn.vol=100   # in final volume of 100uL
# so this gets us 2ng/uL final conc


DNA.needed=total.DNA/DNA[,-c(1,2)]
DNA.needed[DNA.needed>rxn.vol]=rxn.vol
H2O.needed=rxn.vol-DNA.needed

DNA.needed=cbind(DNA[,c(1,2)], DNA.needed)
H2O.needed=cbind(DNA[,c(1,2)], H2O.needed)

x=DNA.needed %>% 
    pivot_longer(cols=starts_with('Plate_'), names_to="DNAsource", values_to="DNAvol")
y=H2O.needed %>% 
    pivot_longer(cols=starts_with('Plate_'), names_to="holder", values_to="H2Ovol")

robot.df=data.frame(sourceplate=x$DNAsource, watersource="watersource", receivingplate=paste0('dil_', x$DNAsource),
            startpos=x$Wells, DNAvol=x$DNAvol, H2Ovol=y$H2Ovol, endpos=x$Wells , stringsAsFactors = FALSE)
#reorder 
robot.df=robot.df[order(robot.df$sourceplate, robot.df$startpos),]

#outfile needs to be csv file name out.csv 

#adding eol argument here should obviate need to run unix2dos
#this will output one file for multiple input plates
write.table(robot.df, outfile, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE, eol = "\r\n")


# this will output one file for each input plate 
rdsf=split(robot.df, robot.df$sourceplate)
for(p in names(rdsf)){ 
    outfile=paste0('/data/yeast/oxford/Oxford_Plate_1234_dil_split_', p, '.csv')
    write.table(rdsf[[p]], outfile, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE, eol = "\r\n")
}

# Old instructions that we used before specifying eol argument in write.table
# Using terminal, drag the folder from Finder so it changes the working directory (cd .....) 
# use command "cat -t 'filename'" to view format (if in unix, doesn't show a line break character)
# use following command "unix2dos 'filename'" to change it to DOS format
# use "cat -t 'filename'" to check it's been changed (should have "^M" at the end of every line)


