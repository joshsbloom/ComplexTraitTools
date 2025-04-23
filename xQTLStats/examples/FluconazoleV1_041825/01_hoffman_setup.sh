#Connect to hoffman2:
ssh username@hoffman2.idre.ucla.edu
#request a compute node
#requests a node for 48 hours with 8G ram and 8 threads
qrsh -l highp,h_rt=48:00:00,h_data=8G -pe shared 8

#on compute node, assuming you have conda in
#assuming you ha
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
module load mamba

#create a conda environment, install a bunch of dependencies, use mamba 
mamba create --name genomics bwa gatk4 bcftools sambamba r-tidyverse r-qs r-vcfr r-alphasimr r-rfast r-qs r-ggpubr r-devtools 

mamba activate genomics 

#start R, in R 
devtools::install_github("joshsbloom/ComplexTraitTools/xQTLStats" ,ref="main")


#for basespace download from personal directory on basespace 
bs auth 
bs list projects
#something like this 
bs download project -i 451502321 -o /u/project/kruglyak/username/projects/FluconazoleV1_041825/fastq/ --extension fastq.gz

