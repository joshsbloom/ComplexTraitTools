power.list=list()

library(MASS)
library(ggplot2)
library(WriteXLS)
install.packages("remotes")
remotes::install_github("grayclhn/dbframe-R-library")
library(dbframe)
# remember dispersion is BCV^2
disps=c(.05, .2,.4,.6,.8,1,2)
#fold changes (added in depletion effects)
fcs=c(1/rev(c(1.01,1.05,1.1,1.2,1.5,2,4)),1.01, 1.05, 1.1, 1.2, 1.5,2,4)
nreps=c(5,25,50,100,200,500, 1000)
# average minimum number of reads per barcode
min.counts=c(10,50,100,1000)


eg=expand.grid(disps,fcs,nreps,min.counts)
colnames(eg)=c('disp', 'fc', 'nbarcode', 'min.count')
eg=data.frame(eg)
eg$power=0

for(i in 1:nrow(eg)) {
    #print(i)
     n=eg[i,'nbarcode']
     disp=eg[i,'disp']
     min.counts=eg[i,'min.count']
     fold=eg[i,'fc']     
     f=as.factor(c(rep('A',n),rep('B',n)))

     psim=RepParallel(1000, {
      # to deal with some stupid random error 
      tryCatch({
        # annoyingly rnegbin and glm define theta differently from edgeR ... as 1/disp
        y=rnegbin(2*n, mu=c(rep(1,n), rep(fold,n))*min.counts, theta=1/disp)

        ## assume pre-estimated dispersion, and estimated without error
         gmodel1=glm(y~f, family=negative.binomial(theta=1/disp)) 
         p=drop1(gmodel1, test='LRT')[[5]][2]
        return(p) 
        }, error=function(e) return(NA) )
     }, simplify='vector', mc.cores=35)
     psim=na.omit(psim)
     # note, this is alpha threshold is just a guess
     eg$power[i]=sum(psim<(.05/5000))/length(psim)
     print(eg[i,])
}
eg$BCV=as.factor(round(sqrt(eg$disp),2))
eg$disp=as.factor(eg$disp) #round(sqrt(eg$disp),2))

library(ggplot)
library(viridis)
ggplot(eg, aes(x=fc, y=power, color=disp))+geom_point()+
    geom_line()+
    facet_grid(rows=vars(min.count), cols=vars(nbarcode))+
   # round(,2)
   # [c(1:5,10:14)]
    xlab("fold change")+ scale_x_continuous(trans='log2', breaks=round(fcs,2) , limits=c(.8, 1.2))  +
    #,2), limits=c(.8,1.2)) +
   #  guide = guide_axis(n.dodge=4))+
    scale_color_viridis(discrete=T) +
     #scale_color_brewer(palette="Set1")+
    theme_bw()+
     theme(axis.text.x = element_text(angle = 45, hjust=1))+
    ggtitle('CRISPR editing design power analysis')

