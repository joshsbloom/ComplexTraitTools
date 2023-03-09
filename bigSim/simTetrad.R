c1=replicate(10000, sample(c(0,0,1,1), 4, replace=F))
c2=replicate(10000, sample(c(0,0,1,1), 4, replace=F))
c3=replicate(10000, sample(c(0,0,1,1), 4, replace=F))
c4=replicate(10000, sample(c(0,0,1,1), 4, replace=F))
c1=as.vector(c1)

c2=as.vector(c2)
c3=as.vector(c3)
c4=as.vector(c4)
tetrad=as.vector(rep(1:10000, each=4))
table(sort(rle(tetrad[which(c1*c2==0)])$lengths))/100
table(sort(rle(tetrad[which(c1*c2*c3==0)])$lengths))/100
table(sort(rle(tetrad[which(c1*c2*c3*c4==0)])$lengths))/100
102-44
table(sort(rle(tetrad[which(c2*c3*c4==0)])$lengths))/100
 
