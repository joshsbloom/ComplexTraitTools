#strain selection
#greedy selection
dists = 1-(abs(cor.mat))
n=100
#subset=g.s

while(nrow(dists)> n) {
    print(nrow(dists))
    cdists=rowSums(dists)
# find the one that's closest to all the others and remove 
    closest=which(cdists == min(cdists))[1]
  #  subset=subset[-closest,]
    dists=dists[-closest,-closest]
}
selected.dists=dists
sum(selected.dists)


xx=rep(NA,500)
dists = 1-(abs(cor.mat))
for(i in 1:500) {
s100=sample(1:nrow(dists), 100)
xx[i]=sum(dists[s100,s100])
}


