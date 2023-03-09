
# number of plasmids
k=1e4
#total diversity
n=1e8


# simulation
x=seq(1,n)
h=c()
for(i in 1:100000) {
    print(i)
    y=rle(sort(sample(x,k,replace=T)))
    h=c(h,sum(y$lengths>1))
    
}
sum(h>0)/length(h)

# analytical
p=as.numeric(k)
for(i in 1:k){
    q=1-(0:(i-1))/n
    p[i]=1-prod(q) 
}
print(p[k])
