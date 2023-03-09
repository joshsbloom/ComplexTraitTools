#testx=svd(t(scale(pg))
#U=testX$u


x.scaled=scale(t(pg))
x.t=t(x.scaled)

y=simBayesC(x.t, target.size=.05, h2=.8/16)

#sxt=svd(x.scaled)
#U=sxt$u
#S=diag(sxt$d)
#Vt=t(sxt$v)

fitted= susie(x.scaled, y[,1], L=200, verbose=T)
b=rep(0, ncol(x.scaled))
b[attr(y, 'causal.indices')]=attr(y, 'causal.betas')
susie_plot(fitted, y='PIP', b=b)
plot(cor(y,x.scaled)[1,])
