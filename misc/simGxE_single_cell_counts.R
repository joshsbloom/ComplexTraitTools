library(MASS)

#number of cells per experiment
ncells.per.exp=500

#simulate some spread in nUMIs per cell
nUMIs=1000+750*rexp(ncells.per.exp*2,1) #meh, something

# biological covariates of interest are two genotypes and two experimental conditions
# initialize
g1=rep(0,ncells.per.exp)
g2=rep(1,ncells.per.exp)
e1=rep(0,2*ncells.per.exp)
e2=rep(1,2*ncells.per.exp)

#one way to specify the experiment design matrix
X=cbind(log(nUMIs),g=c(g1,g2,g1,g2),e=c(e1,e2))
#specify GxE term
int= X[,'g']*X[,'e']
X=cbind(X,int)


#simulate some effects
betas=c(1, -1,2,-1)

#simulate some over-dispersion
BCV=.4
theta=1/BCV^2

#define some expected number of transcripts per total umi per cell
avg.exp.level.per.cell=.001

#get the expected mean of expression  
pMu=(X %*% betas)[,1]

#push into count space with neg bin over-dispersion 
ycounts=rnegbin(length(pMu),exp(pMu)*avg.exp.level.per.cell, theta=theta)

#fit a model of counts given this design matrix with a multiplicative overdispersion factor
m0=glm(ycounts~X, family=quasipoisson())

stripchart(ycounts~paste0(X[,'g'], X[,'e']), vertical=T, method='jitter', pch=19, cex=.2)

#fit a model of counts given this design matrix with quadratic overdispersion factor
m1=glm.nb(ycounts~X)

summary(m0)
summary(m1)

#alternative parameterization, just calculate group means 00,01,10,11
apX=cbind(log(nUMIs), model.matrix(~paste0(X[,'g'], X[,'e'])-1))
colnames(apX)[c(2:5)]=c('00','01','10','11')

m3=glm(ycounts~apX-1, family=quasipoisson())
m4=glm.nb(ycounts~apX-1)



#https://cran.r-project.org/web/packages/contrast/vignettes/contrast.html
#base R with interaction model is doing
#(00-01)-(10-11)
#+00
#-01
#-10
#+11
cVec=c(0,1,-1,-1,1)

Tval=(t(cVec)%*%m3$coefficients)/sqrt(t(cVec)%*%vcov(m3)%*%cVec)
pt(abs(Tval), df=m3$df.resid, lower.tail=F)
#see, same stat as in m0, ~ same p-val, different parameterization

#model.matrix(~paste0(X[,'g'], X[,'e']))
#table(paste0(X[,1],X[,2]))

#I'd probably have done it this way, note it just flips the sign of the T-stat
#interaction contrast
#(01-00)-(11-10)
#+01 -00 -11 +10

#rearrange
#-00
#+01
#+10
#-11
#cVec=c(0,-1,1,1,-1)
#Tval=(t(cVec)%*%m3$coefficients)/sqrt(t(cVec)%*%vcov(m3)%*%cVec)
#note this T is the same as for the interaction term in the model ycounts~X-1 (m0), but sign is flipped

