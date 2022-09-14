# individs in pools
N=1e6
# contrast mode
# resolution
res=100
#centimorgans
cM=2200

p= res/100.0/cM
#peak pos = 226250 
#high.file='/data/multipool/row_1335iteration5N1e+05depth1000effect0.25HGmultiPoolIn.txt'
#low.file= '/data/multipool/row_1335iteration5N1e+05depth1000effect0.25LGmultiPoolIn.txt'

high.file='/data/multipool/row_3518iteration8N1e+06depth1e+06effect1HGmultiPoolIn.txt'
low.file= '/data/multipool/row_3518iteration8N1e+06depth1e+06effect1LGmultiPoolIn.txt'

readf = function(fin,res) {
    dfs    = read.delim(fin, header=F, sep='\t')
    nout   = round(max(dfs[,1])/res)+1
    bmatch = dfs[,1]/res
    bmatch.int =findInterval(bmatch, 1:nout, all.inside=T)

    means  = rep(0, nout)
    means[bmatch.int]=dfs[,2]                # np  (i.e. 1000 * .5)
   
    counts = rep(0, nout)
    counts[bmatch.int]= dfs[,2] + dfs[,3]      # n   (i.e. 1000)
   
    p =rep(0,nout)
    p[bmatch.int] = dfs[,2]/(dfs[,2]+dfs[,3])   # p   (i.e. .5)
   
    variances =  counts * p *( 1-p)  # np(1-p)  (i.e  1000*.5(1-.5)
    variances[variances==0]=NA
    # equal to y, y_var, and d 
    # means = y
    # variances = y_var 
    # counts = d 
    return(list(means=means, variances=variances, counts=counts))
}

hin = readf(high.file,res)
lin = readf(low.file,res)

y=lin$means
y_var=lin$variances
d=lin$counts

y2=hin$means
y_var2=hin$variances
d2=hin$counts

TT=length(hin$means)

kalman=function(y, y_var, d, TT, N, p) {
    mu = rep(0,TT)#numpy.zeros(T)
    V = rep(0,TT)#numpy.zeros(T)
    P = rep(0,TT)# numpy.zeros(T)

    V_pstr = rep(0,TT)# numpy.zeros(T)
    mu_pstr = rep(0,TT) #numpy.zeros(T)

    cc = rep(1,TT)

    mu_initial = 0.5*N # Initial parameters, assumed given (binomial distribution)
    V_initial = 0.25*N

    A = (1.0 - 2.0*p)
    C = 1.0 * d / N
    S = p*(1.0-p)*N
    
    explode=C[1]/(C[1]^2.0*V_initial + ifelse(is.na(y_var[1]),0, y_var[1]))
    
    K = V_initial* ifelse(is.na(explode),0, explode)
    # C[1]/(C[1]^2.0*V_initial + ifelse(is.na(y_var[1]),0, y_var[1]))
    mu[1] = mu_initial + K*(y[1] - C[1]*mu_initial)
    V[1] = (1.0-K*C[1])*V_initial
    if(is.na(y_var[1])) { cc[1]=1 } else {
        cc[1]= dnorm(y[1], C[1]*mu_initial, sqrt(C[1]^2.0*V_initial + y_var[1])) }
    
    # Forward pass:
    for (i in 2:(TT)) {
        if (i == 2) { 
            P[i-1] = A^2.0*V_initial + S }
        else { 
            P[i-1] = A^2.0*V[i-1] + S 
        }

        if ( is.na(y_var[i] )) {
            K = 0
            cc[i] = 1.0
         }
        else{
            K = P[i-1]*C[i]/(C[i]^2.0*P[i-1]+y_var[i])
            cc[i] = dnorm(y[i], C[i]*(A*mu[i-1]+p*N), sqrt(C[i]^2.0*P[i-1] + y_var[i]))
        }
        mu[i] = A * mu[i-1] + N*p + K * (y[i] - C[i]*(A*mu[i-1] + N*p))
        V[i] = (1.0-K*C[i])*P[i-1]
     }

    V_pstr[length(V_pstr)] = V[length(V)]
    mu_pstr[length(V_pstr)] = mu[length(V)]
    logLik = sum(log(cc))

    # Backwards pass:
    for (i in (TT-1):1) {
        J = V[i]*A/P[i]
        mu_pstr[i] = mu[i] + J * (mu_pstr[i+1] - A*(mu[i]) - N*p)
        V_pstr[i] = V[i] + J^2.0 * (V_pstr[i+1] - P[i])
    }
    
    return( list(mu_pstr=mu_pstr, V_pstr=V_pstr, logLik=logLik))

}

k1=kalman(y,y_var,d, TT, N,p)
k2=kalman(y2,y_var2,d2, TT, N,p)
  

lognormpdf=function(x, mu, sigma){-0.5*log(2*pi) - log(sigma) + (-(x-mu)^2.0/2.0/sigma^2.0) }

calcLOD = function(mu_pstr_vec, v_pstr_vec, TT, N) {
  #mu_pstr_vec =  k1$mu_pstr
  #  v_pstr_vec = k1$V_pstr
    LOD = rep(0,TT) #numpy.zeros(T)
    mu_MLE = rep(0,TT) #numpy.zeros(T)

    mu_initial = 0.5*N
    V_initial = 0.25*N
    delta = 0.0025
    x=seq(delta, 1-delta+delta/2, delta)
    p_alt=x

    # row-wise
    log_p_precomp=(sapply(x, function(x) dnorm(N*x, N*p_alt, sqrt(p_alt*(1.0-p_alt)*N),log=T)))
    #p_precomp=(sapply(x, function(x) dnorm(N*x, N*p_alt, sqrt(p_alt*(1.0-p_alt)*N))))
    logreweighter = lognormpdf(N*x, mu_initial, sqrt(V_initial))

    if(is.null(dim(mu_pstr_vec))) {
        for (i in 1:TT) {
            logallsums = rep(0, length(x))
            logtemp = lognormpdf(N*x, mu_pstr_vec[i], sqrt(v_pstr_vec[i])) - logreweighter
            scaler  = max(logtemp)
            #logallsums = logallsums + scaler+ log(p_precomp %*% exp(logtemp - scaler))
            logallsums = logallsums + scaler+ log(jbMult(log_p_precomp,logtemp-scaler  )) 
            p_alt=x[which.max(logallsums)]*N
            mu_MLE[i] = p_alt
            LOD[i] = log10(N) + log10(x[2]-x[1])+ max(logallsums) / log(10.0)
        }
    } else {
        for (i in 1:TT) {
            logallsums = rep(0, length(x))
            for(j in 1:2) {
            logtemp = lognormpdf(N*x, mu_pstr_vec[i,j], sqrt(v_pstr_vec[i,j])) - logreweighter
            scaler  = max(logtemp)
            lst=log(jbMult(log_p_precomp,logtemp-scaler  )) 
            lst[!is.finite(lst)]=-1e-320
            logallsums = logallsums + scaler+ lst
            }
            p_alt=x[which.max(logallsums)]*N
            mu_MLE[i] = p_alt
            LOD[i] = log10(N) + log10(x[2]-x[1])+ max(logallsums) / log(10.0)

        }
    }

    return(list(LOD=LOD, mu_MLE=mu_MLE))
}
jbMult=function(a,b) {apply(a, 1, function(x) { sum(exp((x+b))) }) } 

L1=calcLOD( k1$mu_pstr, k1$V_pstr, TT, N)
L2=calcLOD( k2$mu_pstr, k2$V_pstr, TT, N)
L3= calcLOD( cbind(k1$mu_pstr,k2$mu_pstr), cbind(k1$V_pstr,k2$V_pstr), TT, N)
L=(L1$LOD + L2$LOD) - L3$LOD
L=(L1$LOD + L2$LOD)
plot(L2$mu_MLE, ylim=c(0,1000))
points(y2, col='red')

mu_pstr_vec=cbind(k1$mu_pstr,k2$mu_pstr)
v_pstr_vec =cbind(k1$V_pstr,k2$V_pstr)


L1.1=calcLOD( k1$mu_pstr, k1$V_pstr, TT, N)
L2.1=calcLOD( k2$mu_pstr, k2$V_pstr, TT, N)
L3.1= calcLOD( cbind(k1$mu_pstr,k2$mu_pstr), cbind(k1$V_pstr,k2$V_pstr), TT, N)
L.1=(L1.1$LOD + L2.1$LOD) - L3.1$LOD
plot(L.1)

plot(L2$mu_MLE, ylim=c(0,1000))
points(y2, col='red')
