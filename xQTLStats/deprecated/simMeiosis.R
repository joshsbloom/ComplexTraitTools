library(Meiosis)


set.seed(123L)  ## Seed R's rng
Meiosis::seed_rng(seed = 123L)  ## Seed rng used by Meiosis

n_chr <- 2L     ## number of chromosomes
L <- c(250,51)  #runif(n = n_chr, min = 100, max = 300)  ## sample length of chromosomes in cM
xoparam <- create_xoparam(L,obligate_chiasma=T)  ## no interference, no obligate chiasma
str(xoparam)


n_loci <- c(3000,10) #round(runif(n = n_chr, min = 1500L, max = 2500L))  ## sample number of loci
## sample positions of loci on the chromosome
positions <- lapply(seq_len(n_chr), function(i) sort(runif(n_loci[i], min = 0, max = L[i])))


ind <- replicate(2L, lapply(n_loci, function(n) sample(c(0L, 1L), n, replace = TRUE)), 
    simplify = FALSE)  ## simulate some genotypic data
str(ind)


ind=list(list(rep(0L,n_loci[1]),rep(0L,n_loci[2])),
         list(rep(1L,n_loci[1]),rep(1L,n_loci[2])))

p_geno <- Meiosis::cross_geno(father = ind, mother = ind, positions = positions, 
    xoparam = xoparam)
str(p_geno)

par(mfrow=c(2,1))
plot(p_geno[[1]][[1]])
#plot(p_geno[[1]][[2]])
plot(p_geno[[2]][[1]])
#plot(p_geno[[2]][[2]])

psize=5e5
p_geno <-replicate(psize, Meiosis::cross_geno(father = ind, mother = ind, positions = positions, xoparam=xoparam), 
                   simplify=F)
p_geno_par1=sapply(p_geno,function(x) x$paternal[[1]])

p_geno_par2=sapply(p_geno,function(x) x$paternal[[2]])

blist=list()
blist[[as.character(1)]]=p_geno


for(ngen in 2:48) {
    print(ngen)
    p_geno_h=list()
    for(h in 1:psize) {
        n1=sample(1:psize,1)
        n2=sample(1:psize,1)
        nind=list(list(blist[[as.character(ngen-1)]][[n1]]$paternal[[1]], blist[[as.character(ngen-1)]][[n1]]$paternal[[2]]),
                  list(blist[[as.character(ngen-1)]][[n2]]$maternal[[1]], blist[[as.character(ngen-1)]][[n2]]$maternal[[2]]))
      
        p_geno_h[[as.character(h)]] <- Meiosis::cross_geno(father = nind, mother = nind, positions = positions,  xoparam=xoparam)
    }
    blist[[as.character(ngen)]]=p_geno_h

    idx=sample(1:10000,1)
    par(mfrow=c(2,1))
    plot(blist[[as.character(ngen)]][[idx]][[1]][[1]], main=ngen)
    plot(blist[[as.character(ngen)]][[idx]][[2]][[1]])


    if(!is.null(blist[[as.character(ngen-2)]])){blist[[as.character(ngen-2)]]=NULL}
    print(object.size(blist)/1e6)
}

n_chr=2
f_alleles <- c(0L, 1L)  ## 21 and 65 are arbitrary integers
L=c(250,51)
f <- Meiosis::create_xo_founder(alleles = f_alleles, L = L)
xoparam <- create_xoparam(L,obligate_chiasma=T)  ## no interference, no obligate chiasma

n_loci <- c(3000,10) #round(runif(n = n_chr, min = 1500L, max = 2500L))  ## sample number of loci
## sample positions of loci on the chromosome
positions <- lapply(seq_len(n_chr), function(i) sort(runif(n_loci[i], min = 0, max = L[i])))

ind=list(list(rep(0L,n_loci[1]),rep(0L,n_loci[2])),
         list(rep(1L,n_loci[1]),rep(1L,n_loci[2])))

nc=1e6
mcross=replicate(nc, Meiosis::cross_xo(father = f, mother = f, xoparam = xoparam), simplify=F)

conv <- new(Meiosis::Converter, positions)  ## create a new converter object
conv$insert_founder(f_alleles, ind)  ## insert the (one and only) founder

convme=lapply(mcross, function(x) conv$convert(x))


str(conv$convert(mcross[[1]])) 



p_geno_par48=sapply(blist[["48"]],function(x) x$paternal[[1]])

































                    xoparam = xoparam)

par(mfrow=c(2,1))
plot(p_geno2[[1]][[1]])
plot(p_geno2[[2]][[1]])


p_geno[[500]][2]


n_founder=2
n_ind=1000
n_gen=2
alleles <- lapply(seq_len(n_founder), function(i) c(2L * i - 1L, 2L * i))
founder <- lapply(alleles, create_xo_founder, L = L)
system.time(syn <- make_synthetic(founder, n_ind, n_gen))
