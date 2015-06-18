library(OCplus)
A = MAsim.smyth(ng=15000, n=100, p0=0.2)
A = (A - min(A) + runif(1,0,1))/10

max_iter = 5000
k = 10
nsame = 100
m = nrow(A)
n = ncol(A)
W = as.double(runif(m*k, 0, 1)) # standard random init in (0,1)
H = as.double(runif(k*n, 0 ,1)) # could use libnmf's generatematrix for it

# calls to add: nmf_als, mu, neals, alspg, pg
#dyn.load("/nfs/research2/saezrodriguez/mike/ebits/stats/NMFconsensus/libnmf.so")
dyn.load("/nfs/research2/saezrodriguez/mike/ebits/stats/NMFconsensus/libnmf/nmf_mu.so")
re = .C("nmf_mu", a=as.double(A), w0=as.double(W), h0=as.double(H),
         pm=as.integer(m), pn=as.integer(n), pk=as.integer(k),
         maxiter=as.integer(max_iter), nsame=as.integer(nsame),
         DUP=FALSE, PACKAGE='nmf_mu')

#15kgenes
#fail: 200,175,160,158,156
#success: 100,150,155

#100 genes:
#39+,40x

#10 genes:
# 39+,40x
