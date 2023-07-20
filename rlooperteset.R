Ball1 = Ball
#a = 0

#G = a + B

totalseq = 2000
# m = 500 #500 bp Rloop
#n+m = 1800
Bmid  = sample(Bbin,size=as.integer(0.5+0.3*totalseq),replace=T)
Btop  = sample(Bbin[seq(5,9)],size=as.integer(0.5+0.2*totalseq),replace=T)
Bmid2 = sample(Bbin,size=as.integer(0.5+0.5*totalseq),replace=T)
Bsim1  = c(Btop, Bmid, Bmid2)
Bsim2  = c(Bmid,Btop,Bmid2)
Bsim3 = c(Bmid,Bmid2,Btop)
Bsim4 = rep(0.99,totalseq)
Bsim5 = rep(0.5,totalseq)
Bsim6 = rep(0.25,totalseq)
Bsim7 = rep(0,totalseq)
Bsim8 = rep(-0.25,totalseq)
Bsim9 = rep(-0.5,totalseq)
Bsim10 = rep(-0.99,totalseq)
bf = function(a,B) {return(exp(-1*(a+B)/310))}

sigma = -0.07
N = 1500
A = 1/10.4
C = 1.8
T = 310
K = (2200 * 0.0019858775 * T) / N;
alpha = N*sigma*A


B = data.frame()
Ball0 = 0.5*K*alpha^2
Ball = c(Ball0)


for (n in seq(2,totalseq)) {
  if (n %% 10 == 0) {print(n)}
  if (n == 1) {Ball = c(Ball0)}
  
  set.seed(n)
  
  #m = min(n+1500,totalseq)
  nmax = totalseq
  
  
  B1 = Bsim[seq(n,nmax)]

  m = seq(n-n,nmax-n)+1
  
  top = (2*pi^2*C*K) * (alpha+m*A)^2
  bot = (4*pi^2*C+K*m)
  B2 = (top/bot)
  
  Ball[n] = sum(Bsim[seq(n,nmax)] + ( ((2*pi^2*C*K) * (alpha+m*A)^2) / (4*pi^2*C+K*m) ) )
  # Ball[n] = sum(B1 + B2)
}
Bsim = Bsim4
B = rbind(B,data.frame(x=0.99,Ball=Ball,type='bsim7'))
B$index = rep(seq(1,totalseq),dim(B)[1]/totalseq)
ggplot(B,aes(index,Ball)) + geom_line(color=as.factor(type))

as = c(0,10,10.5,11,1000)
for (i in seq(1,length(as))) {
  if (i == 1) {G = data.frame()}
  a = as[i]
  #G = rbind(G,data.frame(a=as[i],B=B))
  bftop = bf(a=a,B=B$Ball)
  part.func = sum(bftop)
  
  Gcurr = data.frame(a=as[i],Gs = Ball,bftop=bftop,part.func=part.func,p=bftop/part.func,index=seq(1,length(bftop)))
  G = rbind(G,Gcurr)
}
#ggplot(G[G$index <= 2500,],aes(index,p)) + geom_line(aes(color=as.factor(a))) + 
#  coord_cartesian(ylim=c(0,1e-8)) +facet_grid(as.factor(a)~.)
ggplot(G,aes(index,log(1/(1e-23+Gs)))) + geom_line(aes(color=as.factor(a)))#) + facet_grid(as.factor(a)~.)
#  coord_cartesian(ylim=c(0,1e-2))


Bbin = c(-0.36,-0.16,-0.1,-0.06,0.97,0.334,0.45,0.38,-0.12,-0.16,0.6,-0.12,0.45,0.5,0.28,0.8)
nuc = c('G','C','A','T')
B = list()
for (i in seq(1,4)) {
  nuc1 = nuc[i]
  for (j in seq(1,4)) {
    nuc2 = nuc[j]
    B[[nuc1]][[nuc2]] = Bbin[i*j]
  }
}
