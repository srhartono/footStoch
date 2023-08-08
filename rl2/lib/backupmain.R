source('lib/lib9.R')
source('lib/lib1.R')


fa.get_chunk = function(seq1,window=200) {
  if (window > nchar(seq1)) {
    window = nchar(seq1)-1
  }
  df = data.frame(i0 = seq(1,(nchar(seq1)-window)),i1 = seq((1+window),nchar(seq1)))
  df$chunk = apply(df,1,function(x,seq1){substr(seq1,x[1],x[2])},seq1=seq1)
  df2 = t(sapply(df$chunk,fa.calc_CpGprof))
  rownames(df2) = rownames(df)
  colnames(df2) = c('CGdens','GCcont','GCskew')
  df = cbind(df,df2)
  return(df)
#  chunk = sapply(df$i0,substr(),start=i0,stop=)
}
fa.calc_CpGprof = function(seq1) {
  N  = nchar(seq1)
  nA = str_count(seq1,'A')
  nC = str_count(seq1,'C')
  nG = str_count(seq1,'G')
  nT = str_count(seq1,'T')
  nN = str_count(seq1,'N')
  nCG = str_count(seq1,'CG')
  CGdens = 0
  GCcont = 0
  GCskew = 0
  if (N > 0) {
    if (nC/N*nG/N > 0) {
      CGdens = ai((nCG/N)/(nC/N*nG/N)*100)/100
    }
    if (nC+nG > 0) {
      GCskew = ai((nG-nC)/(nG+nC)*100)/100
    }
    GCcont = ai((nC+nG)/N*100)/100
  }
  return(c(CGdens,GCcont,GCskew))
}

#fa.calc_CpGprof(seq1)

# Generate random seq1
set.seed(420)
nuc.bin1 = c('A','C','T','G')
nuc.bin2 = c(rep('G',9),rep('C',3),'A','T')
nuc.bin3 = c(rep('G',6),rep('C',6),'A','T')
seq1a = join(sample(nuc.bin1,size = 200,replace=T))
seq1b = join(sample(nuc.bin1,size = 200,replace=T))
seq1c = join(sample(nuc.bin1,size = 200,replace=T))
seq1d = join(sample(nuc.bin1,size = 200,replace=T))
seq1e = join(sample(nuc.bin1,size = 200,replace=T))
seq1f = join(sample(nuc.bin1,size=200,replace=T))
seq1g = join(sample(nuc.bin1,size=800,replace=T))

seq2a = join(sample(nuc.bin2,size=200,replace=T))
seq2b = join(sample(nuc.bin2,size=200,replace=T))
seq3a = join(sample(nuc.bin3,size=200,replace=T))
seq3b = join(sample(nuc.bin3,size=200,replace=T))

seqT = join(rep('T',2000))
seqG = join(rep('G',2000))
seq1 = join(c(seq1a,seq1b,seq1c,seq1d,seq1e,seq1f,seq1g))
seq2a = join(c(seq2a,seq2b,seq1c,seq1d,seq1e,seq1f,seq1g))
seq3a = join(c(seq3a,seq3b,seq1c,seq1d,seq1e,seq1f,seq1g))
seq2b = join(c(seq1a,seq2a,seq2b,seq1d,seq1e,seq1f,seq1g))
seq3b = join(c(seq1a,seq3a,seq3b,seq1d,seq1e,seq1f,seq1g))
seq2c = join(c(seq1a,seq1b,seq2a,seq2b,seq1e,seq1f,seq1g))
seq3c = join(c(seq1a,seq1b,seq3a,seq3b,seq1e,seq1f,seq1g))
seq2d = join(c(seq1a,seq1b,seq1c,seq2a,seq2b,seq1f,seq1g))
seq3d = join(c(seq1a,seq1b,seq1c,seq3a,seq3b,seq1f,seq1g))
seq2e = join(c(seq1a,seq1b,seq1c,seq1d,seq2a,seq2b,seq1g))
seq3e = join(c(seq1a,seq1b,seq1c,seq1d,seq3a,seq3b,seq1g))
# 
# #seq2 = join(c(seq1f,seq1b,seq1c,seq1d,seq1e))
# #seqT = join(c(seq1a,seq1f,seq1c,seq1d,seq1e))
# #seqG = join(c(seq1a,seq1b,seq1f,seq1d,seq1e))
# seq3 = join(c(seq1a,seq1b,seq1f,seq1g,seq1e,seq1i,seq1h))
# seq4 = join(c(seq1a,seq1b,seq1c,seq1f,seq1g,seq1i,seq1h))
# seq5 = join(c(seq1a,seq1b,seq1c,seq1d,seq1f,seq1g,seq1h))
# 
# seq2c = join(c(seq1a,seq1b,seq1c,seq1d,seq1e,seq1f,seq1g))
# seq3c = join(c(seq1a,seq1b,seq1c,seq1d,seq1e,seq1f,seq1g))

myseqs = c(seq1,seq2a,seq3a,seq2b,seq3b,seq2c,seq3c,seq2d,seq3d,seq2e,seq3e)
seqtypes = c('seq1','seq2a','seq3a','seq2b','seq3b','seq2c','seq3c','seq2d','seq3d','seq2e','seq3e')
dg = file.parse_deltaG()
dg$chunk = pasta(dg$n1,dg$n2);size(dg)

dg.sigma = data.frame()
for (i in 1:length(sigmas)) {
  sigma = sigmas[i]
  sigmaname = as.character(sigma)
  #Calc dg.sigma varied by m.all, without sum(B(n+1))
  alpha = N*sigma*A
  dg.sigma.precalc = a + (2*pi^2*C*K) / (4*pi^2*C + K*m.all) * (alpha+m.all*A)^2
  dg.sigma.precalc[1] = dg.sigma.precalc[1] - a #no rloop means no junction energy (a)
  dg.sigmatemp = data.frame(index=seq(1,length(m.all)),m=m.all)
  dg.sigmatemp[,sigmaname] = dg.sigma.precalc
  if (i == 1) {
    dg.sigma = dg.sigmatemp
  } else {
    dg.sigma[,sigmaname] = dg.sigma.precalc
  }
}

for (i in 1:length(myseqs)) {
  seq1 = myseqs[i]
  seqtype = seqtypes[i]
  window = 200
  print('fa.get_chunk')
  dfseq1 = fa.get_chunk(seq1,window=window);size(dfseq1);
  dfseq1$i0 = dfseq1$i0 + window/2
  print('dg.merge')
  dgseq1 = subset(dg.merge(dg,seq1),select=c(i0,DD,RD,RD.min.DD));size(dgseq1)
  
  

  print('dg.calc')
  dg.all1 = dg.calc(n.all,m.all,dgseq1,dg.sigma,seqtype=seqtype)
  dg.all1$dgmin00 = dg.all1$'0' + dg.all1$RD.min.DD
  dg.all1$dgmin07 = dg.all1$'-0.07' + dg.all1$RD.min.DD
  dg.all1$dgmin14 = dg.all1$'-0.14' + dg.all1$RD.min.DD
  print('save')
  save(dg.all1,file=pasta(seqtype,'_',window,'_.RDS'))
}
files = c()
for (i in 1:length(myseqs)) {
  seq1 = myseqs[i]
  seqtype = seqtypes[i]
  window = 200
  files = c(files,pasta(seqtype,'_',window,'_.RDS'))
}

load(files[1]); dg.seq1.200 = dg.all1
load(files[2]); dg.seq2a.200 = dg.all1
load(files[3]); dg.seq3a.200 = dg.all1
load(files[4]); dg.seq2b.200 = dg.all1
load(files[5]); dg.seq3b.200 = dg.all1
load(files[6]); dg.seq2c.200 = dg.all1
load(files[7]); dg.seq3c.200 = dg.all1
load(files[8]); dg.seq2d.200 = dg.all1
load(files[9]); dg.seq3d.200 = dg.all1
load(files[10]); dg.seq2e.200 = dg.all1
load(files[11]); dg.seq3e.200 = dg.all1
lay(4)
val.all = dg.seq1.200[dg.seq1.200$m == 200,]$dgmin07[1:1200]
val.all[is.na(val.all)] = 0
plot(NA,xlim=c(0,1200),ylim=c(min(val.all),max(val.all)))
lines(dg.seq1.200[dg.seq1.200$m == 200,]$dgmin07[1:1200],col='red',lwd=2)
lines(dg.seq2a.200[dg.seq2a.200$m == 200,]$dgmin07[1:1200],col='purple4',lwd=2)
lines(dg.seq3a.200[dg.seq3a.200$m == 200,]$dgmin07[1:1200],col='orange',lwd=2)
lines(dg.seq2b.200[dg.seq2b.200$m == 200,]$dgmin07[1:1200],col='green4',lwd=2)
lines(dg.seq3b.200[dg.seq3b.200$m == 200,]$dgmin07[1:1200],col='cornflowerblue',lwd=2)
segments(0,0,N1,0,lty=2)
plot(NA,xlim=c(0,1200),ylim=c(min(val.all),max(val.all)))
lines(dg.seq2c.200[dg.seq2c.200$m == 200,]$dgmin07[1:1200],col='blue1',lwd=2)
lines(dg.seq3c.200[dg.seq3c.200$m == 200,]$dgmin07[1:1200],col='blue2',lwd=2)
lines(dg.seq2d.200[dg.seq2d.200$m == 200,]$dgmin07[1:1200],col='blue3',lwd=2)
lines(dg.seq3d.200[dg.seq3d.200$m == 200,]$dgmin07[1:1200],col='blue4',lwd=2)
lines(dg.seq2e.200[dg.seq2e.200$m == 200,]$dgmin07[1:1200],col='purple1',lwd=2)
lines(dg.seq3e.200[dg.seq3e.200$m == 200,]$dgmin07[1:1200],col='purple2',lwd=2)
segments(0,0,N1,0,lty=2)


plot(NA,xlim=c(0,1000),ylim=c(min(val.all),max(val.all)))
lines(dg.all[dg.all$seqtype == '1.T_homopolymer' & dg.all$m == 200,]$dgmin07[1:1000],col='red',lwd=2)
# lines(dg.all[dg.all$seqtype == '2.Homogenized_random' & dg.all$m == 200,]$dgmin07[1:1000],col='purple4


window=200
dfseqT = fa.get_chunk(seqT,window=window)
dfseqT$i0 = dfseqT$i0 + window/2
dfseq1 = fa.get_chunk(seq1,window=window);size(dfseq1);
dfseq1$i0 = dfseq1$i0 + window/2
dfseq2 = fa.get_chunk(seq2,window=window)
dfseq2$i0 = dfseq2$i0 + window/2
dfseqG = fa.get_chunk(seqG,window=window)
dfseqG$i0 = dfseqG$i0 + window/2
dfseq3 = fa.get_chunk(seq3,window=window);size(dfseq3);
dfseq3$i0 = dfseq3$i0 + window/2
dfseq4 = fa.get_chunk(seq4,window=window);size(dfseq4);
dfseq4$i0 = dfseq4$i0 + window/2
dfseq5 = fa.get_chunk(seq5,window=window);size(dfseq5);
dfseq5$i0 = dfseq5$i0 + window/2

#je = -10


dg.merge = function(dg,seq1) {
  dgseq1 = fa.get_chunk(seq1,window=1);size(dgseq1)
  dgseq1 = merge(dgseq1,subset(dg[dg$type == 'DD',],select=c(chunk,dg)),by='chunk');size(dgseq1);colnames(dgseq1)[dim(dgseq1)[2]] = 'DD'
  dgseq1 = merge(dgseq1,subset(dg[dg$type == 'RD',],select=c(chunk,dg)),by='chunk');size(dgseq1);colnames(dgseq1)[dim(dgseq1)[2]] = 'RD'
  dgseq1 = merge(dgseq1,subset(dg[dg$type == 'RD.min.DD',],select=c(chunk,dg)),by='chunk');size(dgseq1);colnames(dgseq1)[dim(dgseq1)[2]] = 'RD.min.DD'
  dgseq1 = dgseq1[order(dgseq1$i0,dgseq1$i1,dgseq1$chunk),]
  return(dgseq1)
}

dg = file.parse_deltaG()
dg$chunk = pasta(dg$n1,dg$n2);size(dg)

dgseq1 = subset(dg.merge(dg,seq1),select=c(i0,DD,RD,RD.min.DD));size(dgseq1)
dgseq2 = subset(dg.merge(dg,seq2),select=c(i0,DD,RD,RD.min.DD));size(dgseq2)
dgseqT = subset(dg.merge(dg,seqT),select=c(i0,DD,RD,RD.min.DD));size(dgseqT)
dgseqG = subset(dg.merge(dg,seqG),select=c(i0,DD,RD,RD.min.DD));size(dgseqG)
dgseq3 = subset(dg.merge(dg,seq3),select=c(i0,DD,RD,RD.min.DD));size(dgseq3)
dgseq4 = subset(dg.merge(dg,seq4),select=c(i0,DD,RD,RD.min.DD));size(dgseq4)
dgseq5 = subset(dg.merge(dg,seq5),select=c(i0,DD,RD,RD.min.DD));size(dgseq5)
#dmseq1 = merge(dfseq1,dgseq1,by=i0);size(dfseq1)

A = 1/10.4
C = 1.8
N = 1500
R = 0.0019858775
T = 310
K = (2200 * R * T) / N
a = 10

#position = n
#m = length of rloop
# dgm = function(n,m.all,dgseq,sigma) {
#   print(n)
#   dgm.res = sapply(m.all,dgn,n=n,dgseq=dgseq,sigma=sigma)
#   return(dgm.res)
# }
# dgn = function(m,dgseq,sigma) {
#   b2 = (2*pi^2*C*K) / (4*pi^2*C + K*m) * (alpha + m*A)^2
#   dgn.res = b1 + b2 + b3
#   #cat(b1,b2,b3,'\n')
#   return(dgn.res)
# }
dgb = function(x,dgseq) {
  b3 = sum(dgseq[dgseq$i0 >= x[1] & dgseq$i0 <= x[1]+x[2],]$RD.min.DD)
  return(b3)
}
#n+1
#n+m
N1 = nchar(seq1)
window = 200
n.all = seq(1,1000)
m.all = seq(0,1000)
A = 1/10.4
C = 1.8
N = 1500
R = 0.0019858775
T = 310
K = (2200 * R * T) / N
a = 10
sigmas = c(0,-0.07,-0.14)

dg = data.frame()
for (i in 1:length(sigmas)) {
  sigma = sigmas[i]
  sigmaname = as.character(sigma)
  #Calc dG varied by m.all, without sum(B(n+1))
  alpha = N*sigma*A
  dg.precalc = a + (2*pi^2*C*K) / (4*pi^2*C + K*m.all) * (alpha+m.all*A)^2
  dg.precalc[1] = dg.precalc[1] - a #no rloop means no junction energy (a)
  dgtemp = data.frame(index=seq(1,length(m.all)),m=m.all)
  dgtemp[,sigmaname] = dg.precalc
  if (i == 1) {
    dg = dgtemp
  } else {
    dg[,sigmaname] = dg.precalc
  }
}

dg.calc = function(n.all,m.all,dgseq,dg,seqtype){
  dg.all = data.frame()
  for (i in 1:length(n.all)) {
    dg.all.temp = dg
    dg.all.temp$n = n.all[i]
    dg.all = rbind(dg.all,dg.all.temp)
  }
  dg.all$RD.min.DD = apply(subset(dg.all,select=c(n,m)),1,dgb,dgseq=dgseq)
  dg.all$seqtype = seqtype
  return(dg.all)
}
dg.allT = dg.calc(n.all,m.all,dgseqT,dg,seqtype='1.T_homopolymer')
dg.all1 = dg.calc(n.all,m.all,dgseq1,dg,seqtype='2.Homogenized_random')
dg.all2 = dg.calc(n.all,m.all,dgseq2,dg,seqtype='3.Homogenized_favorable')
dg.allG = dg.calc(n.all,m.all,dgseqG,dg,seqtype='4.G_homopolymer')
dg.all3 = dg.calc(n.all,m.all,dgseq3,dg,seqtype='5.Homogenized_mid0good')
dg.all4 = dg.calc(n.all,m.all,dgseq4,dg,seqtype='6.Homogenized_mid1good')
dg.all5 = dg.calc(n.all,m.all,dgseq5,dg,seqtype='7.Homogenized_lategood')
dg.all = dg.allT
dg.all = rbind(dg.all,dg.all1)
dg.all = rbind(dg.all,dg.all2)
dg.all = rbind(dg.all,dg.allG)
dg.all = rbind(dg.all,dg.all3)
dg.all = rbind(dg.all,dg.all4)
dg.all = rbind(dg.all,dg.all5)
#saveRDS(dg.all,file='dg.all3.RDS')
dg.all = readRDS('dg.all3.RDS')
dg.all$dgmin00 = dg.all$'0' + dg.all$RD.min.DD
dg.all$dgmin07 = dg.all$'-0.07' + dg.all$RD.min.DD
dg.all$dgmin14 = dg.all$'-0.14' + dg.all$RD.min.DD
reset()
val.all = dg.all$dgmin07
plot(NA,xlim=c(0,1000),ylim=c(min(val.all),max(val.all)))
lines(dg.all[dg.all$seqtype == '1.T_homopolymer' & dg.all$m == 200,]$dgmin07[1:1000],col='red',lwd=2)
lines(dg.all[dg.all$seqtype == '2.Homogenized_random' & dg.all$m == 200,]$dgmin07[1:1000],col='purple4',lwd=2)
lines(dg.all[dg.all$seqtype == '3.Homogenized_favorable' & dg.all$m == 200,]$dgmin07[1:1000],col='orange',lwd=2)
lines(dg.all[dg.all$seqtype == '4.G_homopolymer' & dg.all$m == 200,]$dgmin07[1:1000],col='green4',lwd=2)
lines(dg.all[dg.all$seqtype == '5.Homogenized_mid0good' & dg.all$m == 200,]$dgmin07[1:1000],col='blue4',lwd=3)
lines(dg.all[dg.all$seqtype == '6.Homogenized_mid1good' & dg.all$m == 200,]$dgmin07[1:1000],col='blue3',lwd=3)
lines(dg.all[dg.all$seqtype == '7.Homogenized_lategood' & dg.all$m == 200,]$dgmin07[1:1000],col='cornflowerblue',lwd=3)
segments(0,0,N1,0,lty=2)


plot(NA,xlim=c(0,1000),ylim=c(min(val.all),max(val.all)))
lines(dg.all[dg.all$seqtype == '1.T_homopolymer' & dg.all$m == 200,]$dgmin07[1:1000],col='red',lwd=2)
lines(dg.all[dg.all$seqtype == '2.Homogenized_random' & dg.all$m == 200,]$dgmin07[1:1000],col='purple4',lwd=2)
lines(dg.all[dg.all$seqtype == '3.Homogenized_favorable' & dg.all$m == 200,]$dgmin07[1:1000],col='orange',lwd=2)
lines(dg.all[dg.all$seqtype == '4.G_homopolymer' & dg.all$m == 200,]$dgmin07[1:1000],col='green4',lwd=2)
lines(dg.all[dg.all$seqtype == '5.Homogenized_mid0good' & dg.all$m == 200,]$dgmin07[1:1000],col='blue4',lwd=3)
lines(dg.all[dg.all$seqtype == '6.Homogenized_mid1good' & dg.all$m == 200,]$dgmin07[1:1000],col='blue3',lwd=3)
lines(dg.all[dg.all$seqtype == '7.Homogenized_lategood' & dg.all$m == 200,]$dgmin07[1:1000],col='cornflowerblue',lwd=3)

lay(n1=3,n2=2)
plotthis(dfseqT[dfseqT$i0 <= 1500,],dg.all[dg.all$seqtype == '1.T_homopolymer'         & dg.all$m == 200,]$dgmin07/300,window=200)
plotthis(dfseq1[dfseq1$i0 <= 1500,],dg.all[dg.all$seqtype == '2.Homogenized_random'    & dg.all$m == 500,]$dgmin07/300,window=200)
plotthis(dfseq2[dfseq2$i0 <= 1500,],dg.all[dg.all$seqtype == '3.Homogenized_favorable' & dg.all$m == 500,]$dgmin07/300,window=200)
plotthis(dfseq3[dfseq3$i0 <= 1500,],dg.all[dg.all$seqtype == '5.Homogenized_mid0good'  & dg.all$m == 500,]$dgmin07/300,window=200)
plotthis(dfseq4[dfseq4$i0 <= 1500,],dg.all[dg.all$seqtype == '6.Homogenized_mid1good'  & dg.all$m == 500,]$dgmin07/300,window=200)
plotthis(dfseq5[dfseq5$i0 <= 1500,],dg.all[dg.all$seqtype == '7.Homogenized_lategood'  & dg.all$m == 500,]$dgmin07/300,window=200)
plotthis(dfseqG[dfseqG$i0 <= 1500,],dg.all[dg.all$seqtype == '4.G_homopolymer'         & dg.all$m == 500,]$dgmin07/300,window=200)

plotthis(dfseqT,(dg.allT[dg.allT$m == window,]$dg[1:800]-25)/40,x0=101,x1=900,window=window)
plotthis(dfseq1,(dg.all1[dg.all1$m == window,]$dg[1:800]-25)/40,x0=101,x1=900,window=window)
plotthis(dfseq2,(dg.all2[dg.all2$m == window,]$dg[1:800]-25)/40,x0=101,x1=900,window=window)
plotthis(dfseqT,dg.allG[dg.allG$m == window,]$dg,window=window)

nom = e^(-1*dgm1.n/(R*T))
denom = e^(-1*dgm1.all/(R*T))


#dg.all = dg.all1[dg.all1$n == 1,]$RD.min.DD
plotthis = function(dfseq,dg.all,x0=NA,x1=NA,window) {
  if (is.na(x0)) {
    dg.all.x = seq(window/2,window/2+length(dg.all)-1)
  } else {
    dg.all.x = seq(x0,x1)
  }
  #dg.all.x = seq(window/2,window/2+length(dg.all)-1)
  plot(dg.all.x,dg.all,ylim=c(0,1.5),xlim=c(0,max(dfseq$i0)+window/2),type='l',col='purple')
  lines(dfseq$i0,dfseq$CGdens,col='blue4')
  lines(dfseq$i0,dfseq$GCcont,col='green4')
  lines(dfseq$i0,dfseq$GCskew+0.5,col='red4')
#  lines(seq(window/2,window/2+length(dg.all)-1),(dg.all-50)/20,col='purple')
}
plot(NA,ylim=c(0,1.5),xlim=c(0,size(dfseq2)))
lines(dfseq2$i0,dfseq2$CGdens,col='blue4')
lines(dfseq2$i0,dfseq2$GCcont,col='green4')
lines(dfseq2$i0,dfseq2$GCskew+0.5,col='red4')

plot(NA,ylim=c(0,1.5),xlim=c(0,size(dfseqT)))
lines(dfseqT$i0,dfseqT$CGdens,col='blue4')
lines(dfseqT$i0,dfseqT$GCcont,col='green4')
lines(dfseqT$i0,dfseqT$GCskew+0.5,col='red4')

plot(NA,ylim=c(0,1.5),xlim=c(0,size(dfseqG)))
lines(dfseqG$i0,dfseqG$CGdens,col='blue4')
lines(dfseqG$i0,dfseqG$GCcont,col='green4')
lines(dfseqG$i0,dfseqG$GCskew+0.5,col='red4')

plot(NA,ylim=c(0,1.5),xlim=c(0,size(dfseq5)))
lines(dfseq5$i0,dfseq5$CGdens/3,col='blue4')
lines(dfseq5$i0,dfseq5$GCcont,col='green4')
lines(dfseq5$i0,dfseq5$GCskew+0.5,col='red4')

plot(NA,ylim=c(0,1.5),xlim=c(0,size(dfseq6)))
lines(dfseq6$i0,dfseq6$CGdens/3,col='blue4')
lines(dfseq6$i0,dfseq6$GCcont,col='green4')
lines(dfseq6$i0,dfseq6$GCskew+0.5,col='red4')
