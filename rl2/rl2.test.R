options(max.print = 200)
source("lib/tri.lib.main.R")
source('rl2/lib/rl2.lib9.R')
source('rl2/lib/rl2.lib1.R')
source('lib/tri.lib.misc.R')
#Example
seqr = fa.gen_random(n=4000)
seqT = seqr$seq$seqT
seqG = seqr$seq$seqG
seqgood = seqr$seq$seqgood
seqok = seqr$seq$seqok
# seqr.seqG.CpGprof = fa.get_CpGprof(seqr$seq$seqG,gp=gp,print=TRUE)
# seqr.seqT.CpGprof = fa.get_CpGprof(seqr$seq$seqT,gp=gp,print=TRUE)
# seqr.2c.CpGprof = fa.get_CpGprof(seqr$seq$seq2c,gp=gp,print=TRUE)
# seqr.3c.CpGprof = fa.get_CpGprof(seqr$seq$seq3c,gp=gp,print=TRUE)

cons = rl2.get_constants()
FASTAS.rl2 = rl2.fa.get_FASTAS(force=T)

par = list(
#  n = seq(1,3000),
  m = seq(0,10000),
  sigmas = seq(0,-0.14,by=-0.01),
  window = 200,
  invivo_buffer = 1990,
  invitro_buffer = 10
)
par$len.max = max(par$m)

seqtypeswant = c("pFC9","T7_init_VR_10","T7_init_VR_20","T7_init_VR_1","T7_init_VR_4","T7_init_VR_17","T7_init_VR_18","T7_init_VR_21","T7_init_VR_25","T7_init_VR_31","pFC53_ApaLI","pFC53TVR_10_pair_0_T3TermMix","pFC53TVR_20_pair_0_T3TermMix",'CALM3','BTF3','FUS','MRFAP1L1','CNBP')

FASTASwant = FASTAS.rl2[FASTAS.rl2$gene %in% seqtypeswant,]
myseqs = c(FASTASwant$seq,seqT,seqG)
seqtypes = c(FASTASwant$gene,"seqT","seqG")

length(seqtypes)
length(myseqs)
#myseqs = c(seqVR10,seqVR20,seqT,seqG,seqgood,seqok,seq2a)
#seqtypes = c('T7_init_VR_10','T7_init_VR_20','seqT','seqG','seqgood','seqok','seq2a')
#'seqT','seqG','seq1','seq2a','seq3a','seq2b','seq3b','seq2c','seq3c','seq2d','seq3d','seq2e','seq3e')

#dG = sum(B) + (a + dG.sigma(m,sigma))
B = file.parse_deltaG()
B2 = file.parse_RNAG()
print(paste("- Using max R-loop length (par$m) =",max(par$m)))
dgas = data.frame()
dgas = rbind(dgas,dg.sigma.precalc(m.max=max(par$m),a=10,sigma=-0.00,cons=cons,print=T))
dgas = rbind(dgas,dg.sigma.precalc(m.max=max(par$m),a=10,sigma=-0.01,cons=cons,print=T))
dgas = rbind(dgas,dg.sigma.precalc(m.max=max(par$m),a=10,sigma=-0.04,cons=cons,print=T))
dgas = rbind(dgas,dg.sigma.precalc(m.max=max(par$m),a=10,sigma=-0.07,cons=cons,print=T))
dgas = rbind(dgas,dg.sigma.precalc(m.max=max(par$m),a=10,sigma=-0.11,cons=cons,print=T))
dgas = rbind(dgas,dg.sigma.precalc(m.max=max(par$m),a=10,sigma=-0.14,cons=cons,print=T))
print(dim(dgas)) #10000 x 6



source('rl2/lib/rl2.lib1.R')


dg.all = data.frame()
dgseq.all = data.frame()
dfseq.all = data.frame()
for (i in 1:length(myseqs)) {
  seq1 = myseqs[i]
  N1 = nchar(seq1)
  seqtype = seqtypes[i]
  print(paste('\n',i,'/',length(myseqs),': Doing',seqtypes[i]))
  
  file.dgseq = pasta('resources/',seqtype,'_',par$window,'_dgseq.RDS')
  file.dgall = pasta('resources/',seqtype,'_',par$window,'_dgall.RDS')
  # if (!file.exists(file)) {
  #   next
  # }
  
  print(paste('- rl2.fa.get_chunk',seqtype))
  if (file.exists(file.dfseq)) {
    dfseq1 = readRDS(file.dfseq)
    #dfseq.all = rbind(dfseq.all,dfseq1)
  } else {
    
    dfseq1 = rl2.fa.get_chunk(seq1,window=par$window,xlim=c(0,nchar(seq1)),print=F,my.title=seqtype);size(dfseq1);
    rownames(dfseq1) = dfseq1$i0
    dfseq1$n = dfseq1$i0
    dfseq1$fileInd = i
    dfseq1$seqtype = seqtype
    dfseq1$n = dfseq1$i0

    print(paste('- dg.merge',seqtype))
    dgseq1 = subset(dg.merge(B,seq1),select=c(i0,DD,RD,RD.min.DD));size(dgseq1)
    dgseq1 = merge(dfseq1,dgseq1,by='i0')
    dgseq1$B = dgseq1$RD.min.DD
    dgseq1 = subset(dgseq1,select=-chunk)
    dgseq1 = dgseq1[order(dgseq1$n),]
    saveRDS(dgseq1,file=file.dfseq)
    print('  -> save')
  }  
  ggplot(dgseq1,aes(n,B)) +
    geom_line(aes(y=smooth.spline(-1 * B + 1,spar=0.4)$y),color='purple4') +
    geom_line(aes(y=CGdens),color='cornflowerblue') +
    geom_line(aes(y=GCcont),color='green4') +
    theme_bw() + theme(panel.grid = element_blank()) +
    geom_line(aes(y=GCskew+0.5),color='red2') +
    ggtitle(paste(seqtype,'CpG Profile')) +
    coord_cartesian(xlim=c(0,max(dgseq1$n)), ylim=c(0,1.4)) +
    xlab("Position in plasmid (bp)") +
    ylab("CpG density/GC content/(GC skew +0.5)") +
    annotate(geom='segment',x=0,xend=max(dgseq1$n),y=0.5,yend=0.5,lty=2)
    
  # 
  # dgseq1$B.sum = apply(data.frame(n=dgseq1$n,m=max(dgseq1$n)),1,dgb,dgseq=dgseq1$B)
  # dgseq1$m = max(dgseq1$n) - dgseq1$n + 1
  # mwant = 500
  # nwant = 10
  # Bn = function(x, dfseq) {dfseq[dfseq$n == x[1],]$B.sum - dfseq[dfseq$n == x[1]+x[2],]$B.sum}
  # Bm = function(m, n, dgseq1) {
  #   sapply(n,FUN = Bn,m=m,dgseq1=dgseq1)
  # }
  # test5 = sapply(seq(1,mwant),Bm,n=seq(1,nwant),dgseq1=dgseq1)
  print(paste('- dg.calc',seqtype))
  if (file.exists(file.dgB)) {
    print(paste('- ',file.dgB,'exist!'))
  } else {

    par.n = seq(1,min(nchar(seq1),par$len.max))
    par.m = par$m[par$m <= nchar(seq1)]
    dgB1 = dgB.calcs(par.n,par.m,dgseq1$B)
    dgB.calc = function(par.n,par.m,dgseqB) {
      par.n = seq(1,i)
      dgB1 = as.data.frame(expand.grid(n=par.n,m=par.m))
      dgB1 = dgB1[order(dgB1$n,dgB1$m),]
      start.time <- Sys.time()
      dgB1$B = apply(dgB1,1,dgB1,dg=dgseqB)
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      print(time.taken)  
      return(dgB1)
    }    
    saveRDS(dgB1,file=file.dgB)
    print(paste('- ',file.dgB,'saved!'))
  }
}







#pFC9

# dG final = a + dgall$'0' + dgmin00 + 
seqwant = c('T7_init_VR_10','T7_init_VR_20','pFC9')

dg.all = data.frame()
dgseq.all = data.frame()
dfseq.all = data.frame()
for (i in seq(1,length(seqwant))) {
  seq1 = myseqs[seqtypes == seqwant[i]]
  N1 = nchar(seq1)
  seqtype = seqtypes[seqtypes == seqwant[i]]
  print(paste('\n',i,'/',length(seqwant),': Doing',seqtype))
  
  file.dfseq = pasta('resources/',seqtype,'_',par$window,'_dfseq.RDS')
  dfseq1 = readRDS(file.dfseq); dfseq1$seqtype = seqtype
  dfseq.all = rbind(dfseq.all,dfseq1)
  file.dgseq = pasta('resources/',seqtype,'_',par$window,'_dgseq.RDS')
  dgseq1 = readRDS(file.dgseq); dgseq1$seqtype = seqtype
  dgseq.all = rbind(dgseq.all,dgseq1)
  file.dg = pasta('resources/',seqtype,'_',par$window,'_dgall.RDS')
  dg1 = readRDS(file.dg); dg1$seqtype = seqtype
  dg.all = rbind(dg.all,dg1)
}
saveRDS(dfseq.all,file='resources/dfseq.all.pFC9.RDS')
saveRDS(dgseq.all,file='resources/dgseq.all.pFC9.RDS')
saveRDS(dg.all,file='resources/dg.all.pFC9.RDS')




dfseq.all = readRDS(file='resources/dfseq.all.pFC9.RDS')
dgseq.all = readRDS(file='resources/dgseq.all.pFC9.RDS')
dg.all    = readRDS(file='resources/dg.all.pFC9.RDS')

dg.all.m = melt(dg.all[,c(-3,-4,-5)],id.vars=colnames(dg.all)[c(1,2,6,7,11)])
#dg.all.m = melt(dg.all,id.vars = colnames(dg.all)[c(1,2,6,11)])
dg.all.m$variable = as.character(dg.all.m$variable)
# 
# dg.all.m0 = dg.all.m
# dg.all.m0$a = 0
# 
# dg.all.m.a10 = dg.all.m
# dg.all.m.a10$a = 10
# 
# dg.all.m.a25 = dg.all.m
# dg.all.m.a25$a = 25
# 
# dg.all.m.a50 = dg.all.m
# dg.all.m.a50$a = 50
# 
# dg.all.m.a100 = dg.all.m
# dg.all.m.a100$a = 100
# 
# dg.all.m.final = rbind(dg.all.m0,dg.all.m.a10,dg.all.m.a25,dg.all.m.a50,dg.all.m.a100)
# remove(dg.all.m0)
# 
# remove(dg.all.m10)
# 
# remove(dg.all.m25)
# 
# remove(dg.all.m50)
# 
# remove(dg.all.m100)
# 



#temp = dg.all.m[grep('dgmin',dg.all.m$variable),]




# temp = temp[temp$n <= 2000 & temp$n %% 100 == 0,]
# temp = temp[temp$seqtype != 'seq2a',]
# temp0 = temp0[temp0$seqtype != 'seq2a',]
#temp0 = temp0[temp0$seqtype != 'T7_init_VR_10',]
#temp0 = temp0[temp0$seqtype != 'T7_init_VR_20',]
nwant = 600; 
juncEs = c(0,10,25,50,100)
pdf("variable_vs_juncEs.pdf",width=15,height=5)
for (i in 1:length(juncEs)) {
  temp0 = dg.all.m[dg.all.m$n == nwant,]
  temp0[temp0$m == 0,]$value = temp0[temp0$m == 0,]$value + 10
  
  temp0[temp0$m != 0,]$value = temp0[temp0$m != 0,]$value + juncEs[i]
  p1 = ggplot(temp0,aes(m,value)) +
    geom_line(aes(color=af(seqtype))) + theme_bw() +
    annotate(geom='segment',x=0,xend=1000,y=53.8,yend=53.8,lty=2) +
    theme_bw() + theme(panel.grid=element_blank()) +
    ggtitle(paste('Energy of R-loop length at n =',nwant,'juncE (a)=',juncEs[i])) +
    #scale_color_manual(values=c('seqT' = 'red3','seqok'='purple2','seqgood'='orange3','seqG'='green4','T7_init_VR_20'='blue4')) +
    facet_grid(.~variable) +ylab('Energy (kcal/mol)') + xlab("Nucleotides in R-loop")
  print(p1)
  remove(temp0)
}
dev.off()




#temp$variable2 = gsub('dgmin(.+)','sigma=\\1',temp$variable,perl=T)
pdf("seqtype_vs_juncEs.pdf",width=30,height=5)
#for (seqtype in unique(temp$seqtype)) {
  for (i in 1:length(juncEs)) {
    temp0 = dg.all.m[dg.all.m$n == nwant,]# & dg.all.m$seqtype == seqtype,]
    temp0[temp0$m == 0,]$value = temp0[temp0$m == 0,]$value + 10
    temp0[temp0$m != 0,]$value = temp0[temp0$m != 0,]$value + juncEs[i]

    p2 = ggplot(temp0,aes(m,value)) +
      geom_line(aes(color=af(variable))) + theme_bw() +
      annotate(geom='segment',x=0,xend=1000,y=53.8,yend=53.8,lty=2) +
      annotate(geom='segment',x=0,xend=1000,y=190,yend=190,lty=2) +
      annotate(geom='segment',x=0,xend=1000,y=10.8,yend=10.8,lty=2) +
      theme_bw() + theme(panel.grid=element_blank()) +
      #scale_color_manual(values=c('brown4','orange','cornflowerblue')) +
      ggtitle(paste('Energy of R-loop length at n =',nwant,'juncE (a)=',juncEs[i])) +
      facet_grid(.~seqtype) +
      #'seqT' = 'red3','seqok'='purple3','seqgood'='orange3','seqG'='green4')) +
      ylab('Energy (kcal/mol)') + xlab("Nucleotides in R-loop")
    print(p2)
    remove(temp0)
  }
#}



mwant = 200
pdf("seqtype_vs_position_m200_juncEs.pdf",width=8*length(unique(dg.all$seqtype)),height=4)
#for (seqtype in unique(temp$seqtype)) {
for (i in 1:length(juncEs)) {
  temp0 = dg.all.m[dg.all.m$m == mwant,]# & dg.all.m$seqtype == seqtype,]
  temp0[temp0$m == 0,]$value = temp0[temp0$m == 0,]$value + 10
  temp0[temp0$m != 0,]$value = temp0[temp0$m != 0,]$value + juncEs[i]
  
  p3 = ggplot(dfseq.all,aes(x = dfseq.all$i0,temp4$CGdens)) +
    geom_line(aes(y=CGdens),color='cornflowerblue') +
    geom_line(aes(y=GCcont),color='green4') +
    geom_line(aes(y=GCskew+0.5),color='red4') +
    # geom_line(aes(y=smooth.spline(CGdens,spar=0.1)$y),color='cornflowerblue') +
    # geom_line(aes(y=smooth.spline(GCcont,spar=0.1)$y),color='green4') +
    # geom_line(aes(y=smooth.spline(GCskew+0.5,spar=0.1)$y),color='red4') +
    # ggtitle(paste('Energy of R-loop length at m =',mwant,'juncE (a)=',juncEs[i])) +
    theme_bw() + theme(panel.grid=element_blank()) +
    annotate(geom='rect',xmin=587-292+100,xmax=605-292+100,ymin=-10,ymax=10,fill='black',alpha=0.1,color=NA) +
    annotate(geom='rect',xmin=641-292+100,xmax=841-292+100,ymin=-10,ymax=10,fill='green4',alpha=0.1,color=NA) +
    annotate(geom='rect',xmin=842-292+100,xmax=1317-292+100,ymin=-10,ymax=10,fill='orange2',alpha=0.1,color=NA) +
    geom_line(data=temp0,aes(x=n,y=(value-par$dg.minval)/par$dg.maxval,color=af(variable))) +
    scale_color_manual(values=c('dgmin00'='purple2','dgmin07'='purple3','dgmin14'='purple4')) +
    #geom_line(data=temp0,aes(x=n,y=smooth.spline((value-par$dg.minval)/par$dg.maxval,spar=0.1)$y,color=af(variable))) +
    coord_cartesian(xlim=c(0,4000),ylim=c(0,1.4)) +
    facet_grid(.~seqtype) +
    ylab('Energy (kcal/mol)') + xlab("Posotion at plasmid (bp)") + annotate(geom='segment',x=0,xend=10000,y=0.5,yend=0.5,lty=2)
  # 
  # p2 = ggplot(temp0,aes(m,value)) +
  #   geom_line(aes(color=af(variable))) + theme_bw() +
  #   annotate(geom='segment',x=0,xend=1000,y=53.8,yend=53.8,lty=2) +
  #   annotate(geom='segment',x=0,xend=1000,y=190,yend=190,lty=2) +
  #   annotate(geom='segment',x=0,xend=1000,y=10.8,yend=10.8,lty=2) +
  #   theme_bw() + theme(panel.grid=element_blank()) +
  #   #scale_color_manual(values=c('brown4','orange','cornflowerblue')) +
  #   ggtitle(paste('Energy of R-loop length at n =',nwant,'juncE (a)=',juncEs[i])) +
  #   facet_grid(.~seqtype) +
  #   #'seqT' = 'red3','seqok'='purple3','seqgood'='orange3','seqG'='green4')) +
  #   ylab('Energy (kcal/mol)') + xlab("Nucleotides in R-loop")
  print(p3)
  remove(temp0)
}
dev.off()
#}

#   p2VR20 = ggplot(temp[temp$seqtype == 'T7_init_VR_20' & temp$n == 100,],aes(m,value)) +
#     geom_line(aes(color=af(variable2))) + theme_bw() +
#     annotate(geom='segment',x=0,xend=1000,y=53.8,yend=53.8,lty=2) +
#     annotate(geom='segment',x=0,xend=1000,y=190,yend=190,lty=2) +
#     annotate(geom='segment',x=0,xend=1000,y=10.8,yend=10.8,lty=2) +
#     theme_bw() + theme(panel.grid=element_blank()) +
#     scale_color_manual(values=c('brown4','orange','cornflowerblue')) +
#     #'seqT' = 'red3','seqok'='purple3','seqgood'='orange3','seqG'='green4')) +
#     ylab('Energy (kcal/mol)') + xlab("Nucleotides in R-loop")
#   
# p2VR10
# p2VR20







par$window = 200
#dg.all
par$dg.minval = -55#min(c(dg.all$dgmin07,dg.all$dgmin07,dg.all$dgmin07))
par$dg.maxval = 300#max(c(dg.all$dgmin07,dg.all$dgmin07,dg.all$dgmin07))
plot.this = function(dfseq.all,dg.all,window=200,sigma,myseqtype,m,m1=NA,par,a=0) {
  if (sigma == -0.07) {
    sigmawant = 'dgmin07'
  } else if (sigma == 0) {
    sigmawant = 'dgmin00'
  } else if (sigma == -0.14) {
    sigmawant = 'dgmin14'
  } else {
    sigmawant = 'dgmin07'
  }
  temp3 = dg.all #temp3 = dg.all[dg.all$seqtype == myseqtype,]
  temp3[temp3$m == 0,sigmawant] = temp3[temp3$m == 0,sigmawant] + 10
  temp3[temp3$m != 0,sigmawant] = temp3[temp3$m != 0,sigmawant] + a
  
  temp4 = dfseq.all[dfseq.all$seqtype == myseqtype,]
  p = ggplot(temp4,aes(x = temp4$i0+par$window/2,temp4$CGdens)) +
    geom_line(aes(y=smooth.spline(CGdens,spar=0.2)$y),color='cornflowerblue') +
    geom_line(aes(y=smooth.spline(GCcont,spar=0.2)$y),color='green4') +
    geom_line(aes(y=smooth.spline(GCskew+0.5,spar=0.2)$y),color='red4') +
    ggtitle(paste(myseqtype,sigma,m)) +
    theme_bw() + theme(panel.grid=element_blank()) +
    annotate(geom='rect',xmin=587-292,xmax=605-292,ymin=-10,ymax=10,fill='black',alpha=0.1,color=NA) +
    annotate(geom='rect',xmin=641-292,xmax=841-292,ymin=-10,ymax=10,fill='green4',alpha=0.1,color=NA) +
    annotate(geom='rect',xmin=842-292,xmax=1317-292,ymin=-10,ymax=10,fill='orange2',alpha=0.1,color=NA) +
    facet_grid(.~seqtype) +
    ylab('Energy (kcal/mol)') + xlab("Nucleotides in R-loop") + annotate(geom='segment',x=0,xend=10000,y=0.5,yend=0.5,lty=2)
  
  if (is.na(m1)) {
    temp3b = temp3
    if (defined(temp3b[temp3b$m == m,])) {
      temp3b = temp3b[temp3b$m == m,]
      p2 = p + geom_line(data=temp3b,aes(x=n+par$window/2,y=smooth.spline((temp3b[,sigmawant]-par$dg.minval)/par$dg.maxval)$y),color='purple4',lwd=2) # +    coord_cartesian(xlim=c(0,1000),ylim=c(0,1.4))
      return(p2)
      p3 = ggplot(temp3b,aes(x=n+par$window/2,y=(smooth.spline(temp3b[,sigmawant])$y-par$dg.minval)/par$dg.maxval,group=af(m),color=af(m))) +
        geom_line(aes(color=af(m)),lwd=1) +
        annotate(geom='rect',xmin=587-292,xmax=605-292,ymin=-10,ymax=10,fill='black',alpha=0.1,color=NA) +
        annotate(geom='rect',xmin=641-292,xmax=841-292,ymin=-10,ymax=10,fill='green4',alpha=0.1,color=NA) +
        annotate(geom='rect',xmin=842-292,xmax=1317-292,ymin=-10,ymax=10,fill='orange2',alpha=0.1,color=NA) +
        theme_bw() + theme(panel.grid=element_blank()) +
        ylab('Energy (kcal/mol)') + xlab("Nucleotides in R-loop") + annotate(geom='segment',x=0,xend=10000,y=0.5,yend=0.5,lty=2)
      
    }
  } else {
    temp3b = temp3
    if (defined(temp3b[temp3b$m >= m & temp3b$m <= m1,])) {
      temp3b = temp3b[temp3b$m >= m & temp3b$m <= m1,]
      if (defined(temp3b[temp3b$m %% 50 == 0,])) {
        temp3b = temp3b[temp3b$m %% 50 == 0,]
        print(unique(temp3b$m))
        p2 = p + geom_line(data=temp3b,aes(x=n+par$window/2,y=(temp3b[,sigmawant]-par$dg.minval)/par$dg.maxval,group=af(m),color=af(m)),lwd=0.5) # +    coord_cartesian(xlim=c(0,1000),ylim=c(0,1.4))
        return(p2)
        p3 = ggplot(temp3b,aes(x=n+par$window/2,y=(temp3b[,sigmawant]-par$dg.minval)/par$dg.maxval,group=af(m),color=af(m))) +
          geom_line(aes(color=af(m)),lwd=1) +
          theme_bw() + theme(panel.grid=element_blank()) +
          ylab('Energy (kcal/mol)') + xlab("Nucleotides in R-loop") + annotate(geom='segment',x=0,xend=10000,y=0.5,yend=0.5,lty=2)
        
      }
    }
    else {
      temp3b = temp3
      if (defined(temp3b[temp3b$m == m,])) {
        temp3b = temp3b[temp3b$m == m,]
        p2 = p + geom_line(data=temp3b,aes(x=n+par$window/2,y=smooth.spline((temp3b[,sigmawant]-par$dg.minval)/par$dg.maxval)$y),color='purple4',lwd=2) # +    coord_cartesian(xlim=c(0,1000),ylim=c(0,1.4))
        
        return(p2)
        p3 = ggplot(temp3b,aes(x=n+par$window/2,y=(smooth.spline(temp3b[,sigmawant])$y-par$dg.minval)/par$dg.maxval,group=af(m),color=af(m))) +
          geom_line(aes(color=af(m)),lwd=1) +
          theme_bw() + theme(panel.grid=element_blank()) +
          ylab('Energy (kcal/mol)') + xlab("Nucleotides in R-loop") + annotate(geom='segment',x=0,xend=10000,y=0.5,yend=0.5,lty=2)
        
      }
    }
  }
  return(p2)
  #grid.arrange(p2)# +       coord_cartesian(xlim=c(par$window,000),ylim=c(0,1.4)))#,p3)
  #p2 = ggplot(temp3,aes(x = n,dgmin07)) +
  #  geom_line(color='purple4')
  #print(p2)
}



temp3 = dg.all[dg.all$seqtype == 'T7_init_VR_20',]
temp3b = temp3b[temp3b$m >= 200 & temp3b$m <= 400,]
temp3b = temp3b[temp3b$m %% 50 == 0,]
ggplot(temp3b,aes(x=n+292+par$window/2,y=(temp3b$dgmin07-par$dg.minval)/par$dg.maxval)) +
  geom_line(aes(color=af(m)),lwd=1)






pdf("seqtype_vs_juncEs_CpGprof.pdf",width=30,height=30)
mwant = 200
for (i in 1:length(juncEs)) {
  temp3 = dg.all.m #temp3 = dg.all[dg.all$seqtype == myseqtype,]
  temp3[temp3$m == 0,] = temp3[temp3$m == 0,] + 10
  temp3[temp3$m != 0,] = temp3[temp3$m != 0,] + a
  
  temp4 = dfseq.all #temp4 = dfseq.all[dfseq.all$seqtype == myseqtype,]
  temp3b = temp3
  p = ggplot(temp4,aes(x = temp4$i0+par$window/2,temp4$CGdens)) +
    geom_line(aes(y=smooth.spline(CGdens,spar=0.2)$y),color='cornflowerblue') +
    geom_line(aes(y=smooth.spline(GCcont,spar=0.2)$y),color='green4') +
    geom_line(aes(y=smooth.spline(GCskew+0.5,spar=0.2)$y),color='red4') +
    ggtitle(paste(myseqtype,sigma,m)) +
    theme_bw() + theme(panel.grid=element_blank()) +
    annotate(geom='rect',xmin=587-292,xmax=605-292,ymin=-10,ymax=10,fill='black',alpha=0.1,color=NA) +
    annotate(geom='rect',xmin=641-292,xmax=841-292,ymin=-10,ymax=10,fill='green4',alpha=0.1,color=NA) +
    annotate(geom='rect',xmin=842-292,xmax=1317-292,ymin=-10,ymax=10,fill='orange2',alpha=0.1,color=NA) +
    facet_grid(.~seqtype) +
    ylab('Energy (kcal/mol)') + xlab("Nucleotides in R-loop") + annotate(geom='segment',x=0,xend=10000,y=0.5,yend=0.5,lty=2)
  
  if (defined(temp3b[temp3b$m == mwant,])) {
    temp3b = temp3b[temp3b$m == mwant,]
    p2 = p + geom_line(data=temp3b,aes(x=n+par$window/2,y=smooth.spline((temp3b[,sigmawant]-par$dg.minval)/par$dg.maxval)$y),color='purple4',lwd=2) +    coord_cartesian(xlim=c(0,4000),ylim=c(0,1.4))
    print(p2)
  }
  
  # p2 = plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype=seqtype,m=200,par=par,a=juncEs[i])
  # print(p2 + coord_cartesian(xlim=c(0,4000),ylim=c(0,1.4)))
}
dev.off()
# plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='seqT',m=200,par=par)
# 
# plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='seqok',m=200,par=par)
# 
# 
# 
# plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='T7_init_VR_20',m=200,par=par)
# 
# plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='seqgood',m=200,par=par)
# 
# plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='seqG',m=200,par=par)




plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='seqT',m=200,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='seqok',m=200,par=par)


plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='T7_init_VR_10',m=200,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='T7_init_VR_20',m=200,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='seqgood',m=200,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='seqG',m=200,par=par)






plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='seqT',m=200,m1=400,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='seqok',m=200,m1=400,par=par)


plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='T7_init_VR_10',m=200,m1=400,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='T7_init_VR_20',m=200,m1=400,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='seqgood',m=200,m1=400,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='seqG',m=200,m1=400,par=par)




plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='seqT',m=200,m1=400,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='seqok',m=200,m1=400,par=par)


plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='T7_init_VR_10',m=200,m1=400,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='T7_init_VR_20',m=200,m1=400,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='seqgood',m=200,m1=400,par=par)

plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='seqG',m=200,m1=400,par=par)





plot.this(dfseq.all,dg.all,window=200,sigma= 0   ,myseqtype='T7_init_VR_10',m=200,m1=400,par=par)
plot.this(dfseq.all,dg.all,window=200,sigma= 0   ,myseqtype='T7_init_VR_20',m=200,m1=400,par=par)
plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='T7_init_VR_10',m=200,m1=400,par=par)
plot.this(dfseq.all,dg.all,window=200,sigma=-0.07,myseqtype='T7_init_VR_20',m=200,m1=400,par=par)
plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='T7_init_VR_10',m=200,m1=400,par=par)
plot.this(dfseq.all,dg.all,window=200,sigma=-0.14,myseqtype='T7_init_VR_20',m=200,m1=400,par=par)





dg.calc_prob = function(dfseq.all,dg.all,window=200,sigma,myseqtype,m) {

temp4 = dfseq.all[dfseq.all$i0 <= 2000 & dfseq.all$seqtype == myseqtype,]
temp3 = dg.all[dg.all$seqtype == myseqtype & dg.all$m == 200,]
exp(-1*dg.seq2c.200$dgmin00/(R*T))/sum(exp(-1*dg.seq2c.200$dgmin00/(R*T)))
}
# # # # 
# # # # pdf("Fig1.pdf")
# # # # ggplot(temp,aes(m,value)) +
# # # #   geom_line(aes(color=af(variable))) + theme_bw() +
# # # #   facet_grid(af(n)~seqtype)
# # # # ggplot(temp,aes(m,value)) +
# # # #   geom_line(aes(color=af(seqtype))) + theme_bw() +
# # # #   facet_grid(af(n)~variable)
# # # # dev.off()
# # # # 
# # # # 
# # # 
# # # for (i in 1:length(myseqs)) {
# # #   par$window = 200
# # #   seq1 = myseqs[i]
# # #   N1   = nchar(seq1)
# # #   seqtype = seqtypes[i]
# # #   fileInd = i
# # #   
# # #   dgseq
# # #   file = pasta(seqtype,'_',par$window,'_.RDS')
# # #   dfseq = dfseq
# # # }
# # # lay(4)
# # # reset()
# # # val.all = dg.seq1.200[dg.seq1.200$m == 200,]$dgmin07[1:1200]
# # # val.all[is.na(val.all)] = 0
# # # plot(NA,xlim=c(0,1200),ylim=c(min(val.all),max(val.all)))
# # # lines(dg.seq1.200[dg.seq1.200$m == 200,]$dgmin07[1:1200],col='red',lwd=2)
# # # lines(dg.seq2a.200[dg.seq2a.200$m == 200,]$dgmin07[1:1200],col='purple4',lwd=2)
# # # lines(dg.seq3a.200[dg.seq3a.200$m == 200,]$dgmin07[1:1200],col='orange',lwd=2)
# # # lines(dg.seq2b.200[dg.seq2b.200$m == 200,]$dgmin07[1:1200],col='green4',lwd=2)
# # # lines(dg.seq3b.200[dg.seq3b.200$m == 200,]$dgmin07[1:1200],col='cornflowerblue',lwd=2)
# # # segments(0,0,N1,0,lty=2)
# # # plot(NA,xlim=c(0,1200),ylim=c(min(val.all),max(val.all)))
# # # lines(dg.seq2c.200[dg.seq2c.200$m == 200,]$dgmin07[1:1200],col='blue1',lwd=2)
# # # lines(dg.seq3c.200[dg.seq3c.200$m == 200,]$dgmin07[1:1200],col='blue2',lwd=2)
# # # lines(dg.seq2d.200[dg.seq2d.200$m == 200,]$dgmin07[1:1200],col='blue3',lwd=2)
# # # lines(dg.seq3d.200[dg.seq3d.200$m == 200,]$dgmin07[1:1200],col='blue4',lwd=2)
# # # lines(dg.seq2e.200[dg.seq2e.200$m == 200,]$dgmin07[1:1200],col='purple1',lwd=2)
# # # lines(dg.seq3e.200[dg.seq3e.200$m == 200,]$dgmin07[1:1200],col='purple2',lwd=2)
# # # segments(0,0,N1,0,lty=2)
# # # 
# # 
# # for (i in 1:length(myseqs)) {
# #   seq1 = myseqs[i]
# #   seqtype = seqtypes[i]
# #   par$window = 200
# #   print('fa.get_chunk')
# #   dfseq1 = fa.get_chunk(myseqs[1],window=par$window);size(dfseq1);dfseq1$i0 = dfseq1$i0 + par$window/2
# # }
# # dfseq1 = fa.get_chunk(myseqs[1],window=par$window);size(dfseq1);dfseq1$i0 = dfseq1$i0 + par$window/2
# # dfseq2a = fa.get_chunk(myseqs[2],window=par$window);size(dfseq2a);dfseq2a$i0 = dfseq2a$i0 + par$window/2
# # dfseq3a = fa.get_chunk(myseqs[3],window=par$window);size(dfseq3a);dfseq3a$i0 = dfseq3a$i0 + par$window/2
# # dfseq2b = fa.get_chunk(myseqs[4],window=par$window);size(dfseq2b);dfseq2b$i0 = dfseq2b$i0 + par$window/2
# # dfseq3b = fa.get_chunk(myseqs[5],window=par$window);size(dfseq3b);dfseq3b$i0 = dfseq3b$i0 + par$window/2
# # dfseq2c = fa.get_chunk(myseqs[6],window=par$window);size(dfseq2c);dfseq2c$i0 = dfseq2c$i0 + par$window/2
# # dfseq3c = fa.get_chunk(myseqs[7],window=par$window);size(dfseq3c);dfseq3c$i0 = dfseq3c$i0 + par$window/2
# # dfseq2d = fa.get_chunk(myseqs[8],window=par$window);size(dfseq2d);dfseq2d$i0 = dfseq2d$i0 + par$window/2
# # dfseq3d = fa.get_chunk(myseqs[9],window=par$window);size(dfseq3d);dfseq3d$i0 = dfseq3d$i0 + par$window/2
# # dfseq2e = fa.get_chunk(myseqs[10],window=par$window);size(dfseq2e);dfseq2e$i0 = dfseq2e$i0 + par$window/2
# # dfseq3e = fa.get_chunk(myseqs[11],window=par$window);size(dfseq3e);dfseq3e$i0 = dfseq3e$i0 + par$window/2
# # lay(n1=4,n2=3)
# # plotthis(dfseq1[dfseq1$i0 <= 1500,],dg.seq1.200[dg.seq1.200$m == 1000,],window=200)
# # plotthis(dfseq2a[dfseq2a$i0 <= 1500,],dg.seq2a.200[dg.seq2a.200$m == 1000,],window=200)
# # plotthis(dfseq3a[dfseq3a$i0 <= 1500,],dg.seq3a.200[dg.seq3a.200$m == 1000,],window=200)
# # plotthis(dfseq2b[dfseq2b$i0 <= 1500,],dg.seq2b.200[dg.seq2b.200$m == 1000,],window=200)
# # plotthis(dfseq3b[dfseq3b$i0 <= 1500,],dg.seq3b.200[dg.seq3b.200$m == 1000,],window=200)
# # plotthis(dfseq2c[dfseq2c$i0 <= 1500,],dg.seq2c.200[dg.seq2c.200$m == 1000,],window=200)
# # plotthis(dfseq3c[dfseq3c$i0 <= 1500,],dg.seq3c.200[dg.seq3c.200$m == 1000,],window=200)
# # plotthis(dfseq2d[dfseq2d$i0 <= 1500,],dg.seq2d.200[dg.seq2d.200$m == 1000,],window=200)
# # plotthis(dfseq3d[dfseq3d$i0 <= 1500,],dg.seq3d.200[dg.seq3d.200$m == 1000,],window=200)
# # plotthis(dfseq2e[dfseq2e$i0 <= 1500,],dg.seq2e.200[dg.seq2e.200$m == 1000,],window=200)
# # plotthis(dfseq3e[dfseq3e$i0 <= 1500,],dg.seq3e.200[dg.seq3e.200$m == 1000,],window=200)
# # # 
# # nom = e^(-1*dgm1.n/(R*T))
# # denom = e^(-1*dgm1.all/(R*T))
# dg.seq1.200$probmin00 = exp(-1*dg.seq1.200$dgmin00/(R*T))/sum(exp(-1*dg.seq1.200$dgmin00/(R*T)))
# dg.seq1.200$probmin07 = exp(-1*dg.seq1.200$dgmin07/(R*T))/sum(exp(-1*dg.seq1.200$dgmin07/(R*T)))
# dg.seq1.200$probmin14 = exp(-1*dg.seq1.200$dgmin14/(R*T))/sum(exp(-1*dg.seq1.200$dgmin14/(R*T)))
# 
# dg.seq2a.200$probmin00 = exp(-1*dg.seq2a.200$dgmin00/(R*T))/sum(exp(-1*dg.seq2a.200$dgmin00/(R*T)))
# dg.seq2a.200$probmin07 = exp(-1*dg.seq2a.200$dgmin07/(R*T))/sum(exp(-1*dg.seq2a.200$dgmin07/(R*T)))
# dg.seq2a.200$probmin14 = exp(-1*dg.seq2a.200$dgmin14/(R*T))/sum(exp(-1*dg.seq2a.200$dgmin14/(R*T)))
# dg.seq3a.200$probmin00 = exp(-1*dg.seq3a.200$dgmin00/(R*T))/sum(exp(-1*dg.seq3a.200$dgmin00/(R*T)))
# dg.seq3a.200$probmin07 = exp(-1*dg.seq3a.200$dgmin07/(R*T))/sum(exp(-1*dg.seq3a.200$dgmin07/(R*T)))
# dg.seq3a.200$probmin14 = exp(-1*dg.seq3a.200$dgmin14/(R*T))/sum(exp(-1*dg.seq3a.200$dgmin14/(R*T)))
# 
# dg.seq2b.200$probmin00 = exp(-1*dg.seq2b.200$dgmin00/(R*T))/sum(exp(-1*dg.seq2b.200$dgmin00/(R*T)))
# dg.seq2b.200$probmin07 = exp(-1*dg.seq2b.200$dgmin07/(R*T))/sum(exp(-1*dg.seq2b.200$dgmin07/(R*T)))
# dg.seq2b.200$probmin14 = exp(-1*dg.seq2b.200$dgmin14/(R*T))/sum(exp(-1*dg.seq2b.200$dgmin14/(R*T)))
# dg.seq3b.200$probmin00 = exp(-1*dg.seq3b.200$dgmin00/(R*T))/sum(exp(-1*dg.seq3b.200$dgmin00/(R*T)))
# dg.seq3b.200$probmin07 = exp(-1*dg.seq3b.200$dgmin07/(R*T))/sum(exp(-1*dg.seq3b.200$dgmin07/(R*T)))
# dg.seq3b.200$probmin14 = exp(-1*dg.seq3b.200$dgmin14/(R*T))/sum(exp(-1*dg.seq3b.200$dgmin14/(R*T)))
# 
# dg.seq2c.200$probmin00 = exp(-1*dg.seq2c.200$dgmin00/(R*T))/sum(exp(-1*dg.seq2c.200$dgmin00/(R*T)))
# dg.seq2c.200$probmin07 = exp(-1*dg.seq2c.200$dgmin07/(R*T))/sum(exp(-1*dg.seq2c.200$dgmin07/(R*T)))
# dg.seq2c.200$probmin14 = exp(-1*dg.seq2c.200$dgmin14/(R*T))/sum(exp(-1*dg.seq2c.200$dgmin14/(R*T)))
# dg.seq3c.200$probmin00 = exp(-1*dg.seq3c.200$dgmin00/(R*T))/sum(exp(-1*dg.seq3c.200$dgmin00/(R*T)))
# dg.seq3c.200$probmin07 = exp(-1*dg.seq3c.200$dgmin07/(R*T))/sum(exp(-1*dg.seq3c.200$dgmin07/(R*T)))
# dg.seq3c.200$probmin14 = exp(-1*dg.seq3c.200$dgmin14/(R*T))/sum(exp(-1*dg.seq3c.200$dgmin14/(R*T)))
# 
# dg.seq2d.200$probmin00 = exp(-1*dg.seq2d.200$dgmin00/(R*T))/sum(exp(-1*dg.seq2d.200$dgmin00/(R*T)))
# dg.seq2d.200$probmin07 = exp(-1*dg.seq2d.200$dgmin07/(R*T))/sum(exp(-1*dg.seq2d.200$dgmin07/(R*T)))
# dg.seq2d.200$probmin14 = exp(-1*dg.seq2d.200$dgmin14/(R*T))/sum(exp(-1*dg.seq2d.200$dgmin14/(R*T)))
# dg.seq3d.200$probmin00 = exp(-1*dg.seq3d.200$dgmin00/(R*T))/sum(exp(-1*dg.seq3d.200$dgmin00/(R*T)))
# dg.seq3d.200$probmin07 = exp(-1*dg.seq3d.200$dgmin07/(R*T))/sum(exp(-1*dg.seq3d.200$dgmin07/(R*T)))
# dg.seq3d.200$probmin14 = exp(-1*dg.seq3d.200$dgmin14/(R*T))/sum(exp(-1*dg.seq3d.200$dgmin14/(R*T)))
# 
# dg.seq2e.200$probmin00 = exp(-1*dg.seq2e.200$dgmin00/(R*T))/sum(exp(-1*dg.seq2e.200$dgmin00/(R*T)))
# dg.seq2e.200$probmin07 = exp(-1*dg.seq2e.200$dgmin07/(R*T))/sum(exp(-1*dg.seq2e.200$dgmin07/(R*T)))
# dg.seq2e.200$probmin14 = exp(-1*dg.seq2e.200$dgmin14/(R*T))/sum(exp(-1*dg.seq2e.200$dgmin14/(R*T)))
# dg.seq3e.200$probmin00 = exp(-1*dg.seq3e.200$dgmin00/(R*T))/sum(exp(-1*dg.seq3e.200$dgmin00/(R*T)))
# dg.seq3e.200$probmin07 = exp(-1*dg.seq3e.200$dgmin07/(R*T))/sum(exp(-1*dg.seq3e.200$dgmin07/(R*T)))
# dg.seq3e.200$probmin14 = exp(-1*dg.seq3e.200$dgmin14/(R*T))/sum(exp(-1*dg.seq3e.200$dgmin14/(R*T)))

#dg.all = dg.all1[dg.all1$n == 1,]$RD.min.DD
