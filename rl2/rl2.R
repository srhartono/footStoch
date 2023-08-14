options(max.print = 200)
source('rl2/lib/rl2.lib9.R')
source('rl2/lib/rl2.lib1.R')
source('lib/tri.lib.misc.R')


#gp.get_gp() # get gp.orig into .GlobalEnv
#parseMAINFile(gp=gp.orig,debug=T); # get CLUSTS, PEAKS, BEDS, FASTAS into .GlobalEnv

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
BEDS.rl2 = BEDS
BEDS.rl2[grep('invivo',BEDS.rl2$file),]$BEDS
# Change FASTAS.rl2 from -2000
temp1 = function(fa,bed) {
  fa = fa[fa$type == 'invivo',]
  fa$beg = ''
}


par = list(
#  n = seq(1,3000),
  m = seq(0,2000),
  sigmas = seq(0,-0.14,by=-0.01),
  window = 200,
  invivo_buffer = 1990,
  invitro_buffer = 10
)
par$len.max = max(par$m)

seqtypeswant = c("pFC9","T7_init_VR_10","T7_init_VR_20","T7_init_VR_1","T7_init_VR_4","T7_init_VR_17","T7_init_VR_18","T7_init_VR_21","T7_init_VR_25","T7_init_VR_31","pFC53_ApaLI","pFC53TVR_10_pair_0_T3TermMix","pFC53TVR_20_pair_0_T3TermMix",'CALM3','BTF3','FUS','MRFAP1L1','CNBP',
                 "pFC53TVR_1_pair_0_T3TermMix",
                 "pFC53TVR_4_pair_0_T3TermMix"
                 
)
FASTASwant = FASTAS.rl2[FASTAS.rl2$gene %in% seqtypeswant,]
myseqs = c(FASTASwant$seq,seqT,seqG)
seqtypes = c(FASTASwant$gene,"seqT","seqG")

test1 = FASTAS.rl2[FASTAS.rl2$gene == 'T7_init_VR_10',]
test1 = FASTAS.rl2[FASTAS.rl2$gene == 'CALM3',]
nchar(test1$seq)
nchar(test1$seq.amp)
test1$end - test1$beg

length(seqtypes)
length(myseqs)
#myseqs = c(seqVR10,seqVR20,seqT,seqG,seqgood,seqok,seq2a)
#seqtypes = c('T7_init_VR_10','T7_init_VR_20','seqT','seqG','seqgood','seqok','seq2a')
#'seqT','seqG','seq1','seq2a','seq3a','seq2b','seq3b','seq2c','seq3c','seq2d','seq3d','seq2e','seq3e')

#dG = sum(B) + (a + dG.sigma(m,sigma))
B = file.parse_deltaG()
B2 = file.parse_RNAG()
print(paste("- Using max R-loop length (par$m) =",max(par$m)))
#Example
example.G.a.s()
example.G.a.s = function() {
  dgAS = data.frame()
  dgAS = rbind(dgAS,calc.G_a_sigma(m.max=max(par$m),a=10,sigma=-0.00,cons=cons,print=T))
  dgAS = rbind(dgAS,calc.G_a_sigma(m.max=max(par$m),a=10,sigma=-0.01,cons=cons,print=T))
  dgAS = rbind(dgAS,calc.G_a_sigma(m.max=max(par$m),a=10,sigma=-0.04,cons=cons,print=T))
  dgAS = rbind(dgAS,calc.G_a_sigma(m.max=max(par$m),a=10,sigma=-0.07,cons=cons,print=T))
  dgAS = rbind(dgAS,calc.G_a_sigma(m.max=max(par$m),a=10,sigma=-0.11,cons=cons,print=T))
  dgAS = rbind(dgAS,calc.G_a_sigma(m.max=max(par$m),a=10,sigma=-0.14,cons=cons,print=T))
  print(dim(dgAS)) #10000 x 6
}





dg.all = data.frame()
dgseq.all = data.frame()
dfseq.all = data.frame()

par$force = F
for (i in 5:5) {#:length(myseqs)) {
  seq1 = myseqs[i]
  N1 = nchar(seq1)
  seqtype = seqtypes[i]
  print(paste('\n',i,'/',length(myseqs),': Doing',seqtypes[i]))
  file.dgseq = pasta('rl2/resources/dgtemp/',seqtype,'_',par$window,'_dgseq.RDS')
  file.dgB = pasta('rl2/resources/dgtemp/',seqtype,'_',par$window,'_dgB.RDS')

  print(paste('- rl2.fa.get_chunk',seqtype))
  # if (file.exists(file.dgseq)) { #} & par$force == F) {
  #   dgseq1 = readRDS(file.dgseq)
  # } else {
    dfseq1 = rl2.fa.get_chunk(seq1,window=par$window,xlim=c(0,N1),print=F,my.title=seqtype);size(dfseq1);
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
    saveRDS(dgseq1,file=file.dgseq)
    print('  -> save')
  # }
  p = ggplot(dgseq1,aes(n,B)) +
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

  print(paste('- dg.calc',seqtype))
  print(paste('nchar seq1 =',N1,', len.max =',par$len.max,'but using',min(N1,par$len.max)))
  calc.time(N1,min(N1,par$lenmax))

#  if (N1 > 4000) {next}

  # if (file.exists(file.dgB)) {#} & par$force == F) {
  #   print(paste('- ',file.dgB,'exist!'))
  #   dgB1 = readRDS(file=file.dgB)
  # } else {

    par.n = seq(1,N1)
    par.m = par$m[par$m <= N1]
    dgB1 = dgB.calc(par.n,par.m,dgseq1$B)
    dgB1[dgB1$m == 0,]$B = 0

    saveRDS(dgB1,file=file.dgB)
    print(paste('- ',file.dgB,'saved!'))
  # }
}


