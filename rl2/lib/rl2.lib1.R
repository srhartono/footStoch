# parse deltaG file
file.parse_deltaG = function() {
  B = read.table('rl2/resources/tsv/B.tsv',header=T,sep='\t')
  B = melt(B,id.vars = c('nuc','type'))
  B$variable = ac(B$variable)
  colnames(B) = c('n2','type','n1','B')
  N = B[,c(2,3,1,4)]#subset(B,select=c(n1,n2,type,B))
  sapply(B,class)
  B$chunk = pasta(B$n1,B$n2)
  B0 = dg.get_monobase_deltaG(B)
  B = rbind(B,B0)
  B = B[order(B$type,B$n1,B$n2,B$B),]
  h(B)
  return(B)
}

dg.get_monobase_deltaG = function(B) {
  B0 = B[B$n1 == B$n2,]
  B0$chunk = B0$n1
  B0$B = B0$B / 2
  BN1 = B0
  BN2 = B0
  BN1$chunk = pasta(BN1$chunk, 'N')
  BN1$n2 = 'N'
  BN2$chunk = pasta('N',BN2$chunk)
  BN2$n1 = 'N'
  print(head(BN1))
  print(head(BN2))
  B0 = rbind(B0,BN1)
  B0 = rbind(B0,BN2)
  B0 = rbind(B0,data.frame(n1='N',n2='N',type=unique(B$type),B=0,chunk='NN'))
  return(B0)
}
B = file.parse_deltaG()
file.parse_RNAG = function() {
  B = read.table('rl2/resources/tsv/RNAG.tsv',header=T,sep='\t')
  B = melt(B,id.vars = c('nuc','type'))
  B$variable = ac(B$variable)
  colnames(B) = c('n2','type','n1','B')
  N = B[,c(2,3,1,4)]#subset(B,select=c(n1,n2,type,B))
  sapply(B,class)
  B = B[order(B$type,B$n1,B$n2,B$B),]
  h(B)
  B$chunk = pasta(B$n1,B$n2)
  chunk2 = sapply(B$chunk,revcomp)
  B$n1 = B$chunk
  B$n2 = chunk2
  B$chunk = pasta(B$n1,B$n2)
  return(B)
}

rl2.get_constants  = function(mylist=list()) {
  cons = list(
    a = 10,
    C = 1.8,
    A = 1/10.4,
    N = 1500,
    R = 0.0019858775,
    T = 310,
    K = NA
  ) 
  if (length(mylist) > 0) {
    for (i in 1:length(mylist)) {
      name = names(mylist[i])
      if (name == 'K') {
        message("DONT DIRECTLY MODIFY cons$K!!! Change cons$R, cons$T, and cons$N instead!!!!")
      } else {
        cons[name] = mylist[name]
      }
    }
  }
  cons$K = (2200 * cons$R * cons$T) / cons$N
  return(cons)
}


rl2.fa.get_FASTAS = function(df.fasta=NA,force=F,file=NA,use.default=T) {
  file.default = 'rl2/resources/FASTAS.rl2.RDS'
  
  if (is.na(file)) {
    #print("Using file.default")
    file = file.default
  }
  if (defined(df.fasta) == FALSE) {
    print("Using FASTAS")
    df.fasta = FASTAS
  }
  if (force == F) {
    if (file.exists(file)) {
      FASTAS.rl2 = readRDS(file)
      cat('rl2.fa.get_FASTAS.rl2: LOADED ',file,'\n')
      
      return(FASTAS.rl2)
    }
  }
  #print("Using df.fasta")
  FASTAS.rl2 = df.fasta
  FASTAS.rl2[FASTAS.rl2$strand == '-',]$seq = sapply(FASTAS.rl2[FASTAS.rl2$strand == '-',]$seq, revcomp)
  FASTAS.rl2[FASTAS.rl2$strand == '-',]$seq.amp = sapply(FASTAS.rl2[FASTAS.rl2$strand == '-',]$seq.amp, revcomp)
  FASTAS.rl2$type = 'invivo'
  FASTAS.rl2[grep('invivo',FASTAS.rl2$file,invert=T),]$type = 'invitro'
  FASTAS.rl2 = FASTAS.rl2[order(FASTAS.rl2$type,FASTAS.rl2$gene),]
  saveRDS(FASTAS.rl2,file)
  cat('rl2.fa.get_FASTAS.rl2: SAVED ',file,'; dim=',paste(dim(FASTAS.rl2),collapse=' x '),'\n')
  return(FASTAS.rl2)
}

calc.G_a_sigma = function(m.min=0,m.max,a,sigma,cons,file=NA,use.default=T,force=T,print=F,my.title='') {
  sigma.name = pasta('min',gsub("^\\-","",as.character(sigma)))
  m.precalc = seq(m.min,m.max)
  file.default = pasta('resources/variables/dg.sigma.',sigma.name,'.',max(m.precalc),'.RDS')
  
  #calculate all B sigma for length nmax
  if (use.default == T) {
    file = file.default
    #print(paste('rl2.lib1.R::B.precalc: Using default file:',file))
  }
  if (!is.na(file)) {
    if (force == FALSE) {
      if (file.exists(file)) {
        dg.sigma = readRDS(file)
        #print(paste('rl2.lib1.R::B.precalc: Loaded',file))
        return(dg.sigma)
      }
    # } else {
    #   print('rl2.lib1.R::B.precalc: file.RDS exist but force is TRUE, so overwriting!')
    }
    # if (!file.exists(file)) {
    #   print('rl2.lib1.R::B.precalc: file.RDS does not exist, so redoing!')
    # }
  }
  dg.sigma = data.frame()
  
  #Calc dg.sigma varied by m.precalc, without sum(B(n+1))
  alpha = cons$N * sigma * cons$A
  dg.sigma = data.frame(
    m=m.precalc,
    a = a,
    sigma = sigma,
    dg.sigma = ((2 * pi^2 * cons$C * cons$K) * (alpha + m.precalc * cons$A)^2)/ (4 * pi^2 * cons$C + cons$K * m.precalc) 
  )
  dg.sigma[dg.sigma$m > 0,]$dg.sigma = a + dg.sigma[dg.sigma$m > 0,]$dg.sigma
  saveRDS(dg.sigma,file=file)
  #print(paste('saved as',file,'; dim =',paste(dim(dg.sigma),collapse=' ')))
  cat("\n")
  if (print == T) {
    plot(dg.sigma$m,dg.sigma$dg.sigma,type='l',xlim=c(0,1000),ylim=c(0,500),main=paste(my.title,'\nΔG of R-loop formation minus sum(B), varied by R-loop length\na =',a,'; σ =',sigma),xlab='m\nR-loop length (bp)',ylab='ΔG (kcal/mol)')
    segments(0,46.01797,1000,46.01797,lty=2)
  }
  return(dg.sigma)
}



fa.gen_random = function(n=2000) {
  # Generate random seq1
  mylist = list()
  set.seed(420)
  mylist$nuc.bin1 = c('A','C','T','G')
  mylist$nuc.bin2 = c(rep('G',9),rep('C',3),'A','T') #0.5 GC skew, 0.85 GC cont
  mylist$nuc.bin3 = c(rep('G',6),rep('C',6),'A','T') #0 GC skew, 0.85 GC cont
  mylist$seq.1a = join(sample(mylist$nuc.bin1,size = n/10,replace=T))
  mylist$seq.1b = join(sample(mylist$nuc.bin1,size = n/10,replace=T))
  mylist$seq.1c = join(sample(mylist$nuc.bin1,size = n/10,replace=T))
  mylist$seq.1d = join(sample(mylist$nuc.bin1,size = n/10,replace=T))
  mylist$seq.1e = join(sample(mylist$nuc.bin1,size = n/10,replace=T))
  mylist$seq.1f = join(sample(mylist$nuc.bin1,size=n/10,replace=T))
  mylist$seq.1g = join(sample(mylist$nuc.bin1,size=n/10*4,replace=T))
  
  mylist$seq.2a = join(sample(mylist$nuc.bin2,size=n/10,replace=T))
  mylist$seq.2b = join(sample(mylist$nuc.bin2,size=n/10,replace=T))
  mylist$seq.3a = join(sample(mylist$nuc.bin3,size=n/10,replace=T))
  mylist$seq.3b = join(sample(mylist$nuc.bin3,size=n/10,replace=T))
  mylist$seq = list()
  mylist$seq$seqT = join(rep('T',n))
  mylist$seq$seqG = join(rep('G',n))
  mylist$seq$seqgood = join(sample(mylist$nuc.bin2,size=n,replace=T))
  mylist$seq$seqok = join(sample(mylist$nuc.bin3,size=n,replace=T))
  mylist$seq$seq1a = join(c(mylist$seq.1a,mylist$seq.1b,mylist$seq.1c,mylist$seq.1d,mylist$seq.1e,mylist$seq.1f,mylist$seq.1g))
  mylist$seq$seq1b = join(c(mylist$seq.1a,mylist$seq.1b,mylist$seq.1c,mylist$seq.1d,mylist$seq.1e,mylist$seq.1f,mylist$seq.1g))
  mylist$seq$seq2a = join(c(mylist$seq.2a,mylist$seq.2b,mylist$seq.1c,mylist$seq.1d,mylist$seq.1e,mylist$seq.1f,mylist$seq.1g))
  mylist$seq$seq3a = join(c(mylist$seq.3a,mylist$seq.3b,mylist$seq.1c,mylist$seq.1d,mylist$seq.1e,mylist$seq.1f,mylist$seq.1g))
  mylist$seq$seq2b = join(c(mylist$seq.1a,mylist$seq.2a,mylist$seq.2b,mylist$seq.1d,mylist$seq.1e,mylist$seq.1f,mylist$seq.1g))
  mylist$seq$seq3b = join(c(mylist$seq.1a,mylist$seq.3a,mylist$seq.3b,mylist$seq.1d,mylist$seq.1e,mylist$seq.1f,mylist$seq.1g))
  mylist$seq$seq2c = join(c(mylist$seq.1a,mylist$seq.1b,mylist$seq.2a,mylist$seq.2b,mylist$seq.1e,mylist$seq.1f,mylist$seq.1g))
  mylist$seq$seq3c = join(c(mylist$seq.1a,mylist$seq.1b,mylist$seq.3a,mylist$seq.3b,mylist$seq.1e,mylist$seq.1f,mylist$seq.1g))
  mylist$seq$seq2d = join(c(mylist$seq.1a,mylist$seq.1b,mylist$seq.1c,mylist$seq.2a,mylist$seq.2b,mylist$seq.1f,mylist$seq.1g))
  mylist$seq$seq3d = join(c(mylist$seq.1a,mylist$seq.1b,mylist$seq.1c,mylist$seq.3a,mylist$seq.3b,mylist$seq.1f,mylist$seq.1g))
  mylist$seq$seq2e = join(c(mylist$seq.1a,mylist$seq.1b,mylist$seq.1c,mylist$seq.1d,mylist$seq.2a,mylist$seq.2b,mylist$seq.1g))
  mylist$seq$seq3e = join(c(mylist$seq.1a,mylist$seq.1b,mylist$seq.1c,mylist$seq.1d,mylist$seq.3a,mylist$seq.3b,mylist$seq.1g))
  # return(mylist)  
  
  for (seqnames in names(mylist$seq)) {
    temp.len = nchar(mylist$seq[seqnames])
    if (temp.len < n) {
      print(paste(seqnames,temp.len,n))
      temp.diff = n - temp.len
      mylist$seq$seqnames = pasta(mylist$seq[seqnames],join(sample(mylist$nuc.bin1,temp.diff,replace=T)))
    } else if (temp.len > n) {
      print(paste(seqnames,temp.len,n))
      temp.diff = n - temp.len
      mylist$seq[seqnames] = substr(mylist$seq[seqnames],1,n)
    }
  }
  return(mylist)  
}

rl2.fa.get_chunk = function(seq1,window=200,xlim=NA,ylim=c(0,1.4),print=F,my.title="") {
  if (window > nchar(seq1)) {
    window = nchar(seq1)-1
  }
  my.i0 = seq(1,nchar(seq1))
  my.i1 = sapply(my.i0,function(x,seq1){min(x+window,nchar(seq1))},seq1=seq1)
  df = data.frame(i0=my.i0,i1=my.i1)
  df$chunk = apply(df,1,function(x,seq1){substr(seq1,x[1],x[2])},seq1=seq1)
  df2 = t(sapply(df$chunk,fa.calc_CpGprof))
  rownames(df2) = rownames(df)
  colnames(df2) = c('CGdens','GCcont','GCskew')
  df = cbind(df,df2)
  df$n = df$i0
  if (print == TRUE) {
    if (!defined(xlim)) {xlim = c(min(df$i0),max(df$i0))}
    plot(df$n,df$CGdens,type='l',col='cornflowerblue',main=paste(my.title,'CpG Profile'),
         xlim=xlim, ylim=ylim,
         xlab="Position in plasmid (bp)",
         ylab="CpG density/GC content/(GC skew +0.5)")
    lines(df$n,df$GCcont,type='l',col='green4')
    lines(df$n,df$GCskew+0.5,type='l',col='red2')
    segments(0,0.5,max(df$n),0.5,lty=2)
  }
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

dg.merge = function(dg,seq1) {
  dgseq1 = rl2.fa.get_chunk(seq1,window=1);size(dgseq1)
  dgseq1 = merge(dgseq1,subset(dg[dg$type == 'DD',],select=c(chunk,B)),by='chunk');size(dgseq1);colnames(dgseq1)[dim(dgseq1)[2]] = 'DD'
  dgseq1 = merge(dgseq1,subset(dg[dg$type == 'RD',],select=c(chunk,B)),by='chunk');size(dgseq1);colnames(dgseq1)[dim(dgseq1)[2]] = 'RD'
  dgseq1 = merge(dgseq1,subset(dg[dg$type == 'RD.min.DD',],select=c(chunk,B)),by='chunk');size(dgseq1);colnames(dgseq1)[dim(dgseq1)[2]] = 'RD.min.DD'
  dgseq1 = dgseq1[order(dgseq1$i0,dgseq1$i1,dgseq1$chunk),]
  return(dgseq1)
}

dgB.calc = function(par.n,par.m,dgseqB) {
  dgb = function(x,dg) {
    b3 = sum(dg[x[1]:min(x[1]+(x[2]-1),length(dg))])
    return(b3)
  }
  dgB1 = as.data.frame(expand.grid(n=par.n,m=par.m))
  print(size(dgB1))
  dgB1 = dgB1[order(dgB1$n,dgB1$m),]
  start.time <- Sys.time()
  dgB1$B = apply(dgB1,1,dgb,dg=dgseqB)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)  
  return(dgB1)
}  
# 
# dg.calc = function(n.all,m.all,dgseq,dg,seqtype){
#   # dg.all = data.frame()
#   # for (i in seq(1,length(n.all))) {
#   #   dg.all.temp = dg
#   #   dg.all.temp$n = n.all[i]
#   #   dg.all = rbind(dg.all,dg.all.temp)
#   # }
#   dg.all.temp1 = data.frame(n = rep(seq(1,length(n.all)),dim(dg)[1]*length(n.all)),junk=1)
#   dg.all.temp1 = dg.all.temp1[order(dg.all.temp1$n),]
#   dg.all.temp2 = dg[rep(seq(1,size(dg)),length(n.all)),]
#   dg.all.temp3 = rbind(dg.all.temp1,dg.all.temp2)
#   dg.all = subset(dg.all.temp3,select=-junk)
#   #  dg.all.temp = data.frame(rep(dg,length(n.all)),n=)
#   dg.all$RD.min.DD = apply(subset(dg.all,select=c(n,m)),1,dgb,dgseq=dgseq)
#   dg.all$seqtype = seqtype
#   return(dg.all)
# }
# 
