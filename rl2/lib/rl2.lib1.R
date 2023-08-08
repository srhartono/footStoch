# parse deltaG file
file.parse_deltaG = function() {
  B = read.table('./resources/tsv/B.tsv',header=T,sep='\t')
  B = melt(B,id.vars = c('nuc','type'))
  B$variable = ac(B$variable)
  colnames(B) = c('n2','type','n1','B')
  N = B[,c(2,3,1,4)]#subset(B,select=c(n1,n2,type,B))
  sapply(B,class)
  B = B[order(B$type,B$n1,B$n2,B$B),]
  h(B)
  B$chunk = pasta(B$n1,B$n2)
  return(B)
}

B.precalc = function(file=NA,file.default=T,force=F,par=par) {
  if (file.default == T) {
    print(paste('rl2.lib1.R::B.precalc: Using default file: resources/B.sigma.RDS'))
    file = 'resources/B.sigma.RDS'
  }
  if (!is.na(file)) {
    if (force == FALSE) {
      if (file.exists(file)) {
        B.sigma = readRDS(file)
        print(paste('rl2.lib1.R::B.precalc: Loaded',file))
        return(B.sigma)
      }
    } else {
      print('rl2.lib1.R::B.precalc: file.RDS exist but force is TRUE, so overwriting!')
    }
    if (!file.exists(file)) {
      print('rl2.lib1.R::B.precalc: file.RDS does not exist, so redoing!')
    }
  }
  B.sigma = data.frame()
  for (i in 1:length(par$sigmas)) {
    sigma = par$sigmas[i]
    sigma.name = as.character(sigma)
    #Calc B.sigma varied by par$len, without sum(B(n+1))
    alpha = par$N*sigma*par$A
    B.sigma.precalc    = (2 * pi^2 * par$C * par$K) / (4 * pi^2 * par$C + par$K * par$len) * (alpha + par$len * par$A)^2
    B.sigma.precalc[1] = B.sigma.precalc[1] #no rloop means no junction energy (a)
    B.sigmatemp = data.frame(index=seq(1,length(par$len)),m=par$len)
    B.sigmatemp[,sigma.name] = B.sigma.precalc
    if (i == 1) {
      B.sigma = B.sigmatemp
    } else {
      B.sigma[,sigma.name] = B.sigma.precalc
    }
  }
  saveRDS(B.sigma,file='resources/B.sigma.RDS')
  print('saved as resources/B.sigma.RDS')
  return(B.sigma)
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

rl2.fa.get_chunk = function(seq1,window=200) {
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


dg.merge = function(dg,seq1) {
  dgseq1 = rl2.fa.get_chunk(seq1,window=1);size(dgseq1)
  dgseq1 = merge(dgseq1,subset(dg[dg$type == 'DD',],select=c(chunk,B)),by='chunk');size(dgseq1);colnames(dgseq1)[dim(dgseq1)[2]] = 'DD'
  dgseq1 = merge(dgseq1,subset(dg[dg$type == 'RD',],select=c(chunk,B)),by='chunk');size(dgseq1);colnames(dgseq1)[dim(dgseq1)[2]] = 'RD'
  dgseq1 = merge(dgseq1,subset(dg[dg$type == 'RD.min.DD',],select=c(chunk,B)),by='chunk');size(dgseq1);colnames(dgseq1)[dim(dgseq1)[2]] = 'RD.min.DD'
  dgseq1 = dgseq1[order(dgseq1$i0,dgseq1$i1,dgseq1$chunk),]
  return(dgseq1)
}

dgb = function(x,dgseq) {
  b3 = sum(dgseq[dgseq$i0 >= x[1] & dgseq$i0 <= x[1]+x[2],]$RD.min.DD)
  return(b3)
}

dg.calc = function(n.all,m.all,dgseq,dg,seqtype){
  # dg.all = data.frame()
  # for (i in seq(1,length(n.all))) {
  #   dg.all.temp = dg
  #   dg.all.temp$n = n.all[i]
  #   dg.all = rbind(dg.all,dg.all.temp)
  # }
  dg.all.temp1 = data.frame(n = rep(seq(1,length(n.all)),dim(dg)[1]*length(n.all)),junk=1)
  dg.all.temp1 = dg.all.temp1[order(dg.all.temp1$n),]
  dg.all.temp2 = dg[rep(seq(1,size(dg)),length(n.all)),]
  dg.all.temp3 = rbind(dg.all.temp1,dg.all.temp2)
  dg.all = subset(dg.all.temp3,select=-junk)
  #  dg.all.temp = data.frame(rep(dg,length(n.all)),n=)
  dg.all$RD.min.DD = apply(subset(dg.all,select=c(n,m)),1,dgb,dgseq=dgseq)
  dg.all$seqtype = seqtype
  return(dg.all)
}

