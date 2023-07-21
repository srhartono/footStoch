#R library
library(ggplot2)
library(reshape2)
library(MASS)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(fitdistrplus)
library(digest) #MD5
library(swfscMisc)
library(scater)
source('lib/srhlib.R')
source('lib/9_misc_lib.R')


get_peaktype = function(mystring) {
  if (length(grep("PEAK_TEMP",mystring)) > 0) {
    peaktype = "BOT"
  }
  else {

    peaktype = "TOP"
  }
  return(peaktype)
}

debug_df = function(expected,todebug,type,verbose=F,mydata=data.frame,debug=F) {
  if (expected != todebug) {

    cat(type,'debug FAILED\n')
    cat('- expected: ',expected,'\n')
    cat('- debug   : ',todebug,'\n')
    if (verbose == T)  {

      cat(type,': mydata frame with',dim(mydata)[1],'rows',dim(mydata)[2],'columns\n')
    }
  }
  else {

    cat(type,'debug SUCCESSFUL\n')
    cat('- expected: ',expected,'\n')
    cat('- debug   : ',todebug,'\n')
    if (verbose == T)  {

      cat(type,': mydata frame with',dim(mydata)[1],'rows',dim(mydata)[2],'columns\n')
    }
  }
#  return(invisible())
}

fix_treat_name = function(df) {
  dffix = read.table('./resources/misc/CORRECTED_DESC.tsv',sep='\t',header=T)
  colnames(dffix)[4] = 'treat'
  dffix$treat.long = gsub("^.+desc(.+)$","\\1",dffix$newname,perl=T)
  dffix = subset(dffix, select=c('treat','treat.long'))
  
  df = merge(df,dffix,by='treat.long',all=F)
}

parseCLUSTFile = function(file=NA,divby=-1,debug=F,verbose=F) {
  df0 = readRDS(file)
  unique(df0$treat)
  if (divby != -1) {
    df0$divby = divby
  } else if (defined(df0$divby) == FALSE) {
    df0$divby = 50
  }
  df1 = subset(df0,select=c(x0beg,y0beg,x1beg,y1beg,x0end,y0end,x1end,y1end))
  df2 = subset(df0,select=c(x0beg,y0beg,x1beg,y1beg,x0end,y0end,x1end,y1end)) * df0$divby
  colnames(df1) = paste(colnames(df1),'I',sep='')
  df3 = df0[,!colnames(df0) %in% colnames(df2)]
  df3 = cbind(df2, df1, df3)
  df3.col = colnames(df3)
  df3.col[df3.col == 'treat'] = 'treat.long'
  colnames(df3) = df3.col
  
  df3 = fix_treat_name(df3)
  to_return = df3

  return(to_return)
}

parseFASTAFile = function(file=NA,debug=F,verbose=F) {

  if (debug == F) {

    lines = read.table(file)[,1]
    chrs = lines[grepl("^>(.+)$",lines)]
    seqs = lines[!grepl("^>.+$",lines)]
    chrs = gsub("^>(.+)$","\\1",chrs,perl=T)
    df = data.frame(chr=chrs,gene=chrs,seq=seqs)
    
    if (verbose == T) {

      cat(file,dim(df),'\n')
    }
    
    to_return = df
    return(to_return)
  }
  #---- DEBUG ----
  if (debug == T) {

    parseFASTAFile.test = function() {

      df       = parseFASTAFile(FASTAFILES[1])
      debugmd5 = paste(basename(FASTAFILES[1]),digest(df))
      expected = 'invitro_10buf.fa f29c4bd9229050a9de2e75e0752e5e80'
      
      debug_df(expected=expected,todebug=debugmd5,type='FASTAFile',mydata=df,verbose=T)
      return(invisible())
    }
    return(parseFASTAFile.test())
  }
  return(to_return)
}

parsePEAKFile = function(file,debug=F,verbose=F) {

  if (debug == F) {

    df = read.table(file,sep='\t')
    colnames(df) = c('chr','beg','end','peakname','val','strand','treat','VR','thres','threschar')
    
    df$gene = df$chr
    
    df$read = as.character(df$peakname)
    
    if (length(grepl("^[A-Za-z0-9_]+\\.?(m[0-9]+.+)$",df$read,perl=T)) > 0) {

      df$read = as.character(gsub("^[A-Za-z0-9_]+\\.?(m[0-9]+.+)$","\\1",df$read,perl=T))
    }
    
    if (verbose == T) {

      cat(file,dim(df),'\n')
    }
    
    to_return = df
    return(to_return)
  }
  #---- DEBUG ----
  if (debug == T) {

    parsePEAKFile.test = function() {

      df.PCB0  = parsePEAKFile(PEAKFILES[1])
      df.T7    = parsePEAKFile(PEAKFILES[11])
      df       = rbind(df.PCB0,df.T7)
      debugmd5 = paste(basename(PEAKFILES[1]),basename(PEAKFILES[11]),digest(df.PCB0),digest(df.T7),digest(df))
      expected = 'PCB0_PEAK_C.BED T7_INIT_PEAK_C.BED e90f49e2aeaaf7a92c817ab83e346945 42995643cfbc991594a842ff9f00b9ee f7a010ebeb3fc94f7053bf7215d315e3'
      
      debug_df(expected=expected,todebug=debugmd5,type='PEAKFile',mydata=df,verbose=T)
      return(invisible())
    }
    parsePEAKFile.test()
    return(invisible())
  }
  return(to_return)
}

parseBEDFile = function(file,debug=F,verbose=F) {

  if (debug == F) {
    df = read.table(file,sep='\t')
    colnames(df) = c('chr','beg','end','feature','zero','strand')
    df$gene = df$chr
    
    if (verbose == T) {

      cat(file,dim(df),'\n')
    }   
    
    to_return = df
  }  
  #----- DEBUG -----
  if (debug == T) {

    parseBEDFile.test = function() {

      df.annotation = parseBEDFile(BEDFILES[1])
      expected = 'annotation.bed 1f27dafec79f0221d87c0e4b216838a3'
      debugmd5 = paste(basename(BEDFILES[1]),digest(df.annotation))
      
      debug_df(expected=expected,todebug=debugmd5,type='BEDFile',mydata=df,verbose=T)
      return(invisible())
    }
    parseBEDFile.test()
    return(invisible())
  }
  return(to_return)
}

parseMAINFile = function(gp,debug=F,verbose=F) {

  # Parse CLUSTS
  resourcesDir = './resources/misc/'
  CLUSTFile = paste(resourcesDir,dir(resourcesDir,'final3.all.RDS'),sep='')
  CLUSTS = parseCLUSTFile(CLUSTFile,debug=T,divby=gp$divby)

  if (debug == T) {
    expected = '9b076f1c83a66a57e549a067e72c8194'
    todebug    = digest(CLUSTS)
    debug_df(expected=expected,todebug=todebug,type="CLUSTS",mydata=CLUSTS,verbose=T)
  }  
  
  # Parse PEAKS
  PEAKS  = data.frame()
  
  for (file in PEAKFILES) {

    PEAKS   = rbind(PEAKS,data.frame(parsePEAKFile(file),peaktype=get_peaktype(basename(file))))
  }
  
  if (debug == T) {
    expected = 'ac765900dd2ab39320fabd2c246b812e'
    todebug    = digest(PEAKS)
    debug_df(expected=expected,todebug=todebug,type="PEAKS",mydata=PEAKS,verbose=T)
  }  
  # Parse BEDS
  BEDS  = data.frame()
  
  for (file in BEDFILES) {
    BEDS   = rbind(BEDS,parseBEDFile(file))
  }
  
  if (debug == T) {

    expected = 'b8978ebf63df566c3cdcf064871119f6'
    todebug    = digest(BEDS)
    debug_df(expected=expected,todebug=todebug,type="BEDS",mydata=BEDS,verbose=T)
  }
  # Parse FASTAS
  FASTAS  = data.frame()
  
  for (file in FASTAFILES) {

    FASTAS   = rbind(FASTAS,parseFASTAFile(file))
  }
  if (debug == T) {

    expected = '61edd075305195026718a373e0a3db7e'
    todebug    = digest(FASTAS)
    debug_df(expected=expected,todebug=todebug,type="FASTAS",mydata=FASTAS,verbose=T)
  }  
  mylist=list(PEAKS=PEAKS,BEDS=BEDS,FASTAS=FASTAS,CLUST=CLUSTS)
  return(mylist)
}

get_title = function(mytitle='',params=list(),verbose=F,debug=F) {
  
  if (length(params) == 0) {
    
    return(mytitle)
  }
  for (i in 1:length(params))  {
    
    paramname = names(params)[i]
    paramwant = params[i][[1]]
    if (mytitle != '') {
      
      mytitle = paste(mytitle,',',paramname,'_',paramwant,sep='')
    }
    else {
      
      mytitle = paste(paramname,'_',paramwant,sep='')
    }
  }
  if (mytitle == '') {
    
    mytitle = 'NO_TITLE'
  }
  to_return = mytitle
  
  return(to_return)
}

wrap_title = function(mytitle='',width=10,verbose=F,debug=F) {
  
  to_return = paste(strwrap(gsub(',',' ',mytitle),width=width),collapse='\n')
  
  return(to_return)
}

slice_df = function(mydata,params=list(),params_not_exact=c(),verbose=F,debug=F) {

  #gene.want='any',treat.want='any',peaktype.want='any',VR.want='any',thres.want='any',verbose=F,debug=F)
  if (length(params) == 0) {

    return(mydata)
  }
  for (i in 1:length(params))  {

    orig = mydata
    paramname = names(params)[i]
    paramwant = params[i][[1]]
    not_exact = params_not_exact[i]
    if (not_exact == TRUE) {

      cat(i,'. ',paramname,': ',paramwant,' (not_exact): \t',sep='')
    }
    else {

      cat(i,'. ',paramname,': ',paramwant,'\t',sep='')
    }
    
    if (paramwant == 'any') {

      cat(': From',dim(orig)[1],'x',dim(orig)[2],'to',dim(mydata)[1],'x',dim(mydata)[2],'\n',sep=' ')
      next
    }

    mydata = mydata[!is.na(mydata[,paramname]),]
    if (not_exact == FALSE) {

      mydata = mydata[mydata[,paramname] == paramwant,]
    }
    else {

      mydata = mydata[grep(paramwant,mydata[,paramname]),]
    }
    cat(': From',dim(orig)[1],'x',dim(orig)[2],'to',dim(mydata)[1],'x',dim(mydata)[2],'\n',sep=' ')
  }
  return(mydata)
}

slice_CLUSTS = function(df,CLUSTS) {
  dfclust = CLUSTS[CLUSTS$treat %in% df$treat & CLUSTS$gene %in% df$gene,]
  dfclust = subset(dfclust,select=c('cluster','x0end','y0end','x1end','y1end','x0beg','y0beg','x1beg','y1beg'))
  
  clustmin = data.frame(cluster=0,x0end=0,y0end=0,
                        x0beg=0,y0beg=min(dfclust$x0end)-1,
                        x1end=0,y1end=min(dfclust$x0end)-1,
                        x1beg=min(dfclust$x0end)-1,y1beg=3000)
  clustmax = data.frame(cluster=max(dfclust$cluster)+1,
                        x0end=min(dfclust$x0end)-1,y0end=max(dfclust$y1beg)+1,
                        x0beg=max(dfclust$x0end)+1,y0beg=3000,
                        x1end=max(dfclust$x0end)+1,y1end=3000,
                        x1beg=3000,y1beg=3000)
  dfclust=rbind(clustmin,clustmax,dfclust)
  to_return = dfclust
  return(to_return)
}

get_cluster = function(df, dfclust = data.frame(),gp=list()) {
  df$cluster = 0
  dfclustuniq = unique(dfclust$cluster)
  dfclustuniq = dfclustuniq[order(dfclustuniq)]
  df$RIZ = 0
  df$RTZ = 0
  for (curr_cluster in dfclustuniq) {
    curr = dfclust[dfclust$cluster == curr_cluster,]
    if (defined(df[df$beg >= curr$x0end & df$beg <= curr$x1beg & df$end >= curr$y0end & df$end <= curr$y1beg,]$cluster)) {
      df[df$beg >= curr$x0end & df$beg <= curr$x1beg & df$end >= curr$y0end & df$end <= curr$y1beg,]$cluster = curr_cluster
    }
  }
  # for (i in dfclustuniq) {
  #   curr = dfclust[dfclust$cluster == i,]
  #   print(curr)
  #   if (curr$y0end - curr$x0end <= 5000) {
  #     if (dim(df[df$cluster == i & df$beg >= curr$x0end & df$beg <= curr$x0end+gp$divby,])[1] > 0) {
  #       df[df$cluster == i & df$beg >= curr$x0end & df$beg <= curr$x0end+gp$divby,]$RIZ = 1
  #     }
  #   }
  #   if (curr$y0end - curr$x0end <= 5000) {
  #     if (dim(df[df$cluster == i & df$beg >= curr$y0end-gp$divby & df$beg <= curr$y0end,])[1] > 0) {
  #       df[df$cluster == i & df$beg >= curr$y0end-gp$divby & df$beg <= curr$y0end,]$RIZ = 1
  #     }
  #   }
  #   if (curr$y1beg - curr$x1beg <= 5000) {
  #     if (dim(df[df$cluster == i & df$end >= curr$x1beg & df$end <= curr$x1beg+gp$divby,])[1] > 0) {
  #       df[df$cluster == i & df$end >= curr$x1beg & df$end <= curr$x1beg+gp$divby,]$RTZ = 1
  #     }
  #   }
  #   if (curr$y1beg - curr$x1beg <= 5000) {
  #     if (dim(df[df$cluster == i & df$end >= curr$y1beg-gp$divby & df$end <= curr$y1beg,])[1] > 0) {
  #       df[df$cluster == i & df$end >= curr$y1beg-gp$divby & df$end <= curr$y1beg,]$RTZ = 1
  #     }
  #   }
  # }
  to_return = df
  return(to_return)
}

CLUSTFILES  = myorder(paste('resources/misc/',dir("./resources/misc/","*final*.RDS"),sep=''))
PEAKFILES  = myorder(paste('resources/peaks/',dir("./resources/peaks/","*.BED"),sep=''))
BEDFILES   = myorder(paste('resources/bed/',dir("./resources/bed/","*.bed$"),sep=''))
FASTAFILES = myorder(paste('resources/fa/',dir("./resources/fa/","*.fa"),sep=''))

parseBEDFile(debug=T)
parsePEAKFile(debug=T)
parseFASTAFile(debug=T)
#parseCLUSTFile(debug=T)

MAIN = parseMAINFile(gp=gp,debug=T)

CLUSTS  = MAIN$CLUST
BEDS    = MAIN$BED
PEAKS   = MAIN$PEAK
FASTAS  = MAIN$FASTA

mypar = list(
  'genes'=myorder(unique(PEAKS$gene)),
  'treats'=myorder(unique(PEAKS$treat)),
  'VRs'=myorder(unique(PEAKS$VR)),
  'thres'=myorder(unique(PEAKS$thres)),
  'peaktypes'=myorder(unique(PEAKS$peaktype))
)

