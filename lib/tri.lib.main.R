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
library(pdfCluster)
library(plotly)
library(grid)
library(ggplot2)
library(gridExtra)

source('lib/tri.lib.misc.R')
source('lib/tri.lib.graph_ps.R')

gp.get_gp = function() {
  gp.orig = list(
    minX = 0,
    minY = 0,
    maxX = 3000,
    maxY = 3000,
    windowsmooth=25,
    stepsmooth=1,
    windowdist=50,
    stepdist=1,
    minpval = 0.25,
    minpval2 = 0.25,
    dist.test_min.length = 5,
    dist.test_min.unique = 2,
    dist.test_max.unique.hots = 2,
    divby=25,
    C.transform=FALSE,
    use_clusterFile = FALSE
  )
  .GlobalEnv$gp.orig = gp.orig
}

gp.check_gp = function(gp1,gp2,verbose=F,debug=F) {
  if (digest(gp1) != digest(gp2)) {
    print(pasta('gp and old gp isnt same!\n-',digest(gp1),'\n- ',digest(gp2)))
  } else {
    if (verbose == T) {
      print(pasta('gp and gp.lib are same!\n- ',digest(gp1),'\n- ',digest(gp2)))
    }
  }
  return(gp1)
}

get_peaktype = function(my.string) {
  if (length(grep("PEAK_TEMP",my.string)) > 0) {
    peaktype = "BOT"
  }
  else {
    
    peaktype = "TOP"
  }
  return(peaktype)
}

debug_df = function(expected,todebug,type,verbose=F,my.data=data.frame,debug=F) {
  if (expected != todebug) {
    
    cat(type,'debug FAILED\n')
    cat('- expected: ',expected,'\n')
    cat('- debug   : ',todebug,'\n')
    if (verbose == T)  {
      
      cat(type,': my.data frame with',dim(my.data)[1],'rows',dim(my.data)[2],'columns\n')
    }
  }
  else {
    
    cat(type,'debug SUCCESSFUL\n')
    cat('- expected: ',expected,'\n')
    cat('- debug   : ',todebug,'\n')
    if (verbose == T)  {
      
      cat(type,': my.data frame with',dim(my.data)[1],'rows',dim(my.data)[2],'columns\n')
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
      
      debug_df(expected=expected,todebug=debugmd5,type='FASTAFile',my.data=df,verbose=T)
      return(invisible())
    }
    return(parseFASTAFile.test())
  }
  return(to_return)
}

pos.convert = function(pos.buf1,buf1,buf2) {
  #based on -10 buf, peak at peak.true = 0 -> peak.10 = 10; peak.true 78 = peak.10 at 88.
  # .:. peak.true = peak.buf - buf
  ## pos.true = pos.buf1 - buf1
  ## peak.buf = peak.true + buf
  
  # then to make it based on -2000 buf, peak.2000 = peak.true + 2000
  ## pos.buf2 = pos.true + buf2
  ## return(pos.buf2)
  return(pos.buf1 - buf1 + buf2)
}
    
parsePEAKFile = function(file,debug=F,verbose=F) {
  if (grepl('PCB0',file)) {
    is.invivo = TRUE
  } else {
    is.invivo = FALSE
  }
  if (debug == F) {
    
    df = read.table(file,sep='\t')
    colnames(df) = c('chr','beg','end','peakname','val','strand','treat','VR','thres','threschar')
    
    df$gene = df$chr
    df$read = as.character(df$peakname)
    if (is.invivo == TRUE) { #buf1 = 10, buf2 = 2000
      df$beg = sapply(df$beg,pos.convert,buf1=10,buf2 = 2000)
      df$end = sapply(df$end,pos.convert,buf1=10,buf2 = 2000)
    }
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
      
      debug_df(expected=expected,todebug=debugmd5,type='PEAKFile',my.data=df,verbose=T)
      return(invisible())
    }
    parsePEAKFile.test()
    return(invisible())
  }
  return(to_return)
}

parseBEDFile = function(file,debug=F,verbose=F) {
  if (grepl('invivo',file)) {
    is.invivo = TRUE
  } else {
    is.invivo = FALSE
  }
  
  if (debug == F) {
    df = read.table(file,sep='\t')
    colnames(df) = c('chr','beg','end','feature','zero','strand')
    df$gene = df$chr
    if (defined(df[df$gene == 'RPPH1',])) {
      df = df[df$gene !='RPPH1',]
    }
    if (is.invivo == TRUE) {
      df$beg = sapply(df$beg,pos.convert,10,2000)
      df$end = sapply(df$end,pos.convert,10,2000)
    } 
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
      
      debug_df(expected=expected,todebug=debugmd5,type='BEDFile',my.data=df,verbose=T)
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
    debug_df(expected=expected,todebug=todebug,type="CLUSTS",my.data=CLUSTS,verbose=T)
  }  
  
  # Parse PEAKS
  PEAKS  = data.frame()
  
  for (file in PEAKFILES) {
    
    PEAKS   = rbind(PEAKS,data.frame(parsePEAKFile(file),peaktype=get_peaktype(basename(file))))
  }
  
  if (debug == T) {
    expected = 'ac765900dd2ab39320fabd2c246b812e'
    todebug    = digest(PEAKS)
    debug_df(expected=expected,todebug=todebug,type="PEAKS",my.data=PEAKS,verbose=T)
  }  
  # Parse BEDS
  BEDS  = data.frame()
  
  for (file in BEDFILES) {
    CURR.BEDS = parseBEDFile(file)
    CURR.BEDS$file = file
    # if (grepl('invivo',file)) {
    #   begorig = CURR.BEDS$beg
    #   # begcurr = CURR.BEDS$beg - 10
    #   # endcurr = CURR.BEDS$end + 10
    #   # 
    #   # CURR.BEDS$beg = begcurr - begorig
    #   # CURR.BEDS$end = endcurr - begorig
    #   
    #   CURR.BEDS$beg = (CURR.BEDS$beg - begorig + 2000) #-10 + 2000 = 1990
    #   CURR.BEDS$end = (CURR.BEDS$end - begorig + 2000) #2010 + 2000 = 4010
    #   CURR.BEDS$chr = CURR.BEDS$feature
    #   CURR.BEDS$gene = CURR.BEDS$feature
    #   CURR.BEDS.FW = CURR.BEDS
    #   CURR.BEDS.RV = CURR.BEDS
    #   CURR.BEDS.FW$end = CURR.BEDS.FW$beg + 30
    #   CURR.BEDS.RV$beg = CURR.BEDS.RV$end - 30
    #   CURR.BEDS.FW$feature = "FW_Primer"
    #   CURR.BEDS.RV$feature = "RV_Primer"
    #   CURR.BEDS = rbind(CURR.BEDS.FW,CURR.BEDS.RV)
    #   CURR.BEDS = CURR.BEDS[order(CURR.BEDS$chr,CURR.BEDS$beg,CURR.BEDS$end),]
    # }
    # if (grepl('invivo',file)) {
    #   next
    # } else {
    #   CURR.BEDS.FW = subset(CURR.BEDS[CURR.BEDS$feature == "FW_Primer",],select=-end)
    #   CURR.BEDS.RV = subset(CURR.BEDS[CURR.BEDS$feature == "RV_Primer",],select=c('chr','end'))
    #   CURR.BEDS.FW$feature = "Amplicon"
    #   CURR.BEDS.FW = merge(CURR.BEDS.FW,CURR.BEDS.RV,by='chr')
    #   CURR.BEDS$beg = CURR.BEDS$beg - 10
    #   CURR.BEDS$end = CURR.BEDS$end + 10
    # }
    BEDS   = rbind(BEDS,CURR.BEDS)
  }
  
  BEDS2 = BEDS[grep('invivo',BEDS$file),]
  BEDS = BEDS[grep('invivo',BEDS$file,invert=T),]
  BEDS.FW = subset(BEDS[BEDS$feature == "FW_Primer",],select=-end)
  BEDS.RV = subset(BEDS[BEDS$feature == "RV_Primer",],select=c('chr','end'))
  BEDS.FW$feature = "Amplicon"
  BEDS.FW = merge(BEDS.FW,BEDS.RV,by='chr')
  BEDS.FW$beg = BEDS.FW$beg - 10
  BEDS.FW$end = BEDS.FW$end + 10 
  BEDS = rbind(BEDS,BEDS.FW)
  BEDS = rbind(BEDS,BEDS2)

  if (debug == T) {
    
    expected = 'e75c59cfffd53749823e9a247dc648db  '
    todebug    = digest(BEDS)
    debug_df(expected=expected,todebug=todebug,type="BEDS",my.data=BEDS,verbose=T)
  }
  # Parse FASTAS
  FASTAS  = data.frame()
  
  for (file in FASTAFILES) {
    
    FASTAS   = rbind(FASTAS,parseFASTAFile(file))
  }
  if (debug == T) {
    
    expected = 'ed32b8d93bc599aa11de3003dd1e01d3'
    todebug    = digest(FASTAS)
    debug_df(expected=expected,todebug=todebug,type="FASTAS",my.data=FASTAS,verbose=T)
  }
  
  # Parse BIGFASTAS
  BIGFASTAS  = data.frame()
  
  for (file in BIGFASTAFILES) {
    
    BIGFASTAS   = rbind(BIGFASTAS,parseFASTAFile(file))
  }
  if (debug == T) {
    
    expected = 'dab2095bdbbd4d63ce7428724435644f'
    todebug    = digest(BIGFASTAS)
    debug_df(expected=expected,todebug=todebug,type="BIGFASTAS",my.data=BIGFASTAS,verbose=T)
  }
  
  BIGFASTAS = BIGFASTAS[,-1]
  colnames(FASTAS)[3] = 'seq.amp'
  FASTAS = merge(BIGFASTAS,FASTAS,by='gene')
  FASTAS = subset(FASTAS,select=c('chr','gene','seq','seq.amp'))
  if (debug == T) {
    
    expected = '470c8ca5e36cc56114827887c3adedcc'
    todebug    = digest(FASTAS)
    debug_df(expected=expected,todebug=todebug,type="FASTAS",my.data=FASTAS,verbose=T)
  }
  
  FASTAS = merge(FASTAS,subset(BEDS[BEDS$feature == 'Amplicon',],select=-gene),by='chr')
  # FASTAS[FASTAS$feature == "Amplicon",]$beg = FASTAS[FASTAS$feature == "Amplicon",]$beg - 10
  # FASTAS[FASTAS$feature == "Amplicon",]$end = FASTAS[FASTAS$feature == "Amplicon",]$end + 10
  # 
  # sanity chekc
  FASTAS[!FASTAS$gene %in% BIGFASTAS$gene,]$gene
  BIGFASTAS[!BIGFASTAS$gene %in% FASTAS$gene,]$gene
  (BEDS[!BEDS$gene %in% FASTAS$gene,]$gene)
  (FASTAS[!FASTAS$gene %in% BEDS$gene,]$gene)
  
  #sanity check
  test1 = FASTAS
  test1$len = test1$end - test1$beg
  test1$seq.len = nchar(test1$seq)
  test1$seq.amp.len = nchar(test1$seq.amp)
  test1 = subset(test1,select=c(-seq,-seq.amp))
  head(test1)
  size(test1[test1$len != test1$seq.amp.len,])
  
  my.list=list(PEAKS=PEAKS,BEDS=BEDS,FASTAS=FASTAS,CLUST=CLUSTS)
  .GlobalEnv$BEDS = BEDS
  .GlobalEnv$FASTAS = FASTAS
  .GlobalEnv$CLUSTS = CLUSTS
  .GlobalEnv$PEAKS = PEAKS
}


get_title = function(my.title='',my.params=list(),verbose=F,debug=F) {
  
  if (length(my.params) == 0) {
    
    return(my.title)
  }
  for (i in 1:length(my.params))  {
    
    paramname = names(my.params)[i]
    paramwant = my.params[i][[1]]
    if (my.title != '') {
      
      my.title = paste(my.title,',',paramname,'_',paramwant,sep='')
    }
    else {
      
      my.title = paste(paramname,'_',paramwant,sep='')
    }
  }
  if (my.title == '') {
    
    my.title = 'NO_TITLE'
  }
  to_return = my.title
  
  return(to_return)
}

wrap_title = function(my.title='',width=10,verbose=F,debug=F) {
  
  to_return = paste(strwrap(gsub(',',' ',my.title),width=width),collapse='\n')
  
  return(to_return)
}

slice_bed = function(my.beds=BEDS,my.params=list(),my.params_regex=c(),genewant = NA,verbose=F,debug=F) {
  if (is.na(genewant)) {
    if (defined(my.params_regex$gene)) {
      if (my.params_regex$gene == TRUE) {
        my.bed = my.beds[grep(my.params$gene,my.beds$gene,ignore.case = TRUE),]
      } else {
        my.bed = my.beds[my.beds$gene == my.params$gene,]
      }
    } else {
      my.bed = my.beds[my.beds$gene == my.params$gene,]
    }
    if (defined(my.bed) == FALSE) {
      my.bed = my.beds[grep(paste('^',my.params$gene,'$',sep=''),my.beds$gene,ignore.case = TRUE),]
    }
  } else {
    my.bed = my.beds[grep(paste('^',genewant,'$',sep=''),my.beds$gene),]
    if (defined(my.bed) == FALSE) {
      my.bed = my.beds[grep(paste('^',genewant,'$',sep=''),my.beds$gene,ignore.case = TRUE),]
    }
  }
  return(my.bed)
}

slice_df = function(my.data,my.params=list(),my.params_regex=list(),verbose=F,debug=F) {
  
  #gene.want='any',treat.want='any',peaktype.want='any',VR.want='any',thres.want='any',verbose=F,debug=F)
  if (length(my.params) == 0) {
    
    return(my.data)
  }
  for (i in 1:length(my.params))  {
    
    orig = my.data
    paramname = names(my.params)[i]
    paramwant = my.params[i][[1]]
    isregex = my.params_regex[i]
    if (isregex == TRUE) {
      
      cat(i,'. ',paramname,': ',paramwant,' (regex): \t',sep='')
    }
    else {
      
      cat(i,'. ',paramname,': ',paramwant,'\t',sep='')
    }
    
    if (paramwant == 'any') {
      
      cat(': From',dim(orig)[1],'x',dim(orig)[2],'to',dim(my.data)[1],'x',dim(my.data)[2],'\n',sep=' ')
      next
    }
    
    my.data = my.data[!is.na(my.data[,paramname]),]
    if (isregex == FALSE) {
      
      my.data = my.data[my.data[,paramname] == paramwant,]
    }
    else {
      
      my.data = my.data[grep(paramwant,my.data[,paramname]),]
    }
    cat(': From',dim(orig)[1],'x',dim(orig)[2],'to',dim(my.data)[1],'x',dim(my.data)[2],'\n',sep=' ')
  }
  return(my.data)
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
  df$cluster = -1
  #dfclustuniq = unique(dfclust$cluster)
  #dfclustuniq = dfclustuniq[order(dfclustuniq)]
  dfclustuniq = dfclust$cluster
  df$RIZ = 0
  df$RTZ = 0
  for (clusterInd in seq(1,length(dfclustuniq))) {
    curr_cluster = dfclustuniq[clusterInd]
    curr = dfclust[dfclust$cluster == curr_cluster,]
    if (defined(df[df$cluster == -1 & df$beg >= curr$x0end & df$beg <= curr$x1beg & df$end >= curr$y0end & df$end <= curr$y1beg,]$cluster)) {
      df[df$cluster == -1 & df$beg >= curr$x0end & df$beg <= curr$x1beg & df$end >= curr$y0end & df$end <= curr$y1beg,]$cluster = curr_cluster
    }
  }
  to_return = df
  return(to_return)
}

do_distributionTests = function(df = data.frame(), positionTypes = c('mean','orig'), testTypes=c('unif','norm','hots'), gp=list(), outpdf=F,my.print=F, debug=F, verbose=F) {
  
  #  if (defined(gp$minpval) == F) {
  #    gp$minpval = 0.05
  #  }
  
  do_distributionTests.get_outpdf = function(positionTypes=positionTypes,testTypes=testTypes,gp=gp,outpdf=outpdf) {
    if (outpdf == F) {
      outpdf = paste(paste(positionTypes,sep='_',collapse='_'),'_',paste(testTypes,sep='_',collapse='_'),'.pdf',collapse='',sep='')
      print(outpdf)
      if (defined(gp$my.title)) {
        outpdf = paste(gp$my.title,'_',outpdf,sep='')
      }
      to_return = outpdf
      return(to_return)
    }
  }
  
  outpdf = do_distributionTests.get_outpdf(outpdf=outpdf,gp=gp,positionTypes = positionTypes,testTypes=testTypes)
  
  if (verbose == T) {print(outpdf)}
  
  if (my.print == T) {pdf(outpdf)}
  
  misc = list(
    p.plot = list(),
    p.name = c(),
    p.type = c(),
    sig.coords = data.frame()
  )
  
  # Loop through positionTypes
  for (positionTypesInd in seq(1,length(positionTypes))) {
    positionType      = positionTypes[positionTypesInd]
    
    positionTypePrint = if (positionType == 'orig') {''} else {positionType}
    positionType1     = paste(positionTypePrint,'beg',sep='')
    positionType2     = paste(positionTypePrint,'end',sep='')
    posTypes = c(positionType1,positionType2)
    
    misc[[positionType]] = list(
      positionType1 = positionType1,
      positionType2 = positionType2,
      p.plot = list(),
      p.name = c(),
      sig.coords = data.frame()
    )
    
    
    # Loop through posTypes
    for (posTypesInd in seq(1,length(posTypes))) {
      posType1 = posTypes[posTypesInd]
      posType2 = posType1
      if (grepl("beg",posType2)) {
        posType2 = gsub("beg","end",posType2,perl=T)
      } else if (grepl("end",posType2)) {
        posType2 = gsub("end","beg",posType2,perl=T)
      }
      
      df = df[order(df[,'cluster'],df[,posType1],df[,posType2]),]
      
      # Loop through testTypes
      for (testTypesInd in seq(1,length(testTypes))) {
        testType = testTypes[testTypesInd]
        
        misc[[positionType]][[testType]] = list(
          p.plot = list(),
          p.name = c(),
          sig.coords = data.frame(),
          last.i = 1,
          curr.i = 1
        )
        
        df2 = data.frame()
        
        clusters = my.order(unique(df$cluster))
        misc[[positionType]][[testType]]$clusters = clusters
        misc[[positionType]][[testType]]$p.type = c(misc$p.type,paste(posType1,'.',testType,'.pval',sep=''))
        
        p.type = misc[[positionType]][[testType]]$p.type
        
        df[,p.type] = 0
        
        # Loop through clusters
        for (clustersInd in seq(1,length(clusters))) {
          cluster = as.character(clusters[clustersInd])
          
          misc[[positionType]][[testType]][[cluster]] = list(
            p.plot = list(),
            p.name = c(),
            sig.coords = data.frame()
          )
          
          dm = df[df$cluster == cluster,]
          
          # Get original min and max y-axis for that cluster
          misc[[positionType]][[testType]][[cluster]]$min.y = min(dm$y)
          misc[[positionType]][[testType]][[cluster]]$max.y = max(dm$y)
          
          dm = dm[order(dm$cluster,dm[,posType1],dm[,posType2]),]; dm$y = seq(1,dim(dm)[1]);# = reorder_y(dm)
          
          # Get min and max i for that cluster
          min.i = 1
          max.i = dim(dm)[1]
          
          # Initialize p-value of p.type
          dm[,p.type] = 0
          
          
          # Get min y-axis for that cluster (to be added later)
          
          dm = dm[order(dm$cluster,dm[,posType1],dm[,posType2]),]; dm$y = seq(1,dim(dm)[1]);# = reorder_y(dm)
          
          dm.temp.draw   = dm
          dm.temp.draw$x = dm[,posType1]
          
          p = ggplot(dm.temp.draw,aes(x=x,y=y)) +
            geom_line(col='grey') +
            coord_cartesian(xlim=c(0,3000),ylim=c(0,max(dm$y))) +
            xlab(paste('R-loop',posType1,'position')) + ylab(pasta('Read',posType1,'Number')) +
            ggtitle(paste(testType,'cluster',cluster,'; minpval',gp$minpval))
          
          ##############
          # BIG LOOP 
          # Loop Start
          # - start at i = 1 until i = my.range; my.range = min(i+window or dim(dm)[1])
          # - check p-value by testType (unif, hots, norm, etc)
          #   - if p-value is bigger than pval threshold (gp$minpval), extend and increase pval threshold (+ 1% each iteration) 
          #     until it isn't, then next i jumps to my.range + 1
          #   - else, return and i = i + 1
          
          curr.i = 1
          last.i = 1
          
          while (curr.i < max.i) {
            
            # if verbose: if there was a jump bigger than N in curr.i, then print the positions
            if (verbose == T & last.i != curr.i - 1) {
              print(paste(positionType,testType,'cluster=',cluster,'last.i=',last.i,'curr.i=',curr.i))
            }
            
            last.i = curr.i
            
            my.range = list(
              i0  = curr.i,
              i1  = min(curr.i + gp$windowdist, max.i),
              add = 0,
              p   = 0,
              minpval = gp$minpval, #0.25
              minpval2 = gp$minpval2 #0.5
            )
            my.range$seq  = seq(my.range$i0, my.range$i1)
            my.range$size = length(my.range$seq)
            
            # iterate from i0 to i1 (see my.range above)
            # if p-value is bigger than pval threshold (gp$my.range$minpval), extend my.range and increase pval threshold (+ 1% each iteration)
            # else return curr.i = i + 1
            while (1) {
              if (verbose == T) {cat('  -',my.range$i0,my.range$i1,my.range$add,my.range$p,'\n')}
              my.range$p = 0
              dm2 = dm[my.range$seq,]
              
              
              # check each test
              if (size(dm[my.range$seq,]) >= gp$dist.test_min.length) {
                if (testType == 'unif') {
                  if (length(unique(dm[my.range$seq,][,posType1])) >= gp$dist.test_min.unique) {
                    my.range$hist = hist(dm[my.range$seq,posType1],plot=F,breaks = max(2,sqrt(dim(dm[my.range$seq,][,])[1])))
                    my.range$p = uniform.test(my.range$hist)$p.value
                  }
                } else if(testType == 'norm') {
                  if (length(unique(dm[my.range$seq,][,posType1])) >= gp$dist.test_min.unique) {
                    my.range$p = shapiro.test(dm[my.range$seq,posType1])$p.value
                  }
                }
                else if (testType == 'hots') {
                  if (length(unique(dm[my.range$seq,][,posType1])) < gp$dist.test_max.unique.hots) {
                    my.range$p = 1
                  }
                }
              }
              
              dm[my.range$seq,p.type] = apply(data.frame(dm[my.range$seq, p.type],rep(my.range$p, my.range$size)),1,max)
              
              my.range$minpval = min(1,my.range$minpval + 1/100)
              
              if (my.range$p >= my.range$minpval) {
                my.range$add = my.range$add + 1
              } else {
                break
              }
              
              my.range$i1 = min(my.range$i1 + my.range$add, max.i)
              my.range$seq = seq(my.range$i0, my.range$i1)
              my.range$size = length(my.range$seq)
              if (my.range$i1 >= max.i) {break}
            }
            
            # if the region passed pvalue threshold and extended by "my.range$add", then add that to the curr.i
            if (my.range$add > 0) {
              ibefore = my.range$i0
              curr.i = my.range$i1
              
              misc[[positionType]][[testType]][[cluster]]$sig.coords = rbind (
                misc[[positionType]][[testType]][[cluster]]$sig.coords,
                data.frame(positionType = posType1,
                           testType=testType,
                           cluster=cluster,
                           x0end=min(dm[seq(my.range$i0,my.range$i1),posType1]),
                           y0end=min(dm[seq(my.range$i0,my.range$i1),posType2]),
                           x1beg=max(dm[seq(my.range$i0,my.range$i1),posType1]),
                           y1beg=max(dm[seq(my.range$i0,my.range$i1),posType2]),
                           y0=my.range$i0,
                           y1=my.range$i1,
                           x=100,y=my.range$i0,xmin=100,xmax=150)
              )
              
              # Add to debug drawing
              # get color            
              my.range$col = dist.test_get.color(my.range$p,gp)
              dm.to_draw = data.frame(x=dm[my.range$seq,posType1] + gp$divby,y=dm[my.range$seq,]$y)
              p = p + geom_line(data=dm.to_draw,aes(x=x,y=y),color=my.range$col,lwd=0.5)
              
            }
            
            # if current i is at max i - 1, then break (can't test if i == max.i as there'll only be 1 point)
            if (curr.i >= max.i - 1) {
              
              # add rectangles to p
              p = p + geom_rect(
                data = misc[[positionType]][[testType]][[cluster]]$sig.coords,
                aes(xmin=xmin,xmax=xmax,ymin=y0,ymax=y1),
                fill='red4',color=NA
              )
              
              #append(
              #  misc[[positionType]][[testType]][[cluster]]p.plot,print(p))
              #append(misc[[positionType]][[testType]][[cluster]]$p.plot, print(p))
              misc[[positionType]][[testType]][[cluster]]$p.name = p.type
              
              break
            }
            
            curr.i = curr.i + 1
            
          }
          
          #
          ########################## Loop End
          
          df2 = rbind(df2,dm)
          misc[[positionType]][[testType]][[cluster]]$sig.coords$y0 = misc[[positionType]][[testType]][[cluster]]$sig.coords$y0 + misc[[positionType]][[testType]][[cluster]]$min.y
          misc[[positionType]][[testType]][[cluster]]$sig.coords$y1 = misc[[positionType]][[testType]][[cluster]]$sig.coords$y1 + misc[[positionType]][[testType]][[cluster]]$min.y
          misc[[positionType]][[testType]]$sig.coords = rbind(
            misc[[positionType]][[testType]]$sig.coords,
            misc[[positionType]][[testType]][[cluster]]$sig.coords
          )
        }
        
        if (dim(df)[1] != dim(df2)[1] | dim(df)[2] != dim(df2)[2]) {
          print(paste("Error! Size of df isn't same as size of df2!",paste(dim(df),collapse='_'),paste(dim(df2),collapse='_'),collapse='_'))
        }
        df2 = df2[order(df2[,'cluster'],df2[,posType1],df2[,posType2]),]; df2$y = seq(1,size(df2))
        
        #df[df2[,p.type] != 0,p.type] = df2[df2[,p.type] != 0,p.type]# = cbind(df,df2[,dim(df2)[2]]); colnames(df)[dim(df)[2]] = colnames(df2)[dim(df2)[2]]
        df[,p.type] = apply(data.frame(df[, p.type],df2[,p.type]),1,max)
        misc[[positionType]]$sig.coords = rbind(
          misc[[positionType]]$sig.coords,
          misc[[positionType]][[testType]]$sig.coords
        )
      }
    }
    misc$sig.coords = rbind(
      misc$sig.coords,
      misc[[positionType]]$sig.coords
    )
  }
  to_return = list(df=df, misc=misc)
  return(to_return)
}
#geom_histogram(aes(y=..density..,color=af(meanbeg.unif.sig.coords)),lwd=0.1,fill=NA,binwidth = 5) + facet_grid(cluster~.) +
#  geom_density(aes(color=af(cluster))) + theme_bw() + coord_cartesian(xlim=c(0,3000)) + facet_grid(meanbeg.unif.sig.coords~.) +
#  geom_density(aes(color=af(cluster))) + theme_bw() + coord_cartesian(xlim=c(0,3000)) + facet_grid(meanbeg.unif.sig.coords~.) +
#  geom_density(aes(color=af(cluster),fill=af(meanbeg.unif.sig.coords)),alpha=0.1) + theme_bw() + coord_cartesian(xlim=c(0,3000)) +
#  geom_histogram(aes(y=..density..,color=af(meanbeg.unif.sig.coords)),lwd=0.1,fill=NA,binwidth = 5) +
#  geom_density(aes(color=af(meanbeg.unif.sig.coords))) + theme_bw() + coord_cartesian(xlim=c(0,3000)) +

annot.pvar = function(my.df,my.sc=data.frame()) {
  my.df = my.df[order(my.df$cluster,my.df$meanbeg, my.df$meanend),]; my.df$y = seq(1,dim(my.df)[1])
  for (positionTypesInd in seq(1,length(positionTypes))) {
    positionType      = positionTypes[positionTypesInd]
    
    positionTypePrint = if (positionType == 'orig') {''} else {positionType}
    positionType1     = paste(positionTypePrint,'beg',sep='')
    positionType2     = paste(positionTypePrint,'end',sep='')
    
    posTypes = c(positionType1,positionType2)
    for (posTypesInd in seq(1,length(posTypes))) {
      posType1 = posTypes[posTypesInd]
      posType2 = posType1
      if (grepl("beg",posType2)) {
        posType2 = gsub("beg","end",posType2,perl=T)
      } else if (grepl("end",posType2)) {
        posType2 = gsub("end","beg",posType2,perl=T)
      }
      my.df = my.df[order(my.df[,'cluster'],my.df[,posType1],my.df[,posType2]),]
      my.df$y = seq(1,dim(my.df)[1])
      
      
      for (testTypesInd in seq(1,size(testTypes))) {
        testType = testTypes[testTypesInd]
        pval.varname = paste(posType1,'.',testType,'.pval',sep='')
        varname = paste(posType1,'.',testType,'.sig.coords',sep='')
        print(varname)
        my.df[,varname] = 0
        d1 = print(defined(my.sc[my.sc$positionType == posType1,]$testType))
        d2 = print(defined(my.sc[my.sc$positionType == posType1 & my.sc$testType == testType,]$testType))
        #if (d1 == TRUE & d2 == TRUE) {
        if (defined(my.sc[my.sc$positionType == posType1 & my.sc$testType == testType,]$testType)) {
          my.sc.sub = my.sc[my.sc$positionType == posType1 & my.sc$testType == testType,]
          for (i in seq(1,dim(my.sc.sub)[1])) {
            my.df[my.df$y >= my.sc.sub$y0[i] & my.df$y <= my.sc.sub$y1[i],varname] = 1
            #my.df[my.df[,pval.varname] < 0.5,varname] = 0
            
            if (testType == 'unif' & posType1 == 'meanbeg') {
              print(paste(posType1,testType,i,my.sc.sub$y0[i],my.sc.sub$y1[i]))
            }
          }
        }
      }
    }
  }
  return(my.df)  
}

re.sc = function(my.df,my.sc,positionTypes,testTypes) {#,posTypeWant = 'beg') {
  #my.df = my.list$df
  #my.sc = my.list$misc$sig.coords
  
  #my.df = my.df[order(my.df$cluster,my.df$meanbeg, my.df$meanend),]; my.df$y = seq(1,dim(my.df)[1])
  for (positionTypesInd in seq(1,length(positionTypes))) {
    positionType      = positionTypes[positionTypesInd]
    
    positionTypePrint = if (positionType == 'orig') {''} else {positionType}
    positionType1     = paste(positionTypePrint,'beg',sep='')
    positionType2     = paste(positionTypePrint,'end',sep='')
    
    #    my.df = my.df[order(my.df[,'cluster'],my.df[,positionType1], my.df[,positionType2]),]; my.df$y = seq(1,dim(my.df)[1])
    posTypes = c(positionType1,positionType2)
    for (posTypesInd in seq(1,length(posTypes))) {
      posType1 = posTypes[posTypesInd]
      posType2 = posType1
      if (grepl("beg",posType2)) {
        posType2 = gsub("beg","end",posType2,perl=T)
      } else if (grepl("end",posType2)) {
        posType2 = gsub("end","beg",posType2,perl=T)
      }
      my.df = my.df[order(my.df[,'cluster'],my.df[,posType1],my.df[,posType2]),]
      my.df$y = seq(1,dim(my.df)[1])
      
      
      for (testTypesInd in seq(1,size(testTypes))) {
        testType = testTypes[testTypesInd]
        pval.varname = paste(posType1,'.',testType,'.pval',sep='')
        varname = paste(posType1,'.',testType,'.sig.coords',sep='')
        print(varname)
        my.df[,varname] = 0
        d1 = print(defined(my.sc[my.sc$positionType == posType1,]$testType))
        d2 = print(defined(my.sc[my.sc$positionType == posType1 & my.sc$testType == testType,]$testType))
        #if (d1 == TRUE & d2 == TRUE) {
        if (defined(my.sc[my.sc$positionType == posType1 & my.sc$testType == testType,]$testType)) {
          my.sc.sub = my.sc[my.sc$positionType == posType1 & my.sc$testType == testType,]
          for (i in seq(1,dim(my.sc.sub)[1])) {
            my.df[my.df$y >= my.sc.sub$y0[i] & my.df$y <= my.sc.sub$y1[i],varname] = 1
            #my.df[my.df[,pval.varname] < 0.5,varname] = 0
            
            if (testType == 'unif' & posType1 == 'meanbeg') {
              print(paste(posType1,testType,i,my.sc.sub$y0[i],my.sc.sub$y1[i]))
            }
          }
        }
      }
    }
  }
  
  #  my.df = my.df[order(my.df$cluster,my.df$meanbeg, my.df$meanend),]; my.df$y = seq(1,dim(my.df)[1])
  for (positionTypesInd in seq(1,length(positionTypes))) {
    positionType      = positionTypes[positionTypesInd]
    
    positionTypePrint = if (positionType == 'orig') {''} else {positionType}
    positionType1     = paste(positionTypePrint,'beg',sep='')
    positionType2     = paste(positionTypePrint,'end',sep='')
    
    posTypes = c(positionType1,positionType2)
    for (posTypesInd in seq(1,length(posTypes))) {
      posType1 = posTypes[posTypesInd]
      posType2 = posType1
      if (grepl("beg",posType2)) {
        posType2 = gsub("beg","end",posType2,perl=T)
      } else if (grepl("end",posType2)) {
        posType2 = gsub("end","beg",posType2,perl=T)
      }
      my.df = my.df[order(my.df[,'cluster'],my.df[,posType1],my.df[,posType2]),]
      my.df$y = seq(1,dim(my.df)[1])
      
      
      pval.varname.best = paste(posType1,'.pval.best',sep='')
      varname.best = paste(posType1,'sig.coords.best',sep='')
      testdm = data.frame()
      varnames = c()
      #   testdm.name = data.frame()
      testInd = 1
      for (testTypesInd in seq(1,size(testTypes))) {
        testType = testTypes[testTypesInd]
        pval.varname = paste(posType1,'.',testType,'.pval',sep='')
        varname = paste(posType1,'.',testType,'.sig.coords',sep='')
        print(varname)
        if (defined(my.df[,pval.varname])) {
          varnames = c(varnames,varname)
          if (testInd == 1) {
            testdm = data.frame(varname=my.df[,pval.varname])
          } else {
            testdm = cbind(testdm,data.frame(varname=my.df[,pval.varname]))
          }
          testInd = testInd + 1
          
          #      testdm.name = rbind(testdm.name,varname)
        }
      }
      colnames(testdm) = paste(posType1,'.',testTypes,'.sig.coords',sep='')
      if (defined(testdm)) {
        my.df[,pval.varname.best] = apply(testdm,1,max)
        my.df[,varname.best] = apply(testdm,1,function(x) {varnames[x == max(x)][2]})
      }    
      
    }
  }
  
  for (positionTypesInd in seq(1,length(positionTypes))) {
    positionType      = positionTypes[positionTypesInd]
    
    positionTypePrint = if (positionType == 'orig') {''} else {positionType}
    positionType1     = paste(positionTypePrint,'beg',sep='')
    positionType2     = paste(positionTypePrint,'end',sep='')
    
    posTypes = c(positionType1,positionType2)
    for (posTypesInd in seq(1,length(posTypes))) {
      posType1 = posTypes[posTypesInd]
      posType2 = posType1
      if (grepl("beg",posType2)) {
        posType2 = gsub("beg","end",posType2,perl=T)
      } else if (grepl("end",posType2)) {
        posType2 = gsub("end","beg",posType2,perl=T)
      }
      my.df = my.df[order(my.df[,'cluster'],my.df[,posType1],my.df[,posType2]),]
      my.df$y = seq(1,dim(my.df)[1])
      
      
      pval.varname.best = paste(posType1,'.pval.best',sep='')
      varname.best = paste(posType1,'sig.coords.best',sep='')
      #   testdm.name = data.frame()
      for (testTypesInd in seq(1,size(testTypes))) {
        testType = testTypes[testTypesInd]
        pval.varname = paste(posType1,'.',testType,'.pval',sep='')
        varname = paste(posType1,'.',testType,'.sig.coords',sep='')
        my.df[my.df[,pval.varname] != my.df[,pval.varname.best],varname] = 0
        my.df[my.df[,pval.varname] != my.df[,pval.varname.best],pval.varname] = 0
      }
    }
  }
  
  head(my.df)
  
  #my.scbackup = my.sc
  #my.sc = my.scbackup
  my.sc = data.frame()
  
  #my.df = my.df[order(my.df$cluster,my.df$meanbeg, my.df$meanend),]; my.df$y = seq(1,dim(my.df)[1])
  
  for (positionTypesInd in seq(1,length(positionTypes))) {
    positionType      = positionTypes[positionTypesInd]
    
    positionTypePrint = if (positionType == 'orig') {''} else {positionType}
    positionType1     = paste(positionTypePrint,'beg',sep='')
    positionType2     = paste(positionTypePrint,'end',sep='')
    
    #    my.df = my.df[order(my.df[,'cluster'],my.df[,positionType1], my.df[,positionType2]),]; my.df$y = seq(1,dim(my.df)[1])
    posTypes = c(positionType1,positionType2)
    for (posTypesInd in seq(1,length(posTypes))) {
      posType1 = posTypes[posTypesInd]
      posType2 = posType1
      if (grepl("beg",posType2)) {
        posType2 = gsub("beg","end",posType2,perl=T)
      } else if (grepl("end",posType2)) {
        posType2 = gsub("end","beg",posType2,perl=T)
      }
      my.df = my.df[order(my.df[,'cluster'],my.df[,posType1],my.df[,posType2]),]
      my.df$y = seq(1,dim(my.df)[1])
      
      pval.varname.best = paste(posType1,'.pval.best',sep='')
      varname.best = paste(posType1,'sig.coords.best',sep='')
      #   testdm.name = data.frame()
      my.Xmin = 0
      my.Xmax = 25
      for (testTypesInd in seq(1,size(testTypes))) {
        print(paste(posTypesInd,posTypes[posTypesInd],testTypesInd,testTypes[testTypesInd]))
        testType = testTypes[testTypesInd]
        pval.varname = paste(posType1,'.',testType,'.pval',sep='')
        varname = paste(posType1,'.',testType,'.sig.coords',sep='')
        
        if (defined(my.df[my.df[,varname] != 0,])) {
          test6 = my.df[my.df[,varname] != 0,]
          #test6 = my.df[my.df$meanbeg.unif.sig.coords != 0,]
          #[test6$meanbeg.unif.sig.coords != 0,]
          my.edge6 = findedge(test6$y)
          my.seq6 = data.frame()
          print('here3')
          for (i in 1:(size(my.edge6)-1)) {
            my.seq6 = rbind(
              my.seq6,
              data.frame(
                positionType = posType1,
                testType = testType,
                cluster = test6[my.edge6[i]+1,]$cluster,
                x0end=test6[my.edge6[i]+1,posType1],
                y0end=test6[my.edge6[i]+1,posType2],
                x1beg=test6[my.edge6[i+1],posType1],
                y1beg=test6[my.edge6[i+1],posType2],
                y0 = test6[my.edge6[i]+1,]$y,
                y1 = test6[my.edge6[i+1],]$y,
                xmin = 0, xmax=25, x = 0, y = 0
              )
              
            )
            print('here4')
            i = i + 1
          }
          print('here5')
          my.sc = rbind(my.sc,my.seq6)
        }
        print('here6')
        my.Xmin = my.Xmin + 50
        my.Xmax = my.Xmax + 50
      }
      print('here7')
      
    }
    print('here8')
  }
  
  return(my.sc)
}

CLUSTFILES  = my.order(paste('resources/misc/',dir("./resources/misc/","*final*.RDS"),sep=''))
PEAKFILES  = my.order(paste('resources/peaks/',dir("./resources/peaks/","*.BED"),sep=''))
BEDFILES   = my.order(paste('resources/bed/',dir("./resources/bed/","*\\.bed$"),sep=''))
FASTAFILES = my.order(paste('resources/fa/',dir("./resources/fa/","*.fa"),sep=''))
BIGFASTAFILES = my.order(paste('resources/bigfa/',dir("./resources/bigfa/","bigfa.fa"),sep=''))

#parseBEDFile(debug=T)
#parsePEAKFile(debug=T)
#parseFASTAFile(debug=T)
#parseCLUSTFile(debug=T)


