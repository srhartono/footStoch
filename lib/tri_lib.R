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

gp = list(
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
  divby=25
)


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

get_title = function(mytitle='',myparams=list(),verbose=F,debug=F) {
  
  if (length(myparams) == 0) {
    
    return(mytitle)
  }
  for (i in 1:length(myparams))  {
    
    paramname = names(myparams)[i]
    paramwant = myparams[i][[1]]
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

slice_bed = function(mybeds=BEDS,myparams=list(),myparams_regex=c(),genewant = NA,verbose=F,debug=F) {
  if (is.na(genewant)) {
    if (defined(myparams_regex$gene)) {
      if (myparams_regex$gene == TRUE) {
        mybed = mybeds[grep(myparams$gene,mybeds$gene,ignore.case = TRUE),]
      } else {
        mybed = mybeds[mybeds$gene == myparams$gene,]
      }
    } else {
      mybed = mybeds[mybeds$gene == myparams$gene,]
    }
    if (defined(mybed) == FALSE) {
      mybed = mybeds[grep(paste('^',myparams$gene,'$',sep=''),mybeds$gene,ignore.case = TRUE),]
    }
  } else {
    mybed = mybeds[grep(paste('^',genewant,'$',sep=''),mybeds$gene),]
    if (defined(mybed) == FALSE) {
      mybed = mybeds[grep(paste('^',genewant,'$',sep=''),mybeds$gene,ignore.case = TRUE),]
    }
  }
  return(mybed)
}

slice_df = function(mydata,myparams=list(),myparams_regex=list(),verbose=F,debug=F) {

  #gene.want='any',treat.want='any',peaktype.want='any',VR.want='any',thres.want='any',verbose=F,debug=F)
  if (length(myparams) == 0) {

    return(mydata)
  }
  for (i in 1:length(myparams))  {

    orig = mydata
    paramname = names(myparams)[i]
    paramwant = myparams[i][[1]]
    isregex = myparams_regex[i]
    if (isregex == TRUE) {

      cat(i,'. ',paramname,': ',paramwant,' (regex): \t',sep='')
    }
    else {

      cat(i,'. ',paramname,': ',paramwant,'\t',sep='')
    }
    
    if (paramwant == 'any') {

      cat(': From',dim(orig)[1],'x',dim(orig)[2],'to',dim(mydata)[1],'x',dim(mydata)[2],'\n',sep=' ')
      next
    }

    mydata = mydata[!is.na(mydata[,paramname]),]
    if (isregex == FALSE) {

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

do_distributionTests = function(df = data.frame(), positionTypes = c('mean','orig'), testTypes=c('unif','norm','hots'), gp=list(), outpdf=F,myprint=F, debug=F, verbose=F) {
  
  #  if (defined(gp$minpval) == F) {
  #    gp$minpval = 0.05
  #  }
  
  do_distributionTests.get_outpdf = function(positionTypes=positionTypes,testTypes=testTypes,gp=gp,outpdf=outpdf) {
    if (outpdf == F) {
      outpdf = paste(paste(positionTypes,sep='_',collapse='_'),'_',paste(testTypes,sep='_',collapse='_'),'.pdf',collapse='',sep='')
      print(outpdf)
      if (defined(gp$mytitle)) {
        outpdf = paste(gp$mytitle,'_',outpdf,sep='')
      }
      to_return = outpdf
      return(to_return)
    }
  }
  
  outpdf = do_distributionTests.get_outpdf(outpdf=outpdf,gp=gp,positionTypes = positionTypes,testTypes=testTypes)
  
  if (verbose == T) {print(outpdf)}
  
  if (myprint == T) {pdf(outpdf)}
  
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
        
        clusters = myorder(unique(df$cluster))
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
          # - start at i = 1 until i = myrange; myrange = min(i+window or dim(dm)[1])
          # - check p-value by testType (unif, hots, norm, etc)
          #   - if p-value is bigger than pval threshold (gp$minpval), extend and increase pval threshold (+ 1% each iteration) 
          #     until it isn't, then next i jumps to myrange + 1
          #   - else, return and i = i + 1
          
          curr.i = 1
          last.i = 1
          
          while (curr.i < max.i) {
            
            # if verbose: if there was a jump bigger than N in curr.i, then print the positions
            if (verbose == T & last.i != curr.i - 1) {
              print(paste(positionType,testType,'cluster=',cluster,'last.i=',last.i,'curr.i=',curr.i))
            }
            
            last.i = curr.i
            
            myrange = list(
              i0  = curr.i,
              i1  = min(curr.i + gp$windowdist, max.i),
              add = 0,
              p   = 0,
              minpval = gp$minpval, #0.25
              minpval2 = gp$minpval2 #0.5
            )
            myrange$seq  = seq(myrange$i0, myrange$i1)
            myrange$size = length(myrange$seq)
            
            # iterate from i0 to i1 (see myrange above)
            # if p-value is bigger than pval threshold (gp$myrange$minpval), extend myrange and increase pval threshold (+ 1% each iteration)
            # else return curr.i = i + 1
            while (1) {
              if (verbose == T) {cat('  -',myrange$i0,myrange$i1,myrange$add,myrange$p,'\n')}
              myrange$p = 0
              dm2 = dm[myrange$seq,]
              
              
              # check each test
              if (size(dm[myrange$seq,]) >= gp$dist.test_min.length) {
                if (testType == 'unif') {
                  if (length(unique(dm[myrange$seq,][,posType1])) >= gp$dist.test_min.unique) {
                    myrange$hist = hist(dm[myrange$seq,posType1],plot=F,breaks = max(2,sqrt(dim(dm[myrange$seq,][,])[1])))
                    myrange$p = uniform.test(myrange$hist)$p.value
                  }
                } else if(testType == 'norm') {
                  if (length(unique(dm[myrange$seq,][,posType1])) >= gp$dist.test_min.unique) {
                    myrange$p = shapiro.test(dm[myrange$seq,posType1])$p.value
                  }
                }
                else if (testType == 'hots') {
                  if (length(unique(dm[myrange$seq,][,posType1])) < gp$dist.test_max.unique.hots) {
                    myrange$p = 1
                  }
                }
              }
              
              dm[myrange$seq,p.type] = apply(data.frame(dm[myrange$seq, p.type],rep(myrange$p, myrange$size)),1,max)
              
              myrange$minpval = min(1,myrange$minpval + 1/100)
              
              if (myrange$p >= myrange$minpval) {
                myrange$add = myrange$add + 1
              } else {
                break
              }
              
              myrange$i1 = min(myrange$i1 + myrange$add, max.i)
              myrange$seq = seq(myrange$i0, myrange$i1)
              myrange$size = length(myrange$seq)
              if (myrange$i1 >= max.i) {break}
            }
            
            # if the region passed pvalue threshold and extended by "myrange$add", then add that to the curr.i
            if (myrange$add > 0) {
              ibefore = myrange$i0
              curr.i = myrange$i1
              
              misc[[positionType]][[testType]][[cluster]]$sig.coords = rbind (
                misc[[positionType]][[testType]][[cluster]]$sig.coords,
                data.frame(positionType = posType1,
                           testType=testType,
                           cluster=cluster,
                           x0end=min(dm[seq(myrange$i0,myrange$i1),posType1]),
                           y0end=min(dm[seq(myrange$i0,myrange$i1),posType2]),
                           x1beg=max(dm[seq(myrange$i0,myrange$i1),posType1]),
                           y1beg=max(dm[seq(myrange$i0,myrange$i1),posType2]),
                           y0=myrange$i0,
                           y1=myrange$i1,
                           x=100,y=myrange$i0,xmin=100,xmax=150)
              )

              # Add to debug drawing
              # get color            
              myrange$col = dist.test_get.color(myrange$p,gp)
              dm.to_draw = data.frame(x=dm[myrange$seq,posType1] + gp$divby,y=dm[myrange$seq,]$y)
              p = p + geom_line(data=dm.to_draw,aes(x=x,y=y),color=myrange$col,lwd=0.5)
             
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

annot.pvar = function(mydf,mysc=data.frame()) {
  mydf = mydf[order(mydf$cluster,mydf$meanbeg, mydf$meanend),]; mydf$y = seq(1,dim(mydf)[1])
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
      mydf = mydf[order(mydf[,'cluster'],mydf[,posType1],mydf[,posType2]),]
      mydf$y = seq(1,dim(mydf)[1])
      
      
      for (testTypesInd in seq(1,size(testTypes))) {
        testType = testTypes[testTypesInd]
        pval.varname = paste(posType1,'.',testType,'.pval',sep='')
        varname = paste(posType1,'.',testType,'.sig.coords',sep='')
        print(varname)
        mydf[,varname] = 0
        d1 = print(defined(mysc[mysc$positionType == posType1,]$testType))
        d2 = print(defined(mysc[mysc$positionType == posType1 & mysc$testType == testType,]$testType))
        #if (d1 == TRUE & d2 == TRUE) {
        if (defined(mysc[mysc$positionType == posType1 & mysc$testType == testType,]$testType)) {
          mysc.sub = mysc[mysc$positionType == posType1 & mysc$testType == testType,]
          for (i in seq(1,dim(mysc.sub)[1])) {
            mydf[mydf$y >= mysc.sub$y0[i] & mydf$y <= mysc.sub$y1[i],varname] = 1
            #mydf[mydf[,pval.varname] < 0.5,varname] = 0
            
            if (testType == 'unif' & posType1 == 'meanbeg') {
              print(paste(posType1,testType,i,mysc.sub$y0[i],mysc.sub$y1[i]))
            }
          }
        }
      }
    }
  }
  return(mydf)  
}

re.sc = function(mydf,mysc=data.frame(),positionTypes,testTypes) {#,posTypeWant = 'beg') {
  #mydf = mylist$df
  #mysc = mylist$misc$sig.coords
  
  #mydf = mydf[order(mydf$cluster,mydf$meanbeg, mydf$meanend),]; mydf$y = seq(1,dim(mydf)[1])
  for (positionTypesInd in seq(1,length(positionTypes))) {
    positionType      = positionTypes[positionTypesInd]
    
    positionTypePrint = if (positionType == 'orig') {''} else {positionType}
    positionType1     = paste(positionTypePrint,'beg',sep='')
    positionType2     = paste(positionTypePrint,'end',sep='')
    
#    mydf = mydf[order(mydf[,'cluster'],mydf[,positionType1], mydf[,positionType2]),]; mydf$y = seq(1,dim(mydf)[1])
    posTypes = c(positionType1,positionType2)
    for (posTypesInd in seq(1,length(posTypes))) {
      posType1 = posTypes[posTypesInd]
      posType2 = posType1
      if (grepl("beg",posType2)) {
        posType2 = gsub("beg","end",posType2,perl=T)
      } else if (grepl("end",posType2)) {
        posType2 = gsub("end","beg",posType2,perl=T)
      }
      mydf = mydf[order(mydf[,'cluster'],mydf[,posType1],mydf[,posType2]),]
      mydf$y = seq(1,dim(mydf)[1])
      
      
      for (testTypesInd in seq(1,size(testTypes))) {
        testType = testTypes[testTypesInd]
        pval.varname = paste(posType1,'.',testType,'.pval',sep='')
        varname = paste(posType1,'.',testType,'.sig.coords',sep='')
        print(varname)
        mydf[,varname] = 0
        d1 = print(defined(mysc[mysc$positionType == posType1,]$testType))
        d2 = print(defined(mysc[mysc$positionType == posType1 & mysc$testType == testType,]$testType))
        #if (d1 == TRUE & d2 == TRUE) {
        if (defined(mysc[mysc$positionType == posType1 & mysc$testType == testType,]$testType)) {
          mysc.sub = mysc[mysc$positionType == posType1 & mysc$testType == testType,]
          for (i in seq(1,dim(mysc.sub)[1])) {
            mydf[mydf$y >= mysc.sub$y0[i] & mydf$y <= mysc.sub$y1[i],varname] = 1
            #mydf[mydf[,pval.varname] < 0.5,varname] = 0
            
            if (testType == 'unif' & posType1 == 'meanbeg') {
              print(paste(posType1,testType,i,mysc.sub$y0[i],mysc.sub$y1[i]))
            }
          }
        }
      }
    }
  }
  
#  mydf = mydf[order(mydf$cluster,mydf$meanbeg, mydf$meanend),]; mydf$y = seq(1,dim(mydf)[1])
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
      mydf = mydf[order(mydf[,'cluster'],mydf[,posType1],mydf[,posType2]),]
      mydf$y = seq(1,dim(mydf)[1])
      
      
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
        if (defined(mydf[,pval.varname])) {
          varnames = c(varnames,varname)
          if (testInd == 1) {
            testdm = data.frame(varname=mydf[,pval.varname])
          } else {
            testdm = cbind(testdm,data.frame(varname=mydf[,pval.varname]))
          }
          testInd = testInd + 1

        #      testdm.name = rbind(testdm.name,varname)
        }
      }
      colnames(testdm) = paste(posType1,'.',testTypes,'.sig.coords',sep='')
      if (defined(testdm)) {
        mydf[,pval.varname.best] = apply(testdm,1,max)
        mydf[,varname.best] = apply(testdm,1,function(x) {varnames[x == max(x)][2]})
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
      mydf = mydf[order(mydf[,'cluster'],mydf[,posType1],mydf[,posType2]),]
      mydf$y = seq(1,dim(mydf)[1])
      
      
      pval.varname.best = paste(posType1,'.pval.best',sep='')
      varname.best = paste(posType1,'sig.coords.best',sep='')
      #   testdm.name = data.frame()
      for (testTypesInd in seq(1,size(testTypes))) {
        testType = testTypes[testTypesInd]
        pval.varname = paste(posType1,'.',testType,'.pval',sep='')
        varname = paste(posType1,'.',testType,'.sig.coords',sep='')
        mydf[mydf[,pval.varname] != mydf[,pval.varname.best],varname] = 0
        mydf[mydf[,pval.varname] != mydf[,pval.varname.best],pval.varname] = 0
      }
    }
  }
  
  head(mydf)
  findedge = function(x) {
    
    # Find the edges
    edges <- which(abs(diff(x)) > 1)
    if (!0 %in% edges) {
      edges = c(0,edges)
    }
    if (!size(x) %in% edges) {
      edges = c(edges,size(x))
    }
    return(edges)
  }
  
  #myscbackup = mysc
  #mysc = myscbackup
  mysc = data.frame()
  
  #mydf = mydf[order(mydf$cluster,mydf$meanbeg, mydf$meanend),]; mydf$y = seq(1,dim(mydf)[1])
  
  for (positionTypesInd in seq(1,length(positionTypes))) {
    positionType      = positionTypes[positionTypesInd]
    
    positionTypePrint = if (positionType == 'orig') {''} else {positionType}
    positionType1     = paste(positionTypePrint,'beg',sep='')
    positionType2     = paste(positionTypePrint,'end',sep='')
    
#    mydf = mydf[order(mydf[,'cluster'],mydf[,positionType1], mydf[,positionType2]),]; mydf$y = seq(1,dim(mydf)[1])
    posTypes = c(positionType1,positionType2)
    for (posTypesInd in seq(1,length(posTypes))) {
      posType1 = posTypes[posTypesInd]
      posType2 = posType1
      if (grepl("beg",posType2)) {
        posType2 = gsub("beg","end",posType2,perl=T)
      } else if (grepl("end",posType2)) {
        posType2 = gsub("end","beg",posType2,perl=T)
      }
      mydf = mydf[order(mydf[,'cluster'],mydf[,posType1],mydf[,posType2]),]
      mydf$y = seq(1,dim(mydf)[1])
      
      pval.varname.best = paste(posType1,'.pval.best',sep='')
      varname.best = paste(posType1,'sig.coords.best',sep='')
      #   testdm.name = data.frame()
      myXmin = 0
      myXmax = 25
      for (testTypesInd in seq(1,size(testTypes))) {
        print(paste(posTypesInd,posTypes[posTypesInd],testTypesInd,testTypes[testTypesInd]))
        testType = testTypes[testTypesInd]
        pval.varname = paste(posType1,'.',testType,'.pval',sep='')
        varname = paste(posType1,'.',testType,'.sig.coords',sep='')
        
        if (defined(mydf[mydf[,varname] != 0,])) {
          test6 = mydf[mydf[,varname] != 0,]
        #test6 = mydf[mydf$meanbeg.unif.sig.coords != 0,]
        #[test6$meanbeg.unif.sig.coords != 0,]
          myedge6 = findedge(test6$y)
          myseq6 = data.frame()
          print('here3')
          for (i in 1:(size(myedge6)-1)) {
            myseq6 = rbind(
              myseq6,
              data.frame(
                positionType = posType1,
                testType = testType,
                cluster = test6[myedge6[i]+1,]$cluster,
                x0end=test6[myedge6[i]+1,posType1],
                y0end=test6[myedge6[i]+1,posType2],
                x1beg=test6[myedge6[i+1],posType1],
                y1beg=test6[myedge6[i+1],posType2],
                y0 = test6[myedge6[i]+1,]$y,
                y1 = test6[myedge6[i+1],]$y,
                xmin = 0, xmax=25, x = 0, y = 0
              )
              
            )
            print('here4')
            i = i + 1
          }
          print('here5')
          mysc = rbind(mysc,myseq6)
        }
        print('here6')
        myXmin = myXmin + 50
        myXmax = myXmax + 50
      }
      print('here7')

    }
    print('here8')
  }

  return(mysc)
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

