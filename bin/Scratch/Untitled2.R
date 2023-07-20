# # #    } else {
# # #      lines(dftemp[,vartest1][c(dim(dftemp)[1]-1,dim(dftemp)[1])],dftemp$y[c(dim(dftemp)[1]-1,dim(dftemp)[1])],col='red4',lwd=3)
# # #      rect(xleft = 1,ybottom = dftemp$y[dim(dftemp)[1]-1],xright = 100,ytop = dftemp$y[dim(dftemp)[1]],col='red4',border=NA)
# # #    }
# # 
# # } else {
# #   if (last == 'unif' | last == 'hot')
# #   {
# #     i = i + gp$windowdist
# #   }
# #   last = 'none'
# # }
# 
# if (i > 1) {
#   meanbeggp$xlines = dftemp[,vartest1][seq(dim(dftemp)[1]-1,dim(dftemp)[1])]
#   meanbeggp$ylines = dftemp$y[seq(dim(dftemp)[1]-1,dim(dftemp)[1])]
#   meanbeggp$yrectmins = dftemp$y[dim(dftemp)[1]-1]
#   meanbeggp$yrectmaxs = dftemp$y[dim(dftemp)[1]]
# } else {
#   meanbeggp$xlines = dftemp[,vartest1][seq(1,dim(dftemp)[1])]
#   meanbeggp$ylines = dftemp$y[seq(1,dim(dftemp)[1])]
#   meanbeggp$yrectmins = dftemp$y[1]
#   meanbeggp$yrectmaxs = dftemp$y[dim(dftemp)[1]]
# }



#dm1file = paste(resourcesDir,dir(resourcesDir,'dm1.all.RDS'),sep='')
# for (i in 1:length(RDS)) {
#   temp = readRDS(RDS[i])
#   print(i)
#   if (i == 1) {
#     final0 = temp
#   } else {
#     final0 = rbind(final0,temp)
#   }
# }
# head(final0)
# head(dm1)
# saveRDS(dm1,file="dm1.all.RDS")
# saveRDS(final0,file='final3.all.RDS')


---------
  vartest1 = 'meanbeg'
  myvar = paste(vartest1,'.unif',sep='')
  gp$distbetween = 5
  df[,myvar] = 0
  for (i in seq(1,(max(df$y)),gp$stepdist)) {
    if (i == 1) {
      lastunif = 1
      lasthots = 1
      if (grepl('beg',vartest1)) {vartest2 = gsub('beg','end',vartest1,perl=T)} else {vartest2 = gsub('end','beg',vartest1,perl=T)}
      df = df[order(df$cluster,df[,vartest1],df[,vartest2]),]; df = reorder_y(df)
    }
    myrange = min(i+gp$windowdist,max(df$y))
    seq.myrange  = seq(i,min(i+gp$windowdist,max(df$y)))
    
    
    if (i == 1) {
      plot(df[,vartest1],df$y,type='l',col='grey',xlim = c(0,2500),ylim=c(0,max(df$y)))
      last = 'none'
    }
    
    
    dftemp = df[seq.myrange,]
    meanbeghist = hist(dftemp[,vartest1],plot=F,breaks = max(2,dim(dftemp)[1]/10))
    pvalhots = 0
    if (length(unique(dftemp[,vartest1])) >= 2) {
      pvalmeanbeg = uniform.test(meanbeghist)$p.value
    } else  {
      print(paste(i,'- is a hotspot!'))
      pvalhots = 1
      pvalmeanbeg = 0
    }
    
    unif.bin.cols = brewer.pal(9,"Reds")
    unif.bin.brks = seq(0.05,1,length.out=9); unif.bin.brks.diff = unif.bin.brks[2]-unif.bin.brks[1]
    unif.col = unif.bin.cols[unif.bin.brks < pvalmeanbeg & unif.bin.brks+unif.bin.brks.diff > pvalmeanbeg]
    if (defined(unif.col) == FALSE) {if(pvalmeanbeg > 0.8) {unif.col = unif.bin.cols[length(unif.bin.cols)]} else {unif.col = rgb(0,0,0,0)}}
    unif.col=unif.bin.cols[9]
    #  dftemp.myrange = max(1,dim(dftemp)[1]-gp$windowdist/2)
    dftemp.unif.myrange = max(1,dim(dftemp)[1])#*pvalmeanbeg)
    
    hots.bin.cols = brewer.pal(9,"Greens")
    hots.bin.brks = seq(0.05,1,length.out=9); hots.bin.brks.diff = hots.bin.brks[2]-hots.bin.brks[1]
    hots.col = hots.bin.cols[hots.bin.brks < pvalhots & hots.bin.brks+hots.bin.brks.diff > pvalhots]
    if (defined(hots.col) == FALSE) {if (pvalhots > 0.8) {hots.col = hots.bin.cols[length(hots.bin.cols)]} else {hots.col = rgb(0,0,0,0)}}
    dftemp.hots.myrange = max(1,dim(dftemp)[1]*pvalhots)
    
    meanbeggp = list()
    meanbeggp$xrectadds = 100
    meanbeggp$xrectmins = 1
    meanbeggp$xrectmaxs = meanbeggp$xrectmins + meanbeggp$xrectadds
    if (pvalmeanbeg > 0.05) {
      # if (i == 1 | lastunif >= gp$distbetween | lastunif == -1) {
      #   meanbeggp$xlines = dftemp[,vartest1][seq(1,dftemp.unif.myrange)]
      #   meanbeggp$ylines = dftemp$y[seq(1,dftemp.unif.myrange)]
      #   meanbeggp$yrectmins = dftemp$y[1]
      #   meanbeggp$yrectmaxs = dftemp$y[dftemp.unif.myrange]
      # } else if (lastunif < gp$distbetween) {
      #  meanbeggp$xlines = dftemp[,vartest1][seq(dim(dftemp)[1]- 1 -lastunif,dim(dftemp)[1])]
      #  meanbeggp$ylines = dftemp$y[seq(dim(dftemp)[1]-1-lastunif,dim(dftemp)[1])]
      #  meanbeggp$yrectmins = dftemp$y[dim(dftemp)[1]-1-lastunif]
      #  meanbeggp$yrectmaxs = dftemp$y[dim(dftemp)[1]]
      # }
      meanbeggp$xlines = dftemp[,vartest1][seq(1,2)]
      meanbeggp$ylines = dftemp$y[seq(1,2)]
      meanbeggp$yrectmins = dftemp$y[1]
      meanbeggp$yrectmaxs = dftemp$y[2]
      df[i,myvar] = pvalmeanbeg
      # meanbeggp$xlines = dftemp[,vartest1][seq(dim(dftemp)[1]- 1,dim(dftemp)[1])]
      # meanbeggp$ylines = dftemp$y[seq(dim(dftemp)[1]-1,dim(dftemp)[1])]
      # meanbeggp$yrectmins = dftemp$y[dim(dftemp)[1]-1]
      # meanbeggp$yrectmaxs = dftemp$y[dim(dftemp)[1]]
      
      # meanbeggp$xlines = dftemp[,vartest1][seq(1,dftemp.unif.myrange)]
      #meanbeggp$ylines = dftemp$y[seq(1,dftemp.unif.myrange)]
      #meanbeggp$yrectmins = dftemp$y[1]
      #meanbeggp$yrectmaxs = dftemp$y[dftemp.unif.myrange]
      lines(meanbeggp$xlines+25,meanbeggp$ylines,col=unif.col,lwd=1)
      rect(meanbeggp$xrectmins,meanbeggp$yrectmins,meanbeggp$xrectmaxs,meanbeggp$yrectmaxs,col=unif.col,lwd=1,border=NA)
      lastunif = 0
    } else {
      lastunif = lastunif + 1
    }
    if (pvalhots > 0.05) {
      meanbeggp$xlines = dftemp[,vartest1][seq(1,dftemp.hots.myrange)]
      meanbeggp$ylines = dftemp$y[seq(1,dftemp.hots.myrange)]
      meanbeggp$yrectmins = dftemp$y[1]
      meanbeggp$yrectmaxs = dftemp$y[dftemp.hots.myrange]
      lines(meanbeggp$xlines+25,meanbeggp$ylines,col=hots.col,lwd=4)
      rect(meanbeggp$xrectmins+meanbeggp$xrectadds,meanbeggp$yrectmins,meanbeggp$xrectmaxs+meanbeggp$xrectadds,meanbeggp$yrectmaxs,col=hots.col,lwd=4,border=NA)
      lasthots = 1
    } else {
      lasthots = 0
    }
    
    if (lastunif > gp$distbetween) {
      i = i + gp$windowdist
      lastunif = -1
    }
  }
  
  #-------
  plot(df[,vartest1],df$y,type='l',col='grey',xlim = c(0,2500),ylim=c(0,max(df$y)))
  points(df[df[,myvar] > 0.05,vartest1]+10,df[df[,myvar] > 0.05,]$y,pch='.')
  myvar2 = paste(myvar,2,sep='')
  df[,myvar2] = 0
  for (i in 1:max(df$y)) {
    if (df[i,myvar] > 0.05) {
      pos2 = min(dim(df)[1],i + gp$windowdist)
      df[seq(i,pos2),myvar2] = 1
    }
  }
  #df = df[!is.na(df$beg),]
  for (i in 1:max(df$y)) {
    if (df[i,myvar] > 0.05) {
      if (i > 2) {
        df[seq(i-2,i),myvar2] = 0
      }
    }
  }
  points(df[df[,myvar2] != 0,vartest1]+50,df[df[,myvar2] != 0,]$y,pch='.',col='red4')
  
  
  #   # if (i == 1)
  #   # {
  #   #   #plot(y=df[,vartest1],x=df$y,type='l',col='grey',xlim = c(0,max(df$y)),ylim=c(0,max(df[,vartest1])))
  #   # }
  #   #else
  #   #{
  #   #  #lines(y=dftemp[,vartest1],x=dftemp$y,col='grey')
  #   #}
  #   print(i)
  # #  points(y=dftemp[,vartest1][dim(dftemp)[1]],x=dftemp$y[dim(dftemp)[1]],col='black',pch=15,cex=0.5)
  #   p.hot = 0
  #   p.unif = 0
  #   p.shap = 0
  #   if (length(unique(dftemp[,vartest1])) <= 1)
  #   {
  #     p.unif = 0
  #     p.hot = 1
  #     p.shap = 0
  #   }
  #   else
  #   {
  #     p.unif = uniform.test(hist(dftemp[,vartest1],plot=F))$p.value
  #     p.shap = shapiro.test(dftemp[,vartest1])$p.value
  #     p.hot = 0
  #   }
  #   if (p.hot > 0.05)
  #   {
  #     lines(y=dftemp[,vartest1],x=dftemp$y,col='red',lwd=4)
  #     rect(min(dftemp$y),100,min(dftemp$y)+10,120,col='red',border=NA)
  #   }
  #   else if (p.shap > 0.95)
  #   {
  #     lines(y=dftemp[,vartest1],x=dftemp$y,col='blue4',lwd=4)
  #     rect(min(dftemp$y),150,min(dftemp$y)+10,170,col='blue4',border=NA)
  #   }
  #   else if (p.unif > 0.99)
  #   {
  # #    lines(y=dftemp[,vartest1],x=dftemp$y,col='black',lwd=4)
  #     lines(y=dftemp[,vartest1],x=dftemp$y,col='green3',lwd=4)
  #     rect(min(dftemp$y),50,min(dftemp$y)+10,70,col='green3',border=NA)
  #   }
  #  # lines()
  # }
  
  #shapiro.test(df[,vartest1])
  #uniform.test(hist(df[,vartest1]))
  
  
  
  
  -------------------------------
    source("lib/tri_lib.R")
  #PEAKFILES, BEDFILES, FASTAFILES, 
  #PEAKS, BEDS, FASTAS
  #slice_df
  
  head(PEAKS)
  mypar$treats
  mypar$genes
  mypar$thres
  
  
  # Get BED and FA sequence
  myfa  = FASTAS[FASTAS$chr == params$gene,]
  mybed = BEDS[BEDS$chr == params$gene,]
  myfa$amp.beg = mybed[mybed$feature == "FW_Barcode",]$beg
  myfa$amp.end = mybed[mybed$feature == "RV_Barcode",]$end
  myfa$amp.seq = gsub(paste('^.{',myfa$amp.beg,'}','.(.+)','.{',myfa$amp.end,'}','$'),"\\1",myfa$seq,perl=T)
  
  #myfa$seq = gsub()
  # Get PEAKS
  df = slice_df(PEAKS, params = params,params_not_exact=params_not_exact)
  df$cluster = 1
  df$y = seq(1,dim(df)[1])
  
  df = df[order(df$beg, df$end),]; df$y = seq(1,dim(df)[1])
  df$meanbeg = df$beg
  for (i in seq(1,(max(df$y)),gp$stepsmooth)) {
    myrange = min(i+gp$windowsmooth,max(df$y))
    df$meanbeg[i] = mean(df[seq(i,myrange),]$beg)
  }
  
  df = df[order(df$end, df$beg),]; df$y = seq(1,dim(df)[1])
  df$meanend = df$end
  for (i in seq(1,(max(df$y)),gp$stepsmooth)) {
    myrange = min(i+gp$windowsmooth,max(df$y))
    df$meanend[i] = mean(df[seq(i,myrange),]$end)
  }
  smdf(df)
  
  
  
  df = df[order(df$cluster,df$beg, df$end),]; df$y = seq(1,dim(df)[1])
  df$index = seq(1,dim(df)[1])
  p1.beg = ps(df,'beg',gp=gp,print=F,group='cluster')
  
  df = df[order(df$cluster,df$end, df$beg),]; df$y = seq(1,dim(df)[1])
  p1.end = ps(df,'end',gp=gp,print=F,group='cluster')
  
  
  df = df[order(df$cluster,df$meanbeg, df$meanend),]; df$y = seq(1,dim(df)[1])
  p2.beg = ps(df,by='meanbeg',gp=gp,print=F,group='cluster')
  
  df = df[order(df$cluster,df$meanend, df$meanbeg),]; df$y = seq(1,dim(df)[1])
  p2.end = ps(df,by='meanend',gp=gp,print=F,group='cluster')
  
  grid.arrange(p1.beg,p1.end,p2.beg,p2.end)
  
  
  # Get Cluster
  dfclust = slice_CLUSTS(df,CLUSTS)
  df = get_cluster(df = df, dfclust = dfclust)
  
  df = df[order(df$cluster,df$beg, df$end),]; df$y = seq(1,dim(df)[1])
  p1.beg.c = ps(df,'beg',gp=gp,print=F,group='cluster')
  
  df = df[order(df$cluster,df$end, df$beg),]; df$y = seq(1,dim(df)[1])
  p1.end.c = ps(df,'end',gp=gp,print=F,group='cluster')
  
  
  df = df[order(df$cluster,df$meanbeg, df$meanend),]; df$y = seq(1,dim(df)[1])
  p2.beg.c = ps(df,by='meanbeg',gp=gp,print=F,group='cluster')
  
  df = df[order(df$cluster,df$meanend, df$meanbeg),]; df$y = seq(1,dim(df)[1])
  p2.end.c = ps(df,by='meanend',gp=gp,print=F,group='cluster')
  
  grid.arrange(p1.beg.c,p1.end.c,p2.beg.c,p2.end.c)
  
  ## plot xy
  #p.xy = ggplot.xy(df,dfclust)
  
  # Do tests
  do_distributionTests = function(df = data.frame(), positionTypes = c('mean','orig'), testTypes=c('unif','vhot'), gp=list(), outpdf=F,print=F, debug=F, verbose=F) {
    
    if (defined(gp$minpval) == F) {
      gp$minpval = 0.05
    }
    
    do_distributionTests.get_outpdf = function(positionTypes=positionTypes,testTypes=testTypes,gp=gp,outpdf=outpdf) {
      if (outpdf == F) {
        outpdf = paste(paste(positionTypes,sep='_'),paste(testTypes,sep='_'),'.pdf',sep='')
        if (defined(gp$mytitle)) {
          outpdf = paste(gp$mytitle,'_',outpdf,sep='')
        }
        to_return = outpdf
        return(to_return)
      }
    }
    
    outpdf = do_distributionTests.get_outpdf(print=print,outpdf=outpdf,gp=gp,positionTypes = positionTypes)
    
    if (verbose == T) {print(outpdf)}
    
    return(list())
    
    if (print == T) {pdf(outpdf)}
    
    misc = list(
      p.plot = list(),
      p.name = c(),
      p.type = c(),
      sig.coords = data.frame()
    )
    
    for (positionTypesInd in seq(1,length(positionTypes))) {
      positionType      = positionTypes[positionTypesInd]
      
      positionTypePrint = if (positionType == 'orig') {''} else {positionType}
      positionType1     = paste(positionTypePrint,'beg',sep='')
      positionType2     = paste(positionTypePrint,'end',sep='')
      
      misc[[positionType]] = list(
        positionType1 = positionType1,
        positionType2 = positionType2,
        p.plot = list(),
        p.name = c(),
        sig.coords = data.frame()
      )
      
      df = df[order(df[,cluster],df[,positionType1],df[,positionType2]),]
      
      #ind = 0
      for (testTypesInd in seq(1,length(testTypes))) {
        testType = testTypes[testTypesInd]
        
        misc[[positionType]][[testType]] = list(
          p.plot = list(),
          p.name = c(),
          sig.coords = data.frame(),
          last.i = 1,
          curr.i = 1
        )
        
        clusters = myorder(unique(df$cluster))
        misc[[positionType]][[testType]]$clusters = clusters
        misc[[positionType]][[testType]]$p.type = c(misc$p.type,paste(positionType1,'.',testType,'.pval',sep=''))
        
        p.type = misc[[positionType]][[testType]]$p.type
        
        df[,p.type] = 0
        
        df2 = data.frame()
        
        for (clustersInd in seq(1,length(clusters))) {
          cluster = clusters[clustersInd]
          
          misc[[positionType]][[testType]][[cluster]] = list(
            p.plot = list(),
            p.name = c(),
            sig.coords = data.frame()
          )
          
          dm = df[df$cluster == cluster,]
          
          # Get original min and max y-axis for that cluster
          misc[[positionType]][[testType]][[cluster]]$min.y = min(dm$y)
          misc[[positionType]][[testType]][[cluster]]$max.y = max(dm$y)
          
          dm = dm[order(dm$cluster,dm[,positionType1],dm[,positionType2]),]; dm$y = seq(1,dim(dm)[1]);# = reorder_y(dm)
          
          # Get min and max i for that cluster
          min.i = 1
          max.i = dim(dm)[1]
          
          # Initialize p-value of p.type
          dm[,p.type] = 0
          
          
          # Get min y-axis for that cluster (to be added later)
          
          dm = dm[order(dm$cluster,dm[,positionType1],dm[,positionType2]),]; dm$y = seq(1,dim(dm)[1]);# = reorder_y(dm)
          
          ########################## Loop Start
          # - start at i = 1 until i = myrange; myrange = min(i+window or dim(dm)[1])
          # - check p-value by testType (unif, vhot, hhot, norm, etc)
          #   - if p-value is bigger than pval threshold (gp$minpval), extend and increase pval threshold (+ 1% each iteration) 
          #     until it isn't, then next i jumps to myrange + 1
          #   - else, return and i = i + 1
          
          curr.i = 1
          last.i = 1
          
          while (1) {
            
            if (curr.i > max.i) {
              
              curr.i = 1
              
              # add rectangles to p
              p = p + geom_rect(
                data = misc[[positionType]][[testType]][[cluster]]$sig.coords,
                aes(xmin=xmin,xmax=xmax,ymin=y0,ymax=y1),
                fill='red4',color=NA
              )
              
              misc[[positionType]][[testType]][[cluster]]p.plot = p
              misc[[positionType]][[testType]][[cluster]]p.name = p.type
              
              # if (print == T) {
              #   #p[[clustersInd]][[
              #   myp[[clustersInd]] = local({
              #     i = clustersInd
              #     p1 = p
              #     print(p1)
              #   })
              # }
              
              #mypname = c(mypname,gp$minpval)
              break
            }
            
            #if (misc[[positionType]][[test.type]]$curr.i > max(dm$y)) {
            #  misc[[positionType]][[testType]]$curr.i = 1;
            #  break
            #}
            
            if (verbose == T) {
              if (last.i != curr.i - 1) {
                print(paste(positionType,testType,'cluster=',cluster,'last.i=',last.i,'curr.i=',curr.i))
              }
            }
            
            last.i = curr.i
            
            myrange = min(curr.i + gp$windowdist, max.i)
            seq.myrange  = seq(i,myrange)
            
            
            if (i == 1) {
              
              mydata.ggplot = dm[,positionType1]
              dm.todraw = dm
              dm.todraw$x = dm[,positionType1]
              
              p = ggplot(dm.todraw,aes(x,y)) +
                geom_line(col='grey') +
                coord_cartesian(xlim=c(0,3000),ylim=c(0,max(dm$y))) +
                xlab(paste('R-loop',positionType1,'position')) + ylab('Read Number') +
                ggtitle(paste(testType,'cluster',cluster,'; minpval',gp$minpval))
              last = 'none'
            }
            
            
            misc[[positionType]][[testType]][[cluster]]$seq.myrange.add = 0
            # myrange.temp = min(
            #   misc[[positionType]][[testType]][[cluster]]$curr.i + gp$windowdist + seq.myrange.add,
            #   misc[[positionType]][[testType]][[cluster]]$max.i
            # )
            # seq.myrange.temp = seq(i,myrange.temp)
            # 
            pvalthres = gp$minpval
            
            while (1) {
              #cat('  -',i,myrange.temp,'\n')
              pvalmeanbeg = 0
              misc[[positionType]][[testType]][[cluster]]$myrange = min(
                misc[[positionType]][[testType]][[cluster]]$curr.i + gp$windowdist,
                misc[[positionType]][[testType]][[cluster]]$max.i
              )
              misc[[positionType]][[testType]][[cluster]]$seq.myrange  = seq(misc[[positionType]][[testType]][[cluster]]$curr.i, misc[[positionType]][[testType]][[cluster]]$myrange)
              
              dm.to_test  = dm[seq.myrange.temp,]
              meanbeghist = hist(dm.to_test[,positionType1],plot=F,breaks = max(2,dim(dm.to_test)[1]/10))
              
              if (testType == 'unif') {
                if (length(unique(dm.to_test[,positionType1])) > 1) {
                  pvalmeanbeg = uniform.test(meanbeghist)$p.value
                }
              } else if (testType == 'vhot') {
                if (length(unique(dm.to_test[,positionType1])) < gp$min.len.vhot) {
                  pvalmeanbeg = 1
                }
              }
              dm[seq.myrange.temp,p.type] = as.integer(10000*apply(data.frame(dm[seq.myrange.temp,p.type],rep(pvalmeanbeg,dim(data.frame(dm[seq.myrange.temp,p.type]))[1])),1,max))/10000
              
              if (pvalmeanbeg >= pvalthres) {
                if (seq.myrange.add == 0) {pvalthres = min(1,pvalthres + 1/100)}
                seq.myrange.add = seq.myrange.add + 1
              } else {
                if (seq.myrange.add == 0) {pvalthres = min(1,pvalthres + 1/100)}
                break
              }
              if (myrange.temp == max(dm$y)) {break}
            }
            
            if (seq.myrange.add > 0) {
              ibefore = i
              i = myrange.temp-1
              temp.x0end = min(dm[seq(ibefore,i),positionType1])
              temp.y0end = min(dm[seq(ibefore,i),positionType2])
              temp.x1beg = max(dm[seq(ibefore,i),positionType1])
              temp.y1beg = max(dm[seq(ibefore,i),positionType2])
              
              sig.coords = rbind(sig.coords,data.frame(positionType = positionType1,testType=testType,cluster=cluster,x0end=temp.x0end,y0end=temp.y0end,x1beg=temp.x1beg,y1beg=temp.y1beg,y0=ibefore,y1=i,x=100,y=ibefore,xmin=100,xmax=150))
              #print(paste('ibefore = ',ibefore,'; iafter = ',i))
            }
            
            seq.myrange = seq.myrange.temp
            myrange = myrange.temp
            dmtemp = dm[seq.myrange.temp,]
            dmtemp.unif.myrange = max(1,dim(dmtemp)[1])
            
            unif.col = dist_test.get_color(pvalmeanbeg,gp)
            
            meanbeggp = list()
            meanbeggp$xrectadds = 100
            meanbeggp$xrectmins = 1
            meanbeggp$xrectmaxs = meanbeggp$xrectmins + meanbeggp$xrectadds
            
            
            if (seq.myrange.add > 0) {
              meanbeggp$xlines = dmtemp[,positionType1]
              meanbeggp$ylines = dmtemp$y
              meanbeggp$yrectmins = dmtemp$y[1]
              meanbeggp$yrectmaxs = dmtemp$y[2]
              
              dm.to_draw = data.frame(x=dmtemp[,positionType1]+gp$divby,y=dmtemp$y)
              p = p + geom_line(data=dm.to_draw,aes(x=x,y=y),color=unif.col,lwd=0.5)
            }
            i = i + 1
          }
          
          #
          ########################## Loop End
          
          
          
          df2 = rbind(df2,dm)
          sig.coords$y0 = sig.coords$y0 + min.dm.y
          sig.coords$y1 = sig.coords$y1 + min.dm.y
          misc[[positionType]][[testType]]$sig.coords = rbind(misc[[positionType]][[testType]]$sig.coords,sig.coords)
        }
        df2 = df2[order(df2$cluster,df2[,positionType1],df2[,positionType2]),]
        df2$y = seq(1,dim(df2)[1])
        
        df[df2[,p.type] != 0,p.type] = df2[df2[,p.type] != 0,p.type]# = cbind(df,df2[,dim(df2)[2]]); colnames(df)[dim(df)[2]] = colnames(df2)[dim(df2)[2]]
        sig.coords3 = rbind(sig.coords3,sig.coords2)  
        
      }
    }
    to_return = list(df=df, sig.coords=sig.coords3)
    return(to_return)
  }
  
  positionTypes = c('beg','end')
  positionTypes = c('meanbeg','meanend')
  mylist = do_distributionTests(df,   positionTypes = c('meanbeg','meanend'), gp=gp)
  mydf3 = mylist$df
  mysig.coords3 = mylist$sig.coords
  
  head(df3)
  ggplot.xy(df = df3,dfclust = dfclust)
  df3 = df3[order(df3$cluster,df3$meanbeg, df3$meanend),]; df3$y = seq(1,dim(df3)[1])
  for (testType in testTypes) {
    pval.varname = paste(positionType1,'.',testType,'.pval',sep='')
    varname = paste(positionType1,'.',testType,'.sig.coords',sep='')
    print(varname)
    df3[,varname] = 0
    df3[df3[,pval.varname] > 0.05,varname] = 1 
    #for (i in seq(1,dim(sig.coords3)[1])) {
    #  if (sig.coords3$testType[i] == testType) {
    #    df3[df3$y >= sig.coords3$beg[i] & df3$y <= sig.coords3$end[i],varname] = 1
    #  }
    #}
  }
  positionType1 = 'meanbeg'
  yvar = paste('y.',positionType1,sep='')
  mydf3 = mydf3[order(mydf3$cluster,mydf3$meanbeg,mydf3$meanend),]; mydf3[,yvar] = seq(1,dim(mydf3)[1])
  p2.meanbeg.c.mydf3 = ps(mydf3,by=positionType1,y.var=yvar,gp=gp,print=F,group='cluster') +
    geom_rect(data=mysig.coords3[mysig.coords3$positionType == positionType1 & mysig.coords3$testType == 'unif',],aes(x=0,y=0,xmin=0,ymin=y0,xmax=50,ymax=y1,fill=af(cluster),group=af(cluster))) +
    geom_rect(data=mysig.coords3[mysig.coords3$positionType == positionType1 & mysig.coords3$testType == 'vhot',],aes(x=100,y=100,xmin=100,ymin=y0,xmax=150,ymax=y1,fill=af(cluster),group=af(cluster))) +
    geom_line(data=mydf3,aes(x=meanbeg,y=mydf3[,yvar],group=cluster,alpha=af(meanbeg.vhot.sig.coords)),color='red4',lwd=1) + 
    geom_point(data=mydf3[mydf3$meanend.unif.sig.coords != 0,],aes(x=meanend,y=mydf3[mydf3$meanend.unif.sig.coords != 0,yvar],group=cluster),color='blue2',shape=15,size=0.1) +
    geom_point(data=mydf3[mydf3$meanend.vhot.sig.coords != 0,],aes(x=meanend-1,y=mydf3[mydf3$meanend.vhot.sig.coords != 0,yvar],group=cluster),color='red4',shape=18,size=1) +
    geom_line(data=mydf3,aes(x=meanbeg,y=mydf3[,yvar],group=cluster,alpha=af(meanbeg.unif.sig.coords)),color='black',lwd=0.5) + scale_alpha_manual(values=c('0'=0,'1'=1))
  p2.meanbeg.c.mydf3
  
  positionType1 = 'meanend'
  yvar = paste('y.',positionType1,sep='')
  mydf3 = mydf3[order(mydf3$cluster,mydf3$meanend,mydf3$meanbeg),]; mydf3[,yvar] = seq(1,dim(mydf3)[1])
  p2.meanend.c.mydf3 = ps(mydf3,by=positionType1,y.var=yvar,gp=gp,print=F,group='cluster')
  p2.meanend.c.mydf3
  p2.meanend.c.mydf3 = ps(mydf3,by=positionType1,y.var='y.meanend',gp=gp,print=F,group='cluster') +
    geom_rect(data=mysig.coords3[mysig.coords3$positionType == positionType1 & mysig.coords3$testType == 'unif',],aes(x=0,y=0,xmin=0,ymin=y0,xmax=50,ymax=y1,fill=af(cluster),group=af(cluster))) +
    geom_rect(data=mysig.coords3[mysig.coords3$positionType == positionType1 & mysig.coords3$testType == 'vhot',],aes(x=100,y=100,xmin=100,ymin=y0,xmax=150,ymax=y1,fill=af(cluster),group=af(cluster))) +
    geom_line(data=mydf3,aes(x=meanend,y=mydf3[,yvar],group=cluster,alpha=af(meanend.vhot.sig.coords)),color='red4',lwd=1) +
    geom_point(data=mydf3[mydf3$meanbeg.unif.sig.coords != 0,],aes(x=meanbeg,y=mydf3[mydf3$meanbeg.unif.sig.coords != 0,yvar],group=cluster),color='blue2',shape=15,size=0.1) +
    geom_point(data=mydf3[mydf3$meanbeg.vhot.sig.coords != 0,],aes(x=meanbeg-1,y=mydf3[mydf3$meanbeg.vhot.sig.coords != 0,yvar],group=cluster),color='red4',shape=18,size=1) +
    geom_line(data=mydf3,aes(x=meanend,y=mydf3[,yvar],group=cluster,alpha=af(meanend.unif.sig.coords)),color='black',lwd=0.5) + scale_alpha_manual(values=c('0'=0,'1'=1))
  p2.meanend.c.mydf3
  
  
  #grid.arrange(p2.meanbeg.c.df3,p2.meanend.c.df3,nrow=1,ncol=2)
  
  
  #-------
  plot(df[,positionType1],df$y,type='l',col='grey',xlim = c(0,2500),ylim=c(0,max(df$y)))
  points(df[df[,myvar] > 0.05,positionType1]+10,df[df[,myvar] > 0.05,]$y,pch='.')
  myvar2 = paste(myvar,2,sep='')
  df[,myvar2] = 0
  for (i in 1:max(df$y)) {
    if (df[i,myvar] > 0.05) {
      pos2 = min(dim(df)[1],i + gp$windowdist)
      df[seq(i,pos2),myvar2] = 1
    }
  }
  #df = df[!is.na(df$beg),]
  for (i in 1:max(df$y)) {
    if (df[i,myvar] > 0.05) {
      if (i > 2) {
        df[seq(i-2,i),myvar2] = 0
      }
    }
  }
  points(df[df[,myvar2] != 0,positionType1]+50,df[df[,myvar2] != 0,]$y,pch='.',col='red4')
  