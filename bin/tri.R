#source("lib/tri_lib.R")
#CLUSTFILES, PEAKFILES, BEDFILES, FASTAFILES, 
#CLUSTS, PEAKS, BEDS, FASTAS
#slice_df

head(PEAKS)
mypar$treats
mypar$genes
mypar$thres


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
  divby = 25,
  dist.test_min.length = 5,
  dist.test_min.unique = 2,
  dist.test_max.unique.vhot = 2
)


params = list(
  gene     = 'T7_init_VR_20',
  treat    = '^C$',
  peaktype = 'TOP',
  thres    = 0,
  VR       = 20
)
params_not_exact = c(F,T,F,F,F)

gp$mytitle = get_title(params=params)
gp$mytitle.wrap = wrap_title(gp$mytitle,width = 40)
gp$divby = triclust_get_divby(mytitle=gp$mytitle)$triclust.divby

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

do_distributionTests = function(df = data.frame(), positionTypes = c('mean','orig'), testTypes=c('unif','vhot'), gp=list(), outpdf=F,myprint=F, debug=F, verbose=F) {
  
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

    posTypes = c(positionType1,positionType2)
    for (posTypesInd in seq(1,length(posTypes))) {
      posType1 = posTypes[posTypesInd]
      posType2 = posType1
      if (grepl("beg",posType2)) {
        posType2 = gsub("beg","end",posType2,perl=T)
      } else if (grepl("end",posType2)) {
        posType2 = gsub("end","beg",posType2,perl=T)
      }
      df = df[order(df[,'cluster'],df[,posType1],df[,posType2]),]
      
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
  
        df2 = data.frame()
        
        clusters = myorder(unique(df$cluster))
        misc[[positionType]][[testType]]$clusters = clusters
        misc[[positionType]][[testType]]$p.type = c(misc$p.type,paste(posType1,'.',testType,'.pval',sep=''))
  
        p.type = misc[[positionType]][[testType]]$p.type
        
        df[,p.type] = 0
      
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
          
          ########################## Loop Start
          # - start at i = 1 until i = myrange; myrange = min(i+window or dim(dm)[1])
          # - check p-value by testType (unif, vhot, hhot, norm, etc)
          #   - if p-value is bigger than pval threshold (gp$minpval), extend and increase pval threshold (+ 1% each iteration) 
          #     until it isn't, then next i jumps to myrange + 1
          #   - else, return and i = i + 1
          
          curr.i = 1
          last.i = 1
          
          dm.temp.draw   = dm
          dm.temp.draw$x = dm[,posType1]
          
          p = ggplot(dm.temp.draw,aes(x=x,y=y)) +
            geom_line(col='grey') +
            coord_cartesian(xlim=c(0,3000),ylim=c(0,max(dm$y))) +
            xlab(paste('R-loop',posType1,'position')) + ylab(pasta('Read',posType1,'Number')) +
            ggtitle(paste(testType,'cluster',cluster,'; minpval',gp$minpval))
          
          while (curr.i < max.i) {
            
            # if verbose: if there was a jump bigger than N in curr.i, then print the positions
            if (verbose == T & last.i != curr.i - 1) {
              print(paste(positionType,testType,'cluster=',cluster,'last.i=',last.i,'curr.i=',curr.i))
            }
            
            last.i = curr.i
  
            myrange = list(
              i0  = curr.i,
              i1  = min(curr.i + gp$windowdist, max.i),
              seq = seq(curr.i,min(curr.i + gp$windowdist, max.i)),
              size = length(seq(curr.i,min(curr.i + gp$windowdist, max.i))),
              add = 0,
              p   = 0,
              minpval = gp$minpval,
              minpval2 = gp$minpval2
            )
  
            # iterate from i0 to i1 (see myrange above)
            # if p-value is bigger than pval threshold (gp$myrange$minpval), extend myrange and increase pval threshold (+ 1% each iteration)
            # else return curr.i = i + 1
            while (1) {
              if (verbose == T) {cat('  -',myrange$i0,myrange$i1,myrange$add,myrange$p,'\n')}
              myrange$p = 0
              dm2 = dm[myrange$seq,]
              myrange$hist = hist(dm[myrange$seq,][,posType1],plot=F,breaks = max(2,dim(dm[myrange$seq,][,])[1]/10))
              
              # check each test
              if (testType == 'unif') {
                if (size(dm[myrange$seq,]) >= gp$dist.test_min.length & length(unique(dm[myrange$seq,][,posType1])) >= gp$dist.test_min.unique) {
                  myrange$p = uniform.test(myrange$hist)$p.value
                }
              } else if (testType == 'vhot') {
                if (size(dm[myrange$seq,]) > gp$dist.test_min.length & length(unique(dm[myrange$seq,][,posType1])) < gp$dist.test_max.unique.vhot) {
                  myrange$p = 1
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
            
            if (myrange$add > 0) {
              ibefore = myrange$i0
              curr.i = myrange$i1
              temp.x0end = min(dm[seq(myrange$i0,myrange$i1),posType1])
              temp.y0end = min(dm[seq(myrange$i0,myrange$i1),posType2])
              temp.x1beg = max(dm[seq(myrange$i0,myrange$i1),posType1])
              temp.y1beg = max(dm[seq(myrange$i0,myrange$i1),posType2])
              
              misc[[positionType]][[testType]][[cluster]]$sig.coords = rbind(
                misc[[positionType]][[testType]][[cluster]]$sig.coords,
                data.frame(positionType = posType1,
                           testType=testType,
                           cluster=cluster,
                           x0end=temp.x0end,y0end=temp.y0end,x1beg=temp.x1beg,y1beg=temp.y1beg,
                           y0=myrange$i0,
                           y1=myrange$i1,
                           x=100,y=myrange$i0,xmin=100,xmax=150)
              )
            }
        
            dmtemp.unif.myrange = max(1,dim(dm[myrange$seq,])[1])
            
            myrange$col = dist.test_get.color(myrange$p,gp)
            
         
            if (myrange$add > 0) {
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


positionTypes = c('mean')
mylist = do_distributionTests(df,   positionTypes = positionTypes, gp=gp)
mydf3 = mylist$df
mysig.coords3 = mylist$misc$sig.coords

head(mydf3)
head(mysig.coords3)
#ggplot.xy(df = mydf3,dfclust = dfclust)
mydf3 = mydf3[order(mydf3$cluster,mydf3$meanbeg, mydf3$meanend),]; mydf3$y = seq(1,dim(mydf3)[1])
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
    mydf3 = mydf3[order(mydf3[,'cluster'],mydf3[,posType1],mydf3[,posType2]),]
    mydf3$y = seq(1,dim(mydf3)[1])
 
    
    for (testTypesInd in seq(1,size(testTypes))) {
      testType = testTypes[testTypesInd]
      pval.varname = paste(posType1,'.',testType,'.pval',sep='')
      varname = paste(posType1,'.',testType,'.sig.coords',sep='')
      print(varname)
      mydf3[,varname] = 0
   #   mydf3[mydf3[,pval.varname] > 0.05,varname] = 1
      d1 = print(defined(mysig.coords3[mysig.coords3$positionType == posType1,]$testType))
      d2 = print(defined(mysig.coords3[mysig.coords3$positionType == posType1 & mysig.coords3$testType == testType,]$testType))
      #if (d1 == TRUE & d2 == TRUE) {
      if (defined(mysig.coords3[mysig.coords3$positionType == posType1 & mysig.coords3$testType == testType,]$testType)) {
        mysig.coords3.sub = mysig.coords3[mysig.coords3$positionType == posType1 & mysig.coords3$testType == testType,]
        for (i in seq(1,dim(mysig.coords3.sub)[1])) {
          mydf3[mydf3$y >= mysig.coords3.sub$y0[i] & mydf3$y <= mysig.coords3.sub$y1[i],varname] = 1
          if (testType == 'unif' & posType1 == 'meanbeg') {
            print(paste(posType1,testType,i,mysig.coords3.sub$y0[i],mysig.coords3.sub$y1[i]))
          }
        }
      }
    }
  }
}
positionType1 = 'meanbeg'
yvar = paste('y.',positionType1,sep='')
mydftemp = mydf3[order(mydf3$cluster,mydf3$meanbeg,mydf3$meanend),]; mydftemp[,yvar] = seq(1,dim(mydftemp)[1])
unif1 = mysig.coords3[mysig.coords3$positionType == positionType1 & mysig.coords3$testType == 'unif',]
p2.meanbeg.c.mydf3 = ps(mydftemp,by=positionType1,y.var=yvar,gp=gp,print=F,group='cluster') +
  geom_rect(data=mysig.coords3[mysig.coords3$positionType == positionType1 & mysig.coords3$testType == 'unif',],aes(x=0,y=0,xmin=0,ymin=y0,xmax=50,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
  geom_rect(data=mysig.coords3[mysig.coords3$positionType == positionType1 & mysig.coords3$testType == 'vhot',],aes(x=100,y=100,xmin=100,ymin=y0,xmax=150,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
  geom_line(data=mydftemp,aes(x=meanbeg,y=mydftemp[,yvar],group=cluster,alpha=af(meanbeg.vhot.sig.coords)),color='red4',lwd=1) +
  geom_point(data=mydftemp[mydftemp$meanend.unif.sig.coords != 0,],aes(x=meanend,y=mydftemp[mydftemp$meanend.unif.sig.coords != 0,yvar],group=cluster),color='blue2',shape=15,size=0.1) +
  geom_point(data=mydftemp[mydftemp$meanend.vhot.sig.coords != 0,],aes(x=meanend-1,y=mydftemp[mydftemp$meanend.vhot.sig.coords != 0,yvar],group=cluster),color='red4',shape=18,size=1) +
  geom_line(data=mydftemp,aes(x=meanbeg,y=mydftemp[,yvar],group=cluster,alpha=af(meanbeg.unif.sig.coords)),color='black',lwd=0.5) + scale_alpha_manual(values=c('0'=0,'1'=1))
p2.meanbeg.c.mydf3
p3.meanbeg.c.mydf3 = ggplot(mydftemp,aes(meanbeg,meanend)) +
  geom_point(aes(color=af(cluster)),pch='.') + theme_bw() +
  coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
  geom_rect(data=unif1,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,color=af(cluster)),fill=rgb(0,0,0,0)) +
  geom_rect(data=unif2,aes(x=0,y=0,xmin=y0end,xmax=y1beg,ymin=x0end,ymax=x1beg,fill=af(cluster)),alpha=0.25,color=rgb(0,0,0,0)) +
  geom_rect(data=dfclust,aes(x=x0end,y=y0end,xmin=x0end,ymin=y0end,xmax=x1beg,ymax=y1beg,color=as.factor(cluster)),fill=rgb(0,0,0,0),lwd=1) +
#  geom_rect(data=unif2,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,color=af(cluster)),fill=rgb(0,0,0,0)) +
  geom_text(data=dfclust,aes(x=(x1beg+x0end)/2,y=(y1beg+y0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5) +
  annotate(geom='segment',x=0,y=0,xend=2500,yend=2500) + theme(legend.position='none')


positionType2 = 'meanend'
yvar2 = paste('y.',positionType2,sep='')
mydftemp2 = mydf3[order(mydf3$cluster,mydf3$meanend,mydf3$meanbeg),]; mydftemp2[,yvar2] = seq(1,dim(mydftemp2)[1])
unif2 = mysig.coords3[mysig.coords3$positionType == positionType2 & mysig.coords3$testType == 'unif',]
p2.meanend.c.mydf3 = ps(mydftemp2,by=positionType2,y.var='y.meanend',gp=gp,print=F,group='cluster') +
  geom_rect(data=mysig.coords3[mysig.coords3$positionType == positionType2 & mysig.coords3$testType == 'unif',],aes(x=0,y=0,xmin=0,ymin=y0,xmax=50,ymax=y1,fill=af(cluster),group=af(cluster)),color='black',lwd=0.5) +
  geom_rect(data=mysig.coords3[mysig.coords3$positionType == positionType2 & mysig.coords3$testType == 'vhot',],aes(x=100,y=100,xmin=100,ymin=y0,xmax=150,ymax=y1,fill=af(cluster),group=af(cluster)),color='black',lwd=0.5) +
  geom_line(data=mydftemp2,aes(x=meanend,y=mydftemp2[,yvar2],group=cluster,alpha=af(meanend.vhot.sig.coords)),color='red4',lwd=1) +
  geom_point(data=mydftemp2[mydftemp2$meanbeg.unif.sig.coords != 0,],aes(x=meanbeg,y=mydftemp2[mydftemp2$meanbeg.unif.sig.coords != 0,yvar2],group=cluster),color='blue2',shape=15,size=0.1) +
  geom_point(data=mydftemp2[mydftemp2$meanbeg.vhot.sig.coords != 0,],aes(x=meanbeg-1,y=mydftemp2[mydftemp2$meanbeg.vhot.sig.coords != 0,yvar2],group=cluster),color='red4',shape=18,size=1) +
  geom_line(data=mydftemp2,aes(x=meanend,y=mydftemp2[,yvar2],group=cluster,alpha=af(meanend.unif.sig.coords)),color='black',lwd=0.5) + scale_alpha_manual(values=c('0'=0,'1'=1))
p2.meanend.c.mydf3

p3.meanend.c.mydf3 = ggplot(mydftemp2,aes(meanend,meanbeg)) +
  geom_point(aes(color=af(cluster)),pch='.') + theme_bw() +
  coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
#  geom_rect(data=unif1,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,color=af(cluster)),fill=rgb(0,0,0,0)) +
  geom_rect(data=unif2,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,fill=af(cluster)),alpha=0.25,color=rgb(0,0,0,0)) +
#  geom_rect(data=unif2,aes(x=0,y=0,xmin=y0end,xmax=y1beg,ymin=x0end,ymax=x1beg,color=af(cluster)),fill=rgb(0,0,0,0)) +
  geom_rect(data=unif1,aes(x=0,y=0,xmin=y0end,xmax=y1beg,ymin=x0end,ymax=x1beg,color=af(cluster)),fill=rgb(0,0,0,0)) +
#  geom_rect(data=unif2,aes(x=0,y=0,xmin=y0end,xmax=y1beg,ymin=x0end,ymax=x1beg,color=af(cluster)),fill=rgb(0,0,0,0)) +
  #geom_rect(data=dfclust,aes(x=y0end,y=x0end,xmin=y0end,ymin=x0end,xmax=y1beg,ymax=x1beg,color=as.factor(cluster)),fill=rgb(0,0,0,0),lwd=1) +
  #geom_text(data=dfclust,aes(x=(y1beg+y0end)/2,y=(x1beg+x0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5) +
  annotate(geom='segment',x=0,y=0,xend=2500,yend=2500) + theme(legend.position='none')


pdf(paste(gp$mytitle,'.pdf',sep=''),height=20,width=20)
grid.arrange(p3.meanbeg.c.mydf3,p3.meanend.c.mydf3,p2.meanbeg.c.mydf3,p2.meanend.c.mydf3,nrow=2,ncol=2)
dev.off()