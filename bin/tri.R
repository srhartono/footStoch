source("lib/tri_lib.R")
library(pdfCluster)
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
  dist.test_min.length = 5,
  dist.test_min.unique = 2,
  dist.test_max.unique.hots = 2,
  divby=25
)

myVR = 20

for (myVR in seq(1,31)) {
  print(myVR)
  myparams = list(
    gene     = 'T7_init',#_VR_20',
    treat    = 'C',
    peaktype = 'TOP',
    thres    = 0,
    VR       = myVR
    
  )
  
  myparams$gene = paste('T7_init_VR_',myVR,sep='')
  myparams_regex = list(
    gene = F,
    treat = F,
    peaktype = F,
    thres = F,
    VR = F
  )
  
  gp$mytitle = get_title(myparams=myparams)
  gp$mytitle.wrap = wrap_title(gp$mytitle,width = 40)
  gp$divby = triclust_get_divby(mytitle=gp$mytitle)$triclust.divby
  gp$promoter = mybed[mybed$feature == 'T7_Promoter',]
  # Get BED and FA sequence
  myfa  = FASTAS[FASTAS$chr == myparams$gene,]
  mybed = BEDS[BEDS$chr == myparams$gene,]
  myfa$amp.beg = mybed[mybed$feature == "FW_Primer",]$beg-10
  myfa$amp.end = mybed[mybed$feature == "RV_Primer",]$end+10
  #myfa$amp.seq = gsub(paste('^.{',myfa$amp.beg,'}','.(.+)','.{',myfa$amp.end,'}','$'),"\\1",myfa$seq,perl=T)
  myfa$amp.seq = myfa$seq
  # myfapos = get_fapos(myfa)
  # beg2 = data.frame(beg=myfapos$index,Cpos = myfapos$Cindex,begC = myfapos$Cindex,begdiff = myfapos$diff)
  # end2 = data.frame(end=myfapos$index+1,Cpos = myfapos$Cindex - 1,endC = myfapos$Cindex-1)
  # end2 = merge(end2,beg2,by='Cpos')
  # end2 = rbind(end2,data.frame(Cpos=max(end2$Cpos)+1,end=max(end2$end)+1,endC=max(end2$Cpos)+1,beg=max(beg2$beg),begC=0,begdiff=1))
  # end2 = subset(end2,select=c('end','Cpos','endC','beg','begdiff'));
  # beg2$begpos = beg2$beg
  # end2$endpos = end2$beg
  # end2 = subset(end2,select=-beg)
  # end2 = subset(end2,select=-Cpos)
  # beg2 = subset(beg2,select=-Cpos)
  # colnames(end2)[3] = 'enddiff'
  # head(beg2);head(end2);
  # tail(beg2);tail(end2);
  doclust = function(df) {
    begs = aggregate(df$beg,by=list(df$begcluster),min);colnames(begs) = c('begcluster','begmin')
    ends = aggregate(df$end,by=list(df$begcluster),min);colnames(ends) = c('begcluster','endmin')
    dfclustbeg = merge(begs,ends,by='begcluster')
    dfclustbeg = dfclustbeg[order(dfclustbeg$begmin,dfclustbeg$endmin),]
    dfclustbeg$begcluster2 = seq(0,size(dfclustbeg)-1)
    #dfclustbeg
    df = merge(df,subset(dfclustbeg,select=c('begcluster','begcluster2')),by='begcluster')
    df$begcluster = df$begcluster2
    #df = df[order(df$begcluster,df$beg, df$end),]; df$y = seq(1,dim(df)[1])
    #ps(df,'beg',gp=gp,print=T,group='begcluster')
    

    begs = aggregate(df$beg,by=list(df$endcluster),min);colnames(begs) = c('endcluster','begmin')
    ends = aggregate(df$end,by=list(df$endcluster),min);colnames(ends) = c('endcluster','endmin')
    dfclustend = merge(begs,ends,by='endcluster')
    ends2 = aggregate(df$end,by=list(df$endcluster),max);colnames(ends2) = c('endcluster','endmax')
    dfclustend = merge(dfclustend,ends2,by='endcluster')
    dfclustend = dfclustend[order(dfclustend$endmin,dfclustend$begmin),]
    dfclustend$endcluster2 = seq(0,size(dfclustend)-1)
    #dfclustend
    df = merge(df,subset(dfclustend,select=c('endcluster','endcluster2')),by='endcluster')
    df$endcluster = df$endcluster2
    #df = df[order(df$endcluster,df$end, df$beg),]; df$y = seq(1,dim(df)[1])
    #ps(df,'end',gp=gp,print=T,group='endcluster')
    df$allcluster = paste(df$begcluster,df$endcluster)
    allcluster = data.frame(allcluster = unique(df$allcluster),allcluster.name=1)
    allcluster = allcluster[order(allcluster$allcluster),]
    allcluster$allcluster.name = seq(1,size(allcluster))
    df = merge(df,allcluster,by='allcluster')
    
    begs0 = aggregate(df$beg,by=list(df$allcluster,df$allcluster.name),min);colnames(begs0) = c('allcluster','allcluster.name','x0end')
    begs1 = aggregate(df$beg,by=list(df$allcluster,df$allcluster.name),max);colnames(begs1) = c('allcluster','allcluster.name','x1beg')
    ends0 = aggregate(df$end,by=list(df$allcluster,df$allcluster.name),min);colnames(ends0) = c('allcluster','allcluster.name','y0end')
    ends1 = aggregate(df$end,by=list(df$allcluster,df$allcluster.name),max);colnames(ends1) = c('allcluster','allcluster.name','y1beg')
    allcluster = merge(begs0,ends0,by=c('allcluster','allcluster.name'))
    allcluster = merge(allcluster,begs1,by=c('allcluster','allcluster.name'))
    allcluster = merge(allcluster,ends1,by=c('allcluster','allcluster.name'))
    #allcluster = merge(allcluster,beg1,by='allcluster')
    allcluster$begcluster = gsub('^([0-9]+)\\s+([0-9])+$','\\1',allcluster$allcluster,perl=T)
    allcluster$endcluster = gsub('^([0-9]+)\\s+([0-9]+)$','\\2',allcluster$allcluster,perl=T)
    res = list(df=df,allcluster=allcluster)
    return(res)
  }
  
  #[order(df[df$beg >= 400 & df$beg <= 500,'beg']),] #,n=20)
  #myfa$seq = gsub()
  # Get PEAKS
  df = slice_df(PEAKS, myparams = myparams,myparams_regex=myparams_regex)
  df$index = seq(1,dim(df)[1])
  df$meanbeg = df$beg
  df$meanend = df$end
  df$begI = df$beg #ai(df$beg / gp$divby)
  df$endI = df$end #ai(df$end / gp$divby)
  df$cluster = 0
  df$y = seq(1,dim(df)[1])
  df$begcluster = 0
  df$endcluster =0 
  
  #mycluster
  myclust.annot = function(df,dfclust) {
    if (defined(dfclust[dfclust$type == 'beg',])) {
      dfclustbeg = dfclust[dfclust$type == 'beg',]
      dfclustbeg$beg = dfclustbeg$pos
      dfclustbeg$cluster = seq(1,size(dfclustbeg))
      for (i in 1:size(dfclustbeg)) {
        if (defined(df[ai(df$beg/gp$divby) >= dfclustbeg$beg[i],])) {
          df[ai(df$beg/gp$divby) >= dfclustbeg$beg[i],]$begcluster = dfclustbeg$cluster[i]
        }
      }
      # df = df[order(df$begcluster,df$beg, df$end),]; df$y = seq(1,dim(df)[1])
      # ps(df,'beg',gp=gp,print=T,group='begcluster')
    }
    if (defined(dfclust[dfclust$type == 'end',])) {
      dfclustend = dfclust[dfclust$type == 'end',]
      dfclustend$end = dfclustend$pos+1
      dfclustend$cluster = seq(1,size(dfclustend))
      for (i in 1:size(dfclustend)) {
        if (defined(df[ai(df$end/gp$divby) >= dfclustend$end[i],])) {
          df[ai(df$end/gp$divby) >= dfclustend$end[i],]$endcluster = dfclustend$cluster[i]
        }
      }
    }
    return(df)
  }
  # df = df[order(df$endcluster,df$end, df$end),]; df$y = seq(1,dim(df)[1])
  # ps(df,'end',gp=gp,print=T,group='endcluster')
  mythres = ai(sqrt(size(df)/9230)*100)

  myclust = get_dm3(dm2 = df,divby = gp$divby,threshold = mythres)
  mydf= myclust.annot(df,myclust)
  mres = doclust(mydf)
  dfm = mres$df
  dfmclust = mres$allcluster
  
  
  #pdf cluster
  dfp = df
  dfp$begcluster = groups(pdfCluster(data.frame(beg=dfp$beg)))#,bwtype='adaptive'))
  
  dfp$endcluster = groups(pdfCluster(data.frame(end=dfp$end),n.grid=15))#,bwtype='adaptive'))
  pres = doclust(dfp)
  dfp = pres$df
  dfpclust = pres$allcluster

  #  dfp$begcluster = 
  #kmeans cluster
  dfk = df
  dfk$begcluster = kmeans(dfk$beg,7)$cluster
  dfk$endcluster = kmeans(dfk$end,5)$cluster
  kres = doclust(dfk)
  dfk = kres$df
  dfkclust = kres$allcluster
  
  get_colorgroup = function(df,dfclust) {
    tempz = df
    tempzclust = dfclust
    tempz$cluster = tempz$allcluster.name
    tempz = get_y(tempz[order(tempz$cluster,tempz$beg,tempz$end),])
    
    colzbin = c('Greens','Reds','YlOrBr','Purples','Blues','RdPu')
    colzs = rep(colzbin,10)
    colzs2 = c(rgb(0.25,0.25,0.25))
    for (i in 1:size(colzs)) {
      if (i == 1) {
        plot(x=c(1,100),y=c(1,100),col=colzs2[0],type='l',lwd=5,xlim=c(1,100),ylim=c(1,100));
      }
      colz = colzs[i]
      colzdf = rep(brewer.pal(9,colz)[3:9],3) 
      if (i == 1) {
        colzs2 = c(rgb(0.25,0.25,0.25),colzdf) #as.data.frame(t(colzdf))
      }
      else {
        colzs2 = c(colzs2,colzdf) #rbind(colzs2,t(colzdf))
      }
      lines(x=c(1+i,i+1+100),y=c(1,100),col=colzs2[i],type='l',lwd=5)
    }
    #colzs2 = c(rgb(0.25,0.25,0.25),colzs2)#,rgb(0.5,0.5,0.5))
    
    ind=1
    for (i in seq(2,size(colzs2),21)) {#(22+20)) {
      if (ind == 1) {
        plot(x=c(1,100),y=c(1,100),col=colzs2[1],type='l',lwd=5,xlim=c(1,100),ylim=c(1,100));
      }
      ind = ind + 1
      lines(x=c(1+i,i+1+100),y=c(1,100),col=colzs2[i],type='l',lwd=5)
    }
      lines(x=c(1+24,24+1+100),y=c(1,100),col=colzs2[24],type='l',lwd=5)
    
    colzs2 = c(rgb(0.25,0.25,0.25),colzs2)#,rgb(0.5,0.5,0.5))
    mybreak = c(0,seq(1,size(colzs2)))#,-2)
    tempzclust$begc = ai((tempzclust$x0end-gp$promoter$beg) / 100)
    tempzclust$endc = ai((tempzclust$y1beg-gp$promoter$beg) / 100)
    if (defined(tempzclust[tempzclust$endc < 0,])) {tempzclust[tempzclust$endc < 0,]$endc = 0}
    if (defined(tempzclust[tempzclust$endc > 21,])) {tempzclust[tempzclust$endc > 21,]$endc = 21}
    tempzclust$endc2 = ai(tempzclust$endc) + 1
    tempzclust$colorgroup = tempzclust$begc * 21 + tempzclust$endc2
    tempzclust[tempzclust$colorgroup <= 0,]$colorgroup = 0
    tempz2 = merge(tempz,subset(tempzclust,select=c('allcluster.name','colorgroup')),by='allcluster.name')
    return(list(df=tempz2,dfclust=tempzclust))
  }
  
  dfmlist2 = get_colorgroup(dfm,dfmclust)
  dfklist2 = get_colorgroup(dfk,dfkclust)
  dfplist2 = get_colorgroup(dfp,dfpclust)
  #dfclust = aggregate(dfclust$)
  graph.clust(dfmlist2$df,dfmlist2$dfclust,paste('results/clustby/',gp$mytitle,'_mythres.pdf',sep=''),gp,'colorgroup')
  graph.clust(dfklist2$df,dfklist2$dfclust,paste('results/clustby/',gp$mytitle,'_kmeans5.pdf',sep=''),gp,'colorgroup')
  graph.clust(dfplist2$df,dfplist2$dfclust,paste('results/clustby/',gp$mytitle,'_pdfclust.pdf',sep=''),gp,'colorgroup')
  next
  #segments(0,0,2500,2500)
  
  # get meanbeg and meanend
  df = df[order(df$beg, df$end),]; df$y = seq(1,dim(df)[1])
  meanbeg = get_meanpos(df$beg,gp$windowsmooth,gp$stepsmooth)
  df = df[order(df$end, df$beg),]; df$y = seq(1,dim(df)[1])
  meanend = get_meanpos(df$end,gp$windowsmooth,gp$stepsmooth)
  
  smdf(df)
  
  # print non cluster graph
  df = df[order(df$cluster,df$beg, df$end),]; df$y = seq(1,dim(df)[1])
  p1.beg = ps(df,'beg',gp=gp,print=F,group='cluster')
  
  df = df[order(df$cluster,df$end, df$beg),]; df$y = seq(1,dim(df)[1])
  p1.end = ps(df,'end',gp=gp,print=F,group='cluster')
  
  df = df[order(df$cluster,df$meanbeg, df$meanend),]; df$y = seq(1,dim(df)[1])
  p2.beg = ps(df,by='meanbeg',gp=gp,print=F,group='cluster')
  
  df = df[order(df$cluster,df$meanend, df$meanbeg),]; df$y = seq(1,dim(df)[1])
  p2.end = ps(df,by='meanend',gp=gp,print=F,group='cluster')
  
  p0 = grid.arrange(p1.beg,p1.end,p2.beg,p2.end)
  
  test1 = df
  test1$count = 1
  test2 = aggregate(test1$count,by=list(test1$begI,test1$endI),sum)
  colnames(test2) = c('begI','endI','sum')
  test2 = test2[order(test2$sum,test2$begI,test2$endI),]
  test2$perc = ai(test2$sum/sum(test2$sum)*1000)/10
  
  test2beg = aggregate(test1$count,by=list(test1$begI),sum)
  test2end = aggregate(test1$count,by=list(test1$endI),sum)
  colnames(test2beg) = c('beg','sum')
  colnames(test2end) = c('end','sum')
  allpos = c(test2beg$sum,test2end$sum)
  allperc = ai(allpos/sum(allpos)*1000)/10
  test2.threshold = mean(allperc)+sd(allperc)
  
  test3beg = get_edge_wrapper(test2beg,'beg',test2.threshold)
  test3end = get_edge_wrapper(test2end,'end',test2.threshold)
  plot(test2end[,'end'],test2end$perc,type='l',col='grey',xlab='Rloop End Positions',ylab='Density (%)')
  #points(test2end[,'end'],test2end$perc.smooth,pch=15,cex=0.5,col='grey')
  lines(test2end[,'end'],test2end$perc.smooth,col='red2')
  points(test3end[,'end'],test3end$perc.smooth,col='blue2',pch=15,cex=0.5)
  rect(endz4a,rep(-1,length(endz4a)),endz4b,rep(0,length(endz4b)))
  points(test4end[,'end'],test3beg$perc.smooth,col='blue2',pch=15,cex=0.5)
  get_edge_wrapper = function(test2beg,type,threshold) {
    mysum = sum(test1$count)
    colnames(test2end) = c('end','sum')
    test2end$perc = ai(test2end$sum/sum(test2end$sum)*1000)/10
    test2end$perc.smooth = smooth.spline(test2end$perc,spar=0.4)$y
    test2end = test2end[order(test2end[,type]),]
    test2end$y = seq(1,size(test2end))
    
    test3beg = test2beg[test2beg$perc.smooth > test2.threshold,]# + 1*sd(test2beg$perc),]
    get_edge = function(test2beg,test2.threshold) {
      plot(test2end[,type],test2end$perc,type='l',col='grey')
      #points(test2end[,type],test2end$perc.smooth,pch=15,cex=0.5,col='grey')
      lines(test2end[,type],test2end$perc.smooth,col='red2')
      points(test3end[,type],test3end$perc.smooth,col='blue2',pch=15,cex=0.5)
      #points(test3beg2[,type],rep(4,length(test3beg2$perc.smooth)),col='red4',pch=15,cex=0.5)
      
      begz = findedge(x = test3beg$y); begz = begz[begz != 0]
      begz = test3beg$y[begz]
      #begz0 = test2beg[test2beg$y %in% begz,]
      if (!min(test3beg$y) %in% begz) {begz = c(min(test3beg$y),begz)}
      if (!max(test3beg$y) %in% begz) {begz = c(begz,max(test3beg$y))}
      begz
      begz1 = begz[seq(1,length(begz)-1)] + 1
      begz2 = begz[seq(2,length(begz))]
      begz1
      begz2
      begz1 = test2beg[test2beg$y %in% begz1,][,type]
      begz2 = test2beg[test2beg$y %in% begz2,][,type]
      begz4a = begz1
      begz4b = begz2
      res = data.frame(x=begz4a,y=begz4b)
      rect(begz4a,rep(-1,length(begz4a)),begz4b,rep(0,length(begz4b)))
      return(res)
    }
    test4beg = get_edge(test2beg,test2.threshold)
    
    return(test4beg)
  }
  plot(test2$beg,test2$end,type='p',pch='.',col='grey')
  segments(0,0,2000,2000)
  points(test3beg$x,)
  
  
  
  #
  #1-34
  # 34-59
  # 
  # begz2 = test3beg$begI[begz]
  # if (!min(test3beg$begI) %in% begz2) {
  #   begz2 = c(min(test3beg$begI),begz2)
  # }
  # if (!max(test3beg$begI) %in% begz2) {
  #   begz2 = c(begz2,max(test3beg$begI))
  # }
  # begz2
  #begz2 = test3beg[begz,]$begI
  # begz2 = findedge(begz)
  # begz3 = test3beg$begI[begz[begz2]]
  # begz4 = test3beg[test3beg$begI %in%begz3,]
  
  # begz4a = begz2[c(min(test3beg$begI),seq(1,length(begz)-1))]
  # begz4b = begz3[c(min(test3beg$begI),seq(2,length(begz)))]
  #begz4b = begz3[2:length(begz)]
  #plot(begz3)
  # #rect(begz3,)
  # #ends3 = begz3[begz2 ] + 
  # #y18 <- c(1:3, 5, 4, 7:3, 2*(2:5), rep(10, 4))
  # #xx  <- seq(1, length(y18), length.out = 201)
  # y18 = begz0
  # xx = seq(1,length(y18))
  # #xx = xx
  # (s2   <- smooth.spline(y18)) # GCV
  # (s02  <- smooth.spline(y18, spar = 0.2))
  # (s02. <- smooth.spline(y18, spar = 0.2, cv = NA))
  # plot(xx,y18, main = deparse(s2$call), col.main = 2)
  # lines(s2, col = "gray");
  # s3 = predict(s2, xx2)
  # lines(xx2,s3$y, type='l',col = 2)
  # lines(predict(s02, xx2), col = 3);
  # mtext(deparse(s02$call), col = 3)
  
  
  #points(test3beg$begI,test3beg$perc,col='red4',pch=15)
  points(test3beg$begI,test3beg$perc.smooth,col='red4',pch=15)
  
  test2end$perc = ai(test2end$sum/sum(test2end$sum)*1000)/10
  test3end = test2end[test2end$perc > 3*sd(test2end$perc),]
  
  #test3beg = test2beg[test2beg$sum > 39,]
  #test3end = test2end[test2end$sum > 39,]
  test3 = test2[test2$begI %in% test3beg$begI | test2$endI %in% test3end$endI,]
  #test4 = test3[test3$perc > 2*sd(test3$perc),]
  plot(test1$begI,test1$endI,type='p',col='grey',cex=0.3,pch=15,xlim=c(0,100),ylim=c(0,100))#
  segments(0,0,100,100)
  #points(test4$begI,test4$endI,cex=0.5,pch=15,col='red4')
  points(test3beg$begI,rep(0,size(test3beg)),cex=0.3,pch=15,col='red4')
  points(rep(0,size(test3end)),test3end$endI,cex=0.3,pch=15,col='blue4')
  points(test3$begI,test3$endI,col='black',cex=0.3,pch=15)
  # Get Cluster
  dfclust = slice_CLUSTS(df,CLUSTS)
  #saveRDS(dfclust,file='dfclustVR20.RDS')
  dfclust = readRDS('dfclustVR20.RDS')
  #dfclustVR20 = dfclust
  # Get Cluster from VR20
  # 
  # dfclust[dfclust$cluster == 1,]$x0end = dfclust[dfclust$cluster == 1,]$x0end - gp$divby
  # dfclust[dfclust$cluster == 2,]$x0end = dfclust[dfclust$cluster == 2,]$x0end - gp$divby
  # dfclust[dfclust$cluster == 8,]$x0end = dfclust[dfclust$cluster == 8,]$x0end - gp$divby
  # dfclust[dfclust$cluster == 0,]$x1beg = dfclust[dfclust$cluster == 1,]$x0end - 1
  # dfclust[dfclust$cluster == 2,]$y1beg = 1000
  # dfclust[dfclust$cluster == 4,]$y1beg = 1000
  # dfclust[dfclust$cluster == 6,]$y1beg = dfclust[dfclust$cluster == 6,]$y1beg + gp$divby
  # dfclust[dfclust$cluster == 7,]$y1beg = dfclust[dfclust$cluster == 7,]$y1beg + gp$divby
  # dfclust[dfclust$cluster == 8,]$y0end = dfclust[dfclust$cluster == 8,]$y0end + gp$divby
  # 
  # #dfclust = rbind(dfclust,
  # dfclust3 = data.frame(
  #   cluster=3,
  #   x0end=dfclust[dfclust$cluster == 2,]$x0end,
  #   y0end=dfclust[dfclust$cluster == 2,]$y1beg,
  #   x0beg=dfclust[dfclust$cluster == 2,]$x0end,
  #   y0beg=dfclust[dfclust$cluster == 6,]$y1beg,
  #   x1end=dfclust[dfclust$cluster == 2,]$x0end,
  #   y1end=dfclust[dfclust$cluster == 6,]$y1beg,
  #   x1beg=dfclust[dfclust$cluster == 4,]$x1beg,
  #   y1beg=dfclust[dfclust$cluster == 6,]$y1beg
  # )
  # dfclust[dfclust$cluster == 3,] = dfclust3
  # dfclust[dfclust$cluster == 2,]$x1beg = dfclust[dfclust$cluster == 4,]$x1beg
  # dfclust[dfclust$cluster == 2,]$y1beg = dfclust[dfclust$cluster == 4,]$y1beg
  # 
  # dfclust = dfclust[dfclust$cluster != 4,]
  # 
  # dfclust[dfclust$cluster == 1,]$y1beg = dfclust[dfclust$cluster == 1,]$y1beg + gp$divby
  # dfclust[dfclust$cluster == 1,]$x1beg = dfclust[dfclust$cluster == 1,]$x1beg + gp$divby
  # dfclust[dfclust$cluster == 2,]$y0end = dfclust[dfclust$cluster == 2,]$y0end + gp$divby
  
  # dfclust = dfclustVR20
  #dfclustVR20 = dfclust
  #head(dfclust)
  df = get_cluster(df = df, dfclust = dfclust)
  # 
  # dfclust2 = aggregate(df$beg,by=list(df$cluster),min)
  # colnames(dfclust2) = c('cluster','x0end')
  # dfclust2$x0end = aggregate(df$beg,by=list(df$cluster),min)$x
  # dfclust2$y0end = aggregate(df$end,by=list(df$cluster),min)$x
  # dfclust2$x1beg = aggregate(df$beg,by=list(df$cluster),max)$x
  # dfclust2$y1beg = aggregate(df$end,by=list(df$cluster),max)$x
  # dfclust2$x0beg = dfclust2$x0end
  # dfclust2$y0beg = dfclust2$y1beg
  # dfclust2$x1end = dfclust2$x0end
  # dfclust2$y1end = dfclust2$y1beg
  # 
  # dfclust = dfclust2
  # 
  # #df = df[order(df$cluster,df$beg, df$end),]; df$y = seq(1,dim(df)[1])
  # #test7 = df
  # #test7$beg = test7$begC
  # #test7$beg = test7$begC
  # #p1.beg.c = ps(df,'begC',gp=gp,print=F,group='cluster')
  
  df = df[order(df$cluster,df$beg, df$end),]; df$y = seq(1,dim(df)[1])
  p1.beg.c = ps(df,'beg',gp=gp,print=F,group='cluster')
  
  df = df[order(df$cluster,df$end, df$beg),]; df$y = seq(1,dim(df)[1])
  p1.end.c = ps(df,'end',gp=gp,print=F,group='cluster')
  
  
  df = df[order(df$cluster,df$meanbeg, df$meanend),]; df$y = seq(1,dim(df)[1])
  p2.beg.c = ps(df,by='meanbeg',gp=gp,print=F,group='cluster')
  
  df = df[order(df$cluster,df$meanend, df$meanbeg),]; df$y = seq(1,dim(df)[1])
  p2.end.c = ps(df,by='meanend',gp=gp,print=F,group='cluster')
  
  #grid.arrange(p1.beg.c,p1.end.c,p2.beg.c,p2.end.c)
  df = df[df$cluster != -1,]
  ## plot xy
  #dfclust = data.frame()
  #df$cluster = 1
  #p.xy = ggplot.xy(df,dfclust)
  #plot(df$beg,df$end,pch='.')
  # Do tests
  
  positionTypes = c('mean')#orig') #orig, mean
  testTypes = c('unif','hots','norm')
  mylist = do_distributionTests(df, positionTypes = positionTypes, testTypes = testTypes, gp = gp)
  mydf = mylist$df
  mysc = mylist$misc$sig.coords
  
  mydf = annot.pvar(mydf,mysc)
  mysc = re.sc(mydf,mysc,positionTypes,testTypes)
  mydf = annot.pvar(mydf,mysc)
  head(mydf)
  head(mysc)
  #ggplot.xy(df = mydf,dfclust = dfclust)
  
  
  posType1 = 'meanbeg'
  posType2 = 'meanend'
  positionType1 = 'meanbeg'
  positionType2 = 'meanend'
  
  yvar.meanbeg = paste('y.',positionType1,sep='')
  mydftemp.meanbeg = mydf[order(mydf$cluster,mydf[,positionType1],mydf[,positionType2]),];
  mydftemp.meanbeg[,yvar.meanbeg] = seq(1,dim(mydftemp.meanbeg)[1])
  mysc.meanbeg = mysc[mysc$positionType == positionType1,]
  unif.meanbeg = mysc[mysc$positionType == positionType1 & mysc$testType == 'unif',]
  yvar.meanend = paste('y.',positionType2,sep='')
  mydftemp.meanend = mydf[order(mydf$cluster,mydf[,positionType2],mydf[,positionType1]),];
  mydftemp.meanend[,yvar.meanend] = seq(1,dim(mydftemp.meanend)[1])
  mysc.meanend = mysc[mysc$positionType == positionType2,]
  unif.meanend = mysc[mysc$positionType == positionType2 & mysc$testType == 'unif',]
  
  mycolors = c(brewer.pal(9,"Set1"),brewer.pal(9,"Set3"))
  mybreaks = seq(0,length(mycolors)-1)
  
  
  p0.meanbeg.c.mydf = ggplot(mydftemp.meanbeg,aes(beg,end)) +
    geom_point(aes(color=af(cluster)),pch='.') + theme_bw() +
    #  geom_point(color='black',pch='.') + theme_bw() +
    coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +  
    annotate(geom='segment',x=0,y=0,xend=3000,yend=3000) + theme(legend.position='none')
  
  p1.meanbeg.c.mydf = ggplot(mydftemp.meanbeg,aes(beg,end)) +
    geom_point(aes(color=af(cluster)),pch='.') + theme_bw() +
    coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
    #  geom_rect(data=unif.meanbeg,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,color=af(cluster)),fill=rgb(0,0,0,0)) +
    #  geom_rect(data=unif.meanend,aes(x=0,y=0,xmin=y0end,xmax=y1beg,ymin=x0end,ymax=x1beg,fill=af(cluster)),alpha=0.25,color=rgb(0,0,0,0)) +
    geom_text(data=dfclust,aes(x=(x1beg+x0end)/2,y=(y1beg+y0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5) +
    geom_rect(data=dfclust,aes(x=x0end,y=y0end,xmin=x0end,ymin=y0end,xmax=x1beg,ymax=y1beg,color=as.factor(cluster)),fill=rgb(0,0,0,0),lwd=1) +
    annotate(geom='segment',x=0,y=0,xend=3000,yend=3000) + theme(legend.position='none')
  
  dfclust.meanbeg = aggregate(mydftemp.meanbeg$y.meanbeg,by=list(mydftemp.meanbeg$cluster),min); colnames(dfclust.meanbeg) = c("cluster","y.meanbeg.min")
  dfclust.meanbeg$y.meanbeg.max = aggregate(mydftemp.meanbeg$y.meanbeg,by=list(mydftemp.meanbeg$cluster),max)$x
  dfclust.meanbeg = merge(dfclust.meanbeg,dfclust,by="cluster")
  p2.meanbeg.c.mydf = ps(mydftemp.meanbeg,by=posType1,y.var=yvar.meanbeg,gp=gp,print=F,group='cluster') +
    geom_segment(data=dfclust.meanbeg,aes(x=x0end,xend=x0end,y=y.meanbeg.min,yend=y.meanbeg.max,group=af(cluster)),color='black',lty=1) +
    #  geom_segment(data=dfclust.meanbeg,aes(x=x1beg,xend=x1beg,y=y.meanbeg.min,yend=y.meanbeg.max,group=af(cluster)),color='black',lty=2) +
    #  geom_segment(data=dfclust.meanbeg,aes(x=y0end,xend=y0end,y=y.meanbeg.min,yend=y.meanbeg.max,group=af(cluster)),color='black',lty=1) +
    geom_rect(data = mybed[grep('VR',mybed$feature),],aes(group=0,x=0,y=0,xmin=beg,xmax=end,ymin=0,ymax=max(mydftemp.meanbeg$meanbeg)),color=rgb(0,0,0,0),fill='black',alpha=0.1) +
    annotate(geom='text',x=700,y=max(mydftemp.meanbeg$y.meanbeg)-100,label='VR') +
    geom_segment(data=dfclust.meanbeg,aes(x=y1beg,xend=y1beg,y=y.meanbeg.min,yend=y.meanbeg.max,group=af(cluster)),color='black',lty=2) +
    geom_text(data=dfclust.meanbeg,aes(x=x0end-100,y=(y.meanbeg.min+y.meanbeg.max)/2,label=cluster,group=af(cluster))) +
    geom_rect(data=mysc.meanbeg[mysc.meanbeg$testType == 'unif',],aes(x=0,y=0,xmin=0  ,ymin=y0,xmax=25 ,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
    geom_rect(data=mysc.meanbeg[mysc.meanbeg$testType == 'norm',],aes(x=0,y=0,xmin=50 ,ymin=y0,xmax=75 ,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
    geom_rect(data=mysc.meanbeg[mysc.meanbeg$testType == 'hots',],aes(x=0,y=0,xmin=100,ymin=y0,xmax=125,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
    annotate(geom='text',x=0  ,y=max(mydftemp.meanbeg$y),label='unif',angle=45,hjust=0,vjust=0,size=4) + 
    annotate(geom='text',x=50 ,y=max(mydftemp.meanbeg$y),label='norm',angle=45,hjust=0,vjust=0,size=4) + 
    annotate(geom='text',x=100,y=max(mydftemp.meanbeg$y),label='hots',angle=45,hjust=0,vjust=0,size=4) + 
    #geom_line(data=mydftemp.meanbeg,aes(x=meanbeg,y=mydftemp.meanbeg[,yvar.meanbeg],group=cluster,alpha=af(meanbeg.hots.sig.coords)),color='red4',lwd=1) +
    geom_point(data=mydftemp.meanbeg[mydftemp.meanbeg$meanbeg.norm.sig.coords != 0,],aes(x=meanbeg,y=mydftemp.meanbeg[mydftemp.meanbeg$meanbeg.norm.sig.coords != 0,yvar.meanbeg],group=cluster),color='black',shape=15,size=0.25) +
    geom_point(data=mydftemp.meanbeg[mydftemp.meanbeg$meanbeg.unif.sig.coords != 0,],aes(x=meanbeg,y=mydftemp.meanbeg[mydftemp.meanbeg$meanbeg.unif.sig.coords != 0,yvar.meanbeg],group=cluster),color='blue2',shape=15,size=0.5) +
    geom_point(data=mydftemp.meanbeg[mydftemp.meanbeg$meanbeg.hots.sig.coords != 0,],aes(x=meanbeg-1,y=mydftemp.meanbeg[mydftemp.meanbeg$meanbeg.hots.sig.coords != 0,yvar.meanbeg],group=cluster),color='red4',shape=18,size=0.5) +
    geom_point(data=mydftemp.meanbeg[mydftemp.meanbeg$meanend.norm.sig.coords != 0,],aes(x=meanend,y=mydftemp.meanbeg[mydftemp.meanbeg$meanend.norm.sig.coords != 0,yvar.meanbeg],group=cluster),color='black',shape=15,size=0.25) +
    geom_point(data=mydftemp.meanbeg[mydftemp.meanbeg$meanend.unif.sig.coords != 0,],aes(x=meanend,y=mydftemp.meanbeg[mydftemp.meanbeg$meanend.unif.sig.coords != 0,yvar.meanbeg],group=cluster),color='blue2',shape=15,size=0.5) +
    geom_point(data=mydftemp.meanbeg[mydftemp.meanbeg$meanend.hots.sig.coords != 0,],aes(x=meanend+1,y=mydftemp.meanbeg[mydftemp.meanbeg$meanend.hots.sig.coords != 0,yvar.meanbeg],group=cluster),color='red4',shape=18,size=0.5) +
    #geom_line(data=mydftemp.meanbeg,aes(x=meanbeg,y=mydftemp.meanbeg[,yvar.meanbeg],group=cluster,alpha=af(meanbeg.unif.sig.coords)),color='black',lwd=0.5) + scale_alpha_manual(values=c('0'=0,'1'=1)) +
    scale_x_continuous(breaks=seq(0,3000,100)) + scale_y_continuous(breaks=seq(0,1600,100)) + theme(legend.position='bottom') +
    #scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) + 
    scale_fill_manual(values=mycolors,breaks=mybreaks,label=mybreaks)
  #p2.meanbeg.c.mydf
  
  p0.meanend.c.mydf = ggplot(mydftemp.meanend,aes(end,beg)) +
    geom_point(aes(color=af(cluster)),pch='.') + theme_bw() +
    coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
    annotate(geom='segment',x=0,y=0,xend=3000,yend=3000) + theme(legend.position='none')
  
  p1.meanend.c.mydf = ggplot(mydftemp.meanend,aes(end,beg)) +
    geom_point(aes(color=af(cluster)),pch='.') + theme_bw() +
    coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
    #  geom_rect(data=unif.meanend,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,fill=af(cluster)),alpha=0.25,color=rgb(0,0,0,0)) +
    #  geom_rect(data=unif.meanbeg,aes(x=0,y=0,xmin=y0end,xmax=y1beg,ymin=x0end,ymax=x1beg,color=af(cluster)),fill=rgb(0,0,0,0)) +
    geom_text(data=dfclust,aes(x=(y1beg+y0end)/2,y=(x1beg+x0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5) +
    geom_rect(data=dfclust,aes(x=y0end,y=x0end,xmin=y0end,ymin=x0end,xmax=y1beg,ymax=x1beg,color=as.factor(cluster)),fill=rgb(0,0,0,0),lwd=1) +
    annotate(geom='segment',x=0,y=0,xend=3000,yend=3000) + theme(legend.position='none')
  
  
  dfclust.meanend = aggregate(mydftemp.meanend$y.meanend,by=list(mydftemp.meanend$cluster),min); colnames(dfclust.meanend) = c("cluster","y.meanend.min")
  dfclust.meanend$y.meanend.max = aggregate(mydftemp.meanend$y.meanend,by=list(mydftemp.meanend$cluster),max)$x
  dfclust.meanend = merge(dfclust.meanend,dfclust,by="cluster")
  # dfclust2 = aggregate(df$beg,by=list(df$cluster),min)
  # colnames(dfclust2) = c('cluster','x0end')
  # dfclust2$x0end = aggregate(df$beg,by=list(df$cluster),min)$x
  # dfclust2$y0end = aggregate(df$end,by=list(df$cluster),min)$x
  # dfclust2$x1beg = aggregate(df$beg,by=list(df$cluster),max)$x
  # dfclust2$y1beg = aggregate(df$end,by=list(df$cluster),max)$x
  # dfclust2$x0beg = dfclust2$x0end
  # dfclust2$y0beg = dfclust2$y1beg
  # dfclust2$x1end = dfclust2$x0end
  # dfclust2$y1end = dfclust2$y1beg
  #posType2 = 'end'
  p2.meanend.c.mydf = ps(mydftemp.meanend,by=posType2,y.var='y.meanend',gp=gp,print=F,group='cluster') +
    geom_segment(data=dfclust.meanend,aes(x=x0end,xend=x0end,y=y.meanend.min,yend=y.meanend.max,group=af(cluster)),color='black',lty=1) +
    #  geom_segment(data=dfclust.meanend,aes(x=x1beg,xend=x1beg,y=y.meanend.min,yend=y.meanend.max,group=af(cluster)),color='black',lty=2) +
    #  geom_segment(data=dfclust.meanend,aes(x=y0end,xend=y0end,y=y.meanend.min,yend=y.meanend.max,group=af(cluster)),color='black',lty=1) +
    geom_rect(data = mybed[grep('VR',mybed$feature),],aes(group=0,x=0,y=0,xmin=beg,xmax=end,ymin=0,ymax=max(mydftemp.meanend$meanend)),color=rgb(0,0,0,0),fill='black',alpha=0.1) +
    annotate(geom='text',x=700,y=max(mydftemp.meanend$y.meanend)-100,label='VR') +
    geom_segment(data=dfclust.meanend,aes(x=y1beg,xend=y1beg,y=y.meanend.min,yend=y.meanend.max,group=af(cluster)),color='black',lty=2) +
    geom_rect(data=mysc.meanend[mysc.meanend$testType == 'unif',],aes(x=0,y=0,xmin=0  ,ymin=y0,xmax=25 ,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
    geom_rect(data=mysc.meanend[mysc.meanend$testType == 'norm',],aes(x=0,y=0,xmin=50 ,ymin=y0,xmax=75 ,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
    geom_rect(data=mysc.meanend[mysc.meanend$testType == 'hots',],aes(x=0,y=0,xmin=100,ymin=y0,xmax=125,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
    annotate(geom='text',x=0  ,y=max(mydftemp.meanend$y),label='unif',angle=45,hjust=0,vjust=0,size=4) + 
    annotate(geom='text',x=50 ,y=max(mydftemp.meanend$y),label='norm',angle=45,hjust=0,vjust=0,size=4) + 
    annotate(geom='text',x=100,y=max(mydftemp.meanend$y),label='hots',angle=45,hjust=0,vjust=0,size=4) +
    geom_text(data=dfclust.meanend,aes(x=x0end-100,y=(y.meanend.min+y.meanend.max)/2,label=cluster,group=af(cluster))) +
    geom_point(data=mydftemp.meanend[mydftemp.meanend$meanend.norm.sig.coords != 0,],aes(x=meanend,y=mydftemp.meanend[mydftemp.meanend$meanend.norm.sig.coords != 0,yvar.meanend],group=cluster),color='black',shape=15,size=0.25) +
    geom_point(data=mydftemp.meanend[mydftemp.meanend$meanend.unif.sig.coords != 0,],aes(x=meanend,y=mydftemp.meanend[mydftemp.meanend$meanend.unif.sig.coords != 0,yvar.meanend],group=cluster),color='blue2',shape=15,size=0.5) +
    geom_point(data=mydftemp.meanend[mydftemp.meanend$meanend.hots.sig.coords != 0,],aes(x=meanend+1,y=mydftemp.meanend[mydftemp.meanend$meanend.hots.sig.coords != 0,yvar.meanend],group=cluster),color='red4',shape=18,size=0.5) +
    geom_point(data=mydftemp.meanend[mydftemp.meanend$meanbeg.norm.sig.coords != 0,],aes(x=meanbeg,y=mydftemp.meanend[mydftemp.meanend$meanbeg.norm.sig.coords != 0,yvar.meanend],group=cluster),color='black',shape=15,size=0.25) +
    geom_point(data=mydftemp.meanend[mydftemp.meanend$meanbeg.unif.sig.coords != 0,],aes(x=meanbeg,y=mydftemp.meanend[mydftemp.meanend$meanbeg.unif.sig.coords != 0,yvar.meanend],group=cluster),color='blue2',shape=15,size=0.5) +
    geom_point(data=mydftemp.meanend[mydftemp.meanend$meanbeg.hots.sig.coords != 0,],aes(x=meanbeg-1,y=mydftemp.meanend[mydftemp.meanend$meanbeg.hots.sig.coords != 0,yvar.meanend],group=cluster),color='red4',shape=18,size=0.5) +
    #geom_line(data=mydftemp.meanend,aes(x=meanend,y=mydftemp.meanend[,yvar.meanend],group=cluster,alpha=af(meanend.unif.sig.coords)),color = 'orange3',lwd=0.5) + scale_alpha_manual(values=c('0'=0,'1'=1)) +
    scale_x_continuous(breaks=seq(0,3000,100)) + scale_y_continuous(breaks=seq(0,1600,100)) + theme(legend.position='bottom') +
    #scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) + 
    scale_fill_manual(values=mycolors,breaks=mybreaks,label=mybreaks)
  #p2.meanend.c.mydf
  
  p3.meanbeg.c.mydf = ggplot(mydftemp.meanbeg,aes(beg)) +
    geom_density(aes(color=af(meanbeg.unif.sig.coords))) + theme_bw() + coord_cartesian(xlim=c(0,3000)) + facet_grid(cluster~.) +
    geom_density(data=mydftemp.meanbeg,aes(x = end,color=af(meanbeg.unif.sig.coords+2)),lty=1) +
    #geom_density(data=mydftemp.meanend,aes(x = end,color=af(meanend.unif.sig.coords+2)),lty=1) +
    scale_color_manual(values=c('1'='red4','0'='orange','3'='blue4','2'='cornflowerblue')) +
    theme(strip.text = element_blank()) +
    scale_x_continuous(breaks=seq(0,3000,100)) +
    theme(legend.position = 'none') + ylab('density') + xlab('R-loop End Position (bp)')
  
  p3.meanend.c.mydf = ggplot(mydftemp.meanend,aes(end)) +
    geom_density(aes(color=af(meanend.unif.sig.coords))) + theme_bw() + coord_cartesian(xlim=c(0,3000))+ facet_grid(cluster~.) +
    geom_density(data=mydftemp.meanend,aes(x =beg,color=af(meanend.unif.sig.coords+2)),lty=1) +
    #  geom_density(data=mydftemp.meanbeg,aes(x = end,color=af(meanend.unif.sig.coords+2)),lty=1) +
    scale_color_manual(values=c('1'='red4','0'='orange','3'='blue4','2'='cornflowerblue')) +
    theme(strip.text = element_blank()) +
    scale_x_continuous(breaks=seq(0,3000,100)) +
    theme(legend.position = 'none') + ylab('density') + xlab('R-loop End Position (bp)')
  
  p4.meanbeg.c.mydf = ggplot(mydftemp.meanbeg,aes(beg)) +
    geom_density(aes(color=af(meanbeg.unif.sig.coords))) + theme_bw() + coord_cartesian(xlim=c(0,3000)) +
    geom_density(data=mydftemp.meanend,aes(x = end,color=af(meanend.unif.sig.coords+2)),lty=1) +
    scale_color_manual(values=c('1'='red4','0'='orange','3'='blue4','2'='cornflowerblue')) +
    scale_x_continuous(breaks=seq(0,3000,100)) +
    theme(legend.position = 'bottom') + ylab('density') + xlab('R-loop Start Position (bp)')
  
  p4.meanend.c.mydf = ggplot(mydftemp.meanend,aes(end)) +
    geom_density(aes(color=af(meanend.unif.sig.coords))) + theme_bw() + coord_cartesian(xlim=c(0,3000)) +
    geom_density(data=mydftemp.meanbeg,aes(x = end,color=af(meanend.unif.sig.coords+2)),lty=1) +
    scale_color_manual(values=c('1'='red4','0'='orange','3'='blue4','2'='cornflowerblue')) +
    scale_x_continuous(breaks=seq(0,3000,100)) +
    theme(legend.position = 'bottom') + ylab('density') + xlab('R-loop End Position (bp)')
  # 
  pdf(paste('./results/',gp$mytitle,'.pdf',sep=''),height=50,width=20)
  grid.arrange(
    p0.meanbeg.c.mydf,p0.meanend.c.mydf,
    p1.meanbeg.c.mydf,p1.meanend.c.mydf,
    p2.meanbeg.c.mydf,p2.meanend.c.mydf,
    p3.meanbeg.c.mydf,p3.meanend.c.mydf,
    p4.meanbeg.c.mydf,p4.meanend.c.mydf,
    nrow=5,ncol=2,
    heights=c(1,1,1,1,0.3)
  )
  dev.off()
  
  pdf(paste('./results/',gp$mytitle,'_meanbeg_heatmaponly.pdf',sep=''),height=10,width=10)
  print(p2.meanbeg.c.mydf)
  dev.off()
  pdf(paste('./results/',gp$mytitle,'_meanend_heatmaponly.pdf',sep=''),height=10,width=10)
  print(p2.meanend.c.mydf)
  dev.off()
  # 
  # pdf('VR17_p2.pdf',height=10,width=10)
  # grid.arrange(p2.meanbeg.c.mydf,p2.meanend.c.mydf,nrow=2,ncol=2)
  # dev.off()
  
  if (myVR == 1) {
    myperc = data.frame()
  }
  df$count = 1
  dfcount = aggregate(df$count,by=list(df$cluster),sum)
  colnames(dfcount) = c('cluster','count')
  dfcount$total = sum(dfcount$count)
  dfcount$perc = as.integer(dfcount$count/ dfcount$total*1000+0.5)/10
  dfcount$VR = myparams$VR
  if (defined(myperc[myperc$VR == myparams$VR,])) {
    myperc = myperc[myperc$VR != myparams$VR,]
  }
  myperc = rbind(myperc,dfcount)
}
#saveRDS(myperc,file='resources/myperc.RDS')
myperc = readRDS(file='resources/myperc.RDS')
# myperc = data.frame()
myperc$GCperc = 0
myperc$GCskew = 0
myperc[grep('^(1||5||9||13||28)$',myperc$VR),]$GCperc = 0.4
myperc[grep('^(2||6||10||14||29)$',myperc$VR),]$GCperc = 0.5
myperc[grep('^(3||7||11||15||16||17||18||19||20||21||22||23||24||25||26||30)$',myperc$VR),]$GCperc = 0.6
myperc[grep('^(4||8||12||27||31)$',myperc$VR),]$GCperc = 0.7
myperc[grep('^(1||2||3||4)$',myperc$VR),]$GCskew = 0.0
myperc[grep('^(5||6||7||8)$',myperc$VR),]$GCskew = 0.1
myperc[grep('^(9||10||11||12)$',myperc$VR),]$GCskew = 0.2
myperc[grep('^(13||14||15||16||17||18||19||20||21||22||23||24||25||26||27)$',myperc$VR),]$GCskew = 0.4
myperc[grep('^(28||29||30||31)$',myperc$VR),]$GCskew = 0.6
myperc$Gclust = 0
myperc[grep('^(15||19||23)$',myperc$VR),]$Gclust = 2
myperc[grep('^(16||20||24)$',myperc$VR),]$Gclust = 3
myperc[grep('^(17||21||25)$',myperc$VR),]$Gclust = 4
myperc$ATskew = 0
myperc[grep('^(15||16||17||18)$',myperc$VR),]$ATskew = 0.0
myperc[grep('^(19||20||21||22)$',myperc$VR),]$ATskew = 0.2
myperc[grep('^(23||24||25||26)$',myperc$VR),]$ATskew = 0.4

pdf("cluster.pdf",width=8,height=10)
ggplot(myperc[(myperc$VR < 15 | myperc$VR == 18 | myperc$VR >= 27) & myperc$cluster != -1,],aes(af(cluster),perc)) +
  geom_bar(aes(fill=af(GCskew)),stat='identity',position='dodge',color='black') +
  theme_bw() + xlab('Cluster') + ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) +
  facet_grid(af(GCperc)~.) +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = 'Set3')                                    

ggplot(myperc[(myperc$VR < 15 | myperc$VR == 18 | myperc$VR >= 27) & myperc$cluster != -1,],aes(af(cluster),perc)) +
  geom_bar(aes(fill=af(GCperc)),stat='identity',position='dodge',color='black') +
  theme_bw() + xlab('Cluster') + ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) +
  facet_grid(af(GCskew)~.) +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = 'Set3')                                    

ggplot(myperc[(myperc$VR >= 15 & myperc$VR <= 26) & myperc$cluster != -1,],aes(af(cluster),perc)) +
  geom_bar(aes(fill=af(Gclust)),stat='identity',position='dodge',color='black') +
  theme_bw() + xlab('Cluster') + ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) +
  facet_grid(af(ATskew)~.) +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = 'Set3')                                    

dev.off()

pdf("cluster_By_GCperc_GCskew_Gclustering.pdf",width=5,height=5)
ggplot(myperc[myperc$VR >= 1 & myperc$VR <= 4 & myperc$cluster != -1,],aes(af(cluster),perc)) +
  geom_bar(aes(fill=af(VR)),stat='identity',position='dodge',color='black') +
  theme_bw() + xlab('Cluster') + ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) + ggtitle('% Rloop Peak Distribution\nat VR #1 to 4') +
  theme(legend.position = 'right') + scale_fill_brewer(palette = 'Set3')                                    
ggplot(myperc[(myperc$VR == 2|myperc$VR == 6|myperc$VR == 10|myperc$VR == 14) & myperc$cluster != -1,],aes(af(cluster),perc)) +
  geom_bar(aes(fill=af(VR)),stat='identity',position='dodge',color='black') +
  theme_bw() + xlab('Cluster') + ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) + ggtitle('% Rloop Peak Distribution\nat VR #2, 6, 10, 14') +
  theme(legend.position = 'right') + scale_fill_brewer(palette = 'Set3')                                    
ggplot(myperc[(myperc$VR >= 15 & myperc$VR <= 18) & myperc$cluster != -1,],aes(af(cluster),perc)) +
  geom_bar(aes(fill=af(VR)),stat='identity',position='dodge',color='black') +
  theme_bw() + xlab('Cluster') + ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) + ggtitle('% Rloop Peak Distribution\nat VR #15-21') +
  theme(legend.position = 'right') + scale_fill_brewer(palette = 'Set3')                                    
dev.off()

myperc3 = myperc[(myperc$cluster >= 1 & myperc$cluster <= 3) | myperc$cluster == 6,]
myperc2 = aggregate(myperc3$count,by=list(myperc3$VR),sum)
colnames(myperc2) = c('VR','total.in.VR')
myperc3 = merge(myperc3,myperc2,by=c('VR'),all=T)
myperc3$perc.in.VR = ai(myperc3$count / myperc3$total.in.VR * 1000 + 0.5)/10



myperc4 = myperc[(myperc$cluster < 1 | myperc$cluster > 3) & myperc$cluster != 6,]
myperc2 = aggregate(myperc4$count,by=list(myperc4$VR),sum)
colnames(myperc2) = c('VR','total.out.VR')
myperc4 = merge(myperc4,myperc2,by=c('VR'),all=T)
myperc4$perc.out.VR = ai(myperc4$count / myperc4$total.out.VR * 1000 + 0.5)/10
#myperc3[is.na(myperc3$total.out.VR),]$total.out.VR = 0



pdf("cluster1to4.pdf",width=5,height=5)

p1a = ggplot(myperc3[(myperc3$VR < 15 | myperc3$VR == 18  | myperc3$VR >= 27) & myperc3$VR != 31, ],aes(af(cluster),perc)) +
  geom_bar(aes(fill=af(GCperc)),stat='identity',position='dodge',color='black') +
  geom_text(aes(label=VR,y=0,group=GCperc),stat='identity',position=position_dodge(width=0.9),size=3,vjust=1) +
  theme_bw() + xlab('Cluster') + ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) +
  facet_grid(.~af(GCskew)) +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = 'Set3')                                    
p1b = ggplot(myperc3[(myperc3$VR < 15 | myperc3$VR == 18  | myperc3$VR >= 27) & myperc3$VR != 31,],aes(af(cluster),perc.in.VR)) +
  geom_bar(aes(fill=af(GCperc)),stat='identity',position='dodge',color='black') +
  geom_text(aes(label=VR,y=0,group=GCperc),stat='identity',position=position_dodge(width=0.9),size=3,vjust=1) +
  theme_bw() + ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) +
  facet_grid(.~af(GCskew)) +
  theme(legend.position = 'bottom') + scale_fill_brewer(palette = 'Set3')                                    

p1a = ggplot(myperc3[(myperc3$VR < 15 | myperc3$VR == 18  | myperc3$VR >= 27) & myperc3$VR != 31,],aes(af(GCperc),perc)) +
  geom_bar(aes(fill=af(cluster)),stat='identity',position='dodge',color='black') +
  geom_text(aes(label=VR,y=0,group=GCperc),stat='identity',position=position_dodge(width=0.9),size=3,vjust=1) +
  theme_bw() + ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) +
  facet_grid(.~af(GCskew)) +
  theme(legend.position = 'right') + scale_fill_brewer(palette = 'Set3')                                    
p1b = ggplot(myperc3[(myperc3$VR < 15 | myperc3$VR == 18  | myperc3$VR >= 27) & myperc3$VR != 31,],aes(af(GCperc),perc.in.VR)) +
  geom_bar(aes(fill=af(cluster)),stat='identity',position='dodge',color='black') +
  geom_text(aes(label=VR,y=0,group=GCperc),stat='identity',position=position_dodge(width=0.9),size=3,vjust=1) +
  theme_bw() +  ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) +
  facet_grid(.~af(GCskew)) +
  theme(legend.position = 'right') + scale_fill_brewer(palette = 'Set3')       


p1a = ggplot(myperc3[(myperc3$VR < 15 | myperc3$VR == 18  | myperc3$VR >= 27) & myperc3$VR != 31,],aes(group=af(GCperc),x=GCperc,y=perc)) +
  geom_boxplot(aes(fill=af(GCperc)),outlier.shape = NA) +
  facet_grid(.~af(cluster)) +
  stat_summary(geom='point',fun=mean,aes(x=GCperc,color=af(GCperc))) +
  stat_summary(geom='line',fun=mean,aes(group=1,x=GCperc),color='black') +#,color=af(GCperc))) +
  #  stat_summary(geom='boxplot') +#,aes(fun=mean)) +
  #  geom_bar(aes(fill=af(cluster)),stat='identity',position='dodge',color='black') +
  #  geom_text(aes(label=VR,y=0,group=af(GCperc)),stat='identity',position=position_dodge(width=0.9),size=3,vjust=1) +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) +
  theme(legend.position = 'right') + 
  scale_color_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')                                    

p1b = ggplot(myperc3[(myperc3$VR < 15 | myperc3$VR == 18  | myperc3$VR >= 27) & myperc3$VR != 31,],aes(group=af(GCperc),x=GCperc,y=perc.in.VR)) +
  geom_boxplot(aes(fill=af(GCperc)),outlier.shape = NA) +
  facet_grid(.~af(cluster)) +
  stat_summary(geom='point',fun=mean,aes(x=GCperc,color=af(GCperc))) +
  stat_summary(geom='line',fun=mean,aes(group=1,x=GCperc),color='black') +#,color=af(GCperc))) +
  #  stat_summary(geom='boxplot') +#,aes(fun=mean)) +
  #  geom_bar(aes(fill=af(cluster)),stat='identity',position='dodge',color='black') +
  #  geom_text(aes(label=VR,y=0,group=af(GCperc)),stat='identity',position=position_dodge(width=0.9),size=3,vjust=1) +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) +
  theme(legend.position = 'right') + 
  scale_color_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')  

p1c = ggplot(myperc4[(myperc4$VR < 15 | myperc4$VR == 18  | myperc4$VR >= 27) & myperc4$VR != 31,],aes(group=af(GCperc),x=GCperc,y=perc)) +
  geom_boxplot(aes(fill=af(GCperc)),outlier.shape = NA) +
  facet_grid(.~af(cluster)) +
  stat_summary(geom='point',fun=mean,aes(x=GCperc,color=af(GCperc))) +
  stat_summary(geom='line',fun=mean,aes(group=1,x=GCperc),color='black') +#,color=af(GCperc))) +
  #  stat_summary(geom='boxplot') +#,aes(fun=mean)) +
  #  geom_bar(aes(fill=af(cluster)),stat='identity',position='dodge',color='black') +
  #  geom_text(aes(label=VR,y=0,group=af(GCperc)),stat='identity',position=position_dodge(width=0.9),size=3,vjust=1) +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,50)) +
  theme(legend.position = 'right') + 
  scale_color_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')  

p1d = ggplot(myperc4[(myperc4$VR < 15 | myperc4$VR == 18  | myperc4$VR >= 27) & myperc4$VR != 31,],aes(group=af(GCperc),x=GCperc,y=perc.out.VR)) +
  geom_boxplot(aes(fill=af(GCperc)),outlier.shape = NA) +
  facet_grid(.~af(cluster)) +
  stat_summary(geom='point',fun=mean,aes(x=GCperc,color=af(GCperc))) +
  stat_summary(geom='line',fun=mean,aes(group=1,x=GCperc),color='black') +#,color=af(GCperc))) +
  #  stat_summary(geom='boxplot') +#,aes(fun=mean)) +
  #  geom_bar(aes(fill=af(cluster)),stat='identity',position='dodge',color='black') +
  #  geom_text(aes(label=VR,y=0,group=af(GCperc)),stat='identity',position=position_dodge(width=0.9),size=3,vjust=1) +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent (%)') +
  coord_cartesian(ylim=c(0,100)) +
  theme(legend.position = 'right') + 
  scale_color_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')  
grid.arrange(p1a,p1b,p1c,p1d,ncol=2,nrow=2)

Bprint('DONE')













# + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())
# 
# test1 = mydftemp.meanbeg[mydftemp.meanbeg$meanbeg.unif.sig.coords == 1 & mydftemp.meanbeg$meanbeg.norm.sig.coords == 1 & mydftemp.meanbeg$cluster == 4,]
# plot(test1$y.meanbeg,type='l')
# test1 = test1[test1$y.meanbeg >= 625 & test1$y.meanbeg < 660,]
# # test1 = test1[order(test1$meanbeg),]; test1$y = seq(1,size(test1))
# #test1 = test1[test1$y < 34,]
# uniform.test(hist(test1$meanbeg,breaks = max(1,sqrt(size(test1))),plot = T))
# shapiro.test(test1$meanbeg)$p.value
# plot(y=test1$y.meanbeg,x=test1$meanbeg,type='l')
# lines(x=seq(min(test1$meanbeg),max(test1$meanbeg),length.out=size(test1)),y=seq(min(test1$y),max(test1$y),length.out=size(test1)),col='red4')
# boxplot(mydftemp.meanbeg[mydftemp.meanbeg$meanbeg.unif.sig.coords == 1 & mydftemp.meanbeg$meanbeg.norm.sig.coords == 1 & mydftemp.meanbeg$cluster == 2,]$meanbeg.unif.pval,mydftemp.meanbeg[mydftemp.meanbeg$meanbeg.unif.sig.coords == 1 & mydftemp.meanbeg$meanbeg.norm.sig.coords == 1 & mydftemp.meanbeg$cluster == 2,]$meanbeg.norm.pval)
# 
# p4 = ggplot(mydftemp.meanbeg,aes(x=mydftemp.meanbeg$y.meanbeg,y=1)) +
# #  geom_line(aes(y=mydftemp.meanbeg$meanbeg.unif.pval),color='red4') +
#   stat_smooth(aes(y=mydftemp.meanbeg$meanbeg.unif.pval),color='red4') +
#   stat_smooth(aes(y=mydftemp.meanbeg$meanbeg.hots.pval),color='green4') +
#   stat_smooth(aes(y=mydftemp.meanbeg$meanbeg.unif.pval),color='red4') +
#   #  geom_line(aes(y=mydftemp.meanbeg$meanbeg.norm.pval),color='blue4') +
#   stat_smooth(aes(y=mydftemp.meanbeg$meanbeg.norm.pval),color='blue4') +
#   facet_grid(cluster~.)
# p4

