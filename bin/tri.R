source("lib/tri_lib.R")
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

#[order(df[df$beg >= 400 & df$beg <= 500,'beg']),] #,n=20)
#myfa$seq = gsub()
# Get PEAKS
df = slice_df(PEAKS, myparams = myparams,myparams_regex=myparams_regex)
# dim(df[df$beg %in% beg2$beg,])
# dim(df[df$end %in% end2$end,])
# df$begorig = df$beg
# df$endorig = df$end
# df = merge(df,beg2,by='beg')
# df = merge(df,end2,by='end')
df$cluster = 0
df$y = seq(1,dim(df)[1])


df = df[order(df$beg, df$end),]; df$y = seq(1,dim(df)[1])
df$meanbeg = df$beg
#for (i in seq(1,(max(df$y)),gp$stepsmooth)) {
#  myrange = min(i+gp$windowsmooth,max(df$y))
#  df$meanbeg[i] = mean(df[seq(i,myrange),]$beg)
#}

df = df[order(df$end, df$beg),]; df$y = seq(1,dim(df)[1])
df$meanend = df$end
#for (i in seq(1,(max(df$y)),gp$stepsmooth)) {
#  myrange = min(i+gp$windowsmooth,max(df$y))
#  df$meanend[i] = mean(df[seq(i,myrange),]$end)
#}
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

#grid.arrange(p1.beg,p1.end,p2.beg,p2.end)


# Get Cluster
dfclust = slice_CLUSTS(df,CLUSTS)
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

 dfclust = dfclustVR20
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

pdf('VR17_p2.pdf',height=10,width=10)
grid.arrange(p2.meanbeg.c.mydf,p2.meanend.c.mydf,nrow=2,ncol=2)
dev.off()

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
pdf("cluster1to4.pdf",width=5,height=5)
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

