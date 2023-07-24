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


myparams = list(
  gene     = 'T7_init_VR_20',
  treat    = 'C',
  peaktype = 'TOP',
  thres    = 0,
  VR       = 20
)
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
myfa$amp.beg = mybed[mybed$feature == "FW_Barcode",]$beg
myfa$amp.end = mybed[mybed$feature == "RV_Barcode",]$end
myfa$amp.seq = gsub(paste('^.{',myfa$amp.beg,'}','.(.+)','.{',myfa$amp.end,'}','$'),"\\1",myfa$seq,perl=T)

#myfa$seq = gsub()
# Get PEAKS
df = slice_df(PEAKS, myparams = myparams,myparams_regex=myparams_regex)
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
#dfclustVR20 = dfclust

# Get Cluster from VR20
dfclust = dfclustVR20

head(dfclust)
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

positionTypes = c('mean')
testTypes = c('unif','hots','norm')
mylist = do_distributionTests(df, positionTypes = positionTypes, testTypes = testTypes, gp = gp)
mydf = mylist$df
mysc = mylist$misc$sig.coords

mydf = annot.pvar(mydf,mysc)
mysc = re.sc(mydf,mysc)
mydf = annot.pvar(mydf,mysc)
head(mydf)
head(mysc)
#ggplot.xy(df = mydf,dfclust = dfclust)


positionType1 = 'meanbeg'
yvar.meanbeg = paste('y.',positionType1,sep='')
mydftemp.meanbeg = mydf[order(mydf$cluster,mydf$meanbeg,mydf$meanend),]; mydftemp.meanbeg[,yvar.meanbeg] = seq(1,dim(mydftemp.meanbeg)[1])
mysc.meanbeg = mysc[mysc$positionType == positionType1,]
unif.meanbeg = mysc[mysc$positionType == positionType1 & mysc$testType == 'unif',]

positionType2 = 'meanend'
yvar.meanend = paste('y.',positionType2,sep='')
mydftemp.meanend = mydf[order(mydf$cluster,mydf$meanend,mydf$meanbeg),]; mydftemp.meanend[,yvar.meanend] = seq(1,dim(mydftemp.meanend)[1])
mysc.meanend = mysc[mysc$positionType == positionType2,]
unif.meanend = mysc[mysc$positionType == positionType2 & mysc$testType == 'unif',]

mycolors = c(brewer.pal(9,"Set1"),brewer.pal(9,"Set3"))
mybreaks = seq(0,length(mycolors)-1)
  
p1.meanbeg.c.mydf = ggplot(mydftemp.meanbeg,aes(meanbeg,meanend)) +
  geom_point(aes(color=af(cluster)),pch='.') + theme_bw() +
  coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
  geom_rect(data=unif.meanbeg,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,color=af(cluster)),fill=rgb(0,0,0,0)) +
  geom_rect(data=unif.meanend,aes(x=0,y=0,xmin=y0end,xmax=y1beg,ymin=x0end,ymax=x1beg,fill=af(cluster)),alpha=0.25,color=rgb(0,0,0,0)) +
  geom_text(data=dfclust,aes(x=(x1beg+x0end)/2,y=(y1beg+y0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5) +
  geom_rect(data=dfclust,aes(x=x0end,y=y0end,xmin=x0end,ymin=y0end,xmax=x1beg,ymax=y1beg,color=as.factor(cluster)),fill=rgb(0,0,0,0),lwd=1) +
  annotate(geom='segment',x=0,y=0,xend=2500,yend=2500) + theme(legend.position='none')

p2.meanbeg.c.mydf = ps(mydftemp.meanbeg,by=positionType1,y.var=yvar.meanbeg,gp=gp,print=F,group='cluster') +
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
  scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) + scale_fill_manual(values=mycolors,breaks=mybreaks,label=mybreaks)
p2.meanbeg.c.mydf

p1.meanend.c.mydf = ggplot(mydftemp.meanend,aes(meanend,meanbeg)) +
  geom_point(aes(color=af(cluster)),pch='.') + theme_bw() +
  coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
  geom_rect(data=unif.meanend,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,fill=af(cluster)),alpha=0.25,color=rgb(0,0,0,0)) +
  geom_rect(data=unif.meanbeg,aes(x=0,y=0,xmin=y0end,xmax=y1beg,ymin=x0end,ymax=x1beg,color=af(cluster)),fill=rgb(0,0,0,0)) +
  geom_text(data=dfclust,aes(x=(y1beg+y0end)/2,y=(x1beg+x0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5) +
  geom_rect(data=dfclust,aes(x=y0end,y=x0end,xmin=y0end,ymin=x0end,xmax=y1beg,ymax=x1beg,color=as.factor(cluster)),fill=rgb(0,0,0,0),lwd=1) +
  annotate(geom='segment',x=0,y=0,xend=2500,yend=2500) + theme(legend.position='none')

p2.meanend.c.mydf = ps(mydftemp.meanend,by=positionType2,y.var='y.meanend',gp=gp,print=F,group='cluster') +
  geom_rect(data=mysc.meanend[mysc.meanend$testType == 'unif',],aes(x=0,y=0,xmin=0  ,ymin=y0,xmax=25 ,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
  geom_rect(data=mysc.meanend[mysc.meanend$testType == 'norm',],aes(x=0,y=0,xmin=50 ,ymin=y0,xmax=75 ,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
  geom_rect(data=mysc.meanend[mysc.meanend$testType == 'hots',],aes(x=0,y=0,xmin=100,ymin=y0,xmax=125,ymax=y1,fill=af(cluster),group=af(cluster)),color='black') +
  annotate(geom='text',x=0  ,y=max(mydftemp.meanend$y),label='unif',angle=45,hjust=0,vjust=0,size=4) + 
  annotate(geom='text',x=50 ,y=max(mydftemp.meanend$y),label='norm',angle=45,hjust=0,vjust=0,size=4) + 
  annotate(geom='text',x=100,y=max(mydftemp.meanend$y),label='hots',angle=45,hjust=0,vjust=0,size=4) + 
  #geom_line(data=mydftemp.meanend,aes(x=meanend,y=mydftemp.meanend[,yvar.meanend],group=cluster,alpha=af(meanend.hots.sig.coords)),color='red4',lwd=1) +

  geom_point(data=mydftemp.meanend[mydftemp.meanend$meanend.norm.sig.coords != 0,],aes(x=meanend,y=mydftemp.meanend[mydftemp.meanend$meanend.norm.sig.coords != 0,yvar.meanend],group=cluster),color='black',shape=15,size=0.25) +
  geom_point(data=mydftemp.meanend[mydftemp.meanend$meanend.unif.sig.coords != 0,],aes(x=meanend,y=mydftemp.meanend[mydftemp.meanend$meanend.unif.sig.coords != 0,yvar.meanend],group=cluster),color='blue2',shape=15,size=0.5) +
  geom_point(data=mydftemp.meanend[mydftemp.meanend$meanend.hots.sig.coords != 0,],aes(x=meanend+1,y=mydftemp.meanend[mydftemp.meanend$meanend.hots.sig.coords != 0,yvar.meanend],group=cluster),color='red4',shape=18,size=0.5) +
  geom_point(data=mydftemp.meanend[mydftemp.meanend$meanbeg.norm.sig.coords != 0,],aes(x=meanbeg,y=mydftemp.meanend[mydftemp.meanend$meanbeg.norm.sig.coords != 0,yvar.meanend],group=cluster),color='black',shape=15,size=0.25) +
  geom_point(data=mydftemp.meanend[mydftemp.meanend$meanbeg.unif.sig.coords != 0,],aes(x=meanbeg,y=mydftemp.meanend[mydftemp.meanend$meanbeg.unif.sig.coords != 0,yvar.meanend],group=cluster),color='blue2',shape=15,size=0.5) +
  geom_point(data=mydftemp.meanend[mydftemp.meanend$meanbeg.hots.sig.coords != 0,],aes(x=meanbeg-1,y=mydftemp.meanend[mydftemp.meanend$meanbeg.hots.sig.coords != 0,yvar.meanend],group=cluster),color='red4',shape=18,size=0.5) +
  #geom_line(data=mydftemp.meanend,aes(x=meanend,y=mydftemp.meanend[,yvar.meanend],group=cluster,alpha=af(meanend.unif.sig.coords)),color = 'orange3',lwd=0.5) + scale_alpha_manual(values=c('0'=0,'1'=1)) +
  scale_x_continuous(breaks=seq(0,3000,100)) + scale_y_continuous(breaks=seq(0,1600,100)) + theme(legend.position='bottom') +
  scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) + scale_fill_manual(values=mycolors,breaks=mybreaks,label=mybreaks)
p2.meanend.c.mydf

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

pdf(paste('./results/',gp$mytitle,'.pdf',sep=''),height=40,width=20)
grid.arrange(p1.meanbeg.c.mydf,p1.meanend.c.mydf,
             p2.meanbeg.c.mydf,p2.meanend.c.mydf,
             p3.meanbeg.c.mydf,p3.meanend.c.mydf,
             p4.meanbeg.c.mydf,p4.meanend.c.mydf,
             nrow=4,ncol=2,
             heights=c(1,1,1,0.3)
)
dev.off()
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

