myparams = list(
  gene     = '^T7_init',
  treat    = 'C',
  peaktype = 'TOP',
  thres    = 0,
  VR       = 'any'
  
)

myparams_regex = list(
  gene = T,
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
myfa$amp.seq = myfa$seq

# Get PEAKS
df = slice_df(PEAKS, myparams = myparams,myparams_regex=myparams_regex)

df$cluster = 0
df$y = seq(1,dim(df)[1])


df = df[order(df$beg, df$end),]; df$y = seq(1,dim(df)[1])
df$meanbeg = df$beg


df = df[order(df$end, df$beg),]; df$y = seq(1,dim(df)[1])
df$meanend = df$end

#smdf(df)



df = df[order(df$cluster,df$beg, df$end),]; df$y = seq(1,dim(df)[1])
df$index = seq(1,dim(df)[1])
p1.beg = ps(df,'beg',gp=gp,print=F,group='cluster')

df = df[order(df$cluster,df$end, df$beg),]; df$y = seq(1,dim(df)[1])
p1.end = ps(df,'end',gp=gp,print=F,group='cluster')

#grid.arrange(p1.beg,p1.end)


# Get Cluster
#dfclust = slice_CLUSTS(df,CLUSTS)
dfclust = readRDS('dfclustVR20.RDS')
df = get_cluster(df = df, dfclust = dfclust)

df = df[order(df$cluster,df$beg, df$end),]; df$y = seq(1,dim(df)[1])
p1.beg.c1 = ps(df,'beg',gp=gp,print=F,group='cluster')

df = df[order(df$cluster,df$end, df$beg),]; df$y = seq(1,dim(df)[1])
p1.end.c1 = ps(df,'end',gp=gp,print=F,group='cluster')

df = df[order(df$beg, df$end),]; df$y = seq(1,dim(df)[1])
p1.beg.c2 = ps(df,'beg',gp=gp,print=F,group='cluster')

df = df[order(df$end, df$beg),]; df$y = seq(1,dim(df)[1])
p1.end.c2 = ps(df,'end',gp=gp,print=F,group='cluster')


grid.arrange(p1.beg.c1,p1.end.c1,p1.beg.c2,p1.end.c2)

p1.clust.length = ggplot(df,aes())

annotate_GCperc_GCskew_to_VR = function(df) {
  df$GCperc = 0
  df$GCskew = 0
  df[grep('^(1||5||9||13||28)$',df$VR),]$GCperc = 0.4
  df[grep('^(2||6||10||14||29)$',df$VR),]$GCperc = 0.5
  df[grep('^(3||7||11||15||16||17||18||19||20||21||22||23||24||25||26||30)$',df$VR),]$GCperc = 0.6
  df[grep('^(4||8||12||27||31)$',df$VR),]$GCperc = 0.7
  df[grep('^(1||2||3||4)$',df$VR),]$GCskew = 0.0
  df[grep('^(5||6||7||8)$',df$VR),]$GCskew = 0.1
  df[grep('^(9||10||11||12)$',df$VR),]$GCskew = 0.2
  df[grep('^(13||14||15||16||17||18||19||20||21||22||23||24||25||26||27)$',df$VR),]$GCskew = 0.4
  df[grep('^(28||29||30||31)$',df$VR),]$GCskew = 0.6
  df$Gclust = 0
  df$ATskew = 0
  if (defined(  df[grep('^(15||19||23)$',df$VR),])) {
    df[grep('^(15||19||23)$',df$VR),]$Gclust = 2
    df[grep('^(16||20||24)$',df$VR),]$Gclust = 3
    df[grep('^(17||21||25)$',df$VR),]$Gclust = 4
    df[grep('^(15||16||17||18)$',df$VR),]$ATskew = 0.0
    df[grep('^(19||20||21||22)$',df$VR),]$ATskew = 0.2
    df[grep('^(23||24||25||26)$',df$VR),]$ATskew = 0.4
  }
  
  to_return = df
  return(to_return)
}

df = annotate_GCperc_GCskew_to_VR(df)

temp = df
temp = temp[(temp$VR < 15 | temp$VR == 18  | temp$VR >= 27) & temp$VR != 31,]
temp2 = subset(temp,select=c('cluster','beg','end','VR','GCperc','GCskew'))#,id.vars = c('cluster','VR','beg','end'))
temp2$type1 = temp2$GCperc
temp2$type2 = temp2$GCskew
temp2$type = 'GCperc'
#temp2$type2 = ai(100*(0.4+(0.3*temp2$type2 / 0.6)))/100
temp3 = temp2
temp3$type1 = temp3$GCskew
temp3$type2 = temp3$GCperc
#temp3$type1 = ai(100*(0.4+(0.3*temp3$type1 / 0.6)))/100
temp3$type = 'GCskew'
temp4 = rbind(temp2,temp3)

pdf('testdensity.pdf',height=30,width=30)
p1 = ggplot(temp4[temp4$cluster >= 1 & temp4$cluster <= 6,],aes(beg)) +
  geom_density(data=temp4[temp4$cluster >= 1 & temp4$cluster <= 6 & temp4$type == 'GCperc',],alpha=0.5,aes(color=af(paste(type)))) +
  geom_density(data=temp4[temp4$cluster >= 1 & temp4$cluster <= 6 & temp4$type == 'GCskew',],alpha=0.5,aes(color=af(paste(type)))) +
  facet_grid(paste(cluster,type2)~af(type1),scales = 'free_y')
p2 = ggplot(temp4[temp4$cluster >= 1 & temp4$cluster <= 6,],aes(beg)) +
  geom_density(data=temp4[temp4$cluster >= 1 & temp4$cluster <= 6 & temp4$type == 'GCperc',],alpha=0.5,aes(color=af(paste(type)))) +
  facet_grid(paste(cluster,type2)~af(type1),scales = 'free_y')
p3 = ggplot(temp4[temp4$cluster >= 1 & temp4$cluster <= 6,],aes(beg)) +
  geom_density(data=temp4[temp4$cluster >= 1 & temp4$cluster <= 6 & temp4$type == 'GCskew',],alpha=0.5,aes(color=af(paste(type)))) +
  facet_grid(paste(cluster,type2)~af(type1),scales = 'free_y')
grid.arrange(p1,p2,p3,nrow=1,ncol=3)
dev.off()

p1 = ggplot(temp4[temp4$cluster >= 1 & temp4$cluster <= 6,],aes(beg)) +
  geom_density(data=temp4[temp4$cluster >= 1 & temp4$cluster <= 6 & temp4$type == 'GCperc',],alpha=0.5,aes(color=af(paste(type)))) +
  geom_density(data=temp4[temp4$cluster >= 1 & temp4$cluster <= 6 & temp4$type == 'GCskew',],alpha=0.5,aes(color=af(paste(type)))) +
  facet_grid(paste(cluster,type2)~af(type1),scales = 'free_y')
p2 = ggplot(temp4[temp4$cluster >= 1 & temp4$cluster <= 6,],aes(beg)) +
  geom_density(data=temp4[temp4$cluster >= 1 & temp4$cluster <= 6 & temp4$type == 'GCperc',],alpha=0.5,aes(group=af(type1),color=af(type1))) +
  scale_color_brewer(palette = 'Greens') +
  facet_grid(paste(cluster,type2)~.,scales = 'free_y')
p3 = ggplot(temp4[temp4$cluster >= 1 & temp4$cluster <= 6,],aes(beg)) +
  geom_density(data=temp4[temp4$cluster >= 1 & temp4$cluster <= 6 & temp4$type == 'GCskew',],alpha=0.5,aes(group=af(type1),color=af(type1))) +
  scale_color_brewer(palette = 'Reds') +
  facet_grid(paste(cluster,type2)~.,scales = 'free_y')
pdf('testdensity.pdf',height=30,width=30)
grid.arrange(p1,p2,p3,nrow=1,ncol=3)
dev.off()



p1 = ggplot(temp4[temp4$cluster >= 1 & temp4$cluster <= 6,],aes(beg)) +
  geom_density(data=temp4[temp4$cluster >= 1 & temp4$cluster <= 6 & temp4$type == 'GCperc',],alpha=0.5,aes(color=af(paste(type)))) +
  geom_density(data=temp4[temp4$cluster >= 1 & temp4$cluster <= 6 & temp4$type == 'GCskew',],alpha=0.5,aes(color=af(paste(type)))) +
  facet_grid(paste(cluster,type2)~af(type1),scales = 'free_y')
p2 = ggplot(temp2[temp2$cluster >= 1 & temp2$cluster <= 6,],aes(beg)) +
  geom_density(data=temp2[temp2$cluster >= 1 & temp2$cluster <= 6,],alpha=0.5,aes(group=af(GCperc),color=af(GCperc))) +
  scale_color_manual(values=c('#41ab5d','#238b45','#006d2c','#00441b')) +
#  scale_color_brewer(palette = 'BuGn',direction = 1) +
  facet_grid(paste('Cluster',cluster,'\nGCskew',GCskew)~.,scales = 'free_y')
p3 = ggplot(temp2[temp2$cluster >= 1 & temp2$cluster <= 6,],aes(beg)) +
  geom_density(data=temp2[temp2$cluster >= 1 & temp2$cluster <= 6,],alpha=0.5,aes(group=af(GCskew),color=af(GCskew))) +
#  scale_color_brewer(palette = 'RdBu',direction = -1) +
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
  facet_grid(paste('Cluster',cluster,'\nGCperc',GCperc)~.,scales = 'free_y')
pdf('testdensity.pdf',height=30,width=30)
grid.arrange(p1,p2,p3,nrow=1,ncol=3,widths=c(1,0.2,0.2))
dev.off()


head(temp4)
clusters = unique(temp4$cluster)
VRs = unique(temp4$VR)
temp7 = data.frame()
for (cluster in clusters[order(clusters)]) {
  print(cluster)
  for (VR in VRs[order(VRs)]) {
    print(VR)
    temp5 = temp4[temp4$cluster == cluster & temp4$VR == VR,]
    temp6 = as.data.frame(matrix(apply(temp5,1,myag),byrow = TRUE,nrow=nrow(temp5),ncol=3000))
    colnames(temp6) = seq(1,3000)
    temp6.sum = apply(temp6,2,sum)
    temp6.mean = apply(temp6,2,mean)
    temp6.total = dim(temp6)[1]
    temp6.sum.smooth.4 = smooth.spline(temp6.sum,spar=0.4)
    temp6.sum.smooth.6 = smooth.spline(temp6.sum,spar=0.6)
    temp6.sum.smooth.8 = smooth.spline(temp6.sum,spar=0.8)
    perc.sum=ai(10000*temp6.sum/temp6.total)/100
    perc.sum.smooth.4=smooth.spline(perc.sum,spar=0.4)
    perc.sum.smooth.6=smooth.spline(perc.sum,spar=0.6)
    perc.sum.smooth.8=smooth.spline(perc.sum,spar=0.8)
        # plot(perc.sum,type='l');lines(perc.sum.smooth,col=2,lwd=2)
    temp6 = data.frame(cluster=cluster,
                       VR=VR,
                       pos=colnames(temp6),
                       count=temp6.sum,
                       count.smooth.4=temp6.sum.smooth.4$y,
                       count.smooth.6=temp6.sum.smooth.6$y,
                       count.smooth.8=temp6.sum.smooth.8$y,
                       total=temp6.total,
                       perc.smooth.4=perc.sum.smooth.4$y,
                       perc.smooth.6=perc.sum.smooth.6$y,
                       perc.smooth.8=perc.sum.smooth.8$y,
                       perc=ai(10000*temp6.sum/temp6.total)/100
    )
    temp7 = rbind(temp7,temp6)
  }
}
 
myag = function(x) {
  temp6 = rep(1,3000)
  temp6[seq(x[2],x[3])] = 1
  return(temp6)
}
y <- c(1, 1, 19:1)
plot(y, main = "misbehaviour of \"3RSR\"", col.main = 3)
lines(sm.3RS(y))
lines(smooth(y))
lines(smooth(y, "3RSR"), col = 3, lwd = 2)  # the horror


temp7 = annotate_GCperc_GCskew_to_VR(temp7)
    #temp6.sumsmooth = smooth.spline(temp6.sum,spar=0.8)

# plot(smooth.spline(temp6.sum),type='l')
# lines(smooth.spline(smooth.spline(temp6.sum),spar=0.2),col=3,lwd=2)
# lines(smooth.spline(smooth.spline(temp6.sum),spar=0.8,cv=NA),col=8,lwd=2)
# plot(temp6.sum, main = "misbehaviour of \"3RSR\"", col.main = 3)
# lines(sm.3RS(temp6.sum))
# lines(smooth(temp6.sum))
# lines(smooth(temp6.sum, "3RSR"), col = 3, lwd = 2)  # the horror
#p1 = 

temp7$pos = ai(temp7$pos)
# ggplot(temp7[temp7$GCperc == 0.5 & temp7$cluster == 2,],aes(group=af(GCperc),ai(pos),perc.smooth.8)) +
#   geom_line(aes(group=af(GCperc),color=af(GCperc))) +
#   facet_grid(paste('Cluster',cluster,'\nGCskew',GCskew)~.,scales = 'free_y')
#  geom_smooth(aes(group=af(GCperc),color=af(GCperc))) +


p0.GCperc = ggplot(temp7[temp7$cluster >= 1 & temp7$cluster <= 6 & temp7$pos >= 500 & temp7$pos <= 1200,],aes(group=af(GCperc),x=pos,y=perc.smooth.4)) +
  stat_summary(geom='line',aes(color=af(GCperc)),fun=mean) +
  annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=100,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=95,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=100,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=95,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=100,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=95,label='SNRPN') + theme_bw() +
#  geom_line(aes(group=af(GCperc),color=af(GCperc)))+
  scale_color_manual(values=c('#41ab5d','#238b45','#006d2c','#00441b')) +
  ylab("Average R-loop Peak Density (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("Average R-loop Peak Density\nGrouped by GC content") +
   coord_cartesian(ylim=c(0,100)) +
 facet_grid(paste('Cluster',cluster)~.)


p0.GCskew = ggplot(temp7[temp7$cluster >= 1 & temp7$cluster <= 6 & temp7$pos >= 500 & temp7$pos <= 1200,],aes(group=af(GCskew),x=pos,y=perc.smooth.4)) +
  stat_summary(geom='line',aes(color=af(GCskew)),fun=mean) +
  annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=100,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=95,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=100,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=95,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=100,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=95,label='SNRPN') + theme_bw() +
#  geom_line(aes(group=af(GCskew),color=af(GCskew)))+
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
   coord_cartesian(ylim=c(0,100)) +
  ylab("Average R-loop Peak Density (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("Average R-loop Peak Density\nGrouped by GC Skew") +
 facet_grid(paste('Cluster',cluster)~.)


p1.GCperc = ggplot(temp7[temp7$cluster == 1 & temp7$pos >= 500 & temp7$pos <= 1200,],aes(group=af(GCperc),x=pos,y=perc.smooth.4)) +
  geom_line(aes(group=af(GCperc),color=af(GCperc)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=100,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=95,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=100,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=95,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=100,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=95,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#41ab5d','#238b45','#006d2c','#00441b')) +
  ylab("R-loop Peak Density (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density of Cluster 1\nGrouped by GC content") +  facet_grid(paste('Cluster',cluster,'\nGCskew',GCskew)~.,scales = 'free_y')

p2.GCperc = ggplot(temp7[temp7$cluster == 2 & temp7$pos >= 500 & temp7$pos <= 1200,],aes(group=af(GCperc),x=pos,y=perc.smooth.4)) +
  geom_line(aes(group=af(GCperc),color=af(GCperc)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=100,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=95,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=100,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=95,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=100,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=95,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#41ab5d','#238b45','#006d2c','#00441b')) +
  ylab("R-loop Peak Density (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density of Cluster 2\nGrouped by GC content") +  facet_grid(paste('Cluster',cluster,'\nGCskew',GCskew)~.,scales = 'free_y')

p3.GCperc = ggplot(temp7[temp7$cluster == 3 & temp7$pos >= 500 & temp7$pos <= 1200,],aes(group=af(GCperc),x=pos,y=perc.smooth.4)) +
  geom_line(aes(group=af(GCperc),color=af(GCperc)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=100,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=95,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=100,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=95,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=100,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=95,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#41ab5d','#238b45','#006d2c','#00441b')) +
  ylab("R-loop Peak Density (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density of Cluster 3\nGrouped by GC content") +  facet_grid(paste('Cluster',cluster,'\nGCskew',GCskew)~.,scales = 'free_y')

p6.GCperc = ggplot(temp7[temp7$cluster == 6 & temp7$pos >= 500 & temp7$pos <= 1200,],aes(group=af(GCperc),x=pos,y=perc.smooth.4)) +
  geom_line(aes(group=af(GCperc),color=af(GCperc)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=100,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=95,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=100,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=95,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=100,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=95,label='SNRPN') + theme_bw() +
  ylab("R-loop Peak Density (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density of Cluster 4\nGrouped by GC content") +  facet_grid(paste('Cluster',cluster,'\nGCskew',GCskew)~.,scales = 'free_y')
  # facet_grid(GCskew~cluster,scales='free_y')#paste('Cluster',cluster,'\nGCskew',GCskew)~.,scales = 'free_y')

p1.GCskew = ggplot(temp7[temp7$cluster == 1 & temp7$pos >= 500 & temp7$pos <= 1200,],aes(group=af(GCskew),x=pos,y=perc.smooth.4)) +
  geom_line(aes(group=af(GCskew),color=af(GCskew)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=100,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=95,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=100,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=95,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=100,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=95,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
  ylab("R-loop Peak Density (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density of Cluster 1\nGrouped by GC Skew") + 
  facet_grid(paste('Cluster',cluster,'\nGCperc',GCperc)~.,scales = 'free_y')

p2.GCskew = ggplot(temp7[temp7$cluster == 2 & temp7$pos >= 500 & temp7$pos <= 1200,],aes(group=af(GCskew),x=pos,y=perc.smooth.4)) +
  geom_line(aes(group=af(GCskew),color=af(GCskew)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=100,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=95,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=100,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=95,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=100,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=95,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
  ylab("R-loop Peak Density (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density of Cluster 2\nGrouped by GC Skew") + 
  facet_grid(paste('Cluster',cluster,'\nGCperc',GCperc)~.,scales = 'free_y')

p3.GCskew = ggplot(temp7[temp7$cluster == 3 & temp7$pos >= 500 & temp7$pos <= 1200,],aes(group=af(GCskew),x=pos,y=perc.smooth.4)) +
  geom_line(aes(group=af(GCskew),color=af(GCskew)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=100,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=95,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=100,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=95,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=100,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=95,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
  ylab("R-loop Peak Density (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density of Cluster 3\nGrouped by GC Skew") + 
  facet_grid(paste('Cluster',cluster,'\nGCperc',GCperc)~.,scales = 'free_y')

p6.GCskew = ggplot(temp7[temp7$cluster == 6 & temp7$pos >= 500 & temp7$pos <= 1200,],aes(group=af(GCskew),x=pos,y=perc.smooth.4)) +
  geom_line(aes(group=af(GCskew),color=af(GCskew)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=100,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=95,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=100,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=95,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=100,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=95,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
  ylab("R-loop Peak Density (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density of Cluster 6\nGrouped by GC Skew") + 
  facet_grid(paste('Cluster',cluster,'\nGCperc',GCperc)~.,scales = 'free_y')

pdf('testdensity3.pdf',height=7,width=10)
grid.arrange(p1.GCperc,p1.GCskew,nrow=1,ncol=2,widths=c(0.2,0.2))
grid.arrange(p2.GCperc,p2.GCskew,nrow=1,ncol=2,widths=c(0.2,0.2))
grid.arrange(p3.GCperc,p3.GCskew,nrow=1,ncol=2,widths=c(0.2,0.2))
grid.arrange(p6.GCperc,p6.GCskew,nrow=1,ncol=2,widths=c(0.2,0.2))
grid.arrange(p0.GCperc,p0.GCskew,nrow=1,ncol=2,widths=c(0.2,0.2))
dev.off()

temp7divby = temp7[temp7$GCperc == 0.4,]
temp7divby = subset(temp7divby,select=c('cluster','GCskew','pos','count','perc','perc.smooth.4'))
colnames(temp7divby) = c('cluster','GCskew','pos','countdivby','percdivby','percdivby.smooth.4')

temp8 = merge(temp7,temp7divby,by=c('cluster','GCskew','pos'))
temp8[temp8$perc.smooth.4 < 0,]$perc.smooth.4 = 0
temp8[temp8$percdivby.smooth.4 < 0,]$percdivby.smooth.4 = 0
temp8$perc.smooth.4.GCskewdiff = (abs(temp8$perc.smooth.4)+10) - (abs(temp8$percdivby.smooth.4)+10)


temp7divby = temp7[temp7$GCskew == 0,]
temp7divby = subset(temp7divby,select=c('cluster','GCperc','pos','count','perc','perc.smooth.4'))
colnames(temp7divby) = c('cluster','GCperc','pos','countdivby2','percdivby2','percdivby2.smooth.4')
temp8 = merge(temp8,temp7divby,by=c('cluster','GCperc','pos'))
#temp8[temp8$perc.smooth.4 < 0,]$perc.smooth.4 = 0
temp8[temp8$percdivby2.smooth.4 < 0,]$percdivby2.smooth.4 = 0
temp8$perc.smooth.4.GCpercdiff = (abs(temp8$perc.smooth.4)+10) - (abs(temp8$percdivby2.smooth.4)+10)


p0.GCperc.divby = ggplot(temp8[temp8$cluster >= 1 & temp8$cluster <= 6 & temp8$pos >= 500 & temp8$pos <= 1200,],aes(group=af(GCperc),x=pos,y=perc.smooth.4.GCskewdiff)) +
  stat_summary(geom='line',aes(color=af(GCperc)),fun=mean) +
  annotate(geom='rect',xmin=587,xmax=605,ymin=-70,ymax=70,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=65,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=-70,ymax=70,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=65,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=-70,ymax=70,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=65,label='SNRPN') + theme_bw() +
#  geom_line(aes(group=af(GCperc),color=af(GCperc)))+
  scale_color_manual(values=c('#41ab5d','#238b45','#006d2c','#00441b')) +
  ylab("Average Diff. of R-loop Peak Density\nvs. 40% GC content") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("Average Difference of R-loop Peak Density\nvs. 40% GC content\nGrouped by GC content") +
  #(Density current - Density of R-loop Peak at 40% GC content") +
   coord_cartesian(ylim=c(-70,70)) +
 facet_grid(paste('Cluster',cluster)~.)


p0.GCskew.divby = ggplot(temp8[temp8$cluster >= 1 & temp8$cluster <= 6 &temp8$pos >= 500 & temp8$pos <= 1200,],aes(group=af(GCskew),x=pos,y=perc.smooth.4.GCpercdiff)) +
  stat_summary(geom='line',aes(color=af(GCskew)),fun=mean) +
  annotate(geom='rect',xmin=587,xmax=605,ymin=-70,ymax=70,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=65,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=-70,ymax=70,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=65,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=-70,ymax=70,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=65,label='SNRPN') + theme_bw() +
#  geom_line(aes(group=af(GCskew),color=af(GCskew)))+
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
  ylab("Average Diff. of R-loop Peak Density (%)\nvs. 0 GC skew") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("Average Difference of R-loop Peak Density\nvs. 0 GC skew\nGrouped by GC Skew") +
  #(Density current - Density of R-loop Peak at 0 GC skew") +
   coord_cartesian(ylim=c(-70,70)) +
 facet_grid(paste('Cluster',cluster)~.)

p1.GCperc.divby = ggplot(temp8[temp8$cluster == 1 & temp8$pos >= 500 & temp8$pos <= 1200,],aes(group=af(GCperc),x=pos,y=perc.smooth.4.GCskewdiff)) +
  geom_line(aes(group=af(GCperc),color=af(GCperc)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=-70,ymax=70,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=65,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=-70,ymax=70,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=65,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=-70,ymax=70,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=65,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#41ab5d','#238b45','#006d2c','#00441b')) +
   coord_cartesian(ylim=c(-70,70)) +
  ylab("R-loop Peak Density Difference vs. 40% GC content (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density Difference vs. 40% GC content of Cluster 1\nGrouped by GC Content") + 
 facet_grid(paste('Cluster',cluster,'\nGCskew',GCskew)~.)

p2.GCperc.divby = ggplot(temp8[temp8$cluster == 2 & temp8$pos >= 500 & temp8$pos <= 1200,],aes(group=af(GCperc),x=pos,y=perc.smooth.4.GCskewdiff)) +
  geom_line(aes(group=af(GCperc),color=af(GCperc)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=-70,ymax=70,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=65,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=-70,ymax=70,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=65,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=-70,ymax=70,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=65,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#41ab5d','#238b45','#006d2c','#00441b')) +
   coord_cartesian(ylim=c(-70,70)) +
  ylab("R-loop Peak Density Difference vs. 40% GC content (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density Difference vs. 40% GC content of Cluster 2\nGrouped by GC Content") + 
 facet_grid(paste('Cluster',cluster,'\nGCskew',GCskew)~.)

p3.GCperc.divby = ggplot(temp8[temp8$cluster == 3 & temp8$pos >= 500 & temp8$pos <= 1200,],aes(group=af(GCperc),x=pos,y=perc.smooth.4.GCskewdiff)) +
  geom_line(aes(group=af(GCperc),color=af(GCperc)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=-70,ymax=70,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=65,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=-70,ymax=70,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=65,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=-70,ymax=70,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=65,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#41ab5d','#238b45','#006d2c','#00441b')) +
   coord_cartesian(ylim=c(-70,70)) +
  ylab("R-loop Peak Density Difference vs. 40% GC content (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density Difference vs. 40% GC content of Cluster 3\nGrouped by GC Content") + 
 facet_grid(paste('Cluster',cluster,'\nGCskew',GCskew)~.)

p6.GCperc.divby = ggplot(temp8[temp8$cluster == 6 & temp8$pos >= 500 & temp8$pos <= 1200,],aes(group=af(GCperc),x=pos,y=perc.smooth.4.GCskewdiff)) +
  geom_line(aes(group=af(GCperc),color=af(GCperc)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=-70,ymax=70,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=65,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=-70,ymax=70,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=65,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=-70,ymax=70,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=65,label='SNRPN') + theme_bw() +
   coord_cartesian(ylim=c(-70,70)) +
 scale_color_manual(values=c('#41ab5d','#238b45','#006d2c','#00441b')) +
  ylab("R-loop Peak Density Difference vs. 40% GC content (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density Difference vs. 40% GC content of Cluster 6\nGrouped by GC Content") + 
  facet_grid(paste('Cluster',cluster,'\nGCskew',GCskew)~.)

p1.GCskew.divby = ggplot(temp8[temp8$cluster == 1 & temp8$pos >= 500 & temp8$pos <= 1200,],aes(group=af(GCskew),x=pos,y=perc.smooth.4.GCpercdiff)) +
  geom_line(aes(group=af(GCskew),color=af(GCskew)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=-70,ymax=70,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=65,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=-70,ymax=70,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=65,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=-70,ymax=70,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=65,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
   coord_cartesian(ylim=c(-70,70)) +
  ylab("R-loop Peak Density Difference vs. 0 GC skew (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density Difference vs. 0 GC skew of Cluster 1\nGrouped by GC Skew") + 
 facet_grid(paste('Cluster',cluster,'\nGCperc',GCperc)~.)

p2.GCskew.divby = ggplot(temp8[temp8$cluster == 2 & temp8$pos >= 500 & temp8$pos <= 1200,],aes(group=af(GCskew),x=pos,y=perc.smooth.4.GCpercdiff)) +
  geom_line(aes(group=af(GCskew),color=af(GCskew)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=-70,ymax=70,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=65,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=-70,ymax=70,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=65,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=-70,ymax=70,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=65,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
  ylab("R-loop Peak Density Difference vs. 0 GC skew (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density Difference vs. 0 GC skew of Cluster 2\nGrouped by GC Skew") + 
   coord_cartesian(ylim=c(-70,70)) +
 facet_grid(paste('Cluster',cluster,'\nGCperc',GCperc)~.)

p3.GCskew.divby = ggplot(temp8[temp8$cluster == 3 & temp8$pos >= 500 & temp8$pos <= 1200,],aes(group=af(GCskew),x=pos,y=perc.smooth.4.GCpercdiff)) +
  geom_line(aes(group=af(GCskew),color=af(GCskew)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=-70,ymax=70,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=65,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=-70,ymax=70,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=65,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=-70,ymax=70,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=65,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
  ylab("R-loop Peak Density Difference vs. 0 GC skew (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density Difference vs. 0 GC skew of Cluster 3\nGrouped by GC Skew") + 
   coord_cartesian(ylim=c(-70,70)) +
 facet_grid(paste('Cluster',cluster,'\nGCperc',GCperc)~.)

p6.GCskew.divby = ggplot(temp8[temp8$cluster == 6 & temp8$pos >= 500 & temp8$pos <= 1200,],aes(group=af(GCskew),x=pos,y=perc.smooth.4.GCpercdiff)) +
  geom_line(aes(group=af(GCskew),color=af(GCskew)))+
  annotate(geom='rect',xmin=587,xmax=605,ymin=-70,ymax=70,fill=rgb(0.5,0.5,0.5,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=595,y=65,label='T7') +
  annotate(geom='rect',xmin=641,xmax=841,ymin=-70,ymax=70,fill=rgb(0,0.25,0,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=741,y=65,label='VR') +
  annotate(geom='rect',xmin=842,xmax=1317,ymin=-70,ymax=70,fill=rgb(0,0,0.25,0.1),color=rgb(0,0,0,0.15)) +
  annotate(geom='text',x=1000,y=65,label='SNRPN') + theme_bw() +
  scale_color_manual(values=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')) +
   ylab("R-loop Peak Density Difference vs. 0 GC skew (%)") +
  xlab("Position in Plasmid (bp)") +
  ggtitle("R-loop Peak Density Difference vs. 0 GC skew of Cluster 6\nGrouped by GC Skew") + 
 coord_cartesian(ylim=c(-70,70)) +
  facet_grid(paste('Cluster',cluster,'\nGCperc',GCperc)~.)


pdf('testdensity4.pdf',height=7,width=10)
grid.arrange(p1.GCperc.divby,p1.GCskew.divby,nrow=1,ncol=2,widths=c(0.2,0.2))
grid.arrange(p2.GCperc.divby,p2.GCskew.divby,nrow=1,ncol=2,widths=c(0.2,0.2))
grid.arrange(p3.GCperc.divby,p3.GCskew.divby,nrow=1,ncol=2,widths=c(0.2,0.2))
grid.arrange(p6.GCperc.divby,p6.GCskew.divby,nrow=1,ncol=2,widths=c(0.2,0.2))
grid.arrange(p0.GCperc.divby,p0.GCskew.divby,nrow=1,ncol=2,widths=c(0.2,0.2))
dev.off()


test0 = temp8[temp8$cluster == 2 & temp8$pos >= 700 & temp8$pos <= 710 & temp8$GCskew == 0.0,]$perc.smooth.4
test1 = temp8[temp8$cluster == 2 & temp8$pos >= 700 & temp8$pos <= 710 & temp8$GCskew == 0.1,]$perc.smooth.4
test2 = temp8[temp8$cluster == 2 & temp8$pos >= 700 & temp8$pos <= 710 & temp8$GCskew == 0.2,]$perc.smooth.4
test4 = temp8[temp8$cluster == 2 & temp8$pos >= 700 & temp8$pos <= 710 & temp8$GCskew == 0.4,]$perc.smooth.4
test6 = temp8[temp8$cluster == 2 & temp8$pos >= 700 & temp8$pos <= 710 & temp8$GCskew == 0.6,]$perc.smooth.4
wilcox.test(test0,test1)$p.value/2
wilcox.test(test0,test2)$p.value/2
wilcox.test(test0,test4)$p.value/2
wilcox.test(test0,test6)$p.value/2

#aaa







ggplot(temp4[temp4$cluster >= 1 & temp4$cluster <= 6,],aes(af(type1),beg)) +
  theme_bw() + theme(panel.grid=element_blank()) +
#  geom_boxplot(aes(group=af(paste(type,type1)),fill=af(type)),outlier.shape = NA) +
  stat_smooth(aes(group=af(paste(type)),color=af(paste(type)),fill=af(paste(type))),method = 'lm',formula=y~x) +#log(x+0.1)) +#,color=af(GCskew))) +
  stat_summary(geom='point',fun=mean,aes(group=af(paste(type)),color=af(paste(type)))) +
  stat_smooth(aes(y=end,group=af(paste(type)),color=af(paste(type)),fill=af(paste(type))),method = 'lm',formula=y~x) +#log(x+0.1)) +#,color=af(GCskew))) +
  stat_summary(geom='point',fun=mean,aes(y=end,group=af(paste(type)),color=af(paste(type)))) +
  # stat_smooth(aes(group=af(paste(type,type2)),color=af(paste(type,type2)),fill=af(paste(type,type2))),method = 'lm',formula=y~log(x+0.1)) +#,color=af(GCskew))) +
  # stat_summary(geom='point',fun=mean,aes(group=af(paste(type,type2)),color=af(paste(type,type2)))) +
  facet_grid(af(type2)~cluster) +
  coord_cartesian(ylim=c(600,1200))

p1.clust.beg.length.GCperc = ggplot(
  temp[temp$cluster >= 1 & temp$cluster <= 6,],aes(GCskew,(beg))) +
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_boxplot(aes(fill=af(GCskew))) +
#  geom_density(aes(color=af(GCskew))) +
  stat_smooth(method = 'lm',formula=y~log(x+0.1),color='blue4') +#,color=af(GCskew))) +
  #facet_grid(.~cluster) +
  facet_grid(GCperc~cluster) +
  coord_cartesian(ylim=c(600,1200))

p1.clust.beg.length.GCskew = ggplot(
  temp[temp$cluster >= 1 & temp$cluster <= 6,],aes(GCperc,(beg))) +
  geom_boxplot(aes(fill=af(GCperc))) +
  stat_smooth(method = 'lm',formula=y~log(x+0.1),color='blue4') +#,color=af(GCskew))) +
  theme_bw() + theme(panel.grid=element_blank()) +
  # facet_grid(.~cluster) +
  facet_grid(GCskew~cluster) +
  coord_cartesian(ylim=c(600,1200))


p1.clust.end.length.GCperc = ggplot(
  temp[temp$cluster >= 1 & temp$cluster <= 6,],aes(GCskew,(end))) +
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_boxplot(aes(fill=af(GCskew))) +
  stat_smooth(method = 'lm',formula=y~log(x+0.1),color='blue4') +#,color=af(GCskew))) +
  # facet_grid(.~cluster) +
  facet_grid(GCperc~cluster) +
  coord_cartesian(ylim=c(600,1200))

p1.clust.end.length.GCskew = ggplot(
  temp[temp$cluster >= 1 & temp$cluster <= 6,],aes(GCperc,(end))) +
  geom_boxplot(aes(fill=af(GCperc))) +
  stat_smooth(method = 'lm',formula=y~log(x+0.1),color='blue4') +#,color=af(GCskew))) +
  theme_bw() + theme(panel.grid=element_blank()) +
  # facet_grid(.~cluster) +
  facet_grid(GCskew~cluster) +
  coord_cartesian(ylim=c(600,1200))
grid.arrange(p1.clust.beg.length.GCperc,p1.clust.beg.length.GCskew,
             p1.clust.end.length.GCperc,p1.clust.end.length.GCskew,
             nrow=2,ncol=2)



p1.clust.beg.length.GCperc = ggplot(
  temp[temp$cluster >= 1 & temp$cluster <= 6,],aes(GCskew,(beg/(GCperc)))) +
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_boxplot(aes(fill=af(GCskew))) +
#  geom_density(aes(color=af(GCskew))) +
  stat_smooth(method = 'lm',formula=y~log(x+0.1),color='blue4') +#,color=af(GCskew))) +
  facet_grid(.~cluster) +
  #facet_grid(GCperc~cluster)
  coord_cartesian(ylim=c(600,3000))

p1.clust.beg.length.GCskew = ggplot(
  temp[temp$cluster >= 1 & temp$cluster <= 6,],aes(GCperc,(beg/(GCskew/0.6*0.3+0.4)))) +
  geom_boxplot(aes(fill=af(GCperc))) +
  stat_smooth(method = 'lm',formula=y~log(x+0.1),color='blue4') +#,color=af(GCskew))) +
  theme_bw() + theme(panel.grid=element_blank()) +
  facet_grid(.~cluster) +
  #facet_grid(GCperc~cluster)
  coord_cartesian(ylim=c(600,3000))


p1.clust.end.length.GCperc = ggplot(
  temp[temp$cluster >= 1 & temp$cluster <= 6,],aes(GCskew,(end/(GCperc)))) +
  theme_bw() + theme(panel.grid=element_blank()) +
  geom_boxplot(aes(fill=af(GCskew))) +
  stat_smooth(method = 'lm',formula=y~log(x+0.1),color='blue4') +#,color=af(GCskew))) +
  facet_grid(.~cluster) +
  #facet_grid(GCperc~cluster)
  coord_cartesian(ylim=c(600,3000))

p1.clust.end.length.GCskew = ggplot(
  temp[temp$cluster >= 1 & temp$cluster <= 6,],aes(GCperc,(end/(GCskew/0.6*0.3+0.4)))) +
  geom_boxplot(aes(fill=af(GCperc))) +
  stat_smooth(method = 'lm',formula=y~log(x+0.1),color='blue4') +#,color=af(GCskew))) +
  theme_bw() + theme(panel.grid=element_blank()) +
  facet_grid(.~cluster) +
  #facet_grid(GCperc~cluster)
  coord_cartesian(ylim=c(600,3000))
grid.arrange(p1.clust.beg.length.GCperc,p1.clust.beg.length.GCskew,
             p1.clust.end.length.GCperc,p1.clust.end.length.GCskew,
             nrow=2,ncol=2)
