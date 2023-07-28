
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
mypval = data.frame()
for (clust in unique(myperc$cluster)) {
  temp1 = myperc[myperc$clust == clust,]
  for (GCperc in unique(myperc$GCperc)) {
    temp2 = temp1[temp1$GCperc == GCperc,]
    for (GCperc2 in unique(myperc$GCperc)) {
      temp3 = temp1[temp1$GCperc == GCperc2,]
      p = wilcox.test(temp2$perc,temp3$perc)$p.value/2
      temp4 = data.frame(cluster=clust,GCperc=GCperc,GCperc2=GCperc2,p=p)
      mypval = rbind(mypval,temp4)
    }
  }
}
mypval
mypval$pchar = ''
mypval[mypval$p <= 0.05,]$pchar = "*"
mypval[mypval$p <= 0.01,]$pchar = "*"
mypval[mypval$p <= 0.001,]$pchar = "*"
mypval2 = mypval[mypval$GCperc2 == 0.4,]


mypvalGCskew = data.frame()
for (clust in unique(myperc$cluster)) {
  temp1 = myperc[myperc$clust == clust,]
  for (GCskew in unique(myperc$GCskew)) {
    temp2 = temp1[temp1$GCskew == GCskew,]
    for (GCskew2 in unique(myperc$GCskew)) {
      temp3 = temp1[temp1$GCskew == GCskew2,]
      p = wilcox.test(temp2$perc,temp3$perc)$p.value/2
      temp4 = data.frame(cluster=clust,GCskew=GCskew,GCskew2=GCskew2,p=p)
      mypvalGCskew = rbind(mypvalGCskew,temp4)
    }
  }
}
mypvalGCskew
mypvalGCskew$pchar = ''
mypvalGCskew[mypvalGCskew$p <= 0.05,]$pchar = "*"
mypvalGCskew[mypvalGCskew$p <= 0.01,]$pchar = "*"
mypvalGCskew[mypvalGCskew$p <= 0.001,]$pchar = "*"
mypvalGCskew2 = mypvalGCskew[mypvalGCskew$GCskew2 == 0.0,]

p1aGCperc = ggplot(myperc[(myperc$VR < 15 | myperc$VR == 18  | myperc$VR >= 27) & myperc$VR != 31,],aes(x=GCperc,y=perc)) +
  facet_grid(.~paste('Cluster',af(cluster))) +
  geom_boxplot(aes(fill=af(GCperc)),outlier.shape = NA) +
#  stat_smooth(method = 'lm',formula=y~log(x+0.01),color='black') +#,color=af(GCperc))) +
  stat_smooth(method = 'lm',formula=y~x,color='blue4') +#,color=af(GCskew))) +
  geom_text(data=mypval2,aes(y=70,label=pchar)) +
  geom_point(aes(color = af(GCskew))) +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent R-loop Reads in each cluster (%)') +
  xlab("GC Content") +
  coord_cartesian(ylim=c(0,75)) +
  theme(legend.position = 'right') + 
  scale_x_continuous(breaks = c(0.4,0.5,0.6,0.7)) +
  scale_color_brewer(palette = 'Reds','GC Skew') +
  scale_fill_brewer(palette = 'Greens') +
  ggtitle("% R-loop Distribution in each cluster\nVR init series (1-31)\nGrouped by GC %")

p1aGCskew = ggplot(myperc[(myperc$VR < 15 | myperc$VR == 18  | myperc$VR >= 27) & myperc$VR != 31,],aes(x=GCskew,y=perc)) +
  facet_grid(.~paste('Cluster',af(cluster))) +
  geom_boxplot(aes(fill=af(GCskew)),outlier.shape = NA) +
#  stat_smooth(method = 'lm',formula=y~log(x+0.01),color='black') +#,color=af(GCskew))) +
  stat_smooth(method = 'lm',formula=y~x,color='blue4') +#,color=af(GCskew))) +
  geom_text(data=mypvalGCskew2,aes(y=70,label=pchar)) +
  geom_point(aes(color = af(GCperc))) +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent R-loop Reads in each cluster (%)') +
  xlab("GC Skew") +
  coord_cartesian(ylim=c(0,75)) +
  theme(legend.position = 'right') + 
  scale_color_brewer(palette = 'Reds','GC Cont') +
  scale_fill_brewer(palette = 'Greens') +
  ggtitle("% R-loop Distribution in each cluster\nVR init series (1-31)\nGrouped by GC Skew")

pdf("Rloop_Distr_VR_Init_Series.pdf",width=10,height=10)
grid.arrange(p1aGCperc,p1aGCskew,nrow=2)
dev.off()


p0aGCperc = ggplot(myperc[(myperc$VR < 15 | myperc$VR == 18  | myperc$VR >= 27) & myperc$VR != 31,],aes(x=GCperc,y=perc)) +
  facet_grid(paste('GC Skew',GCskew)~paste('Cluster',af(cluster))) +
  #stat_smooth(method = 'lm',formula=y~log(x+0.01),color='black') +#,color=af(GCperc))) +
  stat_smooth(method = 'lm',formula=y~x,color='blue2',lwd=0.2,fill=rgb(0,0.7,0,0.1)) +#,color=af(GCperc))) +
  geom_point(color='black') +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent R-loop Reads in each cluster (%)') +
  xlab("GC content") +
  coord_cartesian(ylim=c(0,75)) +
  theme(legend.position = 'right') + 
  scale_x_continuous(breaks = c(0.4,0.5,0.6,0.7)) +
  scale_color_brewer(palette = 'Reds','GC Skew') +
  scale_fill_brewer(palette = 'Greens') +
  ggtitle("% R-loop Distribution in each cluster\nVR init series (1-31)")

p0aGCskew = ggplot(myperc[(myperc$VR < 15 | myperc$VR == 18  | myperc$VR >= 27) & myperc$VR != 31,],aes(x=GCskew,y=perc)) +
  facet_grid(paste('GC Cont',GCperc)~paste('Cluster',af(cluster))) +
#  stat_smooth(method = 'lm',formula=y~log(x+0.01),color='black') +#,color=af(GCperc))) +
  stat_smooth(method = 'lm',formula=y~x,color='blue2',lwd=0.2,fill=rgb(0.7,0,0,0.1)) +#,color=af(GCperc))) +
  geom_point(color='black') +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent R-loop Reads in each cluster (%)') +
  xlab("GC skew") +
  coord_cartesian(ylim=c(0,75)) +
  theme(legend.position = 'right') + 
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6)) +
  scale_color_brewer(palette = 'Reds','GC Cont') +
  scale_fill_brewer(palette = 'Greens') +
  ggtitle("% R-loop Distribution in each cluster\nVR init series (1-31)")

pdf("Rloop_Distr_VR_Init_Series_All.pdf",width=10,height=20)
grid.arrange(p0aGCperc,p0aGCskew,nrow=2)
dev.off()


p2aGCperc = ggplot(myperc[(myperc$VR < 15 | myperc$VR == 18  | myperc$VR >= 27) & myperc$VR != 31,],aes(x=GCperc,y=perc)) +
  #facet_grid(paste('GC Skew',GCskew)~paste('Cluster',af(cluster))) +
  facet_grid(.~paste('Cluster',af(cluster))) +
  #stat_smooth(method = 'lm',formula=y~log(x+0.01),color='black') +#,color=af(GCperc))) +
  stat_smooth(aes(group=af(GCskew),color=af(GCskew)),fill=rgb(1,1,1,alpha = 0),method = 'lm',formula=y~x) +#,color='blue2',lwd=0.2,fill=rgb(0,0.7,0,0.1)) +#,color=af(GCperc))) +
  geom_point(aes(color=af(GCskew))) +#color='black') +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent R-loop Reads in each cluster (%)') +
  xlab("GC content") +
  coord_cartesian(ylim=c(0,75)) +
  theme(legend.position = 'right') + 
  scale_x_continuous(breaks = c(0.4,0.5,0.6,0.7)) +
  #scale_color_distiller(palette = 'Reds',direction = 1) +
  scale_color_manual(breaks = c(0,0.1,0.2,0.4,0.6),values = brewer.pal(9,"Reds")[4:9]) +#brewer(palette = 'Reds','GC Skew') +
  scale_fill_brewer(palette = 'Greens') +
  ggtitle("% R-loop Distribution in each cluster\nVR init series (1-31)")


p2aGCskew = ggplot(myperc[(myperc$VR < 15 | myperc$VR == 18  | myperc$VR >= 27) & myperc$VR != 31,],aes(x=GCskew,y=perc)) +
  #facet_grid(paste('GC Skew',GCskew)~paste('Cluster',af(cluster))) +
  facet_grid(.~paste('Cluster',af(cluster))) +
  #stat_smooth(method = 'lm',formula=y~log(x+0.01),color='black') +#,color=af(GCskew))) +
  stat_smooth(aes(group=af(GCperc),color=af(GCperc)),fill=rgb(1,1,1,alpha = 0),method = 'lm',formula=y~x) +#,color='blue2',lwd=0.2,fill=rgb(0,0.7,0,0.1)) +#,color=af(GCskew))) +
  geom_point(aes(color=af(GCperc))) +#color='black') +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent R-loop Reads in each cluster (%)') +
  xlab("GC skew") +
  coord_cartesian(ylim=c(0,75)) +
  theme(legend.position = 'right') + 
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6)) +
  #scale_color_distiller(palette = 'Reds',direction = 1) +
  scale_color_manual(breaks = c(0.4,0.5,0.6,0.7),values = brewer.pal(9,"Greens")[5:9]) +#brewer(palette = 'Reds','GC Skew') +
  scale_fill_brewer(palette = 'Greens') +
  ggtitle("% R-loop Distribution in each cluster\nVR init series (1-31)")

pdf("Rloop_Distr_VR_Init_Series_All2.pdf",width=20,height=10)
grid.arrange(p2aGCperc,p2aGCskew,nrow=2)
dev.off()

myperc5 = myperc
myperc5$GCfactor = (myperc5$GCskew/2 + 0.5)*myperc5$GCperc

p4a = ggplot(myperc5[(myperc5$VR < 15 | myperc5$VR == 18  | myperc5$VR >= 27) & myperc5$VR != 31,],aes(x=GCperc,y=perc)) +
 # facet_grid(paste('GC Cont',GCperc)~paste('Cluster',af(cluster))) +
 #  geom_boxplot(aes(fill=af(GCskew)),outlier.shape = NA) +
 stat_smooth(method = 'lm',formula=y~x,color='blue4') +#,color=af(GCskew))) +

  facet_grid(.~paste('Cluster',af(cluster))) +
  #stat_smooth(method = 'lm',formula=y~log(x+0.01),color='black') +#,color=af(GCfactor))) +
  #stat_smooth(aes(group=af(GCskew),color=af(GCskew)),fill=rgb(1,1,1,alpha = 0),method = 'lm',formula=y~x) +#,color='blue2',lwd=0.2,fill=rgb(0,0.7,0,0.1)) +#,color=af(GCfactor))) +
  geom_point(aes(color=af(GCskew))) +#color='black') +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent R-loop Reads in each cluster (%)') +
  xlab("GC Factor") +
  coord_cartesian(ylim=c(0,75)) +
  # theme(legend.position = 'right') + 
  # scale_x_continuous(breaks = c(0.4,0.5,0.6,0.7)) +
  # #scale_color_distiller(palette = 'Reds',direction = 1) +
  # scale_color_manual(breaks = c(0,0.1,0.2,0.4,0.6),values = brewer.pal(9,"Reds")[4:9]) +#brewer(palette = 'Reds','GC Skew') +
  # scale_fill_brewer(palette = 'Greens') +
  ggtitle("% R-loop Distribution in each cluster\nVR init series (1-31)")

 
p4b =  ggplot(myperc5[(myperc5$VR < 15 | myperc5$VR == 18  | myperc5$VR >= 27) & myperc5$VR != 31,],aes(x=GCskew,y=perc)) +
 # facet_grid(paste('GC Skew',GCskew)~paste('Cluster',af(cluster))) +
 # geom_boxplot(aes(fill=af(GCperc)),outlier.shape = NA) +
  stat_smooth(method = 'lm',formula=y~x,color='blue4') +#,color=af(GCskew))) +
  facet_grid(.~paste('Cluster',af(cluster))) +
  #stat_smooth(method = 'lm',formula=y~log(x+0.01),color='black') +#,color=af(GCfactor))) +
 #stat_smooth(aes(group=af(GCperc),color=af(GCperc)),fill=rgb(1,1,1,alpha = 0),method = 'lm',formula=y~x) +#,color='blue2',lwd=0.2,fill=rgb(0,0.7,0,0.1)) +#,color=af(GCfactor))) +
  geom_point(aes(color=af(GCperc))) +#color='black') +
  theme_bw() + theme(panel.grid = element_blank()) +ylab('Percent R-loop Reads in each cluster (%)') +
  xlab("GC Factor") +
   coord_cartesian(ylim=c(0,75)) +
  # theme(legend.position = 'right') + 
  # scale_x_continuous(breaks = c(0.4,0.5,0.6,0.7)) +
  # #scale_color_distiller(palette = 'Reds',direction = 1) +
  # scale_color_manual(breaks = c(0,0.1,0.2,0.4,0.6),values = brewer.pal(9,"Reds")[4:9]) +#brewer(palette = 'Reds','GC Skew') +
  # scale_fill_brewer(palette = 'Greens') +
  ggtitle("% R-loop Distribution in each cluster\nVR init series (1-31)")

pdf("Rloop_Distr_VR_Init_Series_All4.pdf",width=20,height=10)
grid.arrange(p4a, p4b, nrow=2)
dev.off()
