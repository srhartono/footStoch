#saveRDS(dm1,file="dm1.sgRNA.all.RDS")
mydir = "/Users/mitochy/Work/Project/Ethan//triclust/"
#dm1 = readRDS(file=paste(mydir,"dm1.sgRNA.all.RDS",sep=''))
dm1 = readRDS(file=paste(mydir,"dm1.all.RDS",sep=''))
RDS = dir(paste(mydir,"./final3/",sep=''),"*final3.RDS")
RDS = paste(mydir,"./final3/",RDS,sep='')
#RDS = paste("./final3/",RDS,sep="")
for (i in 1:length(RDS)) {
   temp = readRDS(RDS[i])
   print(i)
   if (i == 1) {
      final0 = temp
   } else {
      final0 = rbind(final0,temp)
   }
}
head(final0)
head(dm1)

get_dmtemp_not100 = function(final0,type='beg') {
   #   dm = data.frame(pos=seq(1,100))
   final0$count = 1
#   dmtemp = aggregate(final0$count,by=list(final0$x0beg,final0$y1beg),sum); colnames(dmtemp) = c("beg","end","count"); dm = dmtemp; #merge(dm,dmtemp,by="pos",all = T)
   if (type == 'beg') {
      dmtemp = aggregate(final0$count,by=list(final0$x0end,final0$y0end),sum); colnames(dmtemp) = c("beg","end","count"); dm = dmtemp; #merge(dm,dmtemp,by="pos",all = T)
   } else if (type == 'end') {
      dmtemp = aggregate(final0$count,by=list(final0$x1beg,final0$y1beg),sum); colnames(dmtemp) = c("beg","end","count"); dm = dmtemp; #merge(dm,dmtemp,by="pos",all = T)
   } else if (type == 'mid') {
      dmtemp = aggregate(final0$count,by=list(final0$x0end,final0$y1beg),sum); colnames(dmtemp) = c("beg","end","count"); dm = dmtemp; #merge(dm,dmtemp,by="pos",all = T)
   } else {
      dmtemp = aggregate(final0$count,by=list(final0$x0end,final0$y0end,final0$x1beg,final0$y1beg),sum); colnames(dmtemp) = c("beg","end",'beg2','end2',"count"); dm = dmtemp; #merge(dm,dmtemp,by="pos",all = T)
   }
   dm[is.na(dm)] = 0
   return(dm)
}

get_dmtemp_want_not100 = function(final0,VRwant='any',treatwant='any',genewant='any',type='beg') {
   #   dm = data.frame(pos=seq(1,100))
   final0$count = 1
   final0$VR = 0
   if (length(grep("VR_",final0$gene)) > 0) {
      final0$VR = as.numeric(as.character(gsub("^.+VR_([0-9]+)$","\\1",final0$gene,perl=T)))
   }
   print(dim(final0)[1])
   print(head(final0$VR))
   if (VRwant != 'any') {
      final0 = final0[final0$VR == VRwant,]
   }
   print(dim(final0)[1])
   if (treatwant != 'any') {
      final0 = final0[grep(paste("^",treatwant,sep=""),final0$treat,perl=T),]
   }
   print(dim(final0)[1])
   if (genewant != 'any') {
      final0 = final0[final0$gene == genewant,]
   }
   print(dim(final0)[1])
   #   dmtemp = aggregate(final0$count,by=list(final0$x0beg,final0$y1beg),sum); colnames(dmtemp) = c("beg","end","count"); dm = dmtemp; #merge(dm,dmtemp,by="pos",all = T)
   if (type == 'beg') {
      dmtemp = aggregate(final0$count,by=list(final0$x0end,final0$y0end),sum); colnames(dmtemp) = c("beg","end","count"); dm = dmtemp; #merge(dm,dmtemp,by="pos",all = T)
   } else if (type == 'end') {
      dmtemp = aggregate(final0$count,by=list(final0$x1beg,final0$y1beg),sum); colnames(dmtemp) = c("beg","end","count"); dm = dmtemp; #merge(dm,dmtemp,by="pos",all = T)
   } else if (type == 'mid') {
      dmtemp = aggregate(final0$count,by=list(final0$x0end,final0$y1beg),sum); colnames(dmtemp) = c("beg","end","count"); dm = dmtemp; #merge(dm,dmtemp,by="pos",all = T)
   } else {
      dmtemp = aggregate(final0$count,by=list(final0$x0end,final0$y0end,final0$x1beg,final0$y1beg),sum); colnames(dmtemp) = c("beg","end",'beg2','end2',"count"); dm = dmtemp; #merge(dm,dmtemp,by="pos",all = T)
   }
   dm[is.na(dm)] = 0
   dm$VR = 1
   return(dm)
}



dmbeg.not100 = get_dmtemp_not100(final0,'beg')
dmend.not100 = get_dmtemp_not100(final0,'end')
dmmid.not100 = get_dmtemp_not100(final0,'mid')
dmall.not100 = get_dmtemp_not100(final0,'all')
dmbeg.not100$perc = as.integer(dmbeg.not100$count/sum(dmbeg.not100$count)*1000)/10
dmend.not100$perc = as.integer(dmend.not100$count/sum(dmend.not100$count)*1000)/10
dmmid.not100$perc = as.integer(dmmid.not100$count/sum(dmmid.not100$count)*1000)/10

dmall.not100 = get_dmtemp_want_not100(final0,VRwant=VRwant,treatwant=treatwant,type='all')
dmall.not100$perc = as.integer(dmall.not100$count/sum(dmall.not100$count)*1000)/10
dmall.not100[dmall.not100$beg == dmall.not100$end,]$end = dmall.not100[dmall.not100$beg == dmall.not100$end,]$end + 1
dmall.not100[dmall.not100$beg2 == dmall.not100$end2,]$beg2 = dmall.not100[dmall.not100$beg2 == dmall.not100$end2,]$beg2 - 1

# ggplot(dmbeg.not100[dmbeg.not100$perc >= 2, ],aes(beg*divby,end*divby)) +
#    geom_point(aes(size=count/sum(count)),color='red4') +
#    geom_point(data=dmend.not100[dmend.not100$perc >= 2,],aes(x=beg*divby,y=end*divby,size=count/sum(count)),color='blue4') +
#    geom_point(data=dmmid.not100[dmmid.not100$perc >= 2,],aes(x=beg*divby,y=end*divby,size=count/sum(count)),color='green4') +
ggplot(dmall.not100[dmall.not100$perc >= 1, ],aes((1+beg)*divby,(1+end)*divby)) +
   geom_point(data=dm1,aes(beg,end),pch=".",color='grey',alpha=0.1) + #aes(size=after_stat(density)),n=20,contour=FALSE)
   geom_point(aes(size=count/sum(count)),color='red4') +
   geom_point(aes(x=beg2*divby,y=end2*divby,size=count/sum(count)),color='blue4') +
   geom_point(aes(x=beg*divby,y=end2*divby,size=count/sum(count)),color='green4') +
   geom_segment(aes(xend=beg*divby,yend=end2*divby),color='red4') +
   geom_segment(aes(x=beg*divby,y=end2*divby,xend=beg2*divby,yend=end2*divby),color='blue4') +
#   geom_segment(aes(x=beg*divby,y=end2*divby,size=count/sum(count)),color='green4') +
   scale_size_continuous(range=c(0.2,1)) +
   coord_cartesian(xlim=c(0,120*divby),ylim=c(0,120*divby)) + 
   annotate(geom = 'segment',x=0,y=0,xend=120*divby,yend=120*divby)
   #geom_density_2d_filled(data=dm1,aes(beg,end),alpha=0.5) +
#x=beg,y=end,),pch=".")

VRs = unique(dm1$VR)[order(unique(dm1$VR))]
for (VRwant in VRs) {
VRwant = 10
treatwant='C'

dm1clust = dmall.not100[dmall.not100$perc >= 1.7,]
if (VRwant != 'any') {
   dm1temp = dm1[dm1$VR == VRwant,]
} else {
   dm1temp = dm1
}
if (treatwant != 'any') {
   dm1temp = dm1temp[dm1temp$treat == treatwant,]
} else {
   dm1temp = dm1temp
}


#df$cluster = 0
# dm1clust$cluster = seq(1,dim(dm1clust)[1])
# dm1clustuniq = unique(dm1clust$cluster)
# dm1clustuniq = dm1clustuniq[order(dm1clustuniq)]
# for (i in dm1clustuniq) {
#    curr = dm1clust[dm1clust$cluster == i,]
#    print(curr)
# #   dm1tmep[dm1tmep$beg >= curr$x0end & dm1tmep$beg <= curr$x1beg & dm1tmep$end >= curr$y0end & dm1tmep$end <= curr$y1beg,]$cluster = i
#    dm1temp[dm1temp$beg >= curr$beg*divby & dm1temp$beg <= curr$beg2*divby & dm1temp$end >= curr$end*divby & dm1temp$end <= curr$end2*divby,]$cluster = i
# }
# dm1temp = dm1temp[order(dm1temp$cluster,dm1temp$beg/50,dm1temp$beg,dm1temp$end/50,dm1temp$end),]
# dm1temp$y = seq(1,dim(dm1temp)[1])

ratio = min(1,as.integer(100*1200*0.2/dim(dm1temp)[1])/100)
if(ratio <= 0.1) {ratio = 0.1}
print(paste(VRwant,dim(dm1temp)[1],ratio))

p2 = ggplot(dm1clust,aes((0+beg)*divby,(1+end)*divby)) +
   annotate(geom = 'segment',x=0,y=0,xend=120*divby,yend=120*divby) +
   annotate(geom='rect',xmin=841,ymin=0*divby,xmax=1317,ymax=120*divby,fill=rgb(0,0,1,0.05),lwd=0) +
   annotate(geom='rect',ymin=841,xmin=0*divby,ymax=1317,xmax=120*divby,fill=rgb(0,0,1,0.05),lwd=0) +
   annotate(geom='rect',xmin=641,ymin=0*divby,xmax=841,ymax=120*divby,fill=rgb(0,1,0,0.1),lwd=0) +
   annotate(geom='rect',ymin=641,xmin=0*divby,ymax=841,xmax=120*divby,fill=rgb(0,1,0,0.1),lwd=0) +
   annotate(geom='rect',ymin=587,xmin=0*divby,ymax=605,xmax=120*divby,fill=rgb(1,0,0,0.1),lwd=0) +
   annotate(geom='rect',xmin=587,ymin=0*divby,xmax=605,ymax=120*divby,fill=rgb(1,0,0,0.1),lwd=0) +
   annotate(geom='text',x=(841+1317)/2,y=120*divby,label='SRNPN') +
   annotate(geom='text',y=(841+1317)/2,x=0,label='SRNPN') +
   annotate(geom='text',x=(641+841)/2,y=120*divby,label='VR') +
   annotate(geom='text',y=(641+841)/2,x=0,label='VR') +
   annotate(geom='text',x=(587+605)/2,y=120*divby,label='T7prom',size=2) +
   annotate(geom='text',y=(587+605)/2,x=0,label='T7prom',size=2) +
   scale_size_continuous(range=c(0.2,1)) +
   coord_cartesian(xlim=c(0,120*divby),ylim=c(0,120*divby)) + 
   geom_point(data=dm1temp,aes(beg,end),pch=".",color='grey',alpha=ratio) + #aes(size=after_stat(density)),n=20,contour=FALSE)
   theme_bw() + theme(legend.position = "none",panel.grid = element_blank()) + 
   ylab('R-loop termination zone (bp)') + xlab('R-loop initiation zone (bp)') +
   ggtitle(paste("T7 init VR",VRwant,'\ntreatment:',treatwant))

p2a = p2 + coord_flip()

p3 = p2 + 
   geom_point(aes(size=count/sum(count)),color='red4') +
   geom_point(aes(x=(0+beg2)*divby,y=(1+end2)*divby,size=count/sum(count)),color='blue4') +
   geom_point(aes(x=(0+beg)*divby,y=(1+end2)*divby,size=count/sum(count)),color='green4') +
   geom_segment(aes(xend=(0+beg)*divby,yend=(1+end2)*divby),color='red4') +
   geom_segment(aes(x=(0+beg)*divby,y=(1+end2)*divby,xend=(0+beg2)*divby,yend=(1+end2)*divby),color='blue4')
   #   geom_segment(aes(x=beg*divby,y=end2*divby,size=count/sum(count)),color='green4') +
   #ggtitle(paste("T7 init VR All treatment All"))
#   pdf(paste("final_T7_init_VR_",VRwant,'treatment_',treatwant,'.pdf',sep=""),width=10,height=5)


dm1temp$cluster = 0
dm1clust$cluster = seq(1,dim(dm1clust)[1])
colnames(dm1clust) = c("x0end","y0end",'x1beg','y1beg','count','perc','cluster')
dm1clustuniq = unique(dm1clust$cluster)
dm1clustuniq = dm1clustuniq[order(dm1clustuniq)]

for (i in dm1clustuniq) {
   curr = dm1clust[dm1clust$cluster == i,]
#   print(curr)
   dm1temp[dm1temp$beg >= curr$x0end*divby & dm1temp$beg <= curr$x1beg*divby & dm1temp$end >= curr$y0end*divby & dm1temp$end <= curr$y1beg*divby,]$cluster = i
}
dm1temp = dm1temp[order(dm1temp$cluster,dm1temp$beg/divby,dm1temp$beg,dm1temp$end/divby,dm1temp$end),]
dm1temp$y = seq(1,dim(dm1temp)[1])

p4 = ggplot(dm1temp,aes(beg,y)) + 
   geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster)))


dm1temp = dm1temp[order(dm1temp$cluster,dm1temp$end/divby,dm1temp$end,dm1temp$beg/divby,dm1temp$beg),]
dm1temp$y = seq(1,dim(dm1temp)[1])

p5 = ggplot(dm1temp,aes(beg,y)) + 
   geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster)))


grid.arrange(p2,p3,p4,p5,ncol=2)
   #   dev.off()
}

#geom_density_2d_filled(data=dm1,aes(beg,end),alpha=0.5) +
#x=beg,y=end,),pch=".")
