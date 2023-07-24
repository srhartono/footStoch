library("fitdistrplus")
library(ggplot2)
library(reshape2)
library(MASS)
library(swfscMisc)
library(RcolorBrewer)
#runif(n,min=0,max=1) = gives uniform dist
#rnorm(n,mean=0,sd=1) = gives random dist
#rexp(n,rate=1) = gives exp dis
#rpois(n,lambda=1)
#rbinom(n,trialsize=trialsize,probsuccess=probsucess)
#rgamma(n,shape=shape,rate=1,scale=1/rate)
#rbeta(n,shape1=shape1,shape2=shape2,ncp=0)
#rweibull(n,shape=shape,scale=1)

params2 = list(
  gene     = 'T7_init_VR_20',
  treat    = 'C_T7Mix_Tx',
  peaktype = 'any',
  thres    = 0,
  VR       = 20
)
params_not_exact = c(F,T,F,F,F)

VRwant = params$VR
treatwant=params$treat
final0$VR = NA
final0[grepl("VR_[0-9]+",final0$gene,perl=T),]$VR = as.numeric(as.character(gsub("^.+VR_([0-9]+).*$","\\1",final0[grepl("VR_[0-9]+",final0$gene,perl=T),]$gene,perl=T)))
final0$thres = 0
dmall.not100 = slice_df(final0,params=params2,params_not_exact = params_not_exact)
dm1$gene = dm1$chr
dm1$thres = dm1$t0
dm1temp = slice_df(dm1,params=params,params_not_exact = params_not_exact)
dm1temp$y = seq(1,dim(dm1temp)[1])
dmall.not100$count = 1
dmall.not100$perc = as.integer(dmall.not100$count/sum(dmall.not100$count)*1000)/10
#dmall.not100[dmall.not100$beg == dmall.not100$end,]$end = dmall.not100[dmall.not100$beg == dmall.not100$end,]$end + 1
#dmall.not100[dmall.not100$beg2 == dmall.not100$end2,]$beg2 = dmall.not100[dmall.not100$beg2 == dmall.not100$end2,]$beg2 - 1


dm1clust = dmall.not100 #[dmall.not100$perc >= 1.7,]
dm1clust$cluster = seq(1,dim(dm1clust)[1])
dm1clustuniq = unique(dm1clust$cluster)
dm1clustuniq = dm1clustuniq[order(dm1clustuniq)]

# __END__
# divby = gp$divby
# # for (i in dm1clustuniq) {
# #   curr = dm1clust[dm1clust$cluster == i,]
# #   if (defined(curr)) {
# #     print(curr)
# # #    if (defined(dm1temp[dm1temp$beg >= curr$beg*divby & dm1temp$beg <= curr$beg2*divby & dm1temp$end >= curr$end*divby & dm1temp$end <= curr$end2*divby,])) {
# #     if (dim(dm1temp[dm1temp$beg >= curr$x0end*divby & dm1temp$beg <= curr$x1beg*divby & dm1temp$end >= curr$y0end*divby & dm1temp$end <= curr$y1beg*divby,])[1] > 0) {
# #       #   dm1tmep[dm1tmep$beg >= curr$x0end & dm1tmep$beg <= curr$x1beg & dm1tmep$end >= curr$y0end & dm1tmep$end <= curr$y1beg,]$cluster = i
# # #      dm1temp[dm1temp$beg >= curr$beg*divby & dm1temp$beg <= curr$beg2*divby & dm1temp$end >= curr$end*divby & dm1temp$end <= curr$end2*divby,]$cluster = i
# #         dm1temp[dm1temp$beg >= curr$x0end*divby & dm1temp$beg <= curr$x1beg*divby & dm1temp$end >= curr$y0end*divby & dm1temp$end <= curr$y1beg*divby,]$cluster = i
# #     }
# #   }
# # }
# # dm1temp = dm1temp[order(dm1temp$cluster,dm1temp$beg/50,dm1temp$beg,dm1temp$end/50,dm1temp$end),]
# # dm1temp$y = seq(1,dim(dm1temp)[1])
# 
# df.testplotstoch = subset(dm1temp,select=c('beg','end','cluster'))
# myclust = dm1clust
# # myclust$x0end = myclust$beg*divby
# # myclust$y0end = myclust$end*divby
# # myclust$x1beg = myclust$beg2*divby
# # myclust$y1beg = myclust$end2*divby
# # myclust$x1end = myclust$x0end
# # myclust$y1end = myclust$y1beg
# # myclust$x0beg = myclust$x0end
# # myclust$y0beg = myclust$y1beg
myclust = subset(myclust,select=c('cluster','x0end','y0end','x1end','y1end','x0beg','y0beg','x1beg','y1beg'))

clustmin = data.frame(cluster=0,x0end=0,y0end=0,
                      x0beg=0,y0beg=min(myclust$x0end)-1,
                      x1end=0,y1end=min(myclust$x0end)-1,
                      x1beg=min(myclust$x0end)-1,y1beg=3000)
clustmax = data.frame(cluster=max(myclust$cluster)+1,
                      x0end=min(myclust$x0end)-1,y0end=max(myclust$y1beg)+1,
                      x0beg=max(myclust$x0end)+1,y0beg=3000,
                      x1end=max(myclust$x0end)+1,y1end=3000,
                      x1beg=3000,y1beg=3000)
myclust=rbind(clustmin,clustmax,myclust)
# 
# 
# df.testplotstoch$cluster = 0
# p.xy = ggplot(df.testplotstoch,aes(beg,end)) +
#   geom_point(pch=".",color='grey') +
#   annotate(geom='segment',x=0,y=0,xend=3000,yend=3000) +
#   theme_bw() +
#   theme(legend.position = 'none') +
#   coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
#   geom_rect(data=myclust,aes(x=x0end,y=y0end,xmin=x0end,ymin=y0end,xmax=x1beg,ymax=y1beg,fill=as.factor(cluster)),color=rgb(0,0,0,0.5),alpha=0.25) +
#   geom_text(data=myclust,aes(x=(x1beg+x0end)/2,y=(y1beg+y0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5) 
# p.xy.clust = ggplot(df.testplotstoch,aes(end,beg)) +
#   geom_point(pch=".",color='grey') +
#   annotate(geom='segment',x=0,y=0,xend=3000,yend=3000) +
#   theme_bw() +
#   theme(legend.position = 'none') +
#   geom_rect(data=myclust,aes(x=y0end,y=x0end,xmin=y0end,ymin=x0end,xmax=y1beg,ymax=x1beg,fill=as.factor(cluster)),color=rgb(0,0,0,0.5),alpha=0.25) +
#   geom_text(data=myclust,aes(x=(y1beg+y0end)/2,y=(x1beg+x0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5) 
# # geom_point(data=myclust,aes(x0end,y0end),color='red2') +
# # geom_point(data=myclust,aes(x1end,y1end),color='red4') +
# # geom_point(data=myclust,aes(x1beg,y1beg),color='blue4') + 
# # geom_segment(data=myclust,aes(x=x0end,y=y0end,xend=x1end,yend=y1end),color='red4') +
# # geom_segment(data=myclust,aes(x=x1end,y=y1end,xend=x1beg,yend=y1beg),color='blue4') +
# # geom_point(data=myclust,aes(x1end,y1end),color='red4') +
# # geom_point(data=myclust,aes(x1beg,y1beg),color='blue4') + 
# 
# 
# df.testplotstoch$cluster = 0
# myclustuniq = unique(myclust$cluster)
# myclustuniq = myclustuniq[order(myclustuniq)]
# df.testplotstoch$RIZ = 0
# df.testplotstoch$RTZ = 0
# for (i in myclustuniq) {
#   curr = myclust[myclust$cluster == i,]
#   print(curr)
#   df.testplotstoch[df.testplotstoch$beg >= curr$x0end & df.testplotstoch$beg <= curr$x1beg & df.testplotstoch$end >= curr$y0end & df.testplotstoch$end <= curr$y1beg,]$cluster = i
# }
# for (i in myclustuniq) {
#   curr = myclust[myclust$cluster == i,]
#   print(curr)
#   if (curr$y0end - curr$x0end <= 5000) {
#     if (dim(df.testplotstoch[df.testplotstoch$cluster == i & df.testplotstoch$beg >= curr$x0end & df.testplotstoch$beg <= curr$x0end+divby,])[1] > 0) {
#       df.testplotstoch[df.testplotstoch$cluster == i & df.testplotstoch$beg >= curr$x0end & df.testplotstoch$beg <= curr$x0end+divby,]$RIZ = 1
#     }
#   }
#   if (curr$y0end - curr$x0end <= 5000) {
#     if (dim(df.testplotstoch[df.testplotstoch$cluster == i & df.testplotstoch$beg >= curr$y0end-divby & df.testplotstoch$beg <= curr$y0end,])[1] > 0) {
#       df.testplotstoch[df.testplotstoch$cluster == i & df.testplotstoch$beg >= curr$y0end-divby & df.testplotstoch$beg <= curr$y0end,]$RIZ = 1
#     }
#   }
#   if (curr$y1beg - curr$x1beg <= 5000) {
#     if (dim(df.testplotstoch[df.testplotstoch$cluster == i & df.testplotstoch$end >= curr$x1beg & df.testplotstoch$end <= curr$x1beg+divby,])[1] > 0) {
#       df.testplotstoch[df.testplotstoch$cluster == i & df.testplotstoch$end >= curr$x1beg & df.testplotstoch$end <= curr$x1beg+divby,]$RTZ = 1
#     }
#   }
#   if (curr$y1beg - curr$x1beg <= 5000) {
#     if (dim(df.testplotstoch[df.testplotstoch$cluster == i & df.testplotstoch$end >= curr$y1beg-divby & df.testplotstoch$end <= curr$y1beg,])[1] > 0) {
#       df.testplotstoch[df.testplotstoch$cluster == i & df.testplotstoch$end >= curr$y1beg-divby & df.testplotstoch$end <= curr$y1beg,]$RTZ = 1
#     }
#   }
# }
# 
# df.testplotstoch = df.testplotstoch[order(df.testplotstoch$cluster,df.testplotstoch$beg/50,df.testplotstoch$beg,df.testplotstoch$end/50,df.testplotstoch$end),]
# df.testplotstoch$y = seq(1,dim(df.testplotstoch)[1])
# window=100
# begtest = test_unif(df.testplotstoch,'beg',divby=divby,window=window)
# begtest$goodup$y = begtest$goodup$index
# goodup=c()
# for (i in seq(1,length(begtest$goodup$y))) {
#   goodup = c(goody,seq(begtest$goodup$y[i],begtest$goodup$y[i]+window))
# }
# goodup = unique(goodup)
# begtestup = df.testplotstoch[df.testplotstoch$y %in% goodup,]
# p.orderbybeg = ggplot(df.testplotstoch,aes(beg,y)) +
#   geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster))) +
#   coord_cartesian(xlim=c(0,3000)) + theme_bw() + theme(legend.position = 'top')
# p.orderbybeg2 = p.orderbybeg + 
#   geom_rect(data=df.testplotstoch[df.testplotstoch$RIZ == 1,],aes(x=0,y=y,xmin=0,xmax=50,ymin=y,ymax=y+1),fill='red4') +
#   geom_rect(data=begtest$goodup,aes(x=0,y=index,xmin=50,xmax=100,ymin=index,ymax=index+window),fill='orange')
# #   geom_rect(data=myclust,aes())
# #p.orderbybeg2 = p.orderbybeg + geom_rect(data=begtest$goodup,aes(x=0,y=index,xmin=0,xmax=3000,ymin=index,ymax=index+1),fill='red4',alpha=0.5)
# #if (dim(begtest$goodsp)[1] > 0) {
# #   p.orderbybeg2 = p.orderbybeg2 + geom_rect(data=begtest$goodsp,aes(x=0,y=index,xmin=2550,xmax=2600,ymin=index,ymax=index+window-50),fill='blue4')
# #}
# #p.orderbybeg2 = p.orderbybeg2 + geom_rect(data=begtest$goodlm,aes(x=0,y=index,xmin=2500,xmax=2600,ymin=index,ymax=index+window-50),fill='green4')
# #p.orderbybeg2 = p.orderbybeg2 + geom_segment(data=begtest$goodup,aes(x=,xend=beg,y=y,yend=y,group=as.factor(cluster)),color='black',shape='.')
# #p.orderbybeg2 = p.orderbybeg + geom_rect(data=begtestup,aes(x=0,y=y,xmin=2600,xmax=2700,ymin=y,ymax=y+1),fill='red4')
# # p.orderbybeg2 = p.orderbybeg + geom_rect(data=begtest$goodup,aes(x=0,y=index,xmin=0,xmax=200,ymin=index,ymax=index+window-50),fill='red4')
# 
# df.testplotstoch = df.testplotstoch[order(df.testplotstoch$cluster,df.testplotstoch$end/50,df.testplotstoch$end,df.testplotstoch$beg/50,df.testplotstoch$beg),]
# df.testplotstoch$y = seq(1,dim(df.testplotstoch)[1])
# endtest = test_unif(df.testplotstoch,'end',divby=divby,window=window)
# endtest$goodup$y = endtest$goodup$index
# goody=c()
# for (i in seq(1,length(endtest$goodup$y))) {
#   goody = c(goody,seq(endtest$goodup$y[i],endtest$goodup$y[i]+window))
# }
# endtest2 = df.testplotstoch[df.testplotstoch$y %in% endtest$goodup$y,]
# p.orderbyend = ggplot(df.testplotstoch,aes(beg,y)) +
#   geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster))) +
#   coord_cartesian(xlim=c(0,3000)) + theme_bw() + theme(legend.position = 'top')
# #p.orderbyend2 = p.orderbyend + geom_rect(data=endtest$goodup,aes(x=0,y=index,xmin=2600,xmax=2700,ymin=index,ymax=index+window),fill='red4')
# p.orderbyend2 = p.orderbyend + 
#   geom_rect(data=df.testplotstoch[df.testplotstoch$RTZ == 1,],aes(x=0,y=y,xmin=2000,xmax=2050,ymin=y,ymax=y+1),fill='red4') +
#   geom_rect(data=endtest$goodup,aes(x=0,y=index,xmin=2050,xmax=2100,ymin=index,ymax=index+window),fill='orange')
# #p.orderbyend2 = p.orderbyend
# #if (dim(endtest$goodps)[1] > 0) {
# #   p.orderbyend2 = p.orderbyend2 + geom_rect(data=endtest$goodlm,aes(x=0,y=index,xmin=2550,xmax=2600,ymin=index,ymax=index+window-50),fill='blue4')
# #}
# #p.orderbyend2 = p.orderbyend2 + geom_rect(data=endtest$goodlm,aes(x=0,y=index,xmin=2500,xmax=2600,ymin=index,ymax=index+window-50),fill='green4')
# #p.orderbyend2 = p.orderbyend2 + geom_point(data=endtest$goodup,aes(x=end,y=y,group=as.factor(cluster)),color='black',shape='.')
# 
# dm2 = rep(0,3000)
# for (i in 1:dim(df.testplotstoch)[1]) {
#   dm2[seq(df.testplotstoch$beg[i],df.testplotstoch$end[i])] = dm2[seq(df.testplotstoch$beg[i],df.testplotstoch$end[i])] + 1
# }
# dm2 = data.frame(beg=seq(1,3000),y=dm2)
# pdens = ggplot(dm2,aes(beg,y)) + geom_line() + theme_bw()
# 
# pdf('tri3_T7_init_VR_20_T7INITVR20APALITRXNDURINGSSB.pdf')
# grid.arrange(p.xy,p.xy.clust,p.orderbybeg2,p.orderbyend2,ncol=2,nrow=2)
# dev.off()
# #head(df.testplotstoch$y)
# #test = ((df.testplotstoch[df.testplotstoch$cluster == 2 & df.testplotstoch$beg >= 700 & df.testplotstoch$beg <= 800,]$beg)); plot(density(test)); uniform.test(hist(test))
# 
# 
