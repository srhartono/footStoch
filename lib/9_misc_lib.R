ps = function(df,by,gp,p=NA,add=F,y.var = NA,myparams=list(),print=F,group='cluster') {
  temp = df
  bybeg = temp$beg
  byend = temp$end
  if (grepl('mean',by))
  {
    bybeg = temp$meanbeg
    byend = temp$meanend
  }
  # } else if (defined(temp[,by])) {
  #   if (grepl('beg',by)) {
  #     bybeg = by
  #     byend = gsub('beg','end',bybeg,perl=T)
  #   } else {
  #     byend = by
  #     bybeg = gsub('beg','end',bybeg,perl=T)
  #   }
  # }
  if (defined(y.var) == FALSE) {
    y.var = 'y'
  }
  
  mycolors = c(brewer.pal(9,"Set1"),brewer.pal(9,"Set3"))
  mybreaks = seq(0,length(mycolors)-1)
  p = ggplot(temp,aes(x=temp[,bybeg],y=temp[,y.var],group=temp[,group])) +
  coord_cartesian(xlim=c(gp$minX,gp$maxX),ylim=c(gp$minY,max(temp[,y.var]))) +#max(500,max(temp$y)))) +
  geom_segment(aes(x=bybeg,xend=byend,y=temp[,y.var],yend=temp[,y.var],group=temp[,group],color=af(temp[,group]))) +#,color=rgb(0,0,0,0)) +
  scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) +
#  geom_rect(aes(xmin=bybeg,xmax=byend,ymin=y-0.5,ymax=y+0.5,fill=cluster),color=rgb(0,0,0,0)) +
#  scale_x_continuous(breaks=seq(0,gp$maxX,500)) +
  xlab('R-loop Position (bp)') + ylab('Read Number') +
  ggtitle(gp$mytitle.wrap) +
  theme_bw() +
  theme(plot.title = element_text(size=5), legend.position='none')
  
  if (add == T) {
    p = p + geom_segment(data=temp,aes(x=bybeg,xend=byend,y=temp[,y.var],yend=temp[,y.var],group=group,color=af(group)))
  }
  if(print == T)
  {
    print(p)
  }
  return(p)
}

ggplot.xy = function(df=data.frame(),dfclust = data.frame(),outpdf = NA) {
  df$cluster = 0
  p.xy.beg = ggplot(df,aes(beg,end)) +
    geom_point(pch=".",color='grey') +
    annotate(geom='segment',x=0,y=0,xend=3000,yend=3000) +
    theme_bw() +
    theme(legend.position = 'none') +
    coord_cartesian(xlim=c(0,3000),ylim=c(0,3000))
  p.xy.clustbeg = p.xy.beg +  geom_rect(data=dfclust,aes(x=x0end,y=y0end,xmin=x0end,ymin=y0end,xmax=x1beg,ymax=y1beg,fill=as.factor(cluster)),color=rgb(0,0,0,0.5),alpha=0.25) +
    geom_text(data=dfclust,aes(x=(x1beg+x0end)/2,y=(y1beg+y0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5)

  p.xy.end = ggplot(df,aes(end,beg)) +
    geom_point(pch=".",color='grey') +
    annotate(geom='segment',x=0,y=0,xend=3000,yend=3000) +
    theme_bw() +
    theme(legend.position = 'none') +
    coord_cartesian(xlim=c(0,3000),ylim=c(0,3000))
    
  p.xy.clustend = p.xy.end + geom_rect(data=dfclust,aes(x=y0end,y=x0end,xmin=y0end,ymin=x0end,xmax=y1beg,ymax=x1beg,fill=as.factor(cluster)),color=rgb(0,0,0,0.5),alpha=0.25) +
    geom_text(data=dfclust,aes(x=(y1beg+y0end)/2,y=(x1beg+x0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5)
  if (!is.na(outpdf)) {
    pdf(outpdf,width=15,height=15)
    print(grid.arrange(p.xy.beg,p.xy.end,p.xy.clustbeg,p.xy.clustend,nrow=2,ncol=2))
    dev.off()
  } else {
    print(grid.arrange(p.xy.beg,p.xy.end,p.xy.clustbeg,p.xy.clustend,nrow=2,ncol=2))
  }
}


ggplot.xy.C = function(df=data.frame(),dfclust = data.frame(),outpdf = NA) {
  df$cluster = 0
  p.xy.begC = ggplot(df,aes(begC,endC)) +
    geom_point(pch=".",color='grey') +
    annotate(geom='segment',x=0,y=0,xend=1000,yend=1000) +
    theme_bw() +
    theme(legend.position = 'none') +
    coord_cartesian(xlim=c(0,1000),ylim=c(0,1000))
  # p.xy.clustbeg = p.xy.beg +  geom_rect(data=dfclust,aes(x=x0end,y=y0end,xmin=x0end,ymin=y0end,xmax=x1beg,ymax=y1beg,fill=as.factor(cluster)),color=rgb(0,0,0,0.5),alpha=0.25) +
    # geom_text(data=dfclust,aes(x=(x1beg+x0end)/2,y=(y1beg+y0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5)

  p.xy.endC = ggplot(df,aes(endC,begC)) +
    geom_point(pch=".",color='grey') +
    annotate(geom='segment',x=0,y=0,xend=1000,yend=1000) +
    theme_bw() +
    theme(legend.position = 'none') +
    coord_cartesian(xlim=c(0,1000),ylim=c(0,1000))
    
  # p.xy.clustend = p.xy.end + geom_rect(data=dfclust,aes(x=y0end,y=x0end,xmin=y0end,ymin=x0end,xmax=y1beg,ymax=x1beg,fill=as.factor(cluster)),color=rgb(0,0,0,0.5),alpha=0.25) +
  #   geom_text(data=dfclust,aes(x=(y1beg+y0end)/2,y=(x1beg+x0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5)
  # if (!is.na(outpdf)) {
  #   pdf(outpdf,width=15,height=15)
  #   print(grid.arrange(p.xy.beg,p.xy.end,p.xy.clustbeg,p.xy.clustend,nrow=2,ncol=2))
  #   dev.off()
  # } else {
  #  print(grid.arrange(p.xy.beg,p.xy.end,p.xy.clustbeg,p.xy.clustend,nrow=2,ncol=2))
  # }
  print(grid.arrange(p.xy.begC,p.xy.endC,nrow=2,ncol=2))
}

dist.test_get.color = function(pval,gp,verbose=F,debug=F) {
  unif.bin.cols = brewer.pal(9,"Reds")
  unif.bin.brks = seq(gp$minpval,1,length.out=9); unif.bin.brks.diff = unif.bin.brks[2]-unif.bin.brks[1]
  unif.col = unif.bin.cols[unif.bin.brks < pval & unif.bin.brks+unif.bin.brks.diff > pval]
  if (defined(unif.col) == FALSE) {
    if(pval > 0.8) {
      unif.col = unif.bin.cols[length(unif.bin.cols)]
    } else {
      unif.col = rgb(0,0,0,0)
    }
  }
  unif.col=unif.bin.cols[9]
  if (verbose == T) {print(unif.col)}
  return(unif.col)
}


get_fapos = function(myfa=list()) {
  
  nuc=as.array(as.character(strsplit(x = as.character(myfa$amp.seq[[1]]),split = '',perl = T)[[1]]))
  test2 = data.frame(nuc=nuc,isC = rep(0,length(nuc)))
  test2[test2$nuc == 'C',]$isC = 1
  test2$index = myfa$amp.beg + seq(1,dim(test2)[1])-1
  head(test2)
  #   nuc isC
  # 1   G   0
  # 2   G   0
  # 3   G   0
  # 4   G   0
  # 5   T   0
  # 6   C   0
  test3 = test2[test2$isC == 1,]
  test3$Cindex = myfa$amp.beg + seq(1,dim(test3)[1])-1
  test3$diff = c(diff(test3$index),0)
  return(test3)
}
#   head(test3[test3$index >= 400 & test3$index <= 500,'index'],n=20)
#   
#   test5 = data.frame(nuc=nuc,isC = rep(0,length(nuc)))
#   test5[test5$nuc == 'G',]$isC = 1
#   test5$index = myfa$amp.beg + seq(1,dim(test5)[1])-1
#   head(test5)
#   #   nuc isC
#   # 1   G   0
#   # 2   G   0
#   # 3   G   0
#   # 4   G   0
#   # 5   T   0
#   # 6   C   0
#   test5 = test5[test5$isC == 1,]
#   test5$Cindex = myfa$amp.beg + seq(1,dim(test5)[1])-1
#   test5$diff = c(diff(test5$index),0)
#   head(test5[test5$index >= 400 & test5$index <= 500,'index'],n=20)
#   
#   test4 = df[df$beg >= 400 & df$beg <= 500,]
#   head(test4[order(test4$beg),]$beg,n=20)
#   head(test3[test3$index >= 400 & test3$index <= 500,'index'],n=20)
#   dim(df[df$beg %in% test3$index,])[1]/dim(df)[1]*100
#   dim(df[df$end+1 %in% test3$index,])[1]/dim(df)[1]*100 
# }