ps = function(df,by,gp,p=NA,add=F,y.var = NA,params=list(),print=F,group='cluster') {
  temp = df
  bybeg = temp$beg
  byend = temp$end
  if (grepl('mean',by))
  {
    bybeg = temp$meanbeg
    byend = temp$meanend
  }
  if (defined(y.var) == FALSE) {
    y.var = 'y'
  }
  p = ggplot(temp,aes(x=temp[,by],y=temp[,y.var],group=temp[,group])) +
  coord_cartesian(xlim=c(gp$minX,gp$maxX),ylim=c(gp$minY,max(temp[,y.var]))) +#max(500,max(temp$y)))) +
  geom_segment(aes(x=bybeg,xend=byend,y=temp[,y.var],yend=temp[,y.var],group=temp[,group],color=af(temp[,group]))) +#,color=rgb(0,0,0,0)) +
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
