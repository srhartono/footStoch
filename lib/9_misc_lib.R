ps = function(df,by,gp,p=NA,add=F,y.var = NA,myparams=list(),print=F,group='cluster',colorgroup=NA) {
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
  
  if (is.na(colorgroup)) {
    colorgroup = group
  }
  
  mycolors = c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Set2"))
  mybreaks = seq(0,length(mycolors)-1)
  #p = 
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
  if (defined(gp)) {
    if (defined(gp$mytitle)) {
      if (grepl('T7',gp$mytitle)) {
        p = p + 
          annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=max(temp[,y.var]),fill=rgb(0.25,0.25,0.25,0.1),color=rgb(0,0,0,0.10)) +
          annotate(geom='text',x=595,y=max(temp[,y.var]),label='T7') +
          annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=max(temp[,y.var]),fill=rgb(0,0.50,0,0.1),color=rgb(0,0,0,0.10)) +
          annotate(geom='text',x=741,y=max(temp[,y.var]),label='VR') +
          annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=max(temp[,y.var]),fill=rgb(0,0,0.50,0.1),color=rgb(0,0,0,0.10)) +
          annotate(geom='text',x=1000,y=max(temp[,y.var]),label='SNRPN') + theme_bw()   
      }
    }
  }
  if (add == T) {
    p = p + geom_segment(data=temp,aes(x=bybeg,xend=byend,y=temp[,y.var],yend=temp[,y.var],group=group,color=af(group)))
  }
  if(print == T)
  {
    print(p)
  }
  return(p)
}

ps2 = function(df,by,gp,p=NA,add=F,y.var = NA,myparams=list(),print=F,group='cluster',colorgroup=NA) {
  temp = df
  if (is.na(colorgroup)) {colorgroup = group}
  bybeg = 'beg'
  byend = 'end'
  if (grepl('mean',by)) {
    bybeg = 'meanbeg'
    byend = 'meanend'
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
  colorgroup='colorgroup'
  mylist = get_color_breaks()
  mycolors = mylist$mycolorgroup
  mybreaks = mylist$mybreaks
  # tempz2 = get_y(tempz2[order(tempz2$colorgroup,tempz2$beg,tempz2$end),])
  # ggplot(tempz2,aes(beg,y,group=colorgroup)) + 
  #   geom_segment(aes(x=beg,xend=end,y=y,yend=y,group=allcluster.name,color=af(colorgroup))) +
  #   scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks)
  # mycolors = c(brewer.pal(9,"Set1"),brewer.pal(9,"Set2"),brewer.pal(9,"Set2"))
  # mybreaks = seq(0,length(mycolors)-1)
  #p = 
  mygroup = af(temp[,group])
  p = ggplot(temp,aes(x=temp[,bybeg],y=temp[,y.var])) +
    coord_cartesian(xlim=c(gp$minX,gp$maxX),ylim=c(gp$minY,max(temp[,y.var]))) +#max(500,max(temp$y)))) +
    geom_segment(aes(x=temp[,bybeg],xend=temp[,byend],y=temp[,y.var],yend=temp[,y.var],color=mygroup)) +
    #,group=temp[,colorgroup],color=af(temp[,colorgroup]))) +#,color=rgb(0,0,0,0)) +
    scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) +
    #  geom_rect(aes(xmin=bybeg,xmax=byend,ymin=y-0.5,ymax=y+0.5,fill=cluster),color=rgb(0,0,0,0)) +
    #  scale_x_continuous(breaks=seq(0,gp$maxX,500)) +
    xlab('R-loop Position (bp)') + ylab('Read Number') +
    ggtitle(gp$mytitle.wrap) +
    theme_bw() +
    theme(plot.title = element_text(size=5), legend.position='none')
  if (defined(gp)) {
    if (defined(gp$mytitle)) {
      if (grepl('T7',gp$mytitle)) {
        p = p + 
          annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=max(temp[,y.var]),fill=rgb(0.25,0.25,0.25,0.1),color=rgb(0,0,0,0.10)) +
          annotate(geom='text',x=595,y=max(temp[,y.var]),label='T7') +
          annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=max(temp[,y.var]),fill=rgb(0,0.50,0,0.1),color=rgb(0,0,0,0.10)) +
          annotate(geom='text',x=741,y=max(temp[,y.var]),label='VR') +
          annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=max(temp[,y.var]),fill=rgb(0,0,0.50,0.1),color=rgb(0,0,0,0.10)) +
          annotate(geom='text',x=1000,y=max(temp[,y.var]),label='SNRPN') + theme_bw()   
      }
    }
  }
  if (add == T) {
    p = p + geom_segment(data=temp,aes(x=bybeg,xend=byend,y=temp[,y.var],yend=temp[,y.var],group=group,color=af(group)))
  }
  if(print == T)
  {
    print(p)
  }
  return(p)
}
get_color_breaks = function() {
   
  colzbin = c('Greens','Reds','YlOrBr','Purples','Blues','RdPu')
  colzs = rep(colzbin,10)
  colzs2 = data.frame()
  for (i in 1:size(colzs)) {
    colz = colzs[i]
    colzdf = brewer.pal(9,colz)[3:9] 
    if (i == 1) {
      colzs2 = c(colzdf) #as.data.frame(t(colzdf))
    }
    else {
      colzs2 = c(colzs2,colzdf) #rbind(colzs2,t(colzdf))
    }
  }
  mybreaks = c(0,seq(1,size(colzs2)))#,-2)
  colzs2 = c(rgb(0.25,0.25,0.25),colzs2)#,rgb(0.5,0.5,0.5))
  mycolors = colzs2
  return(list(mycolorgroup=mycolors,mybreaks=mybreaks))
}
graph.clust = function(df,dfclust,pdfout,gp,group) {
  tempz = df
  tempzclust = dfclust
  tempz$cluster = tempz$allcluster.name
  get_y = function(x) {
    x$y = seq(1,size(x))
    return(x)
  }
  tempz = get_y(tempz[order(tempz$allcluster.name,tempz$beg,tempz$end),])
  mylist = get_color_breaks()
  mycolors = mylist$mycolorgroup
  mybreaks = mylist$mybreaks

  pm2beg = ps2(df = tempz,by = 'beg',gp = gp,group = 'colorgroup')
  # df = tempz
  # by = 'beg'
  # group = 'allcluster.name'
  # colorgruop = 'colorgroup'
  # pm2beg
  
  tempz = get_y(tempz[order(tempz$cluster,tempz$end,tempz$beg),])
  pm2end = ps2(df = tempz,by = 'end',gp = gp,group = 'colorgroup')
  
  pm1 = ggplot(tempz,aes(group=af(allcluster.name),beg,end)) +
    geom_point(pch='.',col='grey') +
    annotate(geom='segment',x = 0,xend=2500,y=0,yend=2500) +
    scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) +
    geom_rect(data=tempzclust,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,color=af(colorgroup)),fill=NA) +
    geom_text(data=tempzclust,aes(x=x0end/2+x1beg/2,y=y0end/2+y1beg/2,label=colorgroup)) +
    theme_bw()
  
  #paste(allcluster.name,'(',allcluster,')',sep='')),fill=NA)
  
  pm1beg = ggplot(tempz,aes(beg)) +#,end)) +
    stat_density(bw = 5,color='black',lwd=0.5,fill=NA) +
    geom_rect(data=tempzclust,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=0,ymax=0.006,fill=begcluster),color=NA,alpha=0.1) + 
    scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) +
    theme_bw() + theme(panel.grid = element_blank())

  pm1end = ggplot(tempz,aes(end)) +#,end)) +
    stat_density(bw = 5,color='black',lwd=0.5,fill=NA) +
    geom_rect(data=tempzclust,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=0,ymax=0.006,fill=endcluster),color=NA,alpha=0.1) + 
    scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) +
    theme_bw() + theme(panel.grid = element_blank())
  
  pdf(pdfout,width=20,height=20)
  print(grid.arrange(pm1,arrangeGrob(pm1beg,pm1end,nrow=2),pm2beg,pm2end,nrow=2,ncol=2))
  dev.off()
  #geom_point(pch='.',col='grey') +
  #annotate(geom='segment',x = 0,xend=2500,y=0,yend=2500) +
  #geom_rect(data=allcluster,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,color=allcluster),fill=NA)
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