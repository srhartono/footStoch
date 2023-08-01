ps = function(df,by,gp,p=NA,add=F,y.var = NA,my.params=list(),print=F,my.group=NA,my.bed=data.frame(),my.title=NA) {
  suppressWarnings({
    if (is.na(my.title)) {my.title = ps.get_my.title.if.NA(df,gp)}
    if (!'C.transform' %in% names(gp)) {
      gp$C.transform = FALSE
    }
    if (!'feature.type' %in% colnames(my.bed)) {
      my.bed$feature.type = my.bed$feature
    }
    if (defined(my.bed[grep('Barcode',my.bed$feature,invert = T),])) {
      my.bed = my.bed[grep('Barcode',my.bed$feature,invert = T),]
    }
    if (defined(my.bed[grep('_diff',my.bed$feature,invert = T),])) {
      my.bed = my.bed[grep('_diff',my.bed$feature,invert = T),]
    }
    temp = df
    bybeg = 'beg'
    byend = 'end'
    if (grepl('mean',by))
    {
      bybeg = 'meanbeg'
      byend = 'meanend'
    }
    
    if (defined(y.var) == FALSE) {
      y.var = 'y'
    }
    xlim0 = gp$minX
    xlim1 = gp$maxX
    
    ylim0 = gp$minY
    ylim1 = max(temp[,y.var])
    
    my.colors = c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Set2"))
    my.breaks = seq(0,length(my.colors)-1)
    if (is.na(my.group)) {
      temp.my.group = af(1)
    } else {
      temp.my.group = af(temp[,my.group])
    }
    p = ggplot(temp,aes(x=temp[,bybeg],y=temp[,y.var],group=temp.my.group)) +
      coord_cartesian(xlim=c(xlim0,xlim1),ylim=c(ylim0,1.2*ylim1)) +#max(500,max(temp$y)))) +
      geom_segment(aes(x=temp[,bybeg],xend=temp[,byend],y=temp[,y.var],yend=temp[,y.var],
                       group=temp.my.group,color=af(temp.my.group))) +#,color=rgb(0,0,0,0)) +
      scale_color_manual(values=my.colors,breaks=my.breaks,label=my.breaks) +
      ggtitle(gp$my.title.wrap) +
      theme_bw() +
      theme(plot.title = element_text(size=10), axis.title.x = element_text(size=10), legend.position='none')
    if (gp$C.transform == TRUE) {
      p = p + xlab(pasta('R-loop Position by Cytosine index (bp)\nSorted by ',by)) + ylab('Read Index')
    } else {
      p = p + xlab(pasta('R-loop Position (bp)\nSorted by ',by)) + ylab('Read Index')
    }
    
    if (defined(gp)) {
      if (defined(gp$my.title)) {
        if (grepl('T7',gp$my.title)) {
          if (size(my.bed) == 0) {
            my.bed = rbind(my.bed,data.frame(beg=587,end=605,feature='T7_Promoter',feature.type='T7_Promoter'))
            my.bed = rbind(my.bed,data.frame(beg=641,end=841,feature='VR',feature.type='VR'))
            my.bed = rbind(my.bed,data.frame(beg=842,end=1317,feature='SNRPN',feature.type='SRNRPN'))
          }
          #          for (i in 1:length(my.bed)) {
          
          p = p + 
            geom_rect(data=my.bed,aes(x=beg,y=end,xmin=beg,xmax=end,ymin=0,ymax=ylim1,group=af(feature.type),fill=af(feature.type)),alpha=0.15,color=rgb(0,0,0,0.15)) +
            geom_text(data=my.bed,aes(x=(beg+end)/2,y=1.02*ylim1,label=feature,group=af(feature.type)),size=3,angle=90,hjust=0)
          #              annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=max(temp[,y.var]),fill=rgb(0.25,0.25,0.25,0.1),color=rgb(0,0,0,0.10)) +
          #              annotate(geom='text',x=595,y=max(temp[,y.var]),label='T7\nPromoter',size=3)
          #          }
          # p = p + 
          #   annotate(geom='rect',xmin=587,xmax=605,ymin=0,ymax=max(temp[,y.var]),fill=rgb(0.25,0.25,0.25,0.1),color=rgb(0,0,0,0.10)) +
          #   annotate(geom='text',x=595,y=max(temp[,y.var]),label='T7\nPromoter') +
          #   annotate(geom='rect',xmin=641,xmax=841,ymin=0,ymax=max(temp[,y.var]),fill=rgb(0,0.50,0,0.1),color=rgb(0,0,0,0.10)) +
          #   annotate(geom='text',x=741,y=max(temp[,y.var]),label='VR') +
          #   annotate(geom='rect',xmin=842,xmax=1317,ymin=0,ymax=max(temp[,y.var]),fill=rgb(0,0,0.50,0.1),color=rgb(0,0,0,0.10)) +
          #   annotate(geom='text',x=1000,y=max(temp[,y.var]),label='SNRPN') + theme_bw()   
          
        }
      }
    }
    if (add == T) {
      p = p + geom_segment(data=temp,aes(x=temp[,bybeg],xend=temp[,byend],y=ylim1,yend=ylim1,group=af(my.group),color=af(my.group)))
    }
    if(print == T)
    {
      print(p)
    }
  })
  return(p)
}

ps.get_my.title.if.NA = function(df=data.frame(),gp=list) {
  my.title = ''
  if ('my.title.wrap' %in% names(gp)) {
    my.title = gp$my.title.wrap
  } else if ('my.title' %in% names(gp)) {
    my.title = gp$my.title
  } else {
    if ('gene' %in% df) {
      my.title = pasta(my.title,'\n',unique(df$chr)[1])
    }
    if ('treat' %in% df) {
      my.title = pasta(my.title,'\n',unique(df$treat)[1])
    }
    if ('threschar' %in% df) {
      my.title = pasta(my.title,'\n',unique(df$threschar)[1])
    }
    if ('peaktype' %in% df) {
      my.title = pasta(my.title,'\n',unique(df$peaktype)[1])
    }
  }
  return(my.title)
}

ps.get_color.breaks = function() {
  
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

get_y = function(x) {
  x$y = seq(1,size(x))
  return(x)
}

graph.clust = function(df,dfclust,pdfout,gp,my.group,my.bed) {
  tempz = df
  tempzclust = dfclust
  beg.density = density(tempz$beg,bw=5)
  end.density = density(tempz$end,bw=5)
  tempz$cluster = tempz$allcluster.name
  mylist = ps.get_color.breaks()
  mycolors = mylist$mycolorgroup
  mybreaks = mylist$mybreaks
  tempz = get_y(tempz[order(tempz$allcluster.name,tempz$beg,tempz$end),])
  pm2beg = ps(df = tempz,by = 'beg',gp = gp,group = 'colorgroup',my.bed=my.bed)
  # df = tempz
  # by = 'beg'
  # group = 'allcluster.name'
  # colorgruop = 'colorgroup'
  # pm2beg
  
  tempz = get_y(tempz[order(tempz$cluster,tempz$end,tempz$beg),])
  pm2end = ps(df = tempz,by = 'end',gp = gp,group = 'colorgroup',my.bed=my.bed)
  
  pm1 = ggplot.xy(df = tempz,dfclust = tempzclust,gp = gp,print = FALSE)
  # #ggplot(tempz,aes(group=af(allcluster.name),beg,end)) +
  #   geom_point(pch='.',col='grey') +
  #   annotate(geom='segment',x = 0,xend=2500,y=0,yend=2500) +
  #   scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) +
  #   geom_rect(data=tempzclust,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,color=af(colorgroup)),fill=NA) +
  #   geom_text(data=tempzclust,aes(x=x0end/2+x1beg/2,y=y0end/2+y1beg/2,label=colorgroup)) +
  #   theme_bw()
  
  #paste(allcluster.name,'(',allcluster,')',sep='')),fill=NA)
  
  tempzclust.beg = aggregate(tempzclust$x0end,by=list(tempzclust$begcluster),min); colnames(tempzclust.beg) = c('begcluster','x0end')
  tempzclust.beg$x1beg = aggregate(tempzclust$x1beg,by=list(tempzclust$begcluster),max)$x
  
  pm1beg = ggplot(tempz,aes(beg)) +#,end)) +
    stat_density(bw = 5,color='black',lwd=0.5,fill=NA) +
    geom_rect(data=tempzclust.beg,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=0,ymax=max(beg.density$y),fill=begcluster),color=NA,alpha=0.2) + 
    scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) +
    theme_bw() + theme(panel.grid = element_blank())
  
  tempzclust.end = aggregate(tempzclust$y0end,by=list(tempzclust$endcluster),min); colnames(tempzclust.end) = c('endcluster','y0end')
  tempzclust.end$y1beg = aggregate(tempzclust$y1beg,by=list(tempzclust$endcluster),max)$x
  
  pm1end = ggplot(tempz,aes(end)) +#,end)) +
    stat_density(bw = 5,color='black',lwd=0.5,fill=NA) +
    geom_rect(data=tempzclust.end,aes(x=0,y=0,xmin=y0end,xmax=y1beg,ymin=0,ymax=max(end.density$y),fill=endcluster),color=NA,alpha=0.2) + 
    scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks) +
    theme_bw() + theme(panel.grid = element_blank())
  
  pdf(pdfout,width=20,height=20)
  # print(plot(NA,xlim=c(0,3000),ylim=c(0,3000)))
  grid.arrange(pm1,arrangeGrob(pm1beg,pm1end,nrow=2),pm2beg,pm2end,nrow=2,ncol=2)
  dev.off()
  print(pdfout)
  #geom_point(pch='.',col='grey') +
  #annotate(geom='segment',x = 0,xend=2500,y=0,yend=2500) +
  #geom_rect(data=allcluster,aes(x=0,y=0,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg,color=allcluster),fill=NA)
}

ps.ggplot.xy = function(df=data.frame(),dfclust = data.frame(),gp,my.by=NA,print=T,outpdf = NA,my.title=NA) {
  if (is.na(my.title)) {ps.get_my.title.if.NA(df,gp)}
  
  df$cluster = 0
  p.xy.beg = ggplot(df,aes(beg,end)) +
    geom_point(pch=".",color='grey') +
    annotate(geom='segment',x=0,y=0,xend=gp$maxX,yend=gp$maxX) +
    theme_bw() +
    theme(legend.position = 'none') +
    coord_cartesian(xlim=c(0,gp$maxX),ylim=c(0,gp$maxX)) + ggtitle(my.title)
  
  p.xy.clustbeg = p.xy.beg +  geom_rect(data=dfclust,aes(x=x0end,y=y0end,xmin=x0end,ymin=y0end,xmax=x1beg,ymax=y1beg,fill=as.factor(cluster)),color=rgb(0,0,0,0.5),alpha=0.25) +
    geom_text(data=dfclust,aes(x=(x1beg+x0end)/2,y=(y1beg+y0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5) + ggtitle(my.title)
  
  p.xy.end = ggplot(df,aes(end,beg)) +
    geom_point(pch=".",color='grey') +
    annotate(geom='segment',x=0,y=0,xend=gp$maxX,yend=gp$maxX) +
    theme_bw() +
    theme(legend.position = 'none') +
    coord_cartesian(xlim=c(0,gp$maxX),ylim=c(0,gp$maxX)) + ggtitle(my.title)
  
  p.xy.clustend = p.xy.end + geom_rect(data=dfclust,aes(x=y0end,y=x0end,xmin=y0end,ymin=x0end,xmax=y1beg,ymax=x1beg,fill=as.factor(cluster)),color=rgb(0,0,0,0.5),alpha=0.25) +
    geom_text(data=dfclust,aes(x=(y1beg+y0end)/2,y=(x1beg+x0end)/2,label=cluster,group=as.factor(cluster)),alpha=0.8,size=5) + ggtitle(my.title)
  
  if (!is.na(outpdf)) {
    pdf(outpdf,width=15,height=15)
    print(grid.arrange(p.xy.beg,p.xy.end,p.xy.clustbeg,p.xy.clustend,nrow=2,ncol=2))
    dev.off()
  } else if (print == T) {
    print(grid.arrange(p.xy.beg,p.xy.end,p.xy.clustbeg,p.xy.clustend,nrow=2,ncol=2))
  } else if (is.na(my.by)) {
    return(p.xy.clustbeg)
  } else if (my.by='beg') {
    return(p.xy.clustbeg)
  } else if (my.by='end') {
    return(p.xy.clustend)
  } else {
    return(p.xy.clustbeg)
  }
}

ps.ggplot.xy.C = function(df=data.frame(),dfclust = data.frame(),outpdf = NA) {
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

dist_test.get.color = function(pval,gp,verbose=F,debug=F) {
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

ps.get_color.group = function(df,dfclust,gp) {
  col.divby = 100
  if ('C.transform' %in% names(gp)) {
    if (gp$C.transform == TRUE) {
      col.divby = 25
    }
  }
  tempz = df
  tempzclust = dfclust
  #  tempz = dfm
  #  tempzclust = dfmclust
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
  my.break = c(0,seq(1,size(colzs2)))#,-2)
  tempzclust$begc = ai((tempzclust$x0end-gp$promoter$beg) / col.divby)
  tempzclust$endc = ai((tempzclust$y1beg-gp$promoter$beg) / col.divby)
  if (defined(tempzclust[tempzclust$endc < 0,])) {tempzclust[tempzclust$endc < 0,]$endc = 0}
  if (defined(tempzclust[tempzclust$endc > 21,])) {tempzclust[tempzclust$endc > 21,]$endc = 21}
  tempzclust$endc2 = ai(tempzclust$endc) + 1
  tempzclust$colorgroup = tempzclust$begc * 21 + tempzclust$endc2
  if (defined(tempzclust[tempzclust$colorgroup <= 0,])) {
    tempzclust[tempzclust$colorgroup <= 0,]$colorgroup = 0
  }
  tempz2 = merge(tempz,subset(tempzclust,select=c('allcluster.name','colorgroup')),by='allcluster.name')
  return(list(df=tempz2,dfclust=tempzclust))
}

ps.CpGprof = function(df,by,gp,my.bed=data.frame(),print=F,my.title=NA) {
  my.by = by
  if (is.na(my.title)) {my.title = ps.get_my.title.if.NA(df,gp)}
  print(my.title)
  temp = df
  xlim0 = gp$minX
  xlim1 = gp$maxX
  
  ylim0 = gp$minY
  ylim1 = 1.5
  p = ggplot(temp,aes(df[,my.by],y=index)) +
    geom_line(aes(df[,my.by],y=CGdens),color='blue2') +
    geom_line(aes(df[,my.by],y=GCcont),color='green4') +
    geom_line(aes(df[,my.by],y=GCskew+0.5),color='red2') +
    geom_line(aes(df[,my.by],y=Gcont),color='red2') +
    theme_bw() +
    theme(plot.title = element_text(size=10), axis.title.x = element_text(size=10), legend.position='none') +
    ggtitle(my.title) +
    coord_cartesian(ylim=c(ylim0,1.2*ylim1),xlim=c(xlim0,xlim1)) +
    annotate(geom='segment',x=gp$minX,xend=gp$maxX,y=0.5,yend=0.5,lty=2) +
    xlab('Position in Plasmid (bp)') + ylab('CG Dens & GC Cont') +
    annotate(geom='segment',x=xlim1*0.9,y=0,xend=xlim1*0.9,yend=ylim1,color='red4') +
    annotate(geom='text',x=xlim1,y=0.5,label='GC Skew',angle=-90,vjust=0,size=4,color='red4')
  for (i in seq(0,ylim1,0.1)) {
    i0 = i
    i1 = i + 0.1
    p = p + annotate(geom='segment',x=xlim1*0.9-xlim1*0.01,y=i,xend=xlim1*0.9+xlim1*0.01,yend=i,color='red4')
    p = p + annotate(geom='text',x=xlim1*0.9+xlim1*0.01*2,y=i,label=i-0.5,hjust=0,size=3,color='red4')
  }  
  
  
  if (defined(my.bed)) {
    if (!'feature.type' %in% colnames(my.bed)) {
      my.bed$feature.type = my.bed$feature
    }
    if (defined(my.bed[grep('Barcode',my.bed$feature,invert = T),])) {
      my.bed = my.bed[grep('Barcode',my.bed$feature,invert = T),]
    }
    if (defined(my.bed[grep('_diff',my.bed$feature,invert = T),])) {
      my.bed = my.bed[grep('_diff',my.bed$feature,invert = T),]
    }
    p = p + 
      geom_rect(data=my.bed,aes(x=beg,y=end,xmin=beg,xmax=end,ymin=0,ymax=ylim1,group=af(feature.type),fill=af(feature.type)),alpha=0.15,color=rgb(0,0,0,0.15)) +
      geom_text(data=my.bed,aes(x=(beg+end)/2,y=1.02*ylim1,label=feature,group=af(feature.type)),size=3,angle=90,hjust=0)
    
    #    p = p + geom_rect(data=my.bed,aes(x=beg,y=end,xmin=beg,xmax=end,ymin=0,ymax=maxY,fill=af()))
  }
  #   plot(NA,xlim=c(0,2500),ylim=c(0,1.4))
  # lines(test1$index,test1$CGdens,col='blue2')
  # lines(test1$index,test1$GCcont,col='green2')
  # lines(test1$index,test1$GCskew+0.5,col='red2')
  # segments(0,0.5,2500,0.5,lty=2)
  if (print == T) {
    print(p)
  } else {
    return(p)
  }
}
