setwd("D:/cygwin64/home/mitochy/INVIVO/PEAKS_LOCAL")

library(grid)
library(ggplot2)

beds = data.frame()
if (file.exists("annot.csv")) {
  beds = read.table("annot.csv",sep=",")
  colnames(beds) = c("gene","beg","end","feature","val","strand")
}

get_mybeds = function(beds, currgene) {
  mybeds = data.frame()
  for (bedsInd in 1:dim(beds)[1]) {
    currbeds = beds[bedsInd,]
    if (length(grep(currbeds$gene,currgene,ignore.case=TRUE)) > 0) {
      if (length(grep("ApaLI",currbeds$gene,ignore.case=TRUE)) > 0 & length(grep("ApaLI",currgene,ignore.case=TRUE)) > 0) {
        mybeds = rbind(mybeds,currbeds)
      } else if (length(grep("ApaLI",currbeds$gene,ignore.case=TRUE)) == 0 & length(grep("ApaLI",currgene,ignore.case=TRUE)) == 0) {
        mybeds = rbind(mybeds,currbeds)
      }
    }
  }
  return(mybeds)
}

maxgenenumber = 99999
files = dir("./",".BED$")
for (filesInd in 1:length(files)) {
  file = files[filesInd]
  filename = basename(file)
  filename = gsub("^(.+).BED$","\\1",filename,perl=T)
  
  df0 = read.table(file,sep="\t")
  #colnames(df0) = c("chr","beg","end","peakz","val","strand","file","mythres0","mythres1","peaktype")
  colnames(df0) = c("chr","beg","end","peakz","val","strand")
  df0 = df0[order(df0$beg, df0$end),]
  df0$file = "DUMMY"
  
  genes = unique(df0$chr)
  treats0 = unique(df0$file)
  treats = c()
  for (treatsInd in 1:length(treats0)) {
    temp = gsub("^.*(PCB[A-Za-z0-9_]+)_Z_10BUF_.+$","\\1",treats0[treatsInd],perl=T)
    treats = c(treats,temp)
  }
  treats = unique(treats)
  treats = treats[order(treats)]
  if (length(grep("sgRNA",file)) > 0) {
    treats = treats[grep("APALI",treats,invert = TRUE)]
    genes = genes[grep("APALI",genes,ignore.case = TRUE,invert = TRUE)]
  }
  for (genesInd in 1:length(genes)) {
    currgene = genes[genesInd]

    mybeds = get_mybeds(beds, currgene)
    
    for (treatsInd in 1:length(treats)) {
      currtreat = treats[treatsInd]
      df = df0[df0$chr == genes[genesInd],]
      df = df[grep(treats[treatsInd],df$file),]
      
      mytitle = paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],sep="")
      print(paste(mytitle,dim(df)[1]))
      cat("\n")
      if (dim(df)[1] < 8) {
        png(paste(mytitle,"_PEAKS.png",sep=""))
        plot(seq(1,10),seq(1,10),title=mytitle)
        text(5,5,paste('There is <8 peaks! (',dim(df)[1],')',sep=""))
        dev.off()
        next
      }
      
      if (dim(df)[1] > maxgenenumber) {
        df = df[sample(seq(1,dim(df)[1]),maxgenenumber,replace=F),]
      }
      
      rownames(df) = paste(seq(1,dim(df)[1]))
      
      dm = df;
      dm$cluster = 1#myclust[,1]
      dm$clusterbeg = 1
      dm$clusterend = 1
      
      dm$mid = as.integer((dm$beg + dm$end)/2/50)*50
      dm$begI = as.integer(dm$beg/2/50)*50
      dm$endI = as.integer(dm$end/2/50)*50
      dm$lenI = dm$endI - dm$begI
      
      dm$midI = as.integer((dm$beg + dm$end)/2/400)*400
      dm$mid2 = as.integer((dm$beg + dm$end)/2)
      dm = dm[order(dm$clusterbeg, dm$cluster,dm$begI,dm$beg,dm$endI,dm$end ),]
      dm = dm[order(dm$clusterend,dm$cluster,dm$midI,dm$mid2),]
      dm$y = seq(1,dim(dm)[1])
      myheight = dim(dm)[1]/400
      mywidth = 3000/400
      
      dm2 = dm
      dm2$cluster2 = 1
      
      p3 = ggplot(dm2,aes(x=beg,y=end)) +
        geom_point(aes(color=as.factor(cluster2)),shape=15,size=2) +
        theme_bw() + coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
        geom_segment(x=0,xend=3000,y=0,yend=3000) +
        ggtitle(paste(mytitle,"\nXY plot with no cluster (n=",dim(dm2)[1],")",sep="")) 
      
      if (dim(mybeds)[1] > 0) {
        p3 = p3 + geom_line(data=mybeds,aes(x=beg,xend=beg,y=0,yend=end),lty=2,color=rgb(1,1,1,0.5)) +
          geom_line(data=mybeds,aes(x=end,xend=end,y=0,yend=end),lty=2,color=rgb(1,1,1,0.5)) +
          geom_line(data=mybeds,aes(y=beg,yend=beg,x=0,xend=end),lty=2,color=rgb(1,1,1,0.5)) +
          geom_line(data=mybeds,aes(y=end,yend=end,x=0,xend=end),lty=2,color=rgb(1,1,1,0.5)) +
          geom_text(data=mybeds,aes(x=(beg+end)/2,y=0,label=feature),size=2,vjust=0) +
          geom_text(data=mybeds,aes(y=(beg+end)/2,x=0,label=feature),size=2,hjust=0)
      }
      
      
      #No cluster
      clusterbeg2 = as.data.frame(aggregate(dm2$beg,by=list(dm2$cluster2),min));colnames(clusterbeg2) = c("cluster2","clusterbeg2")
      clusterend2 = as.data.frame(aggregate(dm2$end,by=list(dm2$cluster2),max));colnames(clusterend2) = c("cluster2","clusterend2")
      
      dm2 = merge(dm2,clusterbeg2,by=c("cluster2"))
      dm2 = merge(dm2,clusterend2,by=c("cluster2"))
      
      dm2 = dm2[order(dm2$clusterbeg2,dm2$cluster2,dm2$begI,dm2$beg,dm2$endI,dm2$end),]
      dm2$y = seq(1,dim(dm2)[1])
      myheight = dim(dm2)[1]/400
      mywidth = 3000/400
      
      p4 = ggplot(dm2,aes(x=beg,y=y)) +
        geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster2))) +
        theme_bw() +
        annotate(geom="rect",xmin=641,xmax=841,ymin=0,ymax=dim(dm2)[1],lty=2,linewidth=1,alpha=0.15) +
        coord_cartesian(xlim=c(0,3000),ylim=c(0,dim(dm2)[1])) +
        geom_point(aes(x=mid,y=y),pch=".") + 
        ggtitle(paste(mytitle,"\nn=",dim(dm2)[1],sep="")) 

      if (dim(mybeds)[1] > 0) {
        p4 = p4 + geom_rect(data=mybeds,aes(x=0,y=0,xmin=beg,xmax=end,ymin=0,ymax=dim(dm2)[1]),fill=rgb(1,1,1,0.1)) +
          geom_text(data=mybeds,aes(x=(beg+end)/2,y=dim(dm2)[1],label=feature),size=2,vjust=1)
      }
      
      
      png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_PEAKSPLOT.png",sep=""),height=myheight*400, width=mywidth*400)
      print(p4)
      dev.off()
      
      png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_XYPLOT.png",sep=""),height=1000, width=1000)
      print(p3)
      dev.off()
      
      
#      png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_XYPLOT_ROTATED.png",sep=""),height=1000, width=1000)
#      pushViewport(viewport(name = "rotate", angle = -45, clip = "off", width = 0.7, height = 0.7))
#      print(p3,vp='rotate')
#      dev.off()
      
      next
      
      
      
    }
  }
}

get_genes = function(files, filesInd) {
  file = files[filesInd]
  filename = basename(file)
  filename = gsub("^(.+).BED$","\\1",filename,perl=T)
  
  df0 = read.table(file,sep="\t")
  #colnames(df0) = c("chr","beg","end","peakz","val","strand","file","mythres0","mythres1","peaktype")
  colnames(df0) = c("chr","beg","end","peakz","val","strand")
  df0 = df0[order(df0$beg, df0$end),]
  df0$file = "DUMMY"
  
  genes = unique(df0$chr)
  treats0 = unique(df0$file)
  treats = c()
  for (treatsInd in 1:length(treats0)) {
    temp = gsub("^.*(PCB[A-Za-z0-9_]+)_Z_10BUF_.+$","\\1",treats0[treatsInd],perl=T)
    treats = c(treats,temp)
  }
  treats = unique(treats)
  treats = treats[order(treats)]
  if (length(grep("sgRNA",file)) > 0) {
    treats = treats[grep("APALI",treats,invert = TRUE)]
    genes = genes[grep("APALI",genes,ignore.case = TRUE,invert = TRUE)]
  }
  return(genes)
}

get_dm2 = function(genes, geneWant) {
  for (i in 1:length(genes)) {
    if (genes[i] == geneWant) {
      genesInd = i
    }
  }
  currgene = genes[genesInd]
  
  mybeds = get_mybeds(beds, currgene)
  
  df = df0[df0$chr == genes[genesInd],]
  
  mytitle = paste(filename,"_gene",genes[genesInd],sep="")
  print(paste(mytitle,dim(df)[1]))
  cat("\n")
  if (dim(df)[1] < 8) {
    png(paste(mytitle,"_PEAKS.png",sep=""))
    plot(seq(1,10),seq(1,10),title=mytitle)
    text(5,5,paste('There is <8 peaks! (',dim(df)[1],')',sep=""))
    dev.off()
    next
  }
  
  if (dim(df)[1] > maxgenenumber) {
    df = df[sample(seq(1,dim(df)[1]),maxgenenumber,replace=F),]
  }
  
  rownames(df) = paste(seq(1,dim(df)[1]))
  
  dm = df;
  dm$cluster = 1#myclust[,1]
  dm$clusterbeg = 1
  dm$clusterend = 1
  
  dm$mid = as.integer((dm$beg + dm$end)/2/50)*50
  dm$begI = as.integer(dm$beg/2/50)*50
  dm$endI = as.integer(dm$end/2/50)*50
  dm$lenI = dm$endI - dm$begI
  
  dm$midI = as.integer((dm$beg + dm$end)/2/400)*400
  dm$mid2 = as.integer((dm$beg + dm$end)/2)
  dm = dm[order(dm$clusterbeg, dm$cluster,dm$begI,dm$beg,dm$endI,dm$end ),]
  dm = dm[order(dm$clusterend,dm$cluster,dm$midI,dm$mid2),]
  dm$y = seq(1,dim(dm)[1])
  myheight = dim(dm)[1]/400
  mywidth = 3000/400
  
  dm2 = dm
  dm2$cluster2 = 1
  
  p3 = ggplot(dm2,aes(x=beg,y=end)) +
    geom_point(aes(color=as.factor(cluster2)),shape=15,size=2) +
    theme_bw() + coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
    geom_segment(x=0,xend=3000,y=0,yend=3000) +
    ggtitle(paste(mytitle,"\nXY plot with no cluster (n=",dim(dm2)[1],")",sep="")) 
  
  if (dim(mybeds)[1] > 0) {
    p3 = p3 + geom_line(data=mybeds,aes(x=beg,xend=beg,y=0,yend=end),lty=2,color=rgb(1,1,1,0.5)) +
      geom_line(data=mybeds,aes(x=end,xend=end,y=0,yend=end),lty=2,color=rgb(1,1,1,0.5)) +
      geom_line(data=mybeds,aes(y=beg,yend=beg,x=0,xend=end),lty=2,color=rgb(1,1,1,0.5)) +
      geom_line(data=mybeds,aes(y=end,yend=end,x=0,xend=end),lty=2,color=rgb(1,1,1,0.5)) +
      geom_text(data=mybeds,aes(x=(beg+end)/2,y=0,label=feature),size=2,vjust=0) +
      geom_text(data=mybeds,aes(y=(beg+end)/2,x=0,label=feature),size=2,hjust=0)
  }
  
  
  #No cluster
  clusterbeg2 = as.data.frame(aggregate(dm2$beg,by=list(dm2$cluster2),min));colnames(clusterbeg2) = c("cluster2","clusterbeg2")
  clusterend2 = as.data.frame(aggregate(dm2$end,by=list(dm2$cluster2),max));colnames(clusterend2) = c("cluster2","clusterend2")
  
  dm2 = merge(dm2,clusterbeg2,by=c("cluster2"))
  dm2 = merge(dm2,clusterend2,by=c("cluster2"))
  
  dm2 = dm2[order(dm2$clusterbeg2,dm2$cluster2,dm2$begI,dm2$beg,dm2$endI,dm2$end),]
  dm2$y = seq(1,dim(dm2)[1])
  myheight = dim(dm2)[1]/400
  mywidth = 3000/400
  return(dm2)
}

genes = get_genes(files,2)
dm2 = get_dm2(genes,"CALM3")
dm4 = dm2[1:4]

get_dm3 = function(dm2, divby,threshold) {
  mybegs = c()
  myends = c()
  for (i in 1:(as.integer(5000/divby))) {
    mybegs[i] = 0
    myends[i] = 0
  }
  a = c()
  for (i in 1:dim(dm2)[1]) {
    begI = as.integer(dm2$beg[i]/divby)
    endI = as.integer(dm2$end[i]/divby)
    a[i] = begI
    mybegs[begI] = mybegs[begI] + 1 
    myends[endI] = myends[endI] + 1
  }
  dm3 = data.frame(x=seq(1,5000/divby),beg=mybegs,end=myends)
  dm3[dm3$beg < 10,]$beg = 0
  dm3[dm3$end < 10,]$end = 0
  return(dm3)
}

divby = 50

dm3 = get_dm3(dm2, divby, 10)
#dm3 = dm3[dm3$beg > 10 | dm3$end > 10,]

plot(NA,xlim=c(0,5000/divby),ylim=c(0,5000/divby))
points(dm2$beg/divby,dm2$end/divby,pch=18,col='grey')
lines(c(0,5000/divby),c(0,5000/divby))
lines(dm3$x,dm3$beg/10,col='blue4')
lines(dm3$end/10,dm3$x,col='red4')

xbegs = dm3[dm3$beg != 0,]$x
for (i in seq(length(xbegs),1,-1)) {
  xbeg0 = xbegs[i-1]
  xbeg1 = xbegs[i]
  print(paste(xbeg0,xbeg1))
  if (dim(dm4[as.integer(dm4$beg/divby) > xbeg0,])[1] == 0) {next}
  if (dim(dm4[as.integer(dm4$beg/divby) <= xbeg1,])[1] == 0) {next}
  if (dim(dm4[as.integer(dm4$beg/divby) > xbeg0 & as.integer(dm4$beg/divby) <= xbeg1,])[1] == 0) {next}
  end2 = dm4[as.integer(dm4$beg/divby) > xbeg0 & as.integer(dm4$beg/divby) <= xbeg1,]
#  if (dim(end2)[1] == 0) {next}
  dm3b = get_dm3(end2, divby, 10)
  dm3c = dm3b
  if (dim(dm3c[dm3c$end > 0,])[1] == 0) {next}
  if (dim(dm3c[dm3c$end > 0,])[1] > 0) {dm3c = dm3c[dm3c$end > 0,]}
  if (dim(dm3c[dm3c$end <= dm3c$x,])[1] > 0) {dm3c = dm3c[dm3c$end <= dm3c$x,]}
  dm3c$end = xbeg0 #dm3c$end/1 + xbeg0
  points(dm3c$end,dm3c$x,col='blue4',pch=15)
#  print(paste(dm3c$end,dm3c$x))
#  lines(dm3c$end,dm3c$x,col='purple4',lwd=1)
}

xends = dm3[dm3$end != 0,]$x
for (i in seq(length(xends),1,-1)) {
  xend0 = xends[i-1]
  xend1 = xends[i]
  print(paste(xend0,xend1,i,"/",length(xends)))
  if (dim(dm4[as.integer(dm4$end/divby) > xend0 & as.integer(dm4$end/divby) <= xend1,])[1] == 0) {next}
  beg2 = dm4[as.integer(dm4$end/divby) > xend0 & as.integer(dm4$end/divby) <= xend1,]
  dm3b = get_dm3(beg2, divby, 10)
  dm3c = dm3b
  if (dim(dm3c[dm3c$beg > 0,])[1] == 0) {next}
  if (dim(dm3c[dm3c$beg > 0,])[1] > 0) {dm3c = dm3c[dm3c$beg > 0,]}
  if (dim(dm3c[dm3c$beg >= dm3c$x,])[1] > 0) {dm3c = dm3c[dm3c$beg >= dm3c$x,]}
  dm3c$beg = xend0 #dm3c$end/1 + xbeg0
  points(dm3c$x,dm3c$beg,col='red4',pch=15)
  lines(c(0,dm3c$x),c(xend0,dm3c$beg),col='purple4',lwd=1)
  #print(paste(dm3c$beg,dm3c$x))
}

# xends = dm3[dm3$beg != 0,]$x
# for (i in 1:length(xends)) {
#   xend = xends[i]
#   beg2 = dm4[as.integer(dm4$end/divby) >= xend,]
#   dm3b = get_dm3(end2, 15, 5)
#   dm3b$beg = dm3b$beg + xend
#   dm3c = dm3b[dm3b$beg >= dm3b$x,]
#   points(dm3c$x,dm3c$beg,col='purple4',shape=15)
# }
dim(dm3[dm3$beg > 50,])
dim(dm3[dm3$end > 50,])
dm3[dm3$beg > 50,]
dm3[dm3$end > 50,]
#for (i in )