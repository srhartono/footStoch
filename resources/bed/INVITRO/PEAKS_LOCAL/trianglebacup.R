#setwd("D:/cygwin64/home/mitochy/INVIVO/PEAKS_LOCAL")
setwd("D:/cygwin64/home/mitochy/INVITRO/PEAKS_LOCAL")

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
get_mytitle = function(genes, geneWant) {
  for (i in 1:length(genes)) {
    if (genes[i] == geneWant) {
      genesInd = i
    }
  }
  currgene = genes[genesInd]
  
  mybeds = get_mybeds(beds, currgene)
  
  df = df0[df0$chr == genes[genesInd],]
  
  mytitle = paste(filename,"_gene",genes[genesInd],sep="")
  return(mytitle)
}
get_dm2 = function(files,filesInd, genes, geneWant) {
  for (i in 1:length(genes)) {
    if (genes[i] == geneWant) {
      genesInd = i
    }
  }
  currgene = genes[genesInd]

  file = files[filesInd]
  filename = basename(file)
  filename = gsub("^(.+).BED$","\\1",filename,perl=T)
  
  df0 = read.table(file,sep="\t")
  #colnames(df0) = c("chr","beg","end","peakz","val","strand","file","mythres0","mythres1","peaktype")
  colnames(df0) = c("chr","beg","end","peakz","val","strand")
  df0 = df0[order(df0$beg, df0$end),]
  df0$file = "DUMMY"
  
  mybeds = get_mybeds(beds, currgene)
  
  df = df0[df0$chr == genes[genesInd],]
  
  mytitle = paste(filename,"_gene",genes[genesInd],sep="")
  print(paste(mytitle,dim(df)[1]))
  # cat("\n")
  if (dim(df)[1] < 8) {
  #   png(paste(mytitle,"_PEAKS.png",sep=""))
  #   plot(seq(1,10),seq(1,10),title=mytitle)
  #   text(5,5,paste('There is <8 peaks! (',dim(df)[1],')',sep=""))
  #   dev.off()
    return(df)
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
  
  # p3 = ggplot(dm2,aes(x=beg,y=end)) +
  #   geom_point(aes(color=as.factor(cluster2)),shape=15,size=2) +
  #   theme_bw() + coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
  #   geom_segment(x=0,xend=3000,y=0,yend=3000) +
  #   ggtitle(paste(mytitle,"\nXY plot with no cluster (n=",dim(dm2)[1],")",sep="")) 
  # 
  # if (dim(mybeds)[1] > 0) {
  #   p3 = p3 + geom_line(data=mybeds,aes(x=beg,xend=beg,y=0,yend=end),lty=2,color=rgb(1,1,1,0.5)) +
  #     geom_line(data=mybeds,aes(x=end,xend=end,y=0,yend=end),lty=2,color=rgb(1,1,1,0.5)) +
  #     geom_line(data=mybeds,aes(y=beg,yend=beg,x=0,xend=end),lty=2,color=rgb(1,1,1,0.5)) +
  #     geom_line(data=mybeds,aes(y=end,yend=end,x=0,xend=end),lty=2,color=rgb(1,1,1,0.5)) +
  #     geom_text(data=mybeds,aes(x=(beg+end)/2,y=0,label=feature),size=2,vjust=0) +
  #     geom_text(data=mybeds,aes(y=(beg+end)/2,x=0,label=feature),size=2,hjust=0)
  # }
  # 
  # 
  #No cluster
  # clusterbeg2 = as.data.frame(aggregate(dm2$beg,by=list(dm2$cluster2),min));colnames(clusterbeg2) = c("cluster2","clusterbeg2")
  # clusterend2 = as.data.frame(aggregate(dm2$end,by=list(dm2$cluster2),max));colnames(clusterend2) = c("cluster2","clusterend2")
  
  # dm2 = merge(dm2,clusterbeg2,by=c("cluster2"))
  # dm2 = merge(dm2,clusterend2,by=c("cluster2"))
  # 
  dm2 = dm2[order(dm2$begI,dm2$beg,dm2$endI,dm2$end),]
  dm2$y = seq(1,dim(dm2)[1])
  myheight = dim(dm2)[1]/400
  mywidth = 3000/400
  return(dm2)
}

get_dm3 = function(dm2, divby,threshold) {
  mybegs = c()
  myends = c()
  for (i in 1:(as.integer(5000/divby))) {
    mybegs[i] = 0
    myends[i] = 0
  }
  for (i in 1:dim(dm2)[1]) {
    begI = as.integer(dm2$beg[i]/divby)
    endI = as.integer(dm2$end[i]/divby)
    mybegs[begI] = mybegs[begI] + 1 
    myends[endI] = myends[endI] + 1
  }
  dm3 = data.frame(x=as.integer(seq(1,5000/divby)),beg=mybegs,end=myends)
  dm3[dm3$beg < threshold,]$beg = 0
  dm3[dm3$end < threshold,]$end = 0
  lastend = 0
  lastendx = -1
  lastbeg = 0
  lastbegx = -1
  for (i in seq(1,dim(dm3)[1],1)) {
  #for (i in seq(dim(dm3)[1],1,-1)) {
    if (lastbegx == -1) {
      lastbeg = dm3$beg[i]
      lastbegx = i
    } else {
      if (lastbeg >= threshold & dm3$beg[i] >= threshold) {
        print(paste("lastbeg=",lastbeg," dm3begi=",dm3$beg[i],sep=""))
        lastbegx = i
        dm3$beg[lastbegx] = dm3$beg[lastbegx] + dm3$beg[i]
        dm3$beg[i] = 0
      } else {
        lastbeg = dm3$beg[i]
        lastbegx = i
      }
    }
  }
  
  for (i in seq(dim(dm3)[1],1,-1)) {
    if (lastendx == -1) {
      lastend = dm3$end[i]
      lastendx = i
#    if (i == dim(dm3)[1]) {
#      lastend = dm3$end[i]
    } else {
      if (lastend >= threshold & dm3$end[i] >= threshold) {
        print(paste("lastend=",lastend," dm3endi=",dm3$end[i],sep=""))
        lastend = dm3$end[i]
        dm3$end[lastendx] = dm3$end[lastendx] + dm3$end[i]
        dm3$end[i] = 0
      } else {
        lastend = dm3$end[i]
        lastendx = i
      }
    }
  }
  return(dm3)
}


initgraph = function(dm2, dm3,divby,mygene="NA") {
  plot(NA,xlim=c(0,5000/divby),ylim=c(0,5000/divby),xlab="beg",ylab="end",main=mygene)
  points(dm2$beg/divby,dm2$end/divby,pch=".",col='grey')
  lines(c(0,5000/divby),c(0,5000/divby))
  lines(dm3$x,log(dm3$beg+1,5),col='red4')
  lines(log(dm3$end+1,5),dm3$x,col='blue4')
  
  dm30 = dm3[as.integer(dm3$beg/5) != 0,]
  text(dm30$x,0,labels = paste(as.integer(dm30$x),sep=""),cex=0.5,adj=c(0.5,1))

  dm30 = dm3[as.integer(dm3$end/5) != 0,]
  text(0,dm30$x,labels = paste(as.integer(dm30$x),sep=""),cex=0.5,adj=c(1,0.5))

#  lines(dm3$x,dm3$beg/(divby/5)10,col='red4')
}

filesInd = 1
genes = get_genes(files,filesInd)
genes = genes[order(genes)]
border0 = 0
border1 = 0
plus = "+"
if (border1 < 0) {plus = "-"}
genewant = "T7_init_VR_20"
mytitle = get_mytitle(genes,genewant)
mytitle = paste(mytitle,"\n(x & y border: ",-1 * border0,"/",plus,border1,")",sep="")
dm2 = get_dm2(files, filesInd, genes,genewant)
dm4 = dm2[1:4]

divby = 16

mythres = 200
mythres2 = 200

dm3 = get_dm3(dm2, divby, mythres)
myx = c()
myy = c()
mytype = c()
initgraph(dm2, dm3, divby,mytitle)

initgraph(dm2, dm3, divby,mytitle)
xbegs = dm3[dm3$beg >= mythres,]$x
xbegs = c(0,xbegs,max(dm3$x))
for (i in seq(length(xbegs),1,-1)) {
  xbeg0 = xbegs[i-1]
  xbeg1 = xbegs[i]
  #if (dim(dm3[dm3$x > xbeg0 & dm3$x <= xbeg1 & dm3$end > mythres,])[1] == 0) {next}
  #end2 = dm3[dm3$x > xbeg0 & dm3$x <= xbeg1 & dm3$end > mythres,]
  if (dim(dm4[as.integer(dm4$beg/divby) > xbeg0 & as.integer(dm4$beg/divby) <= xbeg1,])[1] == 0) {next}
  end2 = dm4[as.integer(dm4$beg/divby) > xbeg0 & as.integer(dm4$beg/divby) <= xbeg1,]
  dm3b = get_dm3(end2, divby, mythres2)
  dm3b = dm3b[dm3b$end > 0,]
  #head(dm3b[dm3b$beg > 0,])
  #dm3b = end2
  print(paste(xbeg0,xbeg1,dm3b[,1], ":", dm3b[,3]))
  # print(paste(xbeg0,xbeg1,",",dm3b[dm3b$end > mythres,]$x,":",dm3b[dm3b$end > mythres,]$end,sep=""))
  dm3c = dm3b
  if (dim(dm3c[dm3c$end >= mythres2,])[1] == 0) {next}
  if (dim(dm3c[dm3c$end >= mythres2,])[1] > 0) {dm3c = dm3c[dm3c$end >= mythres2,]}
  dm3c$end = xbeg0 
 for (j in 1:length(dm3c$x)) {
    curry = dm3c$x[j]
    currx0 = xbeg0
    currx1 = xbeg1
    if (currx1 > as.integer(curry - 100/divby)) {#as.integer(min(beg2$end/divby))) {
      currx1 = as.integer(curry - 100/divby)
    }
    if (currx0 < min(dm2$beg/divby)) {
      currx0 = as.integer(min(dm2$beg/divby))
    }
    if (currx1 < min(dm2$beg/divby)) {
      currx1 = as.integer(min(dm2$beg/divby))
    }
    points(currx0,curry,col=rgb(0,0,1,1),pch=15)
    lines(x=c(currx0,currx1),y=c(curry,curry),col=rgb(0,0,0,0.25))
    points(currx1,curry,col=rgb(0,0,1,1),pch=18)
    myx = c(myx,currx0,currx1)
    myy = c(myy,curry,curry)
    mytype = c(mytype,'beg0','beg1')

  }
  #    points(dm3c$end,dm3c$x,col=rgb(0,0,1,1),pch=15)

}

xends = dm3[dm3$end >= mythres,]$x
xends = c(as.integer(min(dm2$beg/divby)),xends,max(dm3$x))

for (i in seq(length(xends),1,-1)) {
  xend0 = xends[i-1]
  xend1 = xends[i]
  if (dim(dm4[as.integer(dm4$end/divby) > xend0 & as.integer(dm4$end/divby) <= xend1,])[1] == 0) {next}
  beg2 = dm4[as.integer(dm4$end/divby) > xend0 & as.integer(dm4$end/divby) <= xend1,]
  dm3b = get_dm3(beg2, divby, mythres2)
  dm3b = dm3b[dm3b$beg > 0,]
  mymin = as.integer(min(dm3[dm3$end >= mythres2,]$x))
  print(paste("i=",i,",coords=",xend0,xend1,dm3b[,1], ":", dm3b[,3], mymin))
  dm3c = dm3b
  if (dim(dm3c[dm3c$beg >= mythres2,])[1] == 0) {next}
  if (dim(dm3c[dm3c$beg >= mythres2,])[1] > 0) {dm3c = dm3c[dm3c$beg >= mythres2,]}
  dm3c$xend0 = xend0
  for (j in 1:length(dm3c$x)) {
    currx = dm3c$x[j]
    curry0 = xend0
    curry1 = xend1
    if (curry0 < as.integer(currx + 100/divby)) {#as.integer(min(beg2$end/divby))) {
      curry0 = as.integer(currx + 100/divby)
    }
    points(currx,curry0,col=rgb(1,0,0,1),pch=15)
    lines(x=c(currx,currx),y=c(curry0,curry1),col=rgb(0,0,0,0.25))
    points(currx,curry1,col=rgb(1,0,0,1),pch=18)
    myx = c(myx,currx,currx)
    myy = c(myy,curry0,curry1)
    mytype = c(mytype,'end0','end1')
  }
    
  #points(dm3c$x,dm3c$xend1,col=rgb(1,0,0,1),pch=18)
  #lines(c(dm3c$x,dm3c$beg),c(dm3c$x,dm3c$beg2),col=rgb(1,0,0,1))
}

final0 = data.frame(x=myx, y=myy, type=mytype)
final0 = final0[order(final0$x, final0$y),]

myx0 = c()
myy0 = c()
myx1 = c()
myy1 = c()
mypaste = c()
for (i in 1:length(final0$x)) {
  x0 = final0$x[i]
  y0 = final0$y[i]

  for (j in 1:length(final0$y)) {
    x1 = final0$x[j]
    y1 = final0$y[j]
    if (x1 == x0 & y0 == y1) {next}
    if (length(grep(paste(x0,y0,x1,y1),mypaste)) == 0 & x0 >= x1 - border0 & x0 <= x1 + border1) {
      myx0 = c(myx0, x0)
      myy0 = c(myy0, y0)
      myx1 = c(myx1, x1)
      myy1 = c(myy1, y1)
      mypaste = c(mypaste, paste(x0,y0,x1,y1))
    }
  }
}
finalbeg = data.frame(x0=myx0, y0=myy0, x1=myx1,y1=myy1)
finalbeg = finalbeg[order(finalbeg$x0, finalbeg$y0,finalbeg$x1,finalbeg$y1),]

segments(finalbeg$x0,finalbeg$y0,finalbeg$x1,finalbeg$y1,col=rgb(1,0,0,1),pch=18)

myx0 = c()
myy0 = c()
myx1 = c()
myy1 = c()
mypaste = c()
for (i in 1:length(final0$x)) {
  x0 = final0$x[i]
  y0 = final0$y[i]

  for (j in 1:length(final0$y)) {
    x1 = final0$x[j]
    y1 = final0$y[j]
    if (x1 == x0 & y0 == y1) {next}
    if (length(grep(paste(x0,y0,x1,y1),mypaste)) == 0 & y0 >= y1 -  border0 & y0 <= y1 + border1) {
      myx0 = c(myx0, x0)
      myy0 = c(myy0, y0)
      myx1 = c(myx1, x1)
      myy1 = c(myy1, y1)
      mypaste = c(mypaste, paste(x0,y0,x1,y1))
    }
  }
}
finalend = data.frame(x0=myx0, y0=myy0, x1=myx1,y1=myy1)
finalend = finalend[order(finalend$x0, finalend$y0,finalend$x1,finalend$y1),]

segments(finalend$x0,finalend$y0,finalend$x1,finalend$y1,col=rgb(0,0,1,1),pch=18)



# mythres = 10
# xends = dm3[dm3$end >= mythres,]$x
# #-------xend0-----------xend1-------
# for (i in seq(length(xends),1,-1)) {
#   xend0 = xends[i]
# 
#   if (dim(dm4[as.integer(dm4$end/divby) == xend0,])[1] == 0) {next}
#   beg2 = dm4[as.integer(dm4$end/divby) == xend0,]
#   dm3b = get_dm3(beg2, divby, mythres)
#   print(paste(xend0,",",dm3b[dm3b$beg > mythres,]$x,":",dm3b[dm3b$beg > mythres,]$beg,sep=""))
#   dm3c = dm3b
#   if (dim(dm3c[dm3c$beg > mythres,])[1] == 0) {next}
#   if (dim(dm3c[dm3c$beg > mythres,])[1] > 0) {dm3c = dm3c[dm3c$beg > 0,]}
#   dm3c$beg = xend0 
#   points(dm3c$x,dm3c$beg,col=rgb(0,0,1,1),pch=15)
# 
# }
# 
# for (i in seq(length(xbegs),1,-1)) {
#   xbeg0 = xbegs[i-1]
#   xbeg1 = xbegs[i]
#   print(paste(xbeg0,xbeg1))
#   if (dim(dm4[as.integer(dm4$beg/divby) > xbeg0 & as.integer(dm4$beg/divby) <= xbeg1,])[1] == 0) {next}
#   end2 = dm4[as.integer(dm4$beg/divby) >= xbeg0 & as.integer(dm4$beg/divby) <= xbeg1,]
#   dm3b = get_dm3(end2, divby, 10)
#   head(dm3b[dm3b$end > 0,])
#   dm3c = dm3b
#   if (dim(dm3c[dm3c$end > 0,])[1] == 0) {next}
#   if (dim(dm3c[dm3c$end > 0,])[1] > 0) {dm3c = dm3c[dm3c$end > 0,]}
#   if (dim(dm3c[dm3c$end <= dm3c$x,])[1] > 0) {dm3c = dm3c[dm3c$end <= dm3c$x,]}
#   dm3c$end = xbeg1 #dm3c$end/1 + xbeg0
#   points(dm3c$end,dm3c$x,col=rgb(1,0,0,1),pch=15)
#   
#   dm3b = get_dm3(end2, divby, 10)
#   dm3c = dm3b
#   if (dim(dm3c[dm3c$end > 0,])[1] == 0) {next}
#   if (dim(dm3c[dm3c$end > 0,])[1] > 0) {dm3c = dm3c[dm3c$end > 0,]}
#   if (dim(dm3c[dm3c$end <= dm3c$x,])[1] > 0) {dm3c = dm3c[dm3c$end <= dm3c$x,]}
#   dm3c$end = xbeg1 #dm3c$end/1 + xbeg0
#   points(dm3c$end,dm3c$x,col=rgb(0,0,1,1),pch=18)
#   #  print(paste(dm3c$end,dm3c$x))
# #  lines(dm3c$end,dm3c$x,col='purple4',lwd=0.5)
# }
# 
# xends = dm3[dm3$end != 0,]$x
# for (i in seq(length(xends),1,-1)) {
#   xend0 = xends[i-1]
#   if (i == 1) {
#     xend0 = 0
#   }
#   xend1 = xends[i]
#   print(paste(xend0,xend1,i,"/",length(xends)))
#   if (dim(dm4[as.integer(dm4$end/divby) > xend0 & as.integer(dm4$end/divby) <= xend1,])[1] == 0) {next}
#   beg2 = dm4[as.integer(dm4$end/divby) >= xend0 & as.integer(dm4$end/divby) <= xend1,]
#   dm3b = get_dm3(beg2, divby, 10)
#   dm3c = dm3b
#   if (dim(dm3c[dm3c$beg > 0,])[1] == 0) {next}
#   if (dim(dm3c[dm3c$beg > 0,])[1] > 0) {dm3c = dm3c[dm3c$beg > 0,]}
#   if (dim(dm3c[dm3c$beg >= dm3c$x,])[1] > 0) {dm3c = dm3c[dm3c$beg >= dm3c$x,]}
#   dm3c$beg = xend0 #dm3c$beg/1 + xend0
#   if (i == 1) {
#     dm3c$beg = min(beg2$end/divby)
#   }
#   points(dm3c$x,dm3c$beg,col='red4',pch=15)
#   #lines(c(0,dm3c$x),c(xend0,dm3c$beg),col='purple4',lwd=1)
# }
# 

# xends = dm3[dm3$beg != 0,]$x
# for (i in 1:length(xends)) {
#   xend = xends[i]
#   beg2 = dm4[as.integer(dm4$end/divby) >= xend,]
#   dm3b = get_dm3(end2, 15, 5)
#   dm3b$beg = dm3b$beg + xend
#   dm3c = dm3b[dm3b$beg >= dm3b$x,]
# #   points(dm3c$x,dm3c$beg,col='purple4',shape=15)
# # }
# dim(dm3[dm3$beg > 50,])
# dim(dm3[dm3$end > 50,])
# dm3[dm3$beg > 50,]
# dm3[dm3$end > 50,]
# #for (i in )

# colors = c("red4","blue4","green4")
# xbegs = dm3[dm3$beg != 0,]$x
# for (i in seq(length(xbegs),1,-1)) {
#   xbeg0 = xbegs[i]
#   if (dim(dm4[as.integer(dm4$beg/divby) > xbeg0,])[1]== 0) {next}
#   end2 = dm4[as.integer(dm4$beg/divby) >= xbeg0,]
#   
#   temp = kmeans(as.integer(end2$end/divby),3)
#   end2 = cbind(end2, temp$cluster)
#   for (j in 1:max(unique(end2$cluster))) {
#     temp = end2[end2$cluster == j,]
#     points(temp$beg/divby,temp$end/divby,col=colors[j],pch=15)
#   }
#   #end3min = quantile(as.integer(end2$end/divby),probs = 0.01)
#   #end3max = quantile(as.integer(end2$end/divby),probs = 0.99)
#   #end2 = end2[end2$end/divby > end3min & end2$end/divby < end3max,]
# #  points(xbeg0,min(end2$end/divby),col='blue4',pch=15)
# }
# 
# 
# xends = dm3[dm3$end != 0,]$x
# for (i in seq(length(xends),1,-1)) {
#   xend0 = xends[i]
#   if (dim(dm4[as.integer(dm4$end/divby) > xend0,])[1]== 0) {next}
#   beg2 = dm4[as.integer(dm4$end/divby) >= xend0,]
#   beg3min = quantile(as.integer(beg2$beg/divby),probs = 0.05)
#   beg3max = quantile(as.integer(beg2$beg/divby),probs = 0.95)
#   beg2 = beg2[beg2$beg/divby > beg3min & beg2$beg/divby < beg3max,]
#   points(min(beg2$beg/divby),xend0,col='red4',pch=15)
#   print(dim(beg2)[1])
# }
# 
# a = unique(as.integer(end2$end/divby)); a[order(a)]


#dm3 = dm3[dm3$beg > 10 | dm3$end > 10,]

# plot(NA,xlim=c(0,5000/divby),ylim=c(0,5000/divby))
# points(dm2$beg/divby,dm2$end/divby,pch=18,col='grey')
# lines(c(0,5000/divby),c(0,5000/divby))
# #points(dm3[dm3$end != 0,]$end/divby,dm3[dm3$beg != 0,]$x,col='blue4')
# #points(dm3[dm3$beg != 0,]$x,dm3[dm3$beg != 0,]$beg/divby,col='red4')
# lines(dm3$x,dm3$beg/5,col='red4')
# lines(dm3$end/5,dm3$x,col='blue4')



# 
# for (filesInd in 1:length(files)) {
#   file = files[filesInd]
#   filename = basename(file)
#   filename = gsub("^(.+).BED$","\\1",filename,perl=T)
#   
#   df0 = read.table(file,sep="\t")
#   #colnames(df0) = c("chr","beg","end","peakz","val","strand","file","mythres0","mythres1","peaktype")
#   colnames(df0) = c("chr","beg","end","peakz","val","strand")
#   df0 = df0[order(df0$beg, df0$end),]
#   df0$file = "DUMMY"
#   
#   genes = unique(df0$chr)
#   treats0 = unique(df0$file)
#   treats = c()
#   for (treatsInd in 1:length(treats0)) {
#     temp = gsub("^.*(PCB[A-Za-z0-9_]+)_Z_10BUF_.+$","\\1",treats0[treatsInd],perl=T)
#     treats = c(treats,temp)
#   }
#   treats = unique(treats)
#   treats = treats[order(treats)]
#   if (length(grep("sgRNA",file)) > 0) {
#     treats = treats[grep("APALI",treats,invert = TRUE)]
#     genes = genes[grep("APALI",genes,ignore.case = TRUE,invert = TRUE)]
#   }
#   for (genesInd in 1:length(genes)) {
#     currgene = genes[genesInd]
#     
#     mybeds = get_mybeds(beds, currgene)
#     
#     for (treatsInd in 1:length(treats)) {
#       currtreat = treats[treatsInd]
#       df = df0[df0$chr == genes[genesInd],]
#       df = df[grep(treats[treatsInd],df$file),]
#       
#       mytitle = paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],sep="")
#       print(paste(mytitle,dim(df)[1]))
#       cat("\n")
#       if (dim(df)[1] < 8) {
#         png(paste(mytitle,"_PEAKS.png",sep=""))
#         plot(seq(1,10),seq(1,10),title=mytitle)
#         text(5,5,paste('There is <8 peaks! (',dim(df)[1],')',sep=""))
#         dev.off()
#         next
#       }
#       
#       if (dim(df)[1] > maxgenenumber) {
#         df = df[sample(seq(1,dim(df)[1]),maxgenenumber,replace=F),]
#       }
#       
#       rownames(df) = paste(seq(1,dim(df)[1]))
#       
#       dm = df;
#       dm$cluster = 1#myclust[,1]
#       dm$clusterbeg = 1
#       dm$clusterend = 1
#       
#       dm$mid = as.integer((dm$beg + dm$end)/2/50)*50
#       dm$begI = as.integer(dm$beg/2/50)*50
#       dm$endI = as.integer(dm$end/2/50)*50
#       dm$lenI = dm$endI - dm$begI
#       
#       dm$midI = as.integer((dm$beg + dm$end)/2/400)*400
#       dm$mid2 = as.integer((dm$beg + dm$end)/2)
#       dm = dm[order(dm$clusterbeg, dm$cluster,dm$begI,dm$beg,dm$endI,dm$end ),]
#       dm = dm[order(dm$clusterend,dm$cluster,dm$midI,dm$mid2),]
#       dm$y = seq(1,dim(dm)[1])
#       myheight = dim(dm)[1]/400
#       mywidth = 3000/400
#       
#       dm2 = dm
#       dm2$cluster2 = 1
#       
#       p3 = ggplot(dm2,aes(x=beg,y=end)) +
#         geom_point(aes(color=as.factor(cluster2)),shape=15,size=2) +
#         theme_bw() + coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
#         geom_segment(x=0,xend=3000,y=0,yend=3000) +
#         ggtitle(paste(mytitle,"\nXY plot with no cluster (n=",dim(dm2)[1],")",sep="")) 
#       
#       if (dim(mybeds)[1] > 0) {
#         p3 = p3 + geom_line(data=mybeds,aes(x=beg,xend=beg,y=0,yend=end),lty=2,color=rgb(1,1,1,0.5)) +
#           geom_line(data=mybeds,aes(x=end,xend=end,y=0,yend=end),lty=2,color=rgb(1,1,1,0.5)) +
#           geom_line(data=mybeds,aes(y=beg,yend=beg,x=0,xend=end),lty=2,color=rgb(1,1,1,0.5)) +
#           geom_line(data=mybeds,aes(y=end,yend=end,x=0,xend=end),lty=2,color=rgb(1,1,1,0.5)) +
#           geom_text(data=mybeds,aes(x=(beg+end)/2,y=0,label=feature),size=2,vjust=0) +
#           geom_text(data=mybeds,aes(y=(beg+end)/2,x=0,label=feature),size=2,hjust=0)
#       }
#       
#       
#       #No cluster
#       clusterbeg2 = as.data.frame(aggregate(dm2$beg,by=list(dm2$cluster2),min));colnames(clusterbeg2) = c("cluster2","clusterbeg2")
#       clusterend2 = as.data.frame(aggregate(dm2$end,by=list(dm2$cluster2),max));colnames(clusterend2) = c("cluster2","clusterend2")
#       
#       dm2 = merge(dm2,clusterbeg2,by=c("cluster2"))
#       dm2 = merge(dm2,clusterend2,by=c("cluster2"))
#       
#       dm2 = dm2[order(dm2$clusterbeg2,dm2$cluster2,dm2$begI,dm2$beg,dm2$endI,dm2$end),]
#       dm2$y = seq(1,dim(dm2)[1])
#       myheight = dim(dm2)[1]/400
#       mywidth = 3000/400
#       
#       p4 = ggplot(dm2,aes(x=beg,y=y)) +
#         geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster2))) +
#         theme_bw() +
#         annotate(geom="rect",xmin=641,xmax=841,ymin=0,ymax=dim(dm2)[1],lty=2,linewidth=1,alpha=0.15) +
#         coord_cartesian(xlim=c(0,3000),ylim=c(0,dim(dm2)[1])) +
#         geom_point(aes(x=mid,y=y),pch=".") + 
#         ggtitle(paste(mytitle,"\nn=",dim(dm2)[1],sep="")) 
#       
#       if (dim(mybeds)[1] > 0) {
#         p4 = p4 + geom_rect(data=mybeds,aes(x=0,y=0,xmin=beg,xmax=end,ymin=0,ymax=dim(dm2)[1]),fill=rgb(1,1,1,0.1)) +
#           geom_text(data=mybeds,aes(x=(beg+end)/2,y=dim(dm2)[1],label=feature),size=2,vjust=1)
#       }
#       
#       
#       png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_PEAKSPLOT.png",sep=""),height=myheight*400, width=mywidth*400)
#       print(p4)
#       dev.off()
#       
#       png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_XYPLOT.png",sep=""),height=1000, width=1000)
#       print(p3)
#       dev.off()
#       
#       
#       #      png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_XYPLOT_ROTATED.png",sep=""),height=1000, width=1000)
#       #      pushViewport(viewport(name = "rotate", angle = -45, clip = "off", width = 0.7, height = 0.7))
#       #      print(p3,vp='rotate')
#       #      dev.off()
#       
#       next
#       
#       
#       
#     }
#   }
# }



