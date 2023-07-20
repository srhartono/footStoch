setwd("D:/cygwin64/home/mitochy/INVIVO/PEAKS_LOCAL")

library(grid)
library(tidyverse)
library(ggplot2)
#library(dplyr)
library(cluster)
library(factoextra)

beds = read.table("annot.csv",sep=",")
colnames(beds) = c("gene","beg","end","feature","val","strand")


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
    #print(paste(filename,"_gene",genes[genesInd],sep=""))
    #print(mybeds)
    
    for (treatsInd in 1:length(treats)) {
      currtreat = treats[treatsInd]
      df = df0[df0$chr == genes[genesInd],]
      df = df[grep(treats[treatsInd],df$file),]
      
      mytitle = paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],sep="")
      print(paste(mytitle,dim(df)[1]))
      cat("\n")
      if (dim(df)[1] < 8) {
        pdf(paste(mytitle,"_PEAKS.pdf",sep=""))
        plot(seq(1,10),seq(1,10),title=mytitle)
        text(5,5,paste('There is <8 peaks! (',dim(df)[1],')',sep=""))
        dev.off()
        
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
      # 
      # df1k = subset(df,select=c("beg","end"))
      # mypca = prcomp(df1k, center = TRUE, scale = TRUE)
      # summary(mypca)
      # 
      #mytransform = as.data.frame(-mypca$x[,1:2])
      #
      #
      #                       k1 = fviz_nbclust(mytransform, kmeans, method = 'wss')
      #                       k2 = fviz_nbclust(mytransform, kmeans, method = 'silhouette')
      #                       k3 = fviz_nbclust(mytransform, kmeans, method = 'gap_stat')
      #
      #k = 3
      #kmeans_iris = kmeans(mytransform, centers = k, nstart = 50)
      #k3c = fviz_cluster(kmeans_iris, data = mytransform)
      #
      #                       k = 4
      #                       kmeans_iris = kmeans(mytransform, centers = k, nstart = 50)
      #                       k4c = fviz_cluster(kmeans_iris, data = mytransform)
      #
      #                       k = 5
      #                       kmeans_iris = kmeans(mytransform, centers = k, nstart = 50)
      #                       k5c = fviz_cluster(kmeans_iris, data = mytransform)
      #
      #                       k = 6
      #                       kmeans_iris = kmeans(mytransform, centers = k, nstart = 50)
      #                       k6c = fviz_cluster(kmeans_iris, data = mytransform)
      #
      #                       k = 7
      #                       kmeans_iris = kmeans(mytransform, centers = k, nstart = 50)
      #                       k7c = fviz_cluster(kmeans_iris, data = mytransform)
      
      
      #dev.off()
      
      
      #myclust = as.data.frame(kmeans_iris$cluster)
      
      
      dm = df;
      dm$cluster = 1#myclust[,1]
      dm$clusterbeg = 1
      dm$clusterend = 1
      
      #clusterbeg = as.data.frame(aggregate(dm$beg,by=list(dm$cluster),min));colnames(clusterbeg) = c("cluster","clusterbeg")
      #clusterend = as.data.frame(aggregate(dm$end,by=list(dm$cluster),max));colnames(clusterend) = c("cluster","clusterend")
      
      dm = dm #merge(dm,clusterbeg,by=c("cluster"))
      #                       dm = merge(dm,clusterend,by=c("cluster"))
      
      dm$mid = as.integer((dm$beg + dm$end)/2/50)*50
      dm$begI = as.integer(dm$beg/2/50)*50
      dm$endI = as.integer(dm$end/2/50)*50
      dm$lenI = dm$endI - dm$begI
      
      
      dm = dm[order(dm$clusterbeg, dm$cluster,dm$begI,dm$beg,dm$endI,dm$end ),]
      dm$y = seq(1,dim(dm)[1])
      myheight = dim(dm)[1]/400
      mywidth = 3000/400
      
      # p1 = ggplot(dm,aes(x=beg,y=y)) +
      #   geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster))) +
      #   theme_bw() + ggtitle(paste(mytitle,"\nOrdered by beg",sep=""))
      # 
      dm = dm[order(dm$clusterend,dm$cluster,dm$endI,dm$end,dm$begI,dm$beg),]
      dm$y = seq(1,dim(dm)[1])
      myheight = dim(dm)[1]/400
      mywidth = 3000/400
      # 
      # p2 = ggplot(dm,aes(x=beg,y=y)) +
      #   geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster))) +
      #   theme_bw() + ggtitle(paste(mytitle,"\nOrdered by end",sep=""))
      # 
      dm$midI = as.integer((dm$beg + dm$end)/2/400)*400
      dm$mid = as.integer((dm$beg + dm$end)/2)
      dm = dm[order(dm$clusterend,dm$cluster,dm$midI,dm$mid),]
      dm$y = seq(1,dim(dm)[1])
      myheight = 3000/400
      mywidth = 3000/400
      
      
      # p2b = ggplot(dm,aes(x=beg,y=y)) +
      #   geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster))) +
      #   theme_bw() + ggtitle(paste(mytitle,"\nOrdered by mid",sep=""))
      # 
      # 
      # 
      # p2c = ggplot(dm,aes(x=beg,y=end)) +
      #   geom_point(aes(color=as.factor(cluster)),pch=".") +
      #   theme_bw() + coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
      #   geom_segment(x=0,xend=3000,y=0,yend=3000) +
      #   ggtitle(paste(mytitle,"\nXY plot with k=3 clusters",sep=""))
      # 
      
      dm2 = dm
      dm2$cluster2 = 1
      
      p3 = ggplot(dm2,aes(x=beg,y=end)) +
        geom_point(aes(color=as.factor(cluster2)),pch=".") +
        theme_bw() + coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
        geom_segment(x=0,xend=3000,y=0,yend=3000) +
        ggtitle(paste(mytitle,"\nXY plot with no cluster",sep="")) 
        # geom_line(data=mybeds,aes(x=beg,xend=beg,y=0,yend=end),lty=2,color=rgb(1,1,1,0.5)) +
        # geom_line(data=mybeds,aes(x=end,xend=end,y=0,yend=end),lty=2,color=rgb(1,1,1,0.5)) +
        # geom_line(data=mybeds,aes(y=beg,yend=beg,x=0,xend=end),lty=2,color=rgb(1,1,1,0.5)) +
        # geom_line(data=mybeds,aes(y=end,yend=end,x=0,xend=end),lty=2,color=rgb(1,1,1,0.5)) +
        # geom_text(data=mybeds,aes(x=(beg+end)/2,y=0,label=feature),size=2,vjust=0) +
        # geom_text(data=mybeds,aes(y=(beg+end)/2,x=0,label=feature),size=2,hjust=0)
        # 
      
      
      
      #No cluster
      clusterbeg2 = as.data.frame(aggregate(dm2$beg,by=list(dm2$cluster2),min));colnames(clusterbeg2) = c("cluster2","clusterbeg2")
      clusterend2 = as.data.frame(aggregate(dm2$end,by=list(dm2$cluster2),max));colnames(clusterend2) = c("cluster2","clusterend2")
      
      #dm$mid = as.integer((dm$beg + dm$end)/2/50)*50
      #dm$begI = as.integer(dm$beg/2/50)*50
      #dm$endI = as.integer(dm$end/2/50)*50
      #dm$lenI = dm$endI - dm$begI
      
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
        geom_point(aes(x=mid,y=y),pch=".") + ggtitle(paste(filename,"_",genes[genesInd],"_PEAKS",sep=""))
        #geom_rect(data=mybeds,aes(x=0,y=0,xmin=beg,xmax=end,ymin=0,ymax=dim(dm2)[1]),fill=rgb(1,1,1,0.1)) +
        #geom_text(data=mybeds,aes(x=(beg+end)/2,y=dim(dm2)[1],label=feature),size=2,vjust=1)
      
      #       png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_PEAKS.png",sep=""),height=myheight*400, width=mywidth*400)
      #       print(k3c)
      #       dev.off()
      
      png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_PEAKSPLOT.png",sep=""),height=myheight*400, width=mywidth*400)
      print(p4)
      dev.off()
      
      png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_XYPLOT.png",sep=""),height=1000, width=1000)
      print(p3)
      dev.off()
      
      png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_XYPLOT_ROTATED.png",sep=""),height=1000, width=1000)
      pushViewport(viewport(name = "rotate", angle = 45, clip = "off", width = 0.7, height = 0.7))
      print(p3,vp='rotate')
      dev.off()
      
      next
      
      dm2 = dm2[order(dm2$clusterend2,dm2$cluster2,dm2$endI,dm2$end,dm2$begI,dm2$beg),]
      dm2$y = seq(1,dim(dm2)[1])
      myheight = 3000/400
      mywidth = 3000/400
      
      p5 = ggplot(dm2,aes(x=beg,y=y)) +
        geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster2))) +
        theme_bw() +
        annotate(geom="rect",xmin=641,xmax=841,ymin=0,ymax=dim(dm2)[1],lty=2,linewidth=1,alpha=0.15) +
        geom_point(aes(x=mid,y=y),pch=".") + ggtitle(paste(filename,"_",genes[genesInd],"_PEAKS",sep=""))
      
      
      dm2$midI = as.integer((dm2$beg + dm2$end)/2/400)*400
      dm2$mid = as.integer((dm2$beg + dm2$end)/2)
      dm2 = dm2[order(dm2$clusterend2,dm2$cluster2,dm2$midI,dm2$mid),]
      dm2$y = seq(1,dim(dm2)[1])
      myheight = 3000/400
      mywidth = 3000/400
      
      
      p6 = ggplot(dm2,aes(x=beg,y=y)) +
        geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster2))) +
        theme_bw() +
        annotate(geom="rect",xmin=641,xmax=841,ymin=0,ymax=dim(dm2)[1],lty=2,linewidth=1,alpha=0.15) +
        geom_point(aes(x=mid,y=y),pch=".") + ggtitle(paste(filename,"_",genes[genesInd],"_PEAKS",sep=""))
      
      #pdf(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_PEAKS.pdf",sep=""),height=myheight, width=mywidth)
      #print(k1)
      #print(k2)
      #print(k3)
      #print(k3c)
      #print(k4c)
      #print(k5c)
      #print(k6c)
      #print(k7c)
      #print(p1)
      #print(p2)
      #print(p2b)
      #print(p2c)
      #print(p3)
      #print(p4)
      #print(p5)
      #print(p6)
      #dev.off()
      
      #                       png(paste(filename,"_desc",treats[treatsInd],"_gene",genes[genesInd],"_PEAKS.png",sep=""),height=500, width=500)
      #                       print(p3)
      #                       dev.off()
      
      #pushViewport(viewport(name = "rotate", angle = 20, clip = "off", width = 0.7, height = 0.7))
      #print(p3, vp = "rotate"
      
    }
  }
}


