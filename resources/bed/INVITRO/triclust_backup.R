setwd("D:/cygwin64/home/mitochy/INVITRO")
library(plotly)
library(grid)
library(ggplot2)
library(gridExtra)

beds = data.frame()
if (file.exists("annotation.bed")) {
  beds = read.table("annotation.bed",sep="\t")
  colnames(beds)[1:6] = c("gene","beg","end","feature","val","strand")
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
  return(genes)
}
get_mytitle = function(genes, geneWant,mytreatWant,filename) {
  for (i in 1:length(genes)) {
    if (genes[i] == geneWant) {
      genesInd = i
    }
  }
  currgene = genes[genesInd]
  
  mybeds = get_mybeds(beds, currgene)
  
  
  mytitle = paste(filename,"_gene",genes[genesInd],"_desc",mytreatWant,sep="")
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
  colnames(df0)[1:6] = c("chr","beg","end","peakz","val","strand")
  if (dim(df0)[2] > 6) {
    colnames(df0) = c("chr","beg","end","peakz","val","strand","treat","VR","t0","t1")
    df0 = df0[grep("^(0|1|3|4)$",df0$t0),]
  }
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
get_dm3 = function(dm2, divby,threshold,useSD=0) {
  begtotal = c()
  endtotal = c()
  for (i in 1:(as.integer(5000/divby))) {
    begtotal[i] = 0
    endtotal[i] = 0
  }
  for (i in 1:dim(dm2)[1]) {
    begI = as.integer(dm2$beg[i]/divby)
    endI = as.integer(dm2$end[i]/divby)
    begtotal[begI] = begtotal[begI] + 1
    endtotal[endI] = endtotal[endI] + 1
  }
  dm2temp = data.frame(pos=as.integer(seq(1,5000/divby)),begtotal=begtotal,endtotal=endtotal)
  begthreshold = max(1,sd(dm2temp$begtotal))
  endthreshold = max(1,sd(dm2temp$endtotal))
  #plot(density(dm2temp$begtotal));segments(sd(dm2temp$begtotal),0,sd(dm2temp$begtotal),1);segments(sd(dm2temp$begtotal*2),0,sd(dm2temp$begtotal*2),1)
  #dm2temp[dm2temp$begtotal +dm2temp$endtotal > 0,]
  if (useSD == 1) {
    dm2temp[dm2temp$begtotal < begthreshold,]$begtotal = 0
    dm2temp[dm2temp$endtotal < endthreshold,]$endtotal = 0
    dm2temp=dm2temp[dm2temp$begtotal >= begthreshold | dm2temp$endtotal >= endthreshold,]
  } else {
    dm2temp = dm2temp[dm2temp$begtotal >= threshold | dm2temp$endtotal >= threshold,]
  }
  
  #plot(density(dm2temp$begtotal));segments(sd(dm2temp$begtotal),0,sd(dm2temp$begtotal),1,col='orange');segments(sd(dm2temp$begtotal*2),0,sd(dm2temp$begtotal*2),1,col = 'red4')
  # last_endtotal = 0
  # last_endpos = -1
  # last_begtotal = 0
  # last_begpos = -1

  dm3 = dm2temp
  for (i in seq(2,dim(dm3)[1],1)) {
    if (i == 2) {dm3 = dm2temp}
    
    #if (last_begpos == -1) {
    #  last_begtotal = dm3$begtotal[i]
    #  last_begpos = dm3$x[i]
    #} else {

    last_pos = dm3[i-1,]$pos
    curr_pos = dm3[i-0,]$pos
    last_begtotal = dm3[i-1,]$begtotal
    curr_begtotal = dm3[i-0,]$begtotal
    last_endtotal = dm3[i-1,]$endtotal
    curr_endtotal = dm3[i-0,]$endtotal
#    paste("prev: ",last_pos," (",last_begtotal," ",last_endtotal,") vs. curr: ",curr_pos," (",curr_begtotal," ",curr_endtotal,")",sep="")
    paste("prev: ",last_pos," (",last_begtotal,") vs. curr: ",curr_pos," (",curr_begtotal,")",sep="")

    if (last_pos >= curr_pos - 1 & last_pos < curr_pos + 1 & last_begtotal >= threshold & curr_begtotal >= threshold) {
      paste(i,last_begtotal+curr_begtotal)
      #cat(paste(" 1   dm3beg[last_pos]=",dm3$begtotal[last_pos],", last_begtotal=",last_begtotal," dm3begi=",dm3$begtotal[i],"\n",sep=""))
      #print(dm3[last_pos,])
      dm3[i-1,]$begtotal = last_begtotal + curr_begtotal
      dm3[i-0,]$begtotal = 0
      #cat(paste(" --> dm3beg[last_pos]=",dm3$begtotal[last_pos],", last_begtotal=",last_begtotal," dm3begi=",dm3$begtotal[i],"\n",sep=""))
      # } else {
      #   #last_begtotal = dm3$begtotal[i]
      #   #last_pos = i
      # }
    }
    cat(i)
    print(dm3[i,])
    cat("\n\n")
  }
  
  for (i in seq(dim(dm3)[1],1,-1)) {
    if (last_endpos == -1) {
      last_endtotal = dm3$endtotal[i]
      last_endpos = i
    } else {
      if (last_endtotal >= threshold & dm3$endtotal[i] >= threshold) {
        print(paste("last_endtotal=",last_endtotal," dm3endi=",dm3$endtotal[i],sep=""))
        last_endtotal = dm3$endtotal[i]
        dm3$endtotal[last_endpos] = dm3$endtotal[last_endpos] + dm3$endtotal[i]
        dm3$endtotal[i] = 0
      } else {
        last_endtotal = dm3$endtotal[i]
        last_endpos = i
      }
    }
  }
  
  dm3 = dm3[dm3$begtotal > threshold | dm3$endtotal > threshold,]
  dm3beg = subset(dm3[dm3$begtotal > threshold,],select=c("x","begtotal")); colnames(dm3beg) = c("pos","total")
  dm3beg$type = 'beg'
  dm3end = subset(dm3[dm3$endtotal > threshold,],select=c("x","endtotal")); colnames(dm3end) = c("pos","total")
  dm3end$type = 'end'
  dm3 = rbind(dm3beg,dm3end)
  
  initgraph(dm2,dm3,divby,genewant,threshold)
  return(dm3)
}

initgraph = function(dm2, dm3,divby,mygene="NA",threshold=0) {
  dm2$begD = as.integer(dm2$beg/divby)
  dm2$endD = as.integer(dm2$end/divby)
  dm2b = get_dm2b(dm2,divby,threshold)
  xMAX = as.integer(5000/divby)
  yMAX = as.integer(5000/divby)
  dm3beg = dm3[dm3$type == "beg",]
  dm3end = dm3[dm3$type == "end",]
  p = ggplot(dm2,aes(x=begD,y=endD)) +
    theme_bw() + theme(panel.grid=element_blank(),legend.position = "none") +
    annotate(geom="segment",x=0,y=0,xend=xMAX,yend=yMAX) +
    annotate(geom="segment",x=0,y=0,xend=xMAX,yend=0   ,color='grey',alpha=1,lty=1) +
    annotate(geom="segment",x=0,y=0,xend=0,yend=yMAX,color='grey',alpha=1,lty=1) +
    geom_segment(data=dm3beg,aes(x=pos,y=0,xend=pos,yend=-1*log(total+1,5)),color='red4') +
    geom_segment(data=dm3beg,aes(x=pos,y=0,xend=pos,yend=yMAX),color='red4',lty=2,alpha=0.2) +
    geom_segment(data=dm3end,aes(xend=0,y=pos,x=-1 * log(total+1,5),yend=pos),color='blue4') +
    geom_segment(data=dm3end,aes(x=0,y=pos,xend=xMAX,yend=pos),lty=2,col='blue4',alpha=0.2) +
    geom_text(data=dm3beg,aes(x=pos,y=1,label=pos),vjust=0,color='red4') +
    geom_text(data=dm3end,aes(x=1,y=pos,label=pos),hjust=0,color='blue4') +
    #geom_point(shape=15,size=0.25,color='grey') +
    geom_point(data=dm2b,aes(x=beg,y=end,size=sqrt(sum),alpha=1),color='blue4',shape=15) +
    geom_point(data=dm2b,aes(x=beg,y=end,alpha=1),color='purple',shape=15,size=3) +
    #geom_point(aes(x=beg/divby,y=end/divby),color='green4',pch=".",size=0.25,alpha=1) +
    scale_size_continuous(range=c(0.1,1)) +#,size=0.25,alpha=0.25) +
    coord_cartesian(xlim=c(-3,xMAX),ylim=c(-3,yMAX))
  print(p)
}


get_dm2b = function(dm2,divby,threshold=0) {
  temp2 = subset(dm2,select=c("beg","end"))
  temp2$beg = as.integer(temp2$beg / divby)
  temp2$end = as.integer(temp2$end / divby)
  temp2$count = 1
  temp3 = aggregate(temp2$count,by=list(temp2$beg,temp2$end),sum)
  colnames(temp3) = c("beg","end","sum")
  temp3 = temp3[temp3$sum >= threshold,]
  return(temp3)
}

get_dm2c = function(dm2,divby) {
  temp2 = rep(0,5000/divby)
  for (i in 1:dim(dm2)[1]) {
    for (j in as.integer(dm2$beg[i]/divby):as.integer(dm2$end[i]/divby)) {
      temp2[j] = temp2[j] + 1
    }
  }
  temp3 = data.frame(x=seq(length(temp2)),total=temp2)
  return(temp3)
}

dm2_to_matrix = function(dm2,divby) {
  temp2 = subset(dm2,select=c("beg","end"))
  temp2$beg = as.integer(temp2$beg / divby)
  temp2$end = as.integer(temp2$end / divby)
  temp2$count = 1
  temp3 = aggregate(temp2$count,by=list(temp2$beg,temp2$end),sum)
  colnames(temp3) = c("beg","end","sum")
  
  temp4 = matrix(ncol = max(temp3$end)+1,nrow=max(temp3$beg)+1)
  for (i in 1:dim(temp3)[1]) {
    curr.x = temp3$beg[i]
    curr.y = temp3$end[i]
    temp4[curr.x,curr.y] = temp3$sum[i]
  }
  return(temp4)
}

plot_ly_this = function(temp3) {
  p1 = plot_ly(temp3, x = ~-1*beg, y = ~-1*end, z = ~sum,type='scatter3d',size=~sum,
               marker = list(color = ~sum, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
  return(p1)
}


maxgenenumber = 99999


mydir = "ALL"
filesInd = 9

files = dir(mydir,".BED$")
files = paste(normalizePath(dirname(files)),mydir,files,sep="/")


genes = get_genes(files,filesInd)
genes = genes[order(genes)]
border0 = 0
border1 = 0
plus = "+"
if (border1 < 0) {plus = "-"}
#genewant = "T7_init_VR_10"
genewant = genes[1]
#genewant = "pFC53TVR_1_pair_0_T3TermMix"
#genewant = "pFC53_ApaLI"

dm1 = get_dm2(files, filesInd, genes,genewant)

mytreats = unique(dm1$treat)
mytreats =mytreats[order(mytreats)]
mytreatsInd = 9
mytreatwant = mytreats[mytreatsInd]
mytreatwant = 8

dm2 = dm1[dm1$treat == mytreatwant,]

mytitle = get_mytitle(genes,genewant,mytreatwant,basename(files[filesInd]))
mytitle = paste(mytitle,"\n(x & y border: ",-1 * border0,"/",plus,border1,")",sep="")

head(dm2)



dm4 = dm2[1:4]
dim(dm4)


#FUS: 50, 25, 25
divby = 50
mythres = as.integer(sqrt(dim(dm4)[1]/9230)*25)
mythres2 = as.integer(sqrt(dim(dm4)[1]/9230)*25)
mythres

if (length(grep("T7",mytitle,ignore.case=TRUE,perl=TRUE)) > 0) {
  #T7: 25, 100, 100
  divby = 25
  mythres = as.integer(sqrt(dim(dm4)[1]/9230)*100)
  mythres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
  mythres
} else if (length(grep("PFC53.+T3Term",mytitle,ignore.case = TRUE,perl=TRUE)) > 0) {
  #pFC53_T3Term: 15,max(50 etc)
  divby = 15
  mythres = max(50,as.integer(sqrt(dim(dm4)[1]/9230)*100))
  mythres2 = max(50,as.integer(sqrt(dim(dm4)[1]/9230)*100))
  mythres
} else if (length(grep("pFC53.+ApaLI",mytitle,ignore.case = TRUE,perl=TRUE)) > 0) {
  #pFC53_ApaLI: 25,as.integer
  divby = 25
  mythres = as.integer(sqrt(dim(dm4)[1]/9230)*100)
  mythres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
  mythres
} else if (length(grep("pFC9.+ApaLI",mytitle,ignore.case = TRUE,perl=TRUE)) > 0) {
  #pFC9_ApaLI: 25,as.integer
  divby = 25
  mythres = as.integer(sqrt(dim(dm4)[1]/9230)*100)
  mythres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
  mythres
} else if (length(grep("pFC9",mytitle,ignore.case = TRUE,perl=TRUE)) > 0) {
  #pFC9_sgRNA: 15,as.integer
  divby = 50
  mythres = as.integer(sqrt(dim(dm4)[1]/9230)*100)
  mythres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
  mythres
}
#divby=15
mythres=5
mythres2=5
#!-------------------
dm2b = get_dm2b(dm2,divby)
dm2c = get_dm2c(dm2,divby)
dm3 = get_dm3(dm2, divby, mythres)

#dm2b2$beg = dm2b2$beg * divby
#dm2b2$end = dm2b2$end * divby
#initgraph(dm2b2,dm3,divby,genewant)

myx = c()
myy = c()
mytype = c()

initgraph(dm2, dm3, divby,mytitle,mythres)
a = dm3[dm3$begtotal > 0,]
goodbegs = dm3[dm3$begtotal >= mythres,]$x
goodbegs = c(0,goodbegs,max(dm3$x))
for (i in seq(length(goodbegs),1,-1)) {
  xbeg0 = goodbegs[i-1]
  xbeg1 = goodbegs[i]
  #if (dim(dm3[dm3$x > xbeg0 & dm3$x <= xbeg1 & dm3$endtotal > mythres,])[1] == 0) {next}
  #end2 = dm3[dm3$x > xbeg0 & dm3$x <= xbeg1 & dm3$endtotal > mythres,]
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

goodends = dm3[dm3$endtotal >= mythres,]$x
goodends = c(as.integer(min(dm2$beg/divby)),goodends,max(dm3$x))

for (i in seq(length(goodends),1,-1)) {
  xend0 = goodends[i-1]
  xend1 = goodends[i]
  if (dim(dm4[as.integer(dm4$end/divby) > xend0 & as.integer(dm4$end/divby) <= xend1,])[1] == 0) {next}
  beg2 = dm4[as.integer(dm4$end/divby) > xend0 & as.integer(dm4$end/divby) <= xend1,]
  dm3b = get_dm3(beg2, divby, mythres2)
  dm3b = dm3b[dm3b$beg > 0,]
  mymin = as.integer(min(dm3[dm3$endtotal >= mythres2,]$x))
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


#-----------------!

finalbegpaste = c()
finalbegbackup = finalbeg

finalbeg = finalbegbackup
finalbeg2 = data.frame()
for (i in 1:dim(finalbeg)[1]) {
  x0 = finalbeg$x0[i]
  y0 = finalbeg$y0[i]
  x1 = finalbeg$x1[i]
  y1 = finalbeg$y1[i]
  if (y0 > y1) {
    y0 = finalbeg$y1[i]
    y1 = finalbeg$y0[i]
    finalbeg$y0[i] = y0
    finalbeg$y1[i] = y1
  }
  mypaste = paste(x0,y0,x1,y1)
  if (length(finalbegpaste[finalbegpaste == mypaste]) == 0) {
    finalbegpaste = c(finalbegpaste,mypaste)
    finalbeg2 = rbind(finalbeg2,data.frame(x0=x0,y0=y0,x1=x1,y1=y1))
  }
}
finalbeg = finalbeg2
finalbeg$y = seq(1,dim(finalbeg)[1])*10
finalbeg$cluster = seq(1,dim(finalbeg)[1])
finalend$y = seq(1,dim(finalend)[1])*10




# p2 = ggplot(dm2b,aes(beg,end)) +
#    geom_point(aes(alpha=sum,size=sum)) +
#    scale_alpha_continuous(breaks=c(0,5,10,20,40,80,160,320,9999999),range=c(0.1,1)) +#,size=0.05) +#sqrt(sum)/10)) +
#    scale_size_continuous(breaks = c(0,5,10,20,40,80,160,320,9999999),range = c(0.001,1)) +#values=c(0.01,0.05,1)) +
# #   scale_size_continuous(breaks = c(0,1,2,4,8,16,999),range = c(0.01,2)) +#values=c(0.01,0.05,1)) +
#    theme_bw() + coord_cartesian(xlim=c(0,100),ylim=c(0,100)) +    ggtitle(mytitle) +
#    annotate(geom="segment",x=0,y=0,xend=100,yend=100)
# p3 = p2 +
#    geom_segment(data=finalbeg,aes(x=x0,xend=x1,y=y0,yend=y1)) +
#    geom_segment(data=finalend,aes(x=x0,xend=x1,y=y0,yend=y1)) +
#    geom_point(data=finalbeg,aes(x=x0,y=y0),color = "red2",shape=15) +
#    geom_point(data=finalbeg,aes(x=x1,y=y1),color = "red4",shape=18) +
#    geom_point(data=finalend,aes(x=x0,y=y0),color="blue2",shape=15) +
#    geom_point(data=finalend,aes(x=x1,y=y1),color="blue4",shape=18)

p2 = ggplot(dm2b,aes(beg,end)) +
  geom_point(aes(alpha=sqrt(sum),size=sqrt(sum))) +
  geom_point(pch=".",color="grey") +
  scale_size_continuous(range=c(0.1,1)) +
  scale_alpha_continuous(range=c(0.1,1)) +
  #scale_alpha_continuous(breaks=c(0,5,10,20,40,80,160,320,9999999)/100,range=c(0.1,1)) +#,size=0.05) +#sqrt(sum)/10)) +
  #scale_size_continuous(breaks = c(0,5,10,20,40,80,160,320,9999999)/100,range = c(0.001,1)) +#values=c(0.01,0.05,1)) +
  #   scale_size_continuous(breaks = c(0,1,2,4,8,16,999),range = c(0.01,2)) +#values=c(0.01,0.05,1)) +
  theme_bw() + coord_cartesian(xlim=c(0,5000/divby),ylim=c(0,5000/divby)) +    ggtitle(mytitle) +
  annotate(geom="segment",x=0,y=0,xend=5000/divby,yend=5000/divby) + theme(legend.position = "top")
p3 = p2 +
  geom_segment(data=finalbeg,aes(x=x0,xend=x1,y=y0,yend=y1)) +
  #   geom_segment(data=finalend,aes(x=x0,xend=x1,y=y0,yend=y1)) +
  geom_point(data=finalbeg,aes(x=x0,y=y0),color = "red2",shape=15) +
  geom_point(data=finalbeg,aes(x=x1,y=y1),color = "red4",shape=18)
#   geom_point(data=finalend,aes(x=x0,y=y0),color="blue2",shape=15) +
#   geom_point(data=finalend,aes(x=x1,y=y1),color="blue4",shape=18)

# p2
# p3
#pdf("test.pdf")
p4 = ggplot(dm2c,aes(x,total)) +
  geom_line() +
  theme_bw() + coord_cartesian(xlim=c(0,5000/divby),ylim=c(0,max(dm2c$total))) +    ggtitle(mytitle) +
  geom_segment(data=finalbeg,aes(x=x0,y=0,xend=x0,yend=max(dm2c$total)),lty=2,color='red4')
#  geom_segment(data=finalend[!finalend$y1 %in% finalbeg$x0,],aes(x=y1,y=0,xend=y1,yend=max(dm2c$total)),lty=2,color='blue4')


dm2$begD = as.integer(dm2$beg/divby)
dm2$endD = as.integer(dm2$end/divby)
dm2 = dm2[order(dm2$begD, dm2$endD),]
dm2$cluster = 0
for (i in 1:dim(dm2)[1]) {
  begD = dm2[i,]$begD
  endD = dm2[i,]$endD
  for (j in 1:dim(finalbeg)[1]) {
    x0 = finalbeg[j,]$x0
    y0 = finalbeg[j,]$y0
    x1 = finalbeg[j,]$x1
    y1 = finalbeg[j,]$y1
    #if (begD >= x0 - 1 & begD <= x0 + 1 & endD >= y0 & endD <= y1) {
    #if (begD >= x0 - 1 & endD >= y0 & endD <= y1) {
    # if (begD >= x0 - 1 & begD <= y0 & endD >= y0 & endD <= y1) {
    #    if (dm2[i,]$cluster == 0) {
    #       dm2[i,]$cluster = j
    #       print(paste(i,begD,endD,j))
    #       break
    #    }
    # }
    #if (begD >= x0 - 1 & begD <= x0 + 1 & endD >= y1 & endD <= y0) {
    #if (begD >= x0 - 1 & endD >= y1 & endD <= y0) {
    if (dm2[i,]$cluster != 0) {next}
    if (begD >= x0 - 1 & begD <= y0 & endD >= y0 & endD <= y1) {
      dm2[i,]$cluster = j
      print(paste(i,begD,endD,j))
      break
    }
  }
}
for (i in 1:dim(dm2)[1]) {
  begD = dm2[i,]$begD
  endD = dm2[i,]$endD
  for (j in 1:dim(finalbeg)[1]) {
    x0 = finalbeg[j,]$x0
    y0 = finalbeg[j,]$y0
    x1 = finalbeg[j,]$x1
    y1 = finalbeg[j,]$y1
    #if (begD >= x0 - 1 & begD <= x0 + 1 & endD >= y0 & endD <= y1) {
    #if (begD >= x0 - 1 & endD >= y0 & endD <= y1) {
    # if (begD >= x0 - 1 & begD <= y0 & endD >= y0 & endD <= y1) {
    #    if (dm2[i,]$cluster == 0) {
    #       dm2[i,]$cluster = j
    #       print(paste(i,begD,endD,j))
    #       break
    #    }
    # }
    #if (begD >= x0 - 1 & begD <= x0 + 1 & endD >= y1 & endD <= y0) {
    #if (begD >= x0 - 1 & endD >= y1 & endD <= y0) {
    if (dm2[i,]$cluster != 0) {next}
    if (begD >= x0 - 1 & begD <= y1 & endD >= y0 & endD <= y1) {
      dm2[i,]$cluster = j
      print(paste(i,begD,endD,j))
      break
    }
  }
}
#finalbeg=finalbeg[finalbeg$cluster <=8,]
dm2 = dm2[order(dm2$cluster,dm2$begD,dm2$endD),]
dm2$y = seq(1,dim(dm2)[1])
p5 = ggplot(dm2[dm2$cluster!= -1,],aes(begD,endD)) +
  geom_segment(aes(x=begD,xend=endD,y=y,yend=y,color=as.factor(cluster))) +
  theme_bw() + coord_cartesian(xlim=c(0,5000/divby),ylim=c(0,dim(dm2)[1])) +    ggtitle(mytitle) +
  theme(legend.position = "top")
p6 = p5 +
  geom_point(data=finalbeg,aes(x=x0,y=y),color='purple') +
  geom_point(data=finalbeg,aes(x=y0,y=y),color='red4') +
  geom_point(data=finalbeg,aes(x=y1,y=y),color='blue4') +
  geom_segment(data=finalbeg,aes(x=x0,y=y,xend=y0,yend=y),color='red4') +
  geom_segment(data=finalbeg,aes(x=y0,y=y,xend=y1,yend=y),color='blue4') +
  geom_segment(data=finalbeg,aes(x=x0,y=0,xend=x0,yend=dim(dm2)[1]/3),lty=2,lwd=0.5,alpha=0.25,color='purple')+
  geom_segment(data=finalbeg,aes(x=y0,y=dim(dm2)[1]/3,xend=y0,yend=dim(dm2)[1]/3*2),lty=2,lwd=0.5,color='red4')+
  geom_segment(data=finalbeg,aes(x=y1,y=dim(dm2)[1]/3*2,xend=y1,yend=dim(dm2)[1]),lty=2,lwd=0.5,color='blue4')
#   geom_segment(data=finalend[!finalend$y1 %in% finalbeg$x0,],aes(x=y1,y=0,xend=y1,yend=dim(dm2)[1]/2),lwd=0.5,lty=2,color='blue4')

pdf("testpdf.pdf",height=30/2,width=20/2)
grid.arrange(p2,p3,p4,p4,p5,p6,ncol=2,nrow=3,heights=c(2,1,2))
dev.off()#p5
#   geom_segment(data=finalbeg,aes(x=x0,y=dim(dm2)[1]/2,xend=x0,yend=dim(dm2)[1]),lty=2,lwd=0.5,color='red4')+
#plot_ly_this(dm2b)

# Figure 1: QE VR diagram update, init & term series
# Figure 2: VR frequency
# Figure 2b: Distance from TSS initiation, RLoop length
# Figure 3:  something footprint and clusters, triangles
# Figure 4: SSB
# Figure 5: sgRNA
