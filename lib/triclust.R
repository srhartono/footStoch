setwd("D:/cygwin64/home/mitochy/Work/Project/Ethan/footStoch/")
library(plotly)
library(grid)
library(ggplot2)
library(gridExtra)

get_genes = function(file) {
  filename = basename(file)
  filename = gsub("^(.+).BED$","\\1",filename,perl=T)
  
  df = read.table(file,sep="\t")
  colnames(df)[1:6] = c("chr","beg","end","peakz","val","strand")
  if (dim(df)[2] > 6) {
    colnames(df) = c("chr","beg","end","peakz","val","strand","treat","VR","t0","t1")
    df = df[grep("^(0|1|3|4)$",df$t0),]
  }
  df = df[order(df$beg, df$end),]
  
  genes = unique(df[,1])
  return(genes)
}

get_mytitle2 = function(file,gene,treat) {
  bed = BEDS[BEDS$chr == gene,]
  filename = basename(file)
  filename = gsub("^(.+).BED$","\\1",filename,perl=T)
  mytitle = paste(filename,"_gene",gene,"_desc",treat,sep="")
  return(mytitle)
}

get_dm2 = function(file, gene) {
  filename = basename(file)
  filename = gsub("^(.+).BED$","\\1",filename,perl=T)

  df = read.table(file,sep="\t")
  colnames(df)[1:6] = c("chr","beg","end","peakz","val","strand")
  if (dim(df)[2] > 6) {
    colnames(df) = c("chr","beg","end","peakz","val","strand","treat","VR","t0","t1")
    df = df[grep("^(0|1|3|4)$",df$t0),]
  }
  df = df[order(df$beg, df$end),]
  if (!'file' %in% colnames(df)) {
    df$file = "DUMMY"
  }

  if (defined(BEDS[BEDS$chr == gene,])) {
    bed = BEDS[BEDS$chr == curr.gene,]
  }

  if (defined(df[df$chr == gene,])) {
    df = df[df$chr == gene,]
  }

  # if (dim(df)[1] < 8) {
  #   return(df)
  # }
  
  # if (dim(df)[1] > maxgenenumber) {
  #   df = df[sample(seq(1,dim(df)[1]),maxgenenumber,replace=F),]
  # }
  
  rownames(df) = paste(seq(1,dim(df)[1]))
  
  dm = df;
  dm$cluster = 1#myclust[,1]
  dm$clusterbeg = 1
  dm$clusterend = 1
  
  dm$mid = ai((dm$beg + dm$end)/2/50)*50
  dm$begI = ai(dm$beg/2/50)*50
  dm$endI = ai(dm$end/2/50)*50
  dm$lenI = dm$endI - dm$begI
  
  dm$midI = ai((dm$beg + dm$end)/2/400)*400
  dm$mid2 = ai((dm$beg + dm$end)/2)
  dm = dm[order(dm$clusterbeg, dm$cluster,dm$begI,dm$beg,dm$endI,dm$end ),]
  dm = dm[order(dm$clusterend,dm$cluster,dm$midI,dm$mid2),]
  dm$y = seq(1,dim(dm)[1])

  dm$cluster2 = 1
  
  dm = dm[order(dm$begI,dm$beg,dm$endI,dm$end),]
  dm$y = seq(1,dim(dm)[1])
  return(dm)
}

get_dm3 = function(dm2, divby,threshold,useSD=0) {
#  dm2temp = get_dm2b(dm2,divby,threshold)
  dm2 = dm2[order(dm2$beg,dm2$end),]
  begtotal = c()
  endtotal = c()
  for (i in 1:(ai(5000/divby))) {
    begtotal[i] = 0
    endtotal[i] = 0
  }
  for (i in 1:dim(dm2)[1]) {
    begI = ai(dm2$beg[i]/divby)
    endI = ai(dm2$end[i]/divby)
    begtotal[begI] = begtotal[begI] + 1
    endtotal[endI] = endtotal[endI] + 1
  }
  dm2temp = data.frame(pos=ai(seq(1,5000/divby)),begtotal=begtotal,endtotal=endtotal)
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
  for (i in seq(1,dim(dm3)[1],1)) {
    if (i == 1) {dm3$begtotal = dm2temp$begtotal;  best_i = -1; last_pos = -1; last_begtotal = -1;}
    curr_pos        = dm3$pos[i]
    curr_begtotal   = dm3$begtotal[i]
    if (curr_begtotal >= threshold & best_i == -1) {
      best_i = i
      last_pos      = dm3$pos[best_i]
      last_begtotal = dm3$begtotal[best_i]
    } else if (curr_begtotal < threshold) {
      best_i = -1
      last_pos = -1
      last_begtotal = -1
    } else if (curr_begtotal >= threshold & best_i != -1 & last_pos < curr_pos & last_pos >= curr_pos - 2) {
      cat("i=",i,"currpos=",curr_pos,",best_i=",best_i,"bset_pos=",last_pos,"sum=",last_begtotal+curr_begtotal,"\n")
      dm3$begtotal[best_i] = dm3$begtotal[best_i] + curr_begtotal
      dm3$begtotal[i] = 0
      best_i = best_i
      last_pos      = dm3$pos[i]
      last_begtotal = dm3$begtotal[best_i]
    } else if (curr_begtotal >= threshold & best_i != -1 & last_pos < curr_pos - 2) {
      best_i = i
      last_pos      = dm3$pos[best_i]
      last_begtotal = dm3$begtotal[best_i]
    }
  }
  cbind(dm2temp,dm3)

  for (i in seq(1,dim(dm3)[1],1)) {
    if (i == 1) {dm3$endtotal = dm2temp$endtotal;  best_i = -1; last_pos = -1; last_endtotal = -1}
    curr_pos        = dm3$pos[i]
    curr_endtotal   = dm3$endtotal[i]
    if (curr_endtotal >= threshold & best_i == -1) {
      best_i = i
      last_pos      = dm3$pos[best_i]
      last_endtotal = dm3$endtotal[best_i]
    } else if (curr_endtotal < threshold) {
      best_i = -1
      last_pos = -1
      last_endtotal = -1
    } else if (curr_endtotal >= threshold & best_i != -1 & last_pos < curr_pos & last_pos >= curr_pos - 2) {
      cat("i=",i,"currpos=",curr_pos,",best_i=",best_i,"bset_pos=",last_pos,"sum=",last_endtotal+curr_endtotal,"\n")
      dm3$endtotal[best_i] = dm3$endtotal[best_i] + curr_endtotal
      dm3$endtotal[i] = 0
      best_i = best_i
      last_pos      = dm3$pos[i]
      last_endtotal = dm3$endtotal[best_i]
    } else if (curr_endtotal >= threshold & best_i != -1 & last_pos < curr_pos - 2) {
      best_i = i
      last_pos      = dm3$pos[best_i]
      last_endtotal = dm3$endtotal[best_i]
    }
  }
  cbind(dm2temp,dm3)
  dm3[dm3$begtotal < threshold,]$begtotal = 0
  dm3[dm3$endtotal < threshold,]$endtotal = 0
  cbind(dm2temp,dm3)
  # 
  # 
  # for (i in seq(dim(dm3)[1],1,-1)) {
  #   if (last_endpos == -1) {
  #     last_endtotal = dm3$endtotal[i]
  #     last_endpos = i
  #   } else {
  #     if (last_endtotal >= threshold & dm3$endtotal[i] >= threshold) {
  #       print(paste("last_endtotal=",last_endtotal," dm3endi=",dm3$endtotal[i],sep=""))
  #       last_endtotal = dm3$endtotal[i]
  #       dm3$endtotal[last_endpos] = dm3$endtotal[last_endpos] + dm3$endtotal[i]
  #       dm3$endtotal[i] = 0
  #     } else {
  #       last_endtotal = dm3$endtotal[i]
  #       last_endpos = i
  #     }
  #   }
  # }
  # 
  
  dm3 = dm3[dm3$begtotal > threshold | dm3$endtotal > threshold,]
  dm3beg = subset(dm3[dm3$begtotal > threshold,],select=c("pos","begtotal")); colnames(dm3beg) = c("pos","total")
  dm3beg$type = 'beg'
  dm3end = subset(dm3[dm3$endtotal > threshold,],select=c("pos","endtotal")); colnames(dm3end) = c("pos","total")
  dm3end$type = 'end'
  dm3 = rbind(dm3beg,dm3end)
  
  #initgraph(dm2,dm3,divby,genewant,threshold)
  return(dm3)
}

initgraph = function(dm2, dm3,divby,mygene="NA",threshold=0) {
  dm2$begD = ai(dm2$beg/divby)
  dm2$endD = ai(dm2$end/divby)
  dm2b = get_dm2b(dm2,divby,threshold)
  xMAX = ai(5000/divby)
  yMAX = ai(5000/divby)
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
    geom_point(shape=15,size=1,color='black') +
    #geom_point(data=dm2b,aes(x=beg,y=end,size=sqrt(sum),alpha=1),color='blue4',shape=15) +
    geom_point(data=dm2b,aes(x=beg,y=end,alpha=1),color='purple',shape=15,size=1) +
    #geom_point(aes(x=beg/divby,y=end/divby),color='green4',pch=".",size=0.25,alpha=1) +
    scale_size_continuous(range=c(0.1,1)) +#,size=0.25,alpha=0.25) +
    coord_cartesian(xlim=c(-3,xMAX),ylim=c(-3,yMAX))
  return(p)
}

initgraph2 = function(p,goodbegs,goodends) {
  yMIN = 0
  yMAX = 5000/divby
  xMIN = 0
  xMAX = 5000/divby
  p = p + annotate(geom='segment',x=goodbegs[1],y=0,xend=goodbegs[1],yend=yMAX,lty=2,alpha=0.2)
  p = p + annotate(geom='segment',x=goodbegs[1],y=0,xend=goodbegs[1],yend=-2,lty=1,color='red2')
  p = p + annotate(geom='text',x=goodbegs[1],y=1,label=goodbegs[1],vjust=0,color='red2')
  p = p + annotate(geom='segment',x=goodbegs[length(goodbegs)],y=0,xend=goodbegs[length(goodbegs)],yend=yMAX,lty=2,alpha=0.2)
  p = p + annotate(geom='segment',x=goodbegs[length(goodbegs)],y=0,xend=goodbegs[length(goodbegs)],yend=-2,lty=1,color='red2')
  p = p + annotate(geom='text',x=goodbegs[length(goodbegs)],y=1,label=goodbegs[length(goodbegs)],vjust=0,color='red2')
  
  p = p + annotate(geom='segment',x=0,y=goodends[length(goodends)],yend=goodends[length(goodends)],xend=xMAX,lty=2,alpha=0.2)
  p = p + annotate(geom='segment',x=0,y=goodends[1],yend=goodends[1],xend=xMAX,lty=2,alpha=0.2)
  p = p + annotate(geom='segment',y=goodends[length(goodends)],x=0,yend=goodends[length(goodends)],xend=-2,lty=1,color='blue2')
  p = p + annotate(geom='text',y=goodends[length(goodends)],x=1,label=goodends[length(goodends)],hjust=0,color='blue2')
  p = p + annotate(geom='segment',y=goodends[1],x=0,yend=goodends[1],xend=-2,lty=1,color='blue2')
  p = p + annotate(geom='text',y=goodends[1],x=1,label=goodends[1],hjust=0,color='blue2')
  print(p)
}


get_dm2b = function(dm2,divby,threshold=0) {
  temp2 = subset(dm2,select=c("beg","end"))
  temp2$beg = ai(temp2$beg / divby)
  temp2$end = ai(temp2$end / divby)
  temp2$count = 1
  temp3 = aggregate(temp2$count,by=list(temp2$beg,temp2$end),sum)
  colnames(temp3) = c("beg","end","sum")
  temp3 = temp3[temp3$sum >= threshold,]
  return(temp3)
}

get_dm2c = function(dm2,divby) {
  temp2 = rep(0,5000/divby)
  for (i in 1:dim(dm2)[1]) {
    for (j in ai(dm2$beg[i]/divby):ai(dm2$end[i]/divby)) {
      temp2[j] = temp2[j] + 1
    }
  }
  temp3 = data.frame(x=seq(length(temp2)),total=temp2)
  return(temp3)
}

dm2_to_matrix = function(dm2,divby) {
  temp2 = subset(dm2,select=c("beg","end"))
  temp2$beg = ai(temp2$beg / divby)
  temp2$end = ai(temp2$end / divby)
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

mydir = './resources/bed/INVITRO/ALL/'
files = dir(mydir,".BED$")
files = paste(normalizePath(dirname(files)),mydir,files,sep="/")
#filesInd = 9

total_todo = 0
final3 = data.frame()
files2 = c()
genes2 = c()
treat2 = c()

for (filesInd in seq(1,size(files))) {
  if (filesInd == 1) {
    total_todo = 0
    files2 = c()
    genes2 = c()
    treat2 = c()
  }
  curr.file = files[filesInd]
  if (grepl('_TEMP',curr.file)) {
    next
  }

  genes = get_genes(curr.file)
  print(paste(filesInd,curr.file))
  #realgene = gsub("^(.+)_(PEAK_C|PEAK_TEMP_C).BED","\\1",basename(files[filesInd]),perl=T)
  
  genes = genes[order(genes)]

  for (genesInd in seq(1,size(genes))) {
    curr.gene = genes[genesInd]
    dm1 = get_dm2(files[filesInd], curr.gene)
    if (size(dm1) == 0) {
      next
    }
    if (dim(dm1[is.na(dm1$treat),])[1] > 0) {
      dm1[is.na(dm1$treat),]$treat = "NONE"
    }
    if (defined(dm1[dm1$treat == 'NA',])) {
      dm1[is.na(dm1$treat),]$treat = "NONE"
    }

    mytreats = unique(dm1$treat)
    mytreats = mytreats[order(mytreats)]
    for (mytreatsInd in seq(1,size(mytreats))) {
      curr.treat = mytreats[mytreatsInd]
      dm2 = dm1[dm1$treat == curr.treat,]
      mytitle = get_mytitle2(curr.file,curr.gene,curr.treat)
      total_todo = total_todo + 1

      files2 = c(files2,curr.file)
      genes2 = c(genes2,curr.gene)
      treat2 = c(treat2,curr.treat)
      print(paste('      ',curr.gene,mytreatsInd,curr.treat,size(dm2),length(treat2),total_todo))
      
    }
  }
}
total_todo

current_todo = 0
last.files = -1
last.genes = -1
last.treat = -1

filesInd = 1
genesInd = 1
treatInd = 1
lastdone = 1
dm1 = data.frame()
dm2 = data.frame()
loopend = size(files2)
for (currInd in seq(lastdone,loopend)) {
  lastdone = currInd
  
  if (currInd == 1) {
    current_todo = 1
    final3 = data.frame()
  } else {
    current_todo = currInd
  }

  curr.files = files2[currInd]
  curr.genes = genes2[currInd]
  curr.treat = treat2[currInd]
  mytitle = get_mytitle2(curr.files,curr.genes,curr.treat)
  
  if (last.files != curr.files) {
    dm1 = get_dm2(curr.files,curr.genes)
    dm2 = dm1[dm1$treat == curr.treat,]
  } else if (last.genes != curr.genes) {
    dm1 = get_dm2(curr.files,curr.genes)
    dm2 = dm1[dm1$treat == curr.treat,]
  } else if (last.treat != curr.treat) {
    dm2 = dm1[dm1$treat == curr.treat,]
  }
  
  final3.temp = main(dm1,dm2)
  if (defined(final3.temp)) {
    final3.temp$treat = curr.treat
    final3.temp$gene = curr.genes
    final3.temp$file = curr.genes
    final3 = rbind(final3,final3.temp)
  } else {
    print(paste(filesInd,basename(curr.files),genesInd,curr.genes,treatInd,curr.treat,'has no cluster!'))
  }

  
  if (last.files != curr.files) {cat(filesInd,'. ',basename(curr.files),'\n',sep=''); filesInd = filesInd + 1; genesInd = 1; treatInd = 1}
  if (last.genes != curr.genes) {cat('  ',genesInd,'. ',curr.genes,' (n=',size(dm1),')','\n',sep=''); genesInd = genesInd +  1; treatInd = 1}
  if (last.treat != curr.treat) {cat('    ',treatInd,'. ',curr.treat, ' (n=',size(dm2),') (',current_todo,'/',total_todo,')','\n',sep=''); treatInd = treatInd + 1}
  last.files = curr.files
  last.genes = curr.genes
  last.treat = curr.treat
}
# 
# final3 = data.frame()
# for (filesInd in seq(1,size(files))) {
#   curr.file = files[filesInd]
#   if (grepl('_TEMP',curr.file)) {
#     next
#   }
#   
#   genes = get_genes(files,filesInd)
#   print(paste(filesInd,curr.file))
#   realgene = gsub("^(.+)_(PEAK_C|PEAK_TEMP_C).BED","\\1",basename(files[filesInd]),perl=T)
#   
#   genes = genes[order(genes)]
#   print(paste('-',genes))
#   
#   for (genesInd in seq(1,size(genes))) {
#     curr.gene = genes[genesInd]
#     dm1 = get_dm2(files, filesInd, genes,curr.gene)
#     print(paste('   ',genesInd,curr.gene,size(dm1)))
#     
#     if (size(unique(dm1$treat)) == 0) {
#       dm2 = dm1
#       mytitle = get_mytitle2(genes,curr.gene,'NA',basename(files[filesInd]))
#       print(paste('      ',mytreatsInd,'NA',size(dm2)))
#       current_todo = current_todo + 1
#       # final3.temp = main(dm1,dm2)
#       # final3.temp$treat = curr.treat
#       # final3.temp$gene = curr.gene
#       # final3.temp$file = curr.gene
#       # final3 = rbind(final3,final3.temp)
#       if (current_todo %% 100 == 0) {
#         print(paste('Done',current_todo,'/',total_todo))
#       }
#       
#     } else {
#       mytreats = unique(dm1$treat)
#       mytreats = mytreats[order(mytreats)]
#       for (mytreatsInd in seq(1,size(mytreats))) {
#         curr.treat = mytreats[mytreatsInd]
#         dm2 = dm1[dm1$treat == curr.treat,]
#         mytitle = get_mytitle2(genes,curr.gene,curr.treat,basename(files[filesInd]))
#         print(paste('      ',mytreatsInd,curr.treat,size(dm2)))
#         current_todo = current_todo + 1
#         if (current_todo %% 100 == 0) {
#           print(paste('Done',current_todo,'/',total_todo))
#         }
#       #   
#       #   final3.temp = main(dm1,dm2)
#       #   final3.temp$treat = curr.treat
#       #   final3.temp$gene = curr.gene
#       #   final3.temp$file = curr.gene
#       #   final3 = rbind(final3,final3.temp)
#       }
#     }
#   }
# }

main = function(dm1, dm2) {
  dm4 = dm2[1:4]
  dim(dm4)
  
  
  #FUS: 50, 25, 25
  divby = 50
  mythres = ai(sqrt(dim(dm4)[1]/9230)*25)
  mythres2 = ai(sqrt(dim(dm4)[1]/9230)*25)
  mythres
  
  if (length(grep("T7",mytitle,ignore.case=TRUE,perl=TRUE)) > 0) {
    #T7: 25, 100, 100
    divby = 25
    mythres = ai(sqrt(dim(dm4)[1]/9230)*100)
    mythres2 = ai(sqrt(dim(dm4)[1]/9230)*100)
    mythres
  } else if (length(grep("PFC53.+T3Term",mytitle,ignore.case = TRUE,perl=TRUE)) > 0) {
    #pFC53_T3Term: 15,max(50 etc)
    divby = 15
    mythres = max(50,ai(sqrt(dim(dm4)[1]/9230)*100))
    mythres2 = max(50,ai(sqrt(dim(dm4)[1]/9230)*100))
    mythres
  } else if (length(grep("pFC53.+ApaLI",mytitle,ignore.case = TRUE,perl=TRUE)) > 0) {
    #pFC53_ApaLI: 25,ai
    divby = 25
    mythres = ai(sqrt(dim(dm4)[1]/9230)*100)
    mythres2 = ai(sqrt(dim(dm4)[1]/9230)*100)
    mythres
  } else if (length(grep("pFC9.+ApaLI",mytitle,ignore.case = TRUE,perl=TRUE)) > 0) {
    #pFC9_ApaLI: 25,ai
    divby = 25
    mythres = ai(sqrt(dim(dm4)[1]/9230)*100)
    mythres2 = ai(sqrt(dim(dm4)[1]/9230)*100)
    mythres
  } else if (length(grep("pFC9",mytitle,ignore.case = TRUE,perl=TRUE)) > 0) {
    #pFC9_sgRNA: 15,ai
    divby = 50
    mythres = ai(sqrt(dim(dm4)[1]/9230)*100)
    mythres2 = ai(sqrt(dim(dm4)[1]/9230)*100)
    mythres
  }
  
  #!-------------------
  finalbeg = data.frame()
  finalend = data.frame()
  final0 = data.frame()
  
  #!-------------------
  dm2b = get_dm2b(dm2,divby)
  dm2c = get_dm2c(dm2,divby)
  dm3 = get_dm3(dm2, divby, mythres)
  
  yMIN = 0
  yMAX = 5000/divby
  xMIN = 0
  xMAX = 5000/divby
  myx = c()
  myy = c()
  mytype = c()
  
  goodbegs = dm3[dm3$type == 'beg',]$pos
  goodbegs = c(min(dm2b$beg)-1,goodbegs,max(dm2b$beg))
  goodbegs
  goodends = dm3[dm3$type == 'end',]$pos
  goodends = c(min(dm2b$end)-1,goodends,max(dm2b$end))
  goodends
  p = initgraph(dm2, dm3, divby,mytitle,mythres)
  p2 = p
  #initgraph2(p,goodbegs,goodends)
  for (i in seq(length(goodbegs),1,-1)) {
    # if (i == length(goodbegs)) {
    #   final0 = data.frame()
    #   finalbeg = data.frame()
    #   finalend = data.frame()
    # }
    xbeg0 = goodbegs[i-1]
    xbeg1 = goodbegs[i]
    if (dim(dm4[ai(dm4$beg/divby) > xbeg0 & ai(dm4$beg/divby) <= xbeg1,])[1] == 0) {next}
    end2 = dm4[ai(dm4$beg/divby) > xbeg0 & ai(dm4$beg/divby) <= xbeg1,]
    dm3b = get_dm3(end2, divby, mythres2)
    dm3b = dm3b[dm3b$type == 'end' &  dm3b$pos > 0,]
    print(paste(xbeg0,xbeg1,dm3b[,1], ":", dm3b[,3]))
    dm3c = dm3b
    if (dim(dm3c[dm3c$type == 'end' & dm3c$pos  >= mythres2,])[1] == 0) {next}
    if (dim(dm3c[dm3c$type == 'end' & dm3c$pos  >= mythres2,])[1] > 0) {dm3c = dm3c[dm3c$type == 'end' & dm3c$pos  >= mythres2,]}
    for (j in 1:length(dm3c$pos)) {
      curr.y = dm3c$pos[j]
      curr.x0end = xbeg0
      curr.x1beg = xbeg1
      final.temp.beg0 = data.frame(x=curr.x0end,y=curr.y,type='beg0')
      final.temp.beg1 = data.frame(x=curr.x1beg,y=curr.y,type='beg1')
      final.temp = rbind(final.temp.beg0,final.temp.beg1)
      # final0 = rbind(final0,final.temp)
      # 
      # p2 = p2 + 
      #   geom_point(data=final.temp.beg0,aes(x=x,y=y),color=rgb(0,0,1,1),pch=15) +
      #   geom_line( data=final.temp,     aes(x=x,y=y),color=rgb(0,0,1,1)) +
      #   geom_point(data=final.temp.beg1,aes(x=x,y=y),color=rgb(0,0,1,1),pch=15)
      # myx = c(myx,curr.x0end,curr.x1beg)
      # myy = c(myy,curr.y,curr.y)
      # mytype = c(mytype,'beg0','beg1')
      finalbeg = rbind(finalbeg,data.frame(x0end=curr.x0end,y0end=curr.y,x1beg=curr.x1beg,y1beg=curr.y))
    }
  }
  #initgraph2(p2,goodbegs,goodends)
  
  # p2.save = p2
  # myx.save = myx
  # myy.save = myy
  # mytype.save = mytype
  for (i in seq(length(goodends),1,-1)) {
    # if (i == length(goodends)) {
    #   p = p.save
    #   mytype = mytype.save
    #   myx = myx.save
    #   myy = myy.save
    #   finalend=  data.frame()
    # }
    xend0 = goodends[i-1]
    xend1 = goodends[i]
    if (dim(dm4[ai(dm4$end/divby) > xend0 & ai(dm4$end/divby) <= xend1,])[1] == 0) {next}
    beg2 = dm4[ai(dm4$end/divby) > xend0 & ai(dm4$end/divby) <= xend1,]
    dm3b = get_dm3(beg2, divby, mythres2)
    dm3b = dm3b[dm3b$type == 'beg' &  dm3b$pos > 0,]
    dm3c = dm3b
    if (dim(dm3c[dm3c$type == 'beg' & dm3c$pos  >= mythres2,])[1] == 0) {next}
    if (dim(dm3c[dm3c$type == 'beg' & dm3c$pos  >= mythres2,])[1] > 0) {dm3c = dm3c[dm3c$type == 'beg' & dm3c$pos  >= mythres2,]}
    for (j in seq(1:length(dm3c$pos))) {
      cat('\n')
      curr.x = dm3c$pos[j]
      curr.y0end = xend0
      curr.y1beg = xend1
  
      final.temp.end0 = data.frame(x=curr.x,y=curr.y0end,type='end0')
      final.temp.end1 = data.frame(x=curr.x,y=curr.y1beg,type='end1')
      final.temp = rbind(final.temp.end0,final.temp.end1)
      # final0 = rbind(final0,final.temp)
      # 
      # p2 = p2 + 
      #   geom_point(data=final.temp.end0,aes(x=x,y=y),color=rgb(1,0,0,1),pch=15) +
      #   geom_line(data=final.temp, aes(x=x,y=y),color=rgb(1,0,0,1)) +
      #   geom_point(data=final.temp.end1,aes(x=x,y=y),color=rgb(1,0,0,1),pch=18)  
      # myx = c(myx,curr.x,curr.x)
      # myy = c(myy,curr.y0end,curr.y1beg)
      # # mytype = c(mytype,'end0','end1')
      # print(paste('i=',i,',j=',j,': AFTER : ',curr.x,',',curr.y0end,' - ',curr.x,',',curr.y1beg,sep=''))
      finalend = rbind(finalend,data.frame(x0end=curr.x,y0end=curr.y0end,x1beg=curr.x,y1beg=curr.y1beg))
    }
  }
  #p2.save = p2
  #initgraph2(p2,goodbegs,goodends)
  
  for (i in 1:size(finalbeg)) {
    if (i == 1) {
      final0 = data.frame()
    }
    for (j in 1:size(finalend)) {
      if (finalbeg$x0end[i] >= finalend$x0end[j] & finalbeg$x0end[i] <= finalend$x1beg[j]) {
        if (finalbeg$y0end[i] > finalend$y1beg[j]) {next}
        print(paste(i,j,'beg x0end is between end X coords',finalend$x0end[j],'<=',finalbeg$x0end[i],'<=',finalbeg$x1beg[i],'x1beg=',finalbeg$x1beg[i]))
        final0 = rbind(final0,data.frame(x0end=finalbeg$x0end[i],y0end=finalbeg$y0end[i],x1beg=finalbeg$x1beg[i],y1beg=finalend$y1beg[j]))
      } else if (finalbeg$x1beg[i] >= finalend$x0end[j] & finalbeg$x1beg[i] <= finalend$x1beg[j]) {
        if (finalbeg$y0end[i] < finalend$y1beg[j]) {next}
        final0 = rbind(final0,data.frame(x0end=finalbeg$x0end[i],y0end=finalend$y0end[j],x1beg=finalbeg$x1beg[i],y1beg=finalbeg$y1beg[i]))
      }
    }
  }
  
  # initgraph2(p + 
  #              geom_point(data=final0,aes(x=x0end,y=y0end),fill=NA,color='red')+
  #              geom_point(data=final0,aes(x=x0end,y=y1beg),fill=NA,color='red4')+
  #              geom_point(data=final0,aes(x=x0end,y=y1beg),fill=NA,color='blue')+
  #              geom_point(data=final0,aes(x=x1beg,y=y1beg),fill=NA,color='blue4')+
  #              geom_rect(data=final0,aes(x=x0end,y=y0end,xmin=x0end,xmax=x1beg,ymin=y0end,ymax=y1beg),fill=NA,color='black'),
  #            goodbegs,goodends)
  
  #length(finalend)
  #final0 = data.frame(x=myx, y=myy, type=mytype)
  # 
  # final0 = final0[order(final0$x, final0$y),]
  # 
  # myx0end = c()
  # myy0end = c()
  # myx1beg = c()
  # myy1beg = c()
  # mypaste = c()
  # for (i in 1:length(final0$x)) {
  #   x0end = final0$x[i]
  #   y0end = final0$y[i]
  #   
  #   for (j in 1:length(final0$y)) {
  #     x1beg = final0$x[j]
  #     y1beg = final0$y[j]
  #     if (x1beg == x0end & y0end == y1beg) {next}
  #     if (length(grep(paste(x0end,y0end,x1beg,y1beg),mypaste)) == 0 & x0end >= x1beg - border0 & x0end <= x1beg + border1) {
  #       myx0end = c(myx0end, x0end)
  #       myy0end = c(myy0end, y0end)
  #       myx1beg = c(myx1beg, x1beg)
  #       myy1beg = c(myy1beg, y1beg)
  #       mypaste = c(mypaste, paste(x0end,y0end,x1beg,y1beg))
  #     }
  #   }
  # }
  # finalbeg = data.frame(x0end=myx0end, y0end=myy0end, x1beg=myx1beg,y1beg=myy1beg)
  # finalbeg = finalbeg[order(finalbeg$x0end, finalbeg$y0end,finalbeg$x1beg,finalbeg$y1beg),]
  # 
  # p + geom_segment(data=finalbeg,aes(x=x0end,y=y0end,xend=x1beg,yend=y1beg),col=rgb(1,0,0,1))
  # 
  # myx0end = c()
  # myy0end = c()
  # myx1beg = c()
  # myy1beg = c()
  # mypaste = c()
  # for (i in 1:length(final0$x)) {
  #   x0end = final0$x[i]
  #   y0end = final0$y[i]
  #   
  #   for (j in 1:length(final0$y)) {
  #     x1beg = final0$x[j]
  #     y1beg = final0$y[j]
  #     if (x1beg == x0end & y0end == y1beg) {next}
  #     if (length(grep(paste(x0end,y0end,x1beg,y1beg),mypaste)) == 0 & y0end >= y1beg -  border0 & y0end <= y1beg + border1) {
  #       myx0end = c(myx0end, x0end)
  #       myy0end = c(myy0end, y0end)
  #       myx1beg = c(myx1beg, x1beg)
  #       myy1beg = c(myy1beg, y1beg)
  #       mypaste = c(mypaste, paste(x0end,y0end,x1beg,y1beg))
  #     }
  #   }
  # }
  # finalend = data.frame(x0end=myx0end, y0end=myy0end, x1beg=myx1beg,y1beg=myy1beg)
  # finalend = finalend[order(finalend$x0end, finalend$y0end,finalend$x1beg,finalend$y1beg),]
  # 
  # segments(finalend$x0end,finalend$y0end,finalend$x1beg,finalend$y1beg,col=rgb(0,0,1,1),pch=18)
  # 
  # 
  # #-----------------!
  # 
  # finalbegpaste = c()
  # finalbegbackup = finalbeg
  # 
  # finalbeg = finalbegbackup
  # finalbeg2 = data.frame()
  # for (i in 1:dim(finalbeg)[1]) {
  #   x0end = finalbeg$x0end[i]
  #   y0end = finalbeg$y0end[i]
  #   x1beg = finalbeg$x1beg[i]
  #   y1beg = finalbeg$y1beg[i]
  #   if (y0end > y1beg) {
  #     y0end = finalbeg$y1beg[i]
  #     y1beg = finalbeg$y0end[i]
  #     finalbeg$y0end[i] = y0end
  #     finalbeg$y1beg[i] = y1beg
  #   }
  #   mypaste = paste(x0end,y0end,x1beg,y1beg)
  #   if (length(finalbegpaste[finalbegpaste == mypaste]) == 0) {
  #     finalbegpaste = c(finalbegpaste,mypaste)
  #     finalbeg2 = rbind(finalbeg2,data.frame(x0end=x0end,y0end=y0end,x1beg=x1beg,y1beg=y1beg))
  #   }
  # }
  # 
  
  final0$x0beg = final0$x0end
  final0$y0beg = final0$y1beg
  final0$x1end = final0$x0end
  final0$y1end = final0$y1beg
  final0 = final0[order(final0$x0end,final0$y0end,final0$x1beg,final0$y1beg),]
  final0$cluster = seq(1,dim(final0)[1])
  final0$divby = divby
  

  save(final0,paste('resources/CLUSTRDS/',mytitle,'.RDS',sep=''))
  return(final0)
  
  # debug drawing
  p2 = ggplot(dm2b,aes(beg,end)) +
    geom_point(aes(alpha=sqrt(sum),size=sqrt(sum))) +
    geom_point(pch=".",color="grey") +
    scale_size_continuous(range=c(0.1,1)) +
    scale_alpha_continuous(range=c(0.1,1)) +
    theme_bw() + coord_cartesian(xlim=c(0,5000/divby),ylim=c(0,5000/divby)) +    ggtitle(mytitle) +
    annotate(geom="segment",x=0,y=0,xend=5000/divby,yend=5000/divby) + theme(legend.position = "top")
  p3 = p2 +
    geom_point(data=final0,aes(x=x0end,y=y0end),color = "red2",shape=15) +
    geom_point(data=final0,aes(x=x1beg,y=y1beg),color = "red4",shape=18) +
    geom_rect(data=final0,aes(x=x0end,y=y0end,xmin=x0end,ymin=y0end,xmax=x1beg,ymax=y1beg),fill=NA,color='black')
  p4 = ggplot(dm2c,aes(x,total)) +
    geom_line() +
    theme_bw() + coord_cartesian(xlim=c(0,5000/divby),ylim=c(0,max(dm2c$total))) +    ggtitle(mytitle) +
    geom_segment(data=final0,aes(x=x0end,y=0,xend=x0end,yend=max(dm2c$total)),lty=2,color='red4')
  
  dm2$begD = ai(dm2$beg/divby)
  dm2$endD = ai(dm2$end/divby)
  dm2 = dm2[order(dm2$begD, dm2$endD),]
  dm2$cluster = 0
  for (i in 1:dim(dm2)[1]) {
    begD = dm2[i,]$begD
    endD = dm2[i,]$endD
    for (j in 1:dim(final0)[1]) {
      x0end = final0[j,]$x0end
      y0end = final0[j,]$y0end
      x1beg = final0[j,]$x1beg
      y1beg = final0[j,]$y1beg
      if (dm2[i,]$cluster != 0) {next}
      if (begD >= x0end - 1 & begD <= y0end & endD >= y0end & endD <= y1beg) {
        dm2[i,]$cluster = j
        print(paste(i,begD,endD,j))
        break
      }
    }
  }
  for (i in 1:dim(dm2)[1]) {
    begD = dm2[i,]$begD
    endD = dm2[i,]$endD
    for (j in 1:dim(final0)[1]) {
      x0end = final0[j,]$x0end
      y0end = final0[j,]$y0end
      x1beg = final0[j,]$x1beg
      y1beg = final0[j,]$y1beg
    if (dm2[i,]$cluster != 0) {next}
      if (begD >= x0end - 1 & begD <= y1beg & endD >= y0end & endD <= y1beg) {
        dm2[i,]$cluster = j
        print(paste(i,begD,endD,j))
        break
      }
    }
  }
  dm2 = dm2[order(dm2$cluster,dm2$begD,dm2$endD),]
  dm2$y = seq(1,dim(dm2)[1])
  
  mycolors = c(brewer.pal(9,"Set1"),brewer.pal(9,"Set3"))
  mybreaks = seq(0,length(mycolors)-1)
  
  p5 = ggplot(dm2[dm2$cluster!= -1,],aes(begD,endD)) +
    geom_segment(aes(x=begD,xend=endD,y=y,yend=y,color=as.factor(cluster))) +
    theme_bw() + coord_cartesian(xlim=c(0,5000/divby),ylim=c(0,dim(dm2)[1])) +    ggtitle(mytitle) +
    theme(legend.position = "top")
  
  dm2temp = aggregate(dm2$begD,by=list(dm2$cluster),min);colnames(dm2temp) = c('cluster','begDmin')
  dm2temp$begDmin = aggregate(dm2$begD,by=list(dm2$cluster),min)$x
  dm2temp$begDmax = aggregate(dm2$begD,by=list(dm2$cluster),max)$x
  dm2temp$endDmin = aggregate(dm2$endD,by=list(dm2$cluster),min)$x
  dm2temp$endDmax = aggregate(dm2$endD,by=list(dm2$cluster),max)$x
  dm2temp$ymin = aggregate(dm2$y,by=list(dm2$cluster),min)$x
  dm2temp$ymax = aggregate(dm2$y,by=list(dm2$cluster),max)$x
  dm2temp = merge(dm2temp, final0,by='cluster')
  
  p6 = p5 +
    geom_segment(data=dm2temp[dm2temp$cluster != -1,],aes(x=x0end,xend=x0end,y=ymin,yend=ymax,color=af(cluster))) +
    geom_segment(data=dm2temp[dm2temp$cluster != -1,],aes(x=y1beg,xend=y1beg,y=ymin,yend=ymax,color=af(cluster))) +
      scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks)
    #final0,aes(x=x0end,y=y),color='purple') +
    #geom_point(data=final0,aes(x=y0end,y=y),color='red4') +
    #geom_point(data=final0,aes(x=y1beg,y=y),color='blue4') +
    #geom_segment(data=final0,aes(x=x0end,y=y,xend=y0end,yend=y),color='red4') +
    # geom_segment(data=final0,aes(x=x0end,y=y,xend=y0end,yend=y),color='red4') +
    # geom_segment(data=final0,aes(x=y0end,y=y,xend=y1beg,yend=y),color='blue4') +
    # geom_segment(data=final0,aes(x=x0end,y=0,xend=x0end,yend=dim(dm2)[1]/3),lty=2,lwd=0.5,alpha=0.25,color='purple')+
    # geom_segment(data=final0,aes(x=y0end,y=dim(dm2)[1]/3,xend=y0end,yend=dim(dm2)[1]/3*2),lty=2,lwd=0.5,color='red4')+
    # geom_segment(data=final0,aes(x=y1beg,y=dim(dm2)[1]/3*2,xend=y1beg,yend=dim(dm2)[1]),lty=2,lwd=0.5,color='blue4')
  
  dm2 = dm2[order(dm2$cluster,dm2$beg,dm2$end),]
  dm2$y = seq(1,size(dm2))
  p7 = ggplot(dm2[dm2$cluster!= -1,],aes(beg,end)) +
    geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=as.factor(cluster))) +
    theme_bw() + coord_cartesian(xlim=c(0,5000),ylim=c(0,dim(dm2)[1])) +    ggtitle(mytitle) +
    theme(legend.position = "top") +
    geom_segment(data=dm2temp[dm2temp$cluster != -1,],aes(x=x0end*divby,xend=x0end*divby,y=ymin,yend=ymax,color=af(cluster))) +
    geom_segment(data=dm2temp[dm2temp$cluster != -1,],aes(x=y1beg*divby,xend=y1beg*divby,y=ymin,yend=ymax,color=af(cluster))) +
      scale_color_manual(values=mycolors,breaks=mybreaks,label=mybreaks)
    
  
  pdf("testpdf.pdf",height=30/2,width=20/2)
  grid.arrange(p2,p3,p4,p4,p5,p6,p7,p7,ncol=2,nrow=4,heights=c(2,1,2,2))
  dev.off()
}
# Figure 1: QE VR diagram update, init & term series
# Figure 2: VR frequency
# Figure 2b: Distance from TSS initiation, RLoop length
# Figure 3:  something footprint and clusters, triangles
# Figure 4: SSB
# Figure 5: sgRNA



####-end


#!---------
#end