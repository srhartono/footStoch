library(stringr)
library(rlang)
library(gridGraphics)
e = exists
cn = colnames
rn = rownames
aa = as.array
af = as.factor
ac = as.character
an = as.numeric
ai = as.integer
sm = summary

pasta = function(...) {paste(...,sep='')}
formals(pasta)$sep = ''
pasta('a','b')

size = function(df) {
  if (defined(df) & is.data.frame(df)) {
    return(dim(df)[1])
  } else if (defined(df) & is.vector(df)) {
    return(length(df))
  } else if (defined(df)) {
    return(length(df))
  } else {
    return(0)
  }
}

my.order = function(x,orderby=NA) {
  res = tryCatch(
    {
      if (defined(orderby) == FALSE) {
        x[order(x)]
      } else {
        orderbys = data.frame(orderby=orderby,neg=F,type='NORMAL',rownames=F)
        orders = list()
        if (defined(orderbys$orderby[grep('^\\-',orderbys$orderby)])) {
          orderbys$neg[grep('^\\-',orderbys$orderby)] = 'NEG'
          orderbys$orderby[grep('^\\-',orderbys$orderby)] = gsub("^(\\-)(.+)$","\\2",orderbys$orderby[grep('^\\-',orderbys$orderby)])
        }
        if (defined(orderbys$orderby[grep('\\-?rownames\\(',orderbys$orderby)])) {
          orderbys$type[grep('\\-?rownames\\(',orderbys$orderby)] = 'rownames'
        }
        print(orderbys)
        orderInd = seq(1,size(x))
        for (i in seq(size(orderbys),1,-1)) {
          if (orderbys$type[i] == 'rownames') {
            to_order = rownames(x)
          } else {
            to_order = x[,orderbys$orderby[i]]
          }
          if (orderbys$neg[i] == FALSE) {
            x = x[order(to_order,orderInd),]
            orderInd = seq(1,size(x))
          } else {
            if (is.numeric(to_order) == T) {
              to_order = to_order * -1
            } else {
              to_order = as.numeric(to_order) * -1
            }
            print(to_order)
            #to_order = to_order[-to_order]
            x = x[order(to_order,orderInd),]
            orderInd = seq(1,size(x))
          } 
        }
      }
      x
    },
    error=function(err)
    {
      message(paste("!!! srhlib.R::my.order($x,$orderby): Canot reorder this my.data type:>>",class(x)))#,'<<; $orderby is >>',head(orderby),'<<',sep=''))
      message(err)
      return(invisible(x))
    }
  )
  
  to_return = res

  return(to_return)
}

defined = function(my.data,direct=F,debug=F,verbose=F) {
  res = tryCatch (
    {
      if (!any(is.na(my.data)) & !any(is.null(my.data)) & length(my.data)> 0) { 
        # if there is any real value then return TRUE
        if (is.data.frame(my.data)) {
          if (dim(my.data)[1] == 0) {
            FALSE
          } else {
            TRUE
          }
        } else {
          TRUE
        }
      } else if (length(my.data) == 0) {
        FALSE #if (direct == F) {FALSE} else {NULL}
      } else if (length(my.data[is.na(my.data)]) == length(my.data)) {
        FALSE #if (direct == F) {FALSE} else {NULL}
      } else if (length(my.data[is.null(my.data)]) == length(my.data)) {
        FALSE ##if (direct == F) {FALSE} else {NULL}
      } else {#some na, some null, but not all
        TRUE ##if (direct == F) {TRUE} else {my.data}
      }
    },
    error = function(err) {
      if (!e(quo_name(enquo(my.data)), where = .GlobalEnv)) {
        return(invisible(FALSE)) #if (direct == F) {return(invisible(FALSE))} else {return(invisible(NULL))}
      } else {
        return(invisible(TRUE)) ##if (direct == F) {return(TRUE)} else {return(invisible(my.data))}
      }
    },
    warning = function(cond) {
      return(TRUE) #if (direct == F) {return(TRUE)} else {return(invisible(my.data))}
    }
  )
  
  to_return = res
  if (direct == T) {
    if (is_false(res)) {
      to_return = NULL
    } else {
      to_return = my.data
    }
  }

  return(to_return)
}

smdf = function(df) {
  for (i in 1:length(cn(df))) {
    my.cn = cn(df)[i]
    cat(' - ',i,'. ',my.cn,': ',head(unique(df[,cn(df) == my.cn]),n=5),'\n',sep=' ')
  }
  cat(' - ',length(cn(df))+1,'. ','dim(df): ',dim(df),'\n',sep=' ')
  
  to_return=NULL
  return(to_return)
}

reorder_y = function(df) {
  df$y = seq(1,dim(df)[1])

  to_return = df
  return(to_return)
}

#my.dir = "/Users/mitochy/Work/Project/Ethan//triclust/"

triclust_get_divby = function(dm4=data.frame(),my.title=NA,gp=list()) {
  if (defined(gp$my.title)) {
    my.title = gp$my.title
  }
  
  gp$triclust.divby = 50
  gp$triclust.thres1 = max(as.integer(sqrt(dim(dm4)[1]/9230)*25),5)
  gp$triclust.thres2 = max(as.integer(sqrt(dim(dm4)[1]/9230)*25),5)
  gp$triclust.sample = my.title
  if (defined(grep("T7",my.title,ignore.case=TRUE,perl=TRUE))) {
    #T7: 25, 100, 100
    gp$triclust.divby = 25
    gp$triclust.thres1 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.thres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.sample = 'T7_init'
  } else if (defined(grep("PFC53.+T3Term",my.title,ignore.case = TRUE,perl=TRUE)) ) {
    #pFC53_T3Term: 15,max(50 etc)
    gp$triclust.divby = 15
    gp$triclust.thres1 = max(50,as.integer(sqrt(dim(dm4)[1]/9230)*100))
    gp$triclust.thres2 = max(50,as.integer(sqrt(dim(dm4)[1]/9230)*100))
    gp$triclust.sample = 'pFC53_T3Term'
  } else if (defined(grep("pFC53.+ApaLI",my.title,ignore.case = TRUE,perl=TRUE)) ) {
    #pFC53_ApaLI: 25,as.integer
    gp$triclust.divby = 25
    gp$triclust.thres1 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.thres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.sample = 'pFC53_SSB_ApaLI'
  } else if (defined(grep("pFC9.+ApaLI",my.title,ignore.case = TRUE,perl=TRUE)) ) {
    #pFC9_ApaLI: 25,as.integer
    gp$triclust.divby = 25
    gp$triclust.thres1 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.thres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.sample = 'pFC9_SSB_ApaLI'
  } else if (defined(grep("pFC9",my.title,ignore.case = TRUE,perl=TRUE)) ) {
    #pFC9_sgRNA: 15,as.integer
    gp$triclust.divby = 50
    gp$triclust.thres1 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.thres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.sample = 'pFC9_sgRNA'
  }
  
#  print(paste(my.title,gp$triclust.divby,gp$triclust.thres1,gp$triculust.thres2))
  to_return = gp
  return(to_return)
}

findedge = function(x) {
  
  # Find the edges
  edges <- which(abs(diff(x)) > 1)
  if (!0 %in% edges) {
    edges = c(0,edges)
  }
  if (!size(x) %in% edges) {
    edges = c(edges,size(x))
  }
  return(edges)
}

get_meanpos = function(pos,windowsmooth,stepsmooth) {
  my.len = length(pos)
  meanpos = c()
  for (i in seq(1,my.len, stepsmooth)) {
    i0 = i
    i1 = min(i+gp$windowsmooth,my.len)
    meanpos[i] = mean(pos[seq(i0,i1)])
  }
  return(meanpos)
}

get_closest = function(d1,d2) {
  test1 = as.numeric(rep(d1,length(d2)))
  test2 = as.numeric(d2)
  test3 = data.frame(a = test1, b = test2)
  test3$diff = abs(test3$a - test3$b)  
  test3 = test3[order(test3$diff,test3$b,test3$a),]
  return(test3$b[1])
}


fa.get_CytosinePosition = function(my.fa) {
  nuc = as.array(as.character(strsplit(
    x = as.character(my.fa$amp.seq[[1]]),
    split = '',
    perl = T
  )
  [[1]]))
  df = data.frame(nuc=nuc,isC = rep(0,length(nuc)))
  df[df$nuc == 'C',]$isC = 1
  df$index = my.fa$amp.beg + seq(1,dim(df)[1])-1
  head(df)
  #   nuc isC
  # 1   G   0
  # 2   G   0
  # 3   G   0
  # 4   G   0
  # 5   T   0
  # 6   C   0
  
  
  if (defined(df[df$isC == 1,])) {
    df = df[df$isC == 1,]
    df$Cindex = my.fa$amp.beg + seq(1,dim(df)[1])-1
    df$diff = c(diff(df$index),0)
    to_return = df
  } else {
    to_return = data.frame()
  }
  return(to_return)
}
fa.transform_PositionIntoCytosinePosition = function(my_df,my_fa,my_bys,verbose=F) {
  my_fa_c       = fa.get_CytosinePosition(my_fa)
  my.verbose=verbose
  my_df_orig = my_df
  #print(head(my_bys))
  for (my_by in my_bys) {
    my_by_C = pasta(my_by,'.C')
    my_by_C_diff = pasta(my_by,'.C_diff')
    my_by_orig = pasta(my_by,'.orig')
    
    my_fa_c[,my_by] = my_fa_c$index
    my_fa_c[,my_by_C] = my_fa_c$Cindex
    my_fa_c[,my_by_C_diff] = my_fa_c$diff
    my_fa_c.subset = subset(my_fa_c,select=c(my_by,my_by_C,my_by_C_diff))

    my_df[,my_by] = apply(data.frame(my_df[,my_by]),1,get_closest,d2=my_fa_c.subset[,my_by])
    my_df = merge(my_df,my_fa_c.subset,by=my_by)
    my_df[,my_by] = my_df[,my_by_C]
  
    sanity_check(size(my_df_orig)==size(my_df),paste('-> size my_df_orig',size(my_df_orig),'==size my_df_curr',size(my_df),'by',my_by),verbose=my.verbose)
  }  
  to_return = my_df
  return(to_return)
}



fa.calc_dens = function(N,n1,n2,n3) {
  if (N == 0) {return(0)}
  if (n1/N*n2/N == 0) {return(0)}
  return(ai(0.5+((n3*N)/(n1*n2)*1000))/1000)
}
fa.calc_cont = function(N,n1) {
  if (n1 == 0) {return(0)}
  return(ai(0.5+(n1)/N*1000)/1000)
}
fa.calc_skew = function(N,n1,n2) {
  if (n1+n2 == 0) {return(0)}
  if (n1-n2 > 0) {my.sign = 1} else {my.sign = -1}
  return(ai((my.sign * 0.5) + (n1-n2)/(n1+n2)*1000)/1000)
}
fa.get_chunk = function(my.position,my.seq,gp,my.window=NA,my.step=NA) {
  if (is.na(my.window)) {
    if ('CpGprof.window' %in% gp) {
      my.window = gp$GCprof.window
    } else {
      my.window = 200
    }
  }
  if (is.na(my.step)) {
    if ('CpGprof.step' %in% gp) {
      my.step = gp$GCprof.step
    } else {
      my.step = 1
    }
  }
  
  if (my.window > nchar(my.seq)) {my.window = nchar(my.seq)}
  
  iMin = max(1,min(my.position,nchar(my.seq)-1))
  iMax = max(1,min(my.position+my.window,nchar(my.seq)))
  chunk = substr(my.seq,iMin,iMax)
  return(chunk)
}
fa.calculate_CpGprof = function(my.seq) {
  N    = nchar(my.seq)
  nA   = str_count(my.seq,"A")
  nC   = str_count(my.seq,"C")
  nG   = str_count(my.seq,"G")
  nT   = str_count(my.seq,"T")
  nN   = str_count(my.seq,"N")
  nCpG = str_count(my.seq,"CG")
  nG2 = str_count(my.seq,"GG")
  nG3 = str_count(my.seq,"GGG")
  nG4 = str_count(my.seq,"GGGG")
  G2dens = fa.calc_dens(N,nG,nG,nG2)
  G3dens = fa.calc_dens(N,nG2,nG,nG3)
  G4dens = fa.calc_dens(N,nG3,nG,nG4)
  CGdens = fa.calc_dens(N,nC,nG,nCpG)
  GCcont = fa.calc_cont(N,nG+nC)
  GCskew = fa.calc_skew(N,nG,nC)
  GAskew = fa.calc_skew(N,nG,nA)
  GTskew = fa.calc_skew(N,nG,nT)
  ATskew = fa.calc_skew(N,nA,nT)
  Acont  = fa.calc_cont(N,nA)
  Ccont  = fa.calc_cont(N,nC)
  Gcont  = fa.calc_cont(N,nG)
  Tcont  = fa.calc_cont(N,nT)
  G2cont  = fa.calc_cont(N,nG2)
  G3cont  = fa.calc_cont(N,nG3)
  G4cont  = fa.calc_cont(N,nG4)
  return(c(CGdens,GCcont,GCskew,GAskew,GTskew,ATskew,Acont,Ccont,Gcont,Tcont,G2cont,G3cont,G4cont,G2dens,G3dens,G4dens))
}
fa.get_CpGprof = function(my.fa,gp,my.window=NA,my.step=NA,unit.test=FALSE,with.chunk=FALSE) {
  if (is.na(my.window)) {
    if ('CpGprof.window' %in% gp) {
      my.window = gp$GCprof.window
    } else {
      my.window = 200
    }
  }
  if (is.na(my.step)) {
    if ('CpGprof.step' %in% gp) {
      my.step = gp$GCprof.step
    } else {
      my.step = 1
    }
  }
  
  nuc = aa(ac(strsplit(
    x = ac(my.fa$seq[[1]]),
    split = '',
    perl = T
  )[[1]]))
  df = data.frame(nuc=nuc,index = my.fa$amp.beg + seq(1,nchar(my.fa$seq)))
  df$chunk = apply(data.frame(test1$index-my.window/2),1,fa.get_chunk,my.seq=my.fa$seq,gp=gp)
  dfchunk = as.data.frame(t(as.matrix(apply(data.frame(df$chunk),1,fa.calculate_CpGprof),nrow=size(df$chunk),byrow = T)))
  
  colnames(dfchunk) = c('CGdens','GCcont','GCskew','GAskew','GTskew','ATskew','Acont','Ccont','Gcont','Tcont','G2cont','G3cont','G4cont','G2dens','G3dens','G4dens')
  df = cbind(df,dfchunk)
  
  
  #unit test
  if (unit.test == TRUE) {
    ps.CpGprof(df,gp=gp,my.bed=my.bed.orig,my.title=gp$my.title.wrap)
  }
  #df[df$nuc == 'C',]$isC = 
  #  df$chunk = apply(df$index,1,fa.get_chunk(),my.seq=my.fa$seq,gp=gp,my.window=my.window,my.step=my.step)
  if (with.chunk == FALSE) {df = (subset(df,select=c(-chunk)))}
  to_return = list(df=df,names=colnames(dfchunk))
  return(to_return)
}


df.get_density = function(df,by,gp,CpGprof,group=NA) {
  temp = df[df$beg >= 500 & df$beg <= 1300,]
  if (is.na(group)) {
    temp$tempgroup = 1
  } else {
    temp$tempgroup = temp[,group]
  }
  temp$count = 1
  tempcount = aggregate(temp$count,by=list(temp[,by],temp$tempgroup),sum); colnames(tempcount)= c(by,'tempgroup',pasta(by,'.sum'))
  temptotal = aggregate(temp$count,by=list(temp$tempgroup),sum); colnames(temptotal) = c('tempgroup',pasta(by,'.total'))
  temp2 = tempcount
  #temp2 = merge(temp2,tempcount,by=c(by,'tempgroup'))
  temp2 = merge(temp2,temptotal,by=c('tempgroup'))
  temp2[,pasta(by,'.perc')] = ai(temp2[,pasta(by,'.sum')] / temp2[,pasta(by,'.total')] * 10000 + 0.5)/100
  temp3 = CpGprof
  temp3[,by] = temp3$index
  temp3 = merge(temp3,temp2,by=by)
  return(temp3)
}

mysample = function(x) {return(sample(x,size = length(x),replace = FALSE))}

get_cor = function(df,dfd2,by,my.params,gp=list()) {
  bytemp = by
  bytemp.perc = pasta(bytemp,'.perc')
  # dftemp = cbind(subset(df,select=c(bytemp,'end')),dfchunk)
  #dftemp = dfd2[dfd2$beg >= 600 & dfd2$beg <= 821,!colnames(dfd2) %in% colnames(df)]
  dftemp = dfd2[dfd2[,by] >= 641,!colnames(dfd2) %in% colnames(df)]
  colnames(dftemp)[colnames(dftemp) == 'pos'] = bytemp
  colnames(dftemp)
  dftemp.orig = dftemp[,colnames(dftemp) %in% c('beg','end',bytemp.perc,my.CpGprof.name)]
  dftemp.shuf = as.data.frame(apply(dftemp.orig,2,mysample))
  corrplot.col = rev(brewer.pal(n = 11,'RdBu'))
  
  M.orig.s = cor(dftemp.orig,method='spearman')
  M.orig.p = cor(dftemp.orig,method='pearson')
  M.shuf.s = cor(dftemp.shuf,method='spearman')
  M.shuf.p = cor(dftemp.shuf,method='pearson')
  # corrplot(M.orig.s,method='color',type='upper',col=corrplot.col,order = 'hclust'); corrplot(M.orig.s, add = TRUE, type = 'upper', method = 'number', order = 'hclust',number.cex = 0.5, col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n',bg = NA)
  # grid.echo()
  # p1.corr.orig.s = grid.grab()
  # corrplot(M.orig.p,method='color',type='upper',col=corrplot.col,order = 'hclust'); corrplot(M.orig.p, add = TRUE, type = 'upper', method = 'number', order = 'hclust',number.cex = 0.5, col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n',bg = NA)
  # grid.echo()
  # p1.corr.orig.p = grid.grab()
  # corrplot(M.shuf.s,method='color',type='upper',col=corrplot.col,order = 'hclust'); corrplot(M.shuf.s, add = TRUE, type = 'upper', method = 'number', order = 'hclust',number.cex = 0.5, col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n',bg = NA)
  # grid.echo()
  # p1.corr.shuf.s = grid.grab()
  # corrplot(M.shuf.p,method='color',type='upper',col=corrplot.col,order = 'hclust'); corrplot(M.shuf.p, add = TRUE, type = 'upper', method = 'number', order = 'hclust',number.cex = 0.5, col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n',bg = NA)
  # grid.echo()
  # p1.corr.shuf.p = grid.grab()
  # pdf(pasta('results/simple/corrplot_',bytemp,'_VR',my.VR,'.pdf'),width=10,height=10)
  # corrplot(M.orig.s,method='color',type='upper',col=corrplot.col,order = 'hclust'); corrplot(M.orig.s, add = TRUE, type = 'upper', method = 'number', order = 'hclust',number.cex = 0.5, col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n',bg = NA)
  # p2
  # #grid.arrange(arrangeGrob(p1.corr.shuf.s,p1.corr.orig.s,p1.corr.shuf.p,p1.corr.orig.p,nrow=2,ncol=2),p2,nrow=2,ncol=1,heights=c(2,0.5))
  # dev.off()  
  # 
  
  p2 = ggplot(dfd2,aes(dfd2[,bytemp],G2cont)) +#[dfd2$dfd2[,bytemp] >= 620 & dfd2$dfd2[,bytemp] <= 1500,],aes(dfd2[,bytemp],G2cont)) + 
    geom_line(aes(dfd2[,bytemp],smooth.spline(CGdens,spar=1)$y),color='cornflowerblue') +
    geom_line(aes(dfd2[,bytemp],smooth.spline(GCcont,spar=1)$y),color='green4') +
    geom_line(aes(dfd2[,bytemp],smooth.spline(GCskew+0.5,spar=1)$y),color='red4') +
    geom_line(aes(dfd2[,bytemp],smooth.spline(dfd2[,bytemp.perc]/max(dfd2[,bytemp.perc]),spar=1)$y),color='black') +
    geom_line(aes(dfd2[,bytemp],smooth.spline(Gcont,spar=1)$y),lty=2,color='green4') +
    geom_line(aes(dfd2[,bytemp],smooth.spline(G2dens,spar=1)$y),color='purple1') +
    geom_line(aes(dfd2[,bytemp],smooth.spline(G3dens,spar=1)$y),color='purple3') +
    geom_line(aes(dfd2[,bytemp],smooth.spline(G4dens,spar=1)$y),color='purple4') +
    annotate(geom='text',x=300,y=0.60,label='CpGdens _____',color='cornflowerblue',hjust=0) +
    annotate(geom='text',x=300,y=0.55,label='GC % _____',color='green4',hjust=0) +
    annotate(geom='text',x=300,y=0.50,label='GCSkew _____',color='red4',hjust=0) +
    annotate(geom='text',x=300,y=0.45,label=paste('Rloop',bytemp,'density'),color='black',hjust=0) +
    annotate(geom='text',x=300,y=0.40,label='G % - - - -',color='green4',hjust=0) +
    annotate(geom='text',x=300,y=0.35,label='G2dens _____',color='purple1',hjust=0) +
    annotate(geom='text',x=300,y=0.30,label='G3dens _____',color='purple3',hjust=0) +
    annotate(geom='text',x=300,y=0.25,label='G4dens _____',color='purple4',hjust=0) +
    
    ggtitle(pasta(gp$my.title.wrap,'_',by)) +
    
    geom_rect(data=my.bed[grep('Barcode|_diff',my.bed$feature,invert = TRUE),],aes(x=beg,y=0,xmin=beg,xmax=end,ymin=0,ymax=1,fill=af(feature)),color=rgb(0,0,0,0.1),alpha=0.1) +
    geom_text(data=my.bed[grep('Barcode|_diff',my.bed$feature,invert = TRUE),],aes(x=(beg+end)/2,y=1,label=feature,angle=90),color=rgb(0,0,0),size=4,hjust=0) +
    # annotate(geom='rect',xmin=,ymin=0,xmax=841,ymax=1,color=rgb(0,0,0,0.15),fill='green2',alpha=0.15) +
    # annotate(geom='rect',xmin=641,ymin=0,xmax=841,ymax=1,color=rgb(0,0,0,0.15),fill='green2',alpha=0.15) +
    coord_cartesian(xlim=c(300,1500),ylim=c(0,2)) +
    theme_bw() + theme(legend.position = 'none',panel.grid = element_blank())
  
  #dfd2$bytemp.perc.cat = cut(dfd2$bytemp.perc,breaks = c(0,quantile(dfd2$bytemp.perc,probs=seq(0,1,0.25))))#,labels=T)#,max(dfd2$bytemp.perc)+1),labels = F)
  #ggplot(dfd2,aes(bytemp.perc.cat,G2cont)) +
  #  geom_boxplot(aes(y=G2dens,fill=af(bytemp.perc.cat)))
  #ggplot(dfd2,aes(bytemp.perc.cat,G3dens)) +
  #  geom_boxplot(aes(y=G3dens,fill=af(bytemp.perc.cat)))
  
  pdf(pasta('results/simple/corrplot_',gp$my.title,'_',by,'.pdf'),width=10,height=10)
  print(p2)

  M.orig.s = cor(dftemp.orig,method='spearman'); 
  corrplot(M.orig.s,title=pasta('\n',gp$my.title,'_',by,'_spearman'),method='color',type='full',col=corrplot.col,order = 'hclust'); corrplot(M.orig.s, add = TRUE, type = 'full', method = 'number', order = 'hclust',number.cex = 0.5, col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n',bg = NA)
  M.orig.s = as.data.frame(M.orig.s)
  M.orig.s$gene = my.params$gene; M.orig.s$treat = my.params$treat;M.orig.s$thres = my.params$thres; M.orig.s$VR = my.VR;
  M.orig.s$type = 'spearman'
  
  M.orig.s.small = cor(dftemp.orig[dftemp.orig[,by] >= 641 & dftemp.orig[,by] <= 841,],method='spearman'); 
  corrplot(M.orig.s.small,title=pasta('\n',gp$my.title,'_',by,'_spearman.VRonly'),method='color',type='full',col=corrplot.col,order = 'hclust'); corrplot(M.orig.s.small, add = TRUE, type = 'full', method = 'number', order = 'hclust',number.cex = 0.5, col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n',bg = NA)
  M.orig.s.small = as.data.frame(M.orig.s.small)
  M.orig.s.small$gene = my.params$gene; M.orig.s.small$treat = my.params$treat;M.orig.s.small$thres = my.params$thres; M.orig.s.small$VR = my.VR;
  M.orig.s.small$type = 'spearman.VRonly'
  
  M.orig.p = cor(dftemp.orig,method='pearson')
  corrplot(M.orig.p,title=pasta('\n',gp$my.title,'_',by,'_pearson'),method='color',type='full',col=corrplot.col,order = 'hclust'); corrplot(M.orig.p, add = TRUE, type = 'full', method = 'number', order = 'hclust',number.cex = 0.5, col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n',bg = NA)
  M.orig.p = as.data.frame(M.orig.p)
  M.orig.p$gene = my.params$gene; M.orig.p$treat = my.params$treat;M.orig.p$thres = my.params$thres; M.orig.p$VR = my.VR;
  M.orig.p$type = 'pearson'
  
  M.orig.p.small = cor(dftemp.orig[dftemp.orig[,by] >= 641 & dftemp.orig[,by] <= 841,],method='pearson')
  corrplot(M.orig.p.small,title=pasta('\n',gp$my.title,'_',by,'_pearson.VRonly'),method='color',type='full',col=corrplot.col,order = 'hclust'); corrplot(M.orig.p.small, add = TRUE, type = 'full', method = 'number', order = 'hclust',number.cex = 0.5, col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n',bg = NA)
  M.orig.p.small = as.data.frame(M.orig.p.small)
  M.orig.p.small$gene = my.params$gene; M.orig.p.small$treat = my.params$treat;M.orig.p.small$thres = my.params$thres; M.orig.p.small$VR = my.VR;
  M.orig.p.small$type = 'pearson.VRonly'
  
  dev.off() 
  
  dftemp2 = M.orig.s
  dftemp2 = rbind(dftemp2,M.orig.s.small)
  dftemp2 = rbind(dftemp2,M.orig.p)
  dftemp2 = rbind(dftemp2,M.orig.p.small)
  to_return = dftemp2
  
  return(to_return)
}

# fa.calculate_GCprof = function(my_fa,gp,my.window=NA,my.step=NA) {
#   if (is.na(my.window)) {
#     if ('CpGprof.window' %in% gp) {
#       my.window = gp$GCprof.window
#     } else {
#       my.window = 200
#     }
#   }
#   if (is.na(my.step)) {
#     if ('CpGprof.step' %in% gp) {
#       my.step = gp$GCprof.step
#     } else {
#       my.step = 1
#     }
#   }
#   my.step   = gp$GCprof.step
#   my.seq = my.fa$seq[1]
#   iMin = 1
#   iMax = max(my.window,nchar(my.seq)-my.window)
#   if (my.window > nchar(myseq)) {
#     my.window = nchar(myseq)
#   }
#   df = data.frame()
#   for (i in seq(iMin,iMax,my.step)) {
#     chunk = substr(my.seq,i,my.window)
#     nA = str_count(chunk,'A')
#     nC = str_count(chunk,'C')
#     nG = str_count(chunk,'G')
#     nT = str_count(chunk,'T')
#     nN = str_count(chunk,'N')
#     N  = nA + nC + nG + nT + nN
#     # GC Skew
#     GCskew = 0
#     if (nG + nC > 0) {GCskew = ai((nG - nC)/(nG + nC) * 100) / 100}
#     GCcont = 0
#     if (N > 0) {GCcont = ai(GC = (nG + nC) / (N) * 100 + 0.5)/100}
#   }
# }
sanity_check = function(arg,my.message='',verbose=F) {
  if (arg == FALSE) {
    message(my.message,'=',arg)
    return(invisible(FALSE))
  } else {
    if (verbose == T) {
      cat(my.message,'=',arg,'\n')
    }
    return(invisible(TRUE))
  }
}

#test1=fa.get_NucleotidePosition(my.fa)


#dfd.beg3$by='beg'
#dfd.end3$by='end'
#dfd.beg3.orig = dfd.beg3
#dfd.end3.orig = dfd.end3

dfd.beg3.p = mycorplot(dfd.beg3.orig,'beg')
dfd.end3.p = mycorplot(dfd.end3.orig,'end')
grid.arrange(dfd.beg3.p,dfd.end3.p,nrow=1,ncol=2)
mycorplot=  function(dftemp,by) {
  bytemp = by
  bytemp.perc = pasta(bytemp,'.perc')
  #dftemp = melt(dftemp,id.vars = c('gene','treat','thres','VR','type','by',bytemp.perc))
  dftemp = melt(dftemp[grep(bytemp.perc,rownames(dftemp)),],id.vars = c('gene','treat','thres','VR','type','by'))
  head(dftemp)
  dim(dftemp[grep('VRonly',dftemp$type),])
  dftemp$type2 = 'all'
  dftemp[grep('VRonly',dftemp$type),]$type2 = 'VRonly'
  dftemp$type3 = dftemp$type
  dftemp[grep('VRonly',dftemp$type),]$type3 = gsub('^(.+).VRonly$','\\1',dftemp[grep('VRonly',dftemp$type),]$type,perl=T)
  dftemp = dftemp[order(dftemp$variable),]
  p = ggplot(dftemp,aes(variable,value)) +
    geom_boxplot(aes(fill=type2)) + facet_grid(type3~.) +
    annotate(geom='segment',x=0,xend=length(unique(dftemp$variable)),y=0,yend=0,lty=2) +
    theme_bw() + theme(panel.grid = element_blank())
  return(p)
}
