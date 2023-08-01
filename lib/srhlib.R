library(stringr)
library(rlang)
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
  colnames(dfchunk) = c('CGdens','GCcont','GCskew')
  df = cbind(df,dfchunk)
  
  #unit test
  if (unit.test == TRUE) {
    ps.CpGprof(df,gp=gp,my.bed=my.bed.orig,my.title=gp$my.title.wrap)
  }
  #df[df$nuc == 'C',]$isC = 
#  df$chunk = apply(df$index,1,fa.get_chunk(),my.seq=my.fa$seq,gp=gp,my.window=my.window,my.step=my.step)
  if (with.chunk == FALSE) {return(subset(df,select=c(-chunk)))}
  return(df)
}

 # test1 = fa.get_NucleotidePosition(my.fa=my.fa,gp=gp)
# test1$chunk.000 = apply(data.frame(test1$index-000),1,fa.get_chunk,my.seq=my.fa$seq,gp=gp)
# #test1$chunk.050 = apply(data.frame(test1$index-050),1,fa.get_chunk,my.seq=my.fa$seq,gp=gp)
# test1$chunk.100 = apply(data.frame(test1$index-100),1,fa.get_chunk,my.seq=my.fa$seq,gp=gp)
# #test1$chunk.150 = apply(data.frame(test1$index-150),1,fa.get_chunk,my.seq=my.fa$seq,gp=gp)
# test1$chunk.200 = apply(data.frame(test1$index-200),1,fa.get_chunk,my.seq=my.fa$seq,gp=gp)
# 
# test3 = apply(data.frame(test1$chunk.000),1,fa.calculate_CpGprof)
# test2.000 = as.data.frame(t(as.matrix(apply(data.frame(test1$chunk.000),1,fa.calculate_CpGprof),nrow=size(test1$chunk000),byrow = T)))
# colnames(test2.000) = c('CGdens000','GCcont000','GCskew000')
# test1 = cbind(test1,test2.000)
# 
# test2.100 = as.data.frame(t(as.matrix(apply(data.frame(test1$chunk.100),1,fa.calculate_CpGprof),nrow=size(test1$chunk100),byrow = T)))
# colnames(test2.100) = c('CGdens100','GCcont100','GCskew100')
# test1 = cbind(test1,test2.100)
# 
# test2.200 = as.data.frame(t(as.matrix(apply(data.frame(test1$chunk.200),1,fa.calculate_CpGprof),nrow=size(test1$chunk200),byrow = T)))
# colnames(test2.200) = c('CGdens200','GCcont200','GCskew200')
# test1 = cbind(test1,test2.200)
# 

# plot(NA,xlim=c(0,2500),ylim=c(0,1.4))
# lines(test1$index,test1$CGdens000,col='blue2')
# lines(test1$index,test1$CGdens100,col='blue3')
# lines(test1$index,test1$CGdens200,col='blue4')
# lines(test1$index,test1$GCcont000,col='green2')
# lines(test1$index,test1$GCcont100,col='green3')
# lines(test1$index,test1$GCcont200,col='green4')
# lines(test1$index,test1$GCskew000+0.5,col='red2')
# lines(test1$index,test1$GCskew100+0.5,col='red3')
# lines(test1$index,test1$GCskew200+0.5,col='red4')
# segments(0,0.5,2500,0.5,lty=2)

fa.calculate_CpGprof = function(my.seq) {
    N    = nchar(my.seq)
    nA   = str_count(my.seq,"A")
    nC   = str_count(my.seq,"C")
    nG   = str_count(my.seq,"G")
    nT   = str_count(my.seq,"T")
    nN   = str_count(my.seq,"N")
    nCpG = str_count(my.seq,"CG")
    CGdens = fa.calc_dens(N,nC,nG,nCpG)
    GCcont = fa.calc_cont(N,nG,nC)
    GCskew = fa.calc_skew(N,nG,nC)
    return(c(CGdens,GCcont,GCskew))
}

fa.calc_dens = function(N,n1,n2,n3) {
  if (N == 0) {return(0)}
  if (n1/N*n2/N == 0) {return(0)}
  return(ai(0.5+((n3*N)/(n1*n2)*1000))/1000)
}
fa.calc_cont = function(N,n1,n2) {
  if (n1+n2 == 0) {return(0)}
  return(ai(0.5+(n1+n2)/N*1000)/1000)
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

fa.calculate_GCprof = function(my_fa,gp,my.window=NA,my.step=NA) {
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
  my.step   = gp$GCprof.step
  my.seq = my.fa$seq[1]
  iMin = 1
  iMax = max(my.window,nchar(my.seq)-my.window)
  if (my.window > nchar(myseq)) {
    my.window = nchar(myseq)
  }
  df = data.frame()
  for (i in seq(iMin,iMax,my.step)) {
    chunk = substr(my.seq,i,my.window)
    nA = str_count(chunk,'A')
    nC = str_count(chunk,'C')
    nG = str_count(chunk,'G')
    nT = str_count(chunk,'T')
    nN = str_count(chunk,'N')
    N  = nA + nC + nG + nT + nN
    # GC Skew
    GCskew = 0
    if (nG + nC > 0) {GCskew = ai((nG - nC)/(nG + nC) * 100) / 100}
    GCcont = 0
    if (N > 0) {GCcont = ai(GC = (nG + nC) / (N) * 100 + 0.5)/100}
  }
}

fa.calculate_chunk = function(my.seq,nuc1,nuc2,type) {
  if (type == 'skew') {
    n1 = str_count(my.seq,nuc1)
    n1 = str_count(my.seq,nuc1)
  }
}


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

test1=fa.get_NucleotidePosition(my.fa)
