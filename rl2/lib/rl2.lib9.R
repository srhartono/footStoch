library(rlang)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(dplyr)
library(stringr)
ac = as.character
an = as.numeric
ai = as.integer
aa = as.array
af = as.factor
am = as.matrix
ad = as.data.frame
cn = colnames
rn = rownames
e = exists
h  = head
sm = summary
s  = smooth
ss = smooth.spline

h = function(x,n=6) {return(rbind(head(x,n=n),tail(x,n=n)))}
stretch = function(x,n) {
  x = x[ai(seq(0,n-1)/(n-1) * length(x))]
}

reset = function() {
  tryCatch({while (1) {dev.off()}},error=function(cond){return(invisible(NULL))})
}

lay = function(n=NA,n1=NA,n2=NA,byrow=TRUE,verbose=F,debug=F) {
  suppressWarnings(tryCatch({
  if (is.na(n) & is.na(n1) & is.na(n2)) {
    n = 1
    n1 = 1
    n2 = 1
  } else {
    if (!is.na(n)) {
      n1 = ai(sqrt(n))
      n2 = ai(sqrt(n))
      if (n1 < sqrt(n)) {
        if (byrow == TRUE) {
          n1 = (ai(sqrt(n)) + 1)
          n2 = ai(n / n1)+1
        } else {
          n2 = (ai(sqrt(n)) + 1)
          n1 = ai(n / n2)+1
        }
      }
    } else {
      if (is.na(n1)) {
        n1 = 1
      }
      if (is.na(n2)) {
        n2 = 1
      }
#      if (n2 / n1 >= 3) {
#        n1 = n2
#      } else if (n1 / n2 >= 3) {
#          n2 = n1
#      }
      n = n1 * n2
    }
  }
  reset()
  par(mar = c(1, 1, 1, 1))
  if (verbose == TRUE) {
    print(paste('n1',n1,'n2',n2,'n',n))
  }
  layout(matrix(seq(1,n1*n2),n1,n2,byrow=byrow))
  if (debug == TRUE) {
    for (i in 1:(n1*n2)) {
      plot(seq(i-5,i+5),seq(i,i+10),xlim=c(1-10,n1*n2+10),ylim=c(1-10,n1*n2+10),type='p',pch=15,cex=0.2,col='black')
      lines(seq(i-5,i+5),seq(i,i+10),col=i)
      text(i,i,label=paste('debug',i),col=i)
    }
  }
  },
  error = function(cond) {
    n1.orig = n1
    if (is.na(n1.orig)) {n1.orig = 'NA (default)'}
    n2.orig = n2
    if (is.na(n2.orig)) {n2.orig = 'NA (default)'}
    n.orig = n
    if (is.na(n.orig)) {n.orig = 'NA (default)'}
    n1 = 2
    n2 = 2
    n = 4
    reset()
    if (verbose == TRUE) {
      print(paste('n1',n1,'n2',n2,'n',n))
    }

    layout(matrix(1,2,3,4),2,2)
    if (debug == TRUE) {
      for (i in 1:(n1*n2)) {
        plot(seq(i-5,i+5),seq(i,i+10),xlim=c(1-10,n1*n2+10),ylim=c(1-10,n1*n2+10),type='p',pch=15,cex=0.2,col='black')
        lines(seq(i-5,i+5),seq(i,i+10),col=i)
        text(i,i,label=paste('debug','error','\nn=',n.orig,'\nn1=',n1.orig,'\nn2=',n2.orig,'\nso using n1=2, n2=2, n=4'),col=i)
      }
    }
    message(cond)
  }
  ))
}

av = function(x) {
  return(x[[1]])
}

join = function(string,sep='') {
  return(paste(string,collapse=sep))
}

split = function(string,sep='') {
  x = strsplit(string,sep)
  return(x)
}

pasta = paste
formals(pasta)$sep = ''

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

myorder = function(x,orderby=NA) {
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
      message(paste("!!! srhlib.R::myorder($x,$orderby): Canot reorder this mydata type:>>",class(x)))#,'<<; $orderby is >>',head(orderby),'<<',sep=''))
      message(err)
      return(invisible(x))
    }
  )
  
  to_return = res

  return(to_return)
}

defined = function(mydata,direct=F,debug=F,verbose=F) {
  res = tryCatch(
    {
      if (!any(is.na(mydata)) & !any(is.null(mydata)) & length(mydata)> 0) { 
        # if there is any real value then return TRUE
        if (is.data.frame(mydata)) {
          if (dim(mydata)[1] == 0) {
            FALSE
          } else {
            TRUE
          }
        } else {
          TRUE
        }
      } else if (length(mydata) == 0) {
        FALSE #if (direct == F) {FALSE} else {NULL}
      } else if (length(mydata[is.na(mydata)]) == length(mydata)) {
        FALSE #if (direct == F) {FALSE} else {NULL}
      } else if (length(mydata[is.null(mydata)]) == length(mydata)) {
        FALSE ##if (direct == F) {FALSE} else {NULL}
      } else {#some na, some null, but not all
        TRUE ##if (direct == F) {TRUE} else {mydata}
      }
    },
    error = function(err) {
      if (!e(quo_name(enquo(mydata)), where = .GlobalEnv)) {
        return(invisible(FALSE)) #if (direct == F) {return(invisible(FALSE))} else {return(invisible(NULL))}
      } else {
        return(invisible(TRUE)) ##if (direct == F) {return(TRUE)} else {return(invisible(mydata))}
      }
    },
    warning = function(cond) {
      return(TRUE) #if (direct == F) {return(TRUE)} else {return(invisible(mydata))}
    }
  )

  to_return = res
  if (direct == T) {
    if (is_false(res)) {
      to_return = NULL
    } else {
      to_return = mydata
    }
  }

  return(to_return)
}

smdf = function(df) {
  for (i in 1:length(cn(df))) {
    mycn = cn(df)[i]
    cat(' - ',i,'. ',mycn,': ',head(unique(df[,cn(df) == mycn]),n=5),'\n',sep=' ')
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

window_mean = function(vals,window=NA,na.value=NA,na.rm=FALSE,debug=FALSE,verbose=FALSE,offset=0,stretch=FALSE,exact=TRUE,mode=1) {
  if (is.na(window)) {
    window = 200
  }
  window0 = window
  window = window - 1
  if (mode == 1) {
    exact=T
    offset=0
  } else if (mode == 2) {
    exact=F
    offset=-1*ai(window0/2)
  }

  if (!defined(vals)) {
    return(na.value)
  }
  if (window < 2) {
    window = 2
  }
  if (length(vals) < 2) {
    my.mean = mean(vals,na.rm=na.rm)
    if (is.na(my.mean)) {
      return(na.value)
    } else {
      return(my.mean)
    }
  }
  if (window >= length(vals)) {
    window = length(vals)-1
  }
  i0 = 1
  i1 = length(vals) - window
  offsets = c(seq(i0,i0+window-1),seq(i0+window,i1+window+window))#,seq(i1+1,i1+window))
  df0 = c(rep(i0,i0+window-1),seq(i0,i1),seq(i1+1,i1+window));df0
  df1 = c(seq(i0,i0+window-1),seq(i0+window,i1+window),rep(i1+window,window));df1
  df = data.frame(i0 = df0, i1 = df1)

  my.mean = apply(df,1,function(x,vals,na.rm,na.value){mean(vals[x[1]:x[2]],na.rm=na.rm,na.value=na.value)} ,vals=vals,na.rm=na.rm,na.value=na.value)
  #seq(i1+window,i1,-1)))
  dftemp = df
  dftemp$offsets = offsets-window-1
  dftemp$c0 = vals[dftemp$i0]
  dftemp$c1 = vals[dftemp$i1]
  dftemp$my.mean = my.mean 
  min.offset = (min(dftemp$offsets))
  if (offset < min.offset) {
    offset = min.offset
  } else if (offset > 0) {
    offset = 0
  }
  
  # to_return = dftemp[dftemp$offsets >= offset & dftemp$offsets <= offset + length(vals)-1,]$my.mean
  # if (exact == TRUE) {
  #   to_return = dftemp[dftemp$offsets >= offset & dftemp$offsets <= offset + length(seq(i0,i1))-1-1-1,]
  # }
  r0 = window0 + offset
  r1 = window0 + offset + length(vals)-1
  if (exact == TRUE) {
    r0 = window0 + offset+1
    r1 = r0 + length(seq(i0,i1)) - 1-1
    if (window0 %% 2 != 0) {
      r1 = r1 - 1
    }
  }
  to_return = my.mean[r0:r1]#(window0+offset):(length(vals)-1+offset)]
  if (stretch == TRUE) {
    to_return = stretch(to_return,length(vals))
  }
  if (verbose == TRUE) {
    print(head(vals))
    cat('\n- window',window0,'\n- i0',i0,'-',i0+window,'\n- i1',i1,'-',i1+window,'\n\n')
    print(head(dftemp))
    print(tail(dftemp))
    print(r0)
    print(r1)
  }
  if (debug == TRUE) {
    lay(4)
    set.seed(4200)
    vals = rnorm(1000)*1000
    window=20
    windowtemp = ai(window*3/3)
    vals.mean = window_mean(vals,window=windowtemp,verbose=F,mode=1);plot(main='mode 1',vals.mean,col=2,type='l',xlim=c(0,1000),ylim=c(-1000,1000))
    vals.mean = window_mean(vals,window=windowtemp,verbose=F,mode=2);lines(vals.mean,col=4)
    vals.mean = window_mean(vals,window=windowtemp,verbose=F,mode=1,stretch=T);plot(main='mode 1 stretch',vals.mean,col=2,type='l',xlim=c(0,1000),ylim=c(-1000,1000))
    vals.mean = window_mean(vals,window=windowtemp,verbose=F,mode=2);lines(vals.mean,col=4)
    vals.mean = window_mean(vals,window=windowtemp,verbose=F,mode=1);plot(main='mode 1, offset on plotting',seq(windowtemp/2+1,1000-windowtemp/2),vals.mean,col=2,type='l',xlim=c(0,1000),ylim=c(-1000,1000))
    vals.mean = window_mean(vals,window=windowtemp,verbose=F,mode=2);lines(vals.mean,col=4)
    vals.mean = window_mean(vals,window=windowtemp,verbose=F,mode=2);plot(main='mode 2',vals.mean,col=4,type='l',xlim=c(0,1000),ylim=c(-1000,1000))
    return(to_return)
  } else {
    #print(min.offset)
    return(to_return)
  }
}
# Example
# lay(4); vals.mean = window_mean(vals,verbose=TRUE,debug=TRUE)
# lay(4)
# vals.mean = window_mean(vals,window=windowtemp,debug=F,verbose=F,offset=-1*windowtemp/2+1,exact=T,stretch=T);plot(vals.mean,col=1,ylim=c(-200,1000),type='l')
# vals.mean = window_mean(vals,window=windowtemp,debug=F,verbose=F,offset=-1*windowtemp/1+1,exact=T,stretch=T);lines(vals.mean,col=2)
# vals.mean = window_mean(vals,window=windowtemp,debug=F,verbose=F,offset=-0*windowtemp/1+1,exact=T,stretch=T);lines(vals.mean,col=4)


