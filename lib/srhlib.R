library(rlang)
e = exists
cn = colnames
rn = rownames
af = as.factor
ac = as.character
an = as.numeric
ai = as.integer
sm = summary

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
      if (!any(is.na(mydata)) & !any(is.null(mydata)) & length(mydata)> 0) { # if there is any real value then return TRUE
  # if (length(mydata) == 0) {
        #   FALSE
        # } else {
          TRUE
        # }
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

#mydir = "/Users/mitochy/Work/Project/Ethan//triclust/"

triclust_get_divby = function(dm4=data.frame(),mytitle=NA,gp=list()) {
  if (defined(gp$mytitle)) {
    mytitle = gp$mytitle
  }
  
  gp$triclust.divby = 50
  gp$triclust.thres1 = max(as.integer(sqrt(dim(dm4)[1]/9230)*25),5)
  gp$triclust.thres2 = max(as.integer(sqrt(dim(dm4)[1]/9230)*25),5)
  gp$triclust.sample = mytitle
  if (defined(grep("T7",mytitle,ignore.case=TRUE,perl=TRUE))) {
    #T7: 25, 100, 100
    gp$triclust.divby = 25
    gp$triclust.thres1 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.thres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.sample = 'T7_init'
  } else if (defined(grep("PFC53.+T3Term",mytitle,ignore.case = TRUE,perl=TRUE)) ) {
    #pFC53_T3Term: 15,max(50 etc)
    gp$triclust.divby = 15
    gp$triclust.thres1 = max(50,as.integer(sqrt(dim(dm4)[1]/9230)*100))
    gp$triclust.thres2 = max(50,as.integer(sqrt(dim(dm4)[1]/9230)*100))
    gp$triclust.sample = 'pFC53_T3Term'
  } else if (defined(grep("pFC53.+ApaLI",mytitle,ignore.case = TRUE,perl=TRUE)) ) {
    #pFC53_ApaLI: 25,as.integer
    gp$triclust.divby = 25
    gp$triclust.thres1 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.thres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.sample = 'pFC53_SSB_ApaLI'
  } else if (defined(grep("pFC9.+ApaLI",mytitle,ignore.case = TRUE,perl=TRUE)) ) {
    #pFC9_ApaLI: 25,as.integer
    gp$triclust.divby = 25
    gp$triclust.thres1 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.thres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.sample = 'pFC9_SSB_ApaLI'
  } else if (defined(grep("pFC9",mytitle,ignore.case = TRUE,perl=TRUE)) ) {
    #pFC9_sgRNA: 15,as.integer
    gp$triclust.divby = 50
    gp$triclust.thres1 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.thres2 = as.integer(sqrt(dim(dm4)[1]/9230)*100)
    gp$triclust.sample = 'pFC9_sgRNA'
  }
  
#  print(paste(mytitle,gp$triclust.divby,gp$triclust.thres1,gp$triculust.thres2))
  to_return = gp
  return(to_return)
}

