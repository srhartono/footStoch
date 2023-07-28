#finalbeg.save = finalbeg
#finalend.save = finalend
finalbeg = finalbeg.save
finalend = finalend.save
finalbegInd = 1
finalendInd = 1
final0 = data.frame()
breakbeg = 0
iter = 0
finalbeg = finalbeg[order(-finalbeg$x0,-finalbeg$x1,-finalbeg$y0),]
finalend = finalend[order(-finalend$y1,-finalend$y0,-finalend$x0),]
while (breakbeg == 0) {
  if (iter == 0) {
    finalbeg = finalbeg.save
    finalend = finalend.save
    finalbeg = finalbeg[order(-finalbeg$x0,-finalbeg$x1,-finalbeg$y0),]
    #finalbeg = finalbeg[order(finalbeg$x0,finalbeg$x1,finalbeg$y0),]
    finalend = finalend[order(-finalend$y1,-finalend$y0,-finalend$x0),]
    finalbegInd = 1
    finalendInd = 1
    final0 = data.frame()
    breakbeg = 0
  }
  iter = iter + 1
  if (iter > 200) {print('iter break!');break}
  if (finalbegInd > size(finalbeg)) {
    finalbeg = finalbeg[finalbeg$x0 != finalbeg$x1 | finalbeg$y0 != finalbeg$y1,]
    finalend = finalend[finalend$x0 != finalend$x1 | finalend$y0 != finalend$y1,]
    breakbeg = 0
    print('break 1')
    break
  }
  x0a = finalbeg[finalbegInd,]$x0 #3
  y0a = finalbeg[finalbegInd,]$y0 #8
  x1a = finalbeg[finalbegInd,]$x1 #49
  y1a = finalbeg[finalbegInd,]$y1 #8
  if (x0a == x1a & y0a == y1a) {
    breakbeg = 0
    finalbegInd = finalbegInd + 1
    print(paste(finalbegInd,': not using x0a == x1a and y0a == y1a!',x0a,y0a,x1a,y1a))
    next
  }

  breakend = 0
  while (breakend == 0) {
    breakend = 1
    if (finalendInd > size(finalend)) {
      breakend = 0
      finalbeg = finalbeg[finalbeg$x0 != finalbeg$x1 | finalbeg$y0 != finalbeg$y1,]
      finalend = finalend[finalend$x0 != finalend$x1 | finalend$y0 != finalend$y1,]
      break
    }
    x0b = finalend[finalendInd,]$x0 #5
    y0b = finalend[finalendInd,]$y0 #6
    x1b = finalend[finalendInd,]$x1 #5
    y1b = finalend[finalendInd,]$y1 #9
    if (x0b == x1b & y0b == y1b) {
      breakend = 0
      finalendInd = finalendInd + 1
      next
    }
#    print(paste(x0a,y0a,x1a,y1a))
#    print(paste(x0b,y0b,x1b,y1b))
    
    if (defined(myintersect(x0a-1,x1a+1,x0b-1,x1b+1)) & defined(myintersect(y0a-1,y1a+1,y0b-1,y1b+1))) {
      x0end = x0b #5
      y0end = y0b #6
      
      x1beg = x1a #49
      y1beg = y1a #8
      
      final0.temp = data.frame(x0=x0end,y0=y0end,x1=x1beg,y1=y1beg)
      final0 = rbind(final0,final0.temp)
      #3,8 - 49,8
      #5,6 - 5,9
      #5,6 - 49,8
      #3,8-5,8 x2a,
      #5,8-5,9
      finalbeg[finalbegInd,]$x0 = x0a #3
      finalbeg[finalbegInd,]$y0 = y0a #8
      finalbeg[finalbegInd,]$x1 = x1b #5
      finalbeg[finalbegInd,]$y1 = y1a #8

      finalend[finalendInd,]$x0 = x0b #5
      finalend[finalendInd,]$y0 = y0a #8
      finalend[finalendInd,]$x1 = x1b #5
      finalend[finalendInd,]$y1 = y1b #9

      print(paste(iter,',',finalendInd,'  (',x0a,',',y0a,') - (',x1a,',',y1a,')',' & ','(',x0b,',',y0b,') - (',x1b,',',y1b,'): ',x0end,' ',y0end,' ',x1beg,' ',y1beg,' ',sep=''))
      breakend = 1
    } else {
      breakend = 0
      print(paste(iter,',',finalendInd,'  Not Int! (',x0a,',',y0a,') - (',x1a,',',y1a,')',' & ','(',x0b,',',y0b,') - (',x1b,',',y1b,'): ',x0end,' ',y0end,' ',x1beg,' ',y1beg,' ',sep=''))
    }
    if (finalendInd == size(finalend)) {
      breakend = 0
      finalendInd = 1
      break
    }
    if (breakend == 1) {
      finalbegInd = 1
      finalendInd = 1
      print('breakend!')
      break
    }
    finalendInd = finalendInd + 1
  }
  
  if (finalbegInd == size(finalbeg)) {
    breakbeg = 0
    finalbeg = finalbeg[finalbeg$x0 != finalbeg$x1 | finalbeg$y0 != finalbeg$y1,]
    finalend = finalend[finalend$x0 != finalend$x1 | finalend$y0 != finalend$y1,]
    break
  }
  # if (breakbeg == 1) {
  #   finalbeg = finalbeg[finalbeg$x0 != finalbeg$x1 & finalbeg$y0 != finalbeg$y1,]
  #   finalend = finalend[finalend$x0 != finalend$x1 & finalend$y0 != finalend$y1,]
  #   break
  # }
  if (breakend == 0) {
    print(paste(finalbegInd,'beg: breakend 0 so finalbegInd + 1!'))
    finalbegInd = finalbegInd + 1
  } else {
    print(paste(finalbegInd,'/',size(finalbeg),'beg: breakend 1 so finalbegInd is 1 again!'))
  }
}
iter = 0

finalbeg
finalend
final0
p2 + geom_rect(data=final0,aes(xmin=x0,xmax=x1,ymin=y0,ymax=y1),fill=NA,color='black')

plot(NA,xlim=c(0,divby*100),ylim=c(0,divby*100))
points(dm2$beg,dm2$end,pch='.')
segments(0,0,2500,2500)#c(0,0),c(2500,2500))
rect(xleft = final0$x0*divby,ybottom = final0$y0*divby,xright = final0$x1*divby,ytop = (final0$y1+1)*divby,col = rgb(1,1,1,0),border=1)
#rect(xleft = finalend$x0,ybottom = final0$y0,xright = finalend$x1,ytop = finalend$y1,col = rgb(0,0,0,0),border=1)

myintersect = function(x0a,x1a,x0b,x1b) {
  
  test1 = data.frame(x=x0a:x1a)
  test2 = data.frame(x=x0b:x1b)

  #test1 = data.frame(x=c(1:3,10:20))
  #test2 = data.frame(x=c(2:5,12:14))
  test3 = intersect(test1$x,test2$x)
  test3
}
# {
#   test3.edge = findedge(test3)
#   if (defined(test3.edge[test3.edge != 0])) {
#     test3.edge = test3.edge[test3.edge != 0]
#   }
#   if (length(test3.edge) == 1) {
#     
#   }
#   test3.edge.beg = test3.edge[seq(1,length(test3.edge)-1)] + 1
#   test3.edge.end = test3.edge[seq(min(length(test3.edge),1),length(test3.edge))]
#   test3beg = test3[test3.edge.beg]
#   test3end = test3[test3.edge.end]
#   print(paste(test3beg,test3end))
# }
# defined(getinter(x0a,x1a,x0b,x1b))
# defined(getinter(0,0,x0b,x1b))
