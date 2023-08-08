
plotthis = function(dfseq,dg.all.temp,x0=NA,x1=NA,window) {
  if (is.na(x0)) {
    dg.all.temp.x = seq(par$window/2,par$window/2+size(dg.all.temp)-1)
  } else {
    dg.all.temp.x = seq(x0,x1)
  }
  print(length(dg.all.temp.x))
  print(length(dg.all.temp$dgmin00))
  if (defined(dg.all.temp[is.infinite(dg.all.temp$dgmin00),])) {
    dg.all.temp[is.infinite(dg.all.temp$probmin00),]$dgmin00 = 1
  }

  if (defined(dg.all.temp[is.infinite(dg.all.temp$dgmin07),])) {
    dg.all.temp[is.infinite(dg.all.temp$probmin07),]$dgmin07 = 1
  }
  if (defined(dg.all.temp[is.infinite(dg.all.temp$dgmin14),])) {
    dg.all.temp[is.infinite(dg.all.temp$probmin14),]$dgmin14 = 1
  }
  # test = dg.all.temp$probmin07
  # test[is.na(test)]=1
  # test[is.infinite(test)] = 1
  # print(min(test))
  # print(max(test))
  #dg.all.temp.x = seq(par$window/2,par$window/2+length(dg.all.temp)-1)
  minval = min(c(dg.all.temp$dgmin00/500,dg.all.temp$dgmin07/500,dg.all.temp$dgmin14/500,log((dg.all.temp$probmin00)+500)/500,log((dg.all.temp$probmin00)+500)/500,log((dg.all.temp$probmin00)+500)/500))
  maxval = max(c(dg.all.temp$dgmin00/500,dg.all.temp$dgmin07/500,dg.all.temp$dgmin14/500,log((dg.all.temp$probmin00)+500)/500,log((dg.all.temp$probmin00)+500)/500,log((dg.all.temp$probmin00)+500)/500))
  plot(dg.all.temp.x,dg.all.temp$dgmin14/500,ylim=c(minval,maxval),xlim=c(0,max(dfseq$i0)+par$window/2),type='l',col='purple4')
  lines(dg.all.temp.x,dg.all.temp$dgmin00/500,col='purple')
  lines(dg.all.temp.x,dg.all.temp$dgmin07/500,col='purple3')
  lines(dg.all.temp.x,(log(dg.all.temp$probmin00)+500)/500,col='purple',lty=2)
  lines(dg.all.temp.x,(log(dg.all.temp$probmin07)+500)/500,col='purple3',lty=2)
  lines(dg.all.temp.x,(log(dg.all.temp$probmin14)+500)/500,col='purple3',lty=2)
  lines(dfseq$i0,dfseq$CGdens,col='blue4')
  lines(dfseq$i0,dfseq$GCcont,col='green4')
  lines(dfseq$i0,dfseq$GCskew+0.5,col='red4')
#  lines(seq(par$window/2,par$window/2+length(dg.all.temp)-1),(dg.all.temp-50)/20,col='purple')
}
plot(NA,ylim=c(0,1.5),xlim=c(0,size(dfseq2)))
lines(dfseq2$i0,dfseq2$CGdens,col='blue4')
lines(dfseq2$i0,dfseq2$GCcont,col='green4')
lines(dfseq2$i0,dfseq2$GCskew+0.5,col='red4')

plot(NA,ylim=c(0,1.5),xlim=c(0,size(dfseqT)))
lines(dfseqT$i0,dfseqT$CGdens,col='blue4')
lines(dfseqT$i0,dfseqT$GCcont,col='green4')
lines(dfseqT$i0,dfseqT$GCskew+0.5,col='red4')

plot(NA,ylim=c(0,1.5),xlim=c(0,size(dfseqG)))
lines(dfseqG$i0,dfseqG$CGdens,col='blue4')
lines(dfseqG$i0,dfseqG$GCcont,col='green4')
lines(dfseqG$i0,dfseqG$GCskew+0.5,col='red4')

plot(NA,ylim=c(0,1.5),xlim=c(0,size(dfseq5)))
lines(dfseq5$i0,dfseq5$CGdens/3,col='blue4')
lines(dfseq5$i0,dfseq5$GCcont,col='green4')
lines(dfseq5$i0,dfseq5$GCskew+0.5,col='red4')

plot(NA,ylim=c(0,1.5),xlim=c(0,size(dfseq6)))
lines(dfseq6$i0,dfseq6$CGdens/3,col='blue4')
lines(dfseq6$i0,dfseq6$GCcont,col='green4')
lines(dfseq6$i0,dfseq6$GCskew+0.5,col='red4')