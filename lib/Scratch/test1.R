#FASTAS
set.seed(420)
test1 = data.frame(y=seq(1,5000),meanbeg=runif(5000)*5000); test1 = test1[order(test1$meanbeg),]; test1$y = seq(1,size(test1))
test2 = data.frame(y=seq(5001,10000),meanbeg=rnorm(5000)*5000); test2 = test2[order(test2$meanbeg),]; test2$y = seq(5001,10000)
plot(rbind(test1,test2)$y,rbind(test1,test2)$meanbeg,type='l')

test1 = data.frame(y=seq(1,5000),meanbeg=runif(5000)*5000); test1$y = seq(1,size(test1))
test2 = data.frame(y=seq(5001,10000),meanbeg=rnorm(5000)*5000); test2$y = seq(5001,10000)
plot(rbind(test1,test2)$y,rbind(test1,test2)$meanbeg,type='l')

myrange = list()
testTypes = c('unif','norm','hots')

#test1 = test1[order(test1$meanbeg),]; test1$y = seq(1,size(test1))
plot(rbind(test1,test2)$y,rbind(test1,test2)$meanbeg,type='l')

for (testTypeInd in seq(1,length(testTypes))) {
  testType = testTypes[testTypeInd]
  myrange$p = 0
  totest = test1
  plot(totest$y,totest$meanbeg,type='l')
  if (size(totest[,]) >= gp$dist.test_min.length) {
    if (testType == 'unif') {
      if (length(unique(totest[,][,'meanbeg'])) >= gp$dist.test_min.unique) {
        myrange$hist = hist(totest[,'meanbeg'],plot=F,breaks = max(2,sqrt(dim(totest[,][,])[1])))
        myrange$p = uniform.test(myrange$hist)$p.value
      }
    } else if(testType == 'norm') {
      if (length(unique(totest[,][,'meanbeg'])) >= gp$dist.test_min.unique) {
        myrange$p = shapiro.test(totest[,'meanbeg'])$p.value
      }
    } else if (testType == 'hots') {
      if (length(unique(totest[,][,'meanbeg'])) < gp$dist.test_max.unique.vhot) {
        myrange$p = 1
      }
    }
  }
  print(paste(testType,format(myrange$p,digits=3)))
}

for (testTypeInd in seq(1,length(testTypes))) {
  testType = testTypes[testTypeInd]
  myrange$p = 0
  totest = test2
  plot(totest$y,totest$meanbeg,type='l')
  
  if (size(totest[,]) >= gp$dist.test_min.length) {
    if (testType == 'unif') {
      if (length(unique(totest[,][,'meanbeg'])) >= gp$dist.test_min.unique) {
        myrange$hist = hist(totest[,'meanbeg'],plot=F,breaks = max(2,dim(totest[,][,])[1]/10))
        myrange$p = uniform.test(myrange$hist)$p.value
      }
    } else if(testType == 'norm') {
      if (length(unique(totest[,][,'meanbeg'])) >= gp$dist.test_min.unique) {
        myrange$p = shapiro.test(totest[,'meanbeg'])$p.value
      }
    } else if (testType == 'hots') {
      if (length(unique(totest[,][,'meanbeg'])) < gp$dist.test_max.unique.vhot) {
        myrange$p = 1
      }
    }
  }
  print(paste(testType,format(myrange$p,digits=3)))
}
