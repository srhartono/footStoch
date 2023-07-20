mydir = "/Users/mitochy/Work/Project/Ethan//triclust/"
#dm1 = readRDS(file=paste(mydir,"dm1.sgRNA.all.RDS",sep=''))
dm1 = readRDS(file=paste(mydir,"dm1.all.RDS",sep=''))
RDS = dir(paste(mydir,"./final3/",sep=''),"*final3.RDS")
RDS = paste(mydir,"./final3/",RDS,sep='')
#RDS = paste("./final3/",RDS,sep="")
for (i in 1:length(RDS)) {
  temp = readRDS(RDS[i])
  print(i)
  if (i == 1) {
    final0 = temp
  } else {
    final0 = rbind(final0,temp)
  }
}
head(final0)
head(dm1)
saveRDS(dm1,file="dm1.all.RDS")
saveRDS(final0,file='final3.all.RDS')
