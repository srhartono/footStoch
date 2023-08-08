# generate dummy sequence
dummy.gen_seq = function(df) {
  my.seq1 = join(rep('CGCGAAAAAA',10))
  my.seq2 = join(rep('AGAGAAAAAA',10))
  df = data.frame(def=c('def1','def2'),seq=c(my.seq1,my.seq2))
}

df.split = function(df,by,by.split=NA) {
  if (is.na(by.split)) {
    by.split = pasta(by,'.split')
  }
  df[,by.split] = sapply(df[,by],split,'',USE.NAMES = FALSE)
  return(df)
}


