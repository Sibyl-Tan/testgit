gb2fas<-function(path){
  x<-readLines(path)
  grep('^ACCESSION',x)->i
  id1=sub("ACCESSION(\\s*)","",x[i])
  id=paste('>',id1,sep = '')
  
  grep('^ORIGIN',x)->j
  grep('^//',x)->t
  seq=gsub("(\\s*)","",x[(j+1):(t-1)])
  seq=gsub("(\\d*)","",seq)
  
  outfile <- file(description = paste(id1,".fas"), open = 'w')
  writeLines(id,outfile)
  writeLines(seq,outfile)
  close(outfile)
  }

gb2fas("AB115403.gb")