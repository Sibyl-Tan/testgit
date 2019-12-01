dis=0

read_fasta<-function(path){
  x<-readLines(path)
  n=length(x)
  grep('^>',x)->i
  id=sub(">","",x[i])
  start=(i+1)
  end=c(i[-1]-1,n)
  ss<-sapply(seq_along(start), function(i) paste0(x[start[i]:end[i]],collapse = ''))
  df<-data.frame(id=id,seq=ss)
  return(ss)}

x<-read_fasta("AB115403.fas")
y<-read_fasta("AB115404.fas")

seq1<-toupper(c(0,unlist(strsplit(x,""))))
seq2<-toupper(c(0,unlist(strsplit(y,""))))

cal_dis<-function(seq1,seq2){
  if (length(seq1)!=length(seq2)){
    print("ERROR:the length of sequences is not equal!")
  }
  else {
    for (i in 1:length(seq1))
      a<-ifelse(seq1[i]==seq2[i],0,1)
    dis=dis+a
    dis
  }
}

cal_dis(seq1,seq2)