read_fasta<-function(path){
x<-readLines(path)
n=length(x)
grep('^>',x)->i
id=sub(">","",x[i])
start=(i+1)
end=c(i[-1]-1,n)
ss<-sapply(seq_along(start), function(i) paste0(x[start[i]:end[i]],collapse = ''))
df<-data.frame(id=id,seq=ss)
#return(ss)

res=t(sapply(ss,base.freq,USE.NAMES = F))
rownames(res)=df$id
colnames(res)=toupper(colnames(res))

cat("total",length(ss),"sequences")
cat("Labels:\n")
print(df$id)
cat("Base Composition:\n")
print(res)
}

base.freq<-function(inputseq){
  na=nc=nt=ng=0
  l<-unlist(strsplit(inputseq,""))
  for (i in 1:length(l)){
    if (l[i]=='A')
      na=na+1
    else if (l[i]=='C')
      nc=nc+1
    else if (l[i]=='T')
      nt=nt+1
    else
      ng=ng+1
  }
  ratio_a=na/length(l)
  ratio_t=nt/length(l)
  ratio_c=nc/length(l)
  ratio_g=ng/length(l)
  df<-data.frame(a=ratio_a,t=ratio_t,c=ratio_c,g=ratio_g,base_num=length(l))
  return(df)
}

read_fasta("flu_seq.fas")
read_fasta("flu_seq_v2.fas")