dir.create("data")
setwd("data")
acc=paste("AB115",403:433,sep="")
acc1=paste("AJ5345",26:49,sep="")
acc
#can change along the fetch request
file_type="gb"  #gb or fasta
db_type="nucleotide"
list.files()

download_genbank<-function(acc,file_type,db_type)
  {
  for (i in 1:length(acc)){
  url=paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=",db_type,"&rettype=",file_type,"&retmode=text&id=",acc[i],sep="")
  dest=paste(acc[i],".",file_type,sep="")
  download.file(url,destfile =dest)}
}

download_genbank(acc,file_type,db_type)
list.files()