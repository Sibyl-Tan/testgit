##install package"treeio"
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("treeio")

##for using function"read.fasta",package "seqinr" should be loaded
install.packages("seqinr")
install.packages("pheatmap")
library(seqinr)
x=read.fasta("flu_seq.fas")

##a simple function used to calculate the base frequency
base.freq<-function(inputseq){
  count(inputseq,1,freq=T)
}

res=t(sapply(x,base.freq))
rownames(res)=sub("^(\\w+)\\|.*","\\1",names(x))
colnames(res)=toupper(colnames(res))
library(pheatmap)
pheatmap::pheatmap(res)

##for reading flu_seq_v2,package"Biostrings" should be loaded.
BiocManager::install("Biostrings")
require("Biostrings")
y<-read.fasta("flu_seq_v2.fas")
res1=t(sapply(y,base.freq))
rownames(res1)=sub("^(\\w+)\\|.*","\\1",names(y))
colnames(res1)=toupper(colnames(res1))
pheatmap::pheatmap(res1)
