#' Title Batch download genbank/fasta file from URL
#'
#' @param acc accession number or a list of accession number
#' @param file_type the type of file you want to download
#' @param db_type the database you want to search that help narrow down the range
#'
#' @return specified type of file(fasta or genbank or other)
#' @export
#'
#' @examples download_genbank(accn,"gb","nucleotide")
download_genbank<-function(acc,file_type,db_type){
  for (i in 1:length(acc)){
    url=paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=",db_type,"&rettype=",file_type,"&retmode=text&id=",acc[i],sep="")
    dest=paste(acc[i],".",file_type,sep="")
    download.file(url,destfile =dest)}
}



#' Title File format conversion:genbank to fasta
#'
#' @param path the path of the required file
#'
#' @return fasta file
#' @export
#'
#' @examples gb2fas("AB115403.gb")
gb2fas<-function(path){
  x<-readLines(path)
  grep('^ACCESSION',x)->i
  id1=sub("ACCESSION(\\s*)","",x[i])
  id=paste('>',id1,sep = '')

  grep('^ORIGIN',x)->j
  grep('^//',x)->t
  seq=gsub("(\\s*)","",x[(j+1):(t-1)])
  seq=gsub("(\\d*)","",seq)

  outfile <- file(description = paste(id1,".fas",sep=''), open = 'w')
  writeLines(id,outfile)
  writeLines(seq,outfile)
  close(outfile)
}



#' Title Read fasta file
#'
#' @param path the path of the required file
#'
#' @return sequence in fasta file and description
#' @export
#'
#' @examples read_fasta("flu_seq.fas")
read_fasta<-function(path){
  x<-readLines(path)
  n=length(x)
  grep('^>',x)->i
  id=sub(">","",x[i])
  start=(i+1)
  end=c(i[-1]-1,n)
  ss<-sapply(seq_along(start), function(i) paste0(x[start[i]:end[i]],collapse = ''))
  lt<-list(id=id,seq=ss)
  return(lt)
}



#' Title Calculate the base frequency
#'
#' @param inputseq input sequence whose seperation is ''
#'
#' @return a table of each sequence's frequency
#' @export
#'
#' @examples base.freq("AGTCAGTTTAC")
base.freq<-function(inputseq){
  na=nc=nt=ng=0
  l<-toupper(unlist(strsplit(inputseq,"")))
  for (i in 1:length(l)){
    if (l[i]=='A')
      na=na+1
    else if (l[i]=='C')
      nc=nc+1
    else if (l[i]=='T')
      nt=nt+1
    else if (l[i]=='G')
      ng=ng+1
  }
  ratio_a=na/length(l)
  ratio_t=nt/length(l)
  ratio_c=nc/length(l)
  ratio_g=ng/length(l)
  df<-data.frame(a=ratio_a,t=ratio_t,c=ratio_c,g=ratio_g,base_num=length(l))
  return(df)
}



#' Title Parse fasta file(Calculate the base frequency)
#'
#' @param lst a list of id and seq
#'
#' @return res :a table of each sequence's frequency
#' @export
#'
#' @examples pparse_fasta(lt)
parse_fasta<-function(lst){
  res=t(sapply(lst$seq,base.freq,USE.NAMES = F))
  rownames(res)=lst$id
  colnames(res)=toupper(colnames(res))

  cat("total",length(lst$seq),"sequences")
  cat("Labels:\n")
  print(lst$id)
  cat("Base Composition:\n")
  print(res)
  return(res)
}



#' Title Simplified global alignment of 2 sequences
#'
#' @param seq1 input sequence1 for alignment whose seperation is ''(a character string)
#' @param seq2 input sequence2 for alignment whose seperation is ''(a character string)
#' @param score score system
#'
#' @return Align object
#' @export
#'
#' @examples Blast("agctattgca","agtcattgca"),Blast(seq1,seq2)
Blast<-function(seq1, seq2,score = list(match = 5, mismatch = -2, GAP = -6)){
  # simplified blast with scoring system(5,-2,-6)
  # ignoring affine penalty,the gap penalty is same.

  l1 = length(seq1)
  l2 = length(seq2)

  match <- score$match
  mismatch <- score$mismatch
  GAP <- score$GAP

  # initiate the score matrix
  scores <- matrix(NA, l1, l2)
  scores[,1] <- sapply(1:l1-1, function(x) x * GAP)
  scores[1,] <- sapply(1:l2-1, function(x) x * GAP)

  #using point matrix to record the trace path and initiate
  point <- matrix(NA, l1, l2)
  point[,1] <- sapply(1:l1-1, function(x) x='↑')   #seq1 align to '-' all
  point[1,] <- sapply(1:l2-1, function(x) x='←')   #seq2 align to '-' all
  point[1,1]=0

  for (i in 2:l1){                                   #fill in the score matrix
    for (j in 2: l2){
      if ( seq1[i] == seq2[j])
      {diagonal_score <- scores[i-1, j-1] + match}
      else
      {diagonal_score <- scores[i-1, j-1] + mismatch}
      left_score = GAP + scores[i,j - 1]
      up_score = GAP + scores[i - 1,j]
      max_score = max(diagonal_score, left_score, up_score)
      scores[i,j]=max_score

      #record the corresponding directions that obtain the best local score
      if (scores[i,j] == diagonal_score)
        point[i,j]='↖'
      else if (scores[i,j] == left_score)
        point[i,j]='←'
      else
        point[i,j]='↑'
    }
  }

  a1 <- character()   #Trace back
  a2 <- character()
  i = l1
  j = l2
  while (i > 1 && j >1){
    if (point[i,j] == '↖')     #trace back to diagonal
    {a1 <- c(seq1[i],a1)
    a2 <- c(seq2[j],a2)
    i <- i-1
    j <- j-1
    next}
    if (point[i,j] == '←'){     #trace back to left
      a1 <- c("-",a1)
      a2 <- c(seq2[j],a2)
      j <- j-1
      next}
    if (point[i,j] == '↑')       #trace back to upward
    {a1 <- c(seq1[i],a1)
    a2 <- c("-",a2)
    i <- i-1
    next}}

  gap_pos <- which(a1 == '-' | a2 == '-')
  a1 <- paste(a1, collapse='')
  a2 <- paste(a2, collapse='')

  structure(
    list(seq = c(paste(seq1,collapse = ''), paste(seq2,collapse = '')),
         aln = c(a1, a2),
         score = score,
         matrix = scores,
         patrix = point,
         gap_position = gap_pos
    ),
    class = "Align"
  )
}



##' @method print Align
##' @export
print.Align<-function(x,...){
  scores<-x$matrix
  point<-x$patrix
  cat ("Seq1: ", x$seq[1],"\n")
  cat ("Seq2: ", x$seq[2],"\n")
  cat ("Scoring system: ", x$score$match, " for match; ", x$score$mismatch, " for mismatch; ", x$score$GAP, " for gap", "\n\n")
  cat ("Dynamic programming matrix:\n")
  print(scores)
  cat("Trace path matrix:\n")
  print(point)
  cat ("\nAlignment:\n")
  cat (paste(x$aln[1], collapse=''), "\n")
  bar <- rep("|", nchar(x$aln[1]))
  bar[x$gap_position] <- " "
  cat(paste0(bar, collapse = ""), "\n")
  cat (paste(x$aln[2], collapse=''),"\n\n")
  cat ("Optimum alignment score:", scores[length(scores)],"\n")
}
