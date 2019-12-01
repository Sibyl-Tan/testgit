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

#为展示美观，这里分别取两条序列的前9个
seq1<-toupper(c(0,unlist(strsplit(x,""))))[1:10]
seq2<-toupper(c(0,unlist(strsplit(y,""))))[1:10]

#seq1="agtc"
#seq2="tgc"
#seq1<-toupper(c(0,unlist(strsplit(seq1,""))))
#seq2<-toupper(c(0,unlist(strsplit(seq2,""))))

Blast<-function(seq1, seq2){  # 简化的blast
l1 = length(seq1)
l2 = length(seq2)
match=5
mismatch=-2
GAP = -6  # 这里不考虑仿射罚分，任何空位皆罚4分

# 初始化得分矩阵  
scores <- matrix(NA, l1, l2)
scores[,1] <- sapply(1:l1-1, function(x) x * GAP)
scores[1,] <- sapply(1:l2-1, function(x) x * GAP)

#用point矩阵记录回溯路径，并初始化
point <- matrix(NA, l1, l2)    
point[,1] <- sapply(1:l1-1, function(x) x='↑')   #seq1一直匹配空位，往上回溯
point[1,] <- sapply(1:l2-1, function(x) x='←')   #seq2一直匹配空位，往左回溯
point[1,1]=0

  for (i in 2:l1){                                   #填充得分矩阵
      for (j in 2: l2){
        if ( seq1[i] == seq2[j]) 
          {diagonal_score <- scores[i-1, j-1] + match}
        else 
          {diagonal_score <- scores[i-1, j-1] + mismatch}
        left_score = GAP + scores[i,j - 1]
        up_score = GAP + scores[i - 1,j]
        max_score = max(diagonal_score, left_score, up_score)
        scores[i,j]=max_score

        if (scores[i,j] == diagonal_score)   #记录取得局部最大得分对应的方向
          point[i,j]='↖'
        else if (scores[i,j] == left_score)
          point[i,j]='←'
        else
          point[i,j]='↑'
      }
  }

a1 <- character()   #动态规划法回溯
a2 <- character()
i = l1
j = l2
while (i > 1 && j >1){
#  if (point[i,j] == 0)    #异常值处理，break后重新进入while循环
#  break
 if (point[i,j] == '↖')     #往左上回溯
  {a1 <- c(seq1[i],a1)
  a2 <- c(seq2[j],a2)
  i <- i-1
  j <- j-1
  next}
 if (point[i,j] == '←'){     #往左回溯
  a1 <- c("-",a1)
  a2 <- c(seq2[j],a2)
  j <- j-1
  next}
 if (point[i,j] == '↑')                     #往上回溯
  {a1 <- c(seq1[i],a1)
  a2 <- c("-",a2)
  i <- i-1
  next}}

cat ("Seq1: ", seq1,"\n")
cat ("Seq2: ", seq2,"\n")
cat ("Scoring system: ", match, " for match; ", mismatch, " for mismatch; ", GAP, " for gap", "\n\n")
cat ("Dynamic programming matrix:\n")
print(scores)
cat("Trace path matrix:\n")
print(point)
cat ("\nAlignment:\n")
cat (paste(a1, collapse=''), "\n")
cat (paste(a2, collapse=''),"\n\n")
cat ("Optimum alignment score:", scores[l1,l2],"\n")
}

Blast(seq1, seq2)