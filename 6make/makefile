#指定make命令依据makefile文件中的规则，进行构建
#make -f makefile

all:download_genbank.Rout gb2fas.Rout blast.Rout cal_dis.Rout

#makefile默认的target是文件，在运行时，会首先检测当前目录下是否存在这个文件,
#如果存在，则会abort掉make 进程，但目标不是文件的话，则会出现up to date的情况。
#解决的办法就是makefile中的“.PHONY”
.PHONY:download_genbank.Rout
.PHONY:gb2fas.Rout
.PHONY:blast.Rout
.PHONY:cal_dis.Rout

download_genbank.Rout:download_genbank.R
	R CMD BATCH download_genbank.R
  
gb2fas.Rout:download_genbank.Rout gb2fas.R
	R CMD BATCH gb2fas.R
  
blast.Rout:gb2fas.Rout blast.R 
	R CMD BATCH blast.R 
  
cal_dis.Rout:gb2fas.Rout cal_dis.R
	R CMD BATCH cal_dis.R