param <- commandArgs(trailingOnly=T)

list_file = eval(paste(text=param[1]))
outfile   = eval(paste(text=param[2]))

if (length(param)!=2) {
  cat("\nError: two arguments must be supplied.\n")
  cat("ARG1 : the list of annot file (comma delimited)\n")
  cat("ARG2 : the outfile file.\n\n")
  stop("two arguments must be supplied.\n", call.=FALSE)
}

list_file = strsplit(list_file,",")[[1]]

for (chr in 1:22) {
	#read all files in list_files
    data = read.table(gzfile(paste(list_file[1],".",chr,".annot.gz",sep="")),h=T)
    if (length(list_file)>1){
        for (i in 2:length(list_file)){
            temp = read.table(gzfile(paste(list_file[i],".",chr,".annot.gz",sep="")),h=T)
            data = cbind(data,temp[,-(1:4)])
        }
    }
    sumannots_chr = NULL
    annot         = NULL
    for(i in 5:ncol(data)){
        #check if annotation is binary
        if ((sum(data[,i]==1)+sum(data[,i]==0)) == nrow(data)){
            annot         = c(annot,colnames(data)[i])
            sumannots_chr = rbind(sumannots_chr,
                apply(subset(data[,-(1:4)],data[,i]==1),2,sum))
        }
    } 
	if (chr==1) {
        sumannots = sumannots_chr
    } else { 
        sumannots = sumannots + sumannots_chr
    }
}

sumannots = cbind(annot,sumannots)

write.table(sumannots,file=outfile,col.names=T,row.names=F,quote=F,sep="\t")


