#Code from Bruno Santos
tab<-read.table("")
tab1 <- apply(tab,1,function(x) paste(">",x[1],"count=",x[2],"\n",x[1],sep=""))
 write.table(tab1,file=".fa",quote=FALSE,row.names=FALSE) 