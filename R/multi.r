library(mclust)
data<-read.table("all.sp")
args <- commandArgs(trailingOnly = TRUE)
XXX <- args[1]
YYY <- args[2]

out <- paste(paste(as.character(XXX),'-',sep=''),YYY,sep='')

fit<-Mclust(data[,XXX:YYY],G=1:99,modelNames="VVV")
write.table(file=paste(out,".mod",sep=''),fit$parameters$mean)
write.table(file=paste(out,".pro",sep=''),fit$parameters$pro)
dat<-numeric()
for (i in 1:length(fit$parameters$pro)) {
    dat<-cbind(dat,sqrt(diag(fit$parameters$variance$sigma[,,i])))
}
write.table(file=paste(out,".var",sep=''),dat)
write.table(file=paste(out,".bic",sep=''),fit$BIC)
write.table(file=paste(out,".z",sep=''),  fit$z)

