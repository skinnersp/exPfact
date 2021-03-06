library(mclust)
data<-read.table("singlep/all.sp")
fit<-Mclust(data[,XXX:YYY],G=1:99,modelNames="EII")
write.table(file="tmp.mod",fit$parameters$mean)
write.table(file="tmp.pro",fit$parameters$pro)
dat<-numeric()
for (i in 1:length(fit$parameters$pro))
dat<-cbind(dat,sqrt(diag(fit$parameters$variance$sigma[,,i])))
write.table(file="tmp.var",dat)
#write.table(file="tmp.var",fit$parameters$variance)
write.table(file="tmp.bic",fit$BIC)
write.table(file="tmp.z",fit$z)

