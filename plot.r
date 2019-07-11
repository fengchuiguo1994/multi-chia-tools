args=commandArgs(T)
pdf(args[2])
mydat<-read.table(args[1])
names(mydat)<-c("dis",'pet')
plot(c(2.5,8),c(0,1), type = "n", ylab = "distance",xlab = "",xaxt="n")
for(i in 2:15){
    tmp<-mydat[(mydat$pet==i),]
    dd<-density(log10(tmp$dis+1))
    lines(dd, col = rainbow(15)[i], lwd = 2)
}
axis(1,at=c(1,2,3,4,5,6,7,8),labels = c("10","100","1k","10k","100k","1M","10M",'100M'))
dd<-density(log10(mydat$dis+1))
lines(dd, col = 'black', lwd = 2)
id<-paste(2:15)
id[15]<-'all'
colo<-rainbow(15)[2:15]
colo[15]="black"
legend("topright",id,lty=array(1,c(2,2)),col=colo,lwd=2)

plot(c(3,8),c(0,1), type = "n", ylab = "distance",xlab = "",xaxt="n")
for(i in 2:15){
    tmp<-mydat[(mydat$pet>=i),]
    dd<-density(log10(tmp$dis+1))
    lines(dd, col = rainbow(15)[i], lwd = 2)
}
axis(1,at=c(1,2,3,4,5,6,7,8),labels = c("10","100","1k","10k","100k","1M","10M",'100M'))
dd<-density(log10(mydat$dis+1))
lines(dd, col = 'black', lwd = 2)
id<-paste(2:15)
id[15]<-'all'
colo<-rainbow(15)[2:15]
colo[15]="black"
legend("topright",id,lty=array(1,c(2,2)),col=colo,lwd=2)

dev.off()
