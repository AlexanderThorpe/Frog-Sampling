
burned <- round(iter/2):iter ## Adjust for how much burn-in you want to discard.

## Get mean point estimates. Discard first half as burn-in.
tmp_f <- function(x) {
    ## Return some quantiles and also the mean.
    q <- quantile(x,prob=c(.25,.5,.75))
    m <- mean(x)
    return(c(q,mean=m))
}

est_f <- apply(f[,,burned],1:2,tmp_f)
est_d <- apply(d[,burned],1,tmp_f)
est_g <- apply(g[,burned],1,tmp_f)

par(mar=c(2,2,0,0),mfcol=c(4,3))
thin <- round(seq(from=1,to=iter,length.out=100))
## Chains for some random effects.
tmp<-110:121
for (i in 1:4) matplot(t(f[tmp,i,thin]),type="l",lty=1)
matplot(t(d[tmp,thin]),type="l",lty=1)
matplot(t(g[tmp,thin]),type="l",lty=1)

## Chains for group level parameters.
matplot(exp(t(group_level_parameters["loga",,thin])),type="l",lty=1)
matplot(exp(t(group_level_parameters["logb",,thin])),type="l",lty=1)

## Posterior predictive checks.
iter_pp <- round(seq(from=0.75*iter,to=iter,length.out=20))
data_pp <- array(dim=c(dim(data),length(iter_pp)))
for (i in 1:length(iter_pp)) {
    data_pp[,,,i] <- ll(f=f[,,iter_pp[i]],d=d[,iter_pp[i]],g=g[,iter_pp[i]],data=NULL,sample=TRUE)
    data_pp[,,,i][is.na(data)] <- NA ## Make same clips observed as real data.
}

for (i in 1:n_species) {
    plot(apply(data[,,i],1,mean,na.rm=T),apply(data_pp[,,i,],1,mean,na.rm=T),pch=16,xlim=c(0,1),ylim=c(0,1),col=expert[,i]+1)
    mtext(side=3,line=-2,paste(dimnames(data)[[3]][i],"clips"))
    abline(coef=c(0,1))
    plot(apply(data[,,i],2,mean,na.rm=T),apply(data_pp[,,i,],2,mean,na.rm=T),pch=16,xlim=c(0,1),ylim=c(0,1))
    mtext(side=3,line=-2,paste(dimnames(data)[[3]][i],"people"))
    abline(coef=c(0,1))
}


## Acceptance rates
print(mean(apply(d,1,diff)!=0)) ## d and g same.
apply((apply(f,1:2,diff)!=0),3,mean) ## Separate by species.
print( apply(apply(group_level_parameters,1:2,diff)!=0,2:3,mean)[1,])

## Posterior random effect estimates vs. expert.
for (i in 1:n_species) boxplot(est_f["mean",,i]~expert[,i],main=dimnames(data)[[3]][i])

## Check group distributions.
groupmeans <- (apply(exp(group_level_parameters[,,burned]),1:2,mean))
for (i in 1:ncol(groupmeans)) {
    if (i==(ncol(groupmeans)-1)) {
        tmp <- est_d ; nm<-"d"
    } else if (i==ncol(groupmeans)) {
        tmp <- est_g ; nm<-"g"
    } else {
        tmp <- est_f[,,i] ; nm<-dimnames(groupmeans)[[2]][i]
    }
    hist(tmp["mean",],prob=TRUE,main=nm,breaks="FD")
    tmp <- seq(0,1,.001)
    lines(x=tmp,y=dbeta(x=tmp,shape1=groupmeans[1,i],shape2=groupmeans[2,i]),col="red")
}

## Pond occupancy. How to address this? First, get mean of group f distribution:
print(groupmeans[1,]/colSums(groupmeans[,]),3) ## 2x bigger for GBF than LLJ

