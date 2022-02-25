## Load some samples eg:
##load("samples-GBH-2ht.RData")

burned <- round(iter/2):iter ## Adjust for how much burn-in you want to discard.

## Get mean point estimates. Discard first half as burn-in.
est_f_GBF <- (apply(f_GBF[,burned],1,quantile,prob=c(.25,.5,.75)))
est_f_LLJ <- (apply(f_LLJ[,burned],1,quantile,prob=c(.25,.5,.75)))
est_d <- (apply(d[,burned],1,quantile,prob=c(.25,.5,.75)))
est_g <- (apply(g[,burned],1,quantile,prob=c(.25,.5,.75)))
est_f_GBF <- rbind(est_f_GBF,mean=apply(f_GBF[,burned],1,mean))
est_f_LLJ <- rbind(est_f_LLJ,mean=apply(f_LLJ[,burned],1,mean))

par(mar=c(2,2,0,0),mfcol=c(4,3))
thin <- round(seq(from=1,to=iter,length.out=100))
## Chains for some random effects.
##matplot(t(f_GBF[1:10,thin]),type="l",lty=1)
##matplot(t(f_LLJ[1:10,thin]),type="l",lty=1)
##matplot(t(d[1:10,thin]),type="l",lty=1)
##matplot(t(g[1:10,thin]),type="l",lty=1)

## Chains for group level parameters.
matplot(exp(t(group_level_parameters["loga",,thin])),type="l",lty=1)
matplot(exp(t(group_level_parameters["logb",,thin])),type="l",lty=1)

## Posterior predictive checks.
iter_pp <- round(seq(from=0.75*iter,to=iter,length.out=20))
data_GBF_pp <- array(dim=c(dim(data_GBF),length(iter_pp)))
data_LLJ_pp <- array(dim=c(dim(data_LLJ),length(iter_pp)))
for (i in 1:length(iter_pp)) {
    tmp <- ll(f_GBF=f_GBF[,iter_pp[i]],f_LLJ=f_LLJ[,iter_pp[i]],d=d[,iter_pp[i]],g=g[,iter_pp[i]],data_GBF=NULL,data_LLJ=NULL,sample=TRUE,model=model)
    tmp$GBF[is.na(data_GBF)] <- NA ## Make same clips observed as real data.
    tmp$LLJ[is.na(data_LLJ)] <- NA ## Make same clips observed as real data.
    data_GBF_pp[,,i] <- tmp$GBF
    data_LLJ_pp[,,i] <- tmp$LLJ
}

plot(apply(data_GBF,1,mean,na.rm=T),apply(data_GBF_pp,1,mean,na.rm=T),pch=16,xlim=c(0,1),ylim=c(0,1),col=expert_GBF+1)
abline(coef=c(0,1))
plot(apply(data_GBF,2,mean,na.rm=T),apply(data_GBF_pp,2,mean,na.rm=T),pch=16,xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))

plot(apply(data_LLJ,1,mean,na.rm=T),apply(data_LLJ_pp,1,mean,na.rm=T),pch=16,xlim=c(0,1),ylim=c(0,1),col=expert_LLJ+1)
abline(coef=c(0,1))
plot(apply(data_LLJ,2,mean,na.rm=T),apply(data_LLJ_pp,2,mean,na.rm=T),pch=16,xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))

##plot(x=est_f_GBF["50%",],y=est_f_LLJ["50%",],pch=16)


## Acceptance rates
print(mean(apply(d,1,diff)!=0))
print(mean(apply(f_GBF,1,diff)!=0))
print(mean(apply(f_LLJ,1,diff)!=0))
print( apply(apply(group_level_parameters,1:2,diff)!=0,2:3,mean))


## Only for synthetic data.
##plot(x=true_parameters$f_GBF,y=est_f_GBF["50%",],pch=16);abline(coef=c(0,1),col=2)
##plot(x=true_parameters$f_LLJ,y=est_f_LLJ["50%",],pch=16);abline(coef=c(0,1),col=2)
##plot(x=true_parameters$d,y=est_d["50%",],pch=16);abline(coef=c(0,1),col=2)
##plot(x=true_parameters$g,y=est_g["50%",],pch=16);abline(coef=c(0,1),col=2)

## Data row means vs. expert. 
boxplot(apply(data_GBF,1,mean,na.rm=T)~expert_GBF)
##boxplot(apply(data_LLJ,1,mean,na.rm=T)~expert_LLJ)
## Seems like there is signal in there!

## Posterior random effect estimates vs. expert.
boxplot(est_f_GBF["50%",]~expert_GBF)
##boxplot(est_f_LLJ["50%",]~expert_LLJ)


## Check group distributions.

groupmeans <- (apply(exp(group_level_parameters[,,burned]),1:2,mean))
for (i in 1:4) {
    hist(switch(i,est_f_GBF,est_f_LLJ,est_d,est_g)["50%",],prob=TRUE,main="",breaks="FD")
    tmp <- seq(0,1,.001)
    lines(x=tmp,y=dbeta(x=tmp,shape1=groupmeans[1,i],shape2=groupmeans[2,i]),col="red")
}


## Construct ROC for classification accuracy vs. expert. Only works
## for GBF as there is no signal in LLJ, according to expert.
tmp <- seq(0,1,.01) ## Thresholds to test.
ROC_GBF <- array(dim=c(length(tmp),2,2),dimnames=list(paste0(tmp),c("fa","hit"),c("mean","median")))
for (i in 1:length(tmp)) {
    for (j in 1:2) {
        tmp1 <- table(c(est_f_GBF[switch(j,"mean","50%"),]>tmp[i],T,T,F,F) , c(expert_GBF,T,F,T,F))
        ## I added one of each sort and subtracted. Helps maintain shape when there are zero entries in some cells.
        tmp1 <- tmp1-1
        ROC_GBF[i,,j] <- tmp1["TRUE",]/colSums(tmp1)
    }
}
## Plot the ROC. Seems like best discrimination is at threhsold near p=1.
plot(ROC_GBF[,,1],type="l")
lines(ROC_GBF[,,2],col=2)
abline(coef=c(0,1),lty=3)


## Pond occupancy. How to address this? First, get mean of group f distribution:
print(groupmeans[1,1:2]/colSums(groupmeans[,1:2]),3) ## 2x bigger for GBF than LLJ

## Have a look at counting the number of random effects that are biggish.
print(table(est_f_GBF["mean",]>.75))
print(table(est_f_LLJ["mean",]>.75))

