if (FALSE) {
  rm(list=ls())
  load("tree.RData")
  ##load("ground.RData")
  burned <- round(iter/2):iter ## Adjust for how much burn-in you want to discard.
  est_f <- apply(f[,,burned],1:2,mean)
}

## Construct ROC for classification accuracy vs. expert. Does not work
## for species LLJ and WTF as there is no signal, according to expert.
tmp <- seq(0,1,.01) ## Thresholds to test.
ROC <- array(dim=c(length(tmp),2,n_species),dimnames=list(paste0(tmp),c("fa","hit"),dimnames(data)[[3]]))
for (species in 1:n_species) {
    if (dimnames(data)[[3]][species]=="LLJ") next ## Skip LLJ
    if (dimnames(data)[[3]][species]=="WTF") next ## Skip WTF
    for (i in 1:length(tmp)) {
        tmp1 <- table(c(est_f[,species]>tmp[i],T,T,F,F) , c(expert[,species],T,F,T,F))
    ## I added one of each sort and subtracted. Helps maintain shape when there are zero entries in some cells.
        tmp1 <- tmp1-1
        ROC[i,,species] <- tmp1["TRUE",]/colSums(tmp1)
    }
}

if (frogs=="tree") {
    pdf(file="ROCs-tree.pdf",width=6.6,height=3)
    layout(mat=rbind(0,c(0:3,0),0),widths=c(.3,rep(1,3),.05),height=c(.05,1,.3))
} else {
    pdf(file="ROCs-ground.pdf",width=8.5,height=3)
    layout(mat=rbind(0,c(0:4,0),0),widths=c(.3,rep(1,4),.05),height=c(.05,1,.3))
}



par(mar=rep(.5,4))
for (species in 1:n_species) {
    if (dimnames(data)[[3]][species]=="LLJ") next ## Skip LLJ
    if (dimnames(data)[[3]][species]=="WTF") next ## Skip WTF
    plot(ROC[,,species],type="l",lwd=2,axes=F,xlab="",ylab="")
    abline(coef=c(0,1),lty=3)
    axis(side=1,at=seq(0,1,.25),labels=paste0(seq(0,100,25),"%"))
    axis(side=2,at=seq(0,1,.25),labels=FALSE)
    box()
    if (species==3) {
        if (frogs=="tree") {
            mtext(side=1,line=3,"Consensus 'Yes', Expert 'No'")
        } else {
            mtext(side=1,line=3,"Consensus 'Yes', Expert 'No'",at=0)
        }
        
    }
    if (species==1) {
        mtext(side=2,line=3,"Consensus 'Yes', Expert 'Yes'")
        axis(side=2,at=seq(0,1,.25),labels=paste0(seq(0,100,25),"%"))    
    }
    text(x=.7,y=.2,dimnames(data)[[3]][species],cex=2)
    ind <- which.min(abs(tmp-0.675))
##    points(x=ROC[ind,1,species],y=ROC[ind,2,species],col=2,pch=16,cex=2)
}

dev.off()
