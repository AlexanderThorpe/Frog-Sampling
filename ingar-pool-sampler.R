rm(list=ls())

iter= 1e5 ## MCMC length
load("data/processed_data.RData")
model <- "2ht" ## Switch model here. 2HT seems better than SDT.
P <- nrow(data_GBF)
S <- ncol(data_GBF)

ll <- function(f_GBF,f_LLJ,d,g,data_GBF,data_LLJ,sample=FALSE,model) {
    ## Calculate the probabilites of TRUE responses and log-likelihoods.
    ## Uses either model == "sdt" or "2ht".

    if (model=="2ht") {
        ## MPT. Tree starts with frog (yes/no, f), then detect (yes/no, d), then guess
        ## (yes/no, g) if detect fails.
        all_f_GBF <- array(f_GBF,dim=c(length(f_GBF),length(d)))    ## Column-major prob GBF target.
        all_f_LLJ <- array(f_LLJ,dim=c(length(f_LLJ),length(d)))    ## Column-major prob LLJ target.

        all_d <- t(array(d,dim=c(length(d),length(f_GBF)))) ## Row-major prob f.
        all_g <- t(array(g,dim=c(length(d),length(f_GBF)))) ## Row-major prob f.
        probs_GBF <- all_f_GBF*all_d + all_g*(1-all_d) ## Probability of "yes" to GBF question.
        probs_LLJ <- all_f_LLJ*all_d + all_g*(1-all_d) ## Probability of "yes" to LLJ question. 
    }
    if (model=="sdt") {
        audio_f_GBF <- f_GBF ## Actual frog probability of presence GBF.
        audio_f_LLJ <- f_LLJ ## Actual frog probability of presence LLJ.
        person_prob_fa <- pnorm(mean=0,sd=1,q=g,lower.tail=FALSE) ## Probability of "yes" given no frog.
        person_prob_hit<- pnorm(mean=d,sd=1,q=g,lower.tail=FALSE) ## Probability of "yes" given frog.
        probs_GBF <- audio_f_GBF %*% t(person_prob_hit) + (1-audio_f_GBF) %*% t(person_prob_fa) ## Probability of "yes" to GBF.
        probs_LLJ <- audio_f_LLJ %*% t(person_prob_hit) + (1-audio_f_LLJ) %*% t(person_prob_fa) ## Probability of "yes" to LLJ.

    }        
        
    if (sample) return(list(GBF=probs_GBF,LLJ=probs_LLJ)) ## used in sampling synthetic data.
    ll_GBF <- log(data_GBF*probs_GBF + (!data_GBF)*(1-probs_GBF)) ## log likelihoods
    ll_LLJ <- log(data_LLJ*probs_LLJ + (!data_LLJ)*(1-probs_LLJ)) ## log likelihoods
    return(list(
        GBF=list(p=apply(ll_GBF,1,sum,na.rm=TRUE),s=apply(ll_GBF,2,sum,na.rm=TRUE)),
        LLJ=list(p=apply(ll_LLJ,1,sum,na.rm=TRUE),s=apply(ll_LLJ,2,sum,na.rm=TRUE))
    ))
}


## ## Uncomment this section to use synthetic data, for testing.
## true_parameters <- list(
##     f_GBF=rbeta(n=P,shape1=1,shape2=1),
##     f_LLJ=rbeta(n=P,shape1=4,shape2=8),
##     d=rbeta(n=S,shape1=3,shape2=4),
##     g=rbeta(n=S,shape1=1,shape2=2))
## tmp <- ll(f_GBF = true_parameters$f_GBF,f_LLJ = true_parameters$f_LLJ,
##           d = true_parameters$d, g = true_parameters$g, data_GBF = NULL,
##           data_LLJ = NULL, sample = TRUE, model = "2ht")
## data_GBF <- runif(P*S) < tmp$GBF
## ##dimnames(data_GBF) <- dimnames(data_LLJ)
## data_LLJ <- runif(P*S) < tmp$LLJ
## ##dimnames(data_LLJ) <- dimnames(data_GBF)
## tmp <- array(runif(P*S)>.5,dim=dim(data_GBF))
## data_GBF[tmp] <- data_LLJ[tmp] <- NA



## Arrays to store the iterates.
f_GBF <- f_LLJ <- array(dim=c(P,iter))
d <- g <- array(dim=c(S,iter))

group_level_parameters <- array(dim=c(2,4,iter),dimnames=list(c("loga","logb"),c("fGBF","fLLJ","d","g"),NULL))

## Function to estimate start points without edge problems.
mean_start <- function(x) {
    tmp <- mean(x,na.rm=T)
    if (!is.finite(tmp)) tmp <- .5
    if (tmp<.02) tmp<- .02
    if (tmp>.98) tmp<- .98
    return(tmp)
}

## Start points taken from marginal means.
f_GBF[,1] <- apply(data_GBF,1,mean_start)
f_LLJ[,1] <- apply(data_LLJ,1,mean_start)
d[,1] <- g[,1] <- apply(rbind(data_GBF,data_LLJ),2,mean_start)

## Uninformed staart points for group level parameters. 
group_level_parameters[,,1] <- log(c(1,2)) 

## Initial likelihoods
probs <- ll(f_GBF=f_GBF[,1],f_LLJ=f_LLJ[,1],d=d[,1],g=g[,1],data_GBF=data_GBF,data_LLJ=data_LLJ,sample=FALSE,model=model)

for (i in 2:iter) {
    cg <- exp(group_level_parameters[,,i-1]) ## Current group level parameters, back on raw sale. Saves typing.
    ## SAMPLE FOR F.
    ## Propose new f vector and calculate new log_probs from it.
    prop_f_GBF <- f_GBF[,i-1] + rnorm(P,mean=0,sd=0.02) ## Proposal added to last iteration.
    prop_f_LLJ <- f_LLJ[,i-1] + rnorm(P,mean=0,sd=0.02) ## Proposal added to last iteration.
##    tmp <- runif(P)<0.01  ## Add mixture of proposals from prior.
##    prop_f_GBF[tmp] <- runif(sum(tmp))
    prop_f_GBF[prop_f_GBF<=0] <- f_GBF[prop_f_GBF<=0,i-1] ## Any illegal proposals rejected.
    prop_f_GBF[prop_f_GBF>=1] <- f_GBF[prop_f_GBF>=1,i-1] ## Any illegal proposals rejected.
##    tmp <- runif(P)<0.01  ## Add mixture of proposals from prior.
##    prop_f_LLJ[tmp] <- runif(sum(tmp))
    prop_f_LLJ[prop_f_LLJ<=0] <- f_LLJ[prop_f_LLJ<=0,i-1] ## Any illegal proposals rejected.
    prop_f_LLJ[prop_f_LLJ>=1] <- f_LLJ[prop_f_LLJ>=1,i-1] ## Any illegal proposals rejected.
    prop_log_probs <- ll(f_GBF=prop_f_GBF,f_LLJ=prop_f_LLJ,d=d[,i-1],g=g[,i-1],data_GBF=data_GBF,data_LLJ=data_LLJ,sample=FALSE,model=model)

    ## Proposal for f change likelihoods only for one audio file. So we can do
    ## acceptance on those independently.
    tmp <- prop_log_probs$GBF$p - probs$GBF$p ## Difference in log likelihoods
    tmp <- tmp +
        dbeta(prop_f_GBF,shape1=cg["loga","fGBF"],shape2=cg["logb","fGBF"],log=TRUE) -
        dbeta(f_GBF[,i-1],shape1=cg["loga","fGBF"],shape2=cg["logb","fGBF"],log=TRUE) 
    tmp <- runif(P) < exp(tmp) ## M-H accept
    ## Build next vector Gibbs wise.
    f_GBF[,i] <- f_GBF[,i-1]
    f_GBF[tmp,i] <- prop_f_GBF[tmp]
    
    tmp <- prop_log_probs$LLJ$p - probs$LLJ$p ## Difference in log likelihoods
    tmp <- tmp +
        dbeta(prop_f_LLJ,shape1=cg["loga","fLLJ"],shape2=cg["logb","fLLJ"],log=TRUE) -
        dbeta(f_LLJ[,i-1],shape1=cg["loga","fLLJ"],shape2=cg["logb","fLLJ"],log=TRUE)
    tmp <- runif(P) < exp(tmp) ## M-H accept
    ## Build next vector Gibbs wise.
    f_LLJ[,i] <- f_LLJ[,i-1]
    f_LLJ[tmp,i] <- prop_f_LLJ[tmp]
    
    
    ## SAMPLE FOR D & G.
    ## Propose new d & g vectors and calculate new log_probs from them.
    prop_d <- d[,i-1] + rnorm(S,mean=0,sd=0.02) # Proposals based on last iterate.
##    tmp <- runif(S)<0.01  ## Add mixture of proposals from prior.
##    prop_d[tmp] <- runif(sum(tmp))
    prop_d[prop_d<=0] <- d[prop_d<=0,i-1] ## Any illegal proposals rejected.
    prop_d[prop_d>=1] <- d[prop_d>=1,i-1] ## Any illegal proposals rejected.
    prop_g <- g[,i-1] + rnorm(S,mean=0,sd=0.02) # Proposals based on last iterate.
##    tmp <- runif(S)<0.01  ## Add mixture of proposals from prior.
##    prop_g[tmp] <- runif(sum(tmp))
    prop_g[prop_g<=0] <- g[prop_g<=0,i-1] ## Any illegal proposals rejected.
    prop_g[prop_g>=1] <- g[prop_g>=1,i-1] ## Any illegal proposals rejected.
    prop_log_probs <- ll(f_GBF=f_GBF[,i],f_LLJ=f_LLJ[,i],d=prop_d,g=prop_g,data_GBF=data_GBF,data_LLJ=data_LLJ,sample=FALSE,model=model) ## Uses updated f. 

    ## Proposals for d & g change likelihoods only for one subject. So we can do
    ## accpetance on those independently.
    tmp <- prop_log_probs$GBF$s - probs$GBF$s + prop_log_probs$LLJ$s - probs$LLJ$s ## Difference in log likelihoods
    tmp <- tmp +
        dbeta(prop_d,shape1=cg["loga","d"],shape2=cg["logb","d"],log=TRUE) -
        dbeta(d[,i-1],shape1=cg["loga","d"],shape2=cg["logb","d"],log=TRUE)
    tmp <- tmp +
        dbeta(prop_g,shape1=cg["loga","g"],shape2=cg["logb","g"],log=TRUE) -
        dbeta(g[,i-1],shape1=cg["loga","g"],shape2=cg["logb","g"],log=TRUE)
    tmp <- runif(S) < exp(tmp) ## M-H accept

    ## Build next vector Gibbs wise.
    d[,i] <- d[,i-1]
    d[tmp,i] <- prop_d[tmp]
    g[,i] <- g[,i-1]
    g[tmp,i] <- prop_g[tmp]


    ## Proposals for group-level parameters change no data-level likelihoods.
    group_level_parameters[,,i] <- group_level_parameters[,,i-1] ## Copy old ones.
    prop_group <- group_level_parameters[,,i-1]+rnorm(8,mean=0,sd=0.1) ## New proposals
    prop_cg <- exp(prop_group) ## On raw scale. 

    for (j in 1:4) { ## Loop over the parameters, f_GBH, f_LLJ, d, g.
        tmp1 <- switch(j,f_GBF[,i],f_LLJ[,i],d[,i],g[,i])
        tmp <- sum(dbeta(tmp1,shape1=prop_cg["loga",j],shape2=prop_cg["logb",j],log=TRUE) -
                   dbeta(tmp1,shape1=cg["loga",j],shape2=cg["logb",j],log=TRUE)) +
                   sum(dexp(prop_group,rate=1/10,log=TRUE)) - sum(dexp(group_level_parameters[,,i-1],rate=1/10,log=TRUE)) ## Prior.
                  ##sum(dnorm(prop_group,mean=0,sd=3,log=TRUE)) - sum(dnorm(group_level_parameters[,,i-1],mean=0,sd=3,log=TRUE)) ## Prior.
        if (runif(1)<exp(tmp)) group_level_parameters[,j,i]<-prop_group[,j] ## Keep.
    }
 
    ## Progress report.
    if ((i%% 100)==0) cat(" ",i)
}

save.image(file=paste0("samples-",model,".RData"))

