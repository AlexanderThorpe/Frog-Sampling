rm(list=ls())

iter=5e3 ## MCMC length
load("data/processed_data_community.RData")
## Set this to work on the tree or ground frogs:
frogs <- "tree"
##frogs <- "ground"
if (frogs=="tree") {
    data<-data_tree[,,-4]
    expert<-expert_tree[,-4]
} else {
    data<-data_ground[,,-4]
    expert<-expert_ground[,-4]
}

## For now, I am just dropping the "no frog" data. This is because when
## the response is none, it is almost universal that they also avoided
## ticking any of the other frogs. But when the repsonse not-NGF, there
## is some confusion, with quite often other frogs identified (at least
## in the tree frog data, not so much the ground frog data). This is
## also supported by my censoring in the read-data where I removed some
## participants who made too many inconsistent uses of "no frog".

n_species <- dim(data)[3] ## Each person coded for how many frog species?
S <- dim(data)[2] ## How many citizen scientists?
P <- dim(data)[1] ## How many audio clips?

ll <- function(f,d,g,data,sample=FALSE) {
    ## Calculate the probabilites of TRUE responses and log-likelihoods.
    ## Random effects d and g have length dim(data)[2]. Random effects
    ## f is a matrix dim(data)[1],dim(data)[3].
    n_species <- dim(f)[2] ## Each person coded for how many frog species?

    ## MPT. Tree starts with frog (yes/no, f), then detect (yes/no, d),
    ## then guess (yes/no, g) if detect fails.
    all_d <- t(array(d,dim=c(length(d),nrow(f)))) ## Row-major prob d.
    all_g <- t(array(g,dim=c(length(d),nrow(f)))) ## Row-major prob g.
    all_probs <- array(NA,dim=c(nrow(f),length(d),n_species)) ## To store outputs.
    for (species in 1:n_species) {
        all_f <- array(f[,species],dim=c(nrow(f),length(d)))    ## Column-major prob of target species calling.
        all_probs[,,species] <- all_f*all_d + all_g*(1-all_d) ## Probability of "yes" coding.
    }
    if (sample) {
        return(all_probs) ## used in sampling synthetic data.
    } else {
        logprobs <- log(data*all_probs + (!data)*(1-all_probs)) ## log likelihoods
        return(list(
            p=apply(logprobs,c(1,3),sum,na.rm=TRUE), ## Sums for species x clip, used to update f.
            s=apply(logprobs,2,sum,na.rm=TRUE))) ## Sums for person, used to update d&g.
    }
}


## Arrays to store the iterates.
f <- array(dim=c(P,n_species,iter))
d <- g <- array(dim=c(S,iter))
group_level_parameters <- array(dim=c(2,2+n_species,iter),dimnames=list(c("loga","logb"),c(paste0("f_",dimnames(data)[[3]]),"d","g"),NULL))

## Function to estimate start points without edge problems.
mean_start <- function(x) {
    tmp <- mean(x,na.rm=T)
    if (!is.finite(tmp)) tmp <- .5
    if (tmp<.02) tmp<- .02
    if (tmp>.98) tmp<- .98
    return(tmp)
}

## Start points taken from marginal means.
f[,,1] <- apply(data,c(1,3),mean_start)
##d[,1] <- g[,1] <- apply(rbind(data_GBF,data_LLJ),2,mean_start)
d[,1] <- 0.5
g[,1] <- 0.5

## Uninformed start points for group level parameters. 
group_level_parameters[,,1] <- log(c(1,2)) 

## If re-starting from end of last run, do this:
##f[,,1]<-f[,,iter]
##d[,1]<-d[,iter]
##g[,1]<-g[,iter]
##group_level_parameters[,,1]<-group_level_parameters[,,iter]
## Then run from here down. 


## Initial likelihoods
probs <- ll(f=f[,,1],d=d[,1],g=g[,1],data=data,sample=FALSE)

for (i in 2:iter) {
    cg <- exp(group_level_parameters[,,i-1]) ## Current group level parameters, back on raw sale. Saves typing.
    ## SAMPLE FOR F.
    ## Propose new f vector and calculate new log_probs from it.
    prop_f <- f[,,i-1] + rnorm(P*n_species,mean=0,sd=0.02) ## Proposal added to last iteration.
    if (i<(iter/2)) {
        tmp <- runif(P*n_species)<0.01  ## Add 1% mixture of proposals from prior.
        prop_f[tmp] <- runif(sum(tmp))
    }
    prop_f[prop_f<=0] <- f[,,i-1][prop_f<=0] ## Any illegal proposals rejected.
    prop_f[prop_f>=1] <- f[,,i-1][prop_f>=1] ## Any illegal proposals rejected.
    prop_log_probs <- ll(f=prop_f,d=d[,i-1],g=g[,i-1],data=data,sample=FALSE)

    ## Proposal for f change likelihoods only for one audio file, for one species.
    ##. So we can do acceptance on those independently.
    tmp <- prop_log_probs$p - probs$p ## Difference in log likelihoods
    ## To get the difference in group dist likelihood, need to transpose so that
    ## dim(tmp) matches cg.
    tmp_grp <-  dbeta(t(prop_f),shape1=cg["loga",1:n_species],shape2=cg["logb",1:n_species],log=TRUE) -
        dbeta(t(f[,,i-1]),shape1=cg["loga",1:n_species],shape2=cg["logb",1:n_species],log=TRUE)
    tmp <- tmp + t(tmp_grp) ## Back transpose.
    keep <- runif(P*n_species) < exp(tmp) ## M-H accept
    ## Build next vector Gibbs wise.
    f[,,i] <- f[,,i-1]
    f[,,i][keep] <- prop_f[keep]
    
    ## SAMPLE FOR D & G.
    ## Propose new d & g vectors and calculate new log_probs from them.
    prop_d <- d[,i-1] + rnorm(S,mean=0,sd=0.02) # Proposals based on last iterate.
    if (i<(iter/2)) {
        tmp <- runif(S)<0.01  ## Add mixture of proposals from prior.
        prop_d[tmp] <- runif(sum(tmp))
    }
    prop_d[prop_d<=0] <- d[prop_d<=0,i-1] ## Any illegal proposals rejected.
    prop_d[prop_d>=1] <- d[prop_d>=1,i-1] ## Any illegal proposals rejected.
    prop_g <- g[,i-1] + rnorm(S,mean=0,sd=0.02) # Proposals based on last iterate.
    if (i<(iter/2)) {
        tmp <- runif(S)<0.01  ## Add mixture of proposals from prior.
        prop_g[tmp] <- runif(sum(tmp))
    }
    prop_g[prop_g<=0] <- g[prop_g<=0,i-1] ## Any illegal proposals rejected.
    prop_g[prop_g>=1] <- g[prop_g>=1,i-1] ## Any illegal proposals rejected.
    prop_log_probs <- ll(f=f[,,i],d=prop_d,g=prop_g,data=data,sample=FALSE) ## Uses updated f. 

    ## Proposals for d & g change likelihoods only for one subject. So we can do
    ## accpetance on those independently.
    tmp <- prop_log_probs$s - probs$s ## Difference in log likelihoods
    tmp <- tmp + ## Prior likelihood for "d".
        dbeta(prop_d,shape1=cg["loga","d"],shape2=cg["logb","d"],log=TRUE) -
        dbeta(d[,i-1],shape1=cg["loga","d"],shape2=cg["logb","d"],log=TRUE)
    tmp <- tmp + ## Prior likelihood for "g".
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
    prop_group <- group_level_parameters[,,i-1]+rnorm(2*(n_species+2),mean=0,sd=0.1) ## New proposals
    prop_cg <- exp(prop_group) ## On raw scale. 

    for (j in 1:(n_species+2)) { ## Loop over the parameters, f, d, g.
        if (j==(n_species+1)) {
            tmp1 <- d[,i]
        } else if (j==(n_species+2)) {
            tmp1 <- g[,i]
        } else {
            tmp1 <- f[,j,i] ## One species at a time.
        }
        tmp <- sum(dbeta(tmp1,shape1=prop_cg["loga",j],shape2=prop_cg["logb",j],log=TRUE) -
                   dbeta(tmp1,shape1=cg["loga",j],shape2=cg["logb",j],log=TRUE)) +
                  sum(dnorm(prop_group[,j],mean=0,sd=3,log=TRUE)) - sum(dnorm(group_level_parameters[,j,i-1],mean=0,sd=3,log=TRUE)) ## Prior.
        if (runif(1)<exp(tmp)) group_level_parameters[,j,i]<-prop_group[,j] ## Keep.
    }
 
    ## Progress report.
    if ((i%% 100)==0) cat(" ",i)
}

save.image(file=paste0(frogs,".RData"))

