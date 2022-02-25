rm(list=ls())
require(tidyverse)
data <- read_csv(file="Ingar matrix 2_Expert.csv")


## Split into frog species.
data_GBF <- data %>% select(-`Grand Total`)  %>% filter(`Species`=="GBF")
for (i in 3:ncol(data_GBF)) data_GBF[,i] <- as.logical(data_GBF[[i]])
tmp <- data_GBF[[1]] ## Save these for row names later.
data_GBF <- as.matrix(data_GBF[-(1:2)])
row.names(data_GBF) <- tmp

data_LLJ <- data %>% select(-`Grand Total`)  %>% filter(`Species`=="LLJ")
for (i in 3:ncol(data_LLJ)) data_LLJ[,i] <- as.logical(data_LLJ[[i]])
tmp <- data_LLJ[[1]] ## Save these for row names later.
data_LLJ <- as.matrix(data_LLJ[-(1:2)])
row.names(data_LLJ) <- tmp

data_PTF <- data %>% select(-`Grand Total`)  %>% filter(`Species`=="PTF")
for (i in 3:ncol(data_PTF)) data_PTF[,i] <- as.logical(data_PTF[[i]])
tmp <- data_PTF[[1]] ## Save these for row names later.
data_PTF <- as.matrix(data_PTF[-(1:2)])
row.names(data_PTF) <- tmp

data_SMF <- data %>% select(-`Grand Total`)  %>% filter(`Species`=="SMF")
for (i in 3:ncol(data_SMF)) data_SMF[,i] <- as.logical(data_SMF[[i]])
tmp <- data_SMF[[1]] ## Save these for row names later.
data_SMF <- as.matrix(data_SMF[-(1:2)])
row.names(data_PTF) <- tmp

## Check how many classifications per person.
n_class_ground <- apply(data_GBF,2,function(x) sum(!is.na(x))) ## Same for each species.
n_class_tree <- apply(data_LLJ,2,function(x) sum(!is.na(x)))
keep_ground <- names(n_class_ground[n_class_ground>2])
keep_tree <- names(n_class_tree[n_class_tree>2])

keep <- unique(sort(c(keep_ground,keep_tree)))

## Remove people who made 5 or fewer classifications.
data_GBF <- data_GBF[,keep]
data_LLJ <- data_LLJ[,keep]
data_PTF <- data_PTF[,keep]
data_SMF <- data_SMF[,keep]

expert_GBF <- data_GBF[,"olli"]
expert_LLJ <- data_LLJ[,"olli"]
expert_PTF <- data_PTF[,"olli"]
expert_SMF <- data_SMF[,"olli"]

data_GBF <- data_GBF[,-which(colnames(data_GBF)=="olli")]
data_LLJ <- data_LLJ[,-which(colnames(data_LLJ)=="olli")]
data_PTF <- data_PTF[,-which(colnames(data_PTF)=="olli")]
data_SMF <- data_SMF[,-which(colnames(data_SMF)=="olli")]

data <- array(c(data_GBF,data_LLJ,data_PTF,data_SMF), dim=c(length(tmp),length(keep),4), dimnames=list(tmp,keep,list("GBF","LLJ","PTF","SMF")))
expert_data <- array(c(expert_GBF,expert_GBF,expert_PTF,expert_SMF), dim=c(length(tmp), 4), dimnames=list(tmp,list("GBF","LLJ","PTF","SMF")))

save(file="processed_data_4_frogs.RData",list=c("data","expert_data"))
