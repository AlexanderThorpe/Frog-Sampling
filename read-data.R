rm(list=ls())
require(tidyverse)
data <- read_csv(file="Ingar matrix_Expert.csv")


## Split into frog species.
data_GBF <- data %>% select(-`Grand Total`)  %>% filter(`Possible answer`=="GBF")
for (i in 3:ncol(data_GBF)) data_GBF[,i] <- as.logical(data_GBF[[i]])
tmp <- data_GBF[[1]] ## Save these for row names later.
data_GBF <- as.matrix(data_GBF[-(1:2)])
row.names(data_GBF) <- tmp

data_LLJ <- data %>% select(-`Grand Total`)  %>% filter(`Possible answer`=="LLJ")
for (i in 3:ncol(data_LLJ)) data_LLJ[,i] <- as.logical(data_LLJ[[i]])
tmp <- data_LLJ[[1]] ## Save these for row names later.
data_LLJ <- as.matrix(data_LLJ[-(1:2)])
row.names(data_LLJ) <- tmp

## Check how many classificuations per person.
n_class <- apply(data_GBF,2,function(x) sum(!is.na(x))) ## Same for each species.
keep <- names(n_class[n_class>2])
keep <- keep[-grep(x=keep,pattern="kurstinburleson86")] ## He's not good.


## Remove people who made 5 or fewer classifications.
data_GBF <- data_GBF[,keep]
data_LLJ <- data_LLJ[,keep]

## Split out the expert ratings from Oliver.
expert_GBF <- data_GBF[,"Expert"]
expert_LLJ <- data_LLJ[,"Expert"]
data_GBF <- data_GBF[,-which(colnames(data_GBF)=="Expert")]
data_LLJ <- data_LLJ[,-which(colnames(data_LLJ)=="Expert")]


save(file="processed_data.RData",list=c("data_GBF","data_LLJ","expert_GBF","expert_LLJ"))
