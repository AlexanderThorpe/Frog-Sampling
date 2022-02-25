rm(list=ls())
require(tidyverse)
data <- read_csv(file="Ingar Matrix 3 - complete community.csv")

## Split into frog species.
data_list <- expert_list <- list()
for (j in c("BMTF","CEF","EBF","GBF","GSF","LLJ","NGF","NTF","PTF","SMF","WTF")) {
    data_list[[j]] <- data %>% select(-`Grand Total`)  %>% filter(`Row Labels...6`==j)
    for (i in 8:ncol(data_list[[j]])) data_list[[j]][,i] <- as.logical(data_list[[j]][[i]])
    tmp <- data_list[[j]][[5]] ## Save these for row names later.
    data_list[[j]] <- as.matrix(data_list[[j]][-(1:7)])
    row.names(data_list[[j]]) <- tmp
    expert_list[[j]] <- data$Olli[data$`Row Labels...6`==j] ## Expert ratings.
}

data_tree <- array(unlist(data_list[c("BMTF","LLJ","PTF","NTF","GSF","WTF")]),
              dim=c(dim(data_list[[1]]),6), dimnames=c(dimnames(data_list[[1]]),list(c("BMTF","LLJ","PTF","NTF","GSF","WTF"))))
data_ground <- array(unlist(data_list[c("CEF","EBF","GBF","NGF","SMF")]),
                     dim=c(dim(data_list[[1]]),5), dimnames=c(dimnames(data_list[[1]]),list(c("CEF","EBF","GBF","NGF","SMF"))))

expert_tree <- array(unlist(expert_list[c("BMTF","LLJ","PTF","NTF","GSF","WTF")]),
                     dim=c(length(expert_list[[1]]),6),dimnames=dimnames(data_tree)[-2])
expert_ground <- array(unlist(expert_list[c("CEF","EBF","GBF","NGF","SMF")]),
                     dim=c(length(expert_list[[1]]),5),dimnames=dimnames(data_ground)[-2])


## Check how many classifications per person.
n_class_ground <- apply(data_ground,2:3,function(x) sum(!is.na(x)))[,1] ## Same for each species.
n_class_tree <- apply(data_tree,2:3,function(x) sum(!is.na(x)))[,1]
keep_ground <- names(n_class_ground[n_class_ground>2])
keep_tree <- names(n_class_tree[n_class_tree>2])

## Remove people who made 5 or fewer classifications.
data_ground <- data_ground[,keep_ground,]
data_tree <- data_tree[,keep_tree,]

## Does it look reasonable?
##> apply(data_tree,3,table)
##      BMTF  LLJ  PTF  NTF  GSF  WTF
##FALSE 4822 5547 4827 4963 5372 5240
##TRUE  1319  594 1314 1178  769  901
##> apply(data_ground,3,table)
##       CEF  EBF  GBF  NGF  SMF
##FALSE 3577 4210 5236 5422 3578
##TRUE  2554 1921  895  709 2553
##> apply(rbind(expert_tree,1,0),2,table)-1
##  BMTF  LLJ PTF  NTF  GSF  WTF
##0 1214 1260 481 1260 1251 1260
##1   46    0 779    0    9    0
##> apply(rbind(expert_ground,1,0),2,table)-1
##   CEF EBF  GBF  NGF  SMF
##0   92 370 1133 1260   75
##1 1168 890  127    0 1185
 

save(file="processed_data_community.RData",list=c("data_tree","data_ground","expert_tree","expert_ground"))

