##################################################### 
## AIMs for CNP
## Author: Levon Demirdjian
## Last Updated: 04/03/2018 
## Description: Look at results of clustering 
## individuals based on genetic ancestry
##################################################### 

## Load in data
my.data <- read.table("k5.outfile")

## Delete reference data
my.data <- my.data[1044:2486,]
n       <- nrow(my.data)

## Assign each CNP subject their argmax ancestry
pop <- sapply(1:n, function(x) as.numeric(which.max(my.data[x, 2:6])))
hist(pop)