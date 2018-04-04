#################################################
## Description: Create reduced SNP dataset from 
## enormous dataset 
## Author: Levon Demirdjian
## Email: levondem@gmail.com
## Last updated: 4/28/2017
#################################################

## SNP data processing
library(data.table) ## For using fread

## Read in list of AIMS
AIMs <- as.character(read.csv("Data/AIMs.csv", header = FALSE)[,1])

## Read in list of SNPs we're interested in
SNP_list <- sort(as.character(read.table('Data/SNP_List.txt', header = TRUE, sep = "\t")[,1]))
n        <- length(SNP_list)

## Filter out the AIMs
SNP_list <- setdiff(SNP_list, AIMs)
n        <- length(SNP_list)

## Read in list of SNPs in high LD to the SNPs we requested
LD_SNPs  <- read.csv('Data/SNP_Proxy.csv', header = TRUE)

## Read in list of SNPs in SNP dataset
SNP_data   <- fread('../SNP_data/CNP_TEXT.tped', colClasses = c("NULL", "character", rep("NULL", 2888)))

## If a SNP in our list failed QC, use a SNP in high LD instead
SNPs <- rep(0, n)
for(i in 1:n){
  
  ## If SNP is available, use it. Else find a SNP in high LD
  if(SNP_list[i] %in% SNP_data$V2){
    SNPs[i] <- as.character(SNP_list[i])
  }else{
    
    ## There may be multiple SNPs in high LD with the SNP we're interested in
    sub_SNPs <- LD_SNPs[which(LD_SNPs[,1] == SNP_list[i]), 2]
    
    ## If there is more than one such SNP, pick the first one
    if(length(which(sub_SNPs %in% SNP_data$V2)) > 0){
      SNPs[i] <- as.character(sub_SNPs[which(sub_SNPs %in% SNP_data$V2)[1]])
    }else{
      
      ## If SNP and proxy are not available, mark as missing
      SNPs[i] <- NA
    }
  }
}

## Remove missing SNPs and possibly double-counted SNPs
SNPs <- SNPs[!is.na(SNPs)]
SNPs <- unique(SNPs)

## Save the big dataset's indexes of these SNPs
valid_SNPs <- which(SNP_data$V2 %in% SNPs)

## Create matrix to store genotype data into. This may take a while
gene_data  <- matrix(0, nrow = length(valid_SNPs), ncol = 2890)
for(i in 1:length(valid_SNPs)){
  cat("We're on line", i, "\n")
  cur_line <- system(paste("sed '", valid_SNPs[i], "q;d' ../SNP_data/CNP_TEXT.tped", sep = ''), intern = TRUE)
  gene_data[i, ] <- strsplit(cur_line, " ")[[1]]
}

## Eliminate unncessary columns in gene_data
gene_data <- gene_data[,c(-1,-3,-4)]

## Combine each subjects 2 alleles into 1 column
clean_data     <- matrix(0, nrow = nrow(gene_data), ncol = 1444)
clean_data[,1] <- gene_data[,1]
for(j in 1:nrow(gene_data)){
  for(i in 2:((ncol(gene_data)+1)/2)){
    clean_data[j,i] <- paste(gene_data[j,2*(i-1)], ':', gene_data[j,2*(i-1)+1], sep = '')
  }
}

## Create column names (e.g. subject names)
subj_names <- read.table('../SNP_data/CNP_TEXT.tfam')[,1]
colnames(clean_data) <- c("SNP", subj_names)

# Write output
#write.table(clean_data, file = "SNP_Data_new.txt", quote = FALSE,
#            row.names = FALSE)

## Look at function of selected SNPs
# SNP_info <- read.table('SNP_List.txt', header = TRUE, sep = "\t")
# SNP_info <- SNP_info[,c(1,4,5)]
# SNP_info[which(SNP_info[,1] %in% clean_data[,1]), 3]
