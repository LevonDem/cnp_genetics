#############################################
## Description: Create a genetic risk score 
## and explore relations to phenotypes and
## ancestry information
## Author: Levon Demirdjian
## Email: levondem@gmail.com
## Last updated: 04/03/2018
#############################################

library(genetics)    ## For genetic analysis
library(ggplot2)     ## For plotting
library(gridExtra)   ## For plotting
library(grpreg)      ## For penalized regression

## Plotting function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## Load in SNP list and AIMs list
SNP_data       <- as.matrix(read.table('Data/SNP_Data_new.txt', row.names = 1, header = T))
subj_names     <- colnames(SNP_data)
SNP_names      <- rownames(SNP_data)
p              <- length(SNP_names)
ancestry       <- read.table('Data/k5.outfile')[,2:5]
AIM_subj_names <- as.character(read.table('Data/clustering_subject_index.txt')[,2])
AIM_SNPs       <- read.table('Data/AIM_SNPs.txt', header = T)
AIM_SNP_names  <- substr(colnames(AIM_SNPs), 2, 6)
  
######################################## 
## Load in chapman scores (phenotype) ##
######################################## 

## Chapman data
chapman_data   <- read.csv('Data/anhedonia.csv')
ChapPhy        <- chapman_data$CHAPPHY_TOTAL
ChapSoc        <- chapman_data$CHAPSOC_TOTAL
ch_names       <- paste('sub-', chapman_data$PTID, sep = '')
SNP_data_brain <- SNP_data[,which(substr(subj_names, 2, 6) %in% substr(ch_names, 5, 9))]

## Only retain subjects for which we have chapman scores
valid   <- which(substr(ch_names, 5, 9) %in% substr(subj_names, 2, 6))
ChapPhy <- ChapPhy[valid]
ChapSoc <- ChapSoc[valid]

subject_names <- substr(ch_names, 5, 9)[valid]
n             <- length(subject_names)
ancestry      <- ancestry[which(AIM_subj_names %in% subject_names), ]
AIM_SNPs      <- as.matrix(AIM_SNPs[,which(AIM_SNP_names %in% subject_names)])
num_AIMs      <- nrow(AIM_SNPs)

## Separate subjects by illness
Y <- rep(0, n)
Y[substr(subject_names, 1, 1) == 1] <- 0 ## Healthy Controls
Y[substr(subject_names, 1, 1) == 5] <- 1 ## Schizophrenia
Y[substr(subject_names, 1, 1) == 6] <- 2 ## Bipolar Disorder
Y[substr(subject_names, 1, 1) == 7] <- 3 ## ADHD

## Make Y into a factor
Y[Y == 0] <- "control"
Y[Y == 1] <- "schz"
Y[Y == 2] <- "bplr"
Y[Y == 3] <- "adhd"
Y <- factor(Y)
Y <- relevel(Y, ref = "control")
Y <- factor(Y, levels(Y)[c(1,3,4,2)])


##############################
## Impute missing genotypes ##
##############################

SNPs_imputed <- t(SNP_data_brain)
num_missing  <- rep(0, p)
for(i in 1:p){
  
  cur_SNP     <- factor(SNP_data_brain[i,]) 
  SNP_levels  <- sort(levels(cur_SNP))
  
  ## Impute missing SNPs from empirical densities if necessary
  if("0:0" %in% SNP_levels){
    num_missing[i] <- sum(cur_SNP == "0:0")
    emp_dist <- summary(cur_SNP)[2:4]/sum(summary(cur_SNP)[2:4])
    imp_SNPs <- sample(SNP_levels[2:4], sum(cur_SNP == "0:0"), prob = emp_dist, replace = TRUE)
    cur_SNP[cur_SNP == "0:0"] <- imp_SNPs
    cur_SNP          <- factor(cur_SNP, levels(cur_SNP)[2:4])
    SNPs_imputed[,i] <- as.character(cur_SNP)
  }

}

AIMs_imputed <- t(AIM_SNPs)
num_missing  <- rep(0, num_AIMs)
for(i in 1:num_AIMs){
  
  cur_SNP     <- factor(AIM_SNPs[i,]) 
  SNP_levels  <- sort(levels(cur_SNP))
  
  ## Impute missing SNPs from empirical densities if necessary
  if("0:0" %in% SNP_levels){
    num_missing[i] <- sum(cur_SNP == "0:0")
    emp_dist <- summary(cur_SNP)[2:4]/sum(summary(cur_SNP)[2:4])
    imp_SNPs <- sample(SNP_levels[2:4], sum(cur_SNP == "0:0"), prob = emp_dist, replace = TRUE)
    cur_SNP[cur_SNP == "0:0"] <- imp_SNPs
    cur_SNP          <- factor(cur_SNP, levels(cur_SNP)[2:4])
    AIMs_imputed[,i] <- as.character(cur_SNP)
  }
  
}


############################################
## Check LD between AIMs and non-AIM SNPs ##
############################################

# LD_matrix <- matrix(0, nrow = p, ncol = num_AIMs)
# for(i in 1:p){
#   cat("On SNP ", i, ".\n", sep = '')
#   SNP1 <- genotype(factor(SNPs_imputed[,i]), sep = ":")
#   for(j in 1:num_AIMs){
#     AIM1 <- genotype(factor(AIMs_imputed[,j]), sep = ":")
#     LD_results     <- LD(SNP1, AIM1)
#     LD_matrix[i,j] <- LD_results$`P-value`
#   }
# }


########################################################
## Compute additive, recessive, and dominant measures ##
########################################################

## Note: the GRS assumes an additive genetic model
SNPs_additive  <- matrix(0, nrow = nrow(SNPs_imputed), ncol = p)
SNPs_dominant  <- matrix(0, nrow = nrow(SNPs_imputed), ncol = p)
SNPs_recessive <- matrix(0, nrow = nrow(SNPs_imputed), ncol = p)
SNP_weights    <- rep(0, p)
SNP_pvals      <- rep(0, p)

## This is our response phenotype
Chapman <- ChapSoc ## Can change this to ChapSoc or ChapPhy
for(i in 1:p){
  
  cur_SNP     <- factor(SNPs_imputed[,i]) 
  SNP_levels  <- sort(levels(cur_SNP))
  
  minor_allele <- strsplit(SNP_levels[1], ":")[[1]][1]
  major_allele <- strsplit(SNP_levels[3], ":")[[1]][1]
  
  ## See which allele is the risk allele
  homo1 <- paste(minor_allele, ":", minor_allele, sep = '')
  homo2 <- paste(major_allele, ":", major_allele, sep = '')
  
  ## Which group has a higher Chapman score
  score1 <- mean(ChapPhy[cur_SNP == homo1])
  score2 <- mean(ChapPhy[cur_SNP == homo2])
  
  risk_allele <- ""
  safe_allele <- ""
  if(score1 >= score2){
    risk_allele <- minor_allele
    safe_allele <- major_allele
  }else{
    risk_allele <- major_allele
    safe_allele <- minor_allele
  }

  g1 <- paste(safe_allele, ":", safe_allele, sep = '')
  g2 <- paste(safe_allele, ":", risk_allele, sep = '')
  g3 <- paste(risk_allele, ":", safe_allele, sep = '')
  g4 <- paste(risk_allele, ":", risk_allele, sep = '')
  
  ## Save weight of risk allele (use tvalues as weights?)
  ## Here we are assuming an additive model
  risk_vector    <- rep(0, length(cur_SNP))
  risk_vector[cur_SNP == g1] <- 0
  risk_vector[cur_SNP == g2 | cur_SNP == g3] <- 1  
  risk_vector[cur_SNP == g4] <- 2  
    
  results        <- summary(lm(Chapman ~ risk_vector))

  if(results$coefficients[2,1] < 0){
    g1 <- paste(risk_allele, ":", risk_allele, sep = '')
    g4 <- paste(safe_allele, ":", safe_allele, sep = '')
    risk_vector[cur_SNP == g1] <- 0
    risk_vector[cur_SNP == g4] <- 2 
    results        <- summary(lm(Chapman ~ risk_vector))
  }
  
  SNP_weights[i] <- results$coefficients[2,1]
  SNP_pvals[i]   <- results$coefficients[2,4]
  
  ## Additive model
  SNPs_additive[cur_SNP == g1, i] <- 0
  SNPs_additive[cur_SNP == g2, i] <- 1
  SNPs_additive[cur_SNP == g3, i] <- 1
  SNPs_additive[cur_SNP == g4, i] <- 2
  
  ## Dominant model
  SNPs_dominant[cur_SNP == g1, i] <- 0
  SNPs_dominant[cur_SNP == g2, i] <- 1
  SNPs_dominant[cur_SNP == g3, i] <- 1
  SNPs_dominant[cur_SNP == g4, i] <- 1
  
  ## Recessive model
  SNPs_recessive[cur_SNP == g1, i] <- 0
  SNPs_recessive[cur_SNP == g2, i] <- 0
  SNPs_recessive[cur_SNP == g3, i] <- 0
  SNPs_recessive[cur_SNP == g4, i] <- 1
  
}

## Look at relationship between additive model and ancestry
ancestry_pvals <- rep(0, p)
for(i in 1:p)
  ancestry_pvals[i] <- summary(lm(ancestry[,4] ~ SNPs_additive[,i]))$coef[2,4]

white <- ancestry[,4]

## Plot boxplots of pvalues
library(ggplot2)
library(grid)
library(gridExtra)

par(mfrow = c(3,5))
par(mar = 2 * c(2.2,2,1,1))
bp <- list()
k  <- 1
for(i in 1:p){
  if(ancestry_pvals[i] <= 0.0001){
    boxplot(white ~ SNPs_additive[,i],
            ylab = "%_Pop_4_membership",
            xlab = "# of risk alleles",
            main = paste(SNP_names[i]))
    df <- data.frame(Risk_Alleles = factor(SNPs_additive[,i]), Pop4 = white)
    bp[[k]] <- ggplot(df, aes(x = Risk_Alleles, y = Pop4)) +
      #geom_boxplot() + 
      geom_point(alpha = 0.5, size = 0.5, position = position_jitter(width = 0.1))
      
    k <- k + 1
  }
}

grid.arrange(bp[[1]] + coord_flip(), 
             bp[[2]] + coord_flip(), 
             bp[[3]] + coord_flip(), 
             bp[[4]] + coord_flip(), 
             bp[[5]] + coord_flip(), 
             bp[[6]] + coord_flip(), 
             bp[[7]] + coord_flip(), 
             bp[[8]] + coord_flip(), 
             bp[[9]] + coord_flip(), 
             bp[[10]] + coord_flip(), 
             bp[[11]] + coord_flip(),
             bp[[12]] + coord_flip(), 
             bp[[13]] + coord_flip(), 
             bp[[14]] + coord_flip(), 
             widths = c(0.25, 0.25, 0.25, 0.25))

## Weighted GRS (more weight to stronger SNPs)
#SNP_weights[SNP_pvals >= 0.1] <- 0
GRS        <- sapply(1:n, function(i) sum(SNP_weights * SNPs_additive[i,]))
grs_model  <- lm(Chapman ~ Y + GRS * Y - GRS - 1)
grs_model2 <- lm(Chapman ~ Y + white)

summary(grs_model)
summary(grs_model2)

## Put all data into data frame
dat   <- data.frame(GRS, Chapman, Y)

summary(lm(GRS ~  ancestry[,4]))

plot(white, GRS, 
     xlab = "% Population 4 membership", 
     ylab = "Genetic risk score", 
     main = "Ancestry predicts GRS (p-val < 2e-16, R^2 = 17%)")
abline(lm(GRS ~ white), col = 'blue')

#########################
## Plot GRS vs Chapman ##
#########################

# Placeholder plot - prints nothing at all
empty <- ggplot() + geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

# Scatterplot of x and y variables
scatter <- ggplot(dat, aes(x = GRS, y = Chapman, color = Y)) + 
  geom_point(aes(color=Y)) + 
  geom_smooth(method = "lm", fill = NA) +
  scale_color_manual(values = c("steelblue", "purple", "forestgreen", "darkorange")) + 
  theme(legend.position=c(1,1),legend.justification=c(1,1)) 

# Marginal density of x - plot on top
plot_top <- ggplot(dat, aes(GRS, fill=Y)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c("steelblue", "purple", "forestgreen", "darkorange")) + 
  theme(legend.position = "none")

# Marginal density of y - plot on the right
plot_right <- ggplot(dat, aes(Chapman, fill=Y)) + 
  geom_density(alpha=.5) + 
  coord_flip() + 
  scale_fill_manual(values = c("steelblue", "purple", "forestgreen", "darkorange")) + 
  theme(legend.position = "none") 

# Arrange the plots together, with appropriate height and width for each row and column
grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

## Stargazer results
# library(stargazer)
# stargazer(grs_model, intercept.bottom = FALSE, #report = "vcsp*",
#           single.row = TRUE, no.space = TRUE,
#           title = "Genetic Risk Score vs Chapman (physical)")

