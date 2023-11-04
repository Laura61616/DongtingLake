library(piecewiseSEM)
library(nlme)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggpmisc)
library(dplyr)
source('00functions.R')
#----------------------------------Beta---------------------------------#
##########################################################################
### SEM
##########################################################################
data <- readRDS('Data/SEMdata_beta.RDS')$StanNAfree
data[,5:13] <- as.data.frame(scale(data[,5:13]))
data$Site <- paste(data$s1, data$s2)
head(data)
# ---- lme for single path ####
data <- na.omit(data)
summary(lm(TN ~ YRQ + HLU + MAT, data = data))
summary(lm(TP ~ YRQ + HLU + MAT, data = data))
summary(lm(Bloom ~ HLU + MAT + YRQ + TN + TP, data = data))
summary(lm(taxaDiss ~ HLU + MAT+ YRQ + TN + TP + Bloom, data = data))
summary(lm(funcDiss ~ HLU + MAT+ YRQ + TN + TP + Bloom + taxaDiss, data = data))

# ---- full mod ####
beta_pSEM_randomList <- list(
  lm(TN ~ YRQ + HLU + MAT, data = data),
  lm(TP ~ YRQ + HLU + MAT, data = data),
  lm(Bloom ~ HLU + MAT + YRQ + TN + TP, data = data),
  lm(taxaDiss ~ HLU + MAT + YRQ + TN + TP + Bloom, data = data),
  lm(funcDiss ~ HLU + MAT + YRQ + TN + TP + Bloom + taxaDiss, data = data) )
# Run goodness-of-fit tests
beta.psem <- as.psem(beta_pSEM_randomList)

# Evaluate path significance using unstandardized coefficients
coefs(beta_pSEM_randomList, kelp, standardize = "none")

# Obtain standardized regression coefficients
coefs(beta_pSEM_randomList, kelp, standardize = "scale")

# Explore individual model fits
beta.summ <- summary(beta.psem)
beta.summ
# ---- best mod ####
data <- na.omit(data)
beta_pSEM_randomList <- list(
  lm(TN ~ YRQ + HLU + MAT, data = data),
  lm(TP ~ YRQ + HLU + MAT, data = data),
  lm(Bloom ~ HLU + MAT + YRQ + TN + TP, data = data),
  lm(taxaDiss ~ HLU + MAT + YRQ + TN + TP + Bloom, data = data),
  lm(funcDiss ~ HLU + MAT + YRQ + TN + TP + Bloom + taxaDiss, data = data) )
# Run goodness-of-fit tests
beta.psem <- as.psem(beta_pSEM_randomList)
summary(beta.psem)$Cstat
# Obtain standardized regression coefficients
sink('Data/pSEM/BestMod_stanCoeff_beta.txt')
coefs(beta_pSEM_randomList, kelp, standardize = "scale")
sink()
# Explore individual model fits
beta.summ <- summary(beta.psem)
sink('Data/pSEM/BestMod_summary_beta.txt')
beta.summ
sink()
##########################################################################
### selected LME fit plot
##########################################################################
LMPairs <- matrix(c('HLU','TP',
                    'HLU','taxaDiss',
                    'TP','taxaDiss',
                    'taxaDiss','funcDiss'), ncol = 4)
myColor <- c('#FF4F4F','#41AE76','orange','skyblue')

rName <- c('% Human land-use\nintensity', "Total phosphorus",
           'Taxonomic\ndissimilarity','Functional\ndissimilarity')
names(rName) <- c('HLU','TP','taxaDiss','funcDiss')

gg <- NULL
for (i in 1:ncol(LMPairs)) {
  xypair <- LMPairs[,i]
  cat(paste0('Processing ', xypair[1],' verus ', xypair[2],'...\n'))
  
  df <- data.frame(x = data[,xypair[1]], 
                   y = data[,xypair[2]], 
                   Site = data$Site)
  
  fit <- summary(lm(y ~ x, data = na.omit(df)))
  p <- fit$coefficients['x','Pr(>|t|)']
  r <- fit$adj.r.squared
  ind <- case_when(p < 0.001 ~ '***',
                   (p >= 0.001)&(p < 0.01) ~ '**', 
                   (p >= 0.01)&(p < 0.05) ~ '*', 
                   (p >= 0.05)&(p < 0.1) ~ '.', 
                   p >= 0.1 ~ '')
  lab <- paste0('R2 = ', round(r,2), ' p = ', round(p,3), ' ',ind)
  if (p < 0.001) {
    lab <- paste0('R2 = ', round(r,2), ' p < 0.001 ***')
  }
  g <- ggplot(df, aes(x = x, y = y))+
    geom_smooth(method = 'glm', formula = y ~ x, level = 0.95,
                lwd = 1, alpha = 0.1, fill = myColor[i], color = myColor[i])+
    geom_point(alpha = 0.7, size = 1, fill = myColor[i], color = myColor[i])+
    ylab(rName[xypair[2]])+xlab(rName[xypair[1]])+
    annotate("text", label = lab, x = 0, y = 0.8, size = 4) +
    # geom_label(x = 0.1, y = 0.8, label = lab)+
    theme_bw()+
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 45,hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 3) ,
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = 'grey'))
  g
  gg$value <- g; names(gg)[length(gg)] <- i
}
names(gg)
p1 <- ggarrange(gg[[1]],gg[[2]],gg[[3]],gg[[4]], align = 'hv')
p1

pdf('image/02_LME_beta.pdf', width = 6, height = 5)
p1
dev.off()