source('00functions.R')
#----------------------------------Alpha---------------------------------#
##########################################################################
### SEM
##########################################################################
data <- readRDS('Data/SEMdata_alpha.RDS')
data[,1:9] <- scale(data[,1:9])
head(data)
# ---- lme for single path ####
data <- na.omit(data)
summary(lme(TN ~ YRQ + HLU + MAT, random = ~1|Site, data = data))
summary(lme(TP ~ YRQ + HLU + MAT, random = ~1|Site, data = data))
summary(lme(Bloom ~ HLU + MAT+ YRQ + TN + TP, random = ~1|Site, data = data))
summary(lme(TRic ~ HLU + MAT+ YRQ + TN + TP + Bloom, random = ~1|Site, data = data))
summary(lme(Redun ~ HLU + MAT+ YRQ + TN + TP+ Bloom + TRic, random = ~1|Site, data = data))
# ---- full mod ####
alpha_pSEM_randomList <- list(
  lme(TN ~ YRQ + HLU + MAT, random = ~1|Site, data = data),
  lme(TP ~ YRQ + HLU + MAT, random = ~1|Site, data = data),
  lme(Bloom ~ HLU + MAT + YRQ + TN + TP, random = ~1|Site, data = data),
  lme(TRic ~ HLU + MAT + YRQ + TN + TP + Bloom, random = ~1|Site, data = data),
  lme(Redun ~ HLU + MAT + YRQ + TN + TP+ Bloom + TRic, random = ~1|Site, data = data)
)

# Run goodness-of-fit tests
alpha.psem <- as.psem(alpha_pSEM_randomList)

# Evaluate path significance using unstandardized coefficients
coefs(alpha_pSEM_randomList, kelp, standardize = "none")

# Obtain standardized regression coefficients
coefs(alpha_pSEM_randomList, kelp, standardize = "scale")

# Explore individual model fits
alpha.summ <- summary(alpha.psem)
alpha.summ
# ---- best mod ####
data <- na.omit(data)
alpha_pSEM_randomList <- list(
  lme(TN ~ YRQ + HLU + MAT, random = ~1|Site, data = data),
  lme(TP ~ YRQ + HLU + MAT, random = ~1|Site, data = data),
  lme(Bloom ~ HLU + MAT + YRQ + TN + TP, random = ~1|Site, data = data),
  lme(TRic ~ HLU + MAT + YRQ + TN + TP + Bloom, random = ~1|Site, data = data),
  lme(Redun ~ HLU + MAT + YRQ + TN + TP + Bloom + TRic, random = ~1|Site, data = data) )

# Run goodness-of-fit tests
alpha.psem <- as.psem(alpha_pSEM_randomList)
summary(alpha.psem)$Cstat
# Obtain standardized regression coefficients
sink('Data/pSEM/BestMod_stanCoeff_alpha.txt')
coefs(alpha_pSEM_randomList, kelp, standardize = "scale")
sink()
# Explore individual model fits
alpha.summ <- summary(alpha.psem)
sink('Data/pSEM/BestMod_summary_alpha.txt')
alpha.summ
sink()
##########################################################################
### selected LME fit plot
##########################################################################
data <- readRDS('Data/SEMdata_alpha.RDS')
data[,1:9] <- scale(data[,1:9])
head(data)

LMPairs <- matrix(c('HLU','TN',
                    'TN','Bloom',
                    'TN','TRic',
                    'Bloom','TRic',
                    'TRic','Redun'), ncol = 6)
myColor <- c('#FF4F4F','#41AE76','orange','orange','skyblue')

rName <- c('Species Richness', "Redundancy",
           'Human Land Use','Total Nitrogen', "Bloom Intensity")
names(rName) <- c('TRic',"Redun",'HLU','TN','Bloom')

gg <- NULL
for (i in 1:ncol(LMPairs)) {
  xypair <- LMPairs[,i]
  cat(paste0('Processing ', xypair[1],' verus ', xypair[2],'...\n'))
  
  df <- data.frame(x = data[,xypair[1]], y = data[,xypair[2]], Site = data$Site)
  if (i == 2) {
    df <- df[-which.max(df$y),]
  }
  if (i == 4) {
    df <- df[-which.max(df$x),]
  }
  fit <- summary(lme(y ~ x + 1, random = ~1|Site, data = na.omit(df)))
  p <- fit$tTable['x','p-value']
  r <- fit$sigma
  ind <- case_when(p < 0.001 ~ '***',(p >= 0.001)&(p < 0.01) ~ '**', 
                   (p >= 0.01)&(p < 0.05) ~ '*', p >= 0.05 ~ '')
  lab <- paste0('LME R2 = ', round(r,2), ' p = ', round(p,3), ' ',ind)
  if (p < 0.001) {
    lab <- paste0('LME R2 = ', round(r,2), ' p < 0.001 ***')
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
p1 <- ggarrange(gg[[1]],gg[[3]],gg[[4]],gg[[5]], align = 'hv')

pdf('image/02_LME_alpha.pdf', width = 5, height = 4)
p1
dev.off()