source('00functions.R')

##########################################################################
### change point
##########################################################################
# --- BD & --- Env####
BD <- readRDS('Data/UseBDIndices_alpha.RDS')
gamma <- readRDS('Data/UseBDIndices_gamma.RDS')
# --- Env
env <- readRDS('Data/Env_use_alpha.RDS')
UseEnv <- na.omit(env)
# --- combine
ids <- union(rownames(BD), rownames(UseEnv))
data <- cbind.data.frame(BD[ids,],UseEnv[ids,])
rownames(data) <- ids
data$Bloom <- log10(data$Bloom+1)
{
  rName <- c('Species Richness', "Functional Redundancy",colnames(UseEnv)[1:6], 
             'Log10 Bloom Intensity')
  names(rName) <- c('TRic',"Redun",colnames(UseEnv))
  myColor <- c('orange','skyblue',rep('#FF4F4F',ncol(UseEnv)))
}
# --- calculating
paras <- names(rName)
data$year <- as.numeric(gsub('.*\\.','',rownames(data)))

gg2 <- NULL
for (i in 1:length(paras)) {
  pa <- paras[i]
  cat('\r\tprocessing ', pa, '...\r')
  
  add <- na.omit(data[,c('year', pa)]); colnames(add)[2] <- 'value'
  add <- unique(add)
  # ploting
  g2 <- ggplot(add, aes(x = year, y = value))+
    geom_smooth(method = 'gam', formula = y ~ s(x, k=20), level = 0.95,
                lwd = 1, alpha = 0.1, fill = myColor[i], color = myColor[i])+
    stat_poly_eq(aes(
      label =  paste(after_stat(p.value.label), 
                      after_stat(adj.rr.label), 
                      sep = "*\" \"*")),
      formula = y ~ x, parse = TRUE, 
      label.y = "top", label.x = "left",
      p.digits = 3)+
    ylab(rName[pa])+xlab('')+
    theme_bw()+
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 45,hjust=1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = 'grey'))
  g2
  if (i>2) {
    g2 <- g2 + coord_flip()
  }else{
    tmp <- gamma; colnames(tmp)[colnames(tmp)==pa] <- 'value'
    g2 <- g2 + 
      geom_point(alpha = 0.7, 
                 size = 2, 
                 fill = myColor[i], 
                 color = myColor[i])+
      geom_line(aes(x = year, y = value), 
                data = tmp, 
                linetype = 2,
                inherit.aes = F,
                color = 'grey')
  }
  g2
  # keep record
  gg2$value <- g2
  names(gg2)[i] <- pa
}
# --- beta####
beta <- readRDS('Data/SEMdata_beta.RDS')$All
{
  rName <- c('Taxonomic Dissimilarity', "Functional Dissimilarity")
  names(rName) <- c('taxaDiss',"funcDiss")
  myColor <- c('#41AE76','#FED439')
}
paras <- names(rName)

i=1
for (i in 1:length(paras)) {
  pa <- paras[i]
  cat('\r\tprocessing ', pa, '...\r')
  
  add <- na.omit(beta[,c('year', pa)]); colnames(add)[2] <- 'value'
  add <- unique(add)
  # ploting
  g2 <- ggplot(add, aes(x = year, y = value))+
    geom_smooth(method = 'gam', formula = y ~ s(x, k=20), level = 0.95,
                lwd = 1, alpha = 0.1, fill = myColor[i], color = myColor[i])+
    geom_point(alpha = 0.7, size = 2, fill = myColor[i], color = myColor[i])+
    stat_poly_eq(aes(
      label =  paste(after_stat(p.value.label), 
                     after_stat(adj.rr.label), 
                     sep = "*\" \"*")),
      formula = y ~ x, parse = TRUE, 
      label.y = "top", label.x = "left",
      p.digits = 3)+
    ylab(rName[pa])+xlab('')+
    theme_bw()+
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 45,hjust=1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = 'grey'))
  # keep record
  gg2$value <- g2
  names(gg2)[length(gg2)] <- pa
}

# --- ploting####
g1 <- ggarrange(gg2[[1]],gg2[[2]],gg2[[10]],gg2[[11]],align = 'hv', ncol = 4)
g1

g2 <- ggarrange(gg2[[3]],gg2[[4]], gg2[[6]],gg2[[7]],gg2[[8]],gg2[[9]],
                align = 'hv', ncol = 6)
g2

g <- ggarrange(g2, g1,nrow = 2, align = 'hv')
g
pdf('image/01_TrendPlot.pdf', width = 14, height = 6)
g
dev.off()
