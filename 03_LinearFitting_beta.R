source('00functions.R')
##########################################################################
### prepare landuse data for beta div
##########################################################################
Null <- readRDS('landDT/output/NullModelAlphaDat_stan.RDS')
LU <- readRDS('landDT/output/LandUseDat_Alpha.RDS')
##loading BD data
BD <- readRDS("Data/SEMdata_beta.RDS")$StanNAfree
useBD <- BD[,c('s1','s2','y')]
colnames(useBD)[3] <- 'year'
##loading original landuse data
useLU <- LU[,c(1,2,3,6,7,5,8,9)]; colnames(useLU) <- colnames(Null$`1`)
useLU$HLU <- useLU$Cropland + useLU$Impervious

pairs <- unique(useBD[,c('s1','s2','year')])
s = 50; i = 1

final <- NULL
for (s in unique(input$size)) {
  subDT <- subset(useLU, size == s)
  cat(paste0('[processing] ', s, '...\n'))
  for (i in 1:nrow(pairs)) {
    index1 <- which(paste(subDT$site, subDT$year)==paste(pairs[i,'s1'], pairs[i,'year']))
    index2 <- which(paste(subDT$site, subDT$year)==paste(pairs[i,'s2'], pairs[i,'year']))
    
    if ((length(index1)==0)|(length(index2)==0)) {
      next
    }
    add <- data.frame('size' = s, 
                      pairs[i,],
                      'HLU' = abs(subDT[index1,'HLU']-subDT[index2,'HLU']))
    final <- rbind.data.frame(final, add)
  }
}
saveRDS(final, 'landDT/output/LandUseDat_Beta.RDS')

##loading null landuse data
j=1

final <- NULL
for (j in 1:length(Null)) {
  cat(paste0('[processing] ', j, '...\n'))
  useLU <- Null[[j]]
  useLU$HLU <- useLU$Cropland + useLU$Impervious
  
  pairs <- unique(useBD[,c('s1','s2','year')])
  s = 50; i = 1
  
  ff <- NULL
  for (s in unique(input$size)) {
    subDT <- subset(useLU, size == s)
    cat(paste0('\t[processing] ', s, '...\n'))
    for (i in 1:nrow(pairs)) {
      index1 <- which(paste(subDT$site, subDT$year)==paste(pairs[i,'s1'], pairs[i,'year']))
      index2 <- which(paste(subDT$site, subDT$year)==paste(pairs[i,'s2'], pairs[i,'year']))
      
      if ((length(index1)==0)|(length(index2)==0)) {
        next
      }
      add <- data.frame('size' = s, 
                        pairs[i,],
                        'HLU' = abs(subDT[index1,'HLU']-subDT[index2,'HLU']))
      ff <- rbind.data.frame(ff, add)
    }
  }
  final$add <- ff; names(final)[j] <- j
}
saveRDS(final, 'landDT/output/NullModelBetaDat_stan.RDS')
##########################################################################
### Linear Fitting plus Null: beta
##########################################################################
# --- loading landuse data
Null <- readRDS('landDT/output/NullModelBetaDat_stan.RDS')
LU <- readRDS('landDT/output/LandUseDat_Beta.RDS')
# --- loading BD data
BD <- readRDS("Data/SEMdata_beta.RDS")$StanNAfree
useBD <- BD[,c('s1','s2','y','TP','Bloom','taxaDiss','funcDiss')]
colnames(useBD)[3] <- 'year'
# --- LMLU for original data
dat <- useBD
dat$id <- paste(dat$s1, dat$s2, dat$year)
indice <- c('TP','Bloom','taxaDiss','funcDiss')

useLU <- LU
useLU$id <- paste(useLU$s1, useLU$s2, useLU$year)

j = 1
for (j in 1:length(indice)) {
  remove(fig); remove(gg)
  
  ind <- indice[j]
  cat(paste0('\tProcessing ', ind, '...\n'))
  
  data <- merge(dat, useLU, by = 'id', all = T)
  data <- data[, c('size','s1.x','s2.x','year.x',indice,'HLU')]
  colnames(data)[2:4] <- c('s1','s2','year')
  data$commu <- '1'
  
  group <- "commu"; buffer <- "size"
  paras <- 'HLU'
  
  res <- LMLU(data, group, buffer, paras, ind, if.scale = TRUE)
  
  tmp <- res$result
  add <- data.frame('type' = 'original', 'indice' = ind, tmp, 
                    'sig' = case_when(
                      tmp$p.value < 0.05 ~ '1',
                      (tmp$p.value < 0.1)&(tmp$p.value >= 0.05) ~ '2',
                      tmp$p.value >= 0.1 ~ '3'
                    ))
  if(j==1){
    final <- add
  }else{
    final <- rbind.data.frame(final, add)
  }
}
head(final); summary(final$sig)
write.table(final, 'Data/LM/LandLinearFitting_Beta_original.csv', 
            sep = ',', row.names = F,
            quote = F)
# --- LULU for Null model
write.table(t(c("type","indice","group","paras","buffer","site",
                "p.value","r2.adj","estimate",
                "is.positive","sig")), 
            'Data/LM/LandLinearFitting_Beta_NULL.csv', 
            row.names = F, col.names = F,
            sep = ',', quote = F)

i <- 1
for (i in 1:1200) {
  useLU <- Null[[i]]
  useLU$id <- paste(useLU$s1, useLU$s2, useLU$year)
  
  cat(paste0('Processing Null Model ', i, '...\n'))
  
  j <- 1
  for (j in 1:length(indice)) {
    ind <- indice[j]
    cat(paste0('\tProcessing ', ind, '...\n'))
    
    remove(fig); remove(gg)
    
    ind <- indice[j]
    cat(paste0('\tProcessing ', ind, '...\n'))
    
    data <- merge(dat, useLU, by = 'id', all = T)
    data <- data[, c('size','s1.x','s2.x','year.x',indice,'HLU')]
    colnames(data)[2:4] <- c('s1','s2','year')
    data$commu <- '1'
    
    group <- "commu"; buffer <- "size"
    paras <- 'HLU'
    
    res <- LMLU(data, group, buffer, paras, ind, if.scale = T)
    
    tmp <- res$result
    add <- data.frame('type' = as.character(j), 'indice' = ind, tmp, 
                      'sig' = case_when(
                        tmp$p.value < 0.05 ~ '1',
                        (tmp$p.value < 0.1)&(tmp$p.value >= 0.05) ~ '2',
                        tmp$p.value >= 0.1 ~ '3'
                      ))
    
    write.table(add, 'Data/LM/LandLinearFitting_Beta_NULL.csv', 
                row.names = F, col.names = F,
                sep = ',', quote = F, append = T)
  }
}
##########################################################################
### ploting
##########################################################################
library(mgcv)

bd <- read.csv('Data/LM/LandLinearFitting_Beta_original.csv')
bd_null <- read.csv('Data/LM/LandLinearFitting_Beta_NULL.csv')

saveRDS(bd, 'Data/LM/LandLinearFitting_Beta_original.RDS')
saveRDS(bd_null, 'Data/LM/LandLinearFitting_Beta_null.RDS')
# --- beta
my_colors <- c('#FF4F4F','#41AE76','orange','skyblue')

final_a <- readRDS('Data/LM/LandLinearFitting_Beta_original.RDS')
null_a <- readRDS('Data/LM/LandLinearFitting_Beta_null.RDS')
{
  rName <- c('Taxonomic dissimilarity', "Functional dissimilarity",
             "Total Nitrogen","Bloom Intensity")
  names(rName) <- c("TN",'Bloom',"taxaDiss","funcDiss")
}
dat_alpha <- final_a
null_alpha<- null_a

indice <- unique(final_a$indice)

pp <- NULL
xintercept <- c(700,15000,50,50)
i=1
for (i in 1:length(indice)) {
  ind <- indice[i]
  cat(paste0('\r\tProcessing ', ind, '...\r'))
  
  dat <- subset(dat_alpha, indice == ind)
  dat$is.positive <- factor(dat$is.positive, levels = c('TRUE','FALSE'))
  dat$sig <- as.factor(dat$sig)
  
  dat2 <- subset(null_alpha, indice == ind)
  dat2 <- aggregate(dat2$r2.adj,  by = list('buffer' = dat2$buffer), mean)

  xinter <- xintercept[i]
  
  p1 <- ggplot(dat,
              aes(x = log10(buffer), 
                  shape = sig,
                  y = abs(r2.adj)))+
    geom_point(size = 2, color = my_colors[i], fill = my_colors[i]) +
    geom_line(data = dat,
              mapping = aes(x = log10(buffer), 
                            y = abs(r2.adj)),
              inherit.aes = F,
              color = my_colors[i]) +
    # geom_line(mapping = aes(x = log10(buffer),
    #                         y = x),
    #           color = 'grey',
    #           data = dat2, inherit.aes = F) +
    scale_x_continuous(breaks = c(2,3,4), labels = c(100, 1000, 10000)) +
    scale_shape_manual(values = c(16,16,1)) +
    geom_vline(xintercept = log10(xinter), color = my_colors[i], linetype = 2) +
    ylab("Explained variance (R2)") + xlab('Radius of Buffer Area')+
    theme_bw() +
    theme(axis.text.x = element_text( angle = 45, hjust = 1),
          legend.position = 'none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p1
  
  # dat2 <- subset(null_alpha, indice == ind)
  # dat2 <- aggregate(dat2$estimate,  by = list('buffer' = dat2$buffer), mean)
  p2 <- ggplot(dat,
               aes(x = log10(buffer), 
                   # shape = is.positive,
                   y = abs(estimate)))+
    geom_point(size = 2, color = my_colors[i], fill = my_colors[i]) +
    geom_line(data = dat,
              mapping = aes(x = log10(buffer), 
                            y = abs(estimate)),
              inherit.aes = F,
              color = my_colors[i]) +
    # geom_line(mapping = aes(x = log10(buffer),
    #                         y = x),
    #           color = 'grey',
    #           data = dat2, inherit.aes = F) +
    scale_x_continuous(breaks = c(2,3,4), labels = c(100, 1000, 10000)) +
    # scale_shape_manual(values = c(16,17)) +
    geom_vline(xintercept = log10(xinter), color = my_colors[i], linetype = 2) +
    ylab("Human land-use effect (estimate)") + xlab('Radius of Buffer Area')+
    theme_bw() +
    theme(axis.text.x = element_text( angle = 45, hjust = 1),
          legend.position = 'none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p2
  p <- list('r2' = p1, 'estimate' = p2)
  pp$value <- p
  
  names(pp)[length(pp)] <- paste0('beta_',ind)
}

names(pp)
g1 <- ggarrange(pp[[1]][[1]],pp[[2]][[1]], pp[[3]][[1]], pp[[4]][[1]],
                ncol = 1, nrow = 4, align = 'hv')
g1
g2 <- ggarrange(pp[[1]][[2]],pp[[2]][[2]], pp[[3]][[2]], pp[[4]][[2]],
                ncol = 1, nrow = 4, align = 'hv')
g2

pdf('image/03_LinearFitting_beta.pdf', height = 11, width = 4)
g1
dev.off()

saveRDS(pp, 'Data/LM/Fig_beta.RDS')
