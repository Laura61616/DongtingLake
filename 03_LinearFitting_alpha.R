source('00functions.R')
##########################################################################
### plot multiple buffer
##########################################################################
prefix <-'landDT/buffer/point/'
pf <- grep('.shp$',list.files(prefix),value = T)
size <- as.numeric(gsub('.*_','',gsub('\\.shp','',pf)))
pf <- pf[order(size)]; size <- sort(size)
for (i in 1:length(pf)) {
  buf <- st_read(paste0(prefix, pf[i]))
  
  if (i==1) {
    buffer <- list('value' = buf)
  }else{
    buffer$value <- buf
  }
  names(buffer)[i] <- size[i]
}
# --- loading buffers
bf <- 'landDT/buffer/DongtingL.shp'
basin <- st_read(bf)
# --- loading basin
my_color <- my_pal("mixed")(18)

p <- ggplot()+
  labs(x ='Longitude(E)',y="Latitude(N)")+
  geom_sf(aes(geometry = `geometry`),data = basin, 
          fill = 'skyblue', color = 'black', lwd = 0.05) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

for (i in length(buffer):1) {
  p <- p +
    geom_sf(aes(geometry = `geometry`),data = buffer[[i]][1,], 
            color = 'grey', fill = my_color[i], lwd = 0.05, alpha = 0.5)
}

p0 <- p + xlim(112.5, 113.5) + ylim(29, 29.6)
p0
p1 <- p + xlim(112.8, 113) + ylim(29.1, 29.6) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
p1
p2 <- p + xlim(112.98, 113) + ylim(29.32, 29.36)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
p2

pp <- ggarrange(p0, ggarrange(p2, p1, nrow = 2, ncol = 1, align = 'hv'), 
                widths = c(3/4, 1/4))
pp
pdf('image/03_Buffer_point.pdf')
pp
dev.off()
##########################################################################
### Linear Fitting plus Null: env
##########################################################################
# --- loading Null model
Null <- readRDS('landDT/output/NullModelAlphaDat_stan.RDS')
# --- loading original landuse data
LU <- readRDS('landDT/output/LandUseDat_Alpha.RDS')
useLU <- LU[,c(1,2,3,6,7,5,8,9)]; colnames(useLU) <- colnames(Null$`1`)
# --- loading BD data
BD <- readRDS("Data/SEMdata_alpha.RDS")
BD$year <- as.numeric(gsub('.*\\.','',rownames(BD)))
useBD <- BD[,c('year', 'Site','TN','Bloom')]
# --- LMLU for original data
dat <- useBD
indice <- c('TN','Bloom')
for (j in 1:length(indice)) {
  remove(fig); remove(gg)
  
  ind <- indice[j]
  cat(paste0('\tProcessing ', ind, '...\n'))
  
  dat_tmp1 <- dat; dat_tmp1$id <- paste(dat$year, dat$Site)
  useLU_tmp1 <- useLU; useLU_tmp1$id <- paste(useLU$year, useLU$site)
  data <- merge(dat_tmp1[,c("id", "year", "Site", ind)], useLU_tmp1, 
                by = 'id')
  data <- data[,c(5,2,3,4,8:12)]; colnames(data)[2:3] <- c('year','site')
  data$HLU <- data$Cropland + data$Impervious
  
  buffer <- "size"
  paras <- 'HLU'
  
  res <- LMLU_Alpha(data,  buffer, paras, ind, if.scale = TRUE)
  
  tmp <- res$result
  add <- data.frame('type' = 'original', 'indice' = ind, tmp, 
                    'sig' = tmp$p.value < 0.05)
  if(j==1){
    final <- add
  }else{
    final <- rbind.data.frame(final, add)
  }
  
  fig <- res$fig
}
head(final); summary(final$sig)
write.table(final, 'Data/LM/LandLinearFitting_Env_original.csv', sep = ',', row.names = F,
            quote = F)
# --- LULU for Null model
write.table(t(c("type","indice","paras","buffer","site","p.value","r2.adj","estimate",
                "is.positive","sig")), 'Data/LM/LandLinearFitting_Env_NULL.csv', 
            row.names = F, col.names = F,
            sep = ',', quote = F)
for (i in 1:1200) {
  useLU <- Null[[i]]
  cat(paste0('Processing Null Model ', i, '...\n'))
  
  for (j in 1:length(indice)) {
    ind <- indice[j]
    cat(paste0('\tProcessing ', ind, '...\n'))
    
    dat_tmp1 <- dat; dat_tmp1$id <- paste(dat$year, dat$Site)
    useLU_tmp1 <- useLU; useLU_tmp1$id <- paste(useLU$year, useLU$site)
    data <- merge(dat_tmp1[,c("id", "year", "Site", ind)], useLU_tmp1, 
                  by = 'id')
    data <- data[,c(5,2,3,4,8:12)]; colnames(data)[2:3] <- c('year','site')
    data$HLU <- data$Cropland + data$Impervious
    
    buffer <- "size"
    paras <- 'HLU'
    
    res <- LMLU_Alpha(data,  buffer, paras, ind, if.scale = TRUE)
    
    tmp <- res$result
    add <- data.frame('type' = i, 'indice' = ind, tmp, 'sig' = tmp$p.value < 0.05)
    write.table(add, 'Data/LM/LandLinearFitting_Env_NULL.csv', 
                row.names = F, col.names = F,
                sep = ',', quote = F, append = T)
    
    fig <- res$fig
  }
}
##########################################################################
### Linear Fitting plus Null: alpha
##########################################################################
# --- loading Null model
Null <- readRDS('landDT/output/NullModelAlphaDat_stan.RDS')
# --- loading original landuse data
LU <- readRDS('landDT/output/LandUseDat_Alpha.RDS')
useLU <- LU[,c(1,2,3,6,7,5,8,9)]; colnames(useLU) <- colnames(Null$`1`)
# --- loading BD data
BD <- readRDS("Data/AllBDIndices.RDS")
AllBD <- BD$Alpha$Dat
useBD <- AllBD[,c('commu','year', 'site','sr','Redundancy')]
colnames(useBD)[4:ncol(useBD)] <- c('SR','Redun')
# --- LMLU for original data
dat <- subset(useBD, commu == 'BS')
indice <- c('SR','Redun')

for (j in 1:length(indice)) {
  remove(fig); remove(gg)
  
  ind <- indice[j]
  cat(paste0('\tProcessing ', ind, '...\n'))
  
  dat_tmp1 <- dat; dat_tmp1$id <- paste(dat$year, dat$site)
  useLU_tmp1 <- useLU; useLU_tmp1$id <- paste(useLU$year, useLU$site)
  data <- merge(dat_tmp1[,c("id", "commu", "year", "site", ind)], useLU_tmp1, 
                by = 'id')
  data <- data[,c(6,2,3,4,5,9:13)]; colnames(data)[3:4] <- c('year','site')
  data$HLU <- data$Cropland + data$Impervious
  
  buffer <- "size"
  paras <- 'HLU'
  
  res <- LMLU_Alpha(data,  buffer, paras, ind, if.scale = TRUE)
  
  tmp <- res$result
  add <- data.frame('type' = 'original', 'indice' = ind, tmp, 
                    'sig' = tmp$p.value < 0.05)
  if(j==1){
    final <- add
  }else{
    final <- rbind.data.frame(final, add)
  }
}
head(final); summary(final$sig)
write.table(final, 'Data/LM/LandLinearFitting_BD_original.csv', sep = ',', row.names = F,
            quote = F)
# --- LULU for Null model
dat <- subset(useBD, commu == 'BS')
indice <- c('SR','Redun')

write.table(t(c("type","indice","paras","buffer","site","p.value","r2.adj","estimate",
                "is.positive","sig")), 
            'Data/LM/LandLinearFitting_BD_NULL.csv', 
            row.names = F, col.names = F,
            sep = ',', quote = F)
for (i in 1:1200) {
  useLU <- Null[[i]]
  cat(paste0('Processing Null Model ', i, '...\n'))
  
  for (j in 1:length(indice)) {
    ind <- indice[j]
    cat(paste0('\tProcessing ', ind, '...\n'))
    
    dat_tmp1 <- dat; dat_tmp1$id <- paste(dat$year, dat$site)
    useLU_tmp1 <- useLU; useLU_tmp1$id <- paste(useLU$year, useLU$site)
    data <- merge(dat_tmp1[,c("id", "commu", "year", "site", ind)], useLU_tmp1, 
                  by = 'id')
    data <- data[,c(6,2,3,4,5,9:13)]; colnames(data)[3:4] <- c('year','site')
    data$HLU <- data$Cropland + data$Impervious
    
    buffer <- "size"
    paras <- 'HLU'
    
    res <- LMLU_Alpha(data,  buffer, paras, ind, if.scale = TRUE)
    
    tmp <- res$result
    add <- data.frame('type' = i, 'indice' = ind, tmp, 'sig' = tmp$p.value < 0.05)
    write.table(add, 'Data/LM/LandLinearFitting_BD_NULL.csv', 
                row.names = F, col.names = F,
                sep = ',', quote = F, append = T)
    if(j==1){
      final <- add
    }else{
      final <- rbind.data.frame(final, add)
    }
  }
}
##########################################################################
### ploting
##########################################################################
env <- read.csv('Data/LM/LandLinearFitting_Env_original.csv')
env_null <- read.csv('Data/LM/LandLinearFitting_Env_NULL.csv')

bd <- read.csv('Data/LM/LandLinearFitting_BD_original.csv')
bd_null <- read.csv('Data/LM/LandLinearFitting_BD_NULL.csv')


final_a <- rbind.data.frame(env,bd)
null_a <- rbind.data.frame(env_null,bd_null)

saveRDS(final_a, 'Data/LM/LandLinearFitting_Alpha_original2.RDS')
saveRDS(null_a, 'Data/LM/LandLinearFitting_Alpha_null2.RDS')
# --- alpha
my_colors <- c('#FF4F4F','#41AE76','orange','skyblue')

final_a <- readRDS('Data/LM/LandLinearFitting_Alpha_original2.RDS')
null_a <- readRDS('Data/LM/LandLinearFitting_Alpha_null2.RDS')
{
  rName <- c('Species Richness', "Redundancy","Total Nitrogen","Bloom Intensity")
  names(rName) <- c("SR","Redun","TN",'Bloom')
}
dat_alpha <- subset(final_a, paras == 'HLU')
null_alpha<- subset(null_a, paras == 'HLU')

indice <- unique(final_a$indice)

pp <- NULL

xintercept <- data.frame(
  '1' = c(300, 7000),
  '2' = c(700, 3000),
  '3' = c(100, 900),
  '4' = c(NA, 1500))
for (i in 1:length(indice)) {
  ind <- indice[i]
  cat(paste0('\r\tProcessing ', ind, '...\r'))
  
  dat <- subset(dat_alpha, (indice == ind)&(p.value < 0.1))
  dat$is.positive <- factor(dat$is.positive, levels = c('TRUE','FALSE'))
  
  dat2 <- subset(null_alpha, (indice == ind)&(p.value < 0.1))
  dat2 <- aggregate(dat2$r2.adj,  by = list('buffer' = dat2$buffer), mean)

  xinter <- xintercept[,i]
  
  p1 <- ggplot(dat, 
              aes(x = log10(buffer), y = r2.adj))+
    # geom_point(aes(shape = is.positive), color = my_colors[i]) +
    geom_smooth(formula = y ~ s(x, k=5), method = 'gam',
                fill = my_colors[i], color = my_colors[i], alpha = 0.2) +
    geom_line(mapping = aes(x = log10(buffer), y = x), 
              color = 'grey',  linetype = 2,
              data = dat2, inherit.aes = F) +
    # scale_shape_manual(values = c(16,1)) +
    scale_x_continuous(breaks = c(2,3,4), labels = c(100, 1000, 10000)) +
    ylab("% Explained Variance (R2)") + xlab('Radius of Buffer Area')+
    theme_bw()+
    theme(axis.text.x = element_text( angle = 45, hjust = 1),
          legend.position = 'bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  dat3 <- subset(null_alpha, (indice == ind)&(p.value < 0.1))
  dat3 <- aggregate(dat3$estimate,  by = list('buffer' = dat3$buffer), mean)
  p2 <- ggplot(dat, 
               aes(x = log10(buffer), y = abs(estimate)))+
    # geom_point(aes(shape = is.positive), color = my_colors[i]) +
    geom_smooth(formula = y ~ s(x, k=5), method = 'gam',
                fill = my_colors[i], color = my_colors[i], alpha = 0.2) +
    geom_line(mapping = aes(x = log10(buffer), y = abs(x)), 
              color = 'grey',  linetype = 2,
              data = dat3, inherit.aes = F) +
    # scale_shape_manual(values = c(16,1)) +
    scale_x_continuous(breaks = c(2,3,4), labels = c(100, 1000, 10000)) +
    ylab("Human land-use effect (estimate)") + xlab('Radius of Buffer Area')+
    theme_bw()+
    theme(axis.text.x = element_text( angle = 45, hjust = 1),
          legend.position = 'bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p2
  for (j in xinter) {
    p1 <- p1 + geom_vline(xintercept = log10(j), color = my_colors[i], linetype = 2)
    p2 <- p2 + geom_vline(xintercept = log10(j), color = my_colors[i], linetype = 2)
  }
  p <- list('r2' = p1, 'estimate' = p2)
 
  pp$value <- p
  
  names(pp)[length(pp)] <- paste0('alpha_',ind)
}

names(pp)
g1 <- ggarrange(pp[[1]][[1]],pp[[2]][[1]], pp[[3]][[1]], pp[[4]][[1]],
                ncol = 1, nrow = 4, align = 'hv')
g2 <- ggarrange(pp[[1]][[2]],pp[[2]][[2]], pp[[3]][[2]], pp[[4]][[2]],
                ncol = 1, nrow = 4, align = 'hv')

pdf('image/03_LinearFitting_alpha.pdf', height = 11, width = 4)
g1
dev.off()

saveRDS(g, 'Data/LM/Fig_alpha.RDS')
