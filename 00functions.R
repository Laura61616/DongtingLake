library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggpmisc)
library(Rmisc)
library(raster)
library(rgdal)
library(sf)
library(rasterVis)
library(ncdf4)
library(maptools) 
library(mgcv)
library(piecewiseSEM)
library(nlme)
library(vegan)
library(bcp)
library(mvpart)
library(MVPARTwrap)
library(ggalt)
###########mean with nas####
Mean <- function(x){
  return(mean(x, na.rm = T))
}

SD <- function(x){
  return(sd(x, na.rm = T))
}
#########combine samples from replicates####
CombineSample <- function(otu,iden,fun='sum',thred = 0.75){
  ID <- unique(iden)
  result <- matrix(nrow = nrow(otu),ncol = length(ID))
  for (i in 1:length(ID)) {
    otusub <- as.data.frame(otu[,which(iden==ID[i])])
    
    result[,i] <- apply(otusub, 1, function(x){
      if(fun=='sum'){
        return(sum(x))
      }
      if(fun=='average'){
        if (sum(x>0) >0 ) {
          return(sum(x)/sum(x>0))
        }
        else{
          return(0)
        }
      }
      if(fun=='int_ave'){
        if (sum(x>0) >0 ) {
          return(round(sum(x)/sum(x>0),0))
        }
        else{
          return(0)
        }
      }
      if(fun=='select_sum'){
        if(sum(x>0) > thred*length(x)){
          return(sum(x))
        }else{
          return(0)
        }
      }
    })
    
  }
  rownames(result) <- rownames(otusub);colnames(result) <- ID
  return(result)
}
#########combine taxonomy: fun=c('sum','average','count')####
CombineTaxa <- function(otu,cla,fun='sum'){
  ID <- unique(cla)
  table <- as.data.frame(matrix(ncol = ncol(otu),nrow = length(ID), 
                                dimnames = list(ID,colnames(otu))))
  if(fun%in%c('sum','average')){
    for (i in 1:length(ID)) {
      sub <- otu[which(cla==ID[i]),]
      if(length(which(cla==ID[i]))<2)
        table[i,] <- sub
      if(length(which(cla==ID[i]))>1)
        if(fun=='sum')
          table[i,] <- colSums(sub)
        if(fun=='average')
          table[i,] <- colSums(sub)/nrow(sub)
    }
  }
  if(fun=='count'){
    for (i in 1:length(ID)) {
      sub <- otu[which(cla==ID[i]),]
      table[i,] <- colSums(sub>0)
    }
  }
  
  return(table)
}
#########formating beta div####
FormatBeta <- function(dis){
  tmp <- as.matrix(dis)
  tmp[lower.tri(tmp, diag = T)] <- NA
  tmp <- na.omit(melt(tmp))
  tmp$y1 <- as.numeric( gsub('.*\\.','',tmp$Var1))
  tmp$y2 <- as.numeric( gsub('.*\\.','',tmp$Var2))
  
  tmp <- tmp[tmp$y1==tmp$y2,]
  dt <- data.frame('year' = tmp$y1, 'value' = tmp$value)
  
  beta_df <- data.frame(
    s1 = sapply(strsplit(as.character(tmp$Var1), '\\.'), '[', 2),
    y1 = as.numeric(gsub('.*\\.','',tmp$Var1)),
    s2 = sapply(strsplit(as.character(tmp$Var2), '\\.'), '[', 2),
    y2 = as.numeric(gsub('.*\\.','',tmp$Var2)),
    value = tmp$value)
  beta_use <- beta_df[beta_df$y1==beta_df$y2,c(2,1,3,5)]
  colnames(beta_use) <- c('year','s1','s2','value')
  return(beta_use)
}
#########LinearFitting####
LMLU_Alpha <- function(data, buffer, paras, ind, if.scale = FALSE){

  result <- NULL; gg <- NULL
  use <- data
  for (p in paras) {
    dt1 <- use[,c('site', buffer, ind, p)]
    form <- as.formula(paste0(ind, ' ~ ', p))
    sit <- unique(dt1$site)
    for (s in sit) {
      dt <- subset(dt1, site == s)
      buf <- sort(unique(dt[,buffer]))
      
      for (b in buf) {
        dt.sub <- dt[dt[,buffer] == b,]
        dt.sub <- na.omit(dt.sub)
        if (length(unique(dt.sub[,p])) < 2) {
          next
        }
        if (if.scale) {
          dt.sub[,c(p,ind)] <- scale(dt.sub[,c(p,ind)])
        }
        
        fit <- summary(lm(formula = form, data = dt.sub))
        
        if (nrow(fit$coefficients)==1) {
          next
        }
        p.value <- fit$coefficients[p,'Pr(>|t|)']
        r2.adj <- fit$adj.r.squared
        estimate <- fit$coefficients[p,'Estimate']
        dire <- estimate > 0
        
        add <- data.frame('paras' = p, 'buffer' = b, 'site' = s,
                          'p.value' = p.value, 
                          'r2.adj' = r2.adj, 
                          'estimate' = estimate,
                          'is.positive' = dire)
        
        result <- rbind.data.frame(result, add)
      }
    }
    
    res1 <- result[(result$paras==p),]
    res1$sig <- factor(as.character(res1$p.value < 0.05), levels = c('TRUE','FALSE'))
    
    res1 <- res1[res1$sig == 'TRUE',]
    col.pa <- c('red', 'grey'); names(col.pa) <- c('TRUE','FALSE')
    plo1 <- ggplot(res1, 
                  aes(x = buffer, y = r2.adj))+
      geom_smooth(formula = y ~ x, method = 'loess') + 
      scale_color_manual(values = col.pa, limits = c('TRUE','FALSE')) +
      ylab("Adjusted R Square") + xlab('Radius of Buffer Area')+
      theme_bw()+
      theme(axis.text.x = element_text( angle = 45, hjust = 1),
            legend.position = 'none',
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    plo2 <- ggplot(res1, 
                   aes(x = buffer, y = estimate))+
      geom_smooth(formula = y ~ x, method = 'loess') + 
      scale_color_manual(values = col.pa, limits = c('TRUE','FALSE')) +
      ylab("Estimate") + xlab('Radius of Buffer Area')+
      theme_bw()+
      theme(axis.text.x = element_text( angle = 45, hjust = 1),
            legend.position = 'none',
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    plo <- list('r2' = plo1, 'estimate' = plo2)
    gg$value <- plo
    names(gg)[length(gg)] <- p
    
  }
  final <- list('result' = result, 'fig' = gg)
  return(final)
} 
LMLU <- function(data, group = NA, buffer, paras, ind, if.scale = FALSE){
  gro <- unique(data[,group])
  for (g in gro) {
    cat(paste0('\r\t\tProcessing ', g, '...\r'))
    
    use <- data[data[,group]==g,]
    for (p in paras) {
      dt <- use[,c(buffer, ind, p)]
      
      form <- as.formula(paste0(ind, ' ~ ', p))
      
      buf <- sort(unique(dt[,buffer]))
      
      for (b in buf) {
        dt.sub <- dt[dt[,buffer] == b,]
        if (if.scale) {
          dt.sub[,c(p,ind)] <- scale(dt.sub[,c(p,ind)])
        }
        
        fit <- summary(lm(formula = form, data = dt.sub))
        
        if (nrow(fit$coefficients)==1) {
          next
        }
        p.value <- fit$coefficients[p,'Pr(>|t|)']
        r2.adj <- fit$adj.r.squared
        estimate <- fit$coefficients[p,'Estimate']
        dire <- estimate > 0
        
        add <- data.frame('group' = g,
                          'paras' = p, 'buffer' = b, 
                          'p.value' = p.value, 
                          'r2.adj' = r2.adj, 
                          'estimate' = estimate,
                          'is.positive' = dire)
        
        if (exists("result")) {
          result <- rbind.data.frame(result, add)
        }else{
          result <- add
        }
      }
      
      res1 <- result[(result$group==g)&(result$paras==p),]
      res1$sig <- factor(as.character(res1$p.value < 0.05), levels = c('TRUE','FALSE'))
      
      col.pa <- c('red', 'grey'); names(col.pa) <- c('TRUE','FALSE')
      plo1 <- ggplot(res1, 
                    aes(x = buffer, y = r2.adj))+
        geom_point(aes(color = sig), size = 3) +
        geom_line(color = 'grey', linetype = 2) +
        scale_color_manual(values = col.pa, limits = c('TRUE','FALSE')) +
        ylab("Adjusted R Square") + xlab('Radius of Buffer Area')+
        theme_bw()+
        theme(axis.text.x = element_text( angle = 45, hjust = 1),
              legend.position = 'none',
              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      plo2 <- ggplot(res1, 
                     aes(x = buffer, y = estimate))+
        geom_point(aes(color = sig), size = 3) +
        geom_line(color = 'grey', linetype = 2) +
        scale_color_manual(values = col.pa, limits = c('TRUE','FALSE')) +
        ylab("Estimate") + xlab('Radius of Buffer Area')+
        theme_bw()+
        theme(axis.text.x = element_text( angle = 45, hjust = 1),
              legend.position = 'none',
              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      plo <- list('r2' = plo1, 'estimate' = plo2)
      if (exists('gg')) {
        gg$value <- plo
      }else{
        gg <- list('value' = plo)
      }
      
      names(gg)[length(gg)] <- paste(g, p)
      
    }
  }
  
  final <- list('result' = result, 'fig' = gg)
  return(final)
} 
#########Random Forest####
plot_RF <- function(impor, level = rownames(impor), group, group.level = unique(group)){
  impor.scale <- impor[order(impor$'%IncMSE', decreasing = TRUE), ]
  
  #%IncMSE value
  library(ggplot2)
  
  impor.scale$paras <- rownames(impor.scale)
  impor.scale$paras <- factor(impor.scale$paras, levels = level)
  
  impor.scale$group <- group
  impor.scale$group <- factor(group, levels = group.level)
  
  p <- ggplot() +
    geom_col(data = impor.scale, 
             aes(x = paras, y = `%IncMSE`,fill = group), width = 0.5) +
    labs(title = NULL, x = NULL, y = 'Increase in MSE (%)', fill = NULL) +
    scale_fill_npg() +
    theme(panel.grid = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = 'black')) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, 16))
  
  p
  
  # p<0.05  *, p<0.01  **, p<0.001  ***
  for (paras in rownames(impor.scale)) {
    if (impor.scale[paras,'%IncMSE.pval'] >= 0.05) impor.scale[paras,'%IncMSE.sig'] <- ''
    else if (impor.scale[paras,'%IncMSE.pval'] >= 0.01 & impor.scale[paras,'%IncMSE.pval'] < 0.05) impor.scale[paras,'%IncMSE.sig'] <- '*'
    else if (impor.scale[paras,'%IncMSE.pval'] >= 0.001 & impor.scale[paras,'%IncMSE.pval'] < 0.01) impor.scale[paras,'%IncMSE.sig'] <- '**'
    else if (impor.scale[paras,'%IncMSE.pval'] < 0.001) impor.scale[paras,'%IncMSE.sig'] <- '***'
  }
  
  p <- p +
    geom_text(data = impor.scale, aes(x = paras, y = `%IncMSE`, label = `%IncMSE.sig`), 
              nudge_y = 1)
  
  p
  
  return(p)
}
#########determinant####
AutoDeter <- function(data, form_raw){
  require(performance)
  require(effectsize)
  require(relaimpo)
  require(lmerTest)
  library(nlme)
  library(car)
  ## Checking colinearity of variable with a VIF procedure in a simple LM
  data <- as.data.frame(data)
  tmp <- vif(gls(form_raw,data=data))
  tmp <- names(tmp)[which(tmp>10)]
  
  UsePara <- colnames(data)[2:ncol(data)][!(colnames(data)[2:ncol(data)]%in%tmp)]
  ## step AIC with term of spatial autocorrelation
  modclimF.auto <- gls(as.formula(paste0('BD~',
                                         paste0(UsePara, collapse = ' + '))),data = data)
  ## check residuals
  par(mfrow=c(2,2))
  plot(modclimF.auto)
  graphics.off()
  ## confidence intervals
  r_square <- r2(modclimF.auto)[[1]] #adjusted R2
  effSize <- as.data.frame(effectsize(modclimF.auto))
  ## output
  tmp <- summary(modclimF.auto)$tTable
  out1 <- data.frame(effSize[-1,], 'p_value' = tmp[-1,4])
  colnames(out1) <- c('para','size','CI','CI_low','CI_high','p_value')
  
  result <- list('r2' = r_square, 'effSize' = out1)
  return(result)
}

#########double y-axis with ggplot2####
double_y_axis <- function(p1, p2){
  library(gtable)
  library(grid)
  g1 <- ggplot_gtable(ggplot_build(p1))
  g2 <- ggplot_gtable(ggplot_build(p2))
  
  # overlap the panel of 2nd plot on that of 1st plot
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
  
  # axis tweaks
  ia <- which(g2$layout$name == "axis-l")
  ga <- g2$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
  g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
  
  # draw it
  # grid.draw(g)
  return(g)
}
