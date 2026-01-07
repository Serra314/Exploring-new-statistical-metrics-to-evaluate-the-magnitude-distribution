
source('magtest_fun.R')
library(inlabru)
library(doParallel)
library(foreach)
library(dplyr)


  
  #c("Brier.score", "Chi.squared.score", "Chi.squared.score2", "Hellinger.stat",
                  #  "Multinom.loglik.ratio", "TAPGR.Score.ECDF.plot", "TapGRM", "M0") 

TapGRM <- function(n, b, corner_mag, m_min){
  ## Convert Mc to moment and b to beta
  Mt <- 10**((3/2)*m_min + 9.1)
  Mc <- 10**((3/2)*corner_mag + 9.1)
  beta <- (2/3)*b
  
  R1 <- runif(n, 0, 1)
  R2 <- runif(n, 0, 1)
  
  M1 <- Mt*R1**(-1/beta)
  M2 <- Mt - Mc*log(R2)
  
  # Pick minimum value of two sampled options
  mom <- pmin(M1, M2)
  
  ## convert moments back to magnitudes
  ms <- (2/3)*(log10(mom) - 9.1)
  
  return(ms)
}


name.fun.topar <-  c(unlist(strsplit(lsf.str()," *: *")), "TapGRM", "%do%", 'foreach', 'ad.test')


# This is with TAP-GR law
TAPGR.Score.ECDF.plot <- function(b.fore, 
                                  mc.fore, 
                                  type.fore = 'TapGR',
                                  b.obs = 1, 
                                  mc.obs = 6,
                                  type.obs = 'TapGR',
                                  N.o = 100, # expected number of observations
                                  n.fore = 1000, # number of catalogues per period
                                  n.rep = 100, # number of independent periods
                                  M.max = 15,
                                  M0 = 2.5,
                                  return.plot = TRUE){
  
  # create forecast
  m.fore <- foreach(i = 1:n.fore, .export = name.fun.topar) %do% {
    n.p <- rpois(1, N.o)
    M0 <- M0
    if(type.fore == 'TapGR'){
    
      data.frame(m.ij = TapGRM(n = n.p, 
                               b = b.fore, 
                               corner_mag = mc.fore,
                               m_min = M0)) 
    } else if(type.fore == 'GR'){
      beta.fore <- b.fore*log(10)
      data.frame(m.ij = rexp(n.p, beta.fore) + M0)
    } else{
      stop('Uknown type.fore argument, please choose between `TapGR` and `GR`')
    }
    
  }
  # merge forecast
  m.u <- bind_rows(m.fore)
  N.u <- length(m.u$m.ij)
  # create breaks for binning
  m.breaks <- seq((M0 - 0.05), (M.max + 0.05), by = 0.1)
  
  # create bin counts of merged forecast
  h.u <- hist(m.u$m.ij, plot = FALSE, breaks = m.breaks)
  # this calculate the score of a sample from the merged forecast
  boot.score.sample <- foreach(i = 1:n.fore, .combine = rbind, .export = name.fun.topar) %do% {
    M0 <- M0
    m.j <- sample(m.u$m.ij, N.o, replace = TRUE)
    data.frame(chisq.score = Chi.squared.score(m.u$m.ij, m.j, M.max = 12),
               chisq.score_p1 = Chi.squared.score_plus1(m.u$m.ij, m.j, M.max = 12),
               csep.score = m.j.score(m.u$m.ij, m.j, M.max = 12, N.obs = N.o),
               csep.score.un = m.j.score(m.u$m.ij, m.j, M.max = 12, N.obs = N.o, type = 'alt'),
               multinom.score = Multinom.loglik.ratio(m.u$m.ij, m.j),
               hellinger = Hellinger.stat(m.u$m.ij, m.j, m.breaks),
               brier = Brier.score(m.u$m.ij, m.j, m.breaks)
    )
  }
  cat('Finished bootstrap score distribution \n')
  emp.score.distro <- lapply(1:ncol(boot.score.sample), \(x) 
                             ecdf(boot.score.sample[,x]))
  names(emp.score.distro) <- colnames(boot.score.sample)
  # observed score values
  obs.score <- foreach(i = 1:n.rep, .combine = rbind, .export = name.fun.topar) %do% {
    M0 <- M0
    if(type.obs == 'TapGR'){
      m.obs <- TapGRM(n = N.o, 
                      b = b.obs, 
                      corner_mag = mc.obs,
                      m_min = M0) 
     
    } else if(type.obs == 'GR'){
      beta.obs <- b.obs*log(10)
      m.obs <- rexp(N.o, beta.obs) + M0
     
    } else{
      stop('Uknown type.obs argument, please choose between `TapGR` and `GR`')
    }
    h.obs <- hist(m.obs, breaks = m.breaks, plot = FALSE)
    obs.chisq.score = Chi.squared.score(m.u$m.ij, m.obs, M.max = 12)
    obs.chisq.score_p1 = Chi.squared.score_plus1(m.u$m.ij, m.obs, M.max = 12)
    obs.csep.score = m.j.score(m.u$m.ij, m.obs, M.max = 12, N.obs = N.o)
    obs.csep.score.un = m.j.score(m.u$m.ij, m.obs, M.max = 12, N.obs = N.o, type = 'alt')
    obs.multinom.score = Multinom.loglik.ratio(m.u$m.ij, m.obs)
    obs.hellinger = Hellinger.stat(m.u$m.ij, m.obs, m.breaks)
    obs.brier = Brier.score(m.u$m.ij, m.obs, m.breaks)
    data.frame(Chisq = emp.score.distro$chisq.score(obs.chisq.score),
               Chisq.p1 = emp.score.distro$chisq.score_p1(obs.chisq.score_p1),
               Mtest = emp.score.distro$csep.score(obs.csep.score),
               Mtest.un = emp.score.distro$csep.score.un(obs.csep.score.un),
               Multinomial = emp.score.distro$multinom.score(obs.multinom.score),
               Hellinger = emp.score.distro$hellinger(obs.hellinger),
               Brier = emp.score.distro$brier(obs.brier)
    )
  }
  print('iam out')
  xx = seq(0, 1, 0.01)
  df.plot <- foreach(i = 1:ncol(obs.score), .combine = rbind, .export = name.fun.topar) %do% {
    data.frame(value = xx, 
               ECDF = sapply(xx, \(x) mean(obs.score[,i] <= x)), 
               score = colnames(obs.score)[i])
  }
  
  # increasing size of m.obs should help - also I can try it taking a multinomial with prob
  # given by the real theoretical frequencies.
  df.boxplot <- foreach(i = 1:ncol(obs.score), .combine = rbind,.export = name.fun.topar) %do% {
    data.frame(value = obs.score[,i], 
               score = colnames(obs.score)[i])
  }
  
  if(return.plot){
    return(list(boxplot = ggplot(df.boxplot, aes(x = score, y = value)) + 
                  geom_boxplot(),
                ECDF = ggplot(df.plot, aes(x = value, y = ECDF, color = score,
                                           linetype = score)) + geom_line() + 
                  geom_abline(slope = 1, intercept = 0, linetype = 2),
                df.boxplot = df.boxplot,
                df.ecdf = df.plot)
    )  
  }
  else{
    return(list(df.boxplot = df.boxplot,
                df.ecdf = df.plot))
  }
}



get_chisq.pvalue <- function(df){
  qq = seq(0,1, by = 0.1) # set up the break points
  hh = hist(df$value, qq, right = FALSE, plot = FALSE) # tabulate the counts
  cchi <- chisq.test(hh$counts)
  cchi$p.value
}

get_ks.pvalue <- function(df){
  qq = seq(0,1, by = 0.1) # set up the break points
  kks <- ks.test(df$value, punif)
  kks$p.value
}


# plot different TapGR laws
library(PtProcess)

M0 <- 2.5
Mc <- c(5, 6, 7, Inf)
b.v <- c(0.75, 1, 1.25)
b.Mc.grid <- expand.grid(b = b.v, Mc = Mc)

mm.b <- seq(2.5, 10, by = 0.1)
mob <- 10**((3/2)*mm.b + 9.1)

list.plot.tapgr <- lapply(1:nrow(b.Mc.grid), 
                          \(idx) data.frame(mag = mm.b,
                                          freq = log10(vapply(mob, \(x) 1000000*ptappareto(x, (2/3)*b.Mc.grid$b[idx], 
                                                                                           theta = 10**((3/2)*b.Mc.grid$Mc[idx] + 9.1), 
                                                                                           a = 10**((3/2)*M0 + 9.1), 
                                                                                           lower.tail = FALSE), 0)),
                                          Mc = b.Mc.grid$Mc[idx],
                                          b = b.Mc.grid$b[idx])) 

df.plot.tapgr <- Reduce(rbind, list.plot.tapgr)

write.csv(df.plot.tapgr[!is.infinite(df.plot.tapgr$freq),], 'df_tapgr.csv', row.names = FALSE)

pl.dim = 7
pdf('figures/TapGR_vs_GR.pdf', height = pl.dim/2, width = pl.dim)
ggplot(df.plot.tapgr[!is.infinite(df.plot.tapgr$freq),], aes(x = mag, y = freq, 
                                                             color = factor(Mc), 
                                                             linetype = factor(b))) + 
  geom_line() + 
  geom_vline(xintercept = Mc[1:3], alpha = 0.3) + 
  theme_classic() + 
  xlab('Magnitudes') + 
  ylab(~log[10](N[m])) + 
  labs(color = ~m[c], linetype = ~b) + 
  scale_x_continuous(breaks = 3:10) + 
  theme(legend.position = c(0.85, 0.55), 
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = 'white', colour = "grey", size = 1))
dev.off()


# forecast standard GR - observation TapGR 
# same b-value to see difference in the tail

# had to increase number of n.fore to 1000 (from 100) because results too varying if experiment repeated

b.Mc.grid2 <- b.Mc.grid[!is.infinite(b.Mc.grid$Mc),]


# create cluster for parallel computing

###############################
## FORECAST GR vs OBS TAP-GR ##
###############################

# models in this case overestimate the tails of the distribution.
cl <- makePSOCKcluster(detectCores() - 2)
registerDoParallel(cl)
# run the function in parallel
df.q.list.over <- lapply(1:nrow(b.Mc.grid2), \(idx) 
                    TAPGR.Score.ECDF.plot(b.fore = b.Mc.grid2$b[idx], 
                                          mc.fore = NULL, 
                                          type.fore = 'GR', 
                                          b.obs = b.Mc.grid2$b[idx], 
                                          mc.obs = b.Mc.grid2$Mc[idx], 
                                          type.obs = 'TapGR', 
                                          return.plot = FALSE, 
                                          N.o = 500, # expected number of observations
                                          n.fore = 1000, # number of catalogues per period
                                          n.rep = 1000 # number of independent periods
                    )$df.boxplot %>% 
                      mutate(b = b.Mc.grid2$b[idx], 
                             Mc = b.Mc.grid2$Mc[idx])
)
#save(df.q.list.over, file = 'df.q.list.over.Rds')
# stop the cluster
stopCluster(cl)

df.q.over <- Reduce(rbind, df.q.list.over)
df.q.over$panel.name <- paste0('b = ', df.q.over$b, ',Mc = ', df.q.over$Mc)
df.split.over <- split(df.q.over, df.q.over$panel.name)


# calculate prob gamma < 0.05
list.doubleside <- list() 
list.leftside <- list()
list.rightside <- list()
for(idx in 1:length(df.split.over)){
  df_scores <- df.split.over[[idx]]
  df_scores.split <- split(df_scores, df_scores$score)
  list.doubleside[[idx]] <- c(vapply(df_scores.split, \(x) mean(x$value <= 0.025 | x$value >= 0.975), 0),
                               b = unique(df_scores$b), Mc = unique(df_scores$Mc))
  list.leftside[[idx]] <- c(vapply(df_scores.split, \(x) mean(x$value <= 0.05), 0),
                              b = unique(df_scores$b), Mc = unique(df_scores$Mc))
  list.rightside[[idx]] <- c(vapply(df_scores.split, \(x) mean(x$value >= 0.95), 0),
                              b = unique(df_scores$b), Mc = unique(df_scores$Mc))
}

Reduce(rbind, list.doubleside)
Reduce(rbind, list.leftside)
Reduce(rbind, list.rightside)

# calculate p-value of ks uniformity test
list.res.over.ks <- list() 
for(idx in 1:length(df.split.over)){
  df_scores <- df.split.over[[idx]]
  idx.low <- df_scores$value <= 0.5
  df_scores$value[idx.low] <- df_scores$value[idx.low] + (rnorm(sum(idx.low), sd = 0.01))^2
  df_scores$value[!idx.low] <- df_scores$value[!idx.low] - (rexp(nrow(df_scores) - sum(idx.low), 10000))
  df_scores.split <- split(df_scores, df_scores$score)
  list.res.over.ks[[idx]] <- c(vapply(df_scores.split, get_ks.pvalue, 0),
                                  b = unique(df_scores$b), Mc = unique(df_scores$Mc))
}
ks.res.over <- Reduce(rbind, list.res.over.ks)
round(ks.res.over, 3)


library(latex2exp)

df.q.over$panel.name2 <- factor(df.q.over$panel.name, 
                               labels = c('b = 0.75,Mc = 5' = TeX('$b = 0.75, m_c = 5$'),
                                          'b = 0.75,Mc = 6' = TeX('$b = 0.75, m_c = 6$'),
                                          'b = 0.75,Mc = 7' = TeX('$b = 0.75, m_c = 7$'),
                                          'b = 1,Mc = 5' = TeX('$b = 1, m_c = 5$'),
                                          'b = 1,Mc = 6' = TeX('$b = 1, m_c = 6$'),
                                          'b = 1,Mc = 7' = TeX('$b = 1, m_c = 7$'),
                                          'b = 1.25,Mc = 5' = TeX('$b = 1.25, m_c = 5$'),
                                          'b = 1.25,Mc = 6' = TeX('$b = 1.25, m_c = 6$'),
                                          'b = 1.25,Mc = 7' = TeX('$b = 1.25, m_c = 7$')
                                          ))
  
library(RColorBrewer)
dark.pal <- brewer.pal(8, 'Dark2')

pdf('figures/TapGR_overest.pdf', height = pl.dim, width = pl.dim*1.5)
ggplot(df.q.over, aes(x = value, color = score, linetype = score)) + 
  stat_ecdf(geom = 'step', ) + 
  geom_abline(slope = 1, linetype = 2) + 
  labs(color = 'statistics', linetype = 'statistics') + 
  xlab(~gamma[p]) + 
  ylab('ECDF') + 
  facet_wrap(facets = ~panel.name2, labeller = label_parsed) + 
  theme_bw() + 
  theme(legend.box.background = element_rect(color="grey", linewidth = 2)) +
  scale_color_manual(values = dark.pal) 
dev.off()

write.csv(df.q.over, 'df_ecdf_over.csv', row.names = FALSE)


###############################
## FORECAST TAP-GR vs OBS GR ##
###############################

# models in this case underestimate the tails of the distribution.
cl <- makePSOCKcluster(detectCores() - 2)
registerDoParallel(cl)
# run the function in parallel
df.q.list.under <- lapply(1:nrow(b.Mc.grid2), \(idx) 
                         TAPGR.Score.ECDF.plot(b.fore = b.Mc.grid2$b[idx], 
                                               mc.fore = b.Mc.grid2$Mc[idx], 
                                               type.fore = 'TapGR', 
                                               b.obs = b.Mc.grid2$b[idx], 
                                               mc.obs = NULL, 
                                               type.obs = 'GR', 
                                               return.plot = FALSE, 
                                               N.o = 500, # expected number of observations
                                               n.fore = 1000, # number of catalogues per period
                                               n.rep = 1000 # number of independent periods
                         )$df.boxplot %>% 
                           mutate(b = b.Mc.grid2$b[idx], 
                                  Mc = b.Mc.grid2$Mc[idx])
)
# stop the cluster
stopCluster(cl)

df.q.under <- Reduce(rbind, df.q.list.under)
df.q.under$panel.name <- paste0('b = ', df.q.under$b, ',Mc = ', df.q.under$Mc)
df.split.under <- split(df.q.under, df.q.under$panel.name)
df.split.under.pert <- df.split.under

# calculate p-value of chi-squared uniformity test

# calculate prob gamma < 0.05
list.doubleside.under <- list() 
list.leftside.under <- list()
list.rightside.under <- list()
for(idx in 1:length(df.split.under)){
  df_scores <- df.split.under[[idx]]
  df_scores.split <- split(df_scores, df_scores$score)
  list.doubleside.under[[idx]] <- c(vapply(df_scores.split, \(x) mean(x$value <= 0.025 | x$value >= 0.975), 0),
                              b = unique(df_scores$b), Mc = unique(df_scores$Mc))
  list.leftside.under[[idx]] <- c(vapply(df_scores.split, \(x) mean(x$value <= 0.05), 0),
                            b = unique(df_scores$b), Mc = unique(df_scores$Mc))
  list.rightside.under[[idx]] <- c(vapply(df_scores.split, \(x) mean(x$value >= 0.95), 0),
                             b = unique(df_scores$b), Mc = unique(df_scores$Mc))
}

Reduce(rbind, list.doubleside.under)
Reduce(rbind, list.leftside.under)
Reduce(rbind, list.rightside.under)



# calculate p-value of ks uniformity test
list.res.under.ks <- list() 
for(idx in 1:length(df.split.under)){
  df_scores <- df.split.under[[idx]]
  idx.low <- df_scores$value <= 0.5
  df_scores$value[idx.low] <- df_scores$value[idx.low] + (rnorm(sum(idx.low), sd = 0.01))^2
  df_scores$value[!idx.low] <- df_scores$value[!idx.low] - (rnorm(length(idx.low) - sum(idx.low), sd = 0.01))^2
  df_scores.split <- split(df_scores, df_scores$score)
  list.res.under.ks[[idx]] <- c(vapply(df_scores.split, get_ks.pvalue, 0),
                               b = unique(df_scores$b), Mc = unique(df_scores$Mc))
}
ks.res.under <- Reduce(rbind, list.res.under.ks)
round(ks.res.under, 3)


library(latex2exp)

df.q.under$panel.name2 <- factor(df.q.under$panel.name, 
                                labels = c('b = 0.75,Mc = 5' = TeX('$b = 0.75, m_c = 5$'),
                                           'b = 0.75,Mc = 6' = TeX('$b = 0.75, m_c = 6$'),
                                           'b = 0.75,Mc = 7' = TeX('$b = 0.75, m_c = 7$'),
                                           'b = 1,Mc = 5' = TeX('$b = 1, m_c = 5$'),
                                           'b = 1,Mc = 6' = TeX('$b = 1, m_c = 6$'),
                                           'b = 1,Mc = 7' = TeX('$b = 1, m_c = 7$'),
                                           'b = 1.25,Mc = 5' = TeX('$b = 1.25, m_c = 5$'),
                                           'b = 1.25,Mc = 6' = TeX('$b = 1.25, m_c = 6$'),
                                           'b = 1.25,Mc = 7' = TeX('$b = 1.25, m_c = 7$')
                                ))

library(RColorBrewer)
dark.pal <- brewer.pal(8, 'Dark2')


write.csv(df.q.under, 'df_ecdf_under.csv', row.names = FALSE)

pdf('figures/TapGR_underest.pdf', height = pl.dim, width = pl.dim*1.5)
ggplot(df.q.under, aes(x = value, color = score, linetype = score)) + 
  stat_ecdf(geom = 'step', ) + 
  geom_abline(slope = 1, linetype = 2) + 
  labs(color = 'statistics', linetype = 'statistics') + 
  xlab(~gamma[p]) + 
  ylab('ECDF') + 
  facet_wrap(facets = ~panel.name2, labeller = label_parsed) + 
  theme_bw() + 
  theme(legend.box.background = element_rect(color="grey", linewidth = 2)) +
  scale_color_manual(values = dark.pal) 
dev.off()


#######################################
## STUDYING PROB OF GAMMA OUTSIDE CI ##
#######################################
# This is with TAP-GR law
TAPGR.Score.ECDF.plot_red <- function(b.fore, 
                                  mc.fore, 
                                  type.fore = 'TapGR',
                                  b.obs = 1, 
                                  mc.obs = 6,
                                  type.obs = 'TapGR',
                                  N.o = 100, # expected number of observations
                                  n.fore = 1000, # number of catalogues per period
                                  n.rep = 100, # number of independent periods
                                  M.max = 15,
                                  M0 = 2.5,
                                  return.plot = TRUE){
  
  # create forecast
  m.fore <- foreach(i = 1:n.fore, .export = name.fun.topar) %do% {
    n.p <- rpois(1, N.o)
    M0 <- M0
    if(type.fore == 'TapGR'){
      
      data.frame(m.ij = TapGRM(n = n.p, 
                               b = b.fore, 
                               corner_mag = mc.fore,
                               m_min = M0)) 
    } else if(type.fore == 'GR'){
      beta.fore <- b.fore*log(10)
      data.frame(m.ij = rexp(n.p, beta.fore) + M0)
    } else{
      stop('Uknown type.fore argument, please choose between `TapGR` and `GR`')
    }
    
  }
  # merge forecast
  m.u <- bind_rows(m.fore)
  N.u <- length(m.u$m.ij)
  # create breaks for binning
  m.breaks <- seq((M0 - 0.05), (M.max + 0.05), by = 0.1)
  
  # create bin counts of merged forecast
  h.u <- hist(m.u$m.ij, plot = FALSE, breaks = m.breaks)
  # this calculate the score of a sample from the merged forecast
  boot.score.sample <- foreach(i = 1:n.fore, .combine = rbind, .export = name.fun.topar) %do% {
    M0 <- M0
    m.j <- sample(m.u$m.ij, N.o, replace = TRUE)
    data.frame(chisq.score = Chi.squared.score(m.u$m.ij, m.j, M.max = 12),
               chisq.score_p1 = Chi.squared.score_plus1(m.u$m.ij, m.j, M.max = 12),
               csep.score = m.j.score(m.u$m.ij, m.j, M.max = 12, N.obs = N.o),
               multinom.score = Multinom.loglik.ratio(m.u$m.ij, m.j)    )
  }
  cat('Finished bootstrap score distribution \n')
  emp.score.distro <- lapply(1:ncol(boot.score.sample), \(x) 
                             ecdf(boot.score.sample[,x]))
  names(emp.score.distro) <- colnames(boot.score.sample)
  # observed score values
  obs.score <- foreach(i = 1:n.rep, .combine = rbind, .export = name.fun.topar) %do% {
    M0 <- M0
    if(type.obs == 'TapGR'){
      m.obs <- TapGRM(n = N.o, 
                      b = b.obs, 
                      corner_mag = mc.obs,
                      m_min = M0) 
      
    } else if(type.obs == 'GR'){
      beta.obs <- b.obs*log(10)
      m.obs <- rexp(N.o, beta.obs) + M0
      
    } else{
      stop('Uknown type.obs argument, please choose between `TapGR` and `GR`')
    }
    h.obs <- hist(m.obs, breaks = m.breaks, plot = FALSE)
    obs.chisq.score = Chi.squared.score(m.u$m.ij, m.obs, M.max = 12)
    obs.chisq.score_p1 = Chi.squared.score_plus1(m.u$m.ij, m.obs, M.max = 12)
    obs.csep.score = m.j.score(m.u$m.ij, m.obs, M.max = 12, N.obs = N.o)
    obs.multinom.score = Multinom.loglik.ratio(m.u$m.ij, m.obs)
    data.frame(Chisq = emp.score.distro$chisq.score(obs.chisq.score),
               Chisq.p1 = emp.score.distro$chisq.score_p1(obs.chisq.score_p1),
               Mtest = emp.score.distro$csep.score(obs.csep.score),
               Multinomial = emp.score.distro$multinom.score(obs.multinom.score)
               )
  }
  print('iam out')
  xx = seq(0, 1, 0.01)
  df.plot <- foreach(i = 1:ncol(obs.score), .combine = rbind, .export = name.fun.topar) %do% {
    data.frame(value = xx, 
               ECDF = sapply(xx, \(x) mean(obs.score[,i] <= x)), 
               score = colnames(obs.score)[i])
  }
  
  # increasing size of m.obs should help - also I can try it taking a multinomial with prob
  # given by the real theoretical frequencies.
  df.boxplot <- foreach(i = 1:ncol(obs.score), .combine = rbind,.export = name.fun.topar) %do% {
    data.frame(value = obs.score[,i], 
               score = colnames(obs.score)[i])
  }
  
  if(return.plot){
    return(list(boxplot = ggplot(df.boxplot, aes(x = score, y = value)) + 
                  geom_boxplot(),
                ECDF = ggplot(df.plot, aes(x = value, y = ECDF, color = score,
                                           linetype = score)) + geom_line() + 
                  geom_abline(slope = 1, intercept = 0, linetype = 2),
                df.boxplot = df.boxplot,
                df.ecdf = df.plot)
    )  
  }
  else{
    return(list(df.boxplot = df.boxplot,
                df.ecdf = df.plot))
  }
}

# Here we study the probability of taking a decision based on gamma being outside CI-alpha
# this will be considered a function of N number of observations and cutoff magnitude


N.v <- c(10, 25, 50, 100, 200, 350, 500)
M0.v <- c(2.5, 4, 5)

N.M0.grid <- expand.grid(N = N.v, M0 = M0.v)

st = Sys.time()
gammaseq.perN.under. <- lapply(1:nrow(N.M0.grid), \(idx) TAPGR.Score.ECDF.plot(b.fore = 1, 
                                                                               mc.fore = 6, 
                                                                               type.fore = 'TapGR', 
                                                                               b.obs = 1, 
                                                                               mc.obs = NULL, 
                                                                               type.obs = 'GR', 
                                                                               return.plot = FALSE, 
                                                                               N.o = N.M0.grid$N[idx], # expected number of observations
                                                                               n.fore = 1000, # number of catalogues per period
                                                                               n.rep = 1000, # number of independent periods
                                                                               M0 = N.M0.grid$M0[idx])$df.boxplot)

Sys.time() - st
# this takes 31.21199 mins
# but the probabilities are still wonky
# need to try n.rep = 10000 overnight

alpha_ <- 0.05
df.N.under<- foreach(idx = 1:nrow(N.M0.grid), .combine = rbind) %do% {
  gamma.seq <- gamma.seq.perN.under[[idx]]
  gamma.seq.split <- split(gamma.seq, gamma.seq$score)
  p.gamma.out.ci <- vapply(gamma.seq.split, \(x) mean( (x$value <= alpha_/2) | (x$value >= 1 - alpha_/2)), 0 )
  data.frame(t(p.gamma.out.ci), N = N.M0.grid$N[idx], M0 = N.M0.grid$M0[idx])
}

save(df.N.under, file = 'df.M0.N.underb1m6.Rds')
write.csv(df.N.under, 'df_M0_N_under.csv')

pl.score.under <- ggplot(df.N.under) + 
  geom_line(aes(N, Multinomial, color = factor(M0), linetype = 'Multinomial')) +
  geom_line(aes(N, Mtest, color = factor(M0), linetype = 'Mtest')) + 
  geom_line(aes(N, Chisq, color = factor(M0), linetype = 'Chisq')) + 
  geom_hline(yintercept = 0.05) +
  theme_classic() + 
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(fill = 'white', colour = "grey", linewidth = 1)) +
  ylab(~ beta) + 
  labs(colour = ~m[0], linetype = 'statistics')
  
pl.score.under




library(twosamples)





# This is with TAP-GR law
TAPGR.Standard.stat.test <- function(b.fore, 
                                     mc.fore, 
                                     type.fore = 'TapGR',
                                     b.obs = 1, 
                                     mc.obs = 6,
                                     type.obs = 'TapGR',
                                     N.o = 100, # expected number of observations
                                     n.fore = 1000, # number of catalogues per period
                                     n.rep = 100, # number of independent periods
                                     M.max = 15,
                                     M0 = 2.5){
  
  # create forecast
  m.fore <- foreach(i = 1:n.fore, .export = name.fun.topar) %dopar% {
    n.p <- rpois(1, N.o)
    M0 <- M0
    if(type.fore == 'TapGR'){
      
      data.frame(m.ij = TapGRM(n = n.p, 
                               b = b.fore, 
                               corner_mag = mc.fore,
                               m_min = M0)) 
    } else if(type.fore == 'GR'){
      beta.fore <- b.fore*log(10)
      data.frame(m.ij = rexp(n.p, beta.fore) + M0)
    } else{
      stop('Uknown type.fore argument, please choose between `TapGR` and `GR`')
    }
    
  }
  # merge forecast
  m.u <- bind_rows(m.fore)
  N.u <- length(m.u$m.ij)
  
  # observed score values
  obs.and.p.values <- foreach(i = 1:n.rep, .combine = rbind, .export = name.fun.topar) %dopar% {
    M0 <- M0
    if(type.obs == 'TapGR'){
      m.obs <- TapGRM(n = N.o, 
                      b = b.obs, 
                      corner_mag = mc.obs,
                      m_min = M0) 
      
    } else if(type.obs == 'GR'){
      beta.obs <- b.obs*log(10)
      m.obs <- rexp(N.o, beta.obs) + M0
      
    } else{
      stop('Uknown type.obs argument, please choose between `TapGR` and `GR`')
    }
    data.frame(ks = ks.test(m.u$m.ij, m.obs)$p.value,
               wil = wilcox.test(m.u$m.ij, m.obs)$p.value,
               ad = ad.test(m.obs, null = ecdf(m.u$m.ij))$p.value
    )
  }
  
  return(data.frame(ks = mean(obs.and.p.values$ks <= 0.05),
                    wil = mean(obs.and.p.values$wil <= 0.05),
                    ad = mean(obs.and.p.values$ad <= 0.05)))
}

st = Sys.time()
cl <- makePSOCKcluster(detectCores() - 2)
registerDoParallel(cl)
#debugonce(TAPGR.Standard.stat.test)
stand.test.df.under <- foreach(idx = 1:nrow(N.M0.grid), .combine = rbind) %do% {
  print(idx/nrow(N.M0.grid))
  df.p.value <- TAPGR.Standard.stat.test(b.fore = 1, 
                                         mc.fore = 6, 
                                         type.fore = 'TapGR', 
                                         b.obs = 1, 
                                         mc.obs = NULL, 
                                         type.obs = 'GR', 
                                         N.o = N.M0.grid$N[idx], # expected number of observations
                                         n.fore = 1000, # number of catalogues per period
                                         n.rep = 1000, # number of independent periods
                                         M0 = N.M0.grid$M0[idx])
  df.p.value$N <- N.M0.grid$N[idx]
  df.p.value$M0 <- N.M0.grid$M0[idx]
  df.p.value
}
stopCluster(cl)
Sys.time() - st # 2.32877 hours

save(stand.test.df.under, file = 'stats.M0.N.underb1m6.Rds')
load('stats.M0.N.underb1m6.Rds')

ggplot(stand.test.df.under, aes(x = N, linetype = factor(M0))) + 
  geom_line(aes(y = ad, color = 'ad')) + 
  geom_hline(yintercept = 0.05, alpha = 0.5) + 
  theme_bw()


total.df.under <- cbind(stand.test.df.under, 
                  Multinomial = df.N.under$Multinomial,
                  Mtest = df.N.under$Mtest)


ggplot(total.df.under, aes(x = N, linetype = factor(M0))) + 
  geom_line(aes(y = ks, color = 'ks')) + 
  geom_line(aes(y = Mtest, color = 'Mtest')) + 
  geom_hline(yintercept = 0.05, alpha = 0.5) + 
  theme_bw()


################################
## repeating in opposite case ##
################################
## over
N.v <- c(10, 25, 50, 100, 200, 350, 500)
M0.v <- c(2.5, 4, 5)
N.M0.grid <- expand.grid(N = N.v, M0 = M0.v)

st = Sys.time()
set.seed(2)
gamma.seq.perN.over <- lapply(1:nrow(N.M0.grid), \(idx) TAPGR.Score.ECDF.plot(b.fore = 1, 
                                                                         mc.fore = NULL, 
                                                                         type.fore = 'GR', 
                                                                         b.obs = 1, 
                                                                         mc.obs = 6, 
                                                                         type.obs = 'TapGR', 
                                                                         return.plot = FALSE, 
                                                                         N.o = N.M0.grid$N[idx], # expected number of observations
                                                                         n.fore = 1000, # number of catalogues per period
                                                                         n.rep = 1000, # number of independent periods
                                                                         M0 = N.M0.grid$M0[idx]
)$df.boxplot)
Sys.time() - st
# this takes 31.21199 mins
# but the probabilities are still wonky
# need to try n.rep = 10000 overnight

alpha_ <- 0.05
df.N.over <- foreach(idx = 1:nrow(N.M0.grid), .combine = rbind) %do% {
  gamma.seq <- gamma.seq.perN.over[[idx]]
  gamma.seq.split <- split(gamma.seq, gamma.seq$score)
  p.gamma.out.ci <- vapply(gamma.seq.split, \(x) mean( (x$value <= alpha_/2) | (x$value >= 1 - alpha_/2)), 0 )
  data.frame(t(p.gamma.out.ci), N = N.M0.grid$N[idx], M0 = N.M0.grid$M0[idx])
}


save(df.N.over, file = 'df.M0.N.overb1m6.Rds')
write.csv(df.N.over, 'df_M0_N_over.csv')

df.Nbb <- foreach(idx = 1:nrow(N.M0.grid), .combine = rbind) %do% {
  gamma.score <- gamma.seq.perN.over[[idx]]
  data.frame(gamma = gamma.score, N = N.M0.grid$N[idx], M0 = N.M0.grid$M0[idx])
}

write.csv(df.Nbb, 'df_Nbb.csv')
df.Nbb500 <- df.Nbb[df.Nbb$N == 500 & df.Nbb$M0 == 5,]
df.Nbb100 <- df.Nbb[df.Nbb$N == 100 & df.Nbb$M0 == 5,]
df.Nbb = rbind(df.Nbb100, df.Nbb500)
df.Nbb$panel.name = paste0('N = ', df.Nbb$N)

ggplot(df.Nbb, aes(gamma.score, gamma.value)) + 
  geom_boxplot() +
  facet_wrap(facets = ~ panel.name)

png('boxplot_strange.png', res = 300, width = 10, units = "in", height = 5)
ggplot(df.Nbb, aes(gamma.score, gamma.value)) + 
  geom_boxplot() +
  facet_wrap(facets = ~ panel.name) + 
  ylab(expression(gamma)) + 
  xlab('score')
dev.off()

save(df.N.over, file = 'df.M0.N.overb1m6.Rds')

pl.score.over <- ggplot(df.N.over) + 
  geom_line(aes(N, Multinomial, color = factor(M0), linetype = 'Multinomial')) +
  geom_line(aes(N, Mtest, color = factor(M0), linetype = 'Mtest')) + 
  geom_line(aes(N, Chisq, color = factor(M0), linetype = 'Chisq')) + 
  geom_hline(yintercept = 0.05) +
  theme_classic() + 
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(fill = 'white', colour = "grey", size = 1)) +
  ylab(~ beta) + 
  labs(colour = ~m[0], linetype = 'statistics')

pl.score.over


df.N.total <- rbind(cbind(df.N.over, panel.name = 'overestimation'),
                    cbind(df.N.under, panel.name = 'underestimation'))

png('figures/beta_scores.png',  res = 300, width = 10, units = "in", height = 5)
ggplot(df.N.total) + 
  geom_line(aes(N, Multinomial, color = factor(M0), linetype = 'Multinomial')) +
  geom_line(aes(N, Mtest, color = factor(M0), linetype = 'Mtest')) + 
  geom_line(aes(N, Chisq.p1, color = factor(M0), linetype = 'Chisq.p1')) + 
  geom_hline(yintercept = 0.05) +
  facet_wrap(facets = ~ panel.name) + 
  theme_bw() + 
  theme(legend.position = 'bottom',
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = 'white', colour = "grey", size = 1)) +
  ylab(~ beta) + 
  labs(color = ~m[0], linetype = '') 
dev.off()


library(kableExtra)

kbl(df.N.over, booktabs = TRUE) %>%
  save_kable('kable.table.pdf')



st = Sys.time()
cl <- makePSOCKcluster(detectCores() - 2)
registerDoParallel(cl)
#debugonce(TAPGR.Standard.stat.test)
stand.test.df.over <- foreach(idx = 1:nrow(N.M0.grid), .combine = rbind) %do% {
  print(idx/nrow(N.M0.grid))
  df.p.value <- TAPGR.Standard.stat.test(b.fore = 1, 
                                         mc.fore = NULL, 
                                         type.fore = 'GR', 
                                         b.obs = 1, 
                                         mc.obs = 6, 
                                         type.obs = 'TapGR',
                                         N.o = N.M0.grid$N[idx], # expected number of observations
                                         n.fore = 1000, # number of catalogues per period
                                         n.rep = 1000, # number of independent periods
                                         M0 = N.M0.grid$M0[idx])
  df.p.value$N <- N.M0.grid$N[idx]
  df.p.value$M0 <- N.M0.grid$M0[idx]
  df.p.value
}
stopCluster(cl)
Sys.time() - st # 2.32877 hours

save(stand.test.df.over, file = 'stats.M0.N.overb1m6.Rds')
load('stats.M0.N.overb1m6.Rds')

merged.over <- cbind(stand.test.df.over, 
                     Mtest = df.N.over$Mtest,
                     Multinomial = df.N.over$Multinomial)

write.csv(merged.over, file = 'df_M0_N_over_classic_stats.csv')

ggplot(merged.over, aes(x = N, linetype = factor(M0))) + 
  geom_line(aes(y = ks, color = 'KS')) + 
  geom_line(aes(y = wil, color = 'WX')) +
  geom_line(aes(y = ad, color = 'AD')) +
  geom_line(aes(y = Multinomial, color = 'Multinomial')) +
  geom_hline(yintercept = 0.05, alpha = 0.5) + 
  theme_bw()

merged.under <- cbind(stand.test.df.under, 
                      Mtest = df.N.under$Mtest, 
                      Multinomial = df.N.under$Multinomial)
write.csv(merged.under, file = 'df_M0_N_under_classic_stats.csv')

ggplot(merged.under, aes(x = N, linetype = factor(M0))) + 
  geom_line(aes(y = ks, color = 'KS')) + 
  geom_line(aes(y = wil, color = 'WX')) +
  geom_line(aes(y = ad, color = 'AD')) +
  geom_line(aes(y = Multinomial, color = 'Multinomial')) +
  geom_hline(yintercept = 0.05, alpha = 0.5) + 
  theme_bw()






merged.total <- rbind(cbind(merged.under, panel.name = 'underestimation'),
                      cbind(merged.over, panel.name = 'overestimation'))







pdf('figures/beta_stat_test.pdf', height = pl.dim/2, width = pl.dim)
ggplot(merged.total) + 
  geom_line(aes(x = N, y = ks, color = factor(M0), linetype = 'KS')) + 
  geom_line(aes(x = N, y = wil, color = factor(M0), linetype = 'WX')) +
  geom_line(aes(x = N, y = ad, color = factor(M0), linetype = 'AD')) +
  #geom_line(aes(x = N,y = Multinomial, color = factor(M0), linetype ='Multinomial')) +
  geom_hline(yintercept = 0.05, alpha = 0.5) +
  facet_wrap(facets = ~ panel.name) + 
  theme_bw() + 
  theme(legend.position = 'bottom',
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = 'white', colour = "grey", linewidth  = 1)) +
  ylab(~ beta^"'") + 
  labs(color = ~m[0], linetype = '') 
dev.off()

pdf('figures/beta_comparison.pdf', height = pl.dim/2, width = pl.dim)
ggplot(merged.total) + 
  geom_line(aes(x = N, y = ad, color = factor(M0), linetype = 'AD')) + 
  geom_line(aes(x = N,y = Multinomial, color = factor(M0), linetype ='Multinomial')) +
  geom_line(aes(x = N,y = Mtest, color = factor(M0), linetype ='Mtest')) +
  geom_hline(yintercept = 0.05, alpha = 0.5) +
  facet_wrap(facets = ~ panel.name) + 
  theme_bw() + 
  theme(legend.position = 'bottom',
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = 'white', colour = "grey", linewidth  = 1)) +
  ylab(~ beta) + 
  labs(color = ~m[0], linetype = '') 
dev.off()


merged.under.for.tile <- rbind(data.frame(beta = stand.test.df.under$ks,
                                          panel.name = 'KS',
                                          N = stand.test.df.under$N,
                                          M0 = stand.test.df.under$M0),
                               data.frame(beta = stand.test.df.under$wil,
                                          panel.name = 'WX',
                                          N = stand.test.df.under$N,
                                          M0 = stand.test.df.under$M0),
                               data.frame(beta = stand.test.df.under$ad,
                                          panel.name = 'AD',
                                          N = stand.test.df.under$N,
                                          M0 = stand.test.df.under$M0),
                               data.frame(beta = df.N.under$Brier,
                                          panel.name = 'Brier',
                                          N = df.N.under$N,
                                          M0 = df.N.under$M0),
                               data.frame(beta = df.N.under$Chisq,
                                          panel.name = 'Chisq',
                                          N = df.N.under$N,
                                          M0 = df.N.under$M0),
                               data.frame(beta = df.N.under$Chisq.p1,
                                          panel.name = 'Chisq+1',
                                          N = df.N.under$N,
                                          M0 = df.N.under$M0),
                               data.frame(beta = df.N.under$Hellinger,
                                          panel.name = 'Hellinger',
                                          N = df.N.under$N,
                                          M0 = df.N.under$M0),
                               data.frame(beta = df.N.under$Mtest,
                                          panel.name = 'M-test',
                                          N = df.N.under$N,
                                          M0 = df.N.under$M0),
                               data.frame(beta = df.N.under$Mtest.un,
                                          panel.name = 'UM-test',
                                          N = df.N.under$N,
                                          M0 = df.N.under$M0),
                               data.frame(beta = df.N.under$Multinomial,
                                          panel.name = 'MLL',
                                          N = df.N.under$N,
                                          M0 = df.N.under$M0))

pdf('figures/beta_total_under.pdf', height = pl.dim*1.5, width = pl.dim*1.2)
ggplot(merged.under.for.tile, aes(x = factor(N), y = factor(M0), fill = beta)) + 
  geom_tile() + 
  geom_text(aes(label = beta)) + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  scale_fill_gradientn(colours = c('red','yellow', 'green')) + 
  theme_bw() + 
  xlab('N') + 
  ylab(~m[0]) + 
  labs(fill = ~beta)
dev.off()



merged.over.for.tile <- rbind(data.frame(beta = stand.test.df.over$ks,
                                          panel.name = 'KS',
                                          N = stand.test.df.over$N,
                                          M0 = stand.test.df.over$M0),
                               data.frame(beta = stand.test.df.over$wil,
                                          panel.name = 'WX',
                                          N = stand.test.df.over$N,
                                          M0 = stand.test.df.over$M0),
                               data.frame(beta = stand.test.df.over$ad,
                                          panel.name = 'AD',
                                          N = stand.test.df.over$N,
                                          M0 = stand.test.df.over$M0),
                               data.frame(beta = df.N.over$Brier,
                                          panel.name = 'Brier',
                                          N = df.N.over$N,
                                          M0 = df.N.over$M0),
                               data.frame(beta = df.N.over$Chisq,
                                          panel.name = 'Chisq',
                                          N = df.N.over$N,
                                          M0 = df.N.over$M0),
                              data.frame(beta = df.N.over$Chisq.p1,
                                         panel.name = 'Chisq+1',
                                         N = df.N.over$N,
                                         M0 = df.N.over$M0),
                               data.frame(beta = df.N.over$Hellinger,
                                          panel.name = 'Hellinger',
                                          N = df.N.over$N,
                                          M0 = df.N.over$M0),
                               data.frame(beta = df.N.over$Mtest,
                                          panel.name = 'M-test',
                                          N = df.N.over$N,
                                          M0 = df.N.over$M0),
                               data.frame(beta = df.N.over$Mtest.un,
                                          panel.name = 'UM-test',
                                          N = df.N.over$N,
                                          M0 = df.N.over$M0),
                               data.frame(beta = df.N.over$Multinomial,
                                          panel.name = 'MLL',
                                          N = df.N.over$N,
                                          M0 = df.N.over$M0))

pdf('figures/beta_total_over.pdf', height = pl.dim*1.5, width = pl.dim*1.2)
ggplot(merged.over.for.tile, aes(x = factor(N), y = factor(M0), fill = beta)) + 
  geom_tile() + 
  geom_text(aes(label = beta)) + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  scale_fill_gradientn(colours = c('red','yellow', 'green')) + 
  theme_bw() + 
  xlab('N') + 
  ylab(~m[0]) + 
  labs(fill = ~beta)
dev.off()


# small experiment on sample catalogues vs union forecast
samp.fun <- function(x, ns, repl = TRUE){
  if(length(x) == 1){
    return(rep(x, ns))
  } else{
    sample(x, ns, replace = repl)
  }
}


produce.tail.estimates <- function(v.tail = 10, n.cat = 1000, n.obs = 100, size.negbinom = 1){
  n.samp <- rnbinom(n.cat, mu = n.obs, size = size.negbinom) + 1
  res.matrix <- matrix(NA, ncol = 4, nrow = 1000)
  for(i in 1:nrow(res.matrix)){
    l.samp <- lapply(1:n.cat, \(x) rexp(n.samp[x]))
    l.u <- unlist(l.samp)
    l.samp.resamp.cat <- lapply(l.samp, function(x) samp.fun(x, n.obs))
    l.samp.resamp.u <- lapply(l.samp, function(x) samp.fun(l.u, n.obs))
    
    samp.resamp.cat <- unlist(l.samp.resamp.cat)
    samp.resamp.u <- unlist(l.samp.resamp.u)
    
    
    res.matrix[i,] = c(true = pexp(v.tail, lower.tail = FALSE),
                       from.l.u = mean(l.u >= v.tail),
                       from.samp.cat = mean(samp.resamp.cat >= v.tail),
                       from.samp.u = mean(samp.resamp.u >= v.tail) )
    
  }
 return(res.matrix)
}


produce.tail.estimates2 <- function(v.tail = 10, n.cat = 1000, n.obs = 100, size.negbinom = 1){
  n.samp <- rnbinom(n.cat, mu = n.obs, size = size.negbinom) + 1
  res.matrix <- matrix(NA, ncol = 4, nrow = 1000)
  for(i in 1:nrow(res.matrix)){
    l.samp <- lapply(1:n.cat, \(x) rexp(n.samp[x]))
    l.u <- unlist(l.samp)
    l.samp.resamp.cat <- lapply(l.samp, function(x) samp.fun(x, n.obs))
    l.samp.resamp.u <- lapply(l.samp, function(x) samp.fun(l.u, n.obs))
    
    est.resamp.cat <- unlist(lapply(l.samp.resamp.cat, \(x) mean(x >= v.tail)))
    est.resamp.u <- unlist(lapply(l.samp.resamp.u, \(x) mean(x >= v.tail)))
    
    
    res.matrix[i,] = c(true = pexp(v.tail, lower.tail = FALSE),
                       from.l.u = mean(l.u >= v.tail),
                       from.samp.cat = mean(est.resamp.cat),
                       from.samp.u = mean(est.resamp.u))
    
  }
  return(res.matrix)
}


calculate.summaries <- function(res.matrix){
  tt <- rbind(mean = apply(res.matrix, 2, mean),
              sd = apply(res.matrix, 2, sd),
              lower.ci = apply(res.matrix, 2, \(x) quantile(x, 0.025)),
              median = apply(res.matrix, 2, \(x) quantile(x, 0.5)),
              higher.ci = apply(res.matrix, 2, \(x) quantile(x, 0.975)))
  
  colnames(tt) <- c('True', 'estimate L_u', 'samp catalogues', 'samp L_u')
  data.frame(t(tt))
}
res.matrix <- produce.tail.estimates(v.tail = 11)
calculate.summaries(res.matrix)


ggplot(mapping = aes(x = 1:nrow(res.matrix))) + 
  geom_point(aes(y = res.matrix[,1])) + 
  geom_point(aes(y = res.matrix[,2], color = 'L_u')) + 
  geom_point(aes(y = res.matrix[,3], color = 'samples from cat')) + 
  geom_point(aes(y = res.matrix[,4], color = 'samples from L_u'))


v.tails <- seq(6, 16, length.out = 10)

list.of.estimates <- lapply(v.tails, 
                            function(vt) produce.tail.estimates2(v.tail = vt))
list.of.summaries <- lapply(list.of.estimates, 
                            function(rmat) calculate.summaries(res.matrix = rmat))

list.df.est.lu <- lapply(1:length(list.of.summaries), function(sum.idx) 
  data.frame(v.tail = v.tails[sum.idx],
             true = list.of.summaries[[sum.idx]]$mean[1],
             mean = list.of.summaries[[sum.idx]]$mean[2],
             lower.ci = list.of.summaries[[sum.idx]]$lower.ci[2],
             median = list.of.summaries[[sum.idx]]$median[2],
             upper.ci = list.of.summaries[[sum.idx]]$higher.ci[2],
             method = 'estimate from l_u'))
df.est.lu <- Reduce(rbind, list.df.est.lu)

list.df.samp.cat <- lapply(1:length(list.of.summaries), function(sum.idx) 
  data.frame(v.tail = v.tails[sum.idx],
             true = list.of.summaries[[sum.idx]]$mean[1],
             mean = list.of.summaries[[sum.idx]]$mean[3],
             lower.ci = list.of.summaries[[sum.idx]]$lower.ci[3],
             median = list.of.summaries[[sum.idx]]$median[3],
             upper.ci = list.of.summaries[[sum.idx]]$higher.ci[3],
             method = 'samp from catalogues'))
df.samp.cat <- Reduce(rbind, list.df.samp.cat)

list.df.samp.lu <- lapply(1:length(list.of.summaries), function(sum.idx) 
  data.frame(v.tail = v.tails[sum.idx],
             true = list.of.summaries[[sum.idx]]$mean[1],
             mean = list.of.summaries[[sum.idx]]$mean[4],
             lower.ci = list.of.summaries[[sum.idx]]$lower.ci[4],
             median = list.of.summaries[[sum.idx]]$median[4],
             upper.ci = list.of.summaries[[sum.idx]]$higher.ci[4],
             method = 'samp from l_u'))
df.samp.lu <-  Reduce(rbind, list.df.samp.lu)

df.to.plot <- rbind(df.est.lu,
                    df.samp.cat,
                    df.samp.lu)


ggplot(df.to.plot, aes(x = v.tail, y = mean, color = method)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill = method), alpha = 0.2) + 
  geom_line(aes(y = true), color = 'black')

ggplot(df.to.plot[df.to.plot$v.tail > 11,], aes(x = v.tail, y = mean, color = method)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill = method), alpha = 0.2) + 
  geom_line(aes(y = true), color = 'black')





list.of.estimates




