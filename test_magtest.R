source('magtest_fun.R')
library(inlabru)
library(foreach)
library(dplyr)

# set parameters
M0 = 2.5
b.v <- 1

# This is with standard GR law
Score.ECDF.plot <- function(b.fore, b.obs = 1, 
                            N.o = 100, 
                            n.fore = 1000, 
                            n.rep = 100,
                            M.max = 15,
                            return.plot = TRUE){
  beta.obs <- b.obs*log(10)
  beta.fore <- b.fore*log(10)
  # create forecast
  m.fore <- foreach(i = 1:n.fore) %do% {
    #n.p <- rpois(1, N.o)
    data.frame(m.ij = rexp(N.o, beta.fore) + M0) 
  }
  # merge forecast
  m.u <- bind_rows(m.fore)
  N.u <- length(m.u$m.ij)
  # create breaks for binning
  m.breaks <- seq((2.5 - 0.1), (M.max + 0.1 ), by = 0.1)
  
  # create bin counts of merged forecast
  h.u <- hist(m.u$m.ij, plot = FALSE, breaks = m.breaks)
  # this calculate the score of a sample from the merged forecast
  boot.score.sample <- foreach(i = 1:n.fore, .combine = rbind) %do% {
    m.j <- sample(m.u$m.ij, N.o, replace = TRUE)
    h.j <- hist(m.j, breaks = m.breaks, plot = FALSE)
    N.j <- length(m.j)
    data.frame(chisq.score = Chi.squared.score(m.u$m.ij, m.j, M.max = 12),
               chisq.score.p1 = Chi.squared.score_plus1(m.u$m.ij, m.j, M.max = 12),
               chisq.score.mod = Chi.squared.score_mod(m.u$m.ij, m.j, M.max = 12), #chisq.test(aa, p = c.hat.pos)$statistic,
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
  obs.score <- foreach(i = 1:n.rep, .combine = rbind) %do% {
    m.obs <- rexp(N.o, beta.obs) + M0
    h.obs <- hist(m.obs, breaks = m.breaks, plot = FALSE)
    N.obs <- length(m.obs)
   
    obs.chisq.score = Chi.squared.score(m.u$m.ij, m.obs, M.max = 12)
    obs.chisq.score.p1 = Chi.squared.score_plus1(m.u$m.ij, m.obs, M.max = 12)
    obs.chisq.score.mod = Chi.squared.score_mod(m.u$m.ij, m.obs, M.max = 12)
    obs.csep.score = m.j.score(m.u$m.ij, m.obs, M.max = 12, N.obs = N.o)
    obs.csep.score.un = m.j.score(m.u$m.ij, m.obs, M.max = 12, N.obs = N.o, type = 'alt')
    obs.multinom.score = Multinom.loglik.ratio(m.u$m.ij, m.obs)
    obs.hellinger = Hellinger.stat(m.u$m.ij, m.obs, m.breaks)
    obs.brier = Brier.score(m.u$m.ij, m.obs, m.breaks)
    
    data.frame(Chisq = emp.score.distro$chisq.score(obs.chisq.score),
               Chisq.p1 = emp.score.distro$chisq.score.p1(obs.chisq.score.p1),
               Chisq.mod = emp.score.distro$chisq.score.mod(obs.chisq.score.mod),
               `UM-test` = emp.score.distro$csep.score(obs.csep.score),
               `M-test` = emp.score.distro$csep.score.un(obs.csep.score.un),
               MLL = emp.score.distro$multinom.score(obs.multinom.score),
               Hellinger = emp.score.distro$hellinger(obs.hellinger),
               Brier = emp.score.distro$brier(obs.brier)
    )
  }
  
  xx = seq(0, 1, 0.01)
  df.plot <- foreach(i = 1:ncol(obs.score), .combine = rbind) %do% {
    data.frame(value = xx, 
               ECDF = sapply(xx, \(x) mean(obs.score[,i] <= x)), 
               score = colnames(obs.score)[i])
  }
  
  # increasing size of m.obs should help - also I can try it taking a multinomial with prob
  # given by the real theoretical frequencies.
  df.boxplot <- foreach(i = 1:ncol(obs.score), .combine = rbind) %do% {
    data.frame(value = obs.score[,i], 
               score = colnames(obs.score)[i])
  }
  
  if(return.plot){
    return(list(
      boxplot = ggplot(df.boxplot, aes(x = score, y = value)) + 
                  geom_boxplot(),
                
      ECDF = ggplot(df.plot, aes(x = value, y = ECDF, color = score, linetype = score)) + 
        geom_line() + 
        geom_abline(slope = 1, intercept = 0, linetype = 2),
      dfs = list(df.boxplot = df.boxplot,
                 df.ecdf = df.plot) )
    )  
  }
  else{
    return(list(df.boxplot = df.boxplot,
                df.ecdf = df.plot))
  }
}






# 
# gamma.seq <- rep(NA, 10)
# for(fore.period in 1:10){
#   m.u <- rexp(10000, 2) + M0
#   n.obs <- 5
#   omeg <- rexp(n.obs, 2) + M0
#   mlr <- rep(NA, 1000)
#   for(i in 1:1000){
#     mj.samp <- sample(m.u, n.obs, replace = TRUE)
#     mlr[i] <- Multinom.loglik.ratio(m.u_ = m.u, m.j_ = mj.samp)
#   }
#   mlr.obs <- Multinom.loglik.ratio(m.u_ = m.u, m.j_ = omeg)
#   gamma.seq[fore.period] <- mean(mlr <= mlr.obs)
# }
# 
# ecdf.gamma <- ecdf(gamma.seq)
# plot((1:10)/10 , ecdf.gamma((1:10)/10))
# abline(a = 0, b = 1)
# 
# 
# ecdf.gamma <- ecdf(runif(15))
# plot((1:10)/10 , ecdf.gamma((1:10)/10))
# abline(a = 0, b = 1)

# 1000 boot samples from Lambda_U to build the score distribution under the null
# 1000 observed catalogue of 100 events from GR b = 1 

## SEED = 2 AMD 4 ARE THE BEST
set.seed(1)
ECDF.plot.b1 <- Score.ECDF.plot(b.fore = 1, N.o = 100, n.rep = 1000)

write.csv(ECDF.plot.b1$dfs$df.ecdf, 'ecdf_body_1_s1.csv', row.names = FALSE)
write.csv(ECDF.plot.b1$dfs$df.boxplot, 'boxplot_body_1_s1.csv', row.names = FALSE)

set.seed(2)
ECDF.plot.b1 <- Score.ECDF.plot(b.fore = 1, N.o = 100, n.rep = 1000)
write.csv(ECDF.plot.b1$dfs$df.ecdf, 'ecdf_body_1_s2.csv', row.names = FALSE)
write.csv(ECDF.plot.b1$dfs$df.boxplot, 'boxplot_body_1_s2.csv', row.names = FALSE)

set.seed(3)
ECDF.plot.b1 <- Score.ECDF.plot(b.fore = 1, N.o = 100, n.rep = 1000)
write.csv(ECDF.plot.b1$dfs$df.ecdf, 'ecdf_body_1_s3.csv', row.names = FALSE)
write.csv(ECDF.plot.b1$dfs$df.boxplot, 'boxplot_body_1_s3.csv', row.names = FALSE)

set.seed(4)
ECDF.plot.b1 <- Score.ECDF.plot(b.fore = 1, N.o = 100, n.rep = 1000)
write.csv(ECDF.plot.b1$dfs$df.ecdf, 'ecdf_body_1_s4.csv', row.names = FALSE)
write.csv(ECDF.plot.b1$dfs$df.boxplot, 'boxplot_body_1_s4.csv', row.names = FALSE)


ECDF.plot.b1$dfs$df.boxplot

library(RColorBrewer)
dark.pal <- brewer.pal(8, 'Dark2')

ECDF.plot.b1$ECDF.2 <- ECDF.plot.b1$ECDF +
  labs(color = 'statistics',
       linetype = 'statistics',
       title = '') +
  xlab(~gamma[p]) +
  theme_bw() + 
  theme(legend.position.inside  = c(0.3, 0.8),
        legend.box.background = element_rect(color="grey", linewidth = 2),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) + 
  scale_color_manual(values = dark.pal) + 
  annotate('text', x = 0.85, y = 0.15, label = '(a)')
ECDF.plot.b1$boxplot.2 <- ECDF.plot.b1$boxplot + theme_bw() + 
  labs(title = '') + 
  xlab('statistics') + 
  ylab(~gamma[p]) + 
  annotate('text', x = 5.5, y = 0.15, label = '(b)') 
pdf('figures/check_corr_b1.pdf', width = pl.dim*2, height = pl.dim)
multiplot(ECDF.plot.b1$ECDF.2, 
          ECDF.plot.b1$boxplot.2, cols = 2)
dev.off()

nn = 100
ECDF.plot.b1.1 <- Score.ECDF.plot(b.fore = 1.1, return.plot = FALSE, N.o = nn, n.rep = 1000)
ECDF.plot.b1.5 <- Score.ECDF.plot(b.fore = 1.4, return.plot = FALSE, N.o = nn, n.rep = 1000)
ECDF.plot.b09 <- Score.ECDF.plot(b.fore = 0.9, return.plot = FALSE, N.o = nn, n.rep = 1000)
ECDF.plot.b075 <- Score.ECDF.plot(b.fore = 0.6, return.plot = FALSE, N.o = nn, n.rep = 1000)

write.csv(ECDF.plot.b1.1$df.ecdf, 'ecdf_body_11.csv', row.names = FALSE)
write.csv(ECDF.plot.b1.1$df.boxplot, 'boxplot_body_11.csv', row.names = FALSE)

write.csv(ECDF.plot.b1.5$df.ecdf, 'ecdf_body_15.csv', row.names = FALSE)
write.csv(ECDF.plot.b1.5$df.boxplot, 'boxplot_body_15.csv', row.names = FALSE)

write.csv(ECDF.plot.b09$df.ecdf, 'ecdf_body_09.csv', row.names = FALSE)
write.csv(ECDF.plot.b09$df.boxplot, 'boxplot_body_09.csv', row.names = FALSE)

write.csv(ECDF.plot.b075$df.ecdf, 'ecdf_body_075.csv', row.names = FALSE)
write.csv(ECDF.plot.b075$df.boxplot, 'boxplot_body_075.csv', row.names = FALSE)






ecdf.df.toplot <- rbind(ECDF.plot.b1.1$df.ecdf %>% mutate(panel.name = 'b[f] == 1.1'),
                        ECDF.plot.b1.5$df.ecdf %>% mutate(panel.name = 'b[f] == 1.5'),
                        ECDF.plot.b09$df.ecdf %>% mutate(panel.name = 'b[f] == 0.9'),
                        ECDF.plot.b075$df.ecdf %>% mutate(panel.name = 'b[f] == 0.75'))

pdf('figures/check_sens.pdf', width = pl.dim*2, height = pl.dim)
ggplot(ecdf.df.toplot, aes(value, ECDF - value, color = score, linetype = score)) + 
  geom_line() + 
  labs(color = 'statistics',
       linetype = 'statistics') + 
  xlab(~gamma[p]) +
  geom_abline(slope = 0, linetype = 2) + 
  facet_wrap(facets = ~panel.name, labeller = label_parsed) + 
  theme_bw() + 
  theme(legend.box.background = element_rect(color="grey", linewidth = 3)) +
  scale_color_manual(values = dark.pal) 
dev.off()

boxpl.df.toplot <- rbind(ECDF.plot.b1.1$df.boxplot %>% mutate(panel.name = 'b[f] == 1.1'),
                        ECDF.plot.b1.5$df.boxplot %>% mutate(panel.name = 'b[f] == 1.5'),
                        ECDF.plot.b09$df.boxplot %>% mutate(panel.name = 'b[f] == 0.9'),
                        ECDF.plot.b075$df.boxplot %>% mutate(panel.name = 'b[f] == 0.75'))

pdf('figures/check_sens.boxplot.pdf', width = pl.dim*2, height = pl.dim)
ggplot(boxpl.df.toplot, aes(x = score, y = value)) + 
  geom_boxplot() + 
  labs(color = 'statistics') + 
  ylab(~gamma[p]) +
  facet_wrap(facets = ~panel.name, labeller = label_parsed) + 
  theme_bw() 
dev.off()


## code percentage of regection as function of beta.v multiplier
perc.rej <- function(b.fore,  b.obs = 1, thr, n.rep = 1000, N.o = 100, M.max = 20){
  beta.fore <- b.fore*log(10)
  beta.v <- b.obs*log(10)
  # create forecast
  m.fore <- foreach(i = 1:1000) %do% {
    n.p <- rpois(1, N.o)
    data.frame(m.ij = rexp(n.p, beta.fore) + M0) 
  }
  
  # merge forecast
  m.u <- bind_rows(m.fore)
  N.u <- length(m.u$m.ij)
  # create breaks for binning
  m.breaks <- seq((2.5 - 0.1), (M.max + 0.1 ), by = 0.1)
  
  # create bin counts of merged forecast
  h.u <- hist(m.u$m.ij, plot = FALSE, breaks = m.breaks)
  
  boot.score.sample <- foreach(i = 1:1000, .combine = rbind) %do% {
    m.j <- sample(m.u$m.ij, N.o, replace = FALSE)
    h.j <- hist(m.j, breaks = m.breaks, plot = FALSE)
    N.j <- length(m.j)
    #c.hat <- (h.u$counts + h.j$counts)/(N.j + N.u)
    c.hat <- (h.u$counts)/(N.u)
    c.hat.pos <- c.hat[c.hat > 0 | h.j$counts > 0]
    counts.j <- h.j$counts[c.hat > 0 | h.j$counts > 0]
    aa <- counts.j#rmultinom(1, 100, prob = c.hat.pos)#
    data.frame(chisq.score = chisq.test(aa, p = c.hat.pos)$statistic,
               csep.score = m.j.score(m.u$m.ij, m.j, M.max = 12, N.obs = N.o),
               csep.score.un = m.j.score(m.u$m.ij, m.j, M.max = 12, N.obs = N.o, type = 'alt'),
               multinom.score = Multinom.loglik.ratio(m.u$m.ij, m.j),
               hellinger = Hellinger.stat(m.u$m.ij, m.j, m.breaks),
               brier = Brier.score(m.u$m.ij, m.j, m.breaks)
    )
  }
  
  emp.score.distro <- lapply(1:ncol(boot.score.sample), \(x) 
                             ecdf(boot.score.sample[,x]))
  names(emp.score.distro) <- colnames(boot.score.sample)
  cat('Finished to obtain ECDFs of the score under H0', '\n')
  obs.score.samples <- foreach(i = 1:n.rep, .combine = rbind) %do% {
    m.obs <- rexp(N.o, beta.v) + M0
    m.obs <- m.obs[m.obs <= M.max]
    h.obs <- hist(m.obs, breaks = m.breaks, plot = FALSE)
    N.obs <- length(m.obs)
    #c.hat <- (h.u$counts + h.obs$counts)/(N.obs + N.u)
    c.hat <- (h.u$counts)/(N.u)
    c.hat.pos <- c.hat[c.hat > 0 | h.obs$counts > 0]
    counts.obs <- h.obs$counts[c.hat > 0 | h.obs$counts > 0]
    obs.chisq.score = chisq.test(counts.obs, p = c.hat.pos)$statistic
    obs.csep.score = m.j.score(m.u$m.ij, m.obs, M.max = 100, N.obs = N.o)
    obs.csep.score.un = m.j.score(m.u$m.ij, m.obs, M.max = 100, N.obs = N.o, type = 'alt')
    obs.multinom.score = Multinom.loglik.ratio(m.u$m.ij, m.obs)
    obs.hellinger = Hellinger.stat(m.u$m.ij, m.obs, m.breaks)
    obs.brier = Brier.score(m.u$m.ij, m.obs, m.breaks)
    data.frame(chisq.score = emp.score.distro$chisq.score(obs.chisq.score),
               csep.score = emp.score.distro$csep.score(obs.csep.score),
               csep.score.un = emp.score.distro$csep.score.un(obs.csep.score.un),
               multinom.score = emp.score.distro$multinom.score(obs.multinom.score),
               hellinger = emp.score.distro$hellinger(obs.hellinger),
               brier = emp.score.distro$brier(obs.brier)
    )
  }
  
  res.tib <- foreach(i = 1:ncol(obs.score.samples), .combine = rbind) %do% {
    data.frame(value = obs.score.samples[,i], 
               score = colnames(obs.score.samples)[i])
  } %>% 
    mutate(rej = ( value <= thr | value >= 1 - thr ) ) %>%
    group_by(score) %>%
    summarise(perc.rej = mean(rej))
  as.data.frame(res.tib)
}


# rejection if score <= thr or score >= 1 - thr 
# prob of rejection for different b-value
b.vec = seq(0.5, 1.5, 0.1)

res.rej <- foreach(i = 1:length(b.vec), .combine = rbind) %do% {
  cat('b = ', b.vec[i], '\n')
  perc.rej(b.fore = b.vec[i],
           thr = 0.025,
           n.rep = 5000)
  
} %>% mutate(bvalue = rep(b.vec, each = 6))

ggplot(res.rej, aes(x = bvalue, y = perc.rej, color = score)) + geom_line() + 
  geom_hline(yintercept = 0.025, linetype = 2)























