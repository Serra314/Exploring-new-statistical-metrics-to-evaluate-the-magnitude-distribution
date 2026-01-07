library(foreach)
library(ggplot2)
library(dplyr)
library(parallel)
library(goftest)

m.obs.score <- function(m.u_, m.obs_, M.max = 8, M0 = 2.5){
  m.u_ <- m.u_[m.u_ <= M.max]
  m.obs_ <- m.obs_[m.obs_ <= M.max]
  N.obs <- length(m.obs_)
  N.u <- length(m.u_)
  m.breaks <- seq((M0 - 0.1), ( M.max + 0.1 ), by = 0.1)
  h.obs <- hist(m.obs_, plot = FALSE, breaks = m.breaks)
  h.u <- hist(m.u_, plot = FALSE, breaks = m.breaks)
  
  #sum((log(h.u$density + 1) - log(h.obs$density + 1))^2)
  sum((log((N.obs/N.u)*h.u$counts + 1) - log(h.obs$counts + 1))^2)
} 

m.j.score <- function(m.u_, m.j_, M.max = 8, M0 = 2.5, N.obs, type = 'def'){
  m.u_ <- m.u_[m.u_ <= M.max]
  m.j_ <- m.j_[m.j_ <= M.max]
  N.j <- length(m.j_)
  N.u <- length(m.u_)
  m.breaks <- seq((M0 - 0.1), ( M.max + 0.1 ), by = 0.1)
  h.j <- hist(m.j_, plot = FALSE, breaks = m.breaks)
  h.u <- hist(m.u_, plot = FALSE, breaks = m.breaks)
  if(type == 'def'){
    return(sum((log((N.obs/N.u)*h.u$counts + 1) - log((N.obs/N.j)*h.j$counts + 1))^2))  
  }
  else{
    return(return(sum((log((1/N.u)*h.u$counts + 1) - log((1/N.j)*h.j$counts + 1))^2)))
  }
}


Dj.samples <- function(n.fore = 1000, lambda_N = 100, beta.v = log(10), M0 = 2.5, N.obs_ = 100){
  m.fore <- foreach(i = 1:n.fore) %do% {
    n.p <- rpois(1, lambda_N)
    data.frame(m.ij = rexp(n.p, beta.v) + M0)
  }
  m.u <- bind_rows(m.fore)
  list(Dj = unlist(mclapply(m.fore, \(x) m.j.score(m.u$m.ij, x$m.ij, N.obs = N.obs_), mc.cores = 5)),
       mu = m.u)
}

#Chi.squared.score2 <- function(m.u_, m.j_, M.max = 15){
#  m.u_ <- m.u_[m.u_ <= M.max]
#  m.j_ <- m.j_[m.j_ <= M.max]
#  N.j <- length(m.j_)
#  N.u <- length(m.u_)
#  m.breaks <- seq((M0 - 0.1), (M.max + 0.1 ), by = 0.1)
#  h.j <- hist(m.j_, plot = FALSE, breaks = m.breaks)
#  h.u <- hist(m.u_, plot = FALSE, breaks = m.breaks)
#  c.hat <- (h.u$counts + h.j$counts)/(N.j + N.u)
#  idx.nzero <- c.hat > 0 ; nzero <- sum(idx.nzero)
  #chisq.test(cbind(h.j$counts, h.u$counts)[h.j$counts != 0 & h.u$counts != 0, ])$statistic
  # stat = sum(
  #   ( (( h.u$counts[idx.nzero] - c.hat[idx.nzero]*N.u)^2 )/ (N.u*c.hat[idx.nzero]) ) + 
  #     ( (( h.j$counts[idx.nzero] - c.hat[idx.nzero]*N.j)^2) / (N.j*c.hat[idx.nzero]) )
  #   )
  #stat_ = sum( ( ( h.j$counts[idx.nzero] - c.hat[idx.nzero]*N.j)^2 ) / (N.j*c.hat[idx.nzero]) )
 # stat_ <- chisq.test(h.j$counts[idx.nzero], p = c.hat[idx.nzero])$statistic
#  stat_
  #dgf = nzero - 1
  #pchisq(stat_, df = dgf)
  
#}


#Chi.squared.score <- function(m.u_, m.j_, M.max = 15){
#  m.u_ <- m.u_[m.u_ <= M.max]
#  m.j_ <- m.j_[m.j_ <= M.max]
#  N.j <- length(m.j_)
#  N.u <- length(m.u_)
#  m.breaks <- seq((M0 - 0.1), (M.max + 0.1 ), by = 0.1)
#  h.j <- hist(m.j_, plot = FALSE, breaks = m.breaks)
#  h.u <- hist(m.u_, plot = FALSE, breaks = m.breaks)
#  c.hat <- (h.u$counts + h.j$counts + 1)/(N.j + N.u)
#  #idx.nzero <- c.hat > 0 ; nzero <- sum(idx.nzero)
#  out.tosum <- ((h.u$counts - N.u*c.hat)^2)/N.u*c.hat + ((h.j$counts - N.j*c.hat)^2)/N.j*c.hat
#  sum(out.tosum)
#}


Chi.squared.score_plus1 <- function(m.u_, m.j_, M.max = 15){
  m.u_ <- m.u_[m.u_ <= M.max]
  m.j_ <- m.j_[m.j_ <= M.max]
  
  
  m.breaks <- seq((M0 - 0.1), (M.max + 0.1 ), by = 0.1)
  
  h.j <- hist(m.j_, plot = FALSE, breaks = m.breaks)
  h.u <- hist(m.u_, plot = FALSE, breaks = m.breaks)
  
  h.j_1 = h.j$counts + 1
  h.u_1 = h.u$counts + 1
  
  freq.u = (h.u_1/sum(h.u_1))*sum(h.j_1)
  freq.j = h.j_1
  
  chi_stat = ((freq.j - freq.u)^2 )/freq.u
  return(sum(chi_stat))
  
}

Chi.squared.score <- function(m.u_, m.j_, M.max = 15){
  m.u_ <- m.u_[m.u_ <= M.max]
  m.j_ <- m.j_[m.j_ <= M.max]
 
  
  m.breaks <- seq((M0 - 0.1), (M.max + 0.1 ), by = 0.1)
  
  h.j <- hist(m.j_, plot = FALSE, breaks = m.breaks)
  h.u <- hist(m.u_, plot = FALSE, breaks = m.breaks)
  
  h.j_1 = h.j$counts 
  h.u_1 = h.u$counts  
  
  freq.u = (h.u_1/sum(h.u_1))*sum(h.j_1)
  freq.j = h.j_1
  
  idx.nozero = freq.u != 0
  
  chi_stat = ((freq.j[idx.nozero] - freq.u[idx.nozero])^2 )/freq.u[idx.nozero]
  return(sum(chi_stat))

}


Chi.squared.score_mod <- function(m.u_, m.j_, M.max = 15){
  m.u_ <- m.u_[m.u_ <= M.max]
  m.j_ <- m.j_[m.j_ <= M.max]
  
  
  m.breaks <- seq((M0 - 0.1), (M.max + 0.1 ), by = 0.1)
  
  h.j <- hist(m.j_, plot = FALSE, breaks = m.breaks)
  h.u <- hist(m.u_, plot = FALSE, breaks = m.breaks)
  
  h.j_1 = h.j$counts 
  h.u_1 = h.u$counts  
  
  freq.u = (h.u_1/sum(h.u_1))*sum(h.j_1)
  freq.j = h.j_1
  
  chi_stat = ((freq.j - freq.u)^2 )
  return(sum(chi_stat))
  
}


log.d.multin <- function(x, size, prob){
  lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
}

Multinom.loglik.ratio <- function(m.u_, m.j_){
  N.j <- length(m.j_)
  N.u <- length(m.u_)
  m.breaks <- seq((M0 - 0.1), ( max(c(m.u_, m.j_)) + 0.1 ), by = 0.1)
  h.j <- hist(m.j_, plot = FALSE, breaks = m.breaks)
  h.u <- hist(m.u_, plot = FALSE, breaks = m.breaks)
  
  pr.Uj <- ( ( h.u$counts + (N.u/N.j) ) + (h.j$counts + 1) )/sum(( ( h.u$counts + (N.u/N.j) ) + (h.j$counts + 1) ))
  pr.U <- ( h.u$counts + (N.u/N.j)  )/sum( h.u$counts + (N.u/N.j) )
  pr.j <- ( h.j$counts + 1  )/sum( h.j$counts + 1 )
  
  log.lik.Uj <- sum(log.d.multin(x = ( ( h.u$counts + (N.u/N.j) ) + (h.j$counts + 1) ),
                                 size =  sum(( ( h.u$counts + (N.u/N.j) ) + (h.j$counts + 1) )), 
                                 prob = pr.Uj))
  log.lik.U <- sum(log.d.multin(x = ( h.u$counts + (N.u/N.j) ),
                                size =  sum( ( h.u$counts + (N.u/N.j)) ), 
                                prob = pr.U))
  log.lik.j <- sum(log.d.multin(x = ( h.j$counts + 1 ),
                                size =  sum( ( h.j$counts + 1 ) ), 
                                prob = pr.j))
  return(2*(log.lik.Uj - log.lik.U - log.lik.j))
}


Hellinger.stat <- function(m.u_, m.j_, m.breaks_){
  N.j <- length(m.j_)
  N.u <- length(m.u_)
  h.j <- hist(m.j_, plot = FALSE, breaks = m.breaks_)
  h.u <- hist(m.u_, plot = FALSE, breaks = m.breaks_)
  c.hat <- (h.u$counts)/(N.u)
  idx.nzero <- (c.hat > 0 | h.j$counts > 0); nzero <- sum(idx.nzero)
  c.j <- h.j$counts[idx.nzero]/N.j
  c.hat.pos <- c.hat[idx.nzero]
  (1/sqrt(2))*sqrt(sum( (sqrt(c.j) - sqrt(c.hat.pos))^2 ))
}


Brier.score <- function(m.u_, m.j_, m.breaks_){
  dummy.j <- lapply(1:length(m.j_), \(ii)
    sapply(1:(length(m.breaks_) - 1), \(xx) m.j_[ii] > m.breaks_[xx] & m.j_[ii] <= m.breaks_[xx + 1])
  )
  dummy.j <- do.call(rbind, dummy.j)
  with.obs <- colSums(dummy.j) > 0
  N.u <- length(m.u_)
  h.u <- hist(m.u_, plot = FALSE, breaks = m.breaks_)
  c.hat <- (h.u$counts)/(N.u)
  idx.nzero <- which(c.hat > 0 | with.obs)
  c.hat.pos <- c.hat[idx.nzero]
  dummy.j.pos <- dummy.j[,idx.nzero]
  square.diff <- lapply(1:nrow(dummy.j.pos), \(rr)
    (c.hat.pos - dummy.j.pos[rr,])^2
  )
  square.diff <- do.call(c, square.diff)
  sum(square.diff)/length(m.j_)
}



