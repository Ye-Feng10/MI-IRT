library(mirt)
library(MASS)
library(tidyverse)
library(rlist)
#### functions ####
p <- function(t, a, d) 1 / (1 + exp(-outer(t, a, "*") - outer(rep(1, length(t)), d, "*")))

eap <- function(x, a, d, 
                qpts = seq(-4, 4, length = 33), 
                qwts = dnorm(seq(-4, 4, length = 33), 0, 1) / sum(dnorm(seq(-4, 4, length = 33), 0, 1))){
  P <- p(t = qpts, a = a, d = d)
  r <- outer(rep(1, length(qpts)), x, "*")
  L <- apply(P^r * (1 - P)^(1-r), 1, prod)
  est <- sum(qpts * L * qwts) / sum(L * qwts)
  psd <- sqrt(sum((qpts - est)^2 * L * qwts) / sum(L * qwts))
  c(est, psd)
}

eapMI <- function(x, MI, qpts = seq(-4, 4, length = 33), 
                  qwts = dnorm(seq(-4, 4, length = 33), 0, 1) / sum(dnorm(seq(-4, 4, length = 33), 0, 1))){
  eapse <- t(apply(MI, 3, function(mi)
    eap(x = x, a = mi[, 1], d = mi[, 2], qpts = qpts, qwts = qwts)))
  
  eap <- mean(eapse[, 1])
  psd <- sqrt(mean(eapse[, 2]^2) + (1+1/M)*sum((eapse[, 1] - eap)^2) / (M - 1))
  c(eap, psd)
}

inf <- function(t, a, d){
  P <- p(t = t, a = a, d = d)
  outer(rep(1, length(t)), a^2, "*") * P * (1 - P)
}

nextit <- function(idx, t, a, d){
  info <- inf(t = t, a = a, d = d)
  out <- idx[which.max(info)]
  return(out)
}

cat <- function(x, a, d, item1, endtype = "ni", endcrit = 20,
                qpts = seq(-4, 4, length = 33), 
                qwts = dnorm(seq(-4, 4, length = 33), 0, 1) / sum(dnorm(seq(-4, 4, length = 33), 0, 1))){
  i <- 1
  item <- item1
  th_c <- eap(x = x[item], a = a[item], d = d[item], qpts = qpts, qwts = qwts)
  theta <- th_c[1]
  sem <- th_c[2]
  
  if(endtype == "ni"){
    while(i < endcrit){
      i <- i + 1
      item[i] <- nextit(idx = (1:length(a))[-item], t = th_c[1], a = a[-item], d = d[-item])
      th_c <- eap(x = x[item], a = a[item], d = d[item], qpts = qpts, qwts = qwts)
      theta[i] <- th_c[1]
      sem[i] <- th_c[2]
    }
  }
  
  if(endtype == "se"){
    while(th_c[2] > endcrit & (i < I)){
      i <- i + 1
      item[i] <- nextit(idx = (1:length(a))[-item], t = th_c[1], a = a[-item], d = d[-item])
      th_c <- eap(x = x[item], a = a[item], d = d[item], qpts = qpts, qwts = qwts)
      theta[i] <- th_c[1]
      sem[i] <- th_c[2]
    }
  }
  
  list(theta = theta[i], sem = sem[i], thetal = theta, seml = sem, item = item)
}

catMI <- function(x, MI, a, d, item1, endtype = "ni", endcrit = 20,
                  nxi = "fixed", qtl = .5, 
                  qpts = seq(-4, 4, length = 33), 
                  qwts = dnorm(seq(-4, 4, length = 33), 0, 1) / sum(dnorm(seq(-4, 4, length = 33), 0, 1))){
  i <- 1
  item <- item1
  th_c <- eapMI(x = x[item], MI = MI[item, , , drop = FALSE], qpts = qpts, qwts = qwts)
  theta <- th_c[1]
  sem <- th_c[2]
  
  if(endtype == "ni"){
    while(i < endcrit){
      i <- i + 1
      if(nxi == "fixed") item[i] <- nextit(idx = (1:length(a))[-item], t = th_c[1], a = a[-item], d = d[-item])
      if(nxi == "MI") item[i] <- nextitMI(idx = (1:length(a))[-item], t = th_c[1], MI = MI[-item, , , drop = FALSE], qtl = qtl)
      th_c <- eapMI(x = x[item], MI = MI[item, , ], qpts = qpts, qwts = qwts)
      theta[i] <- th_c[1]
      sem[i] <- th_c[2]
    }
  }
  
  if(endtype == "se"){
    while(th_c[2] > endcrit & (i < I)){
      i <- i + 1
      info_i <- inf(t = th_c[1], a = a[-item], d = d[-item])
      item[i] <- (1:length(a))[-item][which.max(info_i)]
      th_c <- eapMI(x = x[item], MI = MI[item, , ], qpts = qpts, qwts = qwts)
      theta[i] <- th_c[1]
      sem[i] <- th_c[2]
    }
  }
  
  list(theta = theta[i], sem = sem[i], thetal = theta, seml = sem, item = item)
}

simfun <- function(seed){
  #### generate calibration data ####
  set.seed(seed)
  res <- NULL
  I <- 100
  N <- 200
  M <- 1000
  a <- rlnorm(I, 0, .5)
  b <- runif(I, -2.5, 2.5)
  d <- -a * b
  poptheta <- rnorm(N)
  popresp <- as.data.frame(t(sapply(poptheta, function(t) as.numeric(p(t, a, d) > runif(I)))))
  
  ## generate population CAT data
  th <- rep(seq(-3, 3, by = .5), each = 200)
  resp <- as.data.frame(t(sapply(th, function(t) as.numeric(p(t, a, d) > runif(I)))))
  
  ## fit model
  mod <- mirt(popresp, 1, "2PL", SE = TRUE, technical = list(NCYCLES = 100000))
  mu <- extract.mirt(mod, "parvec")
  sigma <- extract.mirt(mod, "vcov")
  ah <- mu[seq(1, 2*I - 1, by = 2)]
  dh <- mu[seq(2, 2*I, by = 2)]
  MI <- mvrnorm(n = M, mu = mu, Sigma = sigma)
  MI <- array(t(MI[, c(seq(1, 2*I - 1, by = 2), seq(2, 2*I, by = 2))]), dim = c(I, 2, M))
  
  
  ## run fixed-length cat with 20 items
  
  starttheta <- rnorm(length(th))
  startitem <- sapply(starttheta, function(t) which.max(inf(t = t, a = ah, d = dh)))
  
  pop_cat_20 <- lapply(1:nrow(resp), function(i) 
    cat(x = as.numeric(resp[i, ]), a = a, d = d, item1 = startitem[i]))
  cal_cat_20 <- lapply(1:nrow(resp), function(i) 
    cat(x = as.numeric(resp[i, ]), a = ah, d = dh, item1 = startitem[i]))
  btp_cat_20 <- lapply(1:nrow(resp), function(i) 
    catMI(x = as.numeric(resp[i, ]), MI = MI, a = ah, d = dh, item1 = startitem[i]))
  
  res_cat_20 <- tibble(i = rep(1:2600, 3),
                       type = factor(rep(c("pop", "cal", "btp"), each = 2600)),
                       theta = rep(th, 3),
                       that = c(sapply(cal_cat_20, "[[", "theta"), 
                                sapply(pop_cat_20, "[[", "theta"),
                                sapply(btp_cat_20, "[[", "theta")),
                       sem = c(sapply(cal_cat_20, "[[", "sem"), 
                               sapply(pop_cat_20, "[[", "sem"), 
                               sapply(btp_cat_20, "[[", "sem")))
  res_cat_20 <- res_cat_20 %>% mutate(tbias = theta - that)
  
  res[[1]] <- res_cat_20
  ## run variable-length cat with .4 stopping rule
  
  starttheta <- rnorm(length(th))
  startitem <- sapply(starttheta, function(t) which.max(inf(t = t, a = ah, d = dh)))
  
  pop_cat_4 <- lapply(1:nrow(resp), function(i) 
    cat(x = as.numeric(resp[i, ]), a = a, d = d, item1 = startitem[i], endtype = "se", endcrit = .4))
  cal_cat_4 <- lapply(1:nrow(resp), function(i) 
    cat(x = as.numeric(resp[i, ]), a = ah, d = dh, item1 = startitem[i], endtype = "se", endcrit = .4))
  btp_cat_4 <- lapply(1:nrow(resp), function(i) 
    catMI(x = as.numeric(resp[i, ]), MI = MI, a = ah, d = dh, item1 = startitem[i], endtype = "se", endcrit = .4))
  
  res_cat_4 <- tibble(i = rep(1:2600, 3),
                      type = factor(rep(c("pop", "cal", "btp"), each = 2600)),
                      theta = rep(th, 3),
                      that = c(sapply(cal_cat_4, "[[", "theta"), 
                               sapply(pop_cat_4, "[[", "theta"),
                               sapply(btp_cat_4, "[[", "theta")),
                      sem = c(sapply(cal_cat_4, "[[", "sem"), 
                              sapply(pop_cat_4, "[[", "sem"), 
                              sapply(btp_cat_4, "[[", "sem")),
                      ni = c(sapply(cal_cat_4, function(x) length(x$item)),
                             sapply(pop_cat_4, function(x) length(x$item)),
                             sapply(btp_cat_4, function(x) length(x$item))))
  res[[2]] <- res_cat_4
  
  filename <- paste(seed, ".rds", sep="")
  list.save(res,file = filename)
}

seeds <- sample(100000, 2)
system.time(sapply(seeds, simfun))
