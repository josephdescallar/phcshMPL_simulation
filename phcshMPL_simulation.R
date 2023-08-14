library(parallel)
library(MASS)
library(foreach)
library(doParallel)
library(doSNOW)
library(survival)
library(splines2)
library(gbm)
library(simMSM)
library(Hmisc)

#Function to generate data for simulation
gen.cric.data <- function(n.datasets = 1000, n.obs, prob_event, rho, lambda1, 
                          lambda2, betas, gammaL, gammaR){
  simdata <- list()
  for(sim in 1:n.datasets){
    b11 <- betas[[1]][1]
    b12 <- betas[[1]][2]
    b21 <- betas[[2]][1]
    b22 <- betas[[2]][2]
    seed <- sim
    set.seed(seed)
    #Generate U~unif(0,1)
    u <- runif(n.obs)
    #Generate x1 and x2
    x1 <- rnorm(n.obs, mean=0, sd=1)
    x2 <- rnorm(n.obs, mean=0, sd=1)
    x <- data.matrix(data.frame(x1,x2))
    #Calculate t random variable
    t <- (-log(u) / (lambda1*exp(x1*b11 + x2*b12) + 
                       lambda2*(exp(x1*b21 + x2*b22))))^(1/rho)
    Ul <- runif(n.obs)
    Ur <- runif(n.obs,Ul,1)
    Ue <- runif(n.obs)
    event <- ifelse(Ue < prob_event,1,ifelse(gammaL*Ul <= t & t <= gammaR*Ur,3,
                                             ifelse(t < gammaL*Ul,2,0)))
    #Generate V to determine which risk t belongs
    p1 <- (lambda1*exp(x1*b11 + x2*b12)) / ((lambda1*exp(x1*b11 + x2*b12) + 
                                               lambda2*(exp(x1*b21 + x2*b22))))
    p2 <- (lambda2*exp(x1*b21 + x2*b22)) / ((lambda1*exp(x1*b11 + x2*b12) + 
                                               lambda2*(exp(x1*b21 + x2*b22))))
    v <- runif(n.obs)
    time1 <- ifelse(event==0,gammaR*Ur,ifelse(event==1,t,ifelse(event==2,
                                                                gammaL*Ul,gammaL*Ul)))
    time2 <- ifelse(event==0,NA,ifelse(event==1,t,ifelse(event==2,NA,
                                                         gammaR*Ur)))
    tmid <- ifelse(event==3,rowMeans(cbind(time1,time2)),ifelse(event==2,
                                                                time1/2, time1))
    risk <- ifelse(event==0,NA,ifelse(v <= p1, 1, 2))
    risk_1 <- ifelse((event==1 | event==2 | event==3) & risk == 1, 1, 0)
    risk_2 <- ifelse((event==1 | event==2 | event==3) & risk == 2, 1, 0)
    risk2 <- ifelse(is.na(risk),0,risk)
    x <- data.frame(x1,x2)
    simdata[[sim]] <- data.frame("time"=time1,time2,x,event,risk,tmid,risk2,t,
                                 risk_1,risk_2)
  }
  simdata
}

#Function to display bias, std, cp for cox regression results
sim.results.cox <- function(object, beta, simdata){
  betalist = seBlist = biaslist = mselist = lowerlist = upperlist = list()
  for(sim in 1:simdata){
    betalist[[sim]] = object[[sim]]$coefficients
    biaslist[[sim]] = beta - betalist[[sim]]
    mselist[[sim]] = biaslist[[sim]]^2
    seBlist[[sim]] = sqrt(diag(object[[sim]]$var))
    lowerlist[[sim]] = betalist[[sim]] - 1.96*seBlist[[sim]]
    upperlist[[sim]] = betalist[[sim]] + 1.96*seBlist[[sim]]
  }
  lower = Reduce("rbind", lowerlist)
  upper = Reduce("rbind", upperlist)
  betamat = matrix(beta, nrow = simdata, ncol = length(beta), byrow = TRUE)
  fit = list()
  fit$bias = colMeans(Reduce("rbind",biaslist))
  fit$seB = colMeans(Reduce("rbind",seBlist))
  fit$seBmc = apply(Reduce("rbind",betalist), 2, sd)
  fit$mse = colMeans(Reduce("rbind",mselist))
  fit$cov.prob = colMeans(ifelse(betamat >= lower & betamat <= upper,1,0))
  fit
}

#Function to display bias, std, cp for mpl regression
sim.results.mpl <- function(object, beta, simdata){
  biaslist = valid = seBlist = betalist =  mselist = list()
  lowerlist = upperlist = seBsandlist = lowersandlist = list()
  uppersandlist = list()
  for(sim in 1:simdata){
    valid[[sim]] = object[sim,]$valid
    betalist[[sim]] = unlist(object[sim, ]$"beta")
    biaslist[[sim]] = unlist(beta) - betalist[[sim]]
    seBlist[[sim]] = unlist(object[sim, ]$"seB")
    seBsandlist[[sim]] = unlist(object[sim, ]$"seB.sand")
    mselist[[sim]] = biaslist[[sim]]^2
    lowerlist[[sim]] = betalist[[sim]] - 1.96*seBlist[[sim]]
    upperlist[[sim]] = betalist[[sim]] + 1.96*seBlist[[sim]]
    lowersandlist[[sim]] = betalist[[sim]] - 1.96*seBsandlist[[sim]]
    uppersandlist[[sim]] = betalist[[sim]] + 1.96*seBsandlist[[sim]]
  }
  valid.temp = unlist(valid)
  valid.index = which(valid.temp %in% 1)
  count=1
  biaslist.valid = seBlist.valid = betalist.valid = fit = mselist.valid = list()
  lowerlist.valid = upperlist.valid = seBsandlist.valid = list()
  lowersandlist.valid = uppersandlist.valid = list()
  for(sim in 1:simdata){
    if(valid.temp[sim]==1){
      biaslist.valid[[count]] = biaslist[[sim]]
      seBlist.valid[[count]] = seBlist[[sim]]
      betalist.valid[[count]] = betalist[[sim]]
      mselist.valid[[count]] = mselist[[sim]]
      lowerlist.valid[[count]] = lowerlist[[sim]]
      upperlist.valid[[count]] = upperlist[[sim]]
      seBsandlist.valid[[count]] = seBsandlist[[sim]]
      lowersandlist.valid[[count]] = lowersandlist[[sim]]
      uppersandlist.valid[[count]] = uppersandlist[[sim]]
      count = count + 1
    }
  }
  lower = Reduce("rbind", lowerlist.valid)
  upper = Reduce("rbind", upperlist.valid)
  betamat = matrix(unlist(beta), nrow = count-1, ncol = length(unlist(beta)), 
                   byrow = TRUE)
  lowersand = Reduce("rbind", lowersandlist.valid)
  uppersand = Reduce("rbind", uppersandlist.valid)
  fit$valid = mean(valid.temp)
  fit$bias = colMeans(Reduce("rbind", biaslist.valid))
  fit$seB.asymp = colMeans(Reduce("rbind", seBlist.valid))
  fit$seB.sand = colMeans(Reduce("rbind", seBsandlist.valid))
  fit$seB.mc = apply(Reduce("rbind", betalist.valid), 2, sd)
  fit$mse = colMeans(Reduce("rbind", mselist.valid))
  fit$cov.prob = colMeans(ifelse(betamat >= lower & betamat <= upper,1,0))
  fit$cov.prob.sand = colMeans(ifelse(betamat >= lowersand & betamat <= uppersand,1,0))
  fit
}

#Function to predict baseline hazard from MPL
pred.bh <- function(object, t.points, r, sand = FALSE, rho=0, lambda=0, 
                    var.max=99999, coef){
  psif <- function(x, bknots, iknots){
    mSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
            intercept = object$basis.intercept)}
  pred.psi = psif(t.points, object$b.knots, object$i.knots[[r]])
  pred.h0r = as.vector(pred.psi %*% object$"theta"[[r]]) 
  if(sand == FALSE){
    VarCovMat.theta = object$VarCovMat[object$theta.index[[r]],
                                       object$theta.index[[r]]]
  }
  else if(sand == TRUE){ #up to here
    VarCovMat.theta = object$sand[object$theta.index[[r]], object$theta.index[[r]]]
  }
  pred.h0r.var = diag(pred.psi %*% VarCovMat.theta %*% t(pred.psi))
  if(object$pos.def==1){
    pred.var = pred.h0r.var / pred.h0r^2
    pred.var[pred.var > var.max] <- var.max
    pred.h0r.lower = as.vector(pred.h0r) * exp(-1.96*sqrt(pred.var))
    pred.h0r.upper = as.vector(pred.h0r) * exp(1.96*sqrt(pred.var))
  }
  true.h0r <- lambda*rho*t.points^(rho-1)
  if(object$valid==1){
    rlist <- list("pred.h0r"=pred.h0r, "pred.h0r.lower"=pred.h0r.lower,
                  "pred.h0r.upper"=pred.h0r.upper, "true.h0r"=true.h0r)
  }
  else{
    rlist <- list("pred.h0r"=pred.h0r, "true.h0r"=true.h0r)
  }
  rlist
}

#Function to predic CIF from MPL
pred.CIF <- function(object, t.points, sand=FALSE, rho, lambda, var.max, r){
  psif <- function(x, bknots, iknots){
    mSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
            intercept = object$basis.intercept)}
  PSIf <- function(x, bknots, iknots)
    iSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
            intercept = object$basis.intercept)
  bma = (t.points - object$b.knots[1]) / 2
  bpa = (t.points + object$b.knots[1]) / 2 #up to here 28
  t.gq.change = bma*matrix(object$nodes, nrow=length(bma), 
      ncol = object$gq.points, byrow=TRUE) + bpa
  pred.F0r.psi.gq = pred.F0r.PSI.gq = pred.F0r.h0qt.gq = list()
  pred.F0r.H0qt.gq = list()
  pred.F0r.S0qt.gq = pred.F0r.Integrand.gq = pred.F0r.dhdT.gq = list() 
  pred.F0r.dSdT.gq = pred.F0r.dFdT.gqT = list()
  for(gq in 1:object$gq.points){
    pred.F0r.psi.gq[[gq]] <- psif(t.gq.change[,gq], object$b.knots,
                                  object$i.knots[[r]])
    pred.F0r.PSI.gq[[gq]] <- PSIf(t.gq.change[,gq], object$b.knots,
                                  object$i.knots[[r]])
    pred.F0r.h0qt.gq[[gq]] <- psif(t.gq.change[,gq], object$b.knots,
                                   object$i.knots[[r]]) %*% object$"theta"[[r]]
    pred.F0r.H0qt.gq.r = pred.F0r.S0qt.gq.r = list()
    for(q in 1:object$n.risk){
      pred.F0r.H0qt.gq.r[[q]] = PSIf(t.gq.change[,gq], object$b.knots,
                                    object$i.knots[[q]]) %*% object$"theta"[[q]]
      pred.F0r.S0qt.gq.r[[q]] = exp(-pred.F0r.H0qt.gq.r[[q]])
    }
    pred.F0r.S0qt.gq[[gq]] = Reduce("*", pred.F0r.S0qt.gq.r)
    pred.F0r.Integrand.gq[[gq]] <- object$weights[gq]*pred.F0r.h0qt.gq[[gq]] *
      pred.F0r.S0qt.gq[[gq]]
    pred.F0r.dhdT.gq[[gq]] <- pred.F0r.psi.gq[[gq]]
    pred.F0r.dSdT.gq[[gq]] <- -pred.F0r.PSI.gq[[gq]] *
      as.vector(pred.F0r.S0qt.gq[[gq]])
    pred.F0r.dFdT.gqT[[gq]] <- (as.matrix(pred.F0r.dhdT.gq[[gq]] *
                                as.vector(pred.F0r.S0qt.gq[[gq]]) +
                                as.vector(pred.F0r.h0qt.gq[[gq]]) *
                                as.matrix(pred.F0r.dSdT.gq[[gq]]))) *
      object$weights[gq]
  }
  
  pred.F0r <- as.vector(bma*Reduce("+",pred.F0r.Integrand.gq))
  pred.F0r.dFdT <- bma*Reduce("+",pred.F0r.dFdT.gqT)
  if(sand == FALSE){
    VarCovMat.theta = object$VarCovMat[object$theta.index[[r]],
                                       object$theta.index[[r]]]
  }
  else if(sand == TRUE){
  VarCovMat.theta = object$sand[object$theta.index[[r]], object$theta.index[[r]]]
  }
  pred.F0r.r.var <- diag(pred.F0r.dFdT %*% as.matrix(VarCovMat.theta) %*%
                           t(pred.F0r.dFdT))
  pred.F0r.logOR <- log(pred.F0r / (1-pred.F0r + 1e-12) + 1e-12)
  pred.F0r.logOR.var = ((1/((1-pred.F0r)*pred.F0r))^2)*pred.F0r.r.var
  pred.F0r.log <- log(pred.F0r + 1e-12)
  pred.F0r.log.var <- (1/pred.F0r^2)*pred.F0r.r.var
  if(object$pos.def==1){
    pred.F0r.logOR.lower = pred.F0r.logOR - 1.96*sqrt(pred.F0r.logOR.var)
    pred.F0r.logOR.upper = pred.F0r.logOR + 1.96*sqrt(pred.F0r.logOR.var)
    pred.F0r.lower <- as.vector(exp(pred.F0r.logOR.lower) / 
                      (1 + exp(pred.F0r.logOR.lower)))
    pred.F0r.upper <- as.vector(exp(pred.F0r.logOR.upper) / 
                    (1 + exp(pred.F0r.logOR.upper)))
    
  }
  true.h0r <- lambda*rho*t.points^(rho-1)
  true.H0r <- lambda*(t.points^rho)
  true.S0r <- exp(-true.H0r)
  if(object$valid==1){
    rlist <- list("pred.F0r"=pred.F0r, "pred.F0r.lower"=pred.F0r.lower,
                  "pred.F0r.upper"=pred.F0r.upper)
  }
  else{
    rlist <- list("pred.F0r"=pred.F0r)
  }
  rlist
}

#Scenario 1, n = 200, left/interval = 47.5%, right cens = 47.5%
n.sim=1000
true.beta = list(c(-1,0.5), c(1,-0.5))
ic.scen1.r47.n200 <- gen.cric.data(n.datasets = 1000, n.obs = 200,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 0.89)
table(unlist(lapply(ic.scen1.r47.n200, function(a) a$event))) /
length(unlist(lapply(ic.scen1.r47.n200, function(a) a$event)))
#cox midpoint t with n = 200
ic.scen1.r47.n200.r1 = ic.scen1.r47.n200.r2 = list()
for(sim in 1:n.sim){
  ic.scen1.r47.n200.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
  data = ic.scen1.r47.n200[[sim]])
  ic.scen1.r47.n200.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
  data = ic.scen1.r47.n200[[sim]])
}
sim.results.cox(ic.scen1.r47.n200.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen1.r47.n200.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen1.r47.n200.time <- system.time({ic.scen1.r47.n200.results <- foreach(i = 1:n.sim, 
.combine = rbind, .options.snow = opts) %dopar% { phcsh_mpl3(Surv(time, time2, 
event=event, "interval") ~ x1 + x2, risk = ic.scen1.r47.n200[[i]]$risk, 
data = ic.scen1.r47.n200[[i]], max.outer = 10, max.iter = 3000, 
lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL, aps = TRUE,
knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen1.r47.n200.results, beta = true.beta, simdata = n.sim)


#Scenario 1, n = 200, left/interval = 75%, right cens = 20%
true.beta = list(c(-1,0.5), c(1,-0.5))
ic.scen1.r20.n200 <- gen.cric.data(n.datasets = 1000, n.obs = 200,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.34)
table(unlist(lapply(ic.scen1.r20.n200, function(a) a$event))) /
  length(unlist(lapply(ic.scen1.r20.n200, function(a) a$event)))
#cox midpoint t with n = 200
ic.scen1.r20.n200.r1 = ic.scen1.r20.n200.r2 = list()
for(sim in 1:n.sim){
  ic.scen1.r20.n200.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = ic.scen1.r20.n200[[sim]])
  ic.scen1.r20.n200.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = ic.scen1.r20.n200[[sim]])
}
sim.results.cox(ic.scen1.r20.n200.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen1.r20.n200.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen1.r20.n200.time <- system.time({ic.scen1.r20.n200.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% { phcsh_mpl3(Surv(time, time2, 
  event=event, "interval") ~ x1 + x2, risk = ic.scen1.r20.n200[[i]]$risk, #up to here
  data = ic.scen1.r20.n200[[i]], max.outer = 10, max.iter = 3000, 
  lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL, aps = TRUE,
  knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen1.r20.n200.results, beta = true.beta, simdata = n.sim)


#Scenario 1, n = 1000, left/interval = 47.5%, right cens = 47.5%
n.sim=1000
true.beta = list(c(-1,0.5), c(1,-0.5))
ic.scen1.r47.n1000 <- gen.cric.data(n.datasets = 1000, n.obs = 1000,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 0.89)
table(unlist(lapply(ic.scen1.r47.n1000, function(a) a$event))) /
  length(unlist(lapply(ic.scen1.r47.n1000, function(a) a$event)))
#cox midpoint t with n = 1000
ic.scen1.r47.n1000.r1 = ic.scen1.r47.n1000.r2 = list()
for(sim in 1:n.sim){
  ic.scen1.r47.n1000.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
  data = ic.scen1.r47.n1000[[sim]])
  ic.scen1.r47.n1000.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
  data = ic.scen1.r47.n1000[[sim]])
}
sim.results.cox(ic.scen1.r47.n1000.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen1.r47.n1000.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen1.r47.n1000.time <- system.time({ic.scen1.r47.n1000.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {
  phcsh_mpl3(Surv(time, time2, event=event, "interval") ~ x1 + x2, 
  risk = ic.scen1.r47.n1000[[i]]$risk, data = ic.scen1.r47.n1000[[i]], 
  max.outer = 10, max.iter = 3000, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL, aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, 
  iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen1.r47.n1000.results, beta = true.beta, simdata = n.sim)


#Scenario 1, n = 1000, left/interval = 75%, right cens = 20%
true.beta = list(c(-1,0.5), c(1,-0.5))
ic.scen1.r20.n1000 <- gen.cric.data(n.datasets = 1000, n.obs = 1000,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.34)
table(unlist(lapply(ic.scen1.r20.n1000, function(a) a$event))) /
  length(unlist(lapply(ic.scen1.r20.n1000, function(a) a$event)))
#cox midpoint t with n = 200
ic.scen1.r20.n1000.r1 = ic.scen1.r20.n1000.r2 = list()
for(sim in 1:n.sim){
  ic.scen1.r20.n1000.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = ic.scen1.r20.n1000[[sim]])
  ic.scen1.r20.n1000.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = ic.scen1.r20.n1000[[sim]])
}
sim.results.cox(ic.scen1.r20.n1000.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen1.r20.n1000.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen1.r20.n1000.time <- system.time({ic.scen1.r20.n1000.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {phcsh_mpl3(Surv(time, time2, 
  event=event, "interval") ~ x1 + x2, risk = ic.scen1.r20.n1000[[i]]$risk, 
  data = ic.scen1.r20.n1000[[i]], max.outer = 10, max.iter = 3000,
  lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL, aps = TRUE,
  knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen1.r20.n1000.results, beta = true.beta, simdata = n.sim)


#Scenario 2, n = 200, left/interval = 47.5%, right cens = 47.5%
true.beta = list(c(1,0.5), c(0.5,0.5))
ic.scen2.r47.n200 <- gen.cric.data(n.datasets = 1000, n.obs = 200,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.)
table(unlist(lapply(ic.scen2.r47.n200, function(a) a$event))) /
  length(unlist(lapply(ic.scen2.r47.n200, function(a) a$event)))
#cox midpoint t with n = 200
ic.scen2.r47.n200.r1 = ic.scen2.r47.n200.r2 = list()
for(sim in 1:n.sim){
  ic.scen2.r47.n200.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, data = ic.scen2.r47.n200[[sim]])
  ic.scen2.r47.n200.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, data = ic.scen2.r47.n200[[sim]])
}
sim.results.cox(ic.scen2.r47.n200.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen2.r47.n200.r2, beta = true.beta[[2]], simdata = n.sim)
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen2.r47.n200.time <- system.time({ic.scen2.r47.n200.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {
  phcsh_mpl3(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = ic.scen2.r47.n200[[i]]$risk, data = ic.scen2.r47.n200[[i]], 
  max.outer = 10, max.iter = 3000, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL, aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, 
  iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen2.r47.n200.results, beta = true.beta, simdata = n.sim)


#Scenario 2, n = 200, left/interval = 75%, right cens = 20%
true.beta = list(c(1,0.5), c(0.5,0.5))
ic.scen2.r20.n200 <- gen.cric.data(n.datasets = 1000, n.obs = 200,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.64)
table(unlist(lapply(ic.scen2.r20.n200, function(a) a$event))) /
  length(unlist(lapply(ic.scen2.r20.n200, function(a) a$event)))
#cox midpoint t with n = 200
ic.scen2.r20.n200.r1 = ic.scen2.r20.n200.r2 = list()
for(sim in 1:n.sim){
  ic.scen2.r20.n200.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
  data = ic.scen2.r20.n200[[sim]])
  ic.scen2.r20.n200.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
  data = ic.scen2.r20.n200[[sim]])
}
sim.results.cox(ic.scen2.r20.n200.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen2.r20.n200.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen2.r20.n200.time <- system.time({ic.scen2.r20.n200.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {
  phcsh_mpl3(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = ic.scen2.r20.n200[[i]]$risk, data = ic.scen2.r20.n200[[i]], 
  max.outer = 10, max.iter = 3000, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL, aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, 
  iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen2.r20.n200.results, beta = true.beta, simdata = n.sim)


#Scenario 2, n = 1000, left/interval = 47.5%, right cens = 47.5%
true.beta = list(c(1,0.5), c(0.5,0.5))
ic.scen2.r47.n1000 <- gen.cric.data(n.datasets = 1000, n.obs = 1000,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.)
table(unlist(lapply(ic.scen2.r47.n1000, function(a) a$event))) /
  length(unlist(lapply(ic.scen2.r47.n1000, function(a) a$event)))
#cox midpoint t with n = 1000
ic.scen2.r47.n1000.r1 = ic.scen2.r47.n1000.r2 = list()
for(sim in 1:n.sim){
  ic.scen2.r47.n1000.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
  data = ic.scen2.r47.n1000[[sim]])
  ic.scen2.r47.n1000.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
  data = ic.scen2.r47.n1000[[sim]])
}
sim.results.cox(ic.scen2.r47.n1000.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen2.r47.n1000.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen2.r47.n1000.time <- system.time({ic.scen2.r47.n1000.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {
  phcsh_mpl3(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = ic.scen2.r47.n1000[[i]]$risk, data = ic.scen2.r47.n1000[[i]],
  max.outer = 10, max.iter = 3000, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL, aps = TRUE, knots.perc.limit = c(0.05,0.95), 
  tmid = TRUE, iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen2.r47.n1000.results, beta = true.beta, simdata = n.sim)


#Scenario 2, n = 1000, left/interval = 75%, right cens = 20%
true.beta = list(c(1,0.5), c(0.5,0.5))
ic.scen2.r20.n1000 <- gen.cric.data(n.datasets = 1000, n.obs = 1000,
prob_event = 0.05, rho = 3, lambda1 = 1, lambda2 = 0.5,
betas = true.beta, gammaL = 0.5, gammaR = 1.64)
table(unlist(lapply(ic.scen2.r20.n1000, function(a) a$event))) /
  length(unlist(lapply(ic.scen2.r20.n1000, function(a) a$event)))
#cox midpoint t with n = 1000
ic.scen2.r20.n1000.r1 = ic.scen2.r20.n1000.r2 = list()
for(sim in 1:n.sim){
  ic.scen2.r20.n1000.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
  data = ic.scen2.r20.n1000[[sim]])
  ic.scen2.r20.n1000.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, #up to here
  data = ic.scen2.r20.n1000[[sim]])
}
sim.results.cox(ic.scen2.r20.n1000.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(ic.scen2.r20.n1000.r2, beta = true.beta[[2]], simdata = n.sim)

cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
ic.scen2.r20.n1000.time <- system.time({ic.scen2.r20.n1000.results <- foreach(i = 1:n.sim, 
  .combine = rbind, .options.snow = opts) %dopar% {
  phcsh_mpl3(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = ic.scen2.r20.n1000[[i]]$risk, data = ic.scen2.r20.n1000[[i]], 
  max.outer = 10, max.iter = 3000, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL, aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, 
  iter.disp = FALSE)
}
})
close(pb)
stopCluster(cl)
sim.results.mpl(ic.scen2.r20.n1000.results, beta = true.beta, simdata = n.sim)


#figure 1, baseline hazard and its coverage probabilities for scenario 1,
# int/left cens perc = 47.5%, right cens = 20%, n = 200
#caluclate values for plots
scen1.r47.n200.min <- max(sapply(ic.scen1.r47.n200, function(a)
  min(na.omit(c(a$time, a$time2)))))
scen1.r47.n200.max <- min(sapply(ic.scen1.r47.n200, function(a)
  max(na.omit(c(a$time, a$time2)))))
scen1.plot.t.r47 <- seq(scen1.r47.n200.min, scen1.r47.n200.max, length.out=1000)
scen1.plot.bh.r47 <- scen1.plot.bh.r47.r2 <- list()
for(sim in 1:1000){
  scen1.plot.bh.r47[[sim]] <- pred.bh(ic.scen1.r47.n200.results[sim,],
  scen1.plot.t.r47, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  scen1.plot.bh.r47.r2[[sim]] <- pred.bh(ic.scen1.r47.n200.results[sim,],
   scen1.plot.t.r47, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen1.plot.bh.r47.all <- lapply(scen1.plot.bh.r47, function(a) a$pred.h0r)
scen1.plot.bh.r47.h0r.mean <- colMeans(Reduce("rbind", scen1.plot.bh.r47.all))
scen1.plot.bh.r47.all.lower <- lapply(scen1.plot.bh.r47, 
function(a) a$pred.h0r.lower)
scen1.plot.bh.r47.h0r.mean.lower<- colMeans(Reduce("rbind",
scen1.plot.bh.r47.all.lower))
scen1.plot.bh.r47.all.upper <- lapply(scen1.plot.bh.r47, 
function(a) a$pred.h0r.upper)
scen1.plot.bh.r47.h0r.mean.upper<- colMeans(Reduce("rbind", 
scen1.plot.bh.r47.all.upper))
bh.cp <- colSums(Reduce("rbind", lapply(scen1.plot.bh.r47, 
function(a) ifelse((a$pred.h0r.lower < a$true.h0r)
& (a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000


scen1.plot.bh.r47.r2.all <- lapply(scen1.plot.bh.r47.r2, function(a) a$pred.h0r)
scen1.plot.bh.r47.r2.h0r.mean <- colMeans(Reduce("rbind", 
            scen1.plot.bh.r47.r2.all))
scen1.plot.bh.r47.r2.all.lower <- lapply(scen1.plot.bh.r47.r2, 
function(a) a$pred.h0r.lower)
scen1.plot.bh.r47.r2.h0r.mean.lower<- colMeans(Reduce("rbind", 
scen1.plot.bh.r47.r2.all.lower))
scen1.plot.bh.r47.r2.all.upper <- lapply(scen1.plot.bh.r47.r2, 
function(a) a$pred.h0r.upper)
scen1.plot.bh.r47.r2.h0r.mean.upper<- colMeans(Reduce("rbind", 
scen1.plot.bh.r47.r2.all.upper))
bh.cp.r2 <- lapply(scen1.plot.bh.r47.r2, function(a) 
ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0))
bh.cp.r2 <- colSums(Reduce("rbind", lapply(scen1.plot.bh.r47.r2, 
function(a) ifelse((a$pred.h0r.lower < a$true.h0r) & 
(a$true.h0r < a$pred.h0r.upper),1,0)))) / 1000

par(mfrow=c(2,2))
plot(scen1.plot.t.r47, scen1.plot.bh.r47.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1")
lines(scen1.plot.t.r47, scen1.plot.bh.r47.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47, scen1.plot.bh.r47.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47, scen1.plot.bh.r47[[1]]$true.h0r, col = "black", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')

plot(scen1.plot.t.r47, scen1.plot.bh.r47.r2.h0r.mean, type='l', col='#003366',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2")
lines(scen1.plot.t.r47, scen1.plot.bh.r47.r2.h0r.mean.lower, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47, scen1.plot.bh.r47.r2.h0r.mean.upper, col='#003366',
      lty='dashed')
lines(scen1.plot.t.r47, scen1.plot.bh.r47.r2[[1]]$true.h0r, 
      col = "black",cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=1, bty = 'n')


plot(scen1.plot.t.r47, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0.5,1), main="(C) Risk 1")
lines(scen1.plot.t.r47, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')
plot(scen1.plot.t.r47, bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='#003366', type='l', ylim=c(0.5,1), main="(D) Risk 2")
lines(scen1.plot.t.r47, rep(0.95, 1000), lty=2)
legend("bottomleft", legend=c("MPL", "0.95"), col=c("#003366","black"), 
       lty=c(1,2), bty = 'n')


#figure 2, CIF for scenario 1,
# int/left cens perc = 47.5%, right cens = 20%, n = 200
scen1.plot.cif.r47 <- scen1.plot.r2.r47.r2 <- list()
for(sim in 1:1000){
  scen1.plot.cif.r47[[sim]] <- pred.CIF(ic.scen1.r47.n200.results[sim,],
  scen1.plot.t.r47, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  
  scen1.plot.r2.r47.r2[[sim]] <- pred.CIF(ic.scen1.r47.n200.results[sim,],
  scen1.plot.t.r47, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}
scen1.plot.cif.r47.all <- lapply(scen1.plot.cif.r47, 
function(a) a$pred.F0r)
scen1.plot.cif.r47.F0r.mean <- colMeans(Reduce("rbind",
scen1.plot.cif.r47.all))
scen1.plot.cif.r47.all.lower <- lapply(scen1.plot.cif.r47, 
function(a) a$pred.F0r.lower)
scen1.plot.cif.r47.F0r.mean.lower<- colMeans(Reduce("rbind", 
scen1.plot.cif.r47.all.lower))
scen1.plot.cif.r47.all.upper <- lapply(scen1.plot.cif.r47, 
function(a) a$pred.F0r.upper)
scen1.plot.cif.r47.F0r.mean.upper<- colMeans(Reduce("rbind", 
scen1.plot.cif.r47.all.upper))

scen1.plot.cif.r47.r2.all <- lapply(scen1.plot.r2.r47.r2, 
function(a) a$pred.F0r)
scen1.plot.cif.r47.r2.F0r.mean <- colMeans(Reduce("rbind", 
scen1.plot.cif.r47.r2.all))
scen1.plot.cif.r47.r2.all.lower <- lapply(scen1.plot.r2.r47.r2, 
function(a) a$pred.F0r.lower)
scen1.plot.cif.r47.r2.F0r.mean.lower<- colMeans(Reduce("rbind", 
scen1.plot.cif.r47.r2.all.lower))
scen1.plot.cif.r47.r2.all.upper <- lapply(scen1.plot.r2.r47.r2, 
function(a) a$pred.F0r.upper)
scen1.plot.cif.r47.r2.F0r.mean.upper<- colMeans(Reduce("rbind", 
scen1.plot.cif.r47.r2.all.upper))

integrand.F0r <- function(x, lambda,rho){
  (lambda*rho*x^(rho-1))*exp(-lambda*(x^rho))*exp(-lambda*(x^rho))
}
true.F0r.r1 <- sapply(scen1.plot.t.r47, function(a) integrate(integrand.F0r, 0, 
               a, rho=3, lambda=1)$value)
true.F0r.r2 <- sapply(scen1.plot.t.r47, function(a) integrate(integrand.F0r, 0, 
               a, rho=3, lambda=0.5)$value)
F0r.plot.mat <- sapply(scen1.plot.cif.r47, function(a) a$pred.F0r)
F0r.plot.mean <- apply(F0r.plot.mat, 1, mean)
F0r.plot.l95 <- apply(F0r.plot.mat, 1, quantile, 0.05)
F0r.plot.u95 <- apply(F0r.plot.mat, 1, quantile, 0.95)
F0r.plot.r2.mat <- sapply(scen1.plot.r2.r47.r2, function(a) a$pred.F0r)
F0r.plot.r2.mean <- apply(F0r.plot.r2.mat, 1, mean)
F0r.plot.r2.l95 <- apply(F0r.plot.r2.mat, 1, quantile, 0.05)
F0r.plot.r2.u95 <- apply(F0r.plot.r2.mat, 1, quantile, 0.95)

par(mfrow=c(2,1))
plot(scen1.plot.t.r47, scen1.plot.cif.r47.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(A) Risk 1", lty=1)
lines(scen1.plot.t.r47, scen1.plot.cif.r47.F0r.mean.lower, col='#003366',
      lty=2)
lines(scen1.plot.t.r47, scen1.plot.cif.r47.F0r.mean.upper, col='#003366',
      lty=2)
lines(scen1.plot.t.r47, true.F0r.r1, col = "black", cex=10, lty=1)
lines(scen1.plot.t.r47, F0r.plot.l95, col='black',
      lty=2)
lines(scen1.plot.t.r47, F0r.plot.u95, col='black',
      lty=2)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')

plot(scen1.plot.t.r47, scen1.plot.cif.r47.r2.F0r.mean, type='l', col='#003366',
     xlab='t', ylab='Cumulative Incidence Function', ylim=c(0,1), 
     main="(B) Risk 2", lty=1)
lines(scen1.plot.t.r47, scen1.plot.cif.r47.r2.F0r.mean.lower, col='#003366',
      lty=2)
lines(scen1.plot.t.r47, scen1.plot.cif.r47.r2.F0r.mean.upper, col='#003366',
      lty=2)
lines(scen1.plot.t.r47, true.F0r.r2, col = "black", cex=10, lty=1)
lines(scen1.plot.t.r47, F0r.plot.r2.l95, col='black',
      lty=2)
lines(scen1.plot.t.r47, F0r.plot.r2.u95, col='black',
      lty=2)
legend("topleft", legend=c("MPL", "True"),
       col=c("#003366","black"), lty=c(1,1), bty = 'n')



