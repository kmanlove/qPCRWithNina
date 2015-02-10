#------------------------#
#-- RT PCR detection ----#
#------------------------#

#------------------------#
# 1) sample animals (draw blood) (big sample)
# 2) pick up animals in pipette from sample (small sample -- subject to true 0's)
# 3) run cycles ("false" zeros -- animals in small sample, but no amplification in 40 cycles)
#------------------------#

Ct <- min.exceed + 12

#-- by-cycle expansion function --#

Ct.fun <- function(N0, Eff, reps){
  # build storage objects
  N1 <- N2 <- N3 <- N4 <- N5 <- N6 <- N7 <- N8 <- N9 <- N10 <- N11 <- N12 <- min.exceed <- Nabove <- Nbelow <- Ct <- rep(NA, reps)
  Ngeq13 <- ExceedsFluorescentThreshold <- vector("list", reps)
  # specify exponent for each expansion
  # looping over replicates (not cycles)
  for(i in 1:reps){
  N1[i] <- N0 + rbinom(N0, 1, Eff)
  N2[i] <- N1[i] + rbinom(1, N1[i], Eff)
  N3[i] <- N2[i] + rbinom(1, N2[i], Eff)
  N4[i] <- N3[i] + rbinom(1, N3[i], Eff)
  N5[i] <- N4[i] + rbinom(1, N4[i], Eff)
  N6[i] <- N5[i] + rbinom(1, N5[i], Eff)
  N7[i] <- N6[i] + rbinom(1, N6[i], Eff)
  N8[i] <- N7[i] + rbinom(1, N7[i], Eff)
  N9[i] <- N8[i] + rbinom(1, N8[i], Eff)
  N10[i] <- N9[i] + rbinom(1, N9[i], Eff)
  N11[i] <- N10[i] + rbinom(1, N10[i], Eff)
  N12[i] <- N11[i] + rbinom(1, N11[i], Eff)
  
  # 13 and above, calculate N using exponent.
  Cycles <- seq(1, 33)
  Ngeq13[[i]] <- N12[i] * (1 + Eff) ^ (Cycles)
  ExceedsFluorescentThreshold[[i]] <- ifelse(Ngeq13[[i]] >= (4.9 * (10^11)), "T", "F")
  min.exceed[i] <- min(which(ExceedsFluorescentThreshold[[i]]  == "T")) 
  Nabove[i] <- Ngeq13[[i]][min.exceed[i]]
  Nbelow[i] <-  Ngeq13[[i]][min.exceed[i] - 1]
  Ct[i] <- min.exceed[i] + 12
  }
  Ct.out <- table(factor(Ct, levels = seq(32, 44)))
  Nmat <- cbind(rep(N0, reps), Nabove, Nbelow, Ct, Ct - 1)
  return(list(Ct = Ct[1:reps], min.exceed = min.exceed, Ct.out = Ct.out, Nmat = Nmat))
}


max.No <- 1000
rep.No <- 100
ct.out.tab <- matrix(NA, nrow = 13, ncol = max.No)
nmat.list <- ct.out.obj <- vector("list", max.No)
for(i in 1:max.No){
  ct.out.obj[[i]] <- Ct.fun(N0 = i, Eff = 0.96, reps = rep.No)
  ct.out.tab[ ,i]  <- ct.out.obj[[i]]$Ct.out
  nmat.list[[i]] <- ct.out.obj[[i]]$Nmat
  print(i)
}

nmat.large <- as.data.frame(do.call("rbind", nmat.list))
nmat.large$m <- (nmat.large$Nabove - nmat.large$Nbelow) / 1
nmat.large$b <- (nmat.large$Nbelow) - (nmat.large$m * nmat.large$V5)
nmat.large$CtSoln <- (4.9 * 10 ^{11} - nmat.large$b) / nmat.large$m

# bin over tenths in nmat.large$CtSoln
break.vec <- seq(0, 40, by = .1)
nmat.large$CtSoln.cut <- cut(nmat.large$CtSoln, breaks = break.vec)
bin.025 <- bin.975 <- rep(NA, length(break.vec))
for(i in 1:length(break.vec)){
  k <- subset(nmat.large, CtSoln.cut == as.character(levels(as.factor(nmat.large$CtSoln.cut))[i]))
  bin.025[i] <- quantile(k$V1, 0.025)
  bin.975[i] <- quantile(k$V1, 0.975)
}

plot(x = 0, y = .99, xlim = c(25, 40), ylim = c(.99, 2000), log = "y", xlab = "Cycle", ylab = "Initial copy number")
for(i in 1: length(break.vec)){
  segments(x0 = break.vec[i], x1 = break.vec[i], y0 = bin.025[i], y1 = bin.975[i], )
}


plot(x = 0, y = .99, ylim = c(25, 40), xlim = c(.99, 2000), log = "x", ylab = "Cycle", xlab = "Initial copy number")
for(i in 1: length(break.vec)){
  segments(y0 = break.vec[i], y1 = break.vec[i], x0 = bin.025[i], x1 = bin.975[i], )
}


plot(x = 0, y = .99, ylim = c(25, 40), xlim = c(.99, 2000), log = "x", ylab = "Cycle", xlab = "Initial copy number")
for(i in 1: length(break.vec)){
  segments(y0 = break.vec[i], y1 = break.vec[i], x0 = bin.025[i], x1 = bin.975[i], )
}

out.data <- as.data.frame(cbind(break.vec, bin.025, bin.975))
names(out.data) <- c("break", "Quant2.5", "Quant97.5")

#write.csv(out.data, "~/work/Kezia/Research/ShortTermCollabs/Nina/qPCR/SimData/CtTenthsQuantiles_15Nov2014.csv")
#write.csv(nmat.large, "~/work/Kezia/Research/ShortTermCollabs/Nina/qPCR/SimData/FullSimDataset_15Nov2014.csv")

#-- want a plot of No quantiles over Ct --#
#-- need to simulate fluorescence curves.... --#

#-- plot --#
image(x = 1:501, y = 1:14, z = t(ct.out.tab), col = heat.colors(10), yaxt = "n")
axis(side = 2, at = 1:14, labels = ((1:14) + 31))

image(x = 1:201, y = 1:14, z = t(ct.out.tab[, 1:201]), col = heat.colors(10), yaxt = "n")
axis(side = 2, at = 1:14, labels = ((1:14) + 31))

image(x = 1:101, y = 1:14, z = t(ct.out.tab[, 1:101]), col = heat.colors(10), yaxt = "n")
axis(side = 2, at = 1:14, labels = ((1:14) + 31))

image(x = 1:51, y = 1:14, z = t(ct.out.tab[, 1:51]), col = heat.colors(10), yaxt = "n")
axis(side = 2, at = 1:14, labels = ((1:14) + 31))

Ct.test.1 <- Ct.fun(N0 = 1, Eff = 0.96, reps = 1000)
table(Ct.test.1$Ct.out)

Ct.test.2 <- Ct.fun(N0 = 2, Eff = 0.96, reps = 10000)
table(Ct.test.2)

Ct.test.3 <- Ct.fun(N0 = 3, Eff = 0.96, reps = 10000)
table(Ct.test.3)

Ct.test.4 <- Ct.fun(N0 = 4, Eff = 0.96, reps = 10000)
table(Ct.test.4)

Ct.test.5 <- Ct.fun(N0 = 5, Eff = 0.96, reps = 10000)
table(Ct.test.5)

Ct.test.6 <- Ct.fun(N0 = 6, Eff = 0.96, reps = 10000)
table(Ct.test.6)

Ct.test.7 <- Ct.fun(N0 = 7, Eff = 0.96, reps = 10000)
table(Ct.test.7)

Ct.test.10 <- Ct.fun(N0 = 10, Eff = 0.96, reps = 10000)
table(Ct.test.10)

Ct.test.12 <- Ct.fun(N0 = 12, Eff = 0.96, reps = 10000)
table(Ct.test.12)

Ct.test.13 <- Ct.fun(N0 = 13, Eff = 0.96, reps = 10000)
table(Ct.test.13)

Ct.test.14 <- Ct.fun(N0 = 14, Eff = 0.96, reps = 10000)
table(Ct.test.14)

Ct.test.15 <- Ct.fun(N0 = 15, Eff = 0.96, reps = 10000)
table(Ct.test.15)

Ct.test.20 <- Ct.fun(N0 = 20, Eff = 0.96, reps = 10000)
table(Ct.test.20)

Ct.test.50 <- Ct.fun(N0 = 50, Eff = 0.96, reps = 10000)
table(Ct.test.50)

Ct.test.70 <- Ct.fun(N0 = 70, Eff = 0.96, reps = 10000)
table(Ct.test.70)

Ct.test.90 <- Ct.fun(N0 = 90, Eff = 0.96, reps = 10000)
table(Ct.test.90)

Ct.test.500 <- Ct.fun(N0 = 500, Eff = 0.96, reps = 10000)
table(Ct.test.500)


#-------------------#
#-- OLD JAGS model--#
#-------------------#

# start simple: just estimate fixed effects on lambda--#
sink("qpcr1.bug")
cat("
    model{
    # Priors and constraints
    # for true initial number of gene copies
    for(j in 1:J){
    beta.0[j] ~ dunif(-3, 3)
    }
    beta.1 ~ dunif(-3, 3)
    #       theta.0 ~ dunif(-3, 3)
    #       theta.1 ~ dunif(-3, 3)
    #       
    #       sigma.experimenter ~ dunif(0, 10)
    #       sigma.process ~ dunif(0, 10)
    #       sigma.plate ~ dunif(0, 10)
    #       sigma.obsprocess ~ dunif(0, 10)
    # 
    #       sigma2.experimenter <- pow(sigma.experimenter, 2)
    #       sigma2.process ~ pow(sigma.process, 2)
    #       sigma2.plate ~ pow(sigma.plate, 2)
    #       sigma2.obsprocess ~ pow(sigma.obsprocess, 2)
    # 
    #       tau.experimenter <- pow(sigma.experimenter, -2)
    #       tau.process ~ pow(sigma.process, -2)
    #       tau.plate ~ pow(sigma.plate, -2)
    #       tau.obsprocess ~ pow(sigma.obsprocess, -2)
    
    # likelihood 
    for(j in 1 : J){
    for(i in 1 : N){
    #         Nest.lambda[i] = exp(beta.0 + beta.1 * load[i] + b.experimenter[i]) * exp(theta.0 + theta.1 * Eff[t] + k.plate[t] + epsilon[t]) / (exp (theta.0 + theta.1 * Eff[t] + k.plate[t] + epsilon[t]) + 1)
    log(lambda[i, j]) <- beta.0[j] + beta.1 * cycle[i, j]
    Nest[i, j] ~ dpois(lambda[i, j])
    }
    }
    }
    ")


#----------------------------------------------------#
#-- Code to simulate full qPCR state-space process --#
#----------------------------------------------------#
# 1) Specify variance associated
  # a) experimenter
  # b) pipettings
  # c) plate
  # d) expansion (might not need this -- it might be simulated)
# 2) Calculate N0
# 3) Expand N0 appropriately

# Pipetting process model
Sig_expansion <- 1
Sig_pipette <- 2
Sig_plate <- 0.001

beta_0 <- 1.0
beta_1 <- 1.5
Eff <- 0.96

exp.re <- rnorm(1, 0, Sig_exp)
pipette.resid <- rnorm(1, 0, Sig_pipette)
plate.re <- rnorm(1, 1, Sig_plate)
load.in <- 2

lambda <- exp(beta_0 + beta_1 * load.in + exp.re + pipette.resid)
N0 <- rpois(1, lambda)

# qPCR expansion (e.g., observation) simulation
Ct.expand.fun <- function(N0, Eff, plate.re){
  N <- rep(NA, 43)
  N[1] <- N0 + rbinom(1, N0, Eff)
  for(i in 2 : 15){
#    N[i] <- N[i - 1] + rbinom(1, N[i - 1], Eff)
    N[i] <- N[i - 1] + rbinom(1, N[i - 1], min(Eff * plate.re, 1))
  }
  for(i in 16:43){
    N[i] <-  N[i - 1] * (1 + min(Eff * plate.re, 1))
  }
  return(list(N = N))
}

Ct.expand.test <- Ct.expand.fun(N0, Eff = 0.96, plate.re)

# build dataset
# 1) simulate fluorescence at each cycle
data.array <- array(NA, dim = c(43, 100, 100))
for(i in 1:100){
  for(j in 1:100){
    data.array[, j, i] <- Ct.expand.fun(N0 = i, Eff = 0.96, plate.re)$N
    print(paste(i, "_", j))
  }
}

# 2) extract fluoresence and cycle number for first cycle over 80,000
cycle <- N.data <- matrix(rep(NA, each = 10000), byrow = F, ncol = 100, nrow = 100)
   # cycle is a predictor in bugs model
   # N.data is response in bugs model
for(i in 1:100){
  for(j in 1:100){
    cycle[i, j] <- which(data.array[, j, i] >= 80000)[1]
    N.data[i, j] <- round(data.array[cycle[i, j], j, i], 0)
  }
}

# (random but) low cycle data:
lowcycle.data <- data.array[floor(runif(6, 12)), ,]

# qPCR expansion model
tau <- rnorm(1, 0, Sig_expansion)
plate.re <- rnorm(1, 0, Sig_plate)
n.cycles <- 30

# bugs model
require(rjags)
require(runjags)

#--------------------------------------#
#-- MODEL 0.0: fixed effects on lambda --#
#--------------------------------------#

sink("qpcr0.bug")
cat("
    model{
      # Priors and constraints
      beta ~ dnorm(2, 0.001)
      alpha ~ dnorm(0, 0.001)
    
    # prior test model; no likelihood

    }
    ", fill = T)
sink()

# Bundle data
qpcr0.data <- list(
)

# initial values
qpcr0.inits <- function(){
  list(
    alpha = rnorm(1, 0, .5),
    beta = runif(1, 0, 2)
  )
}


# parameters to monitor
qpcr0.params <- c("alpha", "beta")

# mcmc settings
ni <- 500
nt <- 3
nb <- 250
nc <- 3

# call JAGS from R
qpcr0.call <- jags.model("qpcr0.bug",
                          data = qpcr0.data,
                          inits = qpcr0.inits,
                          n.chains = nc,
                          n.adapt = nb
)

update(qpcr0.call, ni)

qpcr0.coda <- coda.samples(qpcr0.call,
                           qpcr0.params,
                         ni)

summary(qpcr0.coda)
plot(qpcr0.coda)

#----------------------------------------#
#-- MODEL 0.1: fixed effects on lambda --#
#----------------------------------------#

sink("qpcr01.bug")
cat("
    model{
    # Priors and constraints
    beta ~ dnorm(2, 0.001)
    alpha ~ dnorm(0, 0.001)
    
    # likelihood (but no N0 effects yet, so no loop over j)
     for(i in 1 : N){
      log(lambda[i]) <- alpha + beta * cycle[i]
      Nest[i] ~ dpois(lambda[i])
    }
   }
    ", fill = T)
sink()

# Bundle data
qpcr01.data <- list(
    cycle = cycle[1, ],
    Nest = N.data[1, ],
    N = dim(cycle)[1]
)

# initial values
qpcr01.inits <- function(){
  list(
    alpha = rnorm(1, 0, .5),
    beta = runif(1, 0, 2)
  )
}

# parameters to monitor
qpcr01.params <- c("alpha", "beta")

# mcmc settings
ni <- 10000
nt <- 3
nb <- 5000
nc <- 3

# call JAGS from R
qpcr01.call <- jags.model("qpcr01.bug",
                         data = qpcr01.data,
                         inits = qpcr01.inits,
                         n.chains = nc,
                         n.adapt = nb
)

update(qpcr01.call, ni)

qpcr01.coda <- coda.samples(qpcr01.call,
                           qpcr01.params,
                           ni)

summary(qpcr01.coda)
plot(qpcr01.coda)


#--------------------------------------------------#
#-- MODEL 2: separate intercepts for each N0 sim --#
#--------------------------------------------------#

sink("qpcr2.bug")
cat("
    model{
    # Priors and constraints
    beta ~ dnorm(0, tau.beta)
    tau.beta <- pow(sigma.beta, -2)
    sigma.beta ~ dunif(0, 10)
    sigma2.beta <- pow(sigma.beta, 2)

    for(j in 1 : J){
      alpha[j] ~ dnorm(0, 0.00001)
    }

    # likelihood (both over N and over J this time)
    for(j in 1 : J){
    for(i in 1 : N){
    log(lambda[i, j]) <- alpha[j] + beta * cycle[i, j]
    Nest[i, j] ~ dpois(lambda[i, j])
    }
    }
    }
    ", fill = T)
sink()

# Bundle data
qpcr2.data <- list(
  cycle = cycle,
  Nest = N.data,
  N = dim(N.data)[1],
  J = dim(N.data)[2]
)

# initial values
qpcr2.inits <- function(){
  list(
    alpha = rep(0, 100),
    beta = 1.96,
    sigma.beta = runif(1, 0, 10)
  )
}


# parameters to monitor
qpcr2.params <- c("alpha", "beta", "sigma.beta")

# mcmc settings
ni <- 5000
nt <- 3
nb <- 2500
nc <- 3

# call JAGS from R
qpcr2.call <- jags.model("qpcr2.bug",
                         data = qpcr2.data,
                         inits = qpcr2.inits,
                         n.chains = nc,
                         n.adapt = nb
)

update(qpcr2.call, ni)

qpcr2.coda <- coda.samples(qpcr2.call,
                           qpcr2.params,
                           ni)

summary(qpcr2.coda) # misattributing values to alpha instead of beta....
plot(qpcr2.coda)


#-----------------------------------------------------------------#
#-- MODEL 3: Split efficiency and N0 model from expansion model --#
#-----------------------------------------------------------------#

sink("qpcr3.bug")
cat("
    model{
    # Priors and constraints
      alpha.eff ~ dnorm(0, tau.alpha.eff)
     tau.alpha.eff <- pow(sigma.alpha.eff, -2)
     sigma.alpha.eff ~ dunif(0, 10)
     sigma2.alpha.eff <- pow(sigma.alpha.eff, 2)
    
    for(i in 1 : N){
    for(j in 1 : J){
    alpha.expansion[i, j] ~ dnorm(0, 0.00001)
    }
    }

    for(i in 1 : N){
      for(j in 1 : J){
        plate.re[i, j] <- mean.plate
      }
    }

    mean.plate ~ dunif(0, tau.plate)
    tau.plate <- pow(sigma.plate, -2)
    sigma.plate ~ dunif(0, 10)
    sigma2.plate <- pow(sigma.plate, 2)

    # likelihood 
    for(j in 1 : J){
    for(i in 1 : N){
      # efficiency model
      logit(eff[i, j]) <- alpha.eff + plate.re[i, j]
      # expansion model
      log(lambda[i, j]) <- alpha.expansion[i, j] + (eff[i, j] + 1) * cycle[i, j]
      Nest[i, j] ~ dpois(lambda[i, j])
    }
    }
    }
    ", fill = T)
sink()

# Bundle data
qpcr3.data <- list(
  cycle = cycle,
  Nest = N.data,
  N = dim(N.data)[1],
  J = dim(N.data)[2]
)

# initial values
qpcr3.inits <- function(){
  list(
    alpha.eff = rnorm(1, 0, .1),
    alpha.expansion = matrix(runif(dim(N.data)[2] * dim(N.data)[1], 0, 2), nrow = dim(N.data)[1], ncol = dim(N.data)[2], byrow = T)
  )
}


# parameters to monitor
qpcr3.params <- c("alpha.eff", "alpha.expansion", "mean.plate", "sigma.plate", "sigma.alpha.eff")

# mcmc settings
ni <- 3000
nt <- 3
nb <- 1500
nc <- 3

# call JAGS from R
qpcr3.call <- jags.model("qpcr3.bug",
                         data = qpcr3.data,
                         inits = qpcr3.inits,
                         n.chains = nc,
                         n.adapt = nb
)

update(qpcr3.call, ni)

qpcr3.coda <- coda.samples(qpcr3.call,
                           qpcr3.params,
                           ni)

summary(qpcr3.coda) # misattributing values to alpha instead of beta....

plot(qpcr3.coda)