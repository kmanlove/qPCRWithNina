# qPCR data analysis

#-----------------------------------------#
#-- Read in data and required-libraries --#
#-----------------------------------------#
require(rjags)

qpcr <- read.csv("~/work/Kezia/Research/ShortTermCollabs/Nina/qPCR/Data/CompiledData_09Feb2015.csv", header = T)

#-----------------#
#-- Format data --#
#-----------------#
n.treatments <- length(levels(factor(qpcr$sample_id))) # here, calling each sample its own treatment...
n.pipetters <- length(levels(factor(qpcr$user)))
n.assays <- length(levels(factor(qpcr$assay)))
n.plates <- length(levels(factor(qpcr$plate_id)))
n.wells <- dim(qpcr)[1] # assumes each row reflects its own well.
n.cycles <- 40

# Response: Ct (indicator for whether detection has occurred by the cth cycle)
Ct.in <- as.numeric(as.character(qpcr$ct))
Ct.data <- matrix(NA, nrow = dim(qpcr)[1], ncol = 40)
for(i in 1 : dim(qpcr)[1]){
  for(j in 1 : 40){
    Ct.data[i, j] <- ifelse(as.numeric(as.character(qpcr$ct))[i] <= j, 1, 0)
  }
}

Ct.data.vec <- as.vector(t(Ct.data))

# Biological covariate: pipetter
pipetter.ind <- as.numeric(as.factor(qpcr$user))
pipetter.ind.exp <- rep(pipetter.ind, each = 40)

# Biological covariate: treatment
trt.ind <- c(rep(1, 100), rep(2, dim(qpcr)[1] - 100))
trt.ind.exp <- rep(trt.ind, each = 40)

# Expansion covariate: assay
assay.ind <- as.numeric(factor(qpcr$assay))  
assay.ind.exp <- rep(assay.ind, each = 40)

# Expansion covariate: plates
plate.ind <- as.numeric(factor(qpcr$plate_id))
plate.ind.exp <- rep(plate.ind, each = 40)

# Expansion covariate: wells
well.ind <- factor(seq(1 : dim(qpcr)[1]))
well.ind.exp <- rep(well.ind, each = 40)

# Expansion covariate: cycles
cycle.ind <- seq(1:40)
cycle.ind.exp <- rep(cycle.ind, times = dim(qpcr)[1])

# Detection covariate: threshold

#------------------------------------------#
#-- Define objects for state-space model --#
#------------------------------------------#
qpcr.data <- list(
  assay.ind = assay.ind.exp, 
  plate.ind = plate.ind.exp, 
  well.ind = well.ind.exp, 
  treatment.ind = trt.ind.exp, 
  pipetter.ind = pipetter.ind.exp,
  n.assays = length(levels(factor(assay.ind.exp))),
  n.plates = length(levels(factor(plate.ind.exp))) ,
  n.wells = length(levels(factor(well.ind.exp))),
  n.pipetters = length(levels(factor(pipetter.ind.exp))),
  n.samples = dim(qpcr)[1],
  n.treatments = length(levels(factor(trt.ind.exp))),
  n.cycles = 40,
  threshold = 10000,
  Ct.data = Ct.data
  )

qpcr.inits <- function(){
  list(
    sigma.pipetter = runif(1, 0, 10),
    sigma.assay = runif(1, 0, 10), 
    sigma.plate = runif(1, 0, 10),
    sigma.well = runif(1, 0, 10),
    sigma.cycle = runif(1, 0, 10)
  )
}

# parameters to monitor
qpcr.params <- c("beta.treatment", "sigma.pipetter", "sigma.assay", "sigma.plate", "sigma.well", "sigma.cycle")

# mcmc settings
ni <- 500
nt <- 3
nb <- 250
nc <- 3

#-------------------------------------#
#-- Build State-space model in BUGS --#
#-------------------------------------#
-
# bug model
sink("qpcr.bug")
cat("
    model{
    # initial copy-number at each cycle (N.c), and number of gene copies that expand (N.expanding)
    for(s in 1 : n.samples){
     N.expanding[s, 1] ~ dpois(1000)
     N.c[s, 1] ~ dpois(1000)
    } #s
    
    #--------------#
    #--- PRIORS ---#
    #--------------#
    # BIOLOGICAL PROCESS
    # biological error 1: pipetter
    tau.pipetter <- pow(sigma.pipetter, -2)
    sigma.pipetter ~ dunif(0, 10)
    sigma.pipetter2 <- pow(sigma.pipetter, 2)
    
    # Specificy priors on pipetter random effects. 
    for(p in 1:n.pipetters){
    pipetter.re[p] ~ dnorm(0, tau.pipetter) # random pipetter effect
    }
    
    # specify prior on biological fixed effect (.treatment different treatments)
    for(t in 1 : n.treatments){
    beta.treatment[t] ~ dnorm(0, 0.01)T(-10, 10)
    }
    
    # EXPANSION PROCESS
    # observation error 1: between-assay error
    tau.assay <- pow(sigma.assay, -2)
    sigma.assay ~ dunif(0, 50)
    sigma2.assay <- pow(sigma.assay, 2)
    for(a in 1 : n.assays){
    assay.re[a] ~ dnorm(0, tau.assay) # random assay effect
    } #a
    
    # observation error 2: between-plate error
    tau.plate <- pow(sigma.plate, -2)
    sigma.plate ~ dunif(0, 50)
    sigma2.plate <- pow(sigma.plate, 2)
    for(p in 1 : n.plates){
    plate.re[p] ~ dnorm(0, tau.plate) # random plate effect
    } #p
    
    # observation error 3: between-well error (NESTED in plate)
    tau.well <- pow(sigma.well, -2)
    sigma.well ~ dunif(0, 50)
    sigma2.well <- pow(sigma.well, 2)
    for(s in 1 : n.wells){
    well.re[s] ~ dnorm(0, tau.well) # random well effect
    } #w
    
    # observation error 4: between-cycle error (NESTED in well)
    tau.cycle <- pow(sigma.cycle, -2)
    sigma.cycle ~ dunif(0, 50)
    sigma2.cycle <- pow(sigma.cycle, 2)
    for(c in 1 : n.cycles){
    cycle.re[c] ~ dnorm(0, tau.cycle) # random well effect
    } #c
    
    # observation fixed effect: efficiency
    beta.mean.efficiency ~ dnorm(0, 0.01)T(-10, 10)
    
    # Specificy observation effects. 
    for(s in 1 : n.samples){
    for(c in 1 : n.cycles){
      logit(expansion[s, c]) <- beta.mean.efficiency + assay.re[assay.ind[s]] + plate.re[plate.ind[s]] + well.re[well.ind[s]] + cycle.re[c]
    } #c
    } #s
    
    #-------------------#
    #--- LIKELIHOODS ---#
    #-------------------#
    # BIOLOGICAL PROCESS (pipetter + treatment effects)
    # get N_0 estimates for each sample: treatment + pipetter
    for(s in 1 : n.samples){
    log(phi.individ.samples[s]) <- beta.treatment[treatment.ind[s]] + pipetter.re[pipetter.ind[s]]
    } #s
    
    # EXPANSION PROCESS
    for(s in 1 : n.samples){
    for(c in 2 : n.cycles){
    N.expanding[s, c] ~ dbin(expansion[s, c - 1], N.c[s, c - 1]) 
    N.c[s, c] <-  N.c[s, c - 1] + N.expanding[s, c]
    } #c
    } #s 

    # DETECTION PROCESS
    for(s in 1 : n.samples){
    for(c in 1 : n.cycles){
    Ct.data[s, c] ~ dbern(mu.detect[s, c])
    mu.detect[s, c] <- ifelse(N.c[s, c] >= threshold, .9999, 0.0001)
    } #c
    } #s
    
    }", fill = T)
sink()
-
#------------------------------------#
#-- Call and run state-space model --#
#------------------------------------#
qpcr.call <- jags.model("qpcr.bug",
                          data = qpcr.data,
                          inits = qpcr.inits,
                          n.chains = nc,
                          n.adapt = nb
)

update(qpcr.call, ni)

qpcr.coda <- coda.samples(qpcr.call,
                          qpcr.params,
                          ni)

summary(qpcr.coda)
