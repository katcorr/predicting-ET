################################################################
### create_sim_datasets_gen.R                                ###
### data simulated using known data-generating mechanism     ###
### kat correia                                              ###
################################################################

# dat_type = 1 for known data-generating mechanism;
#            2 for resampling of obs. LBs + random censoring applied;
#            3 for resampling all obs + random censoring; 
#            4 for resampling all obs + censoring associated w age

library("tidyverse")
# extraDistr package has function to generate zero-truncated Poisson RV
library("extraDistr")
library("mosaic")
library("xtable")
library("lubridate")

nwomen <- c(500, 1000)
R <- 1000

#####################################################
#### CREATE DATASETS VIA KNOWN MECHANISM         ####
#####################################################

pcens <- c(3.7, 2.3, 1.5)  # about 30%, 45%, 55%

set.seed(6192019)

# initialize data frames

sim.wom.dat <- data.frame(dat_type=NA, dat_pcens=NA, dat_nwomen=NA, dat_num=NA, id=NA, x1=NA, x2.1=NA, x2.2=NA
                          , x3=NA, numET=NA, C=NA, maxET_observed=NA, livebirth=NA
                          , censored=NA, totalETminus1=NA)

sim.emb.dat <- data.frame(dat_type=NA, dat_pcens=NA, dat_nwomen=NA, dat_num=NA, id=NA, x1=NA, x2.1=NA, x2.2=NA
                           , x3=NA, numET=NA, C=NA,                   livebirth=NA
                           , censored=NA, totalETminus1=NA
                           , ETid=NA, livebirthET=NA)

# generate simulated datasets

for (j in pcens){
  for (k in nwomen){
    for (i in 1:R){
      
      # woman-level dataset with woman-level characteristics at first cycle...
      sim.wom.dat0 <- data.frame(dat_type = 1, dat_pcens = j, dat_nwomen = k, dat_num = i
                                , id = c(1:k)
                                , x1 = rnorm(n=k,mean=0,sd=1)
                                , x2 = MASS::mvrnorm(n = k, mu=c(0,0)
                                                   , Sigma=matrix(c(1,0.5,0.5,1)
                                                                  ,nrow=2))
                                , x3 = rbinom(n=k, size=1, prob=0.7)) %>%
        # simulate number of embryos transferred per woman until LB depending on x1, x2, and x3
        mutate(numET = rtpois(k
                              , lambda=exp(log(1.5)+log(1.2)*x2.1 + log(1.5)*x2.2 + log(2)*x3)
                              , a = 0)
               # simulate censoring time
               , C = rtpois(k, lambda = j, a = 0)) %>%
        rowwise() %>%
        mutate(maxET_observed = min(numET, C)
               , livebirth = ifelse(maxET_observed == numET, yes=1, no=0)
               , censored = 1 - livebirth
               , totalETminus1 = maxET_observed - 1) %>%
        ungroup()
      
      # combine 
      sim.wom.dat <- rbind(sim.wom.dat, sim.wom.dat0)
      
      # create embryo-level dataset 
      sim.emb.dat0 <- sim.wom.dat0 %>%
        uncount(maxET_observed) %>%
        mutate(ETid = ave(x1, id,  FUN = seq_along)
               , livebirthET = ifelse(livebirth==1 & ETid==1, yes=1,no=0)) 
      
      # combine
      sim.emb.dat <- rbind(sim.emb.dat, sim.emb.dat0)
      
      remove(sim.wom.dat0, sim.emb.dat0)
    }
  }
}

# save simulated datasets
save(sim.wom.dat, sim.emb.dat, file="data/simdata_knownmech.Rdat")
      

