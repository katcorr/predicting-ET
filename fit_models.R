################################################################
### fit_models.R                                             ###
### 8/9/19                                                   ###
### kat correia                                              ###
################################################################

library("tidyverse")
# extraDistr package has function to generate zero-truncated Poisson RV
library("extraDistr")
library("mosaic")
library("xtable")
# survival approach
library("survival")
library("survminer")
library("riskRegression")
# VGAM package has function to fit zero-truncated OR right-censored (not both) Poisson
library("VGAM")

nwomen <- c(500, 1000)
R <- 1000
up_to <- 50



#####################################################
####  DATASETS VIA KNOWN MECHANISM               ####
#####################################################

# load simulated datasets for known data generating mechanism
load("data/simdata_knownmech.Rdat")

pcens <- c(3.7, 2.3, 1.5)  # about 30%, 45%, 55%

pcens.obs <- array(NA, dim=c(R,length(nwomen),length(pcens))) 
RMSEs <- array(NA, dim=c(R,5,length(nwomen),length(pcens)))
MAEs <-  array(NA, dim=c(R,5,length(nwomen),length(pcens)))
underests <- array(NA, dim=c(R,20,length(nwomen),length(pcens)))
overests <- array(NA, dim=c(R,20,length(nwomen),length(pcens)))
equals <- array(NA, dim=c(R,5,length(nwomen),length(pcens)))


for (j in pcens){
  for (k in nwomen){
    for (i in 1:R){
      
      ck <- which(nwomen==k)
      cj <- which(pcens==j)
      
      use.wom.dat <- sim.wom.dat %>%
        filter(dat_pcens==j, dat_nwomen==k, dat_num==i)
      
      use.emb.dat <- sim.emb.dat %>%
        filter(dat_pcens==j, dat_nwomen==k, dat_num==i)
      
      pcens.obs[i,ck,cj] <- tally(~censored, data=use.wom.dat, format="percent")["1"]
      
      # 1. survival approach where time is #ET
      model1 <- coxph(Surv(time=maxET_observed, event=livebirth) ~ x1 + x2.1 + x2.2 + x3 
                      , data = use.wom.dat
                      , x=TRUE, y=TRUE)
      
      # 2b. right-censored Poisson model (re-characterized)
      model2b <- vglm(SurvS4(totalETminus1,livebirth) ~ x1 + x2.1 + x2.2 + x3
                      , family = cens.poisson()
                      , data = use.wom.dat)
      
      # 3. estimate the probability an individual embryo will yield LB,
      #             then use geometric dist'n (E(X)=1/p)
      model3 <- glm(livebirthET ~ x1 + x2.1 + x2.2 + x3
                    , data = use.emb.dat
                    , family = "binomial")
      
      #####################################################
      #### GET PREDICTIONS                             ####
      #####################################################
      pred0 <- use.wom.dat %>%
        select(id, x1, x2.1, x2.2, x3, numET) 
      
      # calculations for survival approach predictions
      m1.pred0 <- data.table::as.data.table(predictCox(model1, newdata=pred0
                                           , times=c(1:up_to), se=TRUE
                                           , type="survival")) %>%
        mutate(probLB = 1 - survival
               , probLB.lower = 1 - survival.upper
               , probLB.upper = 1 -survival.lower) %>%
        arrange(observation, times)
      
      m1.pred1a <-  m1.pred0[which(m1.pred0$times==up_to & 
                                     m1.pred0$probLB < 0.80)
                             ,"observation"]
      m1.pred1b <- m1.pred0[which(m1.pred0$times==up_to & 
                                    (m1.pred0$probLB >= 0.80 | is.na(m1.pred0$probLB)))
                            ,"observation"]
      m1.pred2 <- data.frame(id = c(m1.pred1a$observation, m1.pred1b$observation)
                             , prediction.model1 = c(rep(up_to,nrow(m1.pred1a))
                                        , min(times ~ observation
                                              , data=filter(m1.pred0,probLB >= 0.80 | is.na(m1.pred0$probLB))))) %>% 
        arrange(id)
      
      pred <- left_join(pred0, m1.pred2, by="id") %>%
        mutate(model2b.pred = as.vector(predict(model2b, newdata=pred0, type="response"))
               , prediction.model2b = model2b.pred + 1
               , predictionC.model2b = ceiling(prediction.model2b)
               , model3.prob = as.vector(predict(model3, newdata=pred0, type="response")) 
               , prediction.model3 = 1/model3.prob
               , predictionC.model3 = ceiling(prediction.model3)
               ##### differences
               , diff.model1 = prediction.model1 - numET
               , diff.model2b = prediction.model2b - numET
               , diffC.model2b = predictionC.model2b - numET
               , diff.model3 = prediction.model3 - numET
               , diffC.model3 = predictionC.model3 - numET
        )
      
        #####################################################
        #### SUMMARIES TO SAVE                           ####
        #####################################################
        
        RMSEs[i,,ck,cj] <- c(sqrt(sum(pred$diff.model1^2)/k)
                             , sqrt(sum(pred$diff.model2b^2)/k)
                             , sqrt(sum(pred$diffC.model2b^2)/k)
                             , sqrt(sum(pred$diff.model3^2)/k)
                             , sqrt(sum(pred$diffC.model3^2)/k))
        
        MAEs[i,,ck,cj] <- c(sum(abs(pred$diff.model1))/k
                            , sum(abs(pred$diff.model2b))/k
                            , sum(abs(pred$diffC.model2b))/k
                            , sum(abs(pred$diff.model3))/k
                            , sum(abs(pred$diffC.model3))/k)
      
        overests[i,,ck,cj] <- c(sum(pred$diff.model1 > 0)/k
                                , quantile(~diff.model1, p=c(0.25,0.50,0.75), data=filter(pred, diff.model1 > 0))
                                , sum(pred$diff.model2b > 0)/k
                                , quantile(~diff.model2b, p=c(0.25,0.50,0.75), data=filter(pred, diff.model2b > 0))
                                , sum(pred$diffC.model2b > 0)/k
                                , quantile(~diffC.model2b, p=c(0.25,0.50,0.75), data=filter(pred, diffC.model2b > 0))
                                , sum(pred$diff.model3 > 0)/k
                                , quantile(~diff.model3, p=c(0.25,0.50,0.75), data=filter(pred, diff.model3 > 0))
                                , sum(pred$diffC.model3 > 0)/k
                                , quantile(~diffC.model3, p=c(0.25,0.50,0.75), data=filter(pred, diffC.model3 > 0))
        )
      
        underests[i,,ck,cj] <- c(sum(pred$diff.model1 < 0)/k
                                 , quantile(~diff.model1, p=c(0.25,0.50,0.75), data=filter(pred, diff.model1 < 0))
                                 , sum(pred$diff.model2b < 0)/k
                                 , quantile(~diff.model2b, p=c(0.25,0.50,0.75), data=filter(pred, diff.model2b < 0))
                                 , sum(pred$diffC.model2b < 0)/k
                                 , quantile(~diffC.model2b, p=c(0.25,0.50,0.75), data=filter(pred, diffC.model2b < 0))
                                 , sum(pred$diff.model3 < 0)/k
                                 , quantile(~diff.model3, p=c(0.25,0.50,0.75), data=filter(pred, diff.model3 < 0))
                                 , sum(pred$diffC.model3 < 0)/k
                                 , quantile(~diffC.model3, p=c(0.25,0.50,0.75), data=filter(pred, diffC.model3 < 0))
        )
      
        equals[i,,ck,cj] <- c(sum(pred$diff.model1 == 0)/k
                              , sum(pred$diff.model2b == 0)/k
                              , sum(pred$diffC.model2b == 0)/k
                              , sum(pred$diff.model3 == 0)/k
                              , sum(pred$diffC.model3 == 0)/k)
        
    }
  }
}

save(pcens.obs, RMSEs, MAEs, underests, overests, equals
     , file="data/results_knownmech.Rda")

