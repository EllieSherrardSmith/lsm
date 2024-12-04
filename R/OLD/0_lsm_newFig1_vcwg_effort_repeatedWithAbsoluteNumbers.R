#################################
##
## VCWG Figure 
##
## sherrar 16 May 2024

## Simulations
library(malariasimulation)
library(malariaEquilibrium)


## 1 To show that residual transmission and resistance agnostic
dt = expand.grid(res = seq(0,1,by=0.1),
                 phiB = seq(0.4,0.99,length=6))

## add in phiI 
dt$phiI = rep(seq(0.45,1,length=6),each=11)

## Add in pyrethroid net estimates
## Lancet PH 2024 Churcher et al. 
# > pyr_only_ITN_parameters[c(1,11,21,31,41,51,61,71,81,91,101),]
# resistance     dn0       rn0    gamman
# 1          0.0 0.33760 0.6396419 3.8087149
# 11         0.1 0.33000 0.6457317 3.6559600
# 21         0.2 0.32070 0.6530264 3.4921798
# 31         0.3 0.31100 0.6604790 3.3166223
# 41         0.4 0.29925 0.6693460 3.1235534
# 51         0.5 0.28500 0.6797429 2.9091222
# 61         0.6 0.26835 0.6914028 2.6716935
# 71         0.7 0.24675 0.7056829 2.3971972
# 81         0.8 0.21750 0.7231990 2.0791716
# 91         0.9 0.16940 0.7467135 1.6698401
# 101        1.0 0.00000 0.7456314 0.7688844
dt$dn0 = rep(c(0.33760,0.33000,0.32070,0.31100,0.29925,0.28500,0.26835,0.24675,0.21750,0.16940,0.00000),6)
dt$rn0 = rep(c(0.6396419,0.6457317,0.6530264,0.6604790,0.6693460,0.6797429,0.6914028,0.7056829,0.7231990,0.7467135,0.7456314),6)
dt$gamman = rep(c(3.8087149,3.6559600,3.4921798,3.3166223,3.1235534,2.9091222,2.6716935,2.3971972,2.0791716,1.6698401,0.7688844),6)

## Add in pyrethroid-PBO
## Add in pyrethroid-pyrrole



###########################################
##
## function for running the model 

lsm_LOOP_f = function(EIR, parms_row,itn_cov){
  ## LSM modelling 
  
  year <- 365
  month <- 30
  sim_length <- 15 * year
  human_population <- 20000
  starting_EIR <- EIR  ## This will change to hit estimated prevalence in U5 yrs of (5%, 10%, 20%, 30%, 40%, 50%, and 60%)
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      
      prevalence_rendering_min_ages = c(0, 0.5, 0) * 365, ## Prev in under 5 years measured
      prevalence_rendering_max_ages = c(5,   5, 100) * 365,
      
      clinical_incidence_rendering_min_ages = c(0, 0.5,  0) * 365, ## All age clin_inc
      clinical_incidence_rendering_max_ages = c(5, 5,  100) * 365,
      
       model_seasonality = TRUE, # seasonal model
      g0 = 0.218087682,
      g = c(-0.005292592, -0.085277539, 0.017356449),
      h = c(0.174739899, -0.082337283, 0.026755121),
   
      individual_mosquitoes = FALSE, ## True by default
      
      bednets = TRUE
      
    )
  )
  
  # set treatment
  # 0.3453279	0.1794239 - from our generic estimates
  simparams <- set_drugs(simparams, list(AL_params,    ## whichever is ACT drug
                                         SP_AQ_params))## whichever is non-ACT drug
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,
                                      time=c(1),
                                      coverage=c(0.4))    # currently these are restricted but it would be a future hope to allow them to change
  simparams <- set_clinical_treatment(simparams, 
                                      drug=2,
                                      time=c(1),
                                      coverage=c(0.1))     # currently these are restricted but it would be a future hope to allow them to change
  
  # set mosquito species
  gam_params <- gamb_params # gambiae

  gam_params['phi_bednets'] <- dt$phiB[parms_row] # proportion biting in bed: update from Nilani (sensitivity analysis?)
  gam_params['phi_indoors'] <- dt$phiI[parms_row] # proportion biting indoors: update from Nilani (sensitivity analysis?)

  simparams <- set_species(simparams,
                           species=list(gam_params),
                           proportions=c(1)
                           
  )
  simparams <- set_equilibrium(simparams, starting_EIR)
  
  bednetparams <- simparams
  
  # set bednets
  bednet_events = data.frame(
    timestep = c(3,6,9,12,15)*365,
    names = rep("itn",5)
  )
  
  ## Bednet 1 corresponds to pyrethroid-only nets no resistance
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    coverages = c(0.6,0.6,itn_cov,itn_cov,itn_cov),
    
    retention = 4 * year, 
    
    dn0 = matrix(dt$dn0[parms_row], nrow=5, ncol=1),
    rn =  matrix(dt$rn[parms_row], nrow=5, ncol=1),
    rnm = matrix(0.24, nrow=5, ncol=1),
    gamman = rep(dt$gamman[parms_row] * 365, 5) 
    
  )
 
  correlationsb1 <- get_correlation_parameters(simparams)
  correlationsb1$inter_round_rho('bednets', 1)
  
  
  ## Specify the LSM coverage
  cc <- get_init_carrying_capacity(bednetparams_1)
  cc
  # 
  # reduction = (cc/human_population)*0.4 * human_population
  # print(reduction)
  reduction = 123220.1
  
  simparams <- bednetparams_1 |>
    set_carrying_capacity(
      carrying_capacity = t(matrix(c(cc-reduction),nrow = 1)),
      timesteps = c(9*365)
    )
 
  output_net1  <- run_simulation(sim_length, simparams,correlationsb1) 
  output_net1$pv_0_1825 = output_net1$n_detect_0_1825/output_net1$n_0_1825
  
  return(output_net1)
  
}

## Working out the absolute numbers of mosquitoes 
# PRINT cc for each species
# >   cc <- get_init_carrying_capacity(bednetparams_1)
# >   cc
# gamb      
# 308050.2

human_population
## This corresponds to 
(308050.2/human_population)
## And a 60% reduction in bites per person gives
## a corresponding
(308050.2/human_population)*0.4 ## bites per person

## this then works out as a absolute reduction in mosquitoes to
((308050.2/human_population)*0.4) * human_population

# aitnlsmabs = list()
# for(i in 1:nrow(dt)){
# 
#     aitnlsmabs[[i]] = lsm_LOOP_f(EIR = 5, ## 5 is ~10%;  20 is ~30%; 90 is ~60%
#                             i+50,
#                             itn_cov = 0.6)
#   print(i)
# }

##saveRDS(aitnlsmabs, "plots manuscript/lowtransmission_newFig1ITNandLSMAbsolutereduction.rds")
## No need to simulate this as it gives same results as the 60% reduction, by design
## I ran a bunch of random to check the reduction estimate was a fixed value

aitnlsmabs = list()
for(i in 1:nrow(dt)){
  ## we alter
  aitnlsmabs[[i]] = lsm_LOOP_f(EIR = 90, ## 5 is ~10%;  20 is ~30%; 90 is ~60%
                               i,
                               itn_cov = 0.6)
  print(i)
}
# saveRDS(aitnlsmabs, "plots manuscript/modtransmission_newFig1ITNandLSMAbsolutereduction.rds")
saveRDS(aitnlsmabs, "plots manuscript/higtransmission_newFig1ITNandLSMAbsolutereduction.rds")

##########################################
##
## Figure 1A 
##
## Show prevalence over time for do nothing, do ITN or do ITN+LSM
## 80% RESISTANCE AND 90% INDOOR BITING


## High transmission first
aa = readRDS("plots manuscript/higtransmission_newFig1comp.rds")
aitn = readRDS("plots manuscript/higtransmission_newFig1ITNonly.rds")
aitnlsm = readRDS("plots manuscript/higtransmission_newFig1ITNandLSM.rds")
aitnlsmabs = readRDS("plots manuscript/higtransmission_newFig1ITNandLSMAbsolutereduction.rds")

for(i in 1:nrow(dt)){
  dt$prevU5_n[i] = mean(aa[[i]]$pv_0_1825[3285:4380])
  dt$prevU5_itn[i] = mean(aitn[[i]]$pv_0_1825[3285:4380])
  dt$prevU5_itnLSM[i] = mean(aitnlsm[[i]]$pv_0_1825[3285:4380])
  dt$rel_reduction_due_to_itn[i] = (dt$prevU5_n[i] - dt$prevU5_itn[i])/dt$prevU5_n[i]
  dt$rel_reduction_due_to_itnlsm[i] = (dt$prevU5_n[i] - dt$prevU5_itnLSM[i])/dt$prevU5_n[i]
  dt$rel_reduction_if_addlsm_to_itn[i] = (dt$prevU5_itn[i] - dt$prevU5_itnLSM[i])/dt$prevU5_itn[i]
  dt$prevU5_itnLSMabsolute[i] = mean(aitnlsmabs[[i]]$pv_0_1825[3285:4380])
  dt$rel_reduction_due_to_itnlsmabs[i] = (dt$prevU5_n[i] - dt$prevU5_itnLSMabsolute[i])/dt$prevU5_n[i]
  dt$rel_reduction_if_addlsm_to_itnabs[i] = (dt$prevU5_itn[i] - dt$prevU5_itnLSMabsolute[i])/dt$prevU5_itn[i]
}

summary(dt$prevU5_n)
summary(dt$rel_reduction_due_to_itn)
summary(dt$rel_reduction_due_to_itnlsm)
summary(dt$rel_reduction_due_to_itnlsmabs)

ddt = subset(dt, dt$res < 0.91)
range(ddt$rel_reduction_due_to_itnlsm-ddt$rel_reduction_due_to_itn)
range(ddt$rel_reduction_due_to_itnlsmabs-ddt$rel_reduction_due_to_itn)

dt$rel_reduction_if_addlsm_to_itn_cum = dt$rel_reduction_if_addlsm_to_itn +dt$rel_reduction_due_to_itn
dt$rel_reduction_if_addlsm_to_itn_cumabs = dt$rel_reduction_if_addlsm_to_itnabs +dt$rel_reduction_due_to_itn

dt$rel_reduction_due_to_itn[which(dt$res == 0.8 &
                                    dt$phiB == 0.872)]
dt$rel_reduction_due_to_itn[53]
par(mfrow = c(3,3))
par(mar = c(4,8,2,1))

plot(aa[[53]]$pv_0_1825 ~ aa[[53]]$timestep,
     type="l",ylim=c(0,0.8),
     ylab = "falciparum prevalence in U5 years (%)",
     yaxt = "n",xaxt = "n",
     xlim = c(8*365,13*365),
     xlab = "Time in years")
text(9*365,0.78,"80% resistance")
text(9*365,0.70,"90% bites indoors")
axis(1, at = c(8,9,10,11,12,13,14,15)*365,
     labels = c(-1,0,1,2,3,4,5,6))
axis(2, las = 2, at = seq(0,0.8,0.1), labels = seq(0,80,10))
mtext("High Transmission", side = 2, line = 5, 
      cex = 1)

no_resistance = c(which(dt$res == 0))
which(dt$phiB == 0.872)

for(i in 53){
  lines(aitn[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="blue",lwd=3)
  lines(aitnlsm[[i]]$pv_0_1825 ~ aitnlsm[[i]]$timestep,col="aquamarine3",lwd=3)
  lines(aa[[i]]$pv_0_1825 ~ aa[[i]]$timestep,col="grey",lwd=3)
  lines(aitnlsmabs[[i]]$pv_0_1825 ~ aitnlsmabs[[i]]$timestep,col="orange",lwd=3)
}

polygon(c(9*365,9*365,12*365,12*365),
        c(0.05,0.15,0.15,0.05),col="grey65",border="white")
polygon(c(10*365,10*365,11*365,11*365),
        c(0.05,0.15,0.15,0.05),col="grey65",border="white")

##########################################
##
## Figure 1B  
##
## Show relative reduction from ITN and additional LSM impact
## against increasing resistance, fixing phiI as 0.9
par(mar = c(4,4,2,1))

plot(dt$rel_reduction_due_to_itn[dt$phiB == 0.872] ~ 
       dt$res[dt$phiB == 0.872],lwd=3,pch="",
     ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
     xlab = "Resistance (%)",xaxt="n",xlim = c(0,0.9)
     )
axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
values = loess((dt$rel_reduction_due_to_itn[dt$phiB == 0.872] ~ 
                  dt$res[dt$phiB == 0.872])) 

values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$phiB == 0.872] ~ 
                   dt$res[dt$phiB == 0.872])) 

values3 = loess((dt$rel_reduction_due_to_itnlsmabs[dt$phiB == 0.872] ~ 
                   dt$res[dt$phiB == 0.872])) 

polygon(c(dt$res[dt$phiB == 0.872],rev(dt$res[dt$phiB == 0.872])),
        c(predict(values2),rep(0,length(predict(values2)))),
        col = "aquamarine3",border=NA)

polygon(c(dt$res[dt$phiB == 0.872],rev(dt$res[dt$phiB == 0.872])),
        c(predict(values3),rep(0,length(predict(values3)))),
        col = "orange",border=NA)


polygon(c(dt$res[dt$phiB == 0.872],rev(dt$res[dt$phiB == 0.872])),
        c(predict(values),rep(0,length(predict(values)))),
        col = "blue",border=NA)

abline(v = 0.8,lty=2, lwd=3, col="grey15")

text(0.2,0.15,"Impact from ITN",col="white")
text(0.28,0.4,"Added impact from LSM",col="white")

##########################################
##
## Figure 1C  
##
## Show relative reduction from ITN and additional LSM impact
## against increasing phiI, resistance at 0.8

plot(dt$rel_reduction_due_to_itn[dt$res == 0.8] ~ 
       dt$phiB[dt$res == 0.8],lwd=3,pch="",
     ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
     xlab = "Proportion of bites in bed (%)",xaxt="n"
)
axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
values = loess((dt$rel_reduction_due_to_itn[dt$res == 0.8] ~ 
                  dt$phiB[dt$res == 0.8])) 

values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$res == 0.8] ~ 
                   dt$phiB[dt$res == 0.8])) 

values3 = loess((dt$rel_reduction_due_to_itnlsmabs[dt$res == 0.8] ~ 
                   dt$phiB[dt$res == 0.8])) 

polygon(c(dt$phiB[dt$res == 0.8],rev(dt$phiB[dt$res == 0.8])),
        c(predict(values2),rep(0,length(predict(values2)))),
        col = "aquamarine3",border=NA)

polygon(c(dt$phiB[dt$res == 0.8],rev(dt$phiB[dt$res == 0.8])),
        c(predict(values3),rep(0,length(predict(values3)))),
        col = "orange",border=NA)

polygon(c(dt$phiB[dt$res == 0.8],rev(dt$phiB[dt$res == 0.8])),
        c(predict(values),rep(0,length(predict(values)))),
        col = "blue",border=NA)

abline(v = 0.9,lty=2, lwd=3, col="grey15")

text(0.8,0.1,"Impact from ITN",col="white")
text(0.55,0.22,"Added impact from LSM",col="white")


####################
###
####
###### Repeat for moderate transmission
####
###
####################

##########################################
##
## Figure 1D 
##
## Show prevalence over time for do nothing, do ITN or do ITN+LSM
## 80% RESISTANCE AND 90% INDOOR BITING

## Moderate transmission first
aa = readRDS("plots manuscript/modtransmission_newFig1comp.rds")
aitn = readRDS("plots manuscript/modtransmission_newFig1ITNonly.rds")
aitnlsm = readRDS("plots manuscript/modtransmission_newFig1ITNandLSM.rds")
aitnlsmabs = readRDS("plots manuscript/modtransmission_newFig1ITNandLSMAbsolutereduction.rds")

for(i in 1:nrow(dt)){
  dt$prevU5_n[i] = mean(aa[[i]]$pv_0_1825[3285:4380])
  dt$prevU5_itn[i] = mean(aitn[[i]]$pv_0_1825[3285:4380])
  dt$prevU5_itnLSM[i] = mean(aitnlsm[[i]]$pv_0_1825[3285:4380])
  dt$rel_reduction_due_to_itn[i] = (dt$prevU5_n[i] - dt$prevU5_itn[i])/dt$prevU5_n[i]
  dt$rel_reduction_due_to_itnlsm[i] = (dt$prevU5_n[i] - dt$prevU5_itnLSM[i])/dt$prevU5_n[i]
  dt$rel_reduction_if_addlsm_to_itn[i] = (dt$prevU5_itn[i] - dt$prevU5_itnLSM[i])/dt$prevU5_itn[i]
  dt$prevU5_itnLSMabsolute[i] = mean(aitnlsmabs[[i]]$pv_0_1825[3285:4380])
  dt$rel_reduction_due_to_itnlsmabs[i] = (dt$prevU5_n[i] - dt$prevU5_itnLSMabsolute[i])/dt$prevU5_n[i]
  dt$rel_reduction_if_addlsm_to_itnabs[i] = (dt$prevU5_itn[i] - dt$prevU5_itnLSMabsolute[i])/dt$prevU5_itn[i]
}

summary(dt$prevU5_n)
summary(dt$rel_reduction_due_to_itn)
summary(dt$rel_reduction_due_to_itnlsm)
summary(dt$rel_reduction_due_to_itnlsmabs)
ddt = subset(dt, dt$res < 0.91)
range(ddt$rel_reduction_due_to_itnlsm-ddt$rel_reduction_due_to_itn)
range(ddt$rel_reduction_due_to_itnlsmabs-ddt$rel_reduction_due_to_itn)

dt$rel_reduction_due_to_itn[which(dt$res == 0.8 &
                                    dt$phiB == 0.872)]
dt$rel_reduction_due_to_itn[53]
par(mar = c(4,8,2,1))

plot(aa[[53]]$pv_0_1825 ~ aa[[53]]$timestep,
     type="l",ylim=c(0,0.8),
     ylab = "falciparum prevalence in U5 years (%)",
     yaxt = "n",xaxt = "n",
     xlim = c(8*365,13*365),
     xlab = "Time in years")
axis(1, at = c(8,9,10,11,12,13,14,15)*365,
     labels = c(-1,0,1,2,3,4,5,6))
axis(2, las = 2, at = seq(0,0.8,0.1), labels = seq(0,80,10))
mtext("Moderate Transmission", side = 2, line = 5, 
      cex = 1)

no_resistance = c(which(dt$res == 0))
which(dt$phiB == 0.872)

for(i in 53){
  lines(aitn[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="blue",lwd=3)
  lines(aitnlsm[[i]]$pv_0_1825 ~ aitnlsm[[i]]$timestep,col="aquamarine3",lwd=3)
  lines(aa[[i]]$pv_0_1825 ~ aa[[i]]$timestep,col="grey",lwd=3)
  lines(aitnlsmabs[[i]]$pv_0_1825 ~ aitnlsmabs[[i]]$timestep,col="orange",lwd=3)
}

##########################################
##
## Figure 1E  
##
## Show relative reduction from ITN and additional LSM impact
## against increasing resistance, fixing phiI as 0.9
par(mar = c(4,4,2,1))

plot(dt$rel_reduction_due_to_itn[dt$phiB == 0.872] ~ 
       dt$res[dt$phiB == 0.872],lwd=3,pch="",
     ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
     xlab = "Resistance (%)",xaxt="n",xlim=c(0,0.9)
)
axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
values = loess((dt$rel_reduction_due_to_itn[dt$phiB == 0.872] ~ 
                  dt$res[dt$phiB == 0.872])) 

values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$phiB == 0.872] ~ 
                   dt$res[dt$phiB == 0.872])) 

values3 = loess((dt$rel_reduction_due_to_itnlsmabs[dt$phiB == 0.872] ~ 
                   dt$res[dt$phiB == 0.872])) 

polygon(c(dt$res[dt$phiB == 0.872],rev(dt$res[dt$phiB == 0.872])),
        c(predict(values2),rep(0,length(predict(values2)))),
        col = "aquamarine3",border=NA)

polygon(c(dt$res[dt$phiB == 0.872],rev(dt$res[dt$phiB == 0.872])),
        c(predict(values3),rep(0,length(predict(values3)))),
        col = "orange",border=NA)


polygon(c(dt$res[dt$phiB == 0.872],rev(dt$res[dt$phiB == 0.872])),
        c(predict(values),rep(0,length(predict(values)))),
        col = "blue",border=NA)

abline(v = 0.8,lty=2, lwd=3, col="grey15")

text(0.2,0.15,"Impact from ITN",col="white")
text(0.28,0.5,"Added impact from LSM",col="white")

##########################################
##
## Figure 1F  
##
## Show relative reduction from ITN and additional LSM impact
## against increasing phiI, resistance at 0.8

plot(dt$rel_reduction_due_to_itn[dt$res == 0.8] ~ 
       dt$phiB[dt$res == 0.8],lwd=3,pch="",
     ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
     xlab = "Proportion of bites in bed (%)",xaxt="n"
)
axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
values = loess((dt$rel_reduction_due_to_itn[dt$res == 0.8] ~ 
                  dt$phiB[dt$res == 0.8])) 

values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$res == 0.8] ~ 
                   dt$phiB[dt$res == 0.8])) 

values3 = loess((dt$rel_reduction_due_to_itnlsmabs[dt$res == 0.8] ~ 
                   dt$phiB[dt$res == 0.8])) 

polygon(c(dt$phiB[dt$res == 0.8],rev(dt$phiB[dt$res == 0.8])),
        c(predict(values2),rep(0,length(predict(values2)))),
        col = "aquamarine3",border=NA)

polygon(c(dt$phiB[dt$res == 0.8],rev(dt$phiB[dt$res == 0.8])),
        c(predict(values3),rep(0,length(predict(values3)))),
        col = "orange",border=NA)

polygon(c(dt$phiB[dt$res == 0.8],rev(dt$phiB[dt$res == 0.8])),
        c(predict(values),rep(0,length(predict(values)))),
        col = "blue",border=NA)

abline(v = 0.9,lty=2, lwd=3, col="grey15")

text(0.8,0.1,"Impact from ITN",col="white")
text(0.55,0.22,"Added impact from LSM",col="white")


##########################################
##
## Figure 1G 
##
## Show prevalence over time for do nothing, do ITN or do ITN+LSM
## 80% RESISTANCE AND 90% INDOOR BITING

## Low transmission 
aa = readRDS("plots manuscript/lowtransmission_newFig1comp.rds")
aitn = readRDS("plots manuscript/lowtransmission_newFig1ITNonly.rds")
aitnlsm = readRDS("plots manuscript/lowtransmission_newFig1ITNandLSM.rds")

for(i in 1:nrow(dt)){
  dt$prevU5_n[i] = mean(aa[[i]]$pv_0_1825[3285:4380])
  dt$prevU5_itn[i] = mean(aitn[[i]]$pv_0_1825[3285:4380])
  dt$prevU5_itnLSM[i] = mean(aitnlsm[[i]]$pv_0_1825[3285:4380])
  dt$rel_reduction_due_to_itn[i] = (dt$prevU5_n[i] - dt$prevU5_itn[i])/dt$prevU5_n[i]
  dt$rel_reduction_due_to_itnlsm[i] = (dt$prevU5_n[i] - dt$prevU5_itnLSM[i])/dt$prevU5_n[i]
  dt$rel_reduction_if_addlsm_to_itn[i] = (dt$prevU5_itn[i] - dt$prevU5_itnLSM[i])/dt$prevU5_itn[i]
}

summary(dt$prevU5_n)
summary(dt$rel_reduction_due_to_itn)
summary(dt$rel_reduction_due_to_itnlsm)
ddt = subset(dt, dt$res < 0.91)
range(ddt$rel_reduction_due_to_itnlsm-ddt$rel_reduction_due_to_itn)

dt$rel_reduction_if_addlsm_to_itn_cum = dt$rel_reduction_if_addlsm_to_itn +dt$rel_reduction_due_to_itn

dt$rel_reduction_due_to_itn[which(dt$res == 0.8 &
                                    dt$phiB == 0.872)]

par(mar = c(4,8,2,1))

plot(aa[[53]]$pv_0_1825 ~ aa[[53]]$timestep,
     type="l",ylim=c(0,0.6),
     ylab = "falciparum prevalence in U5 years (%)",
     yaxt = "n",xaxt = "n",
     xlim = c(8*365,13*365),
     xlab = "Time in years")
text(9*365,0.78,"80% resistance")
text(9*365,0.70,"90% bites indoors")
axis(1, at = c(8,9,10,11,12,13,14,15)*365,
     labels = c(-1,0,1,2,3,4,5,6))
axis(2, las = 2, at = seq(0,0.6,0.1), labels = seq(0,60,10))
mtext("Low transmission",side = 2, cex = 1,line = 5)
no_resistance = c(which(dt$res == 0))
which(dt$phiB == 0.872)

for(i in 53){
  lines(aitn[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="blue",lwd=3)
  lines(aitnlsm[[i]]$pv_0_1825 ~ aitnlsm[[i]]$timestep,col="aquamarine3",lwd=3)
  lines(aa[[i]]$pv_0_1825 ~ aa[[i]]$timestep,col="grey",lwd=3)
  lines(aitnlsm[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="orange",lty = 2,lwd=3)
  
}

##########################################
##
## Figure 1H  
##
## Show relative reduction from ITN and additional LSM impact
## against increasing resistance, fixing phiI as 0.9
par(mar = c(4,4,2,1))

plot(dt$rel_reduction_due_to_itn[dt$phiB == 0.872] ~ 
       dt$res[dt$phiB == 0.872],lwd=3,pch="",
     ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
     xlab = "Resistance (%)",xaxt="n",xlim=c(0,0.9)
)
axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
values = loess((dt$rel_reduction_due_to_itn[dt$phiB == 0.872] ~ 
                  dt$res[dt$phiB == 0.872])) 

values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$phiB == 0.872] ~ 
                   dt$res[dt$phiB == 0.872])) 

polygon(c(dt$res[dt$phiB == 0.872],rev(dt$res[dt$phiB == 0.872])),
        c(predict(values2),rep(0,length(predict(values2)))),
        col = "orange",border=NA)

polygon(c(dt$res[dt$phiB == 0.872],rev(dt$res[dt$phiB == 0.872])),
        c(predict(values),rep(0,length(predict(values)))),
        col = "blue",border=NA)

text(0.2,0.2,"Impact from ITN",col="white")
text(0.28,0.70,"Added impact from LSM",col="white")

abline(v=0.8,col="grey15",lty = 2, lwd = 3)
##########################################
##
## Figure 1I  
##
## Show relative reduction from ITN and additional LSM impact
## against increasing phiI, resistance at 0.8

plot(dt$rel_reduction_due_to_itn[dt$res == 0.8] ~ 
       dt$phiB[dt$res == 0.8],lwd=3,pch="",
     ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
     xlab = "Proportion of bites in bed (%)",xaxt="n"
)
axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
values = loess((dt$rel_reduction_due_to_itn[dt$res == 0] ~ 
                  dt$phiB[dt$res == 0])) 

values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$res == 0.8] ~ 
                   dt$phiB[dt$res == 0.8])) 

polygon(c(dt$phiB[dt$res == 0.8],rev(dt$phiB[dt$res == 0.8])),
        c(predict(values2),rep(0,length(predict(values2)))),
        col = "orange",border=NA)

polygon(c(dt$phiB[dt$res == 0.8],rev(dt$phiB[dt$res == 0.8])),
        c(predict(values),rep(0,length(predict(values)))),
        col = "blue",border=NA)

text(0.8,0.2,"Impact from ITN",col="white")
text(0.55,0.5,"Added impact from LSM",col="white")

abline(v=0.9,col="grey15",lty = 2, lwd = 3)

par(xpd=NA,cex = 1.11)

text(x = -1.27, y = 3.52,"(A)",cex=0.8)
text(x = -0.53, y = 3.52,"(B)",cex=0.8)
text(x = 0.28, y = 3.52,"(C)",cex=0.8)
text(x = -1.27, y = 2.27,"(D)",cex=0.8)
text(x = -0.53, y = 2.27,"(E)",cex=0.8)
text(x = 0.28, y = 2.27,"(F)",cex=0.8)
text(x = -1.27, y = 0.9,"(G)",cex=0.8)
text(x = -0.53, y = 0.9,"(H)",cex=0.8)
text(x = 0.28, y = 0.9,"(I)",cex=0.8)
# #############################################################
# ##
# ## Next simulate increasing cover
# 
# dt2 = seq(0,1,length=11)
# 
# bb = list()
# bitnlsm = list()
# for(i in 1:length(dt2)){
#   bb[[i]] = lsm_LOOP_f(45,         ## this is resistance = 0, phiB = 0.872, pihI = 0.89
#                        #and try 54 as res = 90
#                        itn_cov = 0,
#                        lsm_impact = dt2[i])
#   
#   bitnlsm[[i]] = lsm_LOOP_f(45,
#                             itn_cov = 0.6,
#                             lsm_impact = dt2[i])
#   print(i)
# }
# 
# plot(bitnlsm[[1]]$pv_0_1825 ~ bitnlsm[[1]]$timestep, 
#      type="l",ylim=c(0,0.6),
#      ylab = "falciparum prevalence in U5 years (%)",
#      yaxt = "n",
#      xlim = c(6*365,13*365),
#      xaxt="n",
#      xlab = "Time in days from LSM")
# 
# cols = c("black","grey10","grey20","grey30","grey40","grey50",
#          "lightblue","blue","darkblue","purple","red")
# for(i in 1:11){
#   lines(bitnlsm[[i]]$pv_0_1825 ~ bitnlsm[[i]]$timestep,
#         col = cols[i])
# }
# axis(2, las = 2, at = seq(0,0.6,0.1), labels = seq(0,60,10))
# axis(1, at = c(5,7,9,11,13)*365, labels = c(-4,-2,0,2,4)*365)
# 
# legend("topright",title = "Acheived (and maintained) reduction in carrying capacity",
#        legend = c("None","10%","20%","30%","40%","50%","60%","70%","80%","90%","Total"),
#        lty=1,col=cols,ncol=3,bty="n")


## Repeat Figure 1 as a supplement with no resistance and lower biting

## Show prevalence over time for do nothing, do ITN or do ITN+LSM
## 80% RESISTANCE AND 90% INDOOR BITING

what_resistance_value = 0.0
what_phiB_value = 0.754
which_row = 34

draw_comparison_f = function(what_resistance_value,
                             what_phiB_value,
                             which_row)
{
  
  ##########################################
  ##
  ## Figure 1A 
  ##
  ## Show prevalence over time for do nothing, do ITN or do ITN+LSM
  ## 80% RESISTANCE AND 90% INDOOR BITING
  
  
  ## High transmission first
  aa = readRDS("plots manuscript/higtransmission_newFig1comp.rds")
  aitn = readRDS("plots manuscript/higtransmission_newFig1ITNonly.rds")
  aitnlsm = readRDS("plots manuscript/higtransmission_newFig1ITNandLSM.rds")
  aitnlsmabs = readRDS("plots manuscript/higtransmission_newFig1ITNandLSMAbsolutereduction.rds")
  
  for(i in 1:nrow(dt)){
    dt$prevU5_n[i] = mean(aa[[i]]$pv_0_1825[3285:4380])
    dt$prevU5_itn[i] = mean(aitn[[i]]$pv_0_1825[3285:4380])
    dt$prevU5_itnLSM[i] = mean(aitnlsm[[i]]$pv_0_1825[3285:4380])
    dt$rel_reduction_due_to_itn[i] = (dt$prevU5_n[i] - dt$prevU5_itn[i])/dt$prevU5_n[i]
    dt$rel_reduction_due_to_itnlsm[i] = (dt$prevU5_n[i] - dt$prevU5_itnLSM[i])/dt$prevU5_n[i]
    dt$rel_reduction_if_addlsm_to_itn[i] = (dt$prevU5_itn[i] - dt$prevU5_itnLSM[i])/dt$prevU5_itn[i]
    dt$prevU5_itnLSMabsolute[i] = mean(aitnlsmabs[[i]]$pv_0_1825[3285:4380])
    dt$rel_reduction_due_to_itnlsmabs[i] = (dt$prevU5_n[i] - dt$prevU5_itnLSMabsolute[i])/dt$prevU5_n[i]
    dt$rel_reduction_if_addlsm_to_itnabs[i] = (dt$prevU5_itn[i] - dt$prevU5_itnLSMabsolute[i])/dt$prevU5_itn[i]
  }
  
  summary(dt$prevU5_n)
  summary(dt$rel_reduction_due_to_itn)
  summary(dt$rel_reduction_due_to_itnlsm)
  summary(dt$rel_reduction_due_to_itnlsmabs)
  
  ddt = subset(dt, dt$res < 0.91)
  range(ddt$rel_reduction_due_to_itnlsm-ddt$rel_reduction_due_to_itn)
  range(ddt$rel_reduction_due_to_itnlsmabs-ddt$rel_reduction_due_to_itn)
  
  dt$rel_reduction_if_addlsm_to_itn_cum = dt$rel_reduction_if_addlsm_to_itn +dt$rel_reduction_due_to_itn
  dt$rel_reduction_if_addlsm_to_itn_cumabs = dt$rel_reduction_if_addlsm_to_itnabs +dt$rel_reduction_due_to_itn
  
  dt$rel_reduction_due_to_itn[which(dt$res == what_resistance_value &
                                      dt$phiB == what_phiB_value)]
 
  par(mfrow = c(3,3))
  par(mar = c(4,8,2,1))
  
  plot(aa[[53]]$pv_0_1825 ~ aa[[53]]$timestep,
       type="l",ylim=c(0,0.8),
       ylab = "falciparum prevalence in U5 years (%)",
       yaxt = "n",xaxt = "n",col="white",
       xlim = c(8*365,13*365),
       xlab = "Time in years")
  text(9*365,0.78,paste0(100*what_resistance_value, " % resistance"))
  text(9*365,0.70,paste0(100*what_phiB_value, " % bites in bed"))
  axis(1, at = c(8,9,10,11,12,13,14,15)*365,
       labels = c(-1,0,1,2,3,4,5,6))
  axis(2, las = 2, at = seq(0,0.8,0.1), labels = seq(0,80,10))
  mtext("High Transmission", side = 2, line = 5, 
        cex = 1)
  
  no_resistance = c(which(dt$res == what_resistance_value))
  which(dt$phiB == what_phiB_value)
  
  for(i in which_row){
    lines(aitn[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="blue",lwd=3)
    lines(aitnlsm[[i]]$pv_0_1825 ~ aitnlsm[[i]]$timestep,col="aquamarine3",lwd=3)
    lines(aa[[i]]$pv_0_1825 ~ aa[[i]]$timestep,col="grey",lwd=3)
    lines(aitnlsmabs[[i]]$pv_0_1825 ~ aitnlsmabs[[i]]$timestep,col="orange",lwd=3)
  }
  
  polygon(c(9*365,9*365,12*365,12*365),
          c(0.05,0.15,0.15,0.05),col="grey65",border="white")
  polygon(c(10*365,10*365,11*365,11*365),
          c(0.05,0.15,0.15,0.05),col="grey65",border="white")
  
  ##########################################
  ##
  ## Figure 1B  
  ##
  ## Show relative reduction from ITN and additional LSM impact
  ## against increasing resistance, fixing phiI as 0.9
  par(mar = c(4,4,2,1))
  
  plot(dt$rel_reduction_due_to_itn[dt$phiB == what_phiB_value] ~ 
         dt$res[dt$phiB == what_phiB_value],lwd=3,pch="",
       ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
       xlab = "Resistance (%)",xaxt="n",xlim = c(0,0.9)
  )
  axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
  axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
  values = loess((dt$rel_reduction_due_to_itn[dt$phiB == what_phiB_value] ~ 
                    dt$res[dt$phiB == what_phiB_value])) 
  
  values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$phiB == what_phiB_value] ~ 
                     dt$res[dt$phiB == what_phiB_value])) 
  
  values3 = loess((dt$rel_reduction_due_to_itnlsmabs[dt$phiB == what_phiB_value] ~ 
                     dt$res[dt$phiB == what_phiB_value])) 
  
  polygon(c(dt$res[dt$phiB == what_phiB_value],rev(dt$res[dt$phiB == what_phiB_value])),
          c(predict(values2),rep(0,length(predict(values2)))),
          col = "aquamarine3",border=NA)
  
  polygon(c(dt$res[dt$phiB == what_phiB_value],rev(dt$res[dt$phiB == what_phiB_value])),
          c(predict(values3),rep(0,length(predict(values3)))),
          col = "orange",border=NA)
  
  
  polygon(c(dt$res[dt$phiB == what_phiB_value],rev(dt$res[dt$phiB == what_phiB_value])),
          c(predict(values),rep(0,length(predict(values)))),
          col = "blue",border=NA)
  
  abline(v = what_resistance_value,lty=2, lwd=3, col="grey15")
  
  text(0.2,0.15,"Impact from ITN",col="white")
  text(0.28,0.4,"Added impact from LSM",col="white")
  
  ##########################################
  ##
  ## Figure 1C  
  ##
  ## Show relative reduction from ITN and additional LSM impact
  ## against increasing phiI, resistance at 0.8
  
  plot(dt$rel_reduction_due_to_itn[dt$res == what_resistance_value] ~ 
         dt$phiB[dt$res == 0.8],lwd=3,pch="",
       ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
       xlab = "Proportion of bites in bed (%)",xaxt="n"
  )
  axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
  axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
  values = loess((dt$rel_reduction_due_to_itn[dt$res == what_resistance_value] ~ 
                    dt$phiB[dt$res == what_resistance_value])) 
  
  values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$res == what_resistance_value] ~ 
                     dt$phiB[dt$res == what_resistance_value])) 
  
  values3 = loess((dt$rel_reduction_due_to_itnlsmabs[dt$res == what_resistance_value] ~ 
                     dt$phiB[dt$res == what_resistance_value])) 
  
  polygon(c(dt$phiB[dt$res == what_resistance_value],rev(dt$phiB[dt$res == what_resistance_value])),
          c(predict(values2),rep(0,length(predict(values2)))),
          col = "aquamarine3",border=NA)
  
  polygon(c(dt$phiB[dt$res == what_resistance_value],rev(dt$phiB[dt$res == what_resistance_value])),
          c(predict(values3),rep(0,length(predict(values3)))),
          col = "orange",border=NA)
  
  polygon(c(dt$phiB[dt$res == what_resistance_value],rev(dt$phiB[dt$res == what_resistance_value])),
          c(predict(values),rep(0,length(predict(values)))),
          col = "blue",border=NA)
  
  abline(v = what_phiB_value,lty=2, lwd=3, col="grey15")
  
  text(0.8,0.1,"Impact from ITN",col="white")
  text(0.55,0.22,"Added impact from LSM",col="white")
  
  
  ####################
  ###
  ####
  ###### Repeat for moderate transmission
  ####
  ###
  ####################
  
  ##########################################
  ##
  ## Figure 1D 
  ##
  ## Show prevalence over time for do nothing, do ITN or do ITN+LSM
  ## 80% RESISTANCE AND 90% INDOOR BITING
  
  ## Moderate transmission first
  aa = readRDS("plots manuscript/modtransmission_newFig1comp.rds")
  aitn = readRDS("plots manuscript/modtransmission_newFig1ITNonly.rds")
  aitnlsm = readRDS("plots manuscript/modtransmission_newFig1ITNandLSM.rds")
  aitnlsmabs = readRDS("plots manuscript/modtransmission_newFig1ITNandLSMAbsolutereduction.rds")
  
  for(i in 1:nrow(dt)){
    dt$prevU5_n[i] = mean(aa[[i]]$pv_0_1825[3285:4380])
    dt$prevU5_itn[i] = mean(aitn[[i]]$pv_0_1825[3285:4380])
    dt$prevU5_itnLSM[i] = mean(aitnlsm[[i]]$pv_0_1825[3285:4380])
    dt$rel_reduction_due_to_itn[i] = (dt$prevU5_n[i] - dt$prevU5_itn[i])/dt$prevU5_n[i]
    dt$rel_reduction_due_to_itnlsm[i] = (dt$prevU5_n[i] - dt$prevU5_itnLSM[i])/dt$prevU5_n[i]
    dt$rel_reduction_if_addlsm_to_itn[i] = (dt$prevU5_itn[i] - dt$prevU5_itnLSM[i])/dt$prevU5_itn[i]
    dt$prevU5_itnLSMabsolute[i] = mean(aitnlsmabs[[i]]$pv_0_1825[3285:4380])
    dt$rel_reduction_due_to_itnlsmabs[i] = (dt$prevU5_n[i] - dt$prevU5_itnLSMabsolute[i])/dt$prevU5_n[i]
    dt$rel_reduction_if_addlsm_to_itnabs[i] = (dt$prevU5_itn[i] - dt$prevU5_itnLSMabsolute[i])/dt$prevU5_itn[i]
  }
  
  summary(dt$prevU5_n)
  summary(dt$rel_reduction_due_to_itn)
  summary(dt$rel_reduction_due_to_itnlsm)
  summary(dt$rel_reduction_due_to_itnlsmabs)
  ddt = subset(dt, dt$res < 0.91)
  range(ddt$rel_reduction_due_to_itnlsm-ddt$rel_reduction_due_to_itn)
  range(ddt$rel_reduction_due_to_itnlsmabs-ddt$rel_reduction_due_to_itn)
  
  dt$rel_reduction_due_to_itn[which(dt$res == what_resistance_value &
                                      dt$phiB == what_phiB_value)]
  dt$rel_reduction_due_to_itn[which_row]
  par(mar = c(4,8,2,1))
  
  plot(aa[[which_row]]$pv_0_1825 ~ aa[[which_row]]$timestep,
       type="l",ylim=c(0,0.8),
       ylab = "falciparum prevalence in U5 years (%)",
       yaxt = "n",xaxt = "n",
       xlim = c(8*365,13*365),
       xlab = "Time in years")
  axis(1, at = c(8,9,10,11,12,13,14,15)*365,
       labels = c(-1,0,1,2,3,4,5,6))
  axis(2, las = 2, at = seq(0,0.8,0.1), labels = seq(0,80,10))
  mtext("Moderate Transmission", side = 2, line = 5, 
        cex = 1)
  
  no_resistance = c(which(dt$res == 0))
  which(dt$phiB == what_phiB_value)
  
  for(i in which_row){
    lines(aitn[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="blue",lwd=3)
    lines(aitnlsm[[i]]$pv_0_1825 ~ aitnlsm[[i]]$timestep,col="aquamarine3",lwd=3)
    lines(aa[[i]]$pv_0_1825 ~ aa[[i]]$timestep,col="grey",lwd=3)
    lines(aitnlsmabs[[i]]$pv_0_1825 ~ aitnlsmabs[[i]]$timestep,col="orange",lwd=3)
  }
  
  ##########################################
  ##
  ## Figure 1E  
  ##
  ## Show relative reduction from ITN and additional LSM impact
  ## against increasing resistance, fixing phiI as 0.9
  par(mar = c(4,4,2,1))
  
  plot(dt$rel_reduction_due_to_itn[dt$phiB == what_phiB_value] ~ 
         dt$res[dt$phiB == what_phiB_value],lwd=3,pch="",
       ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
       xlab = "Resistance (%)",xaxt="n",xlim=c(0,0.9)
  )
  axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
  axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
  values = loess((dt$rel_reduction_due_to_itn[dt$phiB == what_phiB_value] ~ 
                    dt$res[dt$phiB == what_phiB_value])) 
  
  values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$phiB == what_phiB_value] ~ 
                     dt$res[dt$phiB == what_phiB_value])) 
  
  values3 = loess((dt$rel_reduction_due_to_itnlsmabs[dt$phiB == what_phiB_value] ~ 
                     dt$res[dt$phiB == what_phiB_value])) 
  
  polygon(c(dt$res[dt$phiB == what_phiB_value],rev(dt$res[dt$phiB == what_phiB_value])),
          c(predict(values2),rep(0,length(predict(values2)))),
          col = "aquamarine3",border=NA)
  
  polygon(c(dt$res[dt$phiB == what_phiB_value],rev(dt$res[dt$phiB == what_phiB_value])),
          c(predict(values3),rep(0,length(predict(values3)))),
          col = "orange",border=NA)
  
  
  polygon(c(dt$res[dt$phiB == what_phiB_value],rev(dt$res[dt$phiB == what_phiB_value])),
          c(predict(values),rep(0,length(predict(values)))),
          col = "blue",border=NA)
  
  abline(v = what_resistance_value,lty=2, lwd=3, col="grey15")
  
  text(0.2,0.15,"Impact from ITN",col="white")
  text(0.28,0.5,"Added impact from LSM",col="white")
  
  ##########################################
  ##
  ## Figure 1F  
  ##
  ## Show relative reduction from ITN and additional LSM impact
  ## against increasing phiI, resistance at 0.8
  
  plot(dt$rel_reduction_due_to_itn[dt$res == what_resistance_value] ~ 
         dt$phiB[dt$res == what_resistance_value],lwd=3,pch="",
       ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
       xlab = "Proportion of bites in bed (%)",xaxt="n"
  )
  axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
  axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
  values = loess((dt$rel_reduction_due_to_itn[dt$res == what_resistance_value] ~ 
                    dt$phiB[dt$res == what_resistance_value])) 
  
  values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$res == what_resistance_value] ~ 
                     dt$phiB[dt$res == what_resistance_value])) 
  
  values3 = loess((dt$rel_reduction_due_to_itnlsmabs[dt$res == what_resistance_value] ~ 
                     dt$phiB[dt$res == what_resistance_value])) 
  
  polygon(c(dt$phiB[dt$res == what_resistance_value],rev(dt$phiB[dt$res == what_resistance_value])),
          c(predict(values2),rep(0,length(predict(values2)))),
          col = "aquamarine3",border=NA)
  
  polygon(c(dt$phiB[dt$res == what_resistance_value],rev(dt$phiB[dt$res == what_resistance_value])),
          c(predict(values3),rep(0,length(predict(values3)))),
          col = "orange",border=NA)
  
  polygon(c(dt$phiB[dt$res == what_resistance_value],rev(dt$phiB[dt$res == what_resistance_value])),
          c(predict(values),rep(0,length(predict(values)))),
          col = "blue",border=NA)
  
  abline(v = what_phiB_value,lty=2, lwd=3, col="grey15")
  
  text(0.8,0.1,"Impact from ITN",col="white")
  text(0.55,0.22,"Added impact from LSM",col="white")
  
  
  ##########################################
  ##
  ## Figure 1G 
  ##
  ## Show prevalence over time for do nothing, do ITN or do ITN+LSM
  ## 80% RESISTANCE AND 90% INDOOR BITING
  
  ## Low transmission 
  aa = readRDS("plots manuscript/lowtransmission_newFig1comp.rds")
  aitn = readRDS("plots manuscript/lowtransmission_newFig1ITNonly.rds")
  aitnlsm = readRDS("plots manuscript/lowtransmission_newFig1ITNandLSM.rds")
  
  for(i in 1:nrow(dt)){
    dt$prevU5_n[i] = mean(aa[[i]]$pv_0_1825[3285:4380])
    dt$prevU5_itn[i] = mean(aitn[[i]]$pv_0_1825[3285:4380])
    dt$prevU5_itnLSM[i] = mean(aitnlsm[[i]]$pv_0_1825[3285:4380])
    dt$rel_reduction_due_to_itn[i] = (dt$prevU5_n[i] - dt$prevU5_itn[i])/dt$prevU5_n[i]
    dt$rel_reduction_due_to_itnlsm[i] = (dt$prevU5_n[i] - dt$prevU5_itnLSM[i])/dt$prevU5_n[i]
    dt$rel_reduction_if_addlsm_to_itn[i] = (dt$prevU5_itn[i] - dt$prevU5_itnLSM[i])/dt$prevU5_itn[i]
  }
  
  summary(dt$prevU5_n)
  summary(dt$rel_reduction_due_to_itn)
  summary(dt$rel_reduction_due_to_itnlsm)
  ddt = subset(dt, dt$res < 0.91)
  range(ddt$rel_reduction_due_to_itnlsm-ddt$rel_reduction_due_to_itn)
  
  dt$rel_reduction_if_addlsm_to_itn_cum = dt$rel_reduction_if_addlsm_to_itn +dt$rel_reduction_due_to_itn
  
  dt$rel_reduction_due_to_itn[which(dt$res == what_resistance_value &
                                      dt$phiB == what_phiB_value)]
  
  par(mar = c(4,8,2,1))
  
  plot(aa[[which_row]]$pv_0_1825 ~ aa[[which_row]]$timestep,
       type="l",ylim=c(0,0.6),
       ylab = "falciparum prevalence in U5 years (%)",
       yaxt = "n",xaxt = "n",
       xlim = c(8*365,13*365),
       xlab = "Time in years")
  text(9*365,0.78,"80% resistance")
  text(9*365,0.70,"90% bites indoors")
  axis(1, at = c(8,9,10,11,12,13,14,15)*365,
       labels = c(-1,0,1,2,3,4,5,6))
  axis(2, las = 2, at = seq(0,0.6,0.1), labels = seq(0,60,10))
  mtext("Low transmission",side = 2, cex = 1,line = 5)
  no_resistance = c(which(dt$res == 0))
  which(dt$phiB == what_phiB_value)
  
  for(i in which_row){
    lines(aitn[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="blue",lwd=3)
    lines(aitnlsm[[i]]$pv_0_1825 ~ aitnlsm[[i]]$timestep,col="aquamarine3",lwd=3)
    lines(aa[[i]]$pv_0_1825 ~ aa[[i]]$timestep,col="grey",lwd=3)
    lines(aitnlsm[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="orange",lty = 2,lwd=3)
    
  }
  
  ##########################################
  ##
  ## Figure 1H  
  ##
  ## Show relative reduction from ITN and additional LSM impact
  ## against increasing resistance, fixing phiI as 0.9
  par(mar = c(4,4,2,1))
  
  plot(dt$rel_reduction_due_to_itn[dt$phiB == what_phiB_value] ~ 
         dt$res[dt$phiB == what_phiB_value],lwd=3,pch="",
       ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
       xlab = "Resistance (%)",xaxt="n",xlim=c(0,0.9)
  )
  axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
  axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
  values = loess((dt$rel_reduction_due_to_itn[dt$phiB == what_phiB_value] ~ 
                    dt$res[dt$phiB == what_phiB_value])) 
  
  values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$phiB == what_phiB_value] ~ 
                     dt$res[dt$phiB == what_phiB_value])) 
  
  polygon(c(dt$res[dt$phiB == what_phiB_value],rev(dt$res[dt$phiB == what_phiB_value])),
          c(predict(values2),rep(0,length(predict(values2)))),
          col = "orange",border=NA)
  
  polygon(c(dt$res[dt$phiB == what_phiB_value],rev(dt$res[dt$phiB == what_phiB_value])),
          c(predict(values),rep(0,length(predict(values)))),
          col = "blue",border=NA)
  
  text(0.2,0.2,"Impact from ITN",col="white")
  text(0.28,0.70,"Added impact from LSM",col="white")
  
  abline(v=what_resistance_value,col="grey15",lty = 2, lwd = 3)
  ##########################################
  ##
  ## Figure 1I  
  ##
  ## Show relative reduction from ITN and additional LSM impact
  ## against increasing phiI, resistance at 0.8
  
  plot(dt$rel_reduction_due_to_itn[dt$res == what_resistance_value] ~ 
         dt$phiB[dt$res == what_resistance_value],lwd=3,pch="",
       ylim=c(0,0.8),yaxt="n",ylab = "Relative Efficacy (%)",
       xlab = "Proportion of bites in bed (%)",xaxt="n"
  )
  axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
  axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
  values = loess((dt$rel_reduction_due_to_itn[dt$res == 0] ~ 
                    dt$phiB[dt$res == 0])) 
  
  values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$res == what_resistance_value] ~ 
                     dt$phiB[dt$res == what_resistance_value])) 
  
  polygon(c(dt$phiB[dt$res == what_resistance_value],rev(dt$phiB[dt$res == what_resistance_value])),
          c(predict(values2),rep(0,length(predict(values2)))),
          col = "orange",border=NA)
  
  polygon(c(dt$phiB[dt$res == what_resistance_value],rev(dt$phiB[dt$res == what_resistance_value])),
          c(predict(values),rep(0,length(predict(values)))),
          col = "blue",border=NA)
  
  text(0.8,0.2,"Impact from ITN",col="white")
  text(0.55,0.5,"Added impact from LSM",col="white")
  
  abline(v=what_phiB_value,col="grey15",lty = 2, lwd = 3)
  
  par(xpd=NA,cex = 1.11)
  
  text(x = -1.27, y = 3.52,"(A)",cex=0.8)
  text(x = -0.53, y = 3.52,"(B)",cex=0.8)
  text(x = 0.28, y = 3.52,"(C)",cex=0.8)
  text(x = -1.27, y = 2.27,"(D)",cex=0.8)
  text(x = -0.53, y = 2.27,"(E)",cex=0.8)
  text(x = 0.28, y = 2.27,"(F)",cex=0.8)
  text(x = -1.27, y = 0.9,"(G)",cex=0.8)
  text(x = -0.53, y = 0.9,"(H)",cex=0.8)
  text(x = 0.28, y = 0.9,"(I)",cex=0.8)
  
}


draw_comparison_f(what_resistance_value = 0.8,
                  what_phiB_value = 0.872,
                  which_row = 54)

draw_comparison_f(what_resistance_value = 0.0,
                  what_phiB_value = 0.754,
                  which_row = 34)
