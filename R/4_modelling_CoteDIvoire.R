########################################
##
## Cote D'Ivoire trial 
## 2024 publication 
## https://malariajournal.biomedcentral.com/articles/10.1186/s12936-024-04953-8
##
########################################


## We will distribute these counts according to the assumed transmission pattern
## in the model over the year. 


## Difference in differences to estimate impacts from LSM in the trial 

## mosquito counts taken from Additional File 1 means... 
##
## b1 Simulating Cote D'Ivoire trial 
## 2024 publication 
## https://malariajournal.biomedcentral.com/articles/10.1186/s12936-024-04953-8
##
########################################


library(malariasimulation)
library(malariaEquilibrium)
library(reshape2)
library(ggplot2)

# LSM impact (reduction in carrying capacity measured as relative reduction in adult densities per trial arm)

## add prevalence data
# data_people = c(0.82,0.82) ## for each arm for calibration 
# incidence = c(88,136) ## per 1000 U5 per year before
# after_incidence = c(96,50)

## Generate a parameter space across which to sample
dtcd = expand.grid(resistance = seq(0.6,0.95,length=3),
                   phiI = c(0.75,0.8,0.85,0.9))
#60% 0.322297065	0.54341427	1.888350601
#78% 0.274512051	0.509450484	1.685005428
#95% 0.163514499	0.522747255	1.100110971

dtcd$dn0 = ifelse(dtcd$resistance == 0.6,0.322297065,
                  ifelse(dtcd$resistance == 0.775, 0.274512051,0.163514499))

dtcd$rn = ifelse(dtcd$resistance == 0.6,0.54341427,
                  ifelse(dtcd$resistance == 0.775, 0.509450484,0.522747255))

dtcd$gamman = ifelse(dtcd$resistance == 0.6,1.888350601,
                  ifelse(dtcd$resistance == 0.775, 1.685005428,1.100110971))

dtcd$phiB = dtcd$phiI-0.05
####################################################
##
## LSM modelling 
##
####################################################

dat_tia = data.frame(
  village = c("Kolekaha","Lofinekaha","kakologo","Nambatiourkaha"),
  eir_est = c(36,12,100,15),
  gamprop = round(c(0.98,0.953,0.880,0.853),2),
  othprop = round(c(0.02,0.007,0.120,0.026),2),
  lsm_mean = c(0,0,0.9,0.5),
  lsm_low = c(0,0,0.7,-0.6),
  lsm_high = c(0,0,0.97,0.85)
)
dat_tia$funprop = round(1 - dat_tia$gamprop - dat_tia$othprop,2)

################################################
##
## Set up to emulate Kenya trial

lsm_cotedivoire_f = function(datarow,
                             res,
                             lsm_range){
  ## Base settings
  year <- 365
  month <- 30
  sim_length <- 10 * year
  human_population <- 10000
  starting_EIR <- dat_tia$eir_est[datarow]  ## This will change to hit estimated prevalence in U5 yrs 
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      
      prevalence_rendering_min_ages = c(0, 5,  15, 0) * 365, ## Prev in under 5 years measured
      prevalence_rendering_max_ages = c(5, 15, 100,100) * 365,
      
      clinical_incidence_rendering_min_ages = c(0, 5, 15,  0) * 365, ## All age clin_inc
      clinical_incidence_rendering_max_ages = c(5, 15,100,  100) * 365,
      
      ## Northern Cote D'Ivoire (Savannes) estimates, **Might update to reflect trends for mosquito densities
      # seasonal_a0	seasonal_a1	seasonal_b1	seasonal_a2	seasonal_b2	seasonal_a3	seasonal_b3
      # 0.2854889	-0.2685007	-0.09308475	-0.05130384	0.06080074	0.05664451	0.007304213
      
      # ssa0 <- 0.218087682 
      # ssa1 <- -0.005292592  
      # ssb1 <- 0.174739899 
      # ssa2 <- -0.085277539 
      # ssb2 <- -0.082337283  
      # ssa3 <- 0.017356449  
      # ssb3 <- 0.026755121
      # 
      model_seasonality = TRUE, # seasonal model
      g0 = 0.2854889,
      g = c(-0.2685007, -0.05130384, 0.05664451),
      h = c(-0.09308475, 0.06080074, 0.007304213),
      
      individual_mosquitoes = FALSE, ## True by default
      
      bednets = TRUE
      
    )
  )
  
  # set species
  gamb_params <- malariasimulation::gamb_params
  gamb_params$phi_indoors <- dtcd$phiI[res]
  gamb_params$phi_bednets <- dtcd$phiB[res]
  # 
  fun_params <- malariasimulation::fun_params
  fun_params$phi_indoors <- dtcd$phiI[res]
  fun_params$phi_bednets <- dtcd$phiB[res]
  
  # other_params <- malariasimulation::arab_params
  # other_params$Q0 <- 0.93
  # other_params$phi_indoors <- dtcd$phiI[datarow]
  # other_params$phi_bednets <- dtcd$phiB[datarow]
  
  # set mosquito species
  simparams <- set_species(simparams,
                           species=list(gamb_params, arab_params, fun_params),
                           proportions=c(dat_tia$gamprop[datarow],
                                         dat_tia$othprop[datarow],
                                         dat_tia$funprop[datarow])
                           
  )
  
  
  
  # set treatment (using National Values from the USA PMI 2023 Report)
  # 0.3453279	0.1794239 - from our generic estimates
  simparams <- set_drugs(simparams, list(AL_params,    ## whichever is ACT drug
                                         SP_AQ_params))## whichever is non-ACT drug
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,
                                      time=c(1),
                                      coverage=c(0.56))    # currently these are restricted but it would be a future hope to allow them to change
  simparams <- set_clinical_treatment(simparams, 
                                      drug=2,
                                      time=c(1),
                                      coverage=c(0.05))     # currently these are restricted but it would be a future hope to allow them to change
  
  simparams <- set_equilibrium(simparams, starting_EIR)
  
  
  bednetparams <- simparams
  # set bednets
  
  bednet_events = data.frame(
    timestep = c(310,year*3+310,year*7+310), ##start 2014 Nov?, 2017 Nov?, [trial occurs March 2019 - Feb 2020] 2020 Nov... 
    names = c("2014","2017","2021")
  )
  
  ## Bednet 1 corresponds to pyrethroid-only nets
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(0.8,0.8,0.8),
    retention = 3.1 * year,
    ## 80% resistance for pyrethroid-ITN 0.22172514	0.720816403	2.051843084
    dn0 = matrix(dtcd$dn0[res], nrow=3, ncol=3),
    rn =  matrix(dtcd$rn[res], nrow=3, ncol=3),
    rnm = matrix(0.24, nrow=3, ncol=3),
    gamman = rep(dtcd$gamman[res] * 365, 3) 
    ## ** will do a sensitivity analysis
  )
  
  ## LSM
  cc <- get_init_carrying_capacity(bednetparams_1)
  cc
  
  ## circle through lsm_mean, lsm_low, lsm_high
  if (lsm_range == "mean"){
    simparams <- bednetparams_1 |>
      set_carrying_capacity(
        carrying_capacity = matrix(c(1 - dat_tia$lsm_mean[datarow],1,
                                     1 - dat_tia$lsm_mean[datarow],1,
                                     1 - dat_tia$lsm_mean[datarow],1),nrow = 2),
        timesteps = c(6*year + 31+30+1,7*year + 31+30+1) ## Starting March 2019 for 2 years
      )
    
  } else if (lsm_range == "low"){
    simparams <- bednetparams_1 |>
      set_carrying_capacity(
        carrying_capacity = matrix(c(1 - dat_tia$lsm_low[datarow],1,
                                     1 - dat_tia$lsm_low[datarow],1,
                                     1 - dat_tia$lsm_low[datarow],1),nrow = 2),
        timesteps = c(6*year + 31+30+1,7*year + 31+30+1) ## Starting March 2019 for 2 years
      )
    
  } else if (lsm_range == "high"){
    simparams <- bednetparams_1 |>
      set_carrying_capacity(
        carrying_capacity = matrix(c(1 - dat_tia$lsm_high[datarow],1,
                                     1 - dat_tia$lsm_high[datarow],1,
                                     1 - dat_tia$lsm_high[datarow],1),nrow = 2),
        timesteps = c(6*year + 31+30+1,7*year + 31+30+1) ## Starting March 2019 for 2 years
      )
    
  } else {
    simparams <- bednetparams_1 |>
      set_carrying_capacity(
        carrying_capacity = matrix(c(1 - dat_tia$counterf[datarow],1,
                                     1 - dat_tia$counterf[datarow],1,
                                     1 - dat_tia$counterf[datarow],1),nrow = 2),
        timesteps = c(6*year + 31+30+1,7*year + 31+30+1) ## Starting March 2019 for 2 years
      )
  }
  
  
  correlationsb1 <- get_correlation_parameters(simparams)
  correlationsb1$inter_round_rho('bednets', 1)
  
  output_net1  <- run_simulation(sim_length, simparams, correlationsb1) 
  output_net1$pv_all_age = output_net1$n_detect_lm_0_36500/output_net1$n_age_0_36500
  output_net1$INC5 = output_net1$n_inc_clinical_0_1825/output_net1$n_age_0_1825
  output_net1$INC15 = output_net1$n_inc_clinical_1825_5475/output_net1$n_age_1825_5475
  output_net1$INColder = output_net1$n_inc_clinical_5475_36500/output_net1$n_age_5475_36500
  output_net1$INCall = output_net1$n_inc_clinical_0_36500/output_net1$n_age_0_36500
 
  return(data.frame(timestep = output_net1$timestep,
                    net_use = output_net1$n_use_net/human_population,
                    prev_all_age = output_net1$pv_all_age,
                    inc_u5 = output_net1$INC5,
                    inc_5_15 = output_net1$INC15,
                    inc_15up = output_net1$INColder,
                    inc_all_Age = output_net1$INCall,
                    ccgam = output_net1$total_M_gamb,
                    ccfun = output_net1$total_M_fun,
                    ccara = output_net1$total_M_arab
                    
                    ) )
  
  
}
dat_tia$counterf = rep(0, nrow(dat_tia))
v1 = v2 = v3 = v4 = list()
v3L = v4L = v3H = v4H = list()
v3_counterf = v4_counterf = list()

for(i in 1:12){
  v3_counterf[[i]] = lsm_cotedivoire_f(datarow = 3,
                                       res = i,
                                       lsm_range = "counterf")
  v4_counterf[[i]] = lsm_cotedivoire_f(datarow = 4,
                                       res = i,
                                       lsm_range = "counterf")
  
}

for(i in 1:nrow(dtcd)){
  v1[[i]] = lsm_cotedivoire_f(datarow = 1,
                              res = i,
                              lsm_range = "mean") ## MEAN
  v2[[i]] = lsm_cotedivoire_f(datarow = 2,
                        res = i,
                        lsm_range = "mean") ## MIN
  v3[[i]] = lsm_cotedivoire_f(datarow = 3,
                        res = i,
                        lsm_range = "mean") ## MAX
  v4[[i]] = lsm_cotedivoire_f(datarow = 4,
                       res = i,
                       lsm_range = "mean") ## MAX

  v3L[[i]] = lsm_cotedivoire_f(datarow = 3,
                       res = i,
                       lsm_range = "low") ## MAX
  v4L[[i]] = lsm_cotedivoire_f(datarow = 4,
                       res = i,
                       lsm_range = "low") ## MAX
  
  v3H[[i]] = lsm_cotedivoire_f(datarow = 3,
                       res = i,
                       lsm_range = "high") ## MAX
  v4H[[i]] = lsm_cotedivoire_f(datarow = 4,
                       res = i,
                       lsm_range = "high") ## MAX
  
  print(i)
}

## Check: Net use when LSM begins (aiming for around 35%)
timestep = c(6*year + 31+30+1,7*year + 31+30+1)
human_population = 10000
v1[[1]]$net_use[timestep]
v2[[1]]$net_use[timestep]
v3[[1]]$net_use[timestep]
v4[[1]]$net_use[timestep]


uncertainty_absolute_num_mos3 = 
  uncertainty_absolute_num_mos4 = array(data = NA,dim = c(12,9))
for(i in 1:12){
  uncertainty_absolute_num_mos3[i,1] = mean(v3[[i]]$ccgam[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3[[i]]$ccgam[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos3[i,2] = mean(v3L[[i]]$ccgam[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3L[[i]]$ccgam[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos3[i,3] = mean(v3H[[i]]$ccgam[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3H[[i]]$ccgam[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  
  uncertainty_absolute_num_mos3[i,4] = mean(v3[[i]]$ccfun[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3[[i]]$ccfun[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos3[i,5] = mean(v3L[[i]]$ccfun[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3L[[i]]$ccfun[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos3[i,6] = mean(v3H[[i]]$ccfun[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3H[[i]]$ccfun[c(6*year + 31+30+1):c(7*year + 31+30+1)])

  uncertainty_absolute_num_mos3[i,7] = mean(v3[[i]]$ccara[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3[[i]]$ccara[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos3[i,8] = mean(v3L[[i]]$ccara[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3L[[i]]$ccara[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos3[i,9] = mean(v3H[[i]]$ccara[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3H[[i]]$ccara[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  
  uncertainty_absolute_num_mos4[i,1] = mean(v4[[i]]$ccgam[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v4[[i]]$ccgam[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos4[i,2] = mean(v4L[[i]]$ccgam[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v4L[[i]]$ccgam[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos4[i,3] = mean(v4H[[i]]$ccgam[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v4H[[i]]$ccgam[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  
  uncertainty_absolute_num_mos4[i,4] = mean(v4[[i]]$ccfun[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v4[[i]]$ccfun[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos4[i,5] = mean(v4L[[i]]$ccfun[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v4L[[i]]$ccfun[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos4[i,6] = mean(v4H[[i]]$ccfun[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v4H[[i]]$ccfun[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  
  uncertainty_absolute_num_mos4[i,7] = mean(v4[[i]]$ccara[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3[[i]]$ccara[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos4[i,8] = mean(v4L[[i]]$ccara[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3L[[i]]$ccara[c(6*year + 31+30+1):c(7*year + 31+30+1)])
  uncertainty_absolute_num_mos4[i,9] = mean(v4H[[i]]$ccara[c(5*year + 31+30+1):c(6*year + 31+30+1)]) - mean(v3H[[i]]$ccara[c(6*year + 31+30+1):c(7*year + 31+30+1)])
}
##GAM
quantile(uncertainty_absolute_num_mos3[,1:3],c(0.05,0.5,0.95))/10000
quantile(uncertainty_absolute_num_mos4[,1:3],c(0.05,0.5,0.95))/10000

##FUN
quantile(uncertainty_absolute_num_mos3[,4:6],c(0.05,0.5,0.95))/10000
quantile(uncertainty_absolute_num_mos4[,4:6],c(0.05,0.5,0.95))/10000

##OTH
quantile(uncertainty_absolute_num_mos3[,7:9],c(0.05,0.5,0.95))/10000
quantile(uncertainty_absolute_num_mos4[,7:9],c(0.05,0.5,0.95))/10000

#######################
##
## Predictions altogether Village level
##
######################

## Create medians and limits files
##MEDIANS
range_estimates_v3countf = array(data=NA, dim=c(nrow(v3_counterf[[1]]),12) )
range_estimates_v4countf = array(data=NA, dim=c(nrow(v4_counterf[[1]]),12) )
range_estimates_v1 = array(data=NA, dim=c(nrow(v1[[1]]),36) )
range_estimates_v2 = array(data=NA, dim=c(nrow(v2[[1]]),36) )
range_estimates_v3 = array(data=NA, dim=c(nrow(v3[[1]]),36) )
range_estimates_v4 = array(data=NA, dim=c(nrow(v4[[1]]),36) )

range_estINC_v3countf = array(data=NA, dim=c(nrow(v3_counterf[[1]]),12) )
range_estINC_v4countf = array(data=NA, dim=c(nrow(v4_counterf[[1]]),12) )
range_estINC_v1 = array(data=NA, dim=c(nrow(v1[[1]]),36) )
range_estINC_v2 = array(data=NA, dim=c(nrow(v2[[1]]),36) )
range_estINC_v3 = array(data=NA, dim=c(nrow(v3[[1]]),36) )
range_estINC_v4 = array(data=NA, dim=c(nrow(v4[[1]]),36) )

range_estINCu5_v3countf = array(data=NA, dim=c(nrow(v3_counterf[[1]]),12) )
range_estINCu5_v4countf = array(data=NA, dim=c(nrow(v4_counterf[[1]]),12) )
range_estINCu5_v1 = array(data=NA, dim=c(nrow(v1[[1]]),36) )
range_estINCu5_v2 = array(data=NA, dim=c(nrow(v2[[1]]),36) )
range_estINCu5_v3 = array(data=NA, dim=c(nrow(v3[[1]]),36) )
range_estINCu5_v4 = array(data=NA, dim=c(nrow(v4[[1]]),36) )


for(i in 1:12){
  range_estimates_v3countf[,i] = v3_counterf[[i]]$prev_all_age
  range_estimates_v4countf[,i] = v4_counterf[[i]]$prev_all_age
  
  range_estimates_v1[,i] = v1[[i]]$prev_all_age
  range_estimates_v1[,i+12] = v1[[i]]$prev_all_age
  range_estimates_v1[,i+24] = v1[[i]]$prev_all_age
  range_estimates_v2[,i] = v2[[i]]$prev_all_age
  range_estimates_v2[,i+12] = v2[[i]]$prev_all_age
  range_estimates_v2[,i+24] = v2[[i]]$prev_all_age
  range_estimates_v3[,i] = v3[[i]]$prev_all_age
  range_estimates_v3[,i+12] = v3L[[i]]$prev_all_age
  range_estimates_v3[,i+24] = v3H[[i]]$prev_all_age
  range_estimates_v4[,i] = v4[[i]]$prev_all_age
  range_estimates_v4[,i+12] = v4L[[i]]$prev_all_age
  range_estimates_v4[,i+24] = v4H[[i]]$prev_all_age

  range_estINC_v3countf[,i] = v3_counterf[[i]]$inc_all_Age
  range_estINC_v4countf[,i] = v4_counterf[[i]]$inc_all_Age
  
  range_estINC_v1[,i] = v1[[i]]$inc_all_Age
  range_estINC_v1[,i+12] = v1[[i]]$inc_all_Age
  range_estINC_v1[,i+24] = v1[[i]]$inc_all_Age
  range_estINC_v2[,i] = v2[[i]]$inc_all_Age
  range_estINC_v2[,i+12] = v2[[i]]$inc_all_Age
  range_estINC_v2[,i+24] = v2[[i]]$inc_all_Age
  range_estINC_v3[,i] = v3[[i]]$inc_all_Age
  range_estINC_v3[,i+12] = v3L[[i]]$inc_all_Age
  range_estINC_v3[,i+24] = v3H[[i]]$inc_all_Age
  range_estINC_v4[,i] = v4[[i]]$inc_all_Age
  range_estINC_v4[,i+12] = v4L[[i]]$inc_all_Age
  range_estINC_v4[,i+24] = v4H[[i]]$inc_all_Age

  range_estINCu5_v3countf[,i] = v3_counterf[[i]]$inc_u5
  range_estINCu5_v4countf[,i] = v4_counterf[[i]]$inc_u5
  
  range_estINCu5_v1[,i] = v1[[i]]$inc_u5
  range_estINCu5_v1[,i+12] = v1[[i]]$inc_u5
  range_estINCu5_v1[,i+24] = v1[[i]]$inc_u5
  range_estINCu5_v2[,i] = v2[[i]]$inc_u5
  range_estINCu5_v2[,i+12] = v2[[i]]$inc_u5
  range_estINCu5_v2[,i+24] = v2[[i]]$inc_u5
  range_estINCu5_v3[,i] = v3[[i]]$inc_u5
  range_estINCu5_v3[,i+12] = v3L[[i]]$inc_u5
  range_estINCu5_v3[,i+24] = v3H[[i]]$inc_u5
  range_estINCu5_v4[,i] = v4[[i]]$inc_u5
  range_estINCu5_v4[,i+12] = v4L[[i]]$inc_u5
  range_estINCu5_v4[,i+24] = v4H[[i]]$inc_u5
}


range_est1 = expand.grid(timestep = v1[[1]]$timestep)
range_est2 = expand.grid(timestep = v2[[1]]$timestep)
range_est3 = expand.grid(timestep = v3[[1]]$timestep)
range_est4 = expand.grid(timestep = v4[[1]]$timestep)

range_inc1 = expand.grid(timestep = v1[[1]]$timestep)
range_inc2 = expand.grid(timestep = v2[[1]]$timestep)
range_inc3 = expand.grid(timestep = v3[[1]]$timestep)
range_inc4 = expand.grid(timestep = v4[[1]]$timestep)

range_u5_inc1 = expand.grid(timestep = v1[[1]]$timestep)
range_u5_inc2 = expand.grid(timestep = v2[[1]]$timestep)
range_u5_inc3 = expand.grid(timestep = v3[[1]]$timestep)
range_u5_inc4 = expand.grid(timestep = v4[[1]]$timestep)

cntf3prev = expand.grid(timestep = v3[[1]]$timestep)
cntf4prev = expand.grid(timestep = v4[[1]]$timestep)
cntf3inc = expand.grid(timestep = v3[[1]]$timestep)
cntf4inc = expand.grid(timestep = v4[[1]]$timestep)
cntf3incU5 = expand.grid(timestep = v3[[1]]$timestep)
cntf4incU5 = expand.grid(timestep = v4[[1]]$timestep)

range_est1[,2] = rowMeans(range_estimates_v1)
range_est2[,2] = rowMeans(range_estimates_v2)
range_est3[,2] = rowMeans(range_estimates_v3)
range_est4[,2] = rowMeans(range_estimates_v4)
cntf3prev[,2] = rowMeans(range_estimates_v3countf)
cntf4prev[,2] = rowMeans(range_estimates_v4countf)

range_inc1[,2] = rowMeans(range_estINC_v1)
range_inc2[,2] = rowMeans(range_estINC_v2)
range_inc3[,2] = rowMeans(range_estINC_v3)
range_inc4[,2] = rowMeans(range_estINC_v4)
cntf3inc[,2] = rowMeans(range_estINC_v3countf)
cntf4inc[,2] = rowMeans(range_estINC_v4countf)

range_u5_inc1[,2] = rowMeans(range_estINCu5_v1)
range_u5_inc2[,2] = rowMeans(range_estINCu5_v2)
range_u5_inc3[,2] = rowMeans(range_estINCu5_v3)
range_u5_inc4[,2] = rowMeans(range_estINCu5_v4)
cntf3incU5[,2] = rowMeans(range_estINCu5_v3countf)
cntf4incU5[,2] = rowMeans(range_estINCu5_v4countf)


# rangegm_est1 = expand.grid(timestep = ALLv1[[1]]$timestep)
# rangegm_est2 = expand.grid(timestep = ALLv2[[1]]$timestep)
# rangegm_est3 = expand.grid(timestep = ALLv3[[1]]$timestep)
# rangegm_est4 = expand.grid(timestep = ALLv4[[1]]$timestep)
# 
# rangegm_est1[,2] = rowMeans(range_estimates_v1)
# rangegm_est2[,2] = rowMeans(range_estimates_v2)
# rangegm_est3[,2] = rowMeans(range_estimates_v3)
# rangegm_est4[,2] = rowMeans(range_estimates_v4)

for(i in 1:nrow(range_est1)){
  range_est1[i,3] = as.numeric(quantile(range_estimates_v1[i,],0.025))
  range_est1[i,4] = as.numeric(quantile(range_estimates_v1[i,],0.975))
  
  range_est2[i,3] = as.numeric(quantile(range_estimates_v2[i,],0.025))
  range_est2[i,4] = as.numeric(quantile(range_estimates_v2[i,],0.975))
  
  range_est3[i,3] = as.numeric(quantile(range_estimates_v3[i,],0.025))
  range_est3[i,4] = as.numeric(quantile(range_estimates_v3[i,],0.975))
  
  range_est4[i,3] = as.numeric(quantile(range_estimates_v4[i,],0.025))
  range_est4[i,4] = as.numeric(quantile(range_estimates_v4[i,],0.975))

  ## Incidence
  range_inc1[i,3] = as.numeric(quantile(range_estINC_v1[i,],0.025))
  range_inc1[i,4] = as.numeric(quantile(range_estINC_v1[i,],0.975))
  
  range_inc2[i,3] = as.numeric(quantile(range_estINC_v2[i,],0.025))
  range_inc2[i,4] = as.numeric(quantile(range_estINC_v2[i,],0.975))
  
  range_inc3[i,3] = as.numeric(quantile(range_estINC_v3[i,],0.025))
  range_inc3[i,4] = as.numeric(quantile(range_estINC_v3[i,],0.975))
  
  range_inc4[i,3] = as.numeric(quantile(range_estINC_v4[i,],0.025))
  range_inc4[i,4] = as.numeric(quantile(range_estINC_v4[i,],0.975))
  
  ## U5 Inc
  range_u5_inc1[i,3] = as.numeric(quantile(range_estINCu5_v1[i,],0.025))
  range_u5_inc1[i,4] = as.numeric(quantile(range_estINCu5_v1[i,],0.975))
  
  range_u5_inc2[i,3] = as.numeric(quantile(range_estINCu5_v2[i,],0.025))
  range_u5_inc2[i,4] = as.numeric(quantile(range_estINCu5_v2[i,],0.975))
  
  range_u5_inc3[i,3] = as.numeric(quantile(range_estINCu5_v3[i,],0.025))
  range_u5_inc3[i,4] = as.numeric(quantile(range_estINCu5_v3[i,],0.975))
  
  range_u5_inc4[i,3] = as.numeric(quantile(range_estINCu5_v4[i,],0.025))
  range_u5_inc4[i,4] = as.numeric(quantile(range_estINCu5_v4[i,],0.975))
  
  # ## gambiae counts
  # rangegm_est1[i,3] = as.numeric(quantile(mosqGam_estimates_v1[i,],0.025))
  # rangegm_est1[i,4] = as.numeric(quantile(mosqGam_estimates_v1[i,],0.975))
  # 
  # rangegm_est2[i,3] = as.numeric(quantile(mosqGam_estimates_v2[i,],0.025))
  # rangegm_est2[i,4] = as.numeric(quantile(mosqGam_estimates_v2[i,],0.975))
  # 
  # rangegm_est3[i,3] = as.numeric(quantile(mosqGam_estimates_v3[i,],0.025))
  # rangegm_est3[i,4] = as.numeric(quantile(mosqGam_estimates_v3[i,],0.975))
  # 
  # rangegm_est4[i,3] = as.numeric(quantile(mosqGam_estimates_v4[i,],0.025))
  # rangegm_est4[i,4] = as.numeric(quantile(mosqGam_estimates_v4[i,],0.975))
  # 
  
}

write.csv(range_inc1,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage1inc.csv")
write.csv(range_inc2,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage2inc.csv")
write.csv(range_inc3,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage3inc.csv")
write.csv(range_inc4,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage4inc.csv")

write.csv(range_u5_inc1,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage1u5_inc.csv")
write.csv(range_u5_inc2,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage2u5_inc.csv")
write.csv(range_u5_inc3,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage3u5_inc.csv")
write.csv(range_u5_inc4,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage4u5_inc.csv")

write.csv(range_est1,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage1.csv")
write.csv(range_est2,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage2.csv")
write.csv(range_est3,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage3.csv")
write.csv(range_est4,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage4.csv")

write.csv(range_estimates_v1,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage1PREV.csv")
write.csv(range_estimates_v2,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage2PREV.csv")
write.csv(range_estimates_v3,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage3PREV.csv")
write.csv(range_estimates_v4,"C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/COTEDIVOIREvillage4PREV.csv")
#########################################
##
## Model estimated reductions
mean(c(cntf3prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2] - range_est3[c(6*year + 31+30+1):c(7*year + 31+30+1),2])/cntf3prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2])
mean(c(cntf4prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2] - range_est4[c(6*year + 31+30+1):c(7*year + 31+30+1),2])/cntf4prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2])

##
mean(c(cntf3prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2] - range_est3[c(6*year + 31+30+1):c(7*year + 31+30+1),3])/cntf3prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2])
mean(c(cntf4prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2] - range_est4[c(6*year + 31+30+1):c(7*year + 31+30+1),3])/cntf4prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2])

##
mean(c(cntf3prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2] - range_est3[c(6*year + 31+30+1):c(7*year + 31+30+1),4])/cntf3prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2])
mean(c(cntf4prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2] - range_est4[c(6*year + 31+30+1):c(7*year + 31+30+1),4])/cntf4prev[c(6*year + 31+30+1):c(7*year + 31+30+1),2])

##########################################
##
## Read in the data
##
range_est1 = read.csv("C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/cdi/COTEDIVOIREvillage1.csv", header = TRUE)
range_est2 = read.csv("C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/cdi/COTEDIVOIREvillage2.csv", header = TRUE)
range_est3 = read.csv("C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/cdi/COTEDIVOIREvillage3.csv", header = TRUE)
range_est4 = read.csv("C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/cdi/COTEDIVOIREvillage4.csv", header = TRUE)

range_estimates_v1 = read.csv("C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/cdi/COTEDIVOIREvillage1PREV.csv", header = TRUE)
range_estimates_v2 = read.csv("C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/cdi/COTEDIVOIREvillage2PREV.csv", header = TRUE)
range_estimates_v3 = read.csv("C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/cdi/COTEDIVOIREvillage3PREV.csv", header = TRUE)
range_estimates_v4 = read.csv("C:/Users/esherrar/OneDrive - LSTM/ICL transfer documents/OneDrive_1_22-09-2024/Documents/lsm/CRITICAL R CODE/model outputs/cdi/COTEDIVOIREvillage4PREV.csv", header = TRUE)

range_est1 = range_est1[,2:5] 
range_est2 = range_est2[,2:5] 
range_est3 = range_est3[,2:5] 
range_est4 = range_est4[,2:5] 

range_estimates_v1 = range_estimates_v1[,2:37]
range_estimates_v2 = range_estimates_v2[,2:37]
range_estimates_v3 = range_estimates_v3[,2:37]
range_estimates_v4 = range_estimates_v4[,2:37]


########################################
##
## Figure 5
year = 365

par(mfrow= c(2,3))
par(mar=c(4,4,1,1))
plot(range_est3[,2] ~ range_est3[,1],
     main = "Kakologo",
     ylab="Prevalence all-ages (%)", yaxt="n",ylim = c(0,1),
     xlab="",pch="",col="aquamarine2",
     xlim = c(5.5*year,7.5*year+30),
     xaxt = "n")
axis(1, at = c(5.5,6,6.5,7)*year,labels = c("Jul 2019","Jan 2020","Jul 2020","Jan 2021"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 6*year + 31+30+1,lty=2)
# mtext("Larvicide village",)

polygon(c(range_est3[,1],rev(range_est3[,1])),
        c(range_est3[,3],rev(range_est3[,4])),
        col=adegenet::transp("aquamarine2",0.4),border=NA)
lines(cntf3prev[,2] ~ cntf3prev[,1],col="grey",lwd=2)
lines(range_est3[,2] ~ range_est3[,1],col="aquamarine2",lwd=1.5)

## add prev data
## methods to generate these prevalence esimtaes
## are in main manuscript - spread cases across seasonal trend
## then selected the matching prevalence estimate for each village
#Kolékaha and Lofinékaha (control villages), and Kakologo and Nambatiourkaha (treatment villages) respectively 
prev_estimates = c(0.34, 0.24, 0.42, 0.27)
prev1 = prev_estimates[3]
inc1 = 339/758
inc2 = 106/758
xx = inc2 * prev1/inc1
SITES3_binomial_mean = c(prev_estimates[3],xx)
SITES3_binomial_low = SITES3_binomial_upp = numeric(2)
SITES3_binomial_low[1] = as.numeric(c(prop.test(x=round(100*prev_estimates[3]), n=100, conf.level=.95, correct=FALSE)$conf.int))[1]
SITES3_binomial_upp[1] = as.numeric(c(prop.test(x=round(100*prev_estimates[3]), n=100, conf.level=.95, correct=FALSE)$conf.int))[2]
SITES3_binomial_low[2] = as.numeric(c(prop.test(x=round(100*xx), n=100, conf.level=.95, correct=FALSE)$conf.int))[1]
SITES3_binomial_upp[2] = as.numeric(c(prop.test(x=round(100*xx), n=100, conf.level=.95, correct=FALSE)$conf.int))[2]

segments(x0 = 2070,
         x1 = 2070,
         y0 = SITES3_binomial_low[1],
         y1 = SITES3_binomial_upp[1], 
         lwd=1.5,
         col="aquamarine2")
segments(x0 = 2070+365,
         x1 = 2070+365,
         y0 = SITES3_binomial_low[2],
         y1 = SITES3_binomial_upp[2], 
         lwd=1.5,
         col="aquamarine2")

points(c(prev_estimates[3],xx) ~ c(2070,2070+365), 
       col="aquamarine2",pch = 19)

##B
plot(range_est4[,2] ~ range_est4[,1],
     main = "Nambatiourkaha",
     ylab="Prevalence all-ages (%)", yaxt="n",ylim = c(0,1),
     xlab="",pch="",col="purple",
     xlim = c(5.5*year,7.5*year+30),
     xaxt = "n")
axis(1, at = c(5.5,6,6.5,7)*year,labels = c("Jul 2019","Jan 2020","Jul 2020","Jan 2021"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 6*year + 31+30+1,lty=2)
# mtext("Larvicide village",)

polygon(c(range_est4[,1],rev(range_est4[,1])),
        c(range_est4[,3],rev(range_est4[,4])),
        col=adegenet::transp("purple",0.4),border=NA)
lines(cntf4prev[,2] ~ cntf4prev[,1],col="grey",lwd=2)
lines(range_est4[,2] ~ range_est4[,1],col="purple",pch=1.5)

## add prev data
prev1 = prev_estimates[4]
inc1 = 339/758
inc2 = 106/758
xx = inc2 * prev1/inc1
SITES4_binomial_mean = c(prev_estimates[4],xx)
SITES4_binomial_low = SITES4_binomial_upp = numeric(2)
SITES4_binomial_low[1] = as.numeric(c(prop.test(x=round(100*prev_estimates[4]), n=100, conf.level=.95, correct=FALSE)$conf.int))[1]
SITES4_binomial_upp[1] = as.numeric(c(prop.test(x=round(100*prev_estimates[4]), n=100, conf.level=.95, correct=FALSE)$conf.int))[2]
SITES4_binomial_low[2] = as.numeric(c(prop.test(x=round(100*xx), n=100, conf.level=.95, correct=FALSE)$conf.int))[1]
SITES4_binomial_upp[2] = as.numeric(c(prop.test(x=round(100*xx), n=100, conf.level=.95, correct=FALSE)$conf.int))[2]

segments(x0 = 2070,
         x1 = 2070,
         y0 = SITES4_binomial_low[1],
         y1 = SITES4_binomial_upp[1], 
         lwd=1.5,
         col="purple")
segments(x0 = 2070+365,
         x1 = 2070+365,
         y0 = SITES4_binomial_low[2],
         y1 = SITES4_binomial_upp[2], 
         lwd=1.5,
         col="purple")

points(c(prev_estimates[4],xx) ~ c(2070,2070+365), 
       col="purple",pch = 19)


## C
##############
##
## bar graph

bars = ## Absolute reduction in mean incidence per 1,000 people per year estimates
  c(sum(range_inc3[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),2]),
    sum(range_inc4[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),2]),
    
    ## Relative reduction in mean incidence estimates
    (sum(range_inc3[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),2]))/sum(cntf3inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]),
    (sum(range_inc4[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),2]))/sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]))
    
sum(range_inc3[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),2])
sum(range_inc3[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),3])
sum(range_inc3[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),4])

100 * (sum(range_inc3[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),2]))/sum(cntf3inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2])
100 * (sum(range_inc3[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),3]))/sum(cntf3inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2])
100 * (sum(range_inc3[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),4]))/sum(cntf3inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2])

## Nam
sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),2])
sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),3])
sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),4])

100 * (sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),2]))/sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2])
100 * (sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),3]))/sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2])
100 * (sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),4]))/sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2])



(291 - 111)/291

##############################
##
############################
barsU = ## Absolute reduction in mean incidence per 1,000 people per year estimates
  c(sum(cntf3inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),3]),
    sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),3]),
    
    ## Relative reduction in mean incidence estimates
    (sum(cntf3inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),3]))/sum(cntf3inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]),
    (sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),3]))/sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2])
  )

barsL = ## Absolute reduction in mean incidence per 1,000 people per year estimates
  c(sum(cntf3inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),4]),
    sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),4]),

    ## Relative reduction in mean incidence estimates
    (sum(cntf3inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc3[c(6*year + 31+30+1):c(7*year + 31+30+1),4]))/sum(cntf3inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]),
    (sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2]) - sum(range_inc4[c(6*year + 31+30+1):c(7*year + 31+30+1),4]))/sum(cntf4inc[c(5*year + 31+30+1):c(6*year + 31+30+1),2])
  )
mat <- matrix(bars, 2, 2)
colnames(mat) <- c("Absolute estimated all-age cases averted",
                   "% reduction incidence (all-age)")
rownames(mat) <- c("Kakologo","Nambatiourkaha")


# plotting settings -------------------------------------------------------
ylim <- c(0,max(c(barsU,barsL,bars)))
angle1 <- rep(c(45,45,135), length.out=2)
angle2 <- rep(c(45,135,135), length.out=2)
density1 <- seq(5,25,length.out=2)
density2 <- seq(5,25,length.out=2)
col <- rep(c("aquamarine2","purple"),2) # rainbow(7)


# plot --------------------------------------------------------------------
par(mar = c(4,5,2,5))
barplot(mat, beside=TRUE, ylim=ylim,
        col=col,angle=angle1,
        density=density1,
        xaxt="n",
        ylab = "",yaxt = "n")
barplot(mat, add=TRUE, beside=TRUE, 
        ylim=ylim, col=col,
        xaxt="n",
        ylab = "",yaxt = "n",
        angle=angle2, density=density2)

x_levels = c(1.5, 2.5,4.5,5.5)
for(i in 1:4){
  segments(x0 = x_levels[i],
           x1 = x_levels[i],
           y0 = barsU[i],
           y1 = barsL[i])
}

axis(2,las = 2, at = seq(0,1,0.2),labels = seq(0,1000,200))
axis(4,las = 2, at = seq(0,1,0.2), labels = seq(0,100,20))
axis(1, at = c(2,5),labels = c("Case counts","Relative reduction"))
mtext("Absolute reduction in",side = 2, line = 4, cex = 0.75) 
mtext("cases averted per 1,000 people",side = 2, line = 3, cex = 0.75)
mtext("Relative reduction in cases (%)",side = 4, line = 2.5,cex = 0.75)
abline(v=6.5,lty=2)
legend("top",legend = c("Kakologo","Nambatiourkaha"),
       col = c("aquamarine2","purple"),pch=15,bty="n")

### D
par(mar=c(4,4,1,1))
plot(range_est1[,2] ~ range_est1[,1],
     main = "Kolekaha",
     ylab="Prevalence all-ages (%)", yaxt="n",ylim = c(0,1),
     xlab="",pch="",col="lightblue",
     xlim = c(5.5*year,7.5*year+30),
     xaxt = "n")
axis(1, at = c(5.5,6,6.5,7)*year,labels = c("Jul 2019","Jan 2020","Jul 2020","Jan 2021"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 6*year + 31+30+1,lty=2)
# mtext("Larvicide village",)

polygon(c(range_est1[,1],rev(range_est1[,1])),
        c(range_est1[,3],rev(range_est1[,4])),
        col=adegenet::transp("lightblue",0.4),border=NA)
lines(range_est1[,2] ~ range_est1[,1],col="lightblue",pch=1.5)

## add prev data
prev1 = prev_estimates[1]
inc1 = 278/496
inc2 = 321/496
xx = inc2 * prev1/inc1
SITES1_binomial_mean = c(prev_estimates[1],xx)
SITES1_binomial_low = SITES1_binomial_upp = numeric(2)
SITES1_binomial_low[1] = as.numeric(c(prop.test(x=round(100*prev_estimates[1]), n=100, conf.level=.95, correct=FALSE)$conf.int))[1]
SITES1_binomial_upp[1] = as.numeric(c(prop.test(x=round(100*prev_estimates[1]), n=100, conf.level=.95, correct=FALSE)$conf.int))[2]
SITES1_binomial_low[2] = as.numeric(c(prop.test(x=round(100*xx), n=100, conf.level=.95, correct=FALSE)$conf.int))[1]
SITES1_binomial_upp[2] = as.numeric(c(prop.test(x=round(100*xx), n=100, conf.level=.95, correct=FALSE)$conf.int))[2]

segments(x0 = 2070,
         x1 = 2070,
         y0 = SITES1_binomial_low[1],
         y1 = SITES1_binomial_upp[1], 
         lwd=1.5,
         col="lightblue")
segments(x0 = 2070+365,
         x1 = 2070+365,
         y0 = SITES1_binomial_low[2],
         y1 = SITES1_binomial_upp[2], 
         lwd=1.5,
         col="lightblue")

points(c(prev_estimates[1],xx) ~ c(2070,2070+365), 
       col="lightblue",pch = 19)

##
par(mar=c(4,4,1,1))
plot(range_est2[,2] ~ range_est2[,1],
     main = "Lofinekaha",
     ylab="Prevalence all-ages (%)", yaxt="n",ylim = c(0,1),
     xlab="",pch="",col="blue",
     xlim = c(5.5*year,7.5*year+30),
     xaxt = "n")
axis(1, at = c(5.5,6,6.5,7)*year,labels = c("Jul 2019","Jan 2020","Jul 2020","Jan 2021"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 6*year + 31+30+1,lty=2)
# mtext("Larvicide village",)

polygon(c(range_est2[,1],rev(range_est2[,1])),
        c(range_est2[,3],rev(range_est2[,4])),
        col=adegenet::transp("blue",0.4),border=NA)
lines(range_est2[,2] ~ range_est2[,1],col="blue",pch=1.5)

## add prev data
prev1 = prev_estimates[2]
inc1 = 215/553
inc2 = 167/553
xx = inc2 * prev1/inc1
SITES2_binomial_mean = c(prev_estimates[2],xx)
SITES2_binomial_low = SITES2_binomial_upp = numeric(2)
SITES2_binomial_low[1] = as.numeric(c(prop.test(x=round(100*prev_estimates[2]), n=100, conf.level=.95, correct=FALSE)$conf.int))[1]
SITES2_binomial_upp[1] = as.numeric(c(prop.test(x=round(100*prev_estimates[2]), n=100, conf.level=.95, correct=FALSE)$conf.int))[2]
SITES2_binomial_low[2] = as.numeric(c(prop.test(x=round(100*xx), n=100, conf.level=.95, correct=FALSE)$conf.int))[1]
SITES2_binomial_upp[2] = as.numeric(c(prop.test(x=round(100*xx), n=100, conf.level=.95, correct=FALSE)$conf.int))[2]

segments(x0 = 2070,
         x1 = 2070,
         y0 = SITES2_binomial_low[1],
         y1 = SITES2_binomial_upp[1], 
         lwd=1.5,
         col="blue")
segments(x0 = 2070+365,
         x1 = 2070+365,
         y0 = SITES2_binomial_low[2],
         y1 = SITES2_binomial_upp[2], 
         lwd=1.5,
         col="blue")

points(c(prev_estimates[2],xx) ~ c(2070,2070+365), 
       col="blue",pch = 19)

#######################################################
##
## Compare the abs versus rel reductions
par(mar = c(5,5,2,5))
## F
vils = c("Musilongo","Musilongo","Kezege","Kezege","Wamondo","Wamondo",
         "Kakologo","Nambatiourkaha")
species = c("gam", "fun", "gam", "fun", "gam", "fun",
            "gam", "gam")
abs_mn = c(41.8, 13.9, 17.2, 9.5, 1.2, 0.9,
           47.2, 3.9)
abs_low = c(39.6, 6.1, 11.1, 0.5, -0.7, 0.7,
            33.2, -5.3)
abs_upp = c(43.3, 15.5, 19.9, 12.7, 1.9, 1.1,
            51.9, 6.8)

rel_mn = c(0.97, 0.87, 0.78, 0.62, 0.47, 0.13,
           0.9, 0.5)
rel_low = c(0.92, 0.25, 0.46, -0.04, -0.42, 0.07,
            0.7, -0.6)
rel_upp = c(0.99, 0.99, 0.92, 0.88, 0.82, 0.21,
            0.97, 0.85)
cols = c("orange","orange","red","red","darkred","darkred","aquamarine2","purple")
spes = c(19,17,19,17,19,17,19,19)
plot(rel_mn ~ abs_mn, col = cols, pch = spes,
     ylab = "", # Relative reduction: vector densities (%)
     xlab = "", # Absolute reduction: vectors per person
     ylim = c(min(rel_low),1),
     yaxt = "n",cex = 2,
     xlim = c(min(abs_low),max(abs_upp)))
axis(2, las = 2, at = c(-0.5,0,0.5,1),labels = c(-50,0,50,100))
mtext("Model derived relative",side = 2, line = 3,cex = 0.8)
mtext("reduction in vector densities (%)",side = 2, line = 2.25,cex = 0.8)

mtext("Model derived absolute",side = 1, line = 2.25,cex = 0.8)
mtext("reduction in vector densities (%)",side = 1, line = 3,cex = 0.8)

for(i in 1:8){
  segments(x0 = abs_upp[i],
           x1 = abs_low[i],
           y0 = rel_mn[i],
           y1 = rel_mn[i],
           col = cols[i],lwd=1)  
  
  segments(x0 = abs_mn[i],
           x1 = abs_mn[i],
           y0 = rel_low[i],
           y1 = rel_upp[i],
           col = cols[i],lwd=1)  
}

legend("right",legend = c("Musilongo","Kezege","Wamondo","Kakologo","Nambatiourkaha"),
       col = c("orange","red","darkred","aquamarine2","purple"),
       pch = 15,bty="n")

legend("bottomright",legend = c("An. gambiae","An. funestus"),
       pch = c(19,17),
       bty="n")

par(xpd=NA,cex = 1.11)

text(x = -225, y = 3.7,"(A)",cex=0.8)
text(x = -128, y = 3.7,"(B)",cex=0.8)
text(x = -20, y = 3.7,"(C)",cex=0.8)
text(x = -225, y = 1.5,"(D)",cex=0.8)
text(x = -128, y = 1.5,"(E)",cex=0.8)
text(x = -20, y = 1.5,"(F)",cex=0.8)
