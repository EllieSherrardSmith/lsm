######################################
##
## RCT Western Kenya, LSM trial 
##
######################################

library(malariasimulation)
library(malariaEquilibrium)
library(reshape2)
library(ggplot2)


## Village specific data
# Proportion An. gambiae s.l.	73.4%	79.7%	59.5%	38.7%	66.3%	76.4%
# Proportion An. funestus s.l.	26.6%	20.3%	40.5%	61.3%	33.7%	23.6%
gamb_prop = c(0.734,0.797,0.595,0.387,0.663,0.764)
fun_prop =  c(0.266,0.203,0.405,0.613,0.337,0.236)

# Historic net use	12.3%	8.1%	20.4%	21.7%	15.8%	9.4%
# July 2005 - 2006	37.8%	15.8%	30.4%	35.4%	25.7%	27.8%
# July 2006 – July 2007	54.3%	36.9%	49.8%	48.4%	43.5%	43.5%
# hist_nets = c(0.123,0.081,0.204,0.217,0.158,0.094)  
# july_2005_nets = c(0.378,0.158,0.304,0.354,0.257,0.278)
# july_2006_nets = c(0.543,0.369,0.498,0.484,0.435,0.435)
nets = read.csv("data/net_use.csv",header=TRUE)

## add prevalence data
data_people = read.csv("data/Fillinger_Kenya highland LSM_MALARIA SURVEY DATA.csv",header=TRUE) 


sites = 1:6
months = c("January","April", "May","June","July","August","November","December")
years = c("2004","2005","2006","2007")

Count_prev_sites = Count_tote_sites = array(dim=c(8,4,6))

for(s in 1:length(sites)){
  for(m in 1:length(months)){
    for(y in 1:length(years)){
      Count_prev_sites[m,y,s] = length(data_people$Microscopy.result.negative..0..or.positive..1.[data_people$Site.code == sites[s] &
                                                                                                    data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                                                                    data_people$Microscopy.result.negative..0..or.positive..1. == 1])      
      Count_tote_sites[m,y,s] = length(data_people$Microscopy.result.negative..0..or.positive..1.[data_people$Site.code == sites[s] &
                                                                                                    data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                                                                    data_people$Microscopy.result.negative..0..or.positive..1. == 0])      
    }
  }
}

Prevalence = ifelse((Count_prev_sites+Count_tote_sites) == 0, -999,Count_prev_sites / (Count_prev_sites+Count_tote_sites))
Prevalence[Prevalence == -999] <- NA

rownames(Prevalence) = c("January","April", "May","June","July","August","November","December")
colnames(Prevalence) = c("2004","2005","2006","2007")


sites1 = as.numeric(c(Prevalence[1,1,1],NA,NA,Prevalence[2:6,1,1],NA,NA,Prevalence[7:8,1,1],
                      Prevalence[1,2,1],NA,NA,Prevalence[2:6,2,1],NA,NA,Prevalence[7:8,2,1],
                      Prevalence[1,3,1],NA,NA,Prevalence[2:6,3,1],NA,NA,Prevalence[7:8,3,1],
                      Prevalence[1,4,1],NA,NA,Prevalence[2:6,4,1],NA,NA,Prevalence[7:8,4,1]))
sites2 = as.numeric(c(Prevalence[1,1,2],NA,NA,Prevalence[2:6,1,2],NA,NA,Prevalence[7:8,1,2],
                      Prevalence[1,2,2],NA,NA,Prevalence[2:6,2,2],NA,NA,Prevalence[7:8,2,2],
                      Prevalence[1,3,2],NA,NA,Prevalence[2:6,3,2],NA,NA,Prevalence[7:8,3,2],
                      Prevalence[1,4,2],NA,NA,Prevalence[2:6,4,2],NA,NA,Prevalence[7:8,4,2]))
sites3 = as.numeric(c(Prevalence[1,1,3],NA,NA,Prevalence[2:6,1,3],NA,NA,Prevalence[7:8,1,3],
                      Prevalence[1,2,3],NA,NA,Prevalence[2:6,2,3],NA,NA,Prevalence[7:8,2,3],
                      Prevalence[1,3,3],NA,NA,Prevalence[2:6,3,3],NA,NA,Prevalence[7:8,3,3],
                      Prevalence[1,4,3],NA,NA,Prevalence[2:6,4,3],NA,NA,Prevalence[7:8,4,3]))
sites4 = as.numeric(c(Prevalence[1,1,4],NA,NA,Prevalence[2:6,1,4],NA,NA,Prevalence[7:8,1,4],
                      Prevalence[1,2,4],NA,NA,Prevalence[2:6,2,4],NA,NA,Prevalence[7:8,2,4],
                      Prevalence[1,3,4],NA,NA,Prevalence[2:6,3,4],NA,NA,Prevalence[7:8,3,4],
                      Prevalence[1,4,4],NA,NA,Prevalence[2:6,4,4],NA,NA,Prevalence[7:8,4,4]))
sites5 = as.numeric(c(Prevalence[1,1,5],NA,NA,Prevalence[2:6,1,5],NA,NA,Prevalence[7:8,1,5],
                      Prevalence[1,2,5],NA,NA,Prevalence[2:6,2,5],NA,NA,Prevalence[7:8,2,5],
                      Prevalence[1,3,5],NA,NA,Prevalence[2:6,3,5],NA,NA,Prevalence[7:8,3,5],
                      Prevalence[1,4,5],NA,NA,Prevalence[2:6,4,5],NA,NA,Prevalence[7:8,4,5]))
sites6 = as.numeric(c(Prevalence[1,1,6],NA,NA,Prevalence[2:6,1,6],NA,NA,Prevalence[7:8,1,6],
                      Prevalence[1,2,6],NA,NA,Prevalence[2:6,2,6],NA,NA,Prevalence[7:8,2,6],
                      Prevalence[1,3,6],NA,NA,Prevalence[2:6,3,6],NA,NA,Prevalence[7:8,3,6],
                      Prevalence[1,4,6],NA,NA,Prevalence[2:6,4,6],NA,NA,Prevalence[7:8,4,6]))

Total_tests = Count_prev_sites+Count_tote_sites
numsites1 = as.numeric(c(Total_tests[1,1,1],NA,NA,Total_tests[2:6,1,1],NA,NA,Total_tests[7:8,1,1],
                         Total_tests[1,2,1],NA,NA,Total_tests[2:6,2,1],NA,NA,Total_tests[7:8,2,1],
                         Total_tests[1,3,1],NA,NA,Total_tests[2:6,3,1],NA,NA,Total_tests[7:8,3,1],
                         Total_tests[1,4,1],NA,NA,Total_tests[2:6,4,1],NA,NA,Total_tests[7:8,4,1]))
numsites2 = as.numeric(c(Total_tests[1,1,2],NA,NA,Total_tests[2:6,1,2],NA,NA,Total_tests[7:8,1,2],
                         Total_tests[1,2,2],NA,NA,Total_tests[2:6,2,2],NA,NA,Total_tests[7:8,2,2],
                         Total_tests[1,3,2],NA,NA,Total_tests[2:6,3,2],NA,NA,Total_tests[7:8,3,2],
                         Total_tests[1,4,2],NA,NA,Total_tests[2:6,4,2],NA,NA,Total_tests[7:8,4,2]))
numsites3 = as.numeric(c(Total_tests[1,1,3],NA,NA,Total_tests[2:6,1,3],NA,NA,Total_tests[7:8,1,3],
                         Total_tests[1,2,3],NA,NA,Total_tests[2:6,2,3],NA,NA,Total_tests[7:8,2,3],
                         Total_tests[1,3,3],NA,NA,Total_tests[2:6,3,3],NA,NA,Total_tests[7:8,3,3],
                         Total_tests[1,4,3],NA,NA,Total_tests[2:6,4,3],NA,NA,Total_tests[7:8,4,3]))
numsites4 = as.numeric(c(Total_tests[1,1,4],NA,NA,Total_tests[2:6,1,4],NA,NA,Total_tests[7:8,1,4],
                         Total_tests[1,2,4],NA,NA,Total_tests[2:6,2,4],NA,NA,Total_tests[7:8,2,4],
                         Total_tests[1,3,4],NA,NA,Total_tests[2:6,3,4],NA,NA,Total_tests[7:8,3,4],
                         Total_tests[1,4,4],NA,NA,Total_tests[2:6,4,4],NA,NA,Total_tests[7:8,4,4]))
numsites5 = as.numeric(c(Total_tests[1,1,5],NA,NA,Total_tests[2:6,1,5],NA,NA,Total_tests[7:8,1,5],
                         Total_tests[1,2,5],NA,NA,Total_tests[2:6,2,5],NA,NA,Total_tests[7:8,2,5],
                         Total_tests[1,3,5],NA,NA,Total_tests[2:6,3,5],NA,NA,Total_tests[7:8,3,5],
                         Total_tests[1,4,5],NA,NA,Total_tests[2:6,4,5],NA,NA,Total_tests[7:8,4,5]))
numsites6 = as.numeric(c(Total_tests[1,1,6],NA,NA,Total_tests[2:6,1,6],NA,NA,Total_tests[7:8,1,6],
                         Total_tests[1,2,6],NA,NA,Total_tests[2:6,2,6],NA,NA,Total_tests[7:8,2,6],
                         Total_tests[1,3,6],NA,NA,Total_tests[2:6,3,6],NA,NA,Total_tests[7:8,3,6],
                         Total_tests[1,4,6],NA,NA,Total_tests[2:6,4,6],NA,NA,Total_tests[7:8,4,6]))

timeline = 1:48 ## JAN 2004 UNTIL DEC 2007

year <- 365

align_timeline = round(seq(2*year+1,2*year+1+48 * 365/12,length=48),0)

## input params
parms = read.csv("R/malariasimulation/params_lsm_range.csv",header=TRUE)

## LSM modelling 

year <- 365
month <- 30
sim_length <- 10 * year
human_population <- 10000
starting_EIR <- 50  ## This will change to hit estimated prevalence in U5 yrs of (5%, 10%, 20%, 30%, 40%, 50%, and 60%)

###########################################
##
## function for running the model 

lsm_LOOP_f = function(parms_row,village,nets_use,root_eir){
  ## LSM modelling 
  
  ## LSM modelling 
  
  year <- 365
  month <- 30
  sim_length <- 10 * year
  human_population <- 10000
  starting_EIR <- root_eir  ## This will change to hit estimated prevalence in U5 yrs of (5%, 10%, 20%, 30%, 40%, 50%, and 60%)
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      
      prevalence_rendering_min_ages = c( 0.5,0, 5,  15, 0) * 365, ## Prev in under 5 years measured
      prevalence_rendering_max_ages = c(10,  5, 15, 100,100) * 365,
      
      clinical_incidence_rendering_min_ages = c(0.5,  0) * 365, ## All age clin_inc
      clinical_incidence_rendering_max_ages = c(5,  100) * 365,
      
      ## Western Kenya estimates, Fitted to RCT trends for mosquito densities
      # ssa0 <- 0.218087682 
      # ssa1 <- -0.005292592  
      # ssb1 <- 0.174739899 
      # ssa2 <- -0.085277539 
      # ssb2 <- -0.082337283  
      # ssa3 <- 0.017356449  
      # ssb3 <- 0.026755121
      # 
      model_seasonality = TRUE, # seasonal model
      g0 = 0.218087682,
      g = c(-0.005292592, -0.085277539, 0.017356449),
      h = c(0.174739899, -0.082337283, 0.026755121),
      
      # model_seasonality = TRUE, # perrennial model
      # g0 = 0.285277,
      # g = c(-0.0248801, -0.0529426, -0.016891),
      # h = c(-0.0216681, -0.0242904, -0.00073646),
      
      individual_mosquitoes = FALSE, ## True by default
      
      bednets = TRUE
      
    )
  )
  
  simparams <- set_equilibrium(simparams, starting_EIR)
  
  # set treatment
  # 0.3453279	0.1794239 - from our generic estimates
  simparams <- set_drugs(simparams, list(AL_params,    ## whichever is ACT drug
                                         SP_AQ_params))## whichever is non-ACT drug
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,
                                      time=c(100,2*year+183),
                                      coverage=c(0.05,0.4))    # currently these are restricted but it would be a future hope to allow them to change
  simparams <- set_clinical_treatment(simparams, 
                                      drug=2,
                                      time=c(100,2*year+183),
                                      coverage=c(0.35,0.1))     # currently these are restricted but it would be a future hope to allow them to change
  
  
  
  
  
  
  
  # set mosquito species
  gam_params <- gamb_params # gambiae
  fun_params <- fun_params  # funestus
  
  gam_params['species'] <- "gambiae sl"
  gam_params['blood_meal_rates'] <- parms$gam_blood_meal_rates[parms_row] # 1/duration of gonothropic cycle
  gam_params['foraging_time'] <- parms$gam_foraging_time[parms_row] # time spent foraging
  gam_params['Q0'] <- parms$gam_Q0[parms_row] # human blood index: update from Nilani (sensitivity analysis?)
  gam_params['phi_bednets'] <- parms$gam_phiB[parms_row] # proportion biting in bed: update from Nilani (sensitivity analysis?)
  gam_params['phi_indoors'] <- parms$gam_phiI[parms_row] # proportion biting indoors: update from Nilani (sensitivity analysis?)
  gam_params['mum'] <- parms$gam_mum[parms_row] # death rate or 1/life expectancy: update from Nilani (sensitivity analysis?)
  
  fun_params['species'] <- "funestus sl"
  fun_params['blood_meal_rates'] <- parms$fun_blood_meal_rates[parms_row] # 1/duration of gonothropic cycle
  fun_params['foraging_time'] <- parms$fun_foraging_time[parms_row] # time spent foraging
  fun_params['Q0'] <- parms$fun_Q0[parms_row] # human blood index: update from Nilani (sensitivity analysis?)
  fun_params['phi_bednets'] <- parms$fun_phiB[parms_row] # proportion biting in bed: update from Nilani (sensitivity analysis?)
  fun_params['phi_indoors'] <- parms$fun_phiI[parms_row] # proportion biting indoors: update from Nilani (sensitivity analysis?)
  fun_params['mum'] <- parms$fun_mum[parms_row] # death rate or 1/life expectancy: update from Nilani (sensitivity analysis?)
  
  
  simparams <- set_species(simparams,
                           species=list(gam_params, fun_params),
                           proportions=c(gamb_prop[village],fun_prop[village])
                           
  )
  
  
  bednetparams <- simparams
  # set bednets
  bednet_events = data.frame(
    timestep = nets_use$timesteps,
    names = rep("use",length(nets_use$timesteps))
  )
  
  ## Bednet 1 corresponds to pyrethroid-only nets no resistance
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    coverages = nets_use$nets,
    
    retention = 5 * year, 
    
    dn0 = matrix(parms$gam_dn0[parms_row], nrow=nrow(nets_use), ncol=2),
    rn =  matrix(parms$gam_rn[parms_row], nrow=nrow(nets_use), ncol=2),
    rnm = matrix(0.24, nrow=nrow(nets_use), ncol=2),
    gamman = rep(parms$gam_gamman[parms_row] * 365, nrow(nets_use)) 
    
  )
  
  # ## LSM
  # lsm_events = data.frame(
  #   habitat_management_timesteps = 3 * year+190, ## july 2005 
  #   name=c("trial-starts")
  # )
  # 
  # habitat_management <- set_habitat_management(
  #   bednetparams_1,
  #   habitat_management_timesteps = lsm_events$habitat_management_timesteps,
  #   lsm_new_eqm = if(village == 1) 
  # matrix(c(parms$gam_rho_1[parms_row],
  #          parms$fun_rho_1[parms_row]),nrow=1,ncol=2) else if(village == 2)
  #            matrix(c(parms$gam_rho_2[parms_row],
  #                     parms$fun_rho_2[parms_row]),nrow=1,ncol=2) else if(village == 3)
  #                       matrix(c(parms$gam_rho_3[parms_row],
  #                                parms$fun_rho_3[parms_row]),nrow=1,ncol=2) else if (village == 4)
  #                                  matrix(c(parms$gam_rho_4[parms_row],
  #                                           parms$fun_rho_4[parms_row]),nrow=1,ncol=2) else if (village == 5)
  #                                             matrix(c(parms$gam_rho_5[parms_row],
  #                                                      parms$fun_rho_5[parms_row]),nrow=1,ncol=2) else if (village == 6)
  #                                                matrix(c(parms$gam_rho_6[parms_row],
  #                                                         parms$fun_rho_6[parms_row]),nrow=1,ncol=2),
  #   lsm_rate_alpha = matrix(-4,nrow=1,ncol=2), # exp(4)
  #   lsm_rate_beta = matrix(0.057,nrow=1,ncol=2)
  # )
  
  ## Specify the LSM coverage
  cc <- get_init_carrying_capacity(bednetparams_1)
  cc
  
  simparams <- bednetparams_1 |>
    set_carrying_capacity(
      carrying_capacity = t(matrix(c(cc * c(1,
                                            1)),nrow = 2)),
      timesteps = c(3*year + 190)
    )
  
  
  correlationsb1 <- get_correlation_parameters(simparams)
  correlationsb1$inter_round_rho('bednets', 1)
  
  output_net1  <- run_simulation(sim_length, simparams,correlationsb1) 
  output_net1$pv_182.5_3650 = output_net1$n_detect_182.5_3650/output_net1$n_182.5_3650
  
  return(data.frame(timestep = output_net1$timestep,
                    net_use = output_net1$n_use_net/human_population,
                    prev_6m_10y = output_net1$pv_182.5_3650))
  
}

## net use estimates  
nets_use1 = data.frame(nets = nets$site1nets,timesteps = nets$timesteps)
nets_use1 = nets_use1[complete.cases(nets_use1$nets),]

nets_use2 = data.frame(nets = nets$site2nets,timesteps = nets$timesteps)
nets_use2 = nets_use2[complete.cases(nets_use2$nets),]

nets_use3 = data.frame(nets = nets$site3nets,timesteps = nets$timesteps)
nets_use3 = nets_use3[complete.cases(nets_use3$nets),]

nets_use4 = data.frame(nets = nets$site4nets,timesteps = nets$timesteps)
nets_use4 = nets_use4[complete.cases(nets_use4$nets),]

nets_use5 = data.frame(nets = nets$site5nets,timesteps = nets$timesteps)
nets_use5 = nets_use5[complete.cases(nets_use5$nets),]

nets_use6 = data.frame(nets = nets$site6nets,timesteps = nets$timesteps)
nets_use6 = nets_use6[complete.cases(nets_use6$nets),]

counterf1 = lsm_LOOP_f(parms_row = 1,
                       village = 1,
                       nets_use = nets_use1, 
                       root_eir = 117.8074)

counterf2 = lsm_LOOP_f(parms_row = 1,
                       village = 2,
                       nets_use = nets_use2, 
                       root_eir = 90.2171)

counterf3 = lsm_LOOP_f(parms_row = 1,
                       village = 3,
                       nets_use = nets_use3, 
                       root_eir = 75.2647)

counterf4 = lsm_LOOP_f(parms_row = 1,
                       village = 4,
                       nets_use = nets_use4, 
                       root_eir = 10.84384)

counterf5 = lsm_LOOP_f(parms_row = 1,
                       village = 5,
                       nets_use = nets_use5, 
                       root_eir = 34.12402)

counterf6 = lsm_LOOP_f(parms_row = 1,
                       village = 6,
                       nets_use = nets_use6, 
                       root_eir = 22.7)

write.csv(counterf1,"R/malariasimulation/model outputs/counterf1.csv")
write.csv(counterf2,"R/malariasimulation/model outputs/counterf2.csv")
write.csv(counterf3,"R/malariasimulation/model outputs/counterf3.csv")
write.csv(counterf4,"R/malariasimulation/model outputs/counterf4.csv")
write.csv(counterf5,"R/malariasimulation/model outputs/counterf5.csv")
write.csv(counterf6,"R/malariasimulation/model outputs/counterf6.csv")
