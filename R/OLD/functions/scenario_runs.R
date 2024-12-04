lsm_sims_f = function(root_eir,
                      phiB,
                      phiI,
                      itn_cov,
                      dn0_med,
                      rn0_med,
                      gamman_med,
                      cc_impact,
                      row_num
                      ){
  ## LSM modelling 
  ## background parameters
  year <- 365
  month <- 30
  sim_length <- 10 * year
  human_population <- 10000
  starting_EIR <- root_eir  ## This will change to hit estimated prevalence in U5 yrs of (5%, 10%, 20%, 30%, 40%, 50%, and 60%)
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      
      prevalence_rendering_min_ages = c(0.5,5, 15, 0) * 365, ## Prev in under 5 years measured
      prevalence_rendering_max_ages = c(5,  15,100,100) * 365,
      
      clinical_incidence_rendering_min_ages =  c(0.5,5, 15, 0) * 365,  ## All age clin_inc
      clinical_incidence_rendering_max_ages =  c(5,  15,100,100) * 365,
      
      severe_incidence_rendering_min_ages =  c(0.5,5, 15, 0) * 365,  ## All age sev_inc
      severe_incidence_rendering_max_ages =  c(5,  15,100,100) * 365,
      
      # mint parameters for two profiles
      # 0.285505	-0.325352	-0.132815	-0.0109352	0.104675	0.0779865	-0.013919
      # seasonal_a0	seasonal_a1	seasonal_b1	seasonal_a2	seasonal_b2	seasonal_a3	seasonal_b3	
      # 0.285277	-0.0248801	-0.0216681	-0.0529426	-0.0242904	-0.016891	-0.0073646
      
      model_seasonality = TRUE, # seasonal model
      g0 = 0.285505,
      g = c(-0.325352,  -0.0109352,  0.0779865),
      h = c(-0.132815,   0.104675,   -0.013919),
      
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
  simparams <- set_drugs(simparams, list(AL_params)) ## ACT drug
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,
                                      time=c(1),
                                      coverage=c(0.5))    
  
  # set mosquito species
  gam_params <- gamb_params # gambiae
  
  gam_params['species'] <- "gambiae sl"
  gam_params['blood_meal_rates'] <- 0.33 # 1/duration of gonothropic cycle
  gam_params['foraging_time'] <- 0.69 # time spent foraging
  gam_params['Q0'] <- 0.92 # human blood index: update from Nilani (sensitivity analysis?)
  gam_params['phi_bednets'] <- phiB # proportion biting in bed: update from Nilani (sensitivity analysis?)
  gam_params['phi_indoors'] <- phiI # proportion biting indoors: update from Nilani (sensitivity analysis?)
  gam_params['mum'] <- 0.132 # death rate or 1/life expectancy: update from Nilani (sensitivity analysis?)
  
  simparams <- set_species(simparams,
                           species=list(gam_params),
                           proportions=c(1)
                           
  )
  
  bednetparams <- simparams
  # set bednets
  bednet_events = data.frame(
    timestep = c(1,4,7)*365, ## every 3 years
    names = rep("ITN campaign",3),
    cover = rep(itn_cov,3)
  )
  
  bednetparams <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    coverages = bednet_events$cover,
    
    retention = 5 * year, 
    
    dn0 = matrix(dn0_med, nrow=nrow(bednet_events), ncol=1),
    rn =  matrix(rn0_med, nrow=nrow(bednet_events), ncol=1),
    rnm = matrix(0.24, nrow=nrow(bednet_events), ncol=1),
    gamman = rep(gamman_med * 365, nrow(bednet_events)) 
    
  )
  
  ## Specify the LSM coverage
  cc <- get_init_carrying_capacity(bednetparams)
  cc
  
  simparams <- bednetparams |>
    set_carrying_capacity(
      carrying_capacity = t(matrix(c(cc * c(cc_impact)),nrow = 1)),
      timesteps = c(7*year)
    )
  
  
  correlationsb1 <- get_correlation_parameters(simparams)
  correlationsb1$inter_round_rho('bednets', 1)
  
  output  <- run_simulation(sim_length, simparams,correlationsb1) 
  output$pv_182.5_1825 = output$n_detect_182.5_1825/output$n_182.5_1825 ## 6-59 months
  output$pv_0_36500 = output$n_detect_0_36500/output$n_0_36500 ## all-age
  
  front_part <- paste0("outputs/sim_", row_num, "_run")
  
  write.csv(output,
            file = paste0(front_part,".csv"), row.names = FALSE)
}
