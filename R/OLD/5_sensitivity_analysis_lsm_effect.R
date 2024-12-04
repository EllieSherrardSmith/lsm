##############################################
##
## Sensitivity analysis of impacts
##
##############################################
library(malariasimulation)

## input params

## Simulations for 2 seasonal profiles, perennial / seasonal
## (annual mean 50% slide positive in 6-59 months) high, (20%) medium, (8%) low transmission profiles
## 40% ITN use, 50% treatment all ACTs
## in bed biting 60% and 90% mosquitoes (single species simulation)
## pyrethroid resistance ranging 0 - 100% in 21 increments
## LSM ranging 0 - 60% in 21 increments

set_up = expand.grid(eir_target = c(0.08,0.2,0.5),
                     phiB = c(0.6,0.9),
                     itn_cov = c(0.4,0.6,0.8),
                     resistance = round(seq(0,1,length=21),2),
                     cc_impact = 1 - seq(0,0.6,length=21)
                     )
netEff = read.csv("data/pyrethroid_only_nets.csv",header=TRUE)

set_up = merge(set_up, netEff,by = "resistance")
set_up$phiI = set_up$phiB + 0.05

dd = set_up[with(set_up, order(-cc_impact, itn_cov, resistance, phiB, eir_target)), ]

dd$root_eir = rep(eirs_calibrated$root_eir,21)

## store the simulation inputs data
write.csv(dd,"R/malariasimulation/scenario_simulations/inputs_set_up.csv")

lsm_sims_f = function(parms_row,root_eir){
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
  gam_params['phi_bednets'] <- set_up$phiB[parms_row] # proportion biting in bed: update from Nilani (sensitivity analysis?)
  gam_params['phi_indoors'] <- set_up$phiI[parms_row] # proportion biting indoors: update from Nilani (sensitivity analysis?)
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
    cover = rep(set_up$itn_cov[parms_row],3)
  )
  
  bednetparams <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    coverages = bednet_events$cover,
    
    retention = 5 * year, 
    
    dn0 = matrix(set_up$dn0_med[parms_row], nrow=nrow(bednet_events), ncol=1),
    rn =  matrix(set_up$rn0_med[parms_row], nrow=nrow(bednet_events), ncol=1),
    rnm = matrix(0.24, nrow=nrow(bednet_events), ncol=1),
    gamman = rep(set_up$gamman_med[parms_row] * 365, nrow(bednet_events)) 
    
  )

  ## Specify the LSM coverage
  cc <- get_init_carrying_capacity(bednetparams)
  cc
  
  simparams <- bednetparams |>
    set_carrying_capacity(
      carrying_capacity = t(matrix(c(cc * c(set_up$cc_impact[parms_row])),nrow = 1)),
      timesteps = c(7*year)
    )
  
  
  correlationsb1 <- get_correlation_parameters(simparams)
  correlationsb1$inter_round_rho('bednets', 1)
  
  output  <- run_simulation(sim_length, simparams,correlationsb1) 
  output$pv_182.5_1825 = output$n_detect_182.5_1825/output$n_182.5_1825 ## 6-59 months
  output$pv_0_36500 = output$n_detect_0_36500/output$n_0_36500 ## all-age

  return(output)
}

## Calibrate for EIR
## Params that do not matter are 
#  cc_impact

set_up = subset(set_up,set_up$cc_impact == 1)





#################################
##
## Calibrate model
##
#################################
library(cali)

## Baseline prevalence targets (which we wish to hit as an "average" for year preceding LSM)
## this is year 6-7, 
eir_store = numeric(length=nrow(set_up))

summary_pfprev_6m_5y = function(x){
  prev_6m_5y <- x$n_detect_182.5_1825/x$n_182.5_1825
  return(prev_6m_5y)
}

for(i in 3:nrow(set_up)){
  ###################################################################
  ## set up parameters
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
  gam_params['phi_bednets'] <- set_up$phiB[parms_row] # proportion biting in bed: update from Nilani (sensitivity analysis?)
  gam_params['phi_indoors'] <- set_up$phiI[parms_row] # proportion biting indoors: update from Nilani (sensitivity analysis?)
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
    cover = rep(set_up$itn_cov[parms_row],3)
  )
  
  bednetparams <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    coverages = bednet_events$cover,
    
    retention = 5 * year, 
    
    dn0 = matrix(set_up$dn0_med[parms_row], nrow=nrow(bednet_events), ncol=1),
    rn =  matrix(set_up$rn0_med[parms_row], nrow=nrow(bednet_events), ncol=1),
    rnm = matrix(0.24, nrow=nrow(bednet_events), ncol=1),
    gamman = rep(set_up$gamman_med[parms_row] * 365, nrow(bednet_events)) 
    
  )
  
  ## Specify the LSM coverage
  cc <- get_init_carrying_capacity(bednetparams)
  cc
  
  simparams <- bednetparams |>
    set_carrying_capacity(
      carrying_capacity = t(matrix(c(cc * c(set_up$cc_impact[parms_row])),nrow = 1)),
      timesteps = c(7*year)
    )
  
  
  correlationsb1 <- get_correlation_parameters(simparams)
  correlationsb1$inter_round_rho('bednets', 1)
  
  #############################################
  # Define target, here two prevalence measures:
  target <- set_up$eir_target[i]
  simparams$timesteps = 6.5*365
  set.seed(123)
  cal1 <- calibrate(target = target,
                    summary_function = summary_pfprev_6m_5y,
                    parameters = simparams,
                    tolerance = 0.01, 
                    low = 3, high = 250)##upper bound needs to be high enough so negative differences are not returned in uniroot
  eir_store[i] = cal1
}


a_root_eir = c(3.15,25,100,3.4,28,120)*(1+set_up$resistance)
set_up$root_eir = rep(a_root_eir[1:42],3)

for(i in 43:nrow(set_up)){
  simsOut[[i]] = lsm_sims_f(parms_row = i, 
                            root_eir = set_up$root_eir[i])
  
}

plot(test$pv_182.5_1825 ~ test$timestep,pch="",ylim=c(0,1),
     yaxt="n",xlab = "Time (days)")
axis(2, las=2, at=seq(0,1,0.2), labels = seq(0,100,20))
abline(h = c(0.08,0.2, 0.5),lty=2)
for(i in 1:nrow(set_up)){
  lines(simsOut[[i]]$pv_182.5_1825 ~ simsOut[[i]]$timestep,col="red")
}
eirs_calibrated = set_up

## Repeat for perennial scenario
## Then align EIR with the follow-up sims for different levels of cc 
## Code drawings for matrix plots

simsOut[[i]] = lsm_sims_f(parms_row = i, 
                          root_eir = set_up$root_eir[i])

##########################################
##
## We will want to draw a matrix of the cases averted / reduction in prev due to lsm
## need to compare like-with-like, the baseline scenario will be no lsm with matched
## itn_cov, resistance, phiB/phiI, eir_target

# Align the set_up row with the model results
# Calculate the total clinical cases over 1, 2 and 3 years given the scenario counterfactual
# Calculate the total clinical cases over 1, 2 and 3 years given the scenario with lsm for n levels
# Calculate the absolute and relative difference in these estimates

## plot the data as matrix of level of resistance versus lsm effect for: 
## each level of eir (low, med, high) with 40% ITN
## each level of eir with 60$ ITN
## each level of eir with 80% ITN
## for the seasonal and perennial settings

zf1 = mean_reduction_n3y_LLIN100 * 100

x = rep(seq(0,1,by=0.05),5)
y = rep(c(0.4,0.5964,0.814,0.8909,0.94,0.9832),each=21)

xcoords = unique(x)
ycoords = unique(y)

surface.mat1 = matrix(zf1,ncol=length(ycoords),nrow=length(xcoords),byrow=F)

## Figure 5
plot.new()
par(mar=c(4,5.5,2,2))
# major tick size and direction, < 0 means outside

matrix_funs = function(surface.matrix,minimum_val,maximum_val,upps,uni,levs,colschoice,cols_conts){
  filled.contour3(xcoords,
                  ycoords,
                  surface.matrix,
                  color=colschoice,
                  plot.axes = { axis(2, at = seq(0.4, 1, by = 0.1), seq(0.4, 1, by = 0.1),las=2)
                    axis(1, at = seq(0, 1, by = 0.2), labels = seq(0, 100, by = 20)) },
                  xlim = c(min(xcoords),max(xcoords)),
                  ylim = c(min(ycoords),max(ycoords)),
                  zlim = c(minimum_val,maximum_val))
  
  # the contour part - draw iso-lines
  contour(xcoords,
          ycoords,
          surface.matrix,
          color=blue2green,
          xlab = "",
          ylab = "",
          nlevels = levs, levels = seq(0,upps,by=uni),
          xlim = c(min(xcoords),max(xcoords)),
          ylim = c(min(ycoords),max(ycoords)),
          zlim = c(minimum_val,maximum_val),
          add=TRUE,                 # add the contour plot to filled-contour,
          #thus making an overlay
          col = cols_conts         # color of overlay-lines
  )
}

##mid left plot: Explanation figure, Pervalence estimate given different phiI and phiB estimates
par(new = "TRUE",  
    plt = c(0.55,0.9,0.1,0.45)  )                 # major tick size and direction, < 0 means outside


matrix_funs(surface.mat1,min(surface.mat1),max(surface.mat1),
            upps=100,uni = 10,levs = 11,colschoice = heat.colors,cols_conts="darkblue")
text(0.95,0.95,"D",cex=1.2,col="darkblue")

######################################################################
#Add a legend:
par(new = "TRUE",
    plt = c(0.93,0.95,0.18,0.4),   # define plot region for legend
    las = 1,
    cex.axis = 1)
#
filled.legend(
  xcoords,
  ycoords,
  surface.mat1,
  color = heat.colors,
  plot.title = "Efficacy (%)",
  xlab = "",
  ylab = "",
  xlim = c(0,100),
  ylim = c(0.4,1),
  zlim = c(min(surface.mat1),max(surface.mat1)))

#Add some figure labels
par(xpd=NA,cex = 1.1)
text(x = -21.6,y = 34,expression(paste("Proportion of mosquito bites taken indoors, ", phi[I])),srt = 90,cex = 0.9)
text(x = -10,y = -26,"Mosquito survival at bioassay (%)",cex = 0.9)

text(x = -16,y = 75,"100% LLIN cover",cex = 0.9)
text(x = -10,y = 82.8,"Efficacy against clinical cases in children under 5 (%)",cex = 0.9)
