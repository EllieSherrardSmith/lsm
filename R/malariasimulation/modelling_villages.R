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


# LSM impact
phi_lsm_gam = c(0, 1, 0.22, 0.88, 1, 1)
phi_lsm_ara = c(1, 1, 1, 1, 1, 1)
phi_lsm_fun = c(0.67, 1, 0.45, 0.11, 1, 1)
rho_alpha = -4
rho_beta = 0.057


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

timeline = 1:48 ## JAN 2004 UNTIL DEC 2007

year <- 365

align_timeline = round(seq(2*year+1,2*year+1+48 * 365/12,length=48),0)



## LSM modelling 

year <- 365
month <- 30
sim_length <- 10 * year
human_population <- 10000
starting_EIR <- 50  ## This will change to hit estimated prevalence in U5 yrs of (5%, 10%, 20%, 30%, 40%, 50%, and 60%)

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
simparams <- set_species(simparams,
                         species=list(gamb_params, arab_params, fun_params),
                         proportions=c(gamb_prop[6],0,fun_prop[6])
                         
)


bednetparams <- simparams
# set bednets

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

bednet_events = data.frame(
  timestep = nets_use6$timesteps,
  names = rep("use",length(nets_use6$timesteps))
  # timestep = c(0, 
  #              2, ## 2004
  #              3, ## 2005 july
  #              4, ## 2006 july
  #              5, 8) * year + c(1,1,190,190,190,1),
  # name=c("Historic1", "Historic2","Y_2005","Y_2006","2007","after")
)

## Bednet 1 corresponds to pyrethroid-only nets
bednetparams_1 <- set_bednets(
  bednetparams,
  
  timesteps = bednet_events$timestep,
  
  coverages = nets_use6$nets,
  
  # coverages = c(hist_nets[6], ## BASELINE 
  #               hist_nets[6], ## BASELINE 
  #               july_2005_nets[6], ## models july 2005
  #               july_2006_nets[6], ## models july 2006
  #               july_2006_nets[6],july_2006_nets[6]),   ## A future estimate (relatively irrelevant so far but likely to want options for i coverages...)
  
  retention = 5 * year, 
  dn0 = matrix(0.387, nrow=nrow(nets_use6), ncol=3),
  ## These all assume no resistance (2001-2011)
  
  rn =  matrix(0.563, nrow=nrow(nets_use6), ncol=3),
  ## 
  
  rnm = matrix(0.24, nrow=nrow(nets_use6), ncol=3),
  ## This is fixed ... something we prob need to investigate statistically at some point!
  
  gamman = rep(2.64 * 365, nrow(nets_use6)) 
  ## These change with the level of resistance
  
)

## LSM
lsm_events = data.frame(
  habitat_management_timesteps = 3 * year+190, ## july 2005 
  name=c("trial-starts")
)

habitat_management <- set_habitat_management(
  bednetparams_1,
  habitat_management_timesteps = lsm_events$habitat_management_timesteps,
  lsm_new_eqm = matrix(c(phi_lsm_gam[6],phi_lsm_ara[6],phi_lsm_fun[6]),nrow=1,ncol=3),
  lsm_rate_alpha = matrix(-4,nrow=1,ncol=3), # exp(4)
  lsm_rate_beta = matrix(0.057,nrow=1,ncol=3)
)

correlationsb1 <- get_correlation_parameters(habitat_management)
correlationsb1$inter_round_rho('bednets', 1)

output_net1  <- run_simulation(sim_length, habitat_management,correlationsb1) 
output_net1$pv_182.5_3650 = output_net1$n_detect_182.5_3650/output_net1$n_182.5_3650

#################################
##
## Calibrate to 2005 
##
#################################


## Baseline prevalence
# 63.3% (50.7 -70.0)	60.1% (42.7 – 72.7)	55.8% (41.7– 77.1)	25.7% (12.6 – 41.4)	44.9% (15.4 – 66.7)	40.2% (18.9 – 63.3)
prev_base = c(0.633,0.601,0.558,0.257,0.449,0.402)

# Define target, here two prevalence measures:
target <- prev_base[6] 
# Time points at which to match target
target_tt <- c(3*year+160)

summary_pfprev_6m_10y = function(x){
  prev_6m_10y <- x$n_detect_182.5_3650/x$n_182.5_3650
  return(prev_6m_10y)
}

set.seed(123)
library(cali)
out_net6 <- calibrate(parameters = habitat_management,
                      target = target,
                      target_tt = target_tt,
                      summary_function = summary_pfprev_6m_10y,
                      tolerance = 0.01, 
                      interval = c(1, 250))##upper bound needs to be high enough so negative differences are not returned in uniroot

out_net6$root

###########################################
##
## function for running the model 

lsm_f = function(i,nets_use,root_eir){
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
  simparams <- set_species(simparams,
                           species=list(gamb_params, arab_params, fun_params),
                           proportions=c(gamb_prop[i],0,fun_prop[i])
                           
  )
  
   
  bednetparams <- simparams
  # set bednets
  bednet_events = data.frame(
    timestep = nets_use$timesteps,
    names = rep("use",length(nets_use$timesteps))
    
    # timestep = c(0, 
    #              3, ## 2004
    #              4, ## 2005 july
    #              5, ## 2006 july
    #              6, 9) * year + c(1,1,190,190,190,1),
    # name=c("Historic1", "Historic2","Y_2005","Y_2006","2007","after")
  )
  
  ## Bednet 1 corresponds to pyrethroid-only nets
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    coverages = nets_use$nets,
    
    # coverages = c(hist_nets[i], ## BASELINE 
    #               hist_nets[i], ## BASELINE 
    #               july_2005_nets[i], ## models july 2005
    #               july_2006_nets[i], ## models july 2006
    #               july_2006_nets[i],july_2006_nets[i]),   ## A future estimate (relatively irrelevant so far but likely to want options for i coverages...)
    
    retention = 5 * year, 
    dn0 = matrix(0.340773639, nrow=nrow(nets_use), ncol=3),
    ## These all assume no resistance (2001-2011)
    
    rn =  matrix(0.637100833, nrow=nrow(nets_use), ncol=3),
    ## 
    
    rnm = matrix(0.24, nrow=nrow(nets_use), ncol=3),
    ## This is fixed ... something we prob need to investigate statistically at some point!
    
    gamman = rep(2.64 * 365, nrow(nets_use)) 
    ## These change with the level of resistance
    
  )
  
  ## LSM
  lsm_events = data.frame(
    habitat_management_timesteps = 3 * year+190, ## july 2005 
    name=c("trial-starts")
  )
  
  habitat_management <- set_habitat_management(
    bednetparams_1,
    habitat_management_timesteps = lsm_events$habitat_management_timesteps,
    lsm_new_eqm = matrix(c(phi_lsm_gam[i],phi_lsm_ara[i],phi_lsm_fun[i]),nrow=1,ncol=3),
    lsm_rate_alpha = matrix(-4,nrow=1,ncol=3), # exp(4)
    lsm_rate_beta = matrix(0.057,nrow=1,ncol=3)
  )
  
  correlationsb1 <- get_correlation_parameters(habitat_management)
  correlationsb1$inter_round_rho('bednets', 1)
  
  output_net1  <- run_simulation(sim_length, habitat_management,correlationsb1) 
  output_net1$pv_182.5_3650 = output_net1$n_detect_182.5_3650/output_net1$n_182.5_3650
  
  return(data.frame(timestep = output_net1$timestep,
                    net_use = output_net1$n_nets/human_population,
                    prev_6m_10y = output_net1$pv_182.5_3650))
  
}

## All come out too high currently
##Options - try Joel's approach?
##Options - Calibrate to just the final data point June 2005

# out_net1$root = 157.8074
# out_net2$root = 100.2171
# out_net3$root = 135.2647
# out_net4$root = 23.84384
# out_net5$root = 54.12402
# out_net6$root = 34.888

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

village_1 = lsm_f(i=1,nets_use = nets_use1, root_eir = 117.8074)
village_2 = lsm_f(i=2,nets_use = nets_use2, root_eir = 90.2171)
village_3 = lsm_f(i=3,nets_use = nets_use3, root_eir = 75.2647)
village_4 = lsm_f(i=4,nets_use = nets_use4, root_eir = 10.84384)
village_5 = lsm_f(i=5,nets_use = nets_use5, root_eir = 34.12402)
village_6 = lsm_f(i=6,nets_use = nets_use6, root_eir = 22.7)

par(mfrow=c(1,2))
plot(village_1$net_use ~ village_1$timestep,
     ylim = c(0,1),ylab = "Proportion of people using nets",
     xlab = "Time in days",xlim = c(2*year,6*year))
abline(v = 3*year+190,lty=2)
text(1850,1,"July 2005")

plot(village_1$prev_6m_10y ~ village_1$timestep,
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year+190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites1[1:18],na.rm=TRUE),18),rev(rep(min(sites1[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites1 ~ align_timeline, col="blue")


par(mfrow=c(1,2))
plot(village_2$net_use ~ village_2$timestep,
     ylim = c(0,1),ylab = "Proportion of people using nets",
     xlab = "Time in days",xlim = c(2*year,6*year))
abline(v = 3*year+190,lty=2)
text(1850,1,"July 2005")

plot(village_2$prev_6m_10y ~ village_2$timestep,
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year+190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites2[1:18],na.rm=TRUE),18),rev(rep(min(sites2[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites2 ~ align_timeline, col="blue")


par(mfrow=c(1,2))
plot(village_3$net_use ~ village_3$timestep,
     ylim = c(0,1),ylab = "Proportion of people using nets",
     xlab = "Time in days",xlim = c(2*year,6*year))
abline(v = 3*year + 190,lty=2)
text(1850,1,"July 2005")

plot(village_3$prev_6m_10y ~ village_3$timestep,
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites3[1:18],na.rm=TRUE),18),rev(rep(min(sites3[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites3 ~ align_timeline, col="blue")


par(mfrow=c(1,2))
plot(village_4$net_use ~ village_4$timestep,
     ylim = c(0,1),ylab = "Proportion of people using nets",
     xlab = "Time in days",xlim = c(2*year,6*year))
abline(v = 3*year + 190,lty=2)
text(1850,1,"July 2005")

plot(village_4$prev_6m_10y ~ village_4$timestep,
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites4[1:18],na.rm=TRUE),18),rev(rep(min(sites4[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites4 ~ align_timeline, col="blue")


par(mfrow=c(1,2))
plot(village_5$net_use ~ village_5$timestep,
     ylim = c(0,1),ylab = "Proportion of people using nets",
     xlab = "Time in days",xlim = c(2*year,6*year))
abline(v = 3*year + 190,lty=2)
text(1850,1,"July 2005")

plot(village_5$prev_6m_10y ~ village_5$timestep,
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites5[1:18],na.rm=TRUE),18),rev(rep(min(sites5[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites5 ~ align_timeline, col="blue")


par(mfrow=c(1,2))
plot(village_6$net_use ~ village_6$timestep,
     ylim = c(0,1),ylab = "Proportion of people using nets",
     xlab = "Time in days",xlim = c(2*year,6*year))
abline(v = 3*year + 190,lty=2)
text(1850,1,"July 2005")

plot(village_6$prev_6m_10y ~ village_6$timestep,
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites6[1:18],na.rm=TRUE),18),rev(rep(min(sites6[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites6 ~ align_timeline, col="blue")

#######################
##
## Predictions altogether Village level
##
######################
par(mfrow = c(2,3))

plot(village_1$prev_6m_10y ~ village_1$timestep,
     main = "Musilongo: LSM",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites1[1:18],na.rm=TRUE),18),rev(rep(min(sites1[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites1 ~ align_timeline, col="blue")



plot(village_2$prev_6m_10y ~ village_2$timestep,
     main = "Kimingini",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites2[1:18],na.rm=TRUE),18),rev(rep(min(sites2[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites2 ~ align_timeline, col="blue")

plot(village_3$prev_6m_10y ~ village_3$timestep,
     main = "Kezege: LSM",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites3[1:18],na.rm=TRUE),18),rev(rep(min(sites3[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites3 ~ align_timeline, col="blue")



plot(village_4$prev_6m_10y ~ village_4$timestep,
     main = "Wamondo: LSM",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites4[1:18],na.rm=TRUE),18),rev(rep(min(sites4[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites4 ~ align_timeline, col="blue")


plot(village_5$prev_6m_10y ~ village_5$timestep,
     main = "Emutete",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites5[1:18],na.rm=TRUE),18),rev(rep(min(sites5[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites5 ~ align_timeline, col="blue")


plot(village_6$prev_6m_10y ~ village_6$timestep,
     main = "Wakikuyu",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",type="l",col="darkblue",
     xlim = c(2*year,6*year))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites6[1:18],na.rm=TRUE),18),rev(rep(min(sites6[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites6 ~ align_timeline, col="blue")


############################
##
## Figure S3 Net use
##
############################
par(mfrow=c(2,3))
plot(village_1$net_use ~ village_1$timestep,
     main = "Musilongo: LSM",type="l",
     ylim = c(0,1),yaxt="n",
     ylab = "Proportion of people using nets",
     xlim = c(2*year,5*year+30),xaxt="n",
     xlab = "Time in days")
axis(1,at=c(2,3,4,5)*year, labels=c("Jan 2004","Jan 2005","Jan 2006","Jan 2007"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year+190,lty=2)
text(1400,1,"July 2005")

plot(village_2$net_use ~ village_2$timestep,
     main = "Kimingini",type="l",
     ylim = c(0,1),yaxt="n",
     ylab = "Proportion of people using nets",
     xlim = c(2*year,5*year+30),xaxt="n",
     xlab = "Time in days")
axis(1,at=c(2,3,4,5)*year, labels=c("Jan 2004","Jan 2005","Jan 2006","Jan 2007"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year+190,lty=2)
text(1400,1,"July 2005")

plot(village_3$net_use ~ village_3$timestep,
     main = "Kezege: LSM",type="l",
     ylim = c(0,1),yaxt="n",
     ylab = "Proportion of people using nets",
     xlim = c(2*year,5*year+30),xaxt="n",
     xlab = "Time in days")
axis(1,at=c(2,3,4,5)*year, labels=c("Jan 2004","Jan 2005","Jan 2006","Jan 2007"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year+190,lty=2)
text(1400,1,"July 2005")

plot(village_4$net_use ~ village_4$timestep,
     main = "Wamondo: LSM",type="l",
     ylim = c(0,1),yaxt="n",
     ylab = "Proportion of people using nets",
     xlim = c(2*year,5*year+30),xaxt="n",
     xlab = "Time in days")
axis(1,at=c(2,3,4,5)*year, labels=c("Jan 2004","Jan 2005","Jan 2006","Jan 2007"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year+190,lty=2)
text(1400,1,"July 2005")

plot(village_5$net_use ~ village_5$timestep,
     main = "Emutete",type="l",
     ylim = c(0,1),yaxt="n",
     ylab = "Proportion of people using nets",
     xlim = c(2*year,5*year+30),xaxt="n",
     xlab = "Time in days")
axis(1,at=c(2,3,4,5)*year, labels=c("Jan 2004","Jan 2005","Jan 2006","Jan 2007"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year+190,lty=2)
text(1400,1,"July 2005")

plot(village_6$net_use ~ village_6$timestep,
     main = "Wakikuyu",type="l",
     ylim = c(0,1),yaxt="n",
     ylab = "Proportion of people using nets",
     xlim = c(2*year,5*year+30),xaxt="n",
     xlab = "Time in days")
axis(1,at=c(2,3,4,5)*year, labels=c("Jan 2004","Jan 2005","Jan 2006","Jan 2007"))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year+190,lty=2)
text(1400,1,"July 2005")
