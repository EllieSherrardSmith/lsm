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
      carrying_capacity = t(matrix(c(cc * c(phi_lsm_gam[village],
                                            phi_lsm_fun[village])),nrow = 2)),
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

v1 = list()

for(i in 1:nrow(parms)){
  v1[[i]] = lsm_LOOP_f(parms_row = i,village = 1,nets_use = nets_use1, root_eir = 117.8074)
}

v2 = list()
v3 = list()
v4 = list()
v5 = list()
v6 = list()

for(i in 1:nrow(parms)){
  v2[[i]] = lsm_LOOP_f(parms_row = i,village = 2,nets_use = nets_use2, root_eir = 90.2171)
  v3[[i]] = lsm_LOOP_f(parms_row = i,village = 3,nets_use = nets_use3, root_eir = 75.2647)
  v4[[i]] = lsm_LOOP_f(parms_row = i,village = 4,nets_use = nets_use4, root_eir = 10.84384)
  v5[[i]] = lsm_LOOP_f(parms_row = i,village = 5,nets_use = nets_use5, root_eir = 34.12402)
  v6[[i]] = lsm_LOOP_f(parms_row = i,village = 6,nets_use = nets_use6, root_eir = 22.7)
}

#######################
##
## Predictions altogether Village level
##
######################

## Create medians and limits files
##MEDIANS
range_estimates_v1 = array(data=NA, dim=c(nrow(v1[[1]]),81) )
range_estimates_v2 = array(data=NA, dim=c(nrow(v2[[1]]),81) )
range_estimates_v3 = array(data=NA, dim=c(nrow(v3[[1]]),81) )
range_estimates_v4 = array(data=NA, dim=c(nrow(v4[[1]]),81) )
range_estimates_v5 = array(data=NA, dim=c(nrow(v5[[1]]),81) )
range_estimates_v6 = array(data=NA, dim=c(nrow(v6[[1]]),81) )
for(i in 1:81){
  range_estimates_v1[,i] = v1[[i]]$prev_6m_10y
  range_estimates_v2[,i] = v2[[i]]$prev_6m_10y
  range_estimates_v3[,i] = v3[[i]]$prev_6m_10y
  range_estimates_v4[,i] = v4[[i]]$prev_6m_10y
  range_estimates_v5[,i] = v5[[i]]$prev_6m_10y
  range_estimates_v6[,i] = v6[[i]]$prev_6m_10y
}
range_est1 = expand.grid(timestep = v1[[1]]$timestep)
range_est2 = expand.grid(timestep = v2[[1]]$timestep)
range_est3 = expand.grid(timestep = v3[[1]]$timestep)
range_est4 = expand.grid(timestep = v4[[1]]$timestep)
range_est5 = expand.grid(timestep = v5[[1]]$timestep)
range_est6 = expand.grid(timestep = v6[[1]]$timestep)

range_est1[,2] = rowMeans(range_estimates_v1)
range_est2[,2] = rowMeans(range_estimates_v2)
range_est3[,2] = rowMeans(range_estimates_v3)
range_est4[,2] = rowMeans(range_estimates_v4)
range_est5[,2] = rowMeans(range_estimates_v5)
range_est6[,2] = rowMeans(range_estimates_v6)

for(i in 1:nrow(range_est1)){
  range_est1[i,3] = as.numeric(quantile(range_estimates_v1[i,],0.025))
  range_est1[i,4] = as.numeric(quantile(range_estimates_v1[i,],0.975))

  range_est2[i,3] = as.numeric(quantile(range_estimates_v2[i,],0.025))
  range_est2[i,4] = as.numeric(quantile(range_estimates_v2[i,],0.975))

  range_est3[i,3] = as.numeric(quantile(range_estimates_v3[i,],0.025))
  range_est3[i,4] = as.numeric(quantile(range_estimates_v3[i,],0.975))

  range_est4[i,3] = as.numeric(quantile(range_estimates_v4[i,],0.025))
  range_est4[i,4] = as.numeric(quantile(range_estimates_v4[i,],0.975))

  range_est5[i,3] = as.numeric(quantile(range_estimates_v5[i,],0.025))
  range_est5[i,4] = as.numeric(quantile(range_estimates_v5[i,],0.975))

  range_est6[i,3] = as.numeric(quantile(range_estimates_v6[i,],0.025))
  range_est6[i,4] = as.numeric(quantile(range_estimates_v6[i,],0.975))
}

write.csv(range_est1,"carrying_cap_model_outputs_village1.csv")
write.csv(range_est2,"carrying_cap_model_outputs_village2.csv")
write.csv(range_est3,"carrying_cap_model_outputs_village3.csv")
write.csv(range_est4,"carrying_cap_model_outputs_village4.csv")
write.csv(range_est5,"carrying_cap_model_outputs_village5.csv")
write.csv(range_est6,"model_outputs_village6.csv")


#####################################################
##
## Figure 2
##
#####################################################
par(mfrow = c(2,3))

plot(range_est1[,2] ~ range_est1[,1],
     main = "Musilongo: LSM",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",pch="",col="darkblue",
     xlim = c(2*year,5*year+30))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

polygon(c(range_est1[,1],rev(range_est1[,1])),
        c(range_est1[,3],rev(range_est1[,4])),
        col=adegenet::transp("blue",0.4),border=NA)
lines(range_est1[,2] ~ range_est1[,1],col="darkblue",pch=1.5)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites1[1:18],na.rm=TRUE),18),rev(rep(min(sites1[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites1 ~ align_timeline, col="blue")


plot(range_est2[,2] ~ range_est2[,1],
     main = "Kimingini",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",pch="",col="darkblue",
     xlim = c(2*year,5*year+30))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

polygon(c(range_est2[,1],rev(range_est2[,1])),
        c(range_est2[,3],rev(range_est2[,4])),
        col=adegenet::transp("blue",0.4),border=NA)
lines(range_est2[,2] ~ range_est2[,1],col="darkblue",pch=1.5)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites2[1:18],na.rm=TRUE),18),rev(rep(min(sites2[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites2 ~ align_timeline, col="blue")

plot(range_est3[,2] ~ range_est3[,1],
     main = "Kezege: LSM",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",pch="",col="darkblue",
     xlim = c(2*year,5*year+30))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

polygon(c(range_est3[,1],rev(range_est3[,1])),
        c(range_est3[,3],rev(range_est3[,4])),
        col=adegenet::transp("blue",0.4),border=NA)
lines(range_est3[,2] ~ range_est3[,1],col="darkblue",pch=1.5)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites3[1:18],na.rm=TRUE),18),rev(rep(min(sites3[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites3 ~ align_timeline, col="blue")

plot(range_est4[,2] ~ range_est4[,1],
     main = "Wamondo: LSM",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",pch="",col="darkblue",
     xlim = c(2*year,5*year+30))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

polygon(c(range_est4[,1],rev(range_est4[,1])),
        c(range_est4[,3],rev(range_est4[,4])),
        col=adegenet::transp("blue",0.4),border=NA)
lines(range_est4[,2] ~ range_est3[,1],col="darkblue",pch=1.5)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites4[1:18],na.rm=TRUE),18),rev(rep(min(sites4[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites4 ~ align_timeline, col="blue")


plot(range_est5[,2] ~ range_est5[,1],
     main = "Emutete",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",pch="",col="darkblue",
     xlim = c(2*year,5*year+30))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

polygon(c(range_est5[,1],rev(range_est5[,1])),
        c(range_est5[,3],rev(range_est5[,4])),
        col=adegenet::transp("blue",0.4),border=NA)
lines(range_est5[,2] ~ range_est3[,1],col="darkblue",pch=1.5)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites5[1:18],na.rm=TRUE),18),rev(rep(min(sites5[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites5 ~ align_timeline, col="blue")


plot(range_est6[,2] ~ range_est6[,1],
     main = "Wakikuyu",
     ylab="Prevalence 6-months to 10-years (%)", yaxt="n",ylim = c(0,1),
     xlab="Time in days",pch="",col="darkblue",
     xlim = c(2*year,5*year+30))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
abline(v = 3*year + 190,lty=2)

polygon(c(range_est6[,1],rev(range_est6[,1])),
        c(range_est6[,3],rev(range_est6[,4])),
        col=adegenet::transp("blue",0.4),border=NA)
lines(range_est6[,2] ~ range_est3[,1],col="darkblue",pch=1.5)

## add prev data
polygon(c(align_timeline[1:18],rev(align_timeline[1:18])),
        c(rep(max(sites6[1:18],na.rm=TRUE),18),rev(rep(min(sites6[1:18],na.rm=TRUE),18))),
        border=NA,col=adegenet::transp("grey",0.2))
points(sites6 ~ align_timeline, col="blue")


#####################################################
##
## Figure 3
##
#####################################################

prev_cr_sect_mod = prev_cr_sect_modL = prev_cr_sect_modU = array(dim=c(length(align_timeline),18))

prev_cr_sect_mod[,1] = range_est1[align_timeline,2]
prev_cr_sect_mod[,2] = range_est2[align_timeline,2]
prev_cr_sect_mod[,3] = range_est3[align_timeline,2]
prev_cr_sect_mod[,4] = range_est4[align_timeline,2]
prev_cr_sect_mod[,5] = range_est5[align_timeline,2]
prev_cr_sect_mod[,6] = range_est6[align_timeline,2]

prev_cr_sect_modL[,1] = range_est1[align_timeline,3]
prev_cr_sect_modL[,2] = range_est2[align_timeline,3]
prev_cr_sect_modL[,3] = range_est3[align_timeline,3]
prev_cr_sect_modL[,4] = range_est4[align_timeline,3]
prev_cr_sect_modL[,5] = range_est5[align_timeline,3]
prev_cr_sect_modL[,6] = range_est6[align_timeline,3]

prev_cr_sect_modU[,1] = range_est1[align_timeline,4]
prev_cr_sect_modU[,2] = range_est2[align_timeline,4]
prev_cr_sect_modU[,3] = range_est3[align_timeline,4]
prev_cr_sect_modU[,4] = range_est4[align_timeline,4]
prev_cr_sect_modU[,5] = range_est5[align_timeline,4]
prev_cr_sect_modU[,6] = range_est6[align_timeline,4]

par(mfrow=c(1,2))
plot(sites1[19:48]~prev_cr_sect_mod[19:48,1],pch=19,col="darkorange",
     ylim=c(0,0.8),xlim=c(0,0.8),yaxt="n",xaxt="n",
     ylab="Trial observed prevalence 6-months to 10-years (%)",
     xlab="Model simulated prevalence 6-months to 10-years (%)")
axis(1,at=seq(0,0.8,0.2),labels=seq(0,80,20))
axis(2,las=2,at=seq(0,0.8,0.2),labels=seq(0,80,20))

abline(a=0,b=1,lty=2,col="grey")

points(sites2[19:48]~prev_cr_sect_mod[19:48,2],pch=19,col="lightblue")
points(sites3[19:48]~prev_cr_sect_mod[19:48,3],pch=19,col="red")
points(sites4[19:48]~prev_cr_sect_mod[19:48,4],pch=19,col="darkred")
points(sites5[19:48]~prev_cr_sect_mod[19:48,5],pch=19,col="blue")
points(sites6[19:48]~prev_cr_sect_mod[19:48,6],pch=19,col="darkblue")


for(i in 19:48){
  segments(x0 = prev_cr_sect_modL[i,1],
           x1 = prev_cr_sect_modU[i,1],
           y0 = sites1[i],y1 = sites1[i],col="darkorange")
  segments(x0 = prev_cr_sect_modL[i,2],
           x1 = prev_cr_sect_modU[i,2],
           y0 = sites2[i],y1 = sites2[i],col="lightblue")
  segments(x0 = prev_cr_sect_modL[i,3],
           x1 = prev_cr_sect_modU[i,3],
           y0 = sites3[i],y1 = sites3[i],col="red")
  segments(x0 = prev_cr_sect_modL[i,4],
           x1 = prev_cr_sect_modU[i,4],
           y0 = sites4[i],y1 = sites4[i],col="darkred")
  segments(x0 = prev_cr_sect_modL[i,5],
           x1 = prev_cr_sect_modU[i,5],
           y0 = sites5[i],y1 = sites5[i],col="blue")
  segments(x0 = prev_cr_sect_modL[i,6],
           x1 = prev_cr_sect_modU[i,6],
           y0 = sites6[i],y1 = sites6[i],col="darkblue")
  
}

##########################
##
## Fig 3B

## trial data
relative_efficacy_from_t0_v1 = (mean(sites1[1:18],na.rm=TRUE) - sites1[19:48])/mean(sites1[1:18],na.rm=TRUE)
relative_efficacy_from_t0_v2 = (mean(sites2[1:18],na.rm=TRUE) - sites2[19:48])/mean(sites2[1:18],na.rm=TRUE)
relative_efficacy_from_t0_v3 = (mean(sites3[1:18],na.rm=TRUE) - sites3[19:48])/mean(sites3[1:18],na.rm=TRUE)
relative_efficacy_from_t0_v4 = (mean(sites4[1:18],na.rm=TRUE) - sites4[19:48])/mean(sites4[1:18],na.rm=TRUE)
relative_efficacy_from_t0_v5 = (mean(sites5[1:18],na.rm=TRUE) - sites5[19:48])/mean(sites5[1:18],na.rm=TRUE)
relative_efficacy_from_t0_v6 = (mean(sites6[1:18],na.rm=TRUE) - sites6[19:48])/mean(sites6[1:18],na.rm=TRUE)

## model outputs
rel_eff_mod_from_t0_v1 = (mean(prev_cr_sect_mod[1:18,1],na.rm=TRUE) - prev_cr_sect_mod[19:48,1])/mean(prev_cr_sect_mod[1:18,1],na.rm=TRUE)
rel_eff_mod_from_t0_v2 = (mean(prev_cr_sect_mod[1:18,2],na.rm=TRUE) - prev_cr_sect_mod[19:48,2])/mean(prev_cr_sect_mod[1:18,2],na.rm=TRUE)
rel_eff_mod_from_t0_v3 = (mean(prev_cr_sect_mod[1:18,3],na.rm=TRUE) - prev_cr_sect_mod[19:48,3])/mean(prev_cr_sect_mod[1:18,3],na.rm=TRUE)
rel_eff_mod_from_t0_v4 = (mean(prev_cr_sect_mod[1:18,4],na.rm=TRUE) - prev_cr_sect_mod[19:48,4])/mean(prev_cr_sect_mod[1:18,4],na.rm=TRUE)
rel_eff_mod_from_t0_v5 = (mean(prev_cr_sect_mod[1:18,5],na.rm=TRUE) - prev_cr_sect_mod[19:48,5])/mean(prev_cr_sect_mod[1:18,5],na.rm=TRUE)
rel_eff_mod_from_t0_v6 = (mean(prev_cr_sect_mod[1:18,6],na.rm=TRUE) - prev_cr_sect_mod[19:48,6])/mean(prev_cr_sect_mod[1:18,6],na.rm=TRUE)

rel_eff_mod_from_t0_v1L = (mean(prev_cr_sect_mod[1:18,1],na.rm=TRUE) - prev_cr_sect_modL[19:48,1])/mean(prev_cr_sect_mod[1:18,1],na.rm=TRUE)
rel_eff_mod_from_t0_v2L = (mean(prev_cr_sect_mod[1:18,2],na.rm=TRUE) - prev_cr_sect_modL[19:48,2])/mean(prev_cr_sect_mod[1:18,2],na.rm=TRUE)
rel_eff_mod_from_t0_v3L = (mean(prev_cr_sect_mod[1:18,3],na.rm=TRUE) - prev_cr_sect_modL[19:48,3])/mean(prev_cr_sect_mod[1:18,3],na.rm=TRUE)
rel_eff_mod_from_t0_v4L = (mean(prev_cr_sect_mod[1:18,4],na.rm=TRUE) - prev_cr_sect_modL[19:48,4])/mean(prev_cr_sect_mod[1:18,4],na.rm=TRUE)
rel_eff_mod_from_t0_v5L = (mean(prev_cr_sect_mod[1:18,5],na.rm=TRUE) - prev_cr_sect_modL[19:48,5])/mean(prev_cr_sect_mod[1:18,5],na.rm=TRUE)
rel_eff_mod_from_t0_v6L = (mean(prev_cr_sect_mod[1:18,6],na.rm=TRUE) - prev_cr_sect_modL[19:48,6])/mean(prev_cr_sect_mod[1:18,6],na.rm=TRUE)

rel_eff_mod_from_t0_v1U = (mean(prev_cr_sect_mod[1:18,1],na.rm=TRUE) - prev_cr_sect_modU[19:48,1])/mean(prev_cr_sect_mod[1:18,1],na.rm=TRUE)
rel_eff_mod_from_t0_v2U = (mean(prev_cr_sect_mod[1:18,2],na.rm=TRUE) - prev_cr_sect_modU[19:48,2])/mean(prev_cr_sect_mod[1:18,2],na.rm=TRUE)
rel_eff_mod_from_t0_v3U = (mean(prev_cr_sect_mod[1:18,3],na.rm=TRUE) - prev_cr_sect_modU[19:48,3])/mean(prev_cr_sect_mod[1:18,3],na.rm=TRUE)
rel_eff_mod_from_t0_v4U = (mean(prev_cr_sect_mod[1:18,4],na.rm=TRUE) - prev_cr_sect_modU[19:48,4])/mean(prev_cr_sect_mod[1:18,4],na.rm=TRUE)
rel_eff_mod_from_t0_v5U = (mean(prev_cr_sect_mod[1:18,5],na.rm=TRUE) - prev_cr_sect_modU[19:48,5])/mean(prev_cr_sect_mod[1:18,5],na.rm=TRUE)
rel_eff_mod_from_t0_v6U = (mean(prev_cr_sect_mod[1:18,6],na.rm=TRUE) - prev_cr_sect_modU[19:48,6])/mean(prev_cr_sect_mod[1:18,6],na.rm=TRUE)

plot(relative_efficacy_from_t0_v1~rel_eff_mod_from_t0_v1,pch=19,col="darkorange",
     ylim=c(0,1),xlim=c(0,1),yaxt="n",xaxt="n",
     ylab="Trial observed change from baseline prevalence (%)",
     xlab="Model simulated change from baseline prevalence (%)")
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))

abline(a=0,b=1,lty=2,col="grey")

points(relative_efficacy_from_t0_v2~rel_eff_mod_from_t0_v2,pch=19,col="lightblue")
points(relative_efficacy_from_t0_v3~rel_eff_mod_from_t0_v3,pch=19,col="red")
points(relative_efficacy_from_t0_v4~rel_eff_mod_from_t0_v4,pch=19,col="darkred")
points(relative_efficacy_from_t0_v5~rel_eff_mod_from_t0_v5,pch=19,col="blue")
points(relative_efficacy_from_t0_v6~rel_eff_mod_from_t0_v6,pch=19,col="darkblue")

for(i in 1:30){
  segments(x0 = rel_eff_mod_from_t0_v1L[i],
           x1 = rel_eff_mod_from_t0_v1U[i],
           y0 = relative_efficacy_from_t0_v1[i],y1 = relative_efficacy_from_t0_v1[i],col="darkorange")
  segments(x0 = rel_eff_mod_from_t0_v2L[i],
           x1 = rel_eff_mod_from_t0_v2U[i],
           y0 = relative_efficacy_from_t0_v2[i],y1 = relative_efficacy_from_t0_v2[i],col="lightblue")
  segments(x0 = rel_eff_mod_from_t0_v3L[i],
           x1 = rel_eff_mod_from_t0_v3U[i],
           y0 = relative_efficacy_from_t0_v3[i],y1 = relative_efficacy_from_t0_v3[i],col="red")
  segments(x0 = rel_eff_mod_from_t0_v4L[i],
           x1 = rel_eff_mod_from_t0_v4U[i],
           y0 = relative_efficacy_from_t0_v4[i],y1 = relative_efficacy_from_t0_v4[i],col="darkred")
  segments(x0 = rel_eff_mod_from_t0_v5L[i],
           x1 = rel_eff_mod_from_t0_v5U[i],
           y0 = relative_efficacy_from_t0_v5[i],y1 = relative_efficacy_from_t0_v5[i],col="blue")
  segments(x0 = rel_eff_mod_from_t0_v6L[i],
           x1 = rel_eff_mod_from_t0_v6U[i],
           y0 = relative_efficacy_from_t0_v6[i],y1 = relative_efficacy_from_t0_v6[i],col="darkblue")
  
}


##########################
##
## Fig 3C
vill1_matched = c(23,25,28,31,35,37)
vill5_matched = c(24,25,28,31,35,37)
vill6_matched =    c(25,28,31,35,37)
range_estimates_v1_matched = range_estimates_v1[align_timeline,]
range_estimates_v1_matched = range_estimates_v1_matched[vill1_matched,]
dim(range_estimates_v1_matched)
vill1_sites_matched = sites1[vill1_matched]

range_estimates_v2_matched = range_estimates_v2[align_timeline,]
range_estimates_v2_matched = range_estimates_v2_matched[vill1_matched,]
vill2_sites_matched = sites2[vill1_matched]

range_estimates_v3_matched = range_estimates_v3[align_timeline,]
range_estimates_v3_matched = range_estimates_v3_matched[vill1_matched,]
vill3_sites_matched = sites3[vill1_matched]

range_estimates_v4_matched = range_estimates_v4[align_timeline,]
range_estimates_v4_matched = range_estimates_v4_matched[vill1_matched,]
vill4_sites_matched = sites4[vill1_matched]

range_estimates_v5_matched = range_estimates_v5[align_timeline,]
range_estimates_v5_matched = range_estimates_v5_matched[vill5_matched,]
vill5_sites_matched = sites5[vill5_matched]

range_estimates_v6_matched = range_estimates_v6[align_timeline,]
range_estimates_v6_matched = range_estimates_v6_matched[vill6_matched,]
vill6_sites_matched = sites6[vill6_matched]

vill1_difference = vill2_difference = 
  vill3_difference = vill4_difference = 
  vill5_difference = array(data=NA,dim=c(81,6))

for(i in 1:6){ ## these are now the matched timesteps for observed and simulated data
  vill1_difference[,i] = vill1_sites_matched[i] - range_estimates_v1_matched[i,]
  vill2_difference[,i] = vill2_sites_matched[i] - range_estimates_v2_matched[i,]
  vill3_difference[,i] = vill3_sites_matched[i] - range_estimates_v3_matched[i,]
  vill4_difference[,i] = vill4_sites_matched[i] - range_estimates_v4_matched[i,]
  vill5_difference[,i] = vill5_sites_matched[i] - range_estimates_v5_matched[i,]
}

vill6_difference = array(data=NA,dim=c(81,6))

for(i in 1:5){ ## these are now the matched timesteps for observed and simulated data
  vill6_difference[,i] = vill6_sites_matched[i] - range_estimates_v6_matched[i,]
}

par(mfrow=c(1,1))
plot(vill1_difference[,1] ~ sample(seq(0.8,1.2,length=80),size=81,replace=TRUE),
     xlim=c(0,41),ylim=c(-0.5,0.5),cex=0.5,col=adegenet::transp("darkorange",0.3),pch=19, 
     xlab = "Timestep matched",xaxt="n",yaxt="n",
     ylab="Absolute difference in observed versus simulated prevalence estimate")
abline(h=0,lty=2)
axis(1,at=c(3,9.5,16.5,23.5,30.5,37.5),labels=c("Nov/Dec 2005","Jan 2006", "Apr 2006", "Jul 2006", "Nov 2006", "Jan 2007"))
axis(2,las=2,seq(-0.4,0.4,0.1),labels = seq(-0.4,0.4,0.1))

points(vill2_difference[,1] ~ sample(seq(2-0.2,2+0.2,length=80),size=81,replace=TRUE),
         cex=0.5,col=adegenet::transp("lightblue",0.5),pch=19)  
  
points(vill3_difference[,1] ~ sample(seq(3-0.2,3+0.2,length=80),size=81,replace=TRUE),
         cex=0.5,col=adegenet::transp("red",0.5),pch=19)  
  
points(vill4_difference[,1] ~ sample(seq(4-0.2,4+0.2,length=80),size=81,replace=TRUE),
         cex=0.5,col=adegenet::transp("darkred",0.5),pch=19)  
  
points(vill5_difference[,1] ~ sample(seq(5-0.2,5+0.2,length=80),size=81,replace=TRUE),
         cex=0.5,col=adegenet::transp("blue",0.5),pch=19)  


## timestep 2
points(vill1_difference[,2] ~ sample(seq(7-0.2,7+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkorange",0.3),pch=19)  

points(vill2_difference[,2] ~ sample(seq(8-0.2,8+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("lightblue",0.5),pch=19)  

points(vill3_difference[,2] ~ sample(seq(9-0.2,9+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("red",0.5),pch=19)  

points(vill4_difference[,2] ~ sample(seq(10-0.2,10+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkred",0.5),pch=19)  

points(vill5_difference[,2] ~ sample(seq(11-0.2,11+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("blue",0.5),pch=19)  

points(vill6_difference[,1] ~ sample(seq(12-0.2,12+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkblue",0.5),pch=19)  


## timestep 3 apr 2006
points(vill1_difference[,3] ~ sample(seq(14-0.2,14+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkorange",0.3),pch=19)  

points(vill2_difference[,3] ~ sample(seq(15-0.2,15+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("lightblue",0.5),pch=19)  

points(vill3_difference[,3] ~ sample(seq(16-0.2,16+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("red",0.5),pch=19)  

points(vill4_difference[,3] ~ sample(seq(17-0.2,17+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkred",0.5),pch=19)  

points(vill5_difference[,3] ~ sample(seq(18-0.2,18+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("blue",0.5),pch=19)  

points(vill6_difference[,2] ~ sample(seq(19-0.2,19+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkblue",0.5),pch=19)  

## timestep 4 jul 2006
points(vill1_difference[,4] ~ sample(seq(21-0.2,21+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkorange",0.3),pch=19)  

points(vill2_difference[,4] ~ sample(seq(22-0.2,22+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("lightblue",0.5),pch=19)  

points(vill3_difference[,4] ~ sample(seq(23-0.2,23+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("red",0.5),pch=19)  

points(vill4_difference[,4] ~ sample(seq(24-0.2,24+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkred",0.5),pch=19)  

points(vill5_difference[,4] ~ sample(seq(25-0.2,25+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("blue",0.5),pch=19)  

points(vill6_difference[,3] ~ sample(seq(26-0.2,26+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkblue",0.5),pch=19)  


## timestep 5 jul 2006
points(vill1_difference[,5] ~ sample(seq(28-0.2,28+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkorange",0.3),pch=19)  

points(vill2_difference[,5] ~ sample(seq(29-0.2,29+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("lightblue",0.5),pch=19)  

points(vill3_difference[,5] ~ sample(seq(30-0.2,30+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("red",0.5),pch=19)  

points(vill4_difference[,5] ~ sample(seq(31-0.2,31+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkred",0.5),pch=19)  

points(vill5_difference[,5] ~ sample(seq(32-0.2,32+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("blue",0.5),pch=19)  

points(vill6_difference[,4] ~ sample(seq(33-0.2,33+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkblue",0.5),pch=19)  


## timestep 5 jul 2006
points(vill1_difference[,6] ~ sample(seq(35-0.2,35+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkorange",0.3),pch=19)  

points(vill2_difference[,6] ~ sample(seq(36-0.2,36+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("lightblue",0.5),pch=19)  

points(vill3_difference[,6] ~ sample(seq(37-0.2,37+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("red",0.5),pch=19)  

points(vill4_difference[,6] ~ sample(seq(38-0.2,38+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkred",0.5),pch=19)  

points(vill5_difference[,6] ~ sample(seq(39-0.2,39+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("blue",0.5),pch=19)  

points(vill6_difference[,5] ~ sample(seq(40-0.2,40+0.2,length=80),size=81,replace=TRUE),
       cex=0.5,col=adegenet::transp("darkblue",0.5),pch=19)  

# for(i in 2:6){
#   points(vill1_difference[,i] ~ sample(seq(i-0.2,i+0.2,length=80),size=81,replace=TRUE),
#          cex=0.5,col=adegenet::transp("darkorange",0.3),pch=19)  
# }
# 
# for(i in 1:6){
#   points(vill2_difference[,i] ~ sample(seq(7+i-0.2,7+i+0.2,length=80),size=81,replace=TRUE),
#          cex=0.5,col=adegenet::transp("lightblue",0.5),pch=19)  
# 
#   points(vill3_difference[,i] ~ sample(seq(14+i-0.2,14+i+0.2,length=80),size=81,replace=TRUE),
#          cex=0.5,col=adegenet::transp("red",0.5),pch=19)  
# 
#   points(vill4_difference[,i] ~ sample(seq(21+i-0.2,21+i+0.2,length=80),size=81,replace=TRUE),
#          cex=0.5,col=adegenet::transp("darkred",0.5),pch=19)  
# 
#   points(vill5_difference[,i] ~ sample(seq(28+i-0.2,28+i+0.2,length=80),size=81,replace=TRUE),
#          cex=0.5,col=adegenet::transp("blue",0.5),pch=19)  
# }
# 
# for(i in 1:5){
#   points(vill6_difference[,i] ~ sample(seq(35+i-0.2,35+i+0.2,length=80),size=81,replace=TRUE),
#          cex=0.5,col=adegenet::transp("darkblue",0.5),pch=19)  
# }

legend("topright",legend = c("Musilongo: LSM",
                             "Kimingini",
                             "Kezege: LSM",
                             "Wamondo: LSM",
                             "Emutete",
                             "Wakikuyu"),
       col = c("darkorange","lightblue","red","darkred",
               "blue","darkblue"),ncol=2,pch=15,bty="n")

text(7,0.45,"Model simulation over estimates effect")
text(7,-0.45,"Model simulation under estimates effect")

# summary(parms[which(vill1_difference[,1] < 0.02),])
# summary(parms[which(vill1_difference[,2] > -0.2),])
# which(vill1_difference[,3] > -0.02);parms[which(vill1_difference[,3] > -0.02),]
# for(i in 1:6){
#   points(vill1_difference[which(parms$gam_Q0 == max(parms$gam_Q0)),i] ~ 
#            sample(seq(i-0.2,i+0.2,length=length(which(parms$gam_Q0 == max(parms$gam_Q0)))),
#                   size=length(which(parms$gam_Q0 == max(parms$gam_Q0))),replace=TRUE),
#          col="darkred",cex=0.5,pch=19)
#   
# }

