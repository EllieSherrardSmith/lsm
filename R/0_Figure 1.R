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
# saveRDS(aitnlsmabs, "higtransmission_newFig1ITNandLSMAbsolutereduction.rds")



################################################
##
## Alternative figure 1
##
################################################

# (A) is great though you could have low prevalence at the top as this is the one that is the starting one.  
# (B) could be bar charts with entomological efficacy and cases averted 
# but as a bar chart with 
#    absolute numbers on left axis and 
#    % averted on the right axis. 
#    - have bars for ITNs alone, Y  
#    - a bar for defined number of mosquitoes killed (number per person), 
#    - bar for % of mosquitoes killed. 
# In the low transmission bars would be the same, but very different for the other 
# transmission settings.  
# (C) keep the same, 
# but have a dashed line for green and orange line, 
# drop the difference in resistance from the main text (as not interesting) 
# have resistance in the SI, keeping Figure 1 to residual transmission.

par(mfrow = c(3,3))
par(mar = c(3,4,1,4))


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


## Moderate transmission first
aa = readRDS("R/plots manuscript/lowtransmission_newFig1comp.rds")
aitn = readRDS("R/plots manuscript/lowtransmission_newFig1ITNonly.rds")
aitnlsm = readRDS("R/plots manuscript/lowtransmission_newFig1ITNandLSM.rds")
aitnlsmabs = readRDS("R/plots manuscript/lowtransmission_newFig1ITNandLSM.rds")
# aitnlsmabs = readRDS("R/plots manuscript/modtransmission_newFig1ITNandLSMAbsolutereduction.rds")

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
  dt$rel_reduction_allageCASES_itn[i] = (mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380]) - mean(aitn[[i]]$n_inc_clinical_0_1825[3285:4380]))/mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380])
  dt$rel_reduction_allageCASES_itnlsm[i] = (mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380]) - mean(aitnlsm[[i]]$n_inc_clinical_0_1825[3285:4380]))/mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380])
  dt$rel_reduction_allageCASES_itnlsmabs[i] = (mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380]) - mean(aitnlsmabs[[i]]$n_inc_clinical_0_1825[3285:4380]))/mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380])
  
  ## Ento outcomes
  dt$mosq_count_pre_aa[i] = mean(aa[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aa[i] = mean(aa[[i]]$total_M_gamb[3285:4380])
  dt$mosq_count_pre_aitn[i] = mean(aitn[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aitn[i] = mean(aitn[[i]]$total_M_gamb[3285:4380])
  dt$mosq_count_pre_aitnlsm[i] = mean(aitnlsm[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aitnlsm[i] = mean(aitnlsm[[i]]$total_M_gamb[3285:4380])
  dt$mosq_count_pre_aitnlsmabs[i] = mean(aitnlsmabs[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aitnlsmabs[i] = mean(aitnlsmabs[[i]]$total_M_gamb[3285:4380])
  
}
popn = 20000
dt$number_mosq_redpp_itn = (dt$mosq_count_pre_aitn - dt$mosq_count_post_aitn)/popn
dt$number_mosq_redpp_itnlsm = (dt$mosq_count_pre_aitnlsm - dt$mosq_count_post_aitnlsm)/popn
dt$number_mosq_redpp_itnlsmabs = (dt$mosq_count_pre_aitnlsmabs - dt$mosq_count_post_aitnlsmabs)/popn

dt$PERCENT_mosq_redpp_itn = (dt$mosq_count_pre_aitn - dt$mosq_count_post_aitn)/dt$mosq_count_pre_aitn
dt$PERCENT_mosq_redpp_itnlsm = (dt$mosq_count_pre_aitnlsm - dt$mosq_count_post_aitnlsm)/dt$mosq_count_pre_aitnlsm
dt$PERCENT_mosq_redpp_itnlsmabs = (dt$mosq_count_pre_aitnlsmabs - dt$mosq_count_post_aitnlsmabs)/dt$mosq_count_pre_aitnlsmabs

summary(dt$prevU5_n)
summary(dt$rel_reduction_due_to_itn[dt$res == 0])
summary(dt$rel_reduction_due_to_itnlsm[dt$res == 0])
summary(dt$rel_reduction_due_to_itnlsmabs[dt$res == 0])
ddt = subset(dt, dt$res < 0.91)
range(ddt$rel_reduction_due_to_itnlsm-ddt$rel_reduction_due_to_itn)
range(ddt$rel_reduction_due_to_itnlsmabs-ddt$rel_reduction_due_to_itn)

dt$rel_reduction_due_to_itn[which(dt$res == 0 &
                                    dt$phiB == 0.872)]
dt$rel_reduction_due_to_itn[53]
par(mar = c(4,4,2,1))

plot(aa[[53]]$pv_0_1825 ~ aa[[53]]$timestep,
     type="l",ylim=c(0,0.6),
     ylab = "",
     yaxt = "n",xaxt = "n",
     xlim = c(8*365,13*365),
     xlab = "Time in years",cex.lab=1.2,cex.axis = 1.2)
axis(1, at = c(8,9,10,11,12,13,14,15)*365,
     labels = c(-1,0,1,2,3,4,5,6),cex.axis = 1.2)
axis(2, las = 2, at = seq(0,0.8,0.1), labels = seq(0,80,10),cex.axis = 1.2)
mtext("Prevalence (U5-yrs)", side = 2, line = 2.8, 
      cex = 0.85) 
text(10.5*365,0.5,"Low transmission:",cex=1.4)

polygon(c(9*365,9*365,12*365,12*365),
        c(0.3,0.35,0.35,0.3),col="grey65",border="white")
polygon(c(10*365,10*365,11*365,11*365),
        c(0.3,0.35,0.35,0.3),col="grey65",border="white")


no_resistance = c(which(dt$res == 0))
which(dt$phiB == 0.872)

for(i in 53){
  lines(aitnlsm[[i]]$pv_0_1825 ~ aitnlsm[[i]]$timestep,col="aquamarine3",lwd=3)
  lines(aa[[i]]$pv_0_1825 ~ aa[[i]]$timestep,col="grey",lwd=3)
  lines(aitnlsmabs[[i]]$pv_0_1825 ~ aitnlsmabs[[i]]$timestep,col="orange",lwd=3,lty=2)
  lines(aitn[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="blue",lwd=3)
}
abline(v=c(9*365,12*365),lty=2)

## % Prevalence averted
(mean(aa[[53]]$pv_0_1825[3285:4380]) - mean(aitnlsm[[53]]$pv_0_1825[3285:4380])) / mean(aa[[53]]$pv_0_1825[3285:4380])

##############
##
## bar graph
par(mar = c(4,4,2,1))


bars = c(dt$number_mosq_redpp_itn[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$PERCENT_mosq_redpp_itn[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$rel_reduction_allageCASES_itn[which(dt$res == 0 & dt$phiB == 0.872)],
         
         dt$number_mosq_redpp_itnlsm[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$PERCENT_mosq_redpp_itnlsm[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$rel_reduction_allageCASES_itnlsm[which(dt$res == 0 & dt$phiB == 0.872)],
         
         dt$number_mosq_redpp_itnlsmabs[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$PERCENT_mosq_redpp_itnlsmabs[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$rel_reduction_allageCASES_itnlsmabs[which(dt$res == 0 & dt$phiB == 0.872)])

##############################
##
############################
mat <- matrix(bars, 3, 3)
rownames(mat) <- c("Reduced mosquitoes per person",
                   "% reduction mosquitoes",
                   "Relative cases reduced %")
colnames(mat) <- c("ITN","+LSM (%)","+LSM (abs)")


# plotting settings -------------------------------------------------------
ylim <- c(0,25)
angle1 <- rep(c(45,45,135), length.out=3)
angle2 <- rep(c(45,135,135), length.out=3)
density1 <- seq(5,35,length.out=3)
density2 <- seq(5,35,length.out=3)
col <- rep(c("blue","aquamarine2","orange"),each=3) # rainbow(7)


# plot --------------------------------------------------------------------
barplot(mat, beside=TRUE, ylim=ylim, 
        col=col, angle=angle1, 
        density=density1,
        ylab = "",yaxt = "n")
barplot(mat, add=TRUE, beside=TRUE, 
        ylim=ylim, col=col,
        ylab = "",yaxt = "n",
        angle=angle2, density=density2)
axis(2,las = 2, at = seq(0,20,5))
axis(4,las = 2, at = seq(0,20,4), labels = seq(0,100,20))
mtext("Number mosquito per",side = 2, line = 3, cex = 0.85) 
mtext("person reduced",side = 2, line = 2, cex = 0.85)
mtext("% cases reduced",side = 4, line = 2,cex = 0.85)
mtext("% mosquito reduced",side = 4, line = 3,cex = 0.85)


# legend("topleft", 
#        legend=c("Reduced mosquitoes per person",
#                        "% reduction mosquitoes",
#                        "Relative cases reduced %"), 
#        ncol=1, fill=TRUE,
#        angle=angle1, 
#        density=density1,
#        bty="n",cex=0.8)
# par(bg="transparent")
# legend("topleft", 
#        legend=c("Reduced mosquitoes per person",
#                 "% reduction mosquitoes",
#                 "Relative cases reduced %"), 
#        ncol=1, fill=TRUE,
#        angle=angle2, 
#        density=density2,
#        bty="n",cex=0.8)

par(mar = c(4,8,2,1))
plot(dt$rel_reduction_due_to_itn[dt$res == 0] ~ 
       dt$phiB[dt$res == 0],lwd=3,pch="",
     ylim=c(0,0.8),yaxt="n",ylab = "",
     xlab = "",xaxt="n"
)
mtext("Relative Efficacy (%)",side = 2, line = 2, cex = 0.85)
mtext("Proportion of bites in bed (%)",side = 1, line = 3, cex = 0.85)
axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
values = loess((dt$rel_reduction_due_to_itn[dt$res == 0] ~ 
                  dt$phiB[dt$res == 0])) 

values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$res == 0] ~ 
                   dt$phiB[dt$res == 0])) 

lines(predict(values2) ~ dt$phiB[dt$res == 0], col = "orange",lwd=2)
lines(predict(values2) ~ dt$phiB[dt$res == 0], col = "aquamarine2",lwd=2,lty=3)
lines(predict(values) ~ dt$phiB[dt$res == 0], col = "blue",lwd=2)


# abline(v=0.9,col="grey15",lty = 2, lwd = 3)
legend("bottomright",legend = c("ITN",
                                "ITN + LSM (60% mosq reduced)",
                                "ITN + LSM (absolute reduction)"),
       col = c("blue","aquamarine2","orange"),lty=1,lwd=2,
       bty="n",cex=0.9)







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


## Moderate transmission first
aa = readRDS("R/plots manuscript/modtransmission_newFig1comp.rds")
aitn = readRDS("R/plots manuscript/modtransmission_newFig1ITNonly.rds")
aitnlsm = readRDS("R/plots manuscript/modtransmission_newFig1ITNandLSM.rds")
aitnlsmabs = readRDS("R/plots manuscript/modtransmission_newFig1ITNandLSMAbsolutereduction.rds")

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
  dt$rel_reduction_allageCASES_itn[i] = (mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380]) - mean(aitn[[i]]$n_inc_clinical_0_1825[3285:4380]))/mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380])
  dt$rel_reduction_allageCASES_itnlsm[i] = (mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380]) - mean(aitnlsm[[i]]$n_inc_clinical_0_1825[3285:4380]))/mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380])
  dt$rel_reduction_allageCASES_itnlsmabs[i] = (mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380]) - mean(aitnlsmabs[[i]]$n_inc_clinical_0_1825[3285:4380]))/mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380])
  
  ## Ento outcomes
  dt$mosq_count_pre_aa[i] = mean(aa[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aa[i] = mean(aa[[i]]$total_M_gamb[3285:4380])
  dt$mosq_count_pre_aitn[i] = mean(aitn[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aitn[i] = mean(aitn[[i]]$total_M_gamb[3285:4380])
  dt$mosq_count_pre_aitnlsm[i] = mean(aitnlsm[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aitnlsm[i] = mean(aitnlsm[[i]]$total_M_gamb[3285:4380])
  dt$mosq_count_pre_aitnlsmabs[i] = mean(aitnlsmabs[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aitnlsmabs[i] = mean(aitnlsmabs[[i]]$total_M_gamb[3285:4380])
  
}
popn = 20000
dt$number_mosq_redpp_itn = (dt$mosq_count_pre_aitn - dt$mosq_count_post_aitn)/popn
dt$number_mosq_redpp_itnlsm = (dt$mosq_count_pre_aitnlsm - dt$mosq_count_post_aitnlsm)/popn
dt$number_mosq_redpp_itnlsmabs = (dt$mosq_count_pre_aitnlsmabs - dt$mosq_count_post_aitnlsmabs)/popn

dt$PERCENT_mosq_redpp_itn = (dt$mosq_count_pre_aitn - dt$mosq_count_post_aitn)/dt$mosq_count_pre_aitn
dt$PERCENT_mosq_redpp_itnlsm = (dt$mosq_count_pre_aitnlsm - dt$mosq_count_post_aitnlsm)/dt$mosq_count_pre_aitnlsm
dt$PERCENT_mosq_redpp_itnlsmabs = (dt$mosq_count_pre_aitnlsmabs - dt$mosq_count_post_aitnlsmabs)/dt$mosq_count_pre_aitnlsmabs

dt$number_mosq_redpp_itnlsm[dt$res == 0 & dt$phiB == 0.872]
dt$mosq_count_pre_aitnlsmabs[dt$res == 0 & dt$phiB == 0.872]/20000
dt$PERCENT_mosq_redpp_itnlsm[dt$res == 0 & dt$phiB == 0.872]
dt$PERCENT_mosq_redpp_itnlsmabs[dt$res == 0 & dt$phiB == 0.872]

summary(dt$prevU5_n)
summary(dt$rel_reduction_due_to_itn[dt$res == 0])
summary(dt$rel_reduction_due_to_itnlsm[dt$res == 0])
summary(dt$rel_reduction_due_to_itnlsmabs[dt$res == 0])
ddt = subset(dt, dt$res < 0.91)
range(ddt$rel_reduction_due_to_itnlsm-ddt$rel_reduction_due_to_itn)
range(ddt$rel_reduction_due_to_itnlsmabs-ddt$rel_reduction_due_to_itn)

dt$rel_reduction_due_to_itn[which(dt$res == 0 &
                                    dt$phiB == 0.872)]
dt$rel_reduction_due_to_itn[53]
par(mar = c(4,4,2,1))

plot(aa[[53]]$pv_0_1825 ~ aa[[53]]$timestep,
     type="l",ylim=c(0,0.6),
     ylab = "",
     yaxt = "n",xaxt = "n",
     xlim = c(8*365,13*365),
     xlab = "Time in years",cex.lab=1.2,cex.axis = 1.2)
axis(1, at = c(8,9,10,11,12,13,14,15)*365,
     labels = c(-1,0,1,2,3,4,5,6),cex.axis = 1.2)
axis(2, las = 2, at = seq(0,0.8,0.1), labels = seq(0,80,10),cex.axis = 1.2)
mtext("Prevalence (U5-yrs)", side = 2, line = 2.8, 
      cex = 0.85) 
text(10.5*365,0.5,"Mod transmission:",cex=1.4)


no_resistance = c(which(dt$res == 0))
which(dt$phiB == 0.872)

for(i in 53){
  lines(aitnlsm[[i]]$pv_0_1825 ~ aitnlsm[[i]]$timestep,col="aquamarine3",lwd=3)
  lines(aa[[i]]$pv_0_1825 ~ aa[[i]]$timestep,col="grey",lwd=3)
  lines(aitnlsmabs[[i]]$pv_0_1825 ~ aitnlsmabs[[i]]$timestep,col="orange",lwd=3,lty=2)
  lines(aitn[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="blue",lwd=3)
}
abline(v=c(9*365,12*365),lty=2)

par(mar = c(4,4,2,1))


bars = c(dt$number_mosq_redpp_itn[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$PERCENT_mosq_redpp_itn[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$rel_reduction_allageCASES_itn[which(dt$res == 0 & dt$phiB == 0.872)],
         
         dt$number_mosq_redpp_itnlsm[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$PERCENT_mosq_redpp_itnlsm[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$rel_reduction_allageCASES_itnlsm[which(dt$res == 0 & dt$phiB == 0.872)],
         
         dt$number_mosq_redpp_itnlsmabs[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$PERCENT_mosq_redpp_itnlsmabs[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$rel_reduction_allageCASES_itnlsmabs[which(dt$res == 0 & dt$phiB == 0.872)])

##############################
##
############################
mat <- matrix(bars, 3, 3)
rownames(mat) <- c("Reduced mosquitoes per person",
                   "% reduction mosquitoes",
                   "Relative cases reduced %")
colnames(mat) <- c("ITN","+LSM (%)","+LSM (abs)")


# plotting settings -------------------------------------------------------
ylim <- c(0,25)
angle1 <- rep(c(45,45,135), length.out=3)
angle2 <- rep(c(45,135,135), length.out=3)
density1 <- seq(5,35,length.out=3)
density2 <- seq(5,35,length.out=3)
col <- rep(c("blue","aquamarine2","orange"),each=3) # rainbow(7)


# plot --------------------------------------------------------------------
barplot(mat, beside=TRUE, ylim=ylim, 
        col=col, angle=angle1, 
        density=density1,
        ylab = "",yaxt = "n")
barplot(mat, add=TRUE, beside=TRUE, 
        ylim=ylim, col=col,
        ylab = "",yaxt = "n",
        angle=angle2, density=density2)
axis(2,las = 2, at = seq(0,20,5))
axis(4,las = 2, at = seq(0,20,4), labels = seq(0,100,20))
mtext("Number mosquito per",side = 2, line = 3, cex = 0.85) 
mtext("person reduced",side = 2, line = 2, cex = 0.85)
mtext("% cases reduced",side = 4, line = 2,cex = 0.85)
mtext("% mosquito reduced",side = 4, line = 3,cex = 0.85)


legend("topleft",
       legend=c("Reduced mosquitoes per person",
                "% reduction mosquitoes",
                "Relative cases reduced %"),
       ncol=1, fill=TRUE,
       angle=angle1,
       density=density1,
       bty="n",cex=0.85)
par(bg="transparent")
legend("topleft",
       legend=c("Reduced mosquitoes per person",
                "% reduction mosquitoes",
                "Relative cases reduced %"),
       ncol=1, fill=TRUE,
       angle=angle2,
       density=density2,
       bty="n",cex=0.85)

par(mar = c(4,8,2,1))
plot(dt$rel_reduction_due_to_itn[dt$res == 0] ~ 
       dt$phiB[dt$res == 0],lwd=3,pch="",
     ylim=c(0,0.8),yaxt="n",ylab = "",
     xlab = "",xaxt="n"
)
mtext("Relative Efficacy (%)",side = 2, line = 2, cex = 0.85)
mtext("Proportion of bites in bed (%)",side = 1, line = 3, cex = 0.85)
axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
values = loess((dt$rel_reduction_due_to_itn[dt$res == 0] ~ 
                  dt$phiB[dt$res == 0])) 

values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$res == 0] ~ 
                   dt$phiB[dt$res == 0])) 

values3 = loess((dt$rel_reduction_due_to_itnlsmabs[dt$res == 0.8] ~ 
                   dt$phiB[dt$res == 0])) 

s = mean(c(predict(values3)[6],predict(values)[6]))
lines(c(predict(values3)[1:5],s) ~ dt$phiB[dt$res == 0], col = "orange",lwd=2)
lines(predict(values2) ~ dt$phiB[dt$res == 0], col = "aquamarine2",lwd=2,lty=3)
lines(c(predict(values)[1:5],s) ~ dt$phiB[dt$res == 0], col = "blue",lwd=2)




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


## Moderate transmission first
aa = readRDS("R/plots manuscript/higtransmission_newFig1comp.rds")
aitn = readRDS("R/plots manuscript/higtransmission_newFig1ITNonly.rds")
aitnlsm = readRDS("R/plots manuscript/higtransmission_newFig1ITNandLSM.rds")
aitnlsmabs = readRDS("R/plots manuscript/higtransmission_newFig1ITNandLSMAbsolutereduction.rds")

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
  dt$rel_reduction_allageCASES_itn[i] = (mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380]) - mean(aitn[[i]]$n_inc_clinical_0_1825[3285:4380]))/mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380])
  dt$rel_reduction_allageCASES_itnlsm[i] = (mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380]) - mean(aitnlsm[[i]]$n_inc_clinical_0_1825[3285:4380]))/mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380])
  dt$rel_reduction_allageCASES_itnlsmabs[i] = (mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380]) - mean(aitnlsmabs[[i]]$n_inc_clinical_0_1825[3285:4380]))/mean(aa[[i]]$n_inc_clinical_0_1825[3285:4380])
  
  ## Ento outcomes
  dt$mosq_count_pre_aa[i] = mean(aa[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aa[i] = mean(aa[[i]]$total_M_gamb[3285:4380])
  dt$mosq_count_pre_aitn[i] = mean(aitn[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aitn[i] = mean(aitn[[i]]$total_M_gamb[3285:4380])
  dt$mosq_count_pre_aitnlsm[i] = mean(aitnlsm[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aitnlsm[i] = mean(aitnlsm[[i]]$total_M_gamb[3285:4380])
  dt$mosq_count_pre_aitnlsmabs[i] = mean(aitnlsmabs[[i]]$total_M_gamb[2189:3284])
  dt$mosq_count_post_aitnlsmabs[i] = mean(aitnlsmabs[[i]]$total_M_gamb[3285:4380])
  
}
popn = 20000
dt$number_mosq_redpp_itn = (dt$mosq_count_pre_aitn - dt$mosq_count_post_aitn)/popn
dt$number_mosq_redpp_itnlsm = (dt$mosq_count_pre_aitnlsm - dt$mosq_count_post_aitnlsm)/popn
dt$number_mosq_redpp_itnlsmabs = (dt$mosq_count_pre_aitnlsmabs - dt$mosq_count_post_aitnlsmabs)/popn

dt$PERCENT_mosq_redpp_itn = (dt$mosq_count_pre_aitn - dt$mosq_count_post_aitn)/dt$mosq_count_pre_aitn
dt$PERCENT_mosq_redpp_itnlsm = (dt$mosq_count_pre_aitnlsm - dt$mosq_count_post_aitnlsm)/dt$mosq_count_pre_aitnlsm
dt$PERCENT_mosq_redpp_itnlsmabs = (dt$mosq_count_pre_aitnlsmabs - dt$mosq_count_post_aitnlsmabs)/dt$mosq_count_pre_aitnlsmabs

dt$PERCENT_mosq_redpp_itnlsmabs[dt$res == 0 & dt$phiB == 0.872]

summary(dt$prevU5_n)
summary(dt$rel_reduction_due_to_itn[dt$res == 0])
summary(dt$rel_reduction_due_to_itnlsm[dt$res == 0])
summary(dt$rel_reduction_due_to_itnlsmabs[dt$res == 0])
ddt = subset(dt, dt$res < 0.91)
range(ddt$rel_reduction_due_to_itnlsm-ddt$rel_reduction_due_to_itn)
range(ddt$rel_reduction_due_to_itnlsmabs-ddt$rel_reduction_due_to_itn)

dt$rel_reduction_due_to_itn[which(dt$res == 0 &
                                    dt$phiB == 0.872)]
dt$rel_reduction_due_to_itn[53]
par(mar = c(4,4,2,1))

plot(aa[[53]]$pv_0_1825 ~ aa[[53]]$timestep,
     type="l",ylim=c(0,0.6),
     ylab = "",
     yaxt = "n",xaxt = "n",
     xlim = c(8*365,13*365),
     xlab = "Time in years",cex.lab=1.2,cex.axis = 1.2)
axis(1, at = c(8,9,10,11,12,13,14,15)*365,
     labels = c(-1,0,1,2,3,4,5,6),cex.axis = 1.2)
axis(2, las = 2, at = seq(0,0.8,0.1), labels = seq(0,80,10),cex.axis = 1.2)
mtext("Prevalence (U5-yrs)", side = 2, line = 2.8, 
      cex = 0.85) 
text(10.5*365,0.1,"High transmission:",cex=1.4)


no_resistance = c(which(dt$res == 0))
which(dt$phiB == 0.872)

for(i in 53){
  lines(aitnlsm[[i]]$pv_0_1825 ~ aitnlsm[[i]]$timestep,col="aquamarine3",lwd=3)
  lines(aa[[i]]$pv_0_1825 ~ aa[[i]]$timestep,col="grey",lwd=3)
  lines(aitnlsmabs[[i]]$pv_0_1825 ~ aitnlsmabs[[i]]$timestep,col="orange",lwd=3,lty=2)
  lines(aitn[[i]]$pv_0_1825 ~ aitn[[i]]$timestep,col="blue",lwd=3)
}
abline(v=c(9*365,12*365),lty=2)

par(mar = c(4,4,2,1))


bars = c(dt$number_mosq_redpp_itn[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$PERCENT_mosq_redpp_itn[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$rel_reduction_allageCASES_itn[which(dt$res == 0 & dt$phiB == 0.872)],
         
         dt$number_mosq_redpp_itnlsm[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$PERCENT_mosq_redpp_itnlsm[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$rel_reduction_allageCASES_itnlsm[which(dt$res == 0 & dt$phiB == 0.872)],
         
         dt$number_mosq_redpp_itnlsmabs[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$PERCENT_mosq_redpp_itnlsmabs[which(dt$res == 0 & dt$phiB == 0.872)],
         20*dt$rel_reduction_allageCASES_itnlsmabs[which(dt$res == 0 & dt$phiB == 0.872)])

##############################
##
############################
mat <- matrix(bars, 3, 3)
rownames(mat) <- c("Reduced mosquitoes per person",
                   "% reduction mosquitoes",
                   "Relative cases reduced %")
colnames(mat) <- c("ITN","+LSM (%)","+LSM (abs)")


# plotting settings -------------------------------------------------------
ylim <- c(0,25)
angle1 <- rep(c(45,45,135), length.out=3)
angle2 <- rep(c(45,135,135), length.out=3)
density1 <- seq(5,35,length.out=3)
density2 <- seq(5,35,length.out=3)
col <- rep(c("blue","aquamarine2","orange"),each=3) # rainbow(7)


# plot --------------------------------------------------------------------
barplot(mat, beside=TRUE, ylim=ylim, 
        col=col, angle=angle1, 
        density=density1,
        ylab = "",yaxt = "n")
barplot(mat, add=TRUE, beside=TRUE, 
        ylim=ylim, col=col,
        ylab = "",yaxt = "n",
        angle=angle2, density=density2)
axis(2,las = 2, at = seq(0,20,5))
axis(4,las = 2, at = seq(0,20,4), labels = seq(0,100,20))
mtext("Number mosquito per",side = 2, line = 3, cex = 0.85) 
mtext("person reduced",side = 2, line = 2, cex = 0.85)
mtext("% cases reduced",side = 4, line = 2,cex = 0.85)
mtext("% mosquito reduced",side = 4, line = 3,cex = 0.85)


# legend("topleft", 
#        legend=c("Reduced mosquitoes per person",
#                 "% reduction mosquitoes",
#                 "Relative cases reduced %"), 
#        ncol=1, fill=TRUE,
#        angle=angle1, 
#        density=density1,
#        bty="n",cex=0.8)
# par(bg="transparent")
# legend("topleft", 
#        legend=c("Reduced mosquitoes per person",
#                 "% reduction mosquitoes",
#                 "Relative cases reduced %"), 
#        ncol=1, fill=TRUE,
#        angle=angle2, 
#        density=density2,
#        bty="n",cex=0.8)

par(mar = c(4,8,2,1))
plot(dt$rel_reduction_due_to_itn[dt$res == 0] ~ 
       dt$phiB[dt$res == 0],lwd=3,pch="",
     ylim=c(0,0.8),yaxt="n",ylab = "",
     xlab = "",xaxt="n"
)
mtext("Relative Efficacy (%)",side = 2, line = 2, cex = 0.85)
mtext("Proportion of bites in bed (%)",side = 1, line = 3, cex = 0.85)
axis(2, las = 2, at = seq(0,0.8,0.2), labels = seq(0,80,20))
axis(1, at = seq(0,1,0.2), labels = seq(0,100,20))
values = loess((dt$rel_reduction_due_to_itn[dt$res == 0] ~ 
                  dt$phiB[dt$res == 0])) 

values2 = loess((dt$rel_reduction_due_to_itnlsm[dt$res == 0] ~ 
                   dt$phiB[dt$res == 0])) 

values3 = loess((dt$rel_reduction_due_to_itnlsmabs[dt$res == 0] ~ 
                   dt$phiB[dt$res == 0])) 

lines(predict(values3) ~ dt$phiB[dt$res == 0], col = "orange",lwd=2)
lines(predict(values2) ~ dt$phiB[dt$res == 0], col = "aquamarine2",lwd=2,lty=3)
lines(predict(values) ~ dt$phiB[dt$res == 0], col = "blue",lwd=2)



par(xpd=NA,cex = 1.11)

text(x = -1.65, y = 3.52,"(A)",cex=0.8)
text(x = -0.7, y = 3.52,"(B)",cex=0.8)
text(x = 0.28, y = 3.52,"(C)",cex=0.8)
text(x = -1.65, y = 2.27,"(D)",cex=0.8)
text(x = -0.7, y = 2.27,"(E)",cex=0.8)
text(x = 0.28, y = 2.27,"(F)",cex=0.8)
text(x = -1.65, y = 0.9,"(G)",cex=0.8)
text(x = -0.7, y = 0.9,"(H)",cex=0.8)
text(x = 0.28, y = 0.9,"(I)",cex=0.8)
