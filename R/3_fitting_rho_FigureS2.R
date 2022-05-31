################################
##
## Estimate the impact of LSM from mosquito adult densities changing

## We have the % reduction given the LSM + LLIN
## and the % reduction with the LLIN

## So any addition reduction for LLIN+LSM arms would be the combination effect

## Reduction in mosquito densities

dat = read.csv("data/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE)


## number sites and which had LSM?
site = 1 #yes
site = 2 #no
site = 3 #yes
site = 4 #yes
site = 5 #no
site = 6 #no

dat$Count_gambiae_females_per_month = dat$GA.TOT - dat$GAMBIAE.males 
dat$Count_funestus_females_per_month = dat$FU.TOT - dat$FUNESTUS.males
dat$Count_culex_females_per_month = dat$CU.TOT - dat$CULEX.males 

years = unique(dat$Year)

agg_dat = function(site,dat,is_this_larvicidal_site){
  
  Count_gambiae_females_per_month = 
    Count_gambiae_total_per_month = 
    Count_funestus_females_per_month = 
    Count_funestus_total_per_month = 
    Count_culex_females_per_month = 
    Count_culex_total_per_month = array(dim=c(12,3))
  
  for(y in 1:3){
    for(m in 1:12){
      Count_gambiae_females_per_month[m,y] = sum(dat$Count_gambiae_females_per_month[dat$Month == m & dat$Year == years[y] & dat$Site.. == site],na.rm=TRUE)  
      Count_gambiae_total_per_month[m,y] =   sum(dat$GA.TOT[dat$Month == m & dat$Year == years[y] & dat$Site.. == site],na.rm=TRUE) 
      Count_funestus_females_per_month[m,y] =sum(dat$Count_funestus_females_per_month[dat$Month == m & dat$Year == years[y] & dat$Site.. == site],na.rm=TRUE) 
      Count_funestus_total_per_month[m,y] =  sum(dat$FU.TOT[dat$Month == m & dat$Year == years[y] & dat$Site.. == site],na.rm=TRUE) 
      Count_culex_females_per_month[m,y] =   sum(dat$Count_culex_females_per_month[dat$Month == m & dat$Year == years[y] & dat$Site.. == site],na.rm=TRUE) 
      Count_culex_total_per_month[m,y] =     sum(dat$CU.TOT[dat$Month == m & dat$Year == years[y] & dat$Site.. == site],na.rm=TRUE) 
    }
  }
  
  d = data.frame(d1=c(Count_gambiae_females_per_month)[1:30], 
                 d2=c(Count_gambiae_total_per_month)[1:30], 
                 d3=c(Count_funestus_females_per_month)[1:30], 
                 d4=c(Count_funestus_total_per_month)[1:30], 
                 d5=c(Count_culex_females_per_month)[1:30], 
                 d6=c(Count_culex_total_per_month)[1:30],
                 TIME = 1:30*30)
  
  yvariab = ifelse(d$d6 > d$d2,d$d6,d$d2)
  plot(yvariab ~ d$TIME,ylab="Mosquito count",xlab="Time in months",
       xaxt="n",bty="n",pch="",yaxt="n",ylim=c(0,max(yvariab)),main = is_this_larvicidal_site)
  axis(1,at=c(1*30,7*30,13*30,19*30,25*30,31*30),labels=c("Jan 2004","Jul","Jan 2005", "Jul", "Jan 2006", "Jul"))
  axis(2,las=2,at=seq(0,max(yvariab),length=6),labels=seq(0,max(yvariab),length=6))
  cols = c("purple","blue","red")
  vec = c(1,3,5)
  vec2 = c(2,4,6)
  
  polygon(c(3*30,5*30,5*30,3*30),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("grey",0.4))
  polygon(c(15*30,17*30,17*30,15*30),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("grey",0.4))
  polygon(c(27*30,29*30,29*30,27*30),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("grey",0.4))
  
  polygon(c(9*30,10*30,10*30,9*30),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("grey",0.4))
  polygon(c(21*30,22*30,22*30,21*30),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("grey",0.4))
  
  legend(22*30,max(d),legend = c("An. gambiae","An. funestus","Culex sp","Total","Females only"),
         col=c(cols,"black","black"),lty=c(1,1,1,1,2))
  
  for(i in 1:3){
    lines(d[,vec[i]] ~ d$TIME,col=cols[i],lty=2)  
    lines(d[,vec2[i]] ~ d$TIME,col=cols[i])  
  }
  
  abline(v=18*30,lty=3)
  return(list(c(Count_gambiae_females_per_month), 
              c(Count_gambiae_total_per_month), 
              c(Count_funestus_females_per_month), 
              c(Count_funestus_total_per_month), 
              c(Count_culex_females_per_month), 
              c(Count_culex_total_per_month)))
  
}

site_1_summary = agg_dat(site = 1,dat = dat,is_this_larvicidal_site = "Larviciding")
site_2_summary = agg_dat(site = 2,dat = dat,is_this_larvicidal_site = "No Larviciding")
site_3_summary = agg_dat(site = 3,dat = dat,is_this_larvicidal_site = "Larviciding")
site_4_summary = agg_dat(site = 4,dat = dat,is_this_larvicidal_site = "Larviciding")
site_5_summary = agg_dat(site = 5,dat = dat,is_this_larvicidal_site = "No Larviciding")
site_6_summary = agg_dat(site = 6,dat = dat,is_this_larvicidal_site = "No Larviciding")

relative_reduction_in_total_gambiae_fem = function(mosq_count_of_interest){
  # difference in number of mosquitoes after larviciding - diff before
  ((sum(site_2_summary[[mosq_count_of_interest]][19:30],site_5_summary[[mosq_count_of_interest]][19:30],site_6_summary[[mosq_count_of_interest]][19:30]) - 
      sum(site_1_summary[[mosq_count_of_interest]][19:30],site_3_summary[[mosq_count_of_interest]][19:30],site_4_summary[[mosq_count_of_interest]][19:30]))/
     sum(site_2_summary[[mosq_count_of_interest]][19:30],site_5_summary[[mosq_count_of_interest]][19:30],site_6_summary[[mosq_count_of_interest]][19:30]) ) -
    ((sum(site_2_summary[[mosq_count_of_interest]][1:18],site_5_summary[[mosq_count_of_interest]][1:18],site_6_summary[[mosq_count_of_interest]][1:18]) - 
        sum(site_1_summary[[mosq_count_of_interest]][1:18],site_3_summary[[mosq_count_of_interest]][1:18],site_4_summary[[mosq_count_of_interest]][1:18]))/
       sum(site_2_summary[[mosq_count_of_interest]][1:18],site_5_summary[[mosq_count_of_interest]][1:18],site_6_summary[[mosq_count_of_interest]][1:18]) ) 
  
}

relative_reduction_in_total_gambiae_fem(mosq_count_of_interest = 1)#gambiae female overall reduction
relative_reduction_in_total_gambiae_fem(mosq_count_of_interest = 2)#gambiae total overall reduction
relative_reduction_in_total_gambiae_fem(mosq_count_of_interest = 3)#funestus female overall reduction
relative_reduction_in_total_gambiae_fem(mosq_count_of_interest = 4)#funestus total overall reduction
relative_reduction_in_total_gambiae_fem(mosq_count_of_interest = 5)#culex female overall reduction
relative_reduction_in_total_gambiae_fem(mosq_count_of_interest = 6)#culex total overall reduction

trial_percentage_reduction = function(control_baseline,lsm_baseline,
                                      control_trial,lsm_year){
  # From Table 3 main publication Fillinger 2009
  100 - ((control_baseline - lsm_baseline) /
           (control_trial - lsm_year))*100 
}
trial_percentage_reduction(control_baseline = mean(c(site_2_summary[[1]][1:18],site_5_summary[[1]][1:18],site_6_summary[[1]][1:18])),
                           lsm_baseline = mean(c(site_1_summary[[1]][1:18],site_3_summary[[1]][1:18],site_4_summary[[1]][1:18])),
                           control_trial = mean(c(site_2_summary[[1]][19:30],site_5_summary[[1]][19:30],site_6_summary[[1]][19:30])),
                           lsm_year = mean(c(site_1_summary[[1]][19:30],site_3_summary[[1]][19:30],site_4_summary[[1]][19:30])))

## 61.1% for gambiae

trial_percentage_reduction(control_baseline = mean(c(site_2_summary[[3]][1:18],site_5_summary[[3]][1:18],site_6_summary[[3]][1:18])),
                           lsm_baseline = mean(c(site_1_summary[[3]][1:18],site_3_summary[[3]][1:18],site_4_summary[[3]][1:18])),
                           control_trial = mean(c(site_2_summary[[3]][19:30],site_5_summary[[3]][19:30],site_6_summary[[3]][19:30])),
                           lsm_year = mean(c(site_1_summary[[3]][19:30],site_3_summary[[3]][19:30],site_4_summary[[3]][19:30])))

## -49.5% for funestus

## This misleads us as we can see that the initial number of funestus was low for LSM arm
## Instead we want to use the functional decay to cap these estimates at 0. 

###################################################
##
## How quick does the mosquito density drop in the lsm arms


Count_gambiae_females_per_month = 
  Count_gambiae_total_per_month = 
  Count_funestus_females_per_month = 
  Count_funestus_total_per_month = 
  Count_culex_females_per_month = 
  Count_culex_total_per_month = array(dim=c(12,3,6))

site = 1:6
for(z in 1:6){
  for(y in 1:3){
    for(m in 1:12){
      Count_gambiae_females_per_month[m,y,z] = sum(dat$Count_gambiae_females_per_month[dat$Month == m & dat$Year == years[y] & dat$Site.. == site[z]],na.rm=TRUE)  
      Count_gambiae_total_per_month[m,y,z] =   sum(dat$GA.TOT[dat$Month == m & dat$Year == years[y] & dat$Site.. == site[z]],na.rm=TRUE) 
      Count_funestus_females_per_month[m,y,z] =sum(dat$Count_funestus_females_per_month[dat$Month == m & dat$Year == years[y] & dat$Site.. == site[z]],na.rm=TRUE) 
      Count_funestus_total_per_month[m,y,z] =  sum(dat$FU.TOT[dat$Month == m & dat$Year == years[y] & dat$Site.. == site[z]],na.rm=TRUE) 
      Count_culex_females_per_month[m,y,z] =   sum(dat$Count_culex_females_per_month[dat$Month == m & dat$Year == years[y] & dat$Site.. == site[z]],na.rm=TRUE) 
      Count_culex_total_per_month[m,y,z] =     sum(dat$CU.TOT[dat$Month == m & dat$Year == years[y] & dat$Site.. == site[z]],na.rm=TRUE) 
    }
  }  
}


#############################################
##
## fitting a function that tracks the reduction in the mosquitoes and estimates
## model parametesr for malariasimulation
##
#############################################

## Look at this for gambiae and funestus separately, 
## Look at this for the trial arms with and without LSM to estimate additional benefit from LSM

param_estimation_f = function(dat1, dat2, dat3, params){
  dat1 = dat1 
  dat2 = dat2
  dat3 = dat3
  
  reduction_f = function(params){ 
    
    phi <- params[1] ## malariasimulation parameter$lsm_new_eqm
    pa  <- params[2] ## malariasimulation parameter$lsm_rate_alpha
    pb  <- params[3] ## malariasimulation parameter$lsm_rate_beta
    
    times = seq(1,14*30,length=15) 
    
    ## function determined in publication
    pred1 = phi + ((1-phi) * (1 + exp(-pa))) / (1 + exp(-pa) * exp(pb * times))
    
    ## trial data from Kenyafor gambiae mosquito
    dat1 = dat1 
    dat2 = dat2
    dat3 = dat3
    
    normalized1 = c((dat1-min(dat1))/(max(dat1)-min(dat1)))
    normalized2 = c((dat2-min(dat2))/(max(dat2)-min(dat2)))
    normalized3 = c((dat3-min(dat3))/(max(dat3)-min(dat3)))
    
    ## maximum likelihood
    output = sum((normalized1-pred1)^2 + (normalized2-pred1)^2 + (normalized3-pred1)^2)
    output
  }
  
  ## parameter estimates prior to learning optimal solution 
  param_guess = c(0.141, -4, 0.1) 
  
  ## optimise the solution for sensible parameter ranges
  ## lsm_new_eqm must be between -1 and 1
  ## lsm_rate_alpha must be positive
  ## lsm_rate_beta must be negative
  satmod = optim(param_guess,lower = c(0,-6,0.001), upper = c(1,5,0.5), reduction_f,method="L-BFGS-B")
  satmod
  
  #
  ## CIs
  #
  
  size.of.grid<-40
  optim.model<-reduction_f(satmod$par)
  phi.range<-seq(0,1,length=size.of.grid)
  pa.range<-seq(-5,5,length=size.of.grid)
  pb.range<-seq(0,0.3,length=size.of.grid)
  
  ci.grid.phi<-array(0, dim = c(size.of.grid,size.of.grid,size.of.grid))
  ci.grid.pa<-array(0, dim = c(size.of.grid,size.of.grid,size.of.grid))
  ci.grid.pb<-array(0, dim = c(size.of.grid,size.of.grid,size.of.grid))
  
  for(a in 1:size.of.grid){
    for(b in 1:size.of.grid){
      for(d in 1:size.of.grid){
        p.vec<-c(phi.range[a],pa.range[b],pb.range[d])
        ci.n.param<-length(satmod$par) 
        ci.fit<-reduction_f(p.vec)     
        ci.grid.phi[a,b,d]<-ifelse(ci.fit<=optim.model+qchisq(0.8, ci.n.param)/2,phi.range[a],NA)
        ci.grid.pa[a,b,d]<-ifelse(ci.fit<=optim.model+qchisq(0.8, ci.n.param)/2,pa.range[b],NA)
        ci.grid.pb[a,b,d]<-ifelse(ci.fit<=optim.model+qchisq(0.8, ci.n.param)/2,pb.range[d],NA)
      }
    }
    print(a)
  }
  
  phi_ci<-range(ci.grid.phi, na.rm=T)##effect size
  pa_ci<-range(ci.grid.pa, na.rm=T)##per bite probability of infection from an infected mosquito to susceptible mouse
  pb_ci<-range(ci.grid.pb, na.rm=T)##per bite probability of infection from an infected mouse to susceptible mosquito 
  uncert_grid = rbind(phi_ci,pa_ci,pb_ci)
  # v<-mean(ci.grid.phi, na.rm=T);r<-mean(ci.grid.pa, na.rm=T);s<-mean(ci.grid.pb, na.rm=T)
  # means_uncert_grid = rbind(v,r,s)
  # 
  
  times2= 1:400
  
  phi <- satmod$par[1]
  pa <- satmod$par[2]
  pb <- satmod$par[3]
  pred1 = phi +  ((1-phi) * (1 + exp(-pa))) / (1 + exp(-pa) * exp(pb * times2) )
  
  pred1low = as.numeric(uncert_grid[1,1]) +  
    ((1-as.numeric(uncert_grid[1,1])) * 
       (1 + exp(-pa))) / (1 + exp(-pa) * exp(pb * times2) )
  pred1upp = as.numeric(uncert_grid[1,2]) +  
    ((1-as.numeric(uncert_grid[1,2])) * 
       (1 + exp(-pa))) / (1 + exp(-pa) * exp(pb * times2) )
  
  return(list(satmod,phi,pa,pb,pred1,
              uncert_grid,
              pred1low,pred1upp))
}
  
## trial arms with LSM
gambiae_mosq_lsm = param_estimation_f(dat1 = c(mean(c(Count_gambiae_females_per_month[1:12,1,1],Count_gambiae_females_per_month[1:5,2,1])),Count_gambiae_females_per_month[6:12,2,1],Count_gambiae_females_per_month[1:7,3,1]), 
                                        dat2 = c(mean(c(Count_gambiae_females_per_month[1:12,1,3],Count_gambiae_females_per_month[1:5,2,3])),Count_gambiae_females_per_month[6:12,2,3],Count_gambiae_females_per_month[1:7,3,3]),
                                        dat3 = c(mean(c(Count_gambiae_females_per_month[1:12,1,4],Count_gambiae_females_per_month[1:5,2,4])),Count_gambiae_females_per_month[6:12,2,4],Count_gambiae_females_per_month[1:7,3,4]),
                                        params = c(0.141, -4, 0.1) )


funestu_mosq_lsm = param_estimation_f(dat1 = c(mean(c(Count_funestus_females_per_month[1:12,1,1],Count_funestus_females_per_month[1:5,2,1])),Count_funestus_females_per_month[6:12,2,1],Count_funestus_females_per_month[1:7,3,1]),
                                        dat2 = c(mean(c(Count_funestus_females_per_month[1:12,1,3],Count_funestus_females_per_month[1:5,2,3])),Count_funestus_females_per_month[6:12,2,3],Count_funestus_females_per_month[1:7,3,3]),
                                        dat3 = c(mean(c(Count_funestus_females_per_month[1:12,1,4],Count_funestus_females_per_month[1:5,2,4])),Count_funestus_females_per_month[6:12,2,4],Count_funestus_females_per_month[1:7,3,4]),
                                        params = c(0.141, -4, 0.1) )

## trial arms without LSM
gambiae_mosq_nolsm = param_estimation_f(dat1 = c(mean(c(Count_gambiae_females_per_month[1:12,1,2],Count_gambiae_females_per_month[1:5,2,2])),Count_gambiae_females_per_month[6:12,2,2],Count_gambiae_females_per_month[1:7,3,2]), 
                                      dat2 = c(mean(c(Count_gambiae_females_per_month[1:12,1,5],Count_gambiae_females_per_month[1:5,2,5])),Count_gambiae_females_per_month[6:12,2,5],Count_gambiae_females_per_month[1:7,3,5]),
                                      dat3 = c(mean(c(Count_gambiae_females_per_month[1:12,1,6],Count_gambiae_females_per_month[1:5,2,6])),Count_gambiae_females_per_month[6:12,2,6],Count_gambiae_females_per_month[1:7,3,6]),
                                      params = c(0.141, -4, 0.1) )

funestu_mosq_nolsm = param_estimation_f(dat1 = c(mean(c(Count_funestus_females_per_month[1:12,1,2],Count_funestus_females_per_month[1:5,2,2])),Count_funestus_females_per_month[6:12,2,2],Count_funestus_females_per_month[1:7,3,2]),
                                      dat2 = c(mean(c(Count_funestus_females_per_month[1:12,1,5],Count_funestus_females_per_month[1:5,2,5])),Count_funestus_females_per_month[6:12,2,5],Count_funestus_females_per_month[1:7,3,5]),
                                      dat3 = c(mean(c(Count_funestus_females_per_month[1:12,1,6],Count_funestus_females_per_month[1:5,2,6])),Count_funestus_females_per_month[6:12,2,6],Count_funestus_females_per_month[1:7,3,6]),
                                      params = c(0.141, -4, 0.1) )

gambiae_mosq_lsm[[2]];gambiae_mosq_lsm[[3]];gambiae_mosq_lsm[[4]]
gambiae_mosq_lsm[[6]]
funestu_mosq_lsm[[2]];funestu_mosq_lsm[[3]];funestu_mosq_lsm[[4]]
funestu_mosq_lsm[[6]]
gambiae_mosq_nolsm[[2]];gambiae_mosq_nolsm[[3]];gambiae_mosq_nolsm[[4]]
gambiae_mosq_nolsm[[6]]
funestu_mosq_nolsm[[2]];funestu_mosq_nolsm[[3]];funestu_mosq_nolsm[[4]]
funestu_mosq_nolsm[[6]]





## create vector of the data
dat1 = c(mean(c(Count_gambiae_females_per_month[1:12,1,1],Count_gambiae_females_per_month[1:5,2,1])),Count_gambiae_females_per_month[6:12,2,1],Count_gambiae_females_per_month[1:7,3,1]) 
dat2 = c(mean(c(Count_gambiae_females_per_month[1:12,1,3],Count_gambiae_females_per_month[1:5,2,3])),Count_gambiae_females_per_month[6:12,2,3],Count_gambiae_females_per_month[1:7,3,3])
dat3 = c(mean(c(Count_gambiae_females_per_month[1:12,1,4],Count_gambiae_females_per_month[1:5,2,4])),Count_gambiae_females_per_month[6:12,2,4],Count_gambiae_females_per_month[1:7,3,4])
dat1b = c(mean(c(Count_funestus_females_per_month[1:12,1,1],Count_funestus_females_per_month[1:5,2,1])),Count_funestus_females_per_month[6:12,2,1],Count_funestus_females_per_month[1:7,3,1]) 
dat2b = c(mean(c(Count_funestus_females_per_month[1:12,1,3],Count_funestus_females_per_month[1:5,2,3])),Count_funestus_females_per_month[6:12,2,3],Count_funestus_females_per_month[1:7,3,3])
dat3b = c(mean(c(Count_funestus_females_per_month[1:12,1,4],Count_funestus_females_per_month[1:5,2,4])),Count_funestus_females_per_month[6:12,2,4],Count_funestus_females_per_month[1:7,3,4])

## mormalize mosquito counts between 0 and 1
normalized1 = c((dat1-min(dat1))/(max(dat1)-min(dat1)))        ## mean of 3:9 is 0.3932584,
normalized2 = c((dat2-min(dat2))/(max(dat2)-min(dat2)))        ## mean of 3:9 is 0.331384,
normalized3 = c((dat3-min(dat3))/(max(dat3)-min(dat3)))        ## mean of 3:9 is 0.2820513,
normalized1b = c((dat1b-min(dat1b))/(max(dat1b)-min(dat1b)))   ## mean of 3:9 is 0.2552,
normalized2b = c((dat2b-min(dat2b))/(max(dat2b)-min(dat2b)))   ## mean of 3:9 is 0.2552,
normalized3b = c((dat3b-min(dat3b))/(max(dat3b)-min(dat3b)))   ## mean of 3:9 is 0.2552,

times = seq(1,14*30,length=15) 

## plot the data and the reduction as observed and estimated
par(mfrow=c(1,1))
plot(normalized1 ~ times,
     ylab=expression(paste(italic("An. gambiae"),"s.l. relative densities")),
     xlab="Time in days since larvicide treatment starts",col="darkred",pch=19)
points(normalized2 ~ times,col="darkred",pch=19)
points(normalized3 ~ times,col="darkred",pch=19)

points(normalized1b ~ times,col="blue",pch=19)
points(normalized2b ~ times,col="blue",pch=19)
points(normalized3b ~ times,col="blue",pch=19)


lines(gambiae_mosq_lsm[[5]] ~ times2,lty=2,lwd=3,col="darkred")
polygon(c(times2,rev(times2)),
  c(gambiae_mosq_lsm[[7]],rev(gambiae_mosq_lsm[[8]])),
        border = NA, col=adegenet::transp("red",0.4))

lines(funestu_mosq_lsm[[5]] ~ times2,lty=2,lwd=3,col="blue")
polygon(c(times2,rev(times2)),
        c(funestu_mosq_lsm[[7]],rev(funestu_mosq_lsm[[8]])),
        border = NA, col=adegenet::transp("blue",0.4))

dat1c = c(mean(c(Count_gambiae_females_per_month[1:12,1,2],Count_gambiae_females_per_month[1:5,2,2])),Count_gambiae_females_per_month[6:12,2,2],Count_gambiae_females_per_month[1:7,3,2]) 
dat2c = c(mean(c(Count_gambiae_females_per_month[1:12,1,5],Count_gambiae_females_per_month[1:5,2,5])),Count_gambiae_females_per_month[6:12,2,5],Count_gambiae_females_per_month[1:7,3,5])
dat3c = c(mean(c(Count_gambiae_females_per_month[1:12,1,6],Count_gambiae_females_per_month[1:5,2,6])),Count_gambiae_females_per_month[6:12,2,6],Count_gambiae_females_per_month[1:7,3,6])
dat1d = c(mean(c(Count_funestus_females_per_month[1:12,1,2],Count_funestus_females_per_month[1:5,2,2])),Count_funestus_females_per_month[6:12,2,2],Count_funestus_females_per_month[1:7,3,2]) 
dat2d = c(mean(c(Count_funestus_females_per_month[1:12,1,5],Count_funestus_females_per_month[1:5,2,5])),Count_funestus_females_per_month[6:12,2,5],Count_funestus_females_per_month[1:7,3,5])
dat3d = c(mean(c(Count_funestus_females_per_month[1:12,1,6],Count_funestus_females_per_month[1:5,2,6])),Count_funestus_females_per_month[6:12,2,6],Count_funestus_females_per_month[1:7,3,6])

normalized1c = c((dat1c-min(dat1c))/(max(dat1c)-min(dat1c)))
normalized2c = c((dat2c-min(dat2c))/(max(dat2c)-min(dat2c)))
normalized3c = c((dat3c-min(dat3c))/(max(dat3c)-min(dat3c)))
normalized1d = c((dat1d-min(dat1d))/(max(dat1d)-min(dat1d)))
normalized2d = c((dat2d-min(dat2d))/(max(dat2d)-min(dat2d)))
normalized3d = c((dat3d-min(dat3d))/(max(dat3d)-min(dat3d)))



points(normalized1c ~ times,col="red",pch=8)
points(normalized2c ~ times,col="red",pch=8)
points(normalized3c ~ times,col="red",pch=8)

points(normalized1d ~ times,col="lightblue",pch=8)
points(normalized2d ~ times,col="lightblue",pch=8)
points(normalized3d ~ times,col="lightblue",pch=8)
lines(pred2gam ~ times2,lty=1,lwd=3,col="red")
lines(pred2fun ~ times2,lty=1,lwd=3,col="lightblue")

## Both trials see a reduction 
## we can work out this reduction in percentage terms
new_eqm_gamb_lsm = gambiae_mosq_lsm[[5]][200]
new_eqm_gamb_nolsm = gambiae_mosq_nolsm[[5]][200]
new_eqm_fun_lsm = funestu_mosq_lsm[[5]][200]
new_eqm_fun_nolsm = funestu_mosq_nolsm[[5]][200]

new_eqm_gamb_nolsm - new_eqm_gamb_lsm 
new_eqm_fun_nolsm - new_eqm_fun_lsm

(new_eqm_gamb_nolsm - new_eqm_gamb_lsm)/new_eqm_gamb_nolsm
(new_eqm_fun_nolsm - new_eqm_fun_lsm)/new_eqm_fun_nolsm

(as.numeric(gambiae_mosq_nolsm[[6]][1,1]) - as.numeric(gambiae_mosq_lsm[[6]][1,1]))
(as.numeric(gambiae_mosq_nolsm[[6]][1,2]) - as.numeric(gambiae_mosq_lsm[[6]][1,1]))

(as.numeric(funestu_mosq_nolsm[[6]][1,1]) - as.numeric(funestu_mosq_lsm[[6]][1,1]))
(as.numeric(funestu_mosq_nolsm[[6]][1,2]) - as.numeric(funestu_mosq_lsm[[6]][1,1]))

########################################
##
## And overall reduction

param_estimation_all_f = function(dat1, dat2, dat3, 
                                  dat1b, dat2b, dat3b, params){
  dat1 = dat1 
  dat2 = dat2
  dat3 = dat3
  dat1b = dat1b 
  dat2b = dat2b
  dat3b = dat3b
  
  reduction_f = function(params){ 
    
    phi <- params[1] ## malariasimulation parameter$lsm_new_eqm
    pa  <- params[2] ## malariasimulation parameter$lsm_rate_alpha
    pb  <- params[3] ## malariasimulation parameter$lsm_rate_beta
    
    times = seq(1,14*30,length=15) 
    
    ## function determined in publication
    pred1 = phi + ((1-phi) * (1 + exp(-pa))) / (1 + exp(-pa) * exp(pb * times))
    
    ## trial data from Kenyafor gambiae mosquito
    dat1 = dat1 
    dat2 = dat2
    dat3 = dat3
    dat1b = dat1b 
    dat2b = dat2b
    dat3b = dat3b
    
    normalized1 = c((dat1-min(dat1))/(max(dat1)-min(dat1)))
    normalized2 = c((dat2-min(dat2))/(max(dat2)-min(dat2)))
    normalized3 = c((dat3-min(dat3))/(max(dat3)-min(dat3)))
    normalized1b = c((dat1b-min(dat1b))/(max(dat1b)-min(dat1b)))
    normalized2b = c((dat2b-min(dat2b))/(max(dat2b)-min(dat2b)))
    normalized3b = c((dat3b-min(dat3b))/(max(dat3b)-min(dat3b)))
    
    ## maximum likelihood
    output = sum((normalized1-pred1)^2 + (normalized2-pred1)^2 + (normalized3-pred1)^2 + 
                   (normalized1b-pred1)^2 + (normalized2b-pred1)^2 + (normalized3b-pred1)^2)
    output
  }
  
  ## parameter estimates prior to learning optimal solution 
  param_guess = c(0.141, -4, 0.1) 
  
  ## optimise the solution for sensible parameter ranges
  ## lsm_new_eqm must be between -1 and 1
  ## lsm_rate_alpha must be positive
  ## lsm_rate_beta must be negative
  satmod = optim(param_guess,lower = c(0,-6,0.001), upper = c(1,5,0.5), reduction_f,method="L-BFGS-B")
  satmod
  
  times2= 1:400
  
  phi <- satmod$par[1]
  pa <- satmod$par[2]
  pb <- satmod$par[3]
  pred1 = phi +  ((1-phi) * (1 + exp(-pa))) / (1 + exp(-pa) * exp(pb * times2)) 
  
  return(list(satmod,phi,pa,pb,pred1))
}


all_arms_reductionLSM = param_estimation_all_f(dat1 = c(mean(c(Count_gambiae_females_per_month[1:12,1,1],Count_gambiae_females_per_month[1:5,2,1])),Count_gambiae_females_per_month[6:12,2,1],Count_gambiae_females_per_month[1:7,3,1]),
                                            dat2 = c(mean(c(Count_gambiae_females_per_month[1:12,1,3],Count_gambiae_females_per_month[1:5,2,3])),Count_gambiae_females_per_month[6:12,2,3],Count_gambiae_females_per_month[1:7,3,3]),
                                            dat3 = c(mean(c(Count_gambiae_females_per_month[1:12,1,4],Count_gambiae_females_per_month[1:5,2,4])),Count_gambiae_females_per_month[6:12,2,4],Count_gambiae_females_per_month[1:7,3,4]),
                                            dat1b = c(mean(c(Count_funestus_females_per_month[1:12,1,1],Count_funestus_females_per_month[1:5,2,1])),Count_funestus_females_per_month[6:12,2,1],Count_funestus_females_per_month[1:7,3,1]),
                                            dat2b = c(mean(c(Count_funestus_females_per_month[1:12,1,3],Count_funestus_females_per_month[1:5,2,3])),Count_funestus_females_per_month[6:12,2,3],Count_funestus_females_per_month[1:7,3,3]),
                                            dat3b = c(mean(c(Count_funestus_females_per_month[1:12,1,4],Count_funestus_females_per_month[1:5,2,4])),Count_funestus_females_per_month[6:12,2,4],Count_funestus_females_per_month[1:7,3,4]),
                                            params = c(0.141, 6, -0.057))

all_arms_reductionNoLSM = param_estimation_all_f(dat1 = c(mean(c(Count_gambiae_females_per_month[1:12,1,2],Count_gambiae_females_per_month[1:5,2,2])),Count_gambiae_females_per_month[6:12,2,2],Count_gambiae_females_per_month[1:7,3,2]), 
                                            dat2 = c(mean(c(Count_gambiae_females_per_month[1:12,1,5],Count_gambiae_females_per_month[1:5,2,5])),Count_gambiae_females_per_month[6:12,2,5],Count_gambiae_females_per_month[1:7,3,5]),
                                            dat3 = c(mean(c(Count_gambiae_females_per_month[1:12,1,6],Count_gambiae_females_per_month[1:5,2,6])),Count_gambiae_females_per_month[6:12,2,6],Count_gambiae_females_per_month[1:7,3,6]),
                                            dat1b = c(mean(c(Count_funestus_females_per_month[1:12,1,2],Count_funestus_females_per_month[1:5,2,2])),Count_funestus_females_per_month[6:12,2,2],Count_funestus_females_per_month[1:7,3,2]), 
                                            dat2b = c(mean(c(Count_funestus_females_per_month[1:12,1,5],Count_funestus_females_per_month[1:5,2,5])),Count_funestus_females_per_month[6:12,2,5],Count_funestus_females_per_month[1:7,3,5]),
                                            dat3b = c(mean(c(Count_funestus_females_per_month[1:12,1,6],Count_funestus_females_per_month[1:5,2,6])),Count_funestus_females_per_month[6:12,2,6],Count_funestus_females_per_month[1:7,3,6]),
                                            params = c(0.141, 6, -0.057))


lines(all_arms_reductionLSM[[5]] ~ times2,lty=2)
lines(all_arms_reductionNoLSM[[5]] ~ times2)

all_arms_reductionNoLSM[[5]][200] - all_arms_reductionLSM[[5]][200]


#####################################################
##
## Estimate reduction of each villages
##
##
VILL_LEVEL_f = function(dat_input){
  
  params = c(0.8, -4, 0.057) 
  param_guess = c(0.8, -4, 0.057) 
  times2= 1:400
  
  reduction_each_f = function(params){
    
    phi <- params[1]
    pa <- -4
    pb <- 0.057
    
   times = seq(1,14*30,length=15) 
    pred1 = phi + ((1-phi) * (1 + exp(-pa))) / (1 + exp(-pa) * exp(pb * times)) 
    
    dat1 = dat_input
    
    normalized1 = c((dat1-min(dat1))/(max(dat1)-min(dat1)))
    
    output = sum((normalized1-pred1)^2)
    output
  }
 
  satmod1 = optim(param_guess,
                  lower = 0,#c(0,-4.1,0), 
                  upper = 1,#c(0.8,-3.9,1), 
                  reduction_each_f,method="L-BFGS-B")
  
  phi <- satmod1$par[1]
  pa <- -4
  pb <- 0.057
  pred1vil = phi +  ((1-phi) * (1 + exp(-pa))) / (1 + exp(-pa) * exp(pb * times2)) 
  
  return(list(satmod1$par[1],
              pred1vil))
  
}


dat1 = c(mean(c(Count_gambiae_females_per_month[1:12,1,1],Count_gambiae_females_per_month[1:5,2,1])),Count_gambiae_females_per_month[6:12,2,1],Count_gambiae_females_per_month[1:7,3,1]) 
dat2 = c(mean(c(Count_gambiae_females_per_month[1:12,1,3],Count_gambiae_females_per_month[1:5,2,3])),Count_gambiae_females_per_month[6:12,2,3],Count_gambiae_females_per_month[1:7,3,3])
dat3 = c(mean(c(Count_gambiae_females_per_month[1:12,1,4],Count_gambiae_females_per_month[1:5,2,4])),Count_gambiae_females_per_month[6:12,2,4],Count_gambiae_females_per_month[1:7,3,4])
dat1b = c(mean(c(Count_funestus_females_per_month[1:12,1,1],Count_funestus_females_per_month[1:5,2,1])),Count_funestus_females_per_month[6:12,2,1],Count_funestus_females_per_month[1:7,3,1]) 
dat2b = c(mean(c(Count_funestus_females_per_month[1:12,1,3],Count_funestus_females_per_month[1:5,2,3])),Count_funestus_females_per_month[6:12,2,3],Count_funestus_females_per_month[1:7,3,3])
dat3b = c(mean(c(Count_funestus_females_per_month[1:12,1,4],Count_funestus_females_per_month[1:5,2,4])),Count_funestus_females_per_month[6:12,2,4],Count_funestus_females_per_month[1:7,3,4])
dat1c = c(mean(c(Count_gambiae_females_per_month[1:12,1,2],Count_gambiae_females_per_month[1:5,2,2])),Count_gambiae_females_per_month[6:12,2,2],Count_gambiae_females_per_month[1:7,3,2]) 
dat2c = c(mean(c(Count_gambiae_females_per_month[1:12,1,5],Count_gambiae_females_per_month[1:5,2,5])),Count_gambiae_females_per_month[6:12,2,5],Count_gambiae_females_per_month[1:7,3,5])
dat3c = c(mean(c(Count_gambiae_females_per_month[1:12,1,6],Count_gambiae_females_per_month[1:5,2,6])),Count_gambiae_females_per_month[6:12,2,6],Count_gambiae_females_per_month[1:7,3,6])
dat1d = c(mean(c(Count_funestus_females_per_month[1:12,1,2],Count_funestus_females_per_month[1:5,2,2])),Count_funestus_females_per_month[6:12,2,2],Count_funestus_females_per_month[1:7,3,2]) 
dat2d = c(mean(c(Count_funestus_females_per_month[1:12,1,5],Count_funestus_females_per_month[1:5,2,5])),Count_funestus_females_per_month[6:12,2,5],Count_funestus_females_per_month[1:7,3,5])
dat3d = c(mean(c(Count_funestus_females_per_month[1:12,1,6],Count_funestus_females_per_month[1:5,2,6])),Count_funestus_females_per_month[6:12,2,6],Count_funestus_females_per_month[1:7,3,6])

normalized1c = c((dat1c-min(dat1c))/(max(dat1c)-min(dat1c)))
normalized2c = c((dat2c-min(dat2c))/(max(dat2c)-min(dat2c)))
normalized3c = c((dat3c-min(dat3c))/(max(dat3c)-min(dat3c)))
normalized1d = c((dat1d-min(dat1d))/(max(dat1d)-min(dat1d)))
normalized2d = c((dat2d-min(dat2d))/(max(dat2d)-min(dat2d)))
normalized3d = c((dat3d-min(dat3d))/(max(dat3d)-min(dat3d)))  

normalized1 = c((dat1-min(dat1))/(max(dat1)-min(dat1)))
normalized2 = c((dat2-min(dat2))/(max(dat2)-min(dat2)))
normalized3 = c((dat3-min(dat3))/(max(dat3)-min(dat3)))
normalized1b = c((dat1b-min(dat1b))/(max(dat1b)-min(dat1b)))
normalized2b = c((dat2b-min(dat2b))/(max(dat2b)-min(dat2b)))
normalized3b = c((dat3b-min(dat3b))/(max(dat3b)-min(dat3b)))


vil1g_estimate = VILL_LEVEL_f(dat_input = dat1)
vil3g_estimate = VILL_LEVEL_f(dat_input = dat2)
vil4g_estimate = VILL_LEVEL_f(dat_input = dat3)
vil2g_estimate = VILL_LEVEL_f(dat_input = dat1c)
vil5g_estimate = VILL_LEVEL_f(dat_input = dat2c)
vil6g_estimate = VILL_LEVEL_f(dat_input = dat3c)


vil1f_estimate = VILL_LEVEL_f(dat_input = dat1b)
vil3f_estimate = VILL_LEVEL_f(dat_input = dat2b)
vil4f_estimate = VILL_LEVEL_f(dat_input = dat3b)
vil2f_estimate = VILL_LEVEL_f(dat_input = dat1d)
vil5f_estimate = VILL_LEVEL_f(dat_input = dat2d)
vil6f_estimate = VILL_LEVEL_f(dat_input = dat3d)


plot(normalized1 ~ times,
     ylab=expression(paste(italic("Anopheles")," relative densities")),
     xlab="Time in days since larvicide treatment starts",col="darkred")
points(normalized2 ~ times,col="darkred",pch=15)
points(normalized3 ~ times,col="darkred",pch=8)

points(normalized1b ~ times,col="blue",pch=0)
points(normalized2b ~ times,col="blue",pch=15)
points(normalized3b ~ times,col="blue",pch=8)

points(normalized1c ~ times,col="red",pch=0)
points(normalized2c ~ times,col="red",pch=15)
points(normalized3c ~ times,col="red",pch=8)

points(normalized1d ~ times,col="lightblue",pch=0)
points(normalized2d ~ times,col="lightblue",pch=15)
points(normalized3d ~ times,col="lightblue",pch=8)

lines(vil1g_estimate[[2]] ~ times2,lty=2,col="red")
lines(vil3g_estimate[[2]] ~ times2,lty=2,col="red")
lines(vil4g_estimate[[2]] ~ times2,lty=2,col="red")
lines(vil2g_estimate[[2]] ~ times2,col="red")
lines(vil5g_estimate[[2]] ~ times2,col="red")
lines(vil6g_estimate[[2]] ~ times2,col="red")


lines(vil1f_estimate[[2]] ~ times2,lty=2,col="blue")
lines(vil3f_estimate[[2]] ~ times2,lty=2,col="blue")
lines(vil4f_estimate[[2]] ~ times2,lty=2,col="blue")
lines(vil2f_estimate[[2]] ~ times2,col="blue")
lines(vil5f_estimate[[2]] ~ times2,col="blue")
lines(vil6f_estimate[[2]] ~ times2,col="blue")

gas = round(c(vil1g_estimate[[1]],vil2g_estimate[[1]],vil3g_estimate[[1]],vil4g_estimate[[1]],vil5g_estimate[[1]],vil6g_estimate[[1]]),2)
fus = round(c(vil1f_estimate[[1]],vil2f_estimate[[1]],vil3f_estimate[[1]],vil4f_estimate[[1]],vil5f_estimate[[1]],vil6f_estimate[[1]]),2)

mean(c(vil2g_estimate[[1]],vil5g_estimate[[1]],vil6g_estimate[[1]]))
mean(c(vil2f_estimate[[1]],vil5f_estimate[[1]],vil6f_estimate[[1]]))

## Global No LSM versus local LSM
## GAMB
gasLSM = round(c(vil1g_estimate[[1]],vil3g_estimate[[1]],vil4g_estimate[[1]]),2)
gasNoLSM = round(c(vil2g_estimate[[1]],vil5g_estimate[[1]],vil6g_estimate[[1]]),2)
## FUNESTUS
fusLSM = round(c(vil1f_estimate[[1]],vil3f_estimate[[1]],vil4f_estimate[[1]]),2)
fusNoLSM = round(c(vil2f_estimate[[1]],vil5f_estimate[[1]],vil6f_estimate[[1]]),2)

(gasNoLSM-gasLSM)/gasNoLSM
(fusNoLSM-fusLSM)/fusNoLSM

#########################################################################
##
## Using the ratio of mosquitoes caught prior to the baseline 
## for species composition
##
#################################################################

species_composition = array(dim=c(6,6))
for(i in 1:6){
  ## the total number of gambiae prior to interventions
  species_composition[i,1] = sum(c(Count_gambiae_females_per_month[1:12,1,i],Count_gambiae_females_per_month[1:5,2,i]))
  
  ## the total number of funestus prior to interventions
  species_composition[i,2] = sum(c(Count_funestus_females_per_month[1:12,1,i],Count_funestus_females_per_month[1:5,2,i]))
  
}

## proportion of mosquito that are gambiae per village (6 villages)
species_composition[,3] = species_composition[,1]/(species_composition[,1] + species_composition[,2])

## of which were arabiensis
species_composition[,4] = species_composition[,3] * 0.03 ## arabiensis
## of which were gambiae
species_composition[,5] = species_composition[,3] * 0.97 ## gambiae ss
## of which are funestus
species_composition[,6] = 1 - species_composition[,3]    ## funestus

village_level_average_gambiae_lsm_reduction = summary(c((vil2g_estimate[[2]][200] - vil1g_estimate[[2]][200])/vil2g_estimate[[2]][200],
       (vil5g_estimate[[2]][200] - vil1g_estimate[[2]][200])/vil5g_estimate[[2]][200],
       (vil6g_estimate[[2]][200] - vil1g_estimate[[2]][200])/vil6g_estimate[[2]][200],
       
       (vil2g_estimate[[2]][200] - vil3g_estimate[[2]][200])/vil2g_estimate[[2]][200],
       (vil5g_estimate[[2]][200] - vil3g_estimate[[2]][200])/vil5g_estimate[[2]][200],
       (vil6g_estimate[[2]][200] - vil3g_estimate[[2]][200])/vil6g_estimate[[2]][200],
       
       (vil2g_estimate[[2]][200] - vil4g_estimate[[2]][200])/vil2g_estimate[[2]][200],
       (vil5g_estimate[[2]][200] - vil4g_estimate[[2]][200])/vil5g_estimate[[2]][200],
       (vil6g_estimate[[2]][200] - vil4g_estimate[[2]][200])/vil6g_estimate[[2]][200]))

village_level_average_funestus_lsm_reduction = summary(c((vil2f_estimate[[2]][200] - vil1f_estimate[[2]][200])/vil2f_estimate[[2]][200],
       (vil5f_estimate[[2]][200] - vil1f_estimate[[2]][200])/vil5f_estimate[[2]][200],
       (vil6f_estimate[[2]][200] - vil1f_estimate[[2]][200])/vil6f_estimate[[2]][200],
       
       (vil2f_estimate[[2]][200] - vil3f_estimate[[2]][200])/vil2f_estimate[[2]][200],
       (vil5f_estimate[[2]][200] - vil3f_estimate[[2]][200])/vil5f_estimate[[2]][200],
       (vil6f_estimate[[2]][200] - vil3f_estimate[[2]][200])/vil6f_estimate[[2]][200],
       
       (vil2f_estimate[[2]][200] - vil4f_estimate[[2]][200])/vil2f_estimate[[2]][200],
       (vil5f_estimate[[2]][200] - vil4f_estimate[[2]][200])/vil5f_estimate[[2]][200],
       (vil6f_estimate[[2]][200] - vil4f_estimate[[2]][200])/vil6f_estimate[[2]][200]))

village_level_average_gambiae_lsm_reduction
village_level_average_funestus_lsm_reduction

######################################
##
## Parameter estimates for the transmission model (malariasimulation)
##
######################################

# LSM vs No LSM (grouping villages)
# this is the reduction on top of mosquito nets
lsm_new_eqm = 