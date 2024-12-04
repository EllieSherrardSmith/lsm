############################################################
##
## Predicting each village in the trial 
## 
##

## model function

## Prevalence data
## microscopy data
data_people = read.csv("data/Fillinger_Kenya highland LSM_MALARIA SURVEY DATA.csv",header=TRUE) 


sites = 1:6
months = c("January","April", "May","June","July","August","November","December")
years = c("2004","2005","2006","2007")

Count_prev_sites = Count_tote_sites = 
  treated = Amodiaquine_used = Chloroquine_used = Quinine_used = ACT_used = SP_used = 
  nets = sprays = coils = array(dim=c(8,4,6))

for(s in 1:length(sites)){
  for(m in 1:length(months)){
    for(y in 1:length(years)){
      Count_prev_sites[m,y,s] = length(data_people$Microscopy.result.negative..0..or.positive..1.[data_people$Site.code == sites[s] &
                                                                                                    data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                                                                    data_people$Microscopy.result.negative..0..or.positive..1. == 1])      
      Count_tote_sites[m,y,s] = length(data_people$Microscopy.result.negative..0..or.positive..1.[data_people$Site.code == sites[s] &
                                                                                                    data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                                                                    data_people$Microscopy.result.negative..0..or.positive..1. == 0])      
      treated[m,y,s] = length(data_people$Malaria.Treatment.in.last.7.days[data_people$Site.code == sites[s] &
                                                                                                    data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                                                                    data_people$Malaria.Treatment.in.last.7.days == 1])      
      ACT_used[m,y,s] = length(data_people$ACT[data_people$Site.code == sites[s] &
                                                 data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                 data_people$ACT == 1])      
      SP_used[m,y,s] = length(data_people$SP[data_people$Site.code == sites[s] &
                                                 data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                 data_people$SP == 1])      
      Amodiaquine_used[m,y,s] = length(data_people$Amodiaquine[data_people$Site.code == sites[s] &
                                                 data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                 data_people$Amodiaquine == 1])      
      Chloroquine_used[m,y,s] = length(data_people$Chloroquine[data_people$Site.code == sites[s] &
                                                                 data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                                 data_people$Chloroquine == 1])      
      Quinine_used[m,y,s] = length(data_people$Qinine[data_people$Site.code == sites[s] &
                                                                 data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                                 data_people$Qinine == 1])      
      nets[m,y,s] = length(data_people$Screened.child.slept.under.ITN.last.night..No..0..or.yes..1.[data_people$Site.code == sites[s] &
                                                        data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                        data_people$Screened.child.slept.under.ITN.last.night..No..0..or.yes..1. == 1])      
      sprays[m,y,s] = length(data_people$Insecticide.Spray.used.last.night..No..0..or.yes..1.[data_people$Site.code == sites[s] &
                                                        data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                        data_people$Insecticide.Spray.used.last.night..No..0..or.yes..1. == 1])      
      coils[m,y,s] = length(data_people$Mosquito.coil.used.last.night..No..0..or.yes..1.[data_people$Site.code == sites[s] &
                                                        data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                        data_people$Mosquito.coil.used.last.night..No..0..or.yes..1. == 1])      
    }
  }
}

total_people = ifelse((Count_prev_sites+Count_tote_sites) == 0, -999, (Count_prev_sites+Count_tote_sites))

net_users = nets/total_people

sites1_nets = as.numeric(c(net_users[1,1,1],NA,NA,net_users[2:6,1,1],NA,NA,net_users[7:8,1,1],
                              net_users[1,2,1],NA,NA,net_users[2:6,2,1],NA,NA,net_users[7:8,2,1],
                              net_users[1,3,1],NA,NA,net_users[2:6,3,1],NA,NA,net_users[7:8,3,1],
                              net_users[1,4,1],NA,NA,net_users[2:6,4,1],NA,NA,net_users[7:8,4,1]))
sites2_nets = as.numeric(c(net_users[1,1,2],NA,NA,net_users[2:6,1,2],NA,NA,net_users[7:8,1,2],
                           net_users[1,2,2],NA,NA,net_users[2:6,2,2],NA,NA,net_users[7:8,2,2],
                           net_users[1,3,2],NA,NA,net_users[2:6,3,2],NA,NA,net_users[7:8,3,2],
                           net_users[1,4,2],NA,NA,net_users[2:6,4,2],NA,NA,net_users[7:8,4,2]))
sites3_nets = as.numeric(c(net_users[1,1,3],NA,NA,net_users[2:6,1,3],NA,NA,net_users[7:8,1,3],
                           net_users[1,2,3],NA,NA,net_users[2:6,2,3],NA,NA,net_users[7:8,2,3],
                           net_users[1,3,3],NA,NA,net_users[2:6,3,3],NA,NA,net_users[7:8,3,3],
                           net_users[1,4,3],NA,NA,net_users[2:6,4,3],NA,NA,net_users[7:8,4,3]))
sites4_nets = as.numeric(c(net_users[1,1,4],NA,NA,net_users[2:6,1,4],NA,NA,net_users[7:8,1,4],
                           net_users[1,2,4],NA,NA,net_users[2:6,2,4],NA,NA,net_users[7:8,2,4],
                           net_users[1,3,4],NA,NA,net_users[2:6,3,4],NA,NA,net_users[7:8,3,4],
                           net_users[1,4,4],NA,NA,net_users[2:6,4,4],NA,NA,net_users[7:8,4,4]))
sites5_nets = as.numeric(c(net_users[1,1,5],NA,NA,net_users[2:6,1,5],NA,NA,net_users[7:8,1,5],
                           net_users[1,2,5],NA,NA,net_users[2:6,2,5],NA,NA,net_users[7:8,2,5],
                           net_users[1,3,5],NA,NA,net_users[2:6,3,5],NA,NA,net_users[7:8,3,5],
                           net_users[1,4,5],NA,NA,net_users[2:6,4,5],NA,NA,net_users[7:8,4,5]))
sites6_nets = as.numeric(c(net_users[1,1,6],NA,NA,net_users[2:6,1,6],NA,NA,net_users[7:8,1,6],
                           net_users[1,2,6],NA,NA,net_users[2:6,2,6],NA,NA,net_users[7:8,2,6],
                           net_users[1,3,6],NA,NA,net_users[2:6,3,6],NA,NA,net_users[7:8,3,6],
                           net_users[1,4,6],NA,NA,net_users[2:6,4,6],NA,NA,net_users[7:8,4,6]))

summarydat = data.frame(months = rep(c("January",NA,NA,"April", "May","June","July","August",NA,NA,"November","December"),4),
                        year = rep(c(2004:2007),each=12),
                        site1nets = sites1_nets,
                        site2nets = sites2_nets,
                        site3nets = sites3_nets,
                        site4nets = sites4_nets,
                        site5nets = sites5_nets,
                        site6nets = sites6_nets,
                        timesteps_list = 1:48,
                        timesteps = round(seq(731,2161,length=48)))
write.csv(summarydat,"data/net_use.csv")

plot(summarydat$site1nets ~ c(1:48))
abline(v=18,lty=2)
plot(summarydat$site2nets ~ c(1:48))
abline(v=18,lty=2)
plot(summarydat$site3nets ~ c(1:48))
abline(v=18,lty=2)
plot(summarydat$site4nets ~ c(1:48))
abline(v=18,lty=2)
plot(summarydat$site5nets ~ c(1:48))
abline(v=18,lty=2)
plot(summarydat$site6nets ~ c(1:48))
abline(v=18,lty=2)

mean(sites1_nets[1:18],na.rm=TRUE)
mean(sites1_nets[19:48],na.rm=TRUE)


total_treated = 100 * treated/Count_prev_sites
total_treated_Amo = 100 * Amodiaquine_used/treated
total_treated_SP = 100 * SP_used/treated
total_treated_Chloro = 100 * Chloroquine_used/treated



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


#Mean baseline prevalence smoothed##
range(c(sites1[1:18]),na.rm=TRUE)
range(c(sites2[1:18]),na.rm=TRUE)
range(c(sites3[1:18]),na.rm=TRUE)
range(c(sites4[1:18]),na.rm=TRUE)
range(c(sites5[1:18]),na.rm=TRUE)
range(c(sites6[1:18]),na.rm=TRUE)

#Mean baseline at jan 2005
mean(c(sites1[1:18]),na.rm=TRUE)
mean(c(sites2[1:18]),na.rm=TRUE)
mean(c(sites3[1:18]),na.rm=TRUE)
mean(c(sites4[1:18]),na.rm=TRUE)
mean(c(sites5[1:18]),na.rm=TRUE)
mean(c(sites6[1:18]),na.rm=TRUE)


## Reduction in mosquito densities

# dat 
dat = read.csv("data/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE)

summary(dat)


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


site_specific_relative_reduction_in_total_anopheles = function(site_summary,mosq_count_of_interest){
  # difference in number of mosquitoes after larviciding - diff before
  ((sum(site_summary[[mosq_count_of_interest]][1:18]) - 
      sum(site_summary[[mosq_count_of_interest]][19:30]))/
     sum(site_summary[[mosq_count_of_interest]][1:18]))
  
}
site_specific_relative_reduction_in_total_anopheles(site_1_summary,1)
site_specific_relative_reduction_in_total_anopheles(site_1_summary,3)

site_specific_relative_reduction_in_total_anopheles(site_2_summary,1)
site_specific_relative_reduction_in_total_anopheles(site_2_summary,3)

site_specific_relative_reduction_in_total_anopheles(site_3_summary,1)
site_specific_relative_reduction_in_total_anopheles(site_3_summary,3)

site_specific_relative_reduction_in_total_anopheles(site_4_summary,1)
site_specific_relative_reduction_in_total_anopheles(site_4_summary,3)

site_specific_relative_reduction_in_total_anopheles(site_5_summary,1)
site_specific_relative_reduction_in_total_anopheles(site_5_summary,3)

site_specific_relative_reduction_in_total_anopheles(site_6_summary,1)
site_specific_relative_reduction_in_total_anopheles(site_6_summary,3)



## match up the seasonality curves
##BIMODAL
TIME = 1:365

theta_c =  0.2866762 
ssa0 = 	0.285277
ssa1 =	0.1
ssb1 =	0.1
ssa2 =	-0.140318 #increases / decreases height
ssb2 = 	-0.0730962
ssa3 =  -0.016891
ssb3 =	0

data = (ssa0+ssa1*cos(2*pi*TIME/365)+ssa2*cos(2*2*pi*TIME/365)+ssa3*cos(3*2*pi*TIME/365)+ssb1*sin(2*pi*TIME/365)+ssb2*sin(2*2*pi*TIME/365)+ ssb3*sin(3*2*pi*TIME/365) )
data2 = rep(data,3)[1:912]
TIME2 = 1:912
lines(data2*200 ~ TIME2)
lines(data2[540:900]*200 ~ TIME2[540:900])



#############################################
##
## Optimising this fit with Optim
## 

## Create strings of yearly data
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
## Where larviciding was completed
#site = 1 #yes (just use year 1)
#site = 2 #no
#site = 3 #yes (just use year 1)
#site = 4 #yes (just use year 1)
#site = 5 #no
#site = 6 #no
plot(Count_gambiae_females_per_month[,1,6]/100) #(using the data with enough counts)


fourier_f<-function(param){
  
  ssa0 <- param[1]
  ssa1 <- param[2]
  ssb1 <- param[3]
  ssa2 <- param[4]
  ssb2 <- param[5]
  ssa3 <- param[6]
  ssb3 <- param[7]
  
  TIME = seq(1,365,length=12)
  pred1<- (ssa0+ssa1*cos(2*pi*TIME/365)+ssa2*cos(2*2*pi*TIME/365)+ssa3*cos(3*2*pi*TIME/365)+ssb1*sin(2*pi*TIME/365)+ssb2*sin(2*2*pi*TIME/365)+ ssb3*sin(3*2*pi*TIME/365) )
  
  data1<- Count_gambiae_females_per_month[,1,1]/100
  data2<- Count_gambiae_females_per_month[,1,3]/100
  data3<- Count_gambiae_females_per_month[,1,4]/100
  
  data4<- Count_gambiae_females_per_month[,1,2]/100
  data5<- Count_gambiae_females_per_month[,1,6]/100
  data6<- Count_gambiae_females_per_month[,2,6]/100
  
  output = sum((data1 - pred1)^2 + (data2 - pred1)^2 + (data3 - pred1)^2 + (data4 - pred1)^2 + (data5 - pred1)^2 + (data6 - pred1)^2)
  output
}

param_guess <- c(0.285277,0.1,0.1,-0.140318,-0.0730962,-0.016891,0)

satmod <- optim(param_guess,fourier_f,method="L-BFGS-B")


ssa0 <- 0.218087682 
ssa1 <- -0.005292592  
ssb1 <- 0.174739899 
ssa2 <- -0.085277539 
ssb2 <- -0.082337283  
ssa3 <- 0.017356449  
ssb3 <- 0.026755121



###################################
## Match up the seasonality with the mosquito density for female An gambiae
##site 1
match_season_function = function(site,larvicide,is_this_larvicidal_site){
  
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
  plot(d$d1 ~ d$TIME,ylab="Mosquito count",xlab="Time in months",
       xaxt="n",bty="n",pch="",yaxt="n",ylim=c(0,max(d$d2)),main = is_this_larvicidal_site)
  axis(1,at=c(1*30,7*30,13*30,19*30,25*30,31*30),labels=c("Jan 2004","Jul","Jan 2005", "Jul", "Jan 2006", "Jul"))
  axis(2,las=2,at=seq(0,max(d$d2),length=6),labels=seq(0,max(d$d2),length=6))
  cols = c("purple","blue","red")
  vec = c(1,3,5)
  vec2 = c(2,4,6)
  
  polygon(c(3*30,5*30,5*30,3*30),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("grey",0.4))
  polygon(c(15*30,17*30,17*30,15*30),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("grey",0.4))
  polygon(c(27*30,29*30,29*30,27*30),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("grey",0.4))
  
  polygon(c(9*30,10*30,10*30,9*30),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("grey",0.4))
  polygon(c(21*30,22*30,22*30,21*30),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("grey",0.4))
  
  #legend(22*30,max(d$d2),legend = c("An. gambiae","An. funestus","Culex sp","Total","Females only"),
  #       col=c(cols,"black","black"),lty=c(1,1,1,1,2),bty="n")
  
  for(i in 1){
    lines(d[,vec[i]] ~ d$TIME,col=cols[i],lty=2)  
    lines(d[,vec2[i]] ~ d$TIME,col=cols[i])  
  }
  
  calibr = quantile(d$d1[1:18],0.9)
  larvicide = larvicide
  calibr2 = ifelse(larvicide == "yes",calibr*0.2,calibr)
  ##BIMODAL
  TIME = 1:365
  
  theta_c =  0.2122221
  ssa0 <- 0.218087682 
  ssa1 <- -0.005292592  
  ssb1 <- 0.174739899 
  ssa2 <- -0.085277539 #increases/decreases height
  ssb2 <- -0.082337283  
  ssa3 <- 0.017356449  
  ssb3 <- 0.026755121
  
  data = (ssa0+ssa1*cos(2*pi*TIME/365)+ssa2*cos(2*2*pi*TIME/365)+ssa3*cos(3*2*pi*TIME/365)+ssb1*sin(2*pi*TIME/365)+ssb2*sin(2*2*pi*TIME/365)+ ssb3*sin(3*2*pi*TIME/365) )
  data2 = rep(data,3)[1:912]
  TIME2 = 1:912
  
  lines(data2[1:540]*calibr ~ TIME2[1:540],lty=2,lwd=2,col="grey20")
  lines(data2[540:900]*calibr2 ~ TIME2[540:900],lty=2,lwd=2,col="grey20") #80% reduction
  abline(v=540,lty=4)
  
}
par(mfrow=c(2,3))
match_season_function(site = 1,larvicide = "yes",is_this_larvicidal_site = "Larviciding")
match_season_function(site = 3,larvicide = "yes",is_this_larvicidal_site = "Larviciding")
match_season_function(site = 4,larvicide = "yes",is_this_larvicidal_site = "Larviciding")
match_season_function(site = 2,larvicide = "no",is_this_larvicidal_site = "No larviciding")
match_season_function(site = 5,larvicide = "no",is_this_larvicidal_site = "No larviciding")
match_season_function(site = 6,larvicide = "no",is_this_larvicidal_site = "No larviciding")
legend(200,66,legend=c("Seasonal pattern"),lwd=2,lty=2,col="grey20",bty="n")
legend(22*30,66,legend = c("An. gambiae total","Females"),
       col=c("purple","purple"),lty=c(1,2),bty="n")


## Net use data

sites = 1:6
months = c("January","April", "May","June","July","August","November","December")
years = c("2004","2005","2006","2007")
ages = 1:4
## by age group toofor(s in 1:length(sites)){
Count_used_sitesA = Count_asked_nets_sitesA = array(dim=c(8,4,6,4))
for(s in 1:length(sites)){
  for(m in 1:length(months)){
    for(y in 1:length(years)){
      for(a in 1:4){
        Count_used_sitesA[m,y,s,a] = length(data_people$Screened.child.slept.under.ITN.last.night..No..0..or.yes..1.[data_people$Age.group.code == ages[a] &
                                                                                                                       data_people$Site.code == sites[s] &
                                                                                                                       data_people$MONTH == months[m] & data_people$YEAR == years[y] &
                                                                                                                       data_people$Screened.child.slept.under.ITN.last.night..No..0..or.yes..1. == 1])      
        
        
        Count_asked_nets_sitesA[m,y,s,a] = length(data_people$Screened.child.slept.under.ITN.last.night..No..0..or.yes..1.[data_people$Age.group.code == ages[a] &
                                                                                                                             data_people$Site.code == sites[s] &
                                                                                                                             data_people$MONTH == months[m] & data_people$YEAR == years[y]])      
        
      }
    }
  }
}
Net_useA = ifelse(Count_asked_nets_sitesA == 0, -999,Count_used_sitesA / Count_asked_nets_sitesA)
Net_useA[Net_useA == -999] <- NA   


sites_age = n_people = array(dim=c(3,6,4))
for(i in 1:6){
  for(a in 1:4){
    for(j in 1:3){
      ## Before
      sites_age[j,i,a] = mean(c(Net_useA[5:8,j,i,a],Net_useA[1:4,j+1,i,a]),na.rm=TRUE) 
      n_people[j,i,a] = sum(c(Count_asked_nets_sitesA[5:8,j,i,a],Count_asked_nets_sitesA[1:4,j+1,i,a]))
    }    
  }
}

## weighted average net use by age
n_people_per_yr = array(dim=c(3,6))
for(i in 1:6){
  for(j in 1:3){
    n_people_per_yr[j,i] = sum(n_people[j,i,1],n_people[j,i,2],n_people[j,i,3],n_people[j,i,4])

  }
}
ratio_people = array(dim=c(3,6,4))
for(i in 1:6){
    for(j in 1:3){
      ratio_people[j,i,] = n_people[j,i,]/n_people_per_yr[j,i]
      
    }
}
weighted_ave_net_use = array(dim=c(3,6))
for(i in 1:6){
  for(j in 1:3){
    weighted_ave_net_use[j,i] = sum(sites_age[j,i,1]*ratio_people[j,i,1] + sites_age[j,i,2]*ratio_people[j,i,2] + 
                                      sites_age[j,i,3]*ratio_people[j,i,3] + sites_age[j,i,4]*ratio_people[j,i,4])
  }  
  }


Count_used_sites = Count_asked_nets_sites = 
  Count_ACT_neg_sites = Count_ACT_pos_sites = array(dim=c(8,4,6))

for(s in 1:length(sites)){
  for(m in 1:length(months)){
    for(y in 1:length(years)){
      Count_used_sites[m,y,s] = length(data_people$Screened.child.slept.under.ITN.last.night..No..0..or.yes..1.[data_people$Site.code == sites[s] &
                                                                                                                  data_people$MONTH == months[m] & data_people$YEAR == years[y] &
                                                                                                               data_people$Screened.child.slept.under.ITN.last.night..No..0..or.yes..1. == 1])      
      Count_asked_nets_sites[m,y,s] = length(data_people$Screened.child.slept.under.ITN.last.night..No..0..or.yes..1.[data_people$Site.code == sites[s] &
                                                                                                                 data_people$MONTH == months[m] & data_people$YEAR == years[y]])      
      Count_ACT_pos_sites[m,y,s] = length(data_people$ACT[data_people$Site.code == sites[s] &
                                                        data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                        data_people$Microscopy.result.negative..0..or.positive..1. == 1])      
      Count_ACT_neg_sites[m,y,s] = length(data_people$ACT[data_people$Site.code == sites[s] &
                                                        data_people$MONTH == months[m] & data_people$YEAR == years[y]  &
                                                        data_people$Microscopy.result.negative..0..or.positive..1. == 0]) 
      }
  }
}


Net_use = ifelse(Count_asked_nets_sites == 0, -999,Count_used_sites / Count_asked_nets_sites)
Net_use[Net_use == -999] <- NA

rownames(Net_use) = c("January","April", "May","June","July","August","November","December")
colnames(Net_use) = c("2004","2005","2006","2007")


sites1 = 100*as.numeric(c(Net_use[1,1,1],NA,NA,Net_use[2:6,1,1],NA,NA,Net_use[7:8,1,1],
                      Net_use[1,2,1],NA,NA,Net_use[2:6,2,1],NA,NA,Net_use[7:8,2,1],
                      Net_use[1,3,1],NA,NA,Net_use[2:6,3,1],NA,NA,Net_use[7:8,3,1],
                      Net_use[1,4,1],NA,NA,Net_use[2:6,4,1],NA,NA,Net_use[7:8,4,1]))
sites2 =  100*as.numeric(c(Net_use[1,1,2],NA,NA,Net_use[2:6,1,2],NA,NA,Net_use[7:8,1,2],
                       Net_use[1,2,2],NA,NA,Net_use[2:6,2,2],NA,NA,Net_use[7:8,2,2],
                       Net_use[1,3,2],NA,NA,Net_use[2:6,3,2],NA,NA,Net_use[7:8,3,2],
                       Net_use[1,4,2],NA,NA,Net_use[2:6,4,2],NA,NA,Net_use[7:8,4,2]))
sites3 =  100*as.numeric(c(Net_use[1,1,3],NA,NA,Net_use[2:6,1,3],NA,NA,Net_use[7:8,1,3],
                       Net_use[1,2,3],NA,NA,Net_use[2:6,2,3],NA,NA,Net_use[7:8,2,3],
                       Net_use[1,3,3],NA,NA,Net_use[2:6,3,3],NA,NA,Net_use[7:8,3,3],
                       Net_use[1,4,3],NA,NA,Net_use[2:6,4,3],NA,NA,Net_use[7:8,4,3]))
sites4 =  100*as.numeric(c(Net_use[1,1,4],NA,NA,Net_use[2:6,1,4],NA,NA,Net_use[7:8,1,4],
                       Net_use[1,2,4],NA,NA,Net_use[2:6,2,4],NA,NA,Net_use[7:8,2,4],
                       Net_use[1,3,4],NA,NA,Net_use[2:6,3,4],NA,NA,Net_use[7:8,3,4],
                       Net_use[1,4,4],NA,NA,Net_use[2:6,4,4],NA,NA,Net_use[7:8,4,4]))
sites5 =  100*as.numeric(c(Net_use[1,1,5],NA,NA,Net_use[2:6,1,5],NA,NA,Net_use[7:8,1,5],
                       Net_use[1,2,5],NA,NA,Net_use[2:6,2,5],NA,NA,Net_use[7:8,2,5],
                       Net_use[1,3,5],NA,NA,Net_use[2:6,3,5],NA,NA,Net_use[7:8,3,5],
                       Net_use[1,4,5],NA,NA,Net_use[2:6,4,5],NA,NA,Net_use[7:8,4,5]))
sites6 =  100*as.numeric(c(Net_use[1,1,6],NA,NA,Net_use[2:6,1,6],NA,NA,Net_use[7:8,1,6],
                       Net_use[1,2,6],NA,NA,Net_use[2:6,2,6],NA,NA,Net_use[7:8,2,6],
                       Net_use[1,3,6],NA,NA,Net_use[2:6,3,6],NA,NA,Net_use[7:8,3,6],
                       Net_use[1,4,6],NA,NA,Net_use[2:6,4,6],NA,NA,Net_use[7:8,4,6]))



timeline = 1:(12*4)
dat_nets = read.csv("C:/Users/esherrar/Documents/larval_source_management/Data resource/mosquito net use per month per arm.csv",header=TRUE) 
head(dat_nets)
par(mfrow=c(1,1))
plot(dat_nets$Site.1,ylim=c(0,100),pch="",col="orange",
     xaxt="n",
     yaxt="n",ylab="Net use reported as screened child slept under net (%)",
     xlab = "Date")
axis(1,at=c(1,7,13,19,25,31),labels=c("Jun 2004","Jan 2005",
                                      "Jun 2005","Jan 2006",
                                      "Jun 2006","Jan 2007"))
axis(2,las=2,at=c(0,20,40,60,80,100))
# lines(dat_nets$Site.3,col="orange")
# lines(dat_nets$Site.4,col="orange")
# lines(dat_nets$Site.2,col="darkgreen")
# lines(dat_nets$Site.5,col="darkgreen")
# lines(dat_nets$Site.6,col="darkgreen")

points(sites1[6:48] ~ timeline[6:48],col="orange",pch=15)
points(sites2[6:48] ~ timeline[6:48],col="darkgreen",pch=15)
points(sites3[6:48] ~ timeline[6:48],col="orange",pch=8)
points(sites4[6:48] ~ timeline[6:48],col="orange",pch=5)
points(sites5[6:48] ~ timeline[6:48],col="darkgreen",pch=5)
points(sites6[6:48] ~ timeline[6:48],col="darkgreen",pch=8)

abline(v=13,lty=2)

mean(sites1[1:12],na.rm=TRUE)
mean(sites2[1:12],na.rm=TRUE)
mean(sites3[1:12],na.rm=TRUE)
mean(sites4[1:12],na.rm=TRUE)
mean(sites5[1:12],na.rm=TRUE)
mean(sites6[1:12],na.rm=TRUE)

mean(sites1[13:24],na.rm=TRUE)
mean(sites2[13:24],na.rm=TRUE)
mean(sites3[13:24],na.rm=TRUE)
mean(sites4[13:24],na.rm=TRUE)
mean(sites5[13:24],na.rm=TRUE)
mean(sites6[13:24],na.rm=TRUE)

mean(sites1[25:48],na.rm=TRUE)
mean(sites2[25:48],na.rm=TRUE)
mean(sites3[25:48],na.rm=TRUE)
mean(sites4[25:48],na.rm=TRUE)
mean(sites5[25:48],na.rm=TRUE)
mean(sites6[25:48],na.rm=TRUE)


##################################################
##
## Using default seasonality to avoid overfitting
##
##################################################

## Need to offset so that nets are introduced in June each year



# Obtain fourier coefficients

## Site is

data = read.csv("C:/Users/esherrar/Documents/larval_source_management/Data resource/Modelling/Site parameters.csv",header=TRUE)

seasonal_fourier_coefs_shift = array(dim=c(nrow(data),7))
offset = 0.5

for(i in 1:nrow(data)){
  
  ##Re-ordering to match the site file we already have
  seasonal_fourier_coefs_shift[i,1] = data$seasonal_a0[i]
  
  seasonal_fourier_coefs_shift[i,2] = data$seasonal_a1[i] * cos(1*2*pi*offset) - data$seasonal_b1[i] * sin(1*2*pi*offset)
  seasonal_fourier_coefs_shift[i,3] = data$seasonal_a1[i] * sin(1*2*pi*offset) + data$seasonal_b1[i] * cos(1*2*pi*offset)
  seasonal_fourier_coefs_shift[i,4] = data$seasonal_a2[i] * cos(2*2*pi*offset) - data$seasonal_b2[i] * sin(2*2*pi*offset)
  seasonal_fourier_coefs_shift[i,5] = data$seasonal_a2[i] * sin(2*2*pi*offset) + data$seasonal_b2[i] * cos(2*2*pi*offset)
  seasonal_fourier_coefs_shift[i,6] = data$seasonal_a3[i] * cos(3*2*pi*offset) - data$seasonal_b3[i] * sin(3*2*pi*offset)
  seasonal_fourier_coefs_shift[i,7] = data$seasonal_a3[i] * sin(3*2*pi*offset) + data$seasonal_b3[i] * cos(3*2*pi*offset)
  
}

colnames(seasonal_fourier_coefs_shift) = names(data[5:11])
head(seasonal_fourier_coefs_shift)
head(data[,5:11])

data[,5:11] = seasonal_fourier_coefs_shift
write.csv(data,"C:/Users/esherrar/Documents/larval_source_management/Data resource/Modelling/Site parameters offset.csv")
##


input_params = expand.grid(uncertainty_ID = 1:100,sites = 1:6)
input_params$prevalence_baseline = c(rnorm(n=100,mean=63.31,sd=5),
                                     rnorm(n=100,mean=60.14,sd=5),
                                     rnorm(n=100,mean=55.81,sd=5),
                                     rnorm(n=100,mean=25.65,sd=5),
                                     rnorm(n=100,mean=44.92,sd=6),
                                     rnorm(n=100,mean=40.16,sd=7))
input_params$net_use_historically = c(rnorm(n=100,mean=0.1234202,sd=0.01),
                                      rnorm(n=100,mean=0.08074346,sd=0.01),
                                      rnorm(n=100,mean=0.2039024,sd=0.01),
                                      rnorm(n=100,mean=0.2169812,sd=0.01),
                                      rnorm(n=100,mean=0.1576931,sd=0.01),
                                      rnorm(n=100,mean=0.0942001,sd=0.01))

input_params$net_use_yr0 = c(rnorm(n=100,mean=0.3781217,sd=0.01),
                             rnorm(n=100,mean=0.15765426,sd=0.01),
                             rnorm(n=100,mean=0.3036978,sd=0.01),
                             rnorm(n=100,mean=0.3541520,sd=0.01),
                             rnorm(n=100,mean=0.2566661,sd=0.01),
                             rnorm(n=100,mean=0.27842454,sd=0.01))
input_params$net_use_yr1 = c(rnorm(n=100,mean=0.5425129,sd=0.01),
                             rnorm(n=100,mean=0.36852716,sd=0.01),
                             rnorm(n=100,mean=0.4981379,sd=0.01),
                             rnorm(n=100,mean=0.483115,sd=0.01),
                             rnorm(n=100,mean=0.4350480,sd=0.01),
                             rnorm(n=100,mean=0.43503708,sd=0.01))

names(data[,5:11])
input_params$seasonal_a0=rep(data$seasonal_a0,100)
input_params$seasonal_a1=rep(data$seasonal_a1,100)
input_params$seasonal_b1=rep(data$seasonal_b1,100)
input_params$seasonal_a2=rep(data$seasonal_a2,100)
input_params$seasonal_b2=rep(data$seasonal_b2,100)
input_params$seasonal_a3=rep(data$seasonal_a3,100)
input_params$seasonal_b3=rep(data$seasonal_b3,100)

write.csv(input_params,"C:/Users/esherrar/Documents/larval_source_management/Data resource/Modelling/input_parameters_for_sites.csv")


##########################################
##
## Ratio Gamb : funestus
site_species_ratios_before = numeric(6)
for(i in 1:6){
  site_species_ratios_before[i] =
    sum(c(Count_gambiae_females_per_month[7:12,1,i],Count_gambiae_females_per_month[1:6,2,i]))/
    (sum(c(Count_gambiae_females_per_month[7:12,1,i],Count_gambiae_females_per_month[1:6,2,i]))+
       sum(c(Count_funestus_females_per_month[7:12,1,i],Count_funestus_females_per_month[1:6,2,i])))
}
