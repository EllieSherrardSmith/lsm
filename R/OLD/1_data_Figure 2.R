##############################################
##
## Figure 1 data figure

## Data used in the transmission modelling

## Prevalence before

## average net use each year

## assumed seasonality

## LSM USED - OR NOT


site = 1 #yes
site = 2 #no
site = 3 #yes
site = 4 #yes
site = 5 #no
site = 6 #no


## Mosquito density reduction


###################################
##
## Figure 1a 

## Demo for site 1
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


plot(sites1[1:43] ~ timeline[1:43],ylab = "",
     xlab="Date",yaxt="n",xlim=c(6,43),
     main="Musilongo: LSM",
     xaxt="n",pch="",type="l",lwd=2,col="darkgreen",ylim=c(0,1))
axis(1,at=c(1,7,13,19,25,31,37,43),labels=c("Jan 2004","Jun 2004","Jan 2005",
                                      "Jun 2005","Jan 2006",
                                      "Jun 2006","Jan 2007","Jun 2007"))
axis(2,las=2,at=seq(0,1,0.2),labels=c(0,20,40,60,80,100))

polygon(c(1:19,19:1),c(rep(0,19),rep(0.073,19)),border=NA,col=adegenet::transp("grey",0.3))
polygon(c(19:31,31:19),c(rep(0,13),rep(0.2097,13)),border=NA,col=adegenet::transp("grey",0.3))
polygon(c(31:43,43:31),c(rep(0,13),rep(0.4923,13)),border=NA,col=adegenet::transp("grey",0.3))


## DIDECODE 101471
data = read.csv("data/site_parameters.csv")

TIME = 1:(365*3.5)
seasonality = (data$seasonal_a0[which(data$DIDE_CODE == 101471)]+
                 -data$seasonal_a1[which(data$DIDE_CODE == 101471)]*cos(2*pi*TIME/365)+
                 data$seasonal_a2[which(data$DIDE_CODE == 101471)]*cos(2*2*pi*TIME/365)+
                 -data$seasonal_a3[which(data$DIDE_CODE == 101471)]*cos(3*2*pi*TIME/365)+
                 -data$seasonal_b1[which(data$DIDE_CODE == 101471)]*sin(2*pi*TIME/365)+
                 data$seasonal_b2[which(data$DIDE_CODE == 101471)]*sin(2*2*pi*TIME/365)+ 
                 -data$seasonal_b3[which(data$DIDE_CODE == 101471)]*sin(3*2*pi*TIME/365) )

seasonality[which(seasonality < 0)] = 0.001

theta_c =  mean(seasonality)
mal1 = seasonality/theta_c

#Normalized Data
normalized = (mal1-min(mal1))/(max(mal1)-min(mal1))

timeline2 = seq(6,48,length=1277)
lines(normalized~timeline2,col=adegenet::transp("blue",0.2),lwd=5)

points(sites1 ~ timeline,col="darkgreen",pch=15,cex=1.4)

# dat
dat = read.csv("data/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE)

summary(dat)

dat$Count_gambiae_females_per_month = dat$GA.TOT - dat$GAMBIAE.males 
dat$Count_funestus_females_per_month = dat$FU.TOT - dat$FUNESTUS.males
dat$Count_culex_females_per_month = dat$CU.TOT - dat$CULEX.males 

years = unique(dat$Year)

site = 1
is_this_larvicidal_site = 1
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

d1=c(Count_gambiae_females_per_month)[1:30] 
d2=c(Count_gambiae_total_per_month)[1:30]
d3=c(Count_funestus_females_per_month)[1:30] 
d4=c(Count_funestus_total_per_month)[1:30]
d5=c(Count_culex_females_per_month)[1:30] 
d6=c(Count_culex_total_per_month)[1:30]
TIME = 1:30

cols = c("purple","red")


polygon(c(3,5,5,3),c(0,0,max(d1),max(d1)),border=FALSE,col=adegenet::transp("blue",0.2))
polygon(c(15,17,17,15),c(0,0,max(d1),max(d1)),border=FALSE,col=adegenet::transp("blue",0.2))
polygon(c(27,29,29,27),c(0,0,max(d1),max(d1)),border=FALSE,col=adegenet::transp("blue",0.2))

polygon(c(9,10,10,9),c(0,0,max(d1),max(d1)),border=FALSE,col=adegenet::transp("blue",0.2))
polygon(c(21,22,22,21),c(0,0,max(d1),max(d1)),border=FALSE,col=adegenet::transp("blue",0.2))

abline(v=19, lty=2,col="grey")

# legend(22*30,max(d),legend = c("An. gambiae","An. funestus","Culex sp","Total","Females only"),
#        col=c(cols,"black","black"),lty=c(1,1,1,1,2))

par(new = T)
plot(c(d1,rep(NA,12)) ~ c(6:47), type="l", axes=F, xlab=NA, ylab=NA, col=cols[2])
par(new = T)
plot(c(d3,rep(NA,12)) ~ c(6:47), type="l", axes=F, xlab=NA, ylab=NA, col=cols[2],lty=2)

legend("topright",legend = c("Prevalence",
                             expression(italic('An. gambiae')),
                             expression(italic('An. funestus')),
                             "Net use (%)",
                             "Assumed seasonal trend"),
       col=c("darkgreen","red","red",adegenet::transp("grey",0.6),adegenet::transp("blue",0.2)),
       pch=15,lty=c(NA,1,2,NA,1),lwd=c(1,1,1,1,10))


plotting_fun_f = function(sites1,
                          print_main_site_and_trt,
                          dat,
                          site_for_lsm,
                          site_base_net_use,
                          site_yr0_net_use,
                          site_yr1_net_use){
  plot(sites1[1:43] ~ timeline[1:43],ylab = "",
       xlab="",yaxt="n",xlim=c(6,43),
       main=print_main_site_and_trt,
       xaxt="n",pch="",type="l",lwd=2,col="darkgreen",ylim=c(0,1))
  axis(1,at=c(1,7,13,19,25,31,37,43),labels=c("Jan 2004","Jun 2004","Jan 2005",
                                              "Jun 2005","Jan 2006",
                                              "Jun 2006","Jan 2007","Jun 2007"))
  axis(2,las=2,at=seq(0,1,0.2),labels=c(0,20,40,60,80,100))
  
  polygon(c(1:19,19:1),c(rep(0,19),rep(site_base_net_use,19)),border=NA,col=adegenet::transp("grey",0.3))
  polygon(c(19:31,31:19),c(rep(0,13),rep(site_yr0_net_use,13)),border=NA,col=adegenet::transp("grey",0.3))
  polygon(c(31:43,43:31),c(rep(0,13),rep(site_yr1_net_use,13)),border=NA,col=adegenet::transp("grey",0.3))
  
  
  
  TIME = 1:(365*3.5)
  # seasonality = (data$seasonal_a0[which(data$DIDE_CODE == 101471)]+
  #                  -data$seasonal_a1[which(data$DIDE_CODE == 101471)]*cos(2*pi*TIME/365)+
  #                  data$seasonal_a2[which(data$DIDE_CODE == 101471)]*cos(2*2*pi*TIME/365)+
  #                  -data$seasonal_a3[which(data$DIDE_CODE == 101471)]*cos(3*2*pi*TIME/365)+
  #                  -data$seasonal_b1[which(data$DIDE_CODE == 101471)]*sin(2*pi*TIME/365)+
  #                  data$seasonal_b2[which(data$DIDE_CODE == 101471)]*sin(2*2*pi*TIME/365)+ 
  #                  -data$seasonal_b3[which(data$DIDE_CODE == 101471)]*sin(3*2*pi*TIME/365) )
  # 
  # seasonality[which(seasonality < 0)] = 0.001
  # 
  # theta_c =  mean(seasonality)
  # mal1 = seasonality/theta_c
  # 
  # #Normalized Data
  # normalized = (mal1-min(mal1))/(max(mal1)-min(mal1))
  # 
  # timeline2 = seq(6,48,length=1277)
  # lines(normalized~timeline2,col=adegenet::transp("blue",0.2),lwd=5)
  # 
  points(sites1 ~ timeline,col="darkgreen",pch=15,cex=1.4)
  
  
  dat$Count_gambiae_females_per_month = dat$GA.TOT - dat$GAMBIAE.males 
  dat$Count_funestus_females_per_month = dat$FU.TOT - dat$FUNESTUS.males
  dat$Count_culex_females_per_month = dat$CU.TOT - dat$CULEX.males 
  
  years = unique(dat$Year)
  
  site = site_for_lsm
  is_this_larvicidal_site = 0
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
  
  d1=c(Count_gambiae_females_per_month)[1:30] 
  d2=c(Count_gambiae_total_per_month)[1:30]
  d3=c(Count_funestus_females_per_month)[1:30] 
  d4=c(Count_funestus_total_per_month)[1:30]
  d5=c(Count_culex_females_per_month)[1:30] 
  d6=c(Count_culex_total_per_month)[1:30]
  TIME = 1:30
  
  cols = c("purple","red")
  
  
  # polygon(c(3,5,5,3),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("blue",0.2))
  # polygon(c(15,17,17,15),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("blue",0.2))
  # polygon(c(27,29,29,27),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("blue",0.2))
  # 
  # polygon(c(9,10,10,9),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("blue",0.2))
  # polygon(c(21,22,22,21),c(0,0,max(d),max(d)),border=FALSE,col=adegenet::transp("blue",0.2))
  # 
  abline(v=19, lty=2,col="grey")
  
  # legend(22*30,max(d),legend = c("An. gambiae","An. funestus","Culex sp","Total","Females only"),
  #        col=c(cols,"black","black"),lty=c(1,1,1,1,2))
  
  par(new = T)
  plot(c(d1,rep(NA,12)) ~ c(6:47), 
       type="l", 
       ylim=c(0,100),
       axes=F, xlab=NA, 
       ylab=NA, col=cols[2])
  par(new = T)
  plot(c(d3,rep(NA,12)) ~ c(6:47), 
       type="l", axes=F, 
       xlab=NA, ylab=NA, 
       ylim=c(0,100), col=cols[2],lty=2)
  
}



##############################################
##
## Potential figure 1

par(mfrow=c(2,3))
par(mar=c(3,5,2,5))

plotting_fun_f(sites1 = sites1,
               print_main_site_and_trt = "Musilongo: Larviciding",
               dat = read.csv("C:/Users/esherrar/Documents/larval_source_management/Data resource/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE),
               site_for_lsm = 1,
               site_base_net_use = 0.073,
               site_yr0_net_use = 0.2097,
               site_yr1_net_use = 0.4923)
axis(4,las=2,at=seq(0,100,20))
mtext(side = 2, line = 2.5,"Prevalence in 0.5 to 10-years (green points, %)",cex = 0.8)
mtext(side = 2, line = 4,"Mosquito net use (grey bars, %)",cex = 0.8)


legend("topright",legend = c("Prevalence",
                             expression(italic('An. gambiae')),
                             expression(italic('An. funestus')),
                             "Net use (%)"),
       col=c("darkgreen","red","red",adegenet::transp("grey",0.6)),#,adegenet::transp("blue",0.2)),
       pch=c(15,15,15,15,NA),lty=c(NA,1,2,NA,1),lwd=c(1,1,1,1,5),cex = 1.1)


plotting_fun_f(sites1 = sites3,
               print_main_site_and_trt = "Kezege: Larviciding",
               dat = read.csv("C:/Users/esherrar/Documents/larval_source_management/Data resource/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE),
               site_for_lsm = 3,
               site_base_net_use = 0.0615,
               site_yr0_net_use = 0.2891,
               site_yr1_net_use = 0.4299)
axis(4,las=2,at=seq(0,100,20))

plotting_fun_f(sites1 = sites4,
               print_main_site_and_trt = "Wamondo: Larviciding",
               dat = read.csv("C:/Users/esherrar/Documents/larval_source_management/Data resource/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE),
               site_for_lsm = 4,
               site_base_net_use = 0.121,
               site_yr0_net_use = 0.2980,
               site_yr1_net_use = 0.4291)
axis(4,las=2,at=seq(0,100,20))
mtext(side = 4, line = 2.5,"Adult mosquito density (count per month)",cex = 0.8)


plotting_fun_f(sites1 = sites2,
               print_main_site_and_trt = "Kimingini: No Larviciding",
               dat = read.csv("C:/Users/esherrar/Documents/larval_source_management/Data resource/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE),
               site_for_lsm = 2,
               site_base_net_use = 0.038,
               site_yr0_net_use = 0.1494,
               site_yr1_net_use = 0.2655)
axis(4,las=2,at=seq(0,100,20))
mtext(side = 2, line = 2.5,"Prevalence in 0.5 to 10-years (green points, %)",cex = 0.8)
mtext(side = 2, line = 4,"Mosquito net use (grey bars, %)",cex = 0.8)

plotting_fun_f(sites1 = sites5,
               print_main_site_and_trt = "Emutete: No Larviciding",
               dat = read.csv("C:/Users/esherrar/Documents/larval_source_management/Data resource/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE),
               site_for_lsm = 5,
               site_base_net_use = 0.087,
               site_yr0_net_use = 0.2019,
               site_yr1_net_use = 0.3634)
axis(4,las=2,at=seq(0,100,20))



plotting_fun_f(sites1 = sites6,
               print_main_site_and_trt = "Wakikuyu: No Larviciding",
               dat = read.csv("C:/Users/esherrar/Documents/larval_source_management/Data resource/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE),
               site_for_lsm = 6,
               site_base_net_use = 0.051,
               site_yr0_net_use = 0.1328,
               site_yr1_net_use = 0.3657)
axis(4,las=2,at=seq(0,100,20))
mtext(side = 4, line = 2.5,"Adult mosquito density (count per month)",cex = 0.8)


########################################################
##
## Next to analysis of reductions

