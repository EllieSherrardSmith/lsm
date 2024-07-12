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

dt = data.frame(time = 1:36, ## months
                cnt = site_1_summary[[1]],
                cnt2 = site_1_summary[[1]]+1)
# M1 <- stats::nls(cnt ~ a/(1 + exp(-b * (time-c)) ), start=list(a=40,b=0.5,c=1))
# Run regression nls with formula SSlogis both from package stats 
fit <- nls(cnt ~ SSlogis(time, a, b, c), data = dt)
summary(fit)
p <- coef(fit)
plot(dt$time,dt$cnt,title(main = "Sigmoid Fitting"),xlab = "Months",
     ylab = 'Counts',pch = 20,ylim=c(0,90))
curve(SSlogis(x,p["a"],p["b"],p["c"]),lwd = 2, col = 'blue',add = TRUE) # Predict
# curve(SSlogis(x,mean(dt$cnt[1:18]),18,-1),lwd = 2, col = 'grey',add = TRUE) # Predict

## look at a logistic curve
d = min(site_1_summary[[1]])
cc = median(site_1_summary[[1]][1:18])
b = -5
e = 18
pred = cc + ((d - cc)/(1 + exp(b * (time - e))))
lines(pred ~ time)

# fit2 <- nls(ctimefit2 <- nls(cnt ~ SSfpl(time, a, b, xmid, scal), data = dt)
# p2 = coef(fit2)
# curve(SSfpl(x,p2["a"],p2["b"],p2["xmid"],p2["scal"]),lwd = 2, col = 'darkred',add = TRUE) # Predict

## LSM sites
dd1 = site_1_summary[[2]]
dd2 = site_3_summary[[2]]
dd3 = site_4_summary[[2]]

## LSM sites
dd1 = site_2_summary[[6]]
dd2 = site_5_summary[[6]]
dd3 = site_6_summary[[6]]

reduction_f = function(params){ 
  
  times = seq(1,36*30,length=36) ## days 

  cc <- params[1] ## prior average mosquito counts
  d  <- params[2] ## post LSM average mosquito counts
  b  <- params[3] ## malariasimulation parameter$lsm_rate_beta
  e  <- times[18] ## malariasimulation parameter$lsm_rate_beta
  
  ## function determined in publication
  pred1 = cc + ((d - cc)/(1 + exp(b * (times - e))))
  
  ## trial data from Kenya for gambiae mosquito
  dat1 = dd1 
  dat2 = dd2
  dat3 = dd3
 
  ## maximum likelihood
  output = sum((dat1-pred1)^2 + (dat2-pred1)^2 + (dat3-pred1)^2)
  output
}

## parameter estimates prior to learning optimal solution 
param_guess = c(mean(c(dd1[1:18],dd2[1:18],dd3[1:18])),
                mean(c(dd1[19:36],dd2[19:36],dd3[19:36])),
                -15) 

## optimise the solution for sensible parameter ranges
## lsm_new_eqm must be between -1 and 1
## lsm_rate_alpha must be positive
## lsm_rate_beta must be negative
satmod = optim(param_guess,lower = c(0,0,-200), upper = c(90,80,-0.001), reduction_f,method="L-BFGS-B")
satmod

## with LSM we go from 20 to 0, absolute reduction of 20 female gam
## from 27 to 0 for any gambiae, absolute reduction of 27
## from 7 to 0 for funestus, absolute reduction of 7
## from 10 to 0 for any funestus, absolute reduction of 10
## from 73 to 28 for culex, absolute reduction of 45 

## with non-LSM we go from 22 to 4, absolute reduction of 18
## from 26 to 5 for any gambiae, absolute reduction of 21
## from 8 to 1 for funestus, absolute reduction of 7
## from 10 to 1 for any funestus, absolute reduction of 9
## from 59 to 22 for culex, absolute reduction of 37
cc = 26.9724649 
d = 0.2332539
b = -15
times = seq(1,36*30,length=36)
pred1 = cc + ((d - cc)/(1 + exp(b * (times - times[18]))))
plot(pred1 ~ times)
