###################################
##
## glm analysis
##
##################################

dat = read.csv("data/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE)

dat$Count_gambiae_females_per_month = dat$GA.TOT - dat$GAMBIAE.males 
dat$Count_funestus_females_per_month = dat$FU.TOT - dat$FUNESTUS.males
dat$Count_culex_females_per_month = dat$CU.TOT - dat$CULEX.males 

years = unique(dat$Year)

#
counts_dat = function(site, dat){
  
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
  
  
  
  return(data.frame(d1=c(Count_gambiae_females_per_month)[1:30], 
                    d2=c(Count_gambiae_total_per_month)[1:30],
                    d3=c(Count_funestus_females_per_month)[1:30], 
                    d4=c(Count_funestus_total_per_month)[1:30],
                    d5=c(Count_culex_females_per_month)[1:30],
                    d6=c(Count_culex_total_per_month)[1:30],
                    TIME = 1:30))
  
  
  
}

#

s1 = counts_dat(site = 1,dat = dat)
s2 = counts_dat(site = 2,dat = dat)
s3 = counts_dat(site = 3,dat = dat)
s4 = counts_dat(site = 4,dat = dat)
s5 = counts_dat(site = 5,dat = dat)
s6 = counts_dat(site = 6,dat = dat)

## create data
## start with vectors, then all Anopheles
dstat = expand.grid(count_vec = c(s1$d1,s1$d3,s2$d1,s2$d3,s3$d1,s3$d3,s4$d1,s4$d3,s5$d1,s5$d3,s6$d1,s6$d3))

dstat$count_mos = c(s1$d2,s1$d4,s2$d2,s2$d4,s3$d2,s3$d4,s4$d2,s4$d4,s5$d2,s5$d4,s6$d2,s6$d4)
dstat$site = rep(1:6, each = 60)
dstat$arm = rep(c(1,0,1,1,0,0),each = 60)          ## LSM or not
dstat$species = rep(rep(c("gam","fun"),each=30),6) ##
dstat$month = rep(1:30,12)
dstat$binary = rep(rep(c(0,1),c(18,12)),12)        ## before and after the intervention
dstat$binary6to12 = rep(rep(c(0,99,1),c(18,6,6)),12)        ## before and after the intervention

head(dstat)


library("rstanarm")
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan1 <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstat,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  

# plot(stan1)

100*(1-exp(stan1$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stan1$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stan1$coefficients[4])) ##difference in difference

100*(1-exp(stan1$coefficients[2]-stan1$ses[2])) ##difference in difference
100*(1-exp(stan1$coefficients[2]+stan1$ses[2])) ##difference in difference
100*(1-exp(stan1$coefficients[3]-stan1$ses[3])) ##difference in difference
100*(1-exp(stan1$coefficients[3]+stan1$ses[3])) ##difference in difference


100*(1-exp(stan1$coefficients[4]-stan1$ses[4])) ##difference in difference
100*(1-exp(stan1$coefficients[4]+stan1$ses[4])) ##difference in difference

par(mfrow = c(2,2))
par(mar = c(4,4,1,1))
boxplot(dstat$count_mos ~ dstat$arm + dstat$binary,
        ylab = "Anopheles mosquitoes counts",
        xlab = "",xaxt="n",
        yaxt = "n") ## 0 no LSM, 1 some, 0 before, 1 after
axis(2, las = 2, seq(0,100,20))
axis(1, las = 1, at=c(1,2,3,4), labels = c("No LSM", "LSM","No LSM", "LSM"))
mtext("Before                                   After", side = 1, line = 2.5)
points(dstat$count_mos[dstat$arm == 0 & dstat$binary == 0] ~ sample(size = length(dstat$count_mos[dstat$arm == 0 & dstat$binary == 0]),
                                                                    x = rnorm(n = 20,mean = 1, sd = 0.15),
                                                                    replace = TRUE),
       pch=19,col = "grey40")
points(dstat$count_mos[dstat$arm == 1 & dstat$binary == 0] ~ sample(size = length(dstat$count_mos[dstat$arm == 0 & dstat$binary == 0]),
                                                                    x = rnorm(n = 20,mean = 2, sd = 0.15),
                                                                    replace = TRUE),
       pch=19,col = "darkgreen")

points(dstat$count_mos[dstat$arm == 0 & dstat$binary == 1] ~ sample(size = length(dstat$count_mos[dstat$arm == 0 & dstat$binary == 1]),
                                                                    x = rnorm(n = 20,mean = 3, sd = 0.15),
                                                                    replace = TRUE),
       pch=19,col = "grey40")
points(dstat$count_mos[dstat$arm == 1 & dstat$binary == 1] ~ sample(size = length(dstat$count_mos[dstat$arm == 0 & dstat$binary == 1]),
                                                                    x = rnorm(n = 20,mean = 4, sd = 0.15),
                                                                    replace = TRUE),
       pch=19,col = "darkgreen")




mid_bound1 = c(exp(posterior_interval(stan1, prob = 0.5)[1,1]),
               exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[3,1]),
               exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[2,1]),
               exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[4,1]))

lower_bound1 = c(exp(posterior_interval(stan1, prob = 0.95)[1,1]),
                 exp(posterior_interval(stan1, prob = 0.95)[1,1] + posterior_interval(stan1, prob = 0.95)[3,1]),
                 exp(posterior_interval(stan1, prob = 0.95)[1,1] + posterior_interval(stan1, prob = 0.95)[2,1]),
                 exp(posterior_interval(stan1, prob = 0.95)[1,1] + posterior_interval(stan1, prob = 0.95)[4,1]))

upper_bound1 = c(exp(posterior_interval(stan1, prob = 0.95)[1,2]),
                 exp(posterior_interval(stan1, prob = 0.95)[1,2] + posterior_interval(stan1, prob = 0.95)[3,2]),
                 exp(posterior_interval(stan1, prob = 0.95)[1,2] + posterior_interval(stan1, prob = 0.95)[2,2]),
                 exp(posterior_interval(stan1, prob = 0.95)[1,2] + posterior_interval(stan1, prob = 0.95)[4,2]))

segs1 = c(0.96,1.96,2.96,3.96)
for(i in 1:4){
  
  segments(x0 = segs1[i],
           x1 = segs1[i],
           y0 = lower_bound1[i],
           y1 = upper_bound1[i],
           lwd=4,col="red")
  
}

points(mid_bound1 ~ segs1,pch=19,cex = 1.2, col = "red")


# pp_check(stan1, plotfun = "hist", nreps = 5)
# y_rep = posterior_predict(stan1)
# check_predictions_f = function(group){
#   y_print = log(dstat$count_mos[which(dstat$binary == group)]+1)
#   y_sample = log(y_rep[sample(1000,100),which(dstat$binary == group)]+1)
#   bayesplot::ppc_dens_overlay(y_print,y_sample)
#   
# }
# check_predictions_f(group = 0)
# check_predictions_f(group = 1)
# 

### And here is th 1 to 6 months difference
dstatEARLY = subset(dstat,dstat$binary6to12 == 99 | dstat$binary6to12 == 0)

stanearly <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstatEARLY,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  

# plot(stanearly)

100*(1-exp(stanearly$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stanearly$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stanearly$coefficients[4])) ##difference in difference

100*(1-exp(stanearly$coefficients[2]-stanearly$ses[2])) ##difference in difference
100*(1-exp(stanearly$coefficients[2]+stanearly$ses[2])) ##difference in difference
100*(1-exp(stanearly$coefficients[3]-stanearly$ses[3])) ##difference in difference
100*(1-exp(stanearly$coefficients[3]+stanearly$ses[3])) ##difference in difference


100*(1-exp(stanearly$coefficients[4]-stanearly$ses[4])) ##difference in difference
100*(1-exp(stanearly$coefficients[4]+stanearly$ses[4])) ##difference in difference


### And now what is the difference in differences for the 7-12 months change
dstatLATE = subset(dstat,dstat$binary6to12 != 99)

stanlate <- rstanarm::stan_glm.nb(
  count_mos ~ binary6to12 * arm, ## species not important so removed
  data = dstatLATE,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  

# plot(stanlate)

100*(1-exp(stanlate$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stanlate$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stanlate$coefficients[4])) ##difference in difference

100*(1-exp(stanlate$coefficients[2]-stanlate$ses[2])) ##difference in difference
100*(1-exp(stanlate$coefficients[2]+stanlate$ses[2])) ##difference in difference
100*(1-exp(stanlate$coefficients[3]-stanlate$ses[3])) ##difference in difference
100*(1-exp(stanlate$coefficients[3]+stanlate$ses[3])) ##difference in difference


100*(1-exp(stanlate$coefficients[4]-stanlate$ses[4])) ##difference in difference
100*(1-exp(stanlate$coefficients[4]+stanlate$ses[4])) ##difference in difference

boxplot(log(dstat$count_mos+1) ~ dstat$arm + dstat$binary,
        ylab = "Anopheles mosquitoes counts (log scale)",
        xlab = "",xaxt="n",
        yaxt = "n") ## 0 no LSM, 1 some, 0 before, 1 after
axis(2, las = 2, at = log(seq(0,100,20)+1), labels=seq(0,100,20))
axis(1, las = 1, at=c(1,2,3,4), labels = c("No LSM", "LSM","No LSM", "LSM"))
mtext("Before                                   After", side = 1, line = 2.5)

points(log(dstat$count_mos[dstat$arm == 0 & dstat$binary == 0]+1) ~ sample(size = length(dstat$count_mos[dstat$arm == 0 & dstat$binary == 0]),
                                                                    x = rnorm(n = 20,mean = 1, sd = 0.15),
                                                                    replace = TRUE),
       pch=19,col = "grey40")

points(log(dstat$count_mos[dstat$arm == 1 & dstat$binary == 0]+1) ~ sample(size = length(dstat$count_mos[dstat$arm == 0 & dstat$binary == 0]),
                                                                    x = rnorm(n = 20,mean = 2, sd = 0.15),
                                                                    replace = TRUE),
       pch=19,col = "darkgreen")

points(log(dstat$count_mos[dstat$arm == 0 & dstat$binary6to12 == 99]+1) ~ sample(size = length(dstat$count_mos[dstat$arm == 0 & dstat$binary6to12 == 99]),
                                                                          x = rnorm(n = 20,mean = 2.8, sd = 0.05),
                                                                          replace = TRUE),
       pch=19,col = adegenet::transp("darkblue",0.8))

points(log(dstat$count_mos[dstat$arm == 0 & dstat$binary6to12 == 1]+1) ~ sample(size = length(dstat$count_mos[dstat$arm == 0 & dstat$binary6to12 == 1]),
                                                                         x = rnorm(n = 20,mean = 3.2, sd = 0.05),
                                                                         replace = TRUE),
       pch=19,col = adegenet::transp("blue",0.4))


points(log(dstat$count_mos[dstat$arm == 1 & dstat$binary6to12 == 99]+1) ~ sample(size = length(dstat$count_mos[dstat$arm == 0 & dstat$binary6to12 == 99]),
                                                                          x = rnorm(n = 20,mean = 3.8, sd = 0.05),
                                                                          replace = TRUE),
       pch=19,col = adegenet::transp("aquamarine2",0.4))

points(log(dstat$count_mos[dstat$arm == 1 & dstat$binary6to12 == 1]+1) ~ sample(size = length(dstat$count_mos[dstat$arm == 0 & dstat$binary6to12 == 1]),
                                                                         x = rnorm(n = 20,mean = 4.2, sd = 0.05),
                                                                         replace = TRUE),
       pch=19,col = adegenet::transp("darkgreen",0.4))



mid_bound1late = c(exp(posterior_interval(stanlate, prob = 0.5)[1,1]),
                   exp(posterior_interval(stanlate, prob = 0.5)[1,1] + posterior_interval(stanlate, prob = 0.5)[3,1]),
                   exp(posterior_interval(stanlate, prob = 0.5)[1,1] + posterior_interval(stanlate, prob = 0.5)[2,1]),
                   exp(posterior_interval(stanlate, prob = 0.5)[1,1] + posterior_interval(stanlate, prob = 0.5)[4,1]))

lower_bound1late = c(exp(posterior_interval(stanlate, prob = 0.95)[1,1]),
                     exp(posterior_interval(stanlate, prob = 0.95)[1,1] + posterior_interval(stanlate, prob = 0.95)[3,1]),
                     exp(posterior_interval(stanlate, prob = 0.95)[1,1] + posterior_interval(stanlate, prob = 0.95)[2,1]),
                     exp(posterior_interval(stanlate, prob = 0.95)[1,1] + posterior_interval(stanlate, prob = 0.95)[4,1]))

upper_bound1late = c(exp(posterior_interval(stanlate, prob = 0.95)[1,2]),
                     exp(posterior_interval(stanlate, prob = 0.95)[1,2] + posterior_interval(stanlate, prob = 0.95)[3,2]),
                     exp(posterior_interval(stanlate, prob = 0.95)[1,2] + posterior_interval(stanlate, prob = 0.95)[2,2]),
                     exp(posterior_interval(stanlate, prob = 0.95)[1,2] + posterior_interval(stanlate, prob = 0.95)[4,2]))

segs1late = c(1.1,2.1,3.1,4.1)
for(i in 1:4){
  
  segments(x0 = segs1late[i],
           x1 = segs1late[i],
           y0 = log(lower_bound1late[i]+1),
           y1 = log(upper_bound1late[i]+1),
           lwd=4,col="blue")
  
}

points(log(mid_bound1late+1) ~ segs1late,pch=19,cex = 1.2, col = "blue")



# segs1 = c(0.96,1.96,2.96,3.96)
# for(i in 1:4){
#   
#   segments(x0 = segs1[i],
#            x1 = segs1[i],
#            y0 = lower_bound1[i],
#            y1 = upper_bound1[i],
#            lwd=4,col="red")
#   
# }


mid_bound1early = c(exp(posterior_interval(stanearly, prob = 0.5)[1,1]),
                    exp(posterior_interval(stanearly, prob = 0.5)[1,1] + posterior_interval(stanearly, prob = 0.5)[3,1]),
                    exp(posterior_interval(stanearly, prob = 0.5)[1,1] + posterior_interval(stanearly, prob = 0.5)[2,1]),
                    exp(posterior_interval(stanearly, prob = 0.5)[1,1] + posterior_interval(stanearly, prob = 0.5)[4,1]))

lower_bound1early = c(exp(posterior_interval(stanearly, prob = 0.95)[1,1]),
                      exp(posterior_interval(stanearly, prob = 0.95)[1,1] + posterior_interval(stanearly, prob = 0.95)[3,1]),
                      exp(posterior_interval(stanearly, prob = 0.95)[1,1] + posterior_interval(stanearly, prob = 0.95)[2,1]),
                      exp(posterior_interval(stanearly, prob = 0.95)[1,1] + posterior_interval(stanearly, prob = 0.95)[4,1]))

upper_bound1early = c(exp(posterior_interval(stanearly, prob = 0.95)[1,2]),
                      exp(posterior_interval(stanearly, prob = 0.95)[1,2] + posterior_interval(stanearly, prob = 0.95)[3,2]),
                      exp(posterior_interval(stanearly, prob = 0.95)[1,2] + posterior_interval(stanearly, prob = 0.95)[2,2]),
                      exp(posterior_interval(stanearly, prob = 0.95)[1,2] + posterior_interval(stanearly, prob = 0.95)[4,2]))

segs1early = c(0.88,1.88,2.88,3.88)
for(i in 1:4){
  
  segments(x0 = segs1early[i],
           x1 = segs1early[i],
           y0 = log(lower_bound1early[i]+1),
           y1 = log(upper_bound1early[i]+1),
           lwd=4,col="black")
  
}

points(log(mid_bound1late+1) ~ segs1early,pch=19,cex = 1.2, col = "black")

legend("topright",legend = c("0-6 months",
                             "7-12 months"),
       col = c("black","blue"),
       pch=19,border=NA,bty="n")

####################################################
##
## And village specific estimates
##
####################################################

dat = read.csv("data/Fillinger_Kenya highland LSM_ADULT MOSQUITO_SPRAY CATCHES_RAW DATA.csv",header=TRUE)

dat$Count_gambiae_females_per_month = dat$GA.TOT - dat$GAMBIAE.males 
dat$Count_funestus_females_per_month = dat$FU.TOT - dat$FUNESTUS.males
dat$Count_culex_females_per_month = dat$CU.TOT - dat$CULEX.males 

years = unique(dat$Year)

s1 = counts_dat(site = 1,dat = dat)
s2 = counts_dat(site = 2,dat = dat)
s3 = counts_dat(site = 3,dat = dat)
s4 = counts_dat(site = 4,dat = dat)
s5 = counts_dat(site = 5,dat = dat)
s6 = counts_dat(site = 6,dat = dat)

## create data
## start with vectors, then all Anopheles
dstat = expand.grid(count_vec = c(s1$d1,s1$d3,s2$d1,s2$d3,s3$d1,s3$d3,s4$d1,s4$d3,s5$d1,s5$d3,s6$d1,s6$d3))

dstat$count_mos = c(s1$d2,s1$d4,s2$d2,s2$d4,s3$d2,s3$d4,s4$d2,s4$d4,s5$d2,s5$d4,s6$d2,s6$d4)
dstat$site = rep(1:6, each = 60)
dstat$arm = rep(c(1,0,1,1,0,0),each = 60)          ## LSM or not
dstat$species = rep(rep(c("gam","fun"),each=30),6) ##
dstat$month = rep(1:30,12)
dstat$binary = rep(rep(c(0,1),c(18,12)),12)        ## before and after the intervention
dstat$binary6to12 = rep(rep(c(0,99,1),c(18,6,6)),12)        ## before and after the intervention

head(dstat)

gam = subset(dstat, dstat$species == "gam")
boxplot(gam$count_mos ~ gam$arm + gam$binary,
        ylab = "gambiae counts",
        xlab = "",xaxt = "n",ylim = c(0,160), 
        yaxt = "n") ## 0 no LSM, 1 some, 0 before, 1 after
axis(2, las = 2, at = seq(0,100,20))
axis(1, at = c(1,2,3,4), labels = c("No LSM", "LSM",  "No LSM",  "LSM"))
mtext("Before                                           After" ,side = 1, line = 2.5)


## Village 1
## subset for site 1
dstat1 = subset(dstat, dstat$site != 3 & dstat$site != 4 & dstat$species == "gam")
stan1 <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstat1,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  
100*(1-exp(stan1$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stan1$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stan1$coefficients[4])) ##difference in difference

100*(1-exp(stan1$coefficients[2]-stan1$ses[2])) ##difference in difference
100*(1-exp(stan1$coefficients[2]+stan1$ses[2])) ##difference in difference
100*(1-exp(stan1$coefficients[4]-stan1$ses[4])) ##difference in difference
100*(1-exp(stan1$coefficients[4]+stan1$ses[4])) ##difference in difference

mid_b1 = c(exp(posterior_interval(stan1, prob = 0.5)[1,1]),
           exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[3,1]),
           exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[2,1]),
           exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[4,1]))

lower_b1 = c(exp(posterior_interval(stan1, prob = 0.95)[1,1]),
             exp(posterior_interval(stan1, prob = 0.95)[1,1] + posterior_interval(stan1, prob = 0.95)[3,1]),
             exp(posterior_interval(stan1, prob = 0.95)[1,1] + posterior_interval(stan1, prob = 0.95)[2,1]),
             exp(posterior_interval(stan1, prob = 0.95)[1,1] + posterior_interval(stan1, prob = 0.95)[4,1]))

upper_b1 = c(exp(posterior_interval(stan1, prob = 0.95)[1,2]),
             exp(posterior_interval(stan1, prob = 0.95)[1,2] + posterior_interval(stan1, prob = 0.95)[3,2]),
             exp(posterior_interval(stan1, prob = 0.95)[1,2] + posterior_interval(stan1, prob = 0.95)[2,2]),
             exp(posterior_interval(stan1, prob = 0.95)[1,2] + posterior_interval(stan1, prob = 0.95)[4,2]))

segsv1 = c(0.7,1.7,2.7,3.7)
for(i in c(1:4)){
  
  segments(x0 = segsv1[i],
           x1 = segsv1[i],
           y0 = lower_b1[i],
           y1 = upper_b1[i],
           lwd=4,col="orange")
  
}

points(mid_b1 ~ segsv1,pch=19,cex = 1.2, col = "orange")



points(dstat1$count_mos[dstat1$arm == 0 & dstat1$binary == 0] ~ sample(size = length(dstat1$count_mos[dstat1$arm == 0 & dstat1$binary == 0]),
                                                              x = rnorm(n = 20,mean = 0.7, sd = 0.04),
                                                              replace = TRUE),
       pch=19,col = "grey40")

points(dstat1$count_mos[dstat1$arm == 1 & dstat1$binary == 0] ~ sample(size = length(dstat1$count_mos[dstat1$arm == 1 & dstat1$binary == 0]),
                                                                              x = rnorm(n = 20,mean = 1.7, sd = 0.04),
                                                                              replace = TRUE),
       pch=19,col = "darkgreen")

points(dstat1$count_mos[dstat1$arm == 0 & dstat1$binary == 1] ~ sample(size = length(dstat1$count_mos[dstat1$arm == 0 & dstat1$binary == 1]),
                                                                              x = rnorm(n = 20,mean = 2.7, sd = 0.04),
                                                                              replace = TRUE),
       pch=19,col = adegenet::transp("darkblue",0.8))


points(dstat1$count_mos[dstat1$arm == 1 & dstat1$binary == 1] ~ sample(size = length(dstat1$count_mos[dstat1$arm == 1 & dstat1$binary == 1]),
                                                                              x = rnorm(n = 20,mean = 3.7, sd = 0.04),
                                                                              replace = TRUE),
       pch=19,col = adegenet::transp("aquamarine2",0.4))


## Village 3
## subset for site 3
dstat3 = subset(dstat, dstat$site != 1 & dstat$site != 4 & dstat$species == "gam")
stan3 <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstat3,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  
100*(1-exp(stan3$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stan3$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stan3$coefficients[4])) ##difference in difference

100*(1-exp(stan3$coefficients[2]-stan3$ses[2])) ##difference in difference
100*(1-exp(stan3$coefficients[2]+stan3$ses[2])) ##difference in difference
100*(1-exp(stan3$coefficients[4]-stan3$ses[4])) ##difference in difference
100*(1-exp(stan3$coefficients[4]+stan3$ses[4])) ##difference in difference

mid_b1 = c(exp(posterior_interval(stan3, prob = 0.5)[1,1]),
           exp(posterior_interval(stan3, prob = 0.5)[1,1] + posterior_interval(stan3, prob = 0.5)[3,1]),
           exp(posterior_interval(stan3, prob = 0.5)[1,1] + posterior_interval(stan3, prob = 0.5)[2,1]),
           exp(posterior_interval(stan3, prob = 0.5)[1,1] + posterior_interval(stan3, prob = 0.5)[4,1]))

lower_b1 = c(exp(posterior_interval(stan3, prob = 0.95)[1,1]),
             exp(posterior_interval(stan3, prob = 0.95)[1,1] + posterior_interval(stan3, prob = 0.95)[3,1]),
             exp(posterior_interval(stan3, prob = 0.95)[1,1] + posterior_interval(stan3, prob = 0.95)[2,1]),
             exp(posterior_interval(stan3, prob = 0.95)[1,1] + posterior_interval(stan3, prob = 0.95)[4,1]))

upper_b1 = c(exp(posterior_interval(stan3, prob = 0.95)[1,2]),
             exp(posterior_interval(stan3, prob = 0.95)[1,2] + posterior_interval(stan3, prob = 0.95)[3,2]),
             exp(posterior_interval(stan3, prob = 0.95)[1,2] + posterior_interval(stan3, prob = 0.95)[2,2]),
             exp(posterior_interval(stan3, prob = 0.95)[1,2] + posterior_interval(stan3, prob = 0.95)[4,2]))

segsv1 = c(1,2,3,4)
for(i in 1:4){
  
  segments(x0 = segsv1[i],
           x1 = segsv1[i],
           y0 = lower_b1[i],
           y1 = upper_b1[i],
           lwd=4,col="red")
  
}

points(mid_b1 ~ segsv1,pch=19,cex = 1.2, col = "red")

points(dstat3$count_mos[dstat3$arm == 0 & dstat3$binary == 0] ~ sample(size = length(dstat3$count_mos[dstat3$arm == 0 & dstat3$binary == 0]),
                                                                       x = rnorm(n = 20,mean = 1, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("grey40",0.6))

points(dstat3$count_mos[dstat3$arm == 1 & dstat3$binary == 0] ~ sample(size = length(dstat3$count_mos[dstat3$arm == 1 & dstat3$binary == 0]),
                                                                       x = rnorm(n = 20,mean = 2, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("darkgreen",0.6))

points(dstat3$count_mos[dstat3$arm == 0 & dstat3$binary == 1] ~ sample(size = length(dstat3$count_mos[dstat3$arm == 0 & dstat3$binary == 1]),
                                                                       x = rnorm(n = 20,mean = 3, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("darkblue",0.5))


points(dstat3$count_mos[dstat3$arm == 1 & dstat3$binary == 1] ~ sample(size = length(dstat3$count_mos[dstat3$arm == 1 & dstat1$binary == 1]),
                                                                       x = rnorm(n = 20,mean = 4, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("aquamarine3",0.4))


## Village 4
## subset for site 4
dstat4 = subset(dstat, dstat$site != 1 & dstat$site != 3 & dstat$species == "gam")
stan4 <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstat4,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  
100*(1-exp(stan4$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stan4$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stan4$coefficients[4])) ##difference in difference

100*(1-exp(stan4$coefficients[2]-stan4$ses[2])) ##difference in difference
100*(1-exp(stan4$coefficients[2]+stan4$ses[2])) ##difference in difference
100*(1-exp(stan4$coefficients[4]-stan4$ses[4])) ##difference in difference
100*(1-exp(stan4$coefficients[4]+stan4$ses[4])) ##difference in difference

mid_b1 = c(exp(posterior_interval(stan4, prob = 0.5)[1,1]),
           exp(posterior_interval(stan4, prob = 0.5)[1,1] + posterior_interval(stan4, prob = 0.5)[3,1]),
           exp(posterior_interval(stan4, prob = 0.5)[1,1] + posterior_interval(stan4, prob = 0.5)[2,1]),
           exp(posterior_interval(stan4, prob = 0.5)[1,1] + posterior_interval(stan4, prob = 0.5)[4,1]))

lower_b1 = c(exp(posterior_interval(stan4, prob = 0.95)[1,1]),
             exp(posterior_interval(stan4, prob = 0.95)[1,1] + posterior_interval(stan4, prob = 0.95)[3,1]),
             exp(posterior_interval(stan4, prob = 0.95)[1,1] + posterior_interval(stan4, prob = 0.95)[2,1]),
             exp(posterior_interval(stan4, prob = 0.95)[1,1] + posterior_interval(stan4, prob = 0.95)[4,1]))

upper_b1 = c(exp(posterior_interval(stan4, prob = 0.95)[1,2]),
             exp(posterior_interval(stan4, prob = 0.95)[1,2] + posterior_interval(stan4, prob = 0.95)[3,2]),
             exp(posterior_interval(stan4, prob = 0.95)[1,2] + posterior_interval(stan4, prob = 0.95)[2,2]),
             exp(posterior_interval(stan4, prob = 0.95)[1,2] + posterior_interval(stan4, prob = 0.95)[4,2]))

segsv1 = c(1.3,2.3,3.3,4.3)
for(i in 1:4){
  
  segments(x0 = segsv1[i],
           x1 = segsv1[i],
           y0 = lower_b1[i],
           y1 = upper_b1[i],
           lwd=4,col="darkred")
  
}

points(mid_b1 ~ segsv1,pch=19,cex = 1.2, col = "darkred")

points(dstat4$count_mos[dstat4$arm == 0 & dstat4$binary == 0] ~ sample(size = length(dstat4$count_mos[dstat4$arm == 0 & dstat4$binary == 0]),
                                                                       x = rnorm(n = 20,mean = 1.3, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("grey40",0.3))

points(dstat4$count_mos[dstat4$arm == 1 & dstat4$binary == 0] ~ sample(size = length(dstat4$count_mos[dstat4$arm == 1 & dstat4$binary == 0]),
                                                                       x = rnorm(n = 20,mean = 2.3, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("darkgreen",0.3))

points(dstat4$count_mos[dstat4$arm == 0 & dstat4$binary == 1] ~ sample(size = length(dstat4$count_mos[dstat4$arm == 0 & dstat4$binary == 1]),
                                                                       x = rnorm(n = 20,mean = 3.3, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("darkblue",0.2))


points(dstat4$count_mos[dstat4$arm == 1 & dstat4$binary == 1] ~ sample(size = length(dstat4$count_mos[dstat4$arm == 1 & dstat1$binary == 1]),
                                                                       x = rnorm(n = 20,mean = 4.3, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("aquamarine3",0.8))


legend("topright",legend = c("Musilonge",
                             "Kezege",
                             "Wamondo"),
       col = c("orange","red","darkred"),
       pch=19,border = NA, bty = "n")



fun = subset(dstat, dstat$species == "fun")
boxplot(fun$count_mos ~ fun$arm + fun$binary,
        ylab = "funestus counts",
        xlab = "",xaxt = "n",ylim = c(0,50), 
        yaxt = "n") ## 0 no LSM, 1 some, 0 before, 1 after
axis(2, las = 2, at = seq(0,100,20))
axis(1, at = c(1,2,3,4), labels = c("No LSM", "LSM",  "No LSM",  "LSM"))
mtext("Before                                           After" ,side = 1, line = 2.5)

## Village 1
## subset for site 1
dstat1 = subset(dstat, dstat$site != 3 & dstat$site != 4 & dstat$species == "fun")
stan1 <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstat1,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  
100*(1-exp(stan1$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stan1$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stan1$coefficients[4])) ##difference in difference

100*(1-exp(stan1$coefficients[2]-stan1$ses[2])) ##difference in difference
100*(1-exp(stan1$coefficients[2]+stan1$ses[2])) ##difference in difference
100*(1-exp(stan1$coefficients[4]-stan1$ses[4])) ##difference in difference
100*(1-exp(stan1$coefficients[4]+stan1$ses[4])) ##difference in difference

mid_b1 = c(exp(posterior_interval(stan1, prob = 0.5)[1,1]),
           exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[3,1]),
           exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[2,1]),
           exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[4,1]))

lower_b1 = c(exp(posterior_interval(stan1, prob = 0.95)[1,1]),
             exp(posterior_interval(stan1, prob = 0.95)[1,1] + posterior_interval(stan1, prob = 0.95)[3,1]),
             exp(posterior_interval(stan1, prob = 0.95)[1,1] + posterior_interval(stan1, prob = 0.95)[2,1]),
             exp(posterior_interval(stan1, prob = 0.95)[1,1] + posterior_interval(stan1, prob = 0.95)[4,1]))

upper_b1 = c(exp(posterior_interval(stan1, prob = 0.95)[1,2]),
             exp(posterior_interval(stan1, prob = 0.95)[1,2] + posterior_interval(stan1, prob = 0.95)[3,2]),
             exp(posterior_interval(stan1, prob = 0.95)[1,2] + posterior_interval(stan1, prob = 0.95)[2,2]),
             exp(posterior_interval(stan1, prob = 0.95)[1,2] + posterior_interval(stan1, prob = 0.95)[4,2]))

segsv1 = c(0.7,1.7,2.7,3.7)
for(i in 1:4){
  
  segments(x0 = segsv1[i],
           x1 = segsv1[i],
           y0 = lower_b1[i],
           y1 = upper_b1[i],
           lwd=4,col="orange")
  
}

points(mid_b1 ~ segsv1,pch=19,cex = 1.2, col = "orange")


points(dstat1$count_mos[dstat1$arm == 0 & dstat1$binary == 0] ~ sample(size = length(dstat1$count_mos[dstat1$arm == 0 & dstat1$binary == 0]),
                                                                       x = rnorm(n = 20,mean = 0.7, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = "grey40")

points(dstat1$count_mos[dstat1$arm == 1 & dstat1$binary == 0] ~ sample(size = length(dstat1$count_mos[dstat1$arm == 1 & dstat1$binary == 0]),
                                                                       x = rnorm(n = 20,mean = 1.7, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = "darkgreen")

points(dstat1$count_mos[dstat1$arm == 0 & dstat1$binary == 1] ~ sample(size = length(dstat1$count_mos[dstat1$arm == 0 & dstat1$binary == 1]),
                                                                       x = rnorm(n = 20,mean = 2.7, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("darkblue",0.8))


points(dstat1$count_mos[dstat1$arm == 1 & dstat1$binary == 1] ~ sample(size = length(dstat1$count_mos[dstat1$arm == 1 & dstat1$binary == 1]),
                                                                       x = rnorm(n = 20,mean = 3.7, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("aquamarine2",0.4))


## Village 3
## subset for site 3
dstat3 = subset(dstat, dstat$site != 1 & dstat$site != 4 & dstat$species == "fun")
stan3 <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstat3,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  
100*(1-exp(stan3$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stan3$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stan3$coefficients[4])) ##difference in difference

100*(1-exp(stan3$coefficients[2]-stan3$ses[2])) ##difference in difference
100*(1-exp(stan3$coefficients[2]+stan3$ses[2])) ##difference in difference
100*(1-exp(stan3$coefficients[4]-stan3$ses[4])) ##difference in difference
100*(1-exp(stan3$coefficients[4]+stan3$ses[4])) ##difference in difference

mid_b1 = c(exp(posterior_interval(stan3, prob = 0.5)[1,1]),
           exp(posterior_interval(stan3, prob = 0.5)[1,1] + posterior_interval(stan3, prob = 0.5)[3,1]),
           exp(posterior_interval(stan3, prob = 0.5)[1,1] + posterior_interval(stan3, prob = 0.5)[2,1]),
           exp(posterior_interval(stan3, prob = 0.5)[1,1] + posterior_interval(stan3, prob = 0.5)[4,1]))

lower_b1 = c(exp(posterior_interval(stan3, prob = 0.95)[1,1]),
             exp(posterior_interval(stan3, prob = 0.95)[1,1] + posterior_interval(stan3, prob = 0.95)[3,1]),
             exp(posterior_interval(stan3, prob = 0.95)[1,1] + posterior_interval(stan3, prob = 0.95)[2,1]),
             exp(posterior_interval(stan3, prob = 0.95)[1,1] + posterior_interval(stan3, prob = 0.95)[4,1]))

upper_b1 = c(exp(posterior_interval(stan3, prob = 0.95)[1,2]),
             exp(posterior_interval(stan3, prob = 0.95)[1,2] + posterior_interval(stan3, prob = 0.95)[3,2]),
             exp(posterior_interval(stan3, prob = 0.95)[1,2] + posterior_interval(stan3, prob = 0.95)[2,2]),
             exp(posterior_interval(stan3, prob = 0.95)[1,2] + posterior_interval(stan3, prob = 0.95)[4,2]))

segsv1 = c(1,2,3,4)
for(i in 1:4){
  
  segments(x0 = segsv1[i],
           x1 = segsv1[i],
           y0 = lower_b1[i],
           y1 = upper_b1[i],
           lwd=4,col="red")
  
}

points(mid_b1 ~ segsv1,pch=19,cex = 1.2, col = "red")

points(dstat3$count_mos[dstat3$arm == 0 & dstat3$binary == 0] ~ sample(size = length(dstat3$count_mos[dstat3$arm == 0 & dstat3$binary == 0]),
                                                                       x = rnorm(n = 20,mean = 1, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("grey40",0.6))

points(dstat3$count_mos[dstat3$arm == 1 & dstat3$binary == 0] ~ sample(size = length(dstat3$count_mos[dstat3$arm == 1 & dstat3$binary == 0]),
                                                                       x = rnorm(n = 20,mean = 2, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("darkgreen",0.6))

points(dstat3$count_mos[dstat3$arm == 0 & dstat3$binary == 1] ~ sample(size = length(dstat3$count_mos[dstat3$arm == 0 & dstat3$binary == 1]),
                                                                       x = rnorm(n = 20,mean = 3, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("darkblue",0.5))


points(dstat3$count_mos[dstat3$arm == 1 & dstat3$binary == 1] ~ sample(size = length(dstat3$count_mos[dstat3$arm == 1 & dstat1$binary == 1]),
                                                                       x = rnorm(n = 20,mean = 4, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("aquamarine3",0.4))


## Village 4
## subset for site 4
dstat4 = subset(dstat, dstat$site != 1 & dstat$site != 3 & dstat$species == "fun")
stan4 <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstat4,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  
100*(1-exp(stan4$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stan4$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stan4$coefficients[4])) ##difference in difference

100*(1-exp(stan4$coefficients[2]-stan4$ses[2])) ##difference in difference
100*(1-exp(stan4$coefficients[2]+stan4$ses[2])) ##difference in difference
100*(1-exp(stan4$coefficients[4]-stan4$ses[4])) ##difference in difference
100*(1-exp(stan4$coefficients[4]+stan4$ses[4])) ##difference in difference

mid_b1 = c(exp(posterior_interval(stan4, prob = 0.5)[1,1]),
           exp(posterior_interval(stan4, prob = 0.5)[1,1] + posterior_interval(stan4, prob = 0.5)[3,1]),
           exp(posterior_interval(stan4, prob = 0.5)[1,1] + posterior_interval(stan4, prob = 0.5)[2,1]),
           exp(posterior_interval(stan4, prob = 0.5)[1,1] + posterior_interval(stan4, prob = 0.5)[4,1]))

lower_b1 = c(exp(posterior_interval(stan4, prob = 0.95)[1,1]),
             exp(posterior_interval(stan4, prob = 0.95)[1,1] + posterior_interval(stan4, prob = 0.95)[3,1]),
             exp(posterior_interval(stan4, prob = 0.95)[1,1] + posterior_interval(stan4, prob = 0.95)[2,1]),
             exp(posterior_interval(stan4, prob = 0.95)[1,1] + posterior_interval(stan4, prob = 0.95)[4,1]))

upper_b1 = c(exp(posterior_interval(stan4, prob = 0.95)[1,2]),
             exp(posterior_interval(stan4, prob = 0.95)[1,2] + posterior_interval(stan4, prob = 0.95)[3,2]),
             exp(posterior_interval(stan4, prob = 0.95)[1,2] + posterior_interval(stan4, prob = 0.95)[2,2]),
             exp(posterior_interval(stan4, prob = 0.95)[1,2] + posterior_interval(stan4, prob = 0.95)[4,2]))

segsv1 = c(1.3,2.3,3.3,4.3)
for(i in 1:4){
  
  segments(x0 = segsv1[i],
           x1 = segsv1[i],
           y0 = lower_b1[i],
           y1 = upper_b1[i],
           lwd=4,col="darkred")
  
}

points(mid_b1 ~ segsv1,pch=19,cex = 1.2, col = "darkred")

points(dstat4$count_mos[dstat4$arm == 0 & dstat4$binary == 0] ~ sample(size = length(dstat4$count_mos[dstat4$arm == 0 & dstat4$binary == 0]),
                                                                       x = rnorm(n = 20,mean = 1.3, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("grey40",0.3))

points(dstat4$count_mos[dstat4$arm == 1 & dstat4$binary == 0] ~ sample(size = length(dstat4$count_mos[dstat4$arm == 1 & dstat4$binary == 0]),
                                                                       x = rnorm(n = 20,mean = 2.3, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("darkgreen",0.3))

points(dstat4$count_mos[dstat4$arm == 0 & dstat4$binary == 1] ~ sample(size = length(dstat4$count_mos[dstat4$arm == 0 & dstat4$binary == 1]),
                                                                       x = rnorm(n = 20,mean = 3.3, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("darkblue",0.2))


points(dstat4$count_mos[dstat4$arm == 1 & dstat4$binary == 1] ~ sample(size = length(dstat4$count_mos[dstat4$arm == 1 & dstat1$binary == 1]),
                                                                       x = rnorm(n = 20,mean = 4.3, sd = 0.04),
                                                                       replace = TRUE),
       pch=19,col = adegenet::transp("aquamarine3",0.8))

legend("topright",legend = c("Musilonge",
                             "Kezege",
                             "Wamondo"),
       col = c("orange","red","darkred"),
       pch=19,border = NA, bty = "n")


par(xpd=NA,cex = 1.11)

text(x = -5, y = 125,"(A)",cex=0.8)
text(x = -0.05, y = 125,"(B)",cex=0.8)
text(x = -5, y = 55,"(C)",cex=0.8)
text(x = -0.05, y = 55,"(D)",cex=0.8)










## Village 1
## subset for site 1
dstat1f = subset(dstat, dstat$site != 3 & dstat$site != 4 & dstat$species == "fun")
stan1f <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstat1f,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  
100*(1-exp(stan1f$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stan1f$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stan1f$coefficients[4])) ##difference in difference

100*(1-exp(stan1f$coefficients[2]-stan1f$ses[2])) ##difference in difference
100*(1-exp(stan1f$coefficients[2]+stan1f$ses[2])) ##difference in difference
100*(1-exp(stan1f$coefficients[4]-stan1f$ses[4])) ##difference in difference
100*(1-exp(stan1f$coefficients[4]+stan1f$ses[4])) ##difference in difference

## Village 3
## subset for site 3
dstat3f = subset(dstat, dstat$site != 1 & dstat$site != 4 & dstat$species == "fun")
stan3f <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstat3f,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  
100*(1-exp(stan3f$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stan3f$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stan3f$coefficients[4])) ##difference in difference

100*(1-exp(stan3f$coefficients[2]-stan3f$ses[2])) ##difference in difference
100*(1-exp(stan3f$coefficients[2]+stan3f$ses[2])) ##difference in difference
100*(1-exp(stan3f$coefficients[4]-stan3f$ses[4])) ##difference in difference
100*(1-exp(stan3f$coefficients[4]+stan3f$ses[4])) ##difference in difference

## Village 4
## subset for site 4
dstat4f = subset(dstat, dstat$site != 1 & dstat$site != 3 & dstat$species == "fun")
stan4f <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ## species not important so removed
  data = dstat4f,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  
100*(1-exp(stan4f$coefficients[2])) ##median_change_after_nonLSM
100*(1-exp(stan4f$coefficients[3])) ##median_change_after_LSM
100*(1-exp(stan4f$coefficients[4])) ##difference in difference

100*(1-exp(stan4f$coefficients[2]-stan4f$ses[2])) ##difference in difference
100*(1-exp(stan4f$coefficients[2]+stan4f$ses[2])) ##difference in difference
100*(1-exp(stan4f$coefficients[4]-stan4f$ses[4])) ##difference in difference
100*(1-exp(stan4f$coefficients[4]+stan4f$ses[4])) ##difference in difference

phi_lsm_fun = c(as.numeric(exp(stan1f$coefficients[4])),
                1, 
                as.numeric(exp(stan3f$coefficients[4])), 
                as.numeric(exp(stan4f$coefficients[4])),
                1, 
                1)
rho_alpha = -4
rho_beta = 0.057

## lower bound
phi_lsm_gam_low = c(as.numeric(exp(stan1$coefficients[4]-stan1$ses[4])),
                    1, 
                    as.numeric(exp(stan3$coefficients[4]-stan3$ses[4])), 
                    as.numeric(exp(stan4$coefficients[4]-stan4$ses[4])),
                    1, 
                    1)
phi_lsm_ara_low = c(1, 1, 1, 1, 1, 1)
phi_lsm_fun_low = c(as.numeric(exp(stan1f$coefficients[4]-stan1f$ses[4])),
                    1, 
                    as.numeric(exp(stan3f$coefficients[4]-stan3f$ses[4])), 
                    as.numeric(exp(stan4f$coefficients[4]-stan4f$ses[4])),
                    1, 
                    1)


## upper bound
phi_lsm_gam_upp = c(as.numeric(exp(stan1$coefficients[4]+stan1$ses[4])),
                    1, 
                    as.numeric(exp(stan3$coefficients[4]+stan3$ses[4])), 
                    1,#as.numeric(exp(stan4$coefficients[4]+stan4$ses[4])),
                    1, 
                    1)
phi_lsm_ara_upp = c(1, 1, 1, 1, 1, 1)
phi_lsm_fun_upp = c(as.numeric(exp(stan1f$coefficients[4]+stan1f$ses[4])),
                    1, 
                    as.numeric(exp(stan3f$coefficients[4]+stan3f$ses[4])), 
                    as.numeric(exp(stan4f$coefficients[4]+stan4f$ses[4])),
                    1, 
                    1)

