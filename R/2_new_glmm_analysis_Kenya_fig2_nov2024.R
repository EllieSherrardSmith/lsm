###################################
##
## glm analysis
##
##################################
# install.packages(c("rstanarm",
#                    "bayesplot", "ggplot2", "broom"))
library("rstanarm")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library("bayesplot")
library("ggplot2")
library("broom")

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
dstat$site_score = rep(c("b","a","c","d","a","a"), each = 60)
dstat$sitea = rep(c("a","b","c","d","e","f"), each = 60)
dstat$arm = rep(c(1,0,1,1,0,0),each = 60)          ## LSM or not
dstat$species = rep(rep(c("gam","fun"),each=30),6) ##
dstat$month = rep(1:30,12)
dstat$binary = rep(rep(c(0,1),c(18,12)),12)        ## before and after the intervention
dstat$binary6to12 = rep(rep(c(0,99,1),c(18,6,6)),12)        ## before and after the intervention

head(dstat)


## In the absence of covariates >
stan0 <- rstanarm::stan_glm.nb(
  count_mos ~ binary * arm, ##
  data = dstat,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
) 
summary(stan0)
plot(stan0)
plot(stan0, plotfun = "trace")
prior_summary(stan0)
pp_check(stan0, plotfun = "stat", stat = "mean")

##################################################
##
## THIS MODEL IS THE ONE

stan4 <- rstanarm::stan_glmer.nb(
  count_mos ~ binary * arm * species + (1|sitea), ##
  data = dstat,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 2000 # for speed of example only
)
print(stan4)
plot(stan4)
plot(stan4, plotfun = "trace")
# bayesplot::color_scheme_set("pink")
(trace <- plot(stan4, "trace", pars = "(Intercept)",ylab = "Intercept"))
(trace <- plot(stan4, "trace", pars = "binary"))
(trace <- plot(stan4, "trace", pars = "arm"))
(trace <- plot(stan4, "trace", pars = "speciesgam"))
(trace <- plot(stan4, "trace", pars = "binary:arm"))
(trace <- plot(stan4, "trace", pars = "binary:speciesgam"))
(trace <- plot(stan4, "trace", pars = "arm:speciesgam"))
(trace <- plot(stan4, "trace", pars = "binary:arm:speciesgam"))
(trace <- plot(stan4, "trace", pars = c("b[(Intercept) sitea:a]",
                                        "b[(Intercept) sitea:b]",
                                        "b[(Intercept) sitea:c]",
                                        "b[(Intercept) sitea:d]",
                                        "b[(Intercept) sitea:e]",
                                        "b[(Intercept) sitea:f]")))
(trace <- plot(stan4, "trace", pars = c("reciprocal_dispersion",
                                        "Sigma[sitea:(Intercept),(Intercept)]")))

summary(stan4)
priors = prior_summary(stan4)
names(priors)
priors$prior$scale
priors$prior$adjusted_scale
pp_check(stan4, plotfun = "stat", stat = "mean")
pp_check(stan4, plotfun = "dens_overlay")
yy_pp = posterior_predict(stan4)
n_sims = nrow(yy_pp)
subset = sample(n_sims, 100)
ppc_dens_overlay(log10(dstat$count_mos+1), log10(yy_pp[subset,]+1))

set00 = which(dstat$arm == 0 & dstat$binary == 0 & dstat$species == "fun")
set01 = which(dstat$arm == 0 & dstat$binary == 1 & dstat$species == "fun")
set10 = which(dstat$arm == 1 & dstat$binary == 0 & dstat$species == "fun")
set11 = which(dstat$arm == 1 & dstat$binary == 1 & dstat$species == "fun")
n_sims = nrow(yy_pp)
subset = sample(n_sims, 100)
ppc_dens_overlay(log10(dstat$count_mos[set00]+1), log10(yy_pp[subset,set00]+1))
ppc_dens_overlay(log10(dstat$count_mos[set01]+1), log10(yy_pp[subset,set01]+1))
ppc_dens_overlay(log10(dstat$count_mos[set10]+1), log10(yy_pp[subset,set10]+1))
ppc_dens_overlay(log10(dstat$count_mos[set11]+1), log10(yy_pp[subset,set11]+1))
ppc_dens_overlay(log10(dstat$count_mos+1), log10(yy_pp[subset,]+1))

set00 = which(dstat$arm == 0 & dstat$binary == 0 & dstat$species == "gam")
set01 = which(dstat$arm == 0 & dstat$binary == 1 & dstat$species == "gam")
set10 = which(dstat$arm == 1 & dstat$binary == 0 & dstat$species == "gam")
set11 = which(dstat$arm == 1 & dstat$binary == 1 & dstat$species == "gam")
ppc_dens_overlay(log10(dstat$count_mos[set00]+1), log10(yy_pp[subset,set00]+1))
ppc_dens_overlay(log10(dstat$count_mos[set01]+1), log10(yy_pp[subset,set01]+1))
ppc_dens_overlay(log10(dstat$count_mos[set10]+1), log10(yy_pp[subset,set10]+1))
ppc_dens_overlay(log10(dstat$count_mos[set11]+1), log10(yy_pp[subset,set11]+1))

# library(shinystan)

## equation
pre_or_post = 0 ## 0 or 1
arm_val = 0 ## 0 controls and 1 LSM
which_species = 0 ## 0 or 1 where 0 is funestus, 1 is gambiae
translate_f = function(pre_or_post,
                       arm_val,
                       which_species){
  cnt = (stan4$coefficients[1] + 
              stan4$coefficients[2] * pre_or_post +
              stan4$coefficients[3] * arm_val +
              stan4$coefficients[4] * which_species +
              stan4$coefficients[5] * pre_or_post * arm_val +
              stan4$coefficients[6] * pre_or_post * which_species +
              stan4$coefficients[7] * arm_val * which_species +
              stan4$coefficients[8] * pre_or_post * arm_val * which_species)
  return(as.numeric(cnt))
}
B_pre_con_fun = translate_f(pre_or_post = 0,arm_val = 0,which_species = 0 )
B_pre_con_gam = translate_f(pre_or_post = 0,arm_val = 0,which_species = 1 )
A_pre_trt_fun = translate_f(pre_or_post = 0,arm_val = 1,which_species = 0 )
A_pre_trt_gam = translate_f(pre_or_post = 0,arm_val = 1,which_species = 1 )

D_post_con_fun = translate_f(pre_or_post = 1,arm_val = 0,which_species = 0 )
D_post_con_gam = translate_f(pre_or_post = 1,arm_val = 0,which_species = 1 )
C_post_trt_fun = translate_f(pre_or_post = 1,arm_val = 1,which_species = 0 )
C_post_trt_gam = translate_f(pre_or_post = 1,arm_val = 1,which_species = 1 )

(exp(B_pre_con_fun) - exp(D_post_con_fun))/exp(B_pre_con_fun)  ## 100*(1-exp(stan4$coefficients[2]))
(exp(B_pre_con_fun) - exp(A_pre_trt_fun))/exp(B_pre_con_fun)  ## 100*(1-exp(stan4$coefficients[3]))

cc2 = exp((C_post_trt_fun-A_pre_trt_fun)-(D_post_con_fun-B_pre_con_fun))
100*(1-cc2)
## 100*(1-exp(stan4$coefficients[5]))

(exp(B_pre_con_gam) - exp(D_post_con_gam))/exp(B_pre_con_gam)   ## 100*(1-exp(stan4$coefficients[2]+stan4$coefficients[6]))
(exp(B_pre_con_gam) - exp(A_pre_trt_gam))/exp(B_pre_con_gam) ## 100*(1-exp(stan4$coefficients[3]+stan4$coefficients[7]))

cc2 = exp((C_post_trt_gam-A_pre_trt_gam)-(D_post_con_gam-B_pre_con_gam))
100*(1-cc2)

### Generate uncertainty
translate_uncert_f = function(num_reps){
  B_pre_con_fun = B_pre_con_gam = A_pre_trt_fun = A_pre_trt_gam = 
    D_post_con_fun = D_post_con_gam = C_post_trt_fun = C_post_trt_gam = numeric(num_reps)
  
  translate_f = function(pre_or_post,
                         arm_val,
                         which_species,rr){
    cnt = (sims4[rr,1] + 
                sims4[rr,2] * pre_or_post +
                sims4[rr,3] * arm_val +
                sims4[rr,4] * which_species +
                sims4[rr,5] * pre_or_post * arm_val +
                sims4[rr,6] * pre_or_post * which_species +
                sims4[rr,7] * arm_val * which_species +
                sims4[rr,8] * pre_or_post * arm_val * which_species)
    
    return(as.numeric(cnt))
  }
  sims4 = as.matrix(stan4)
  n_sims = nrow(sims4)
  random_draw = sample(n_sims,num_reps)
  for(j in 1:num_reps){
    B_pre_con_fun[j] = translate_f(pre_or_post = 0,arm_val = 0,which_species = 0, rr=random_draw[j])
    B_pre_con_gam[j] = translate_f(pre_or_post = 0,arm_val = 0,which_species = 1, rr=random_draw[j])
    A_pre_trt_fun[j] = translate_f(pre_or_post = 0,arm_val = 1,which_species = 0, rr=random_draw[j])
    A_pre_trt_gam[j] = translate_f(pre_or_post = 0,arm_val = 1,which_species = 1, rr=random_draw[j])
    
    D_post_con_fun[j] = translate_f(pre_or_post = 1,arm_val = 0,which_species = 0, rr=random_draw[j])
    D_post_con_gam[j] = translate_f(pre_or_post = 1,arm_val = 0,which_species = 1, rr=random_draw[j])
    C_post_trt_fun[j] = translate_f(pre_or_post = 1,arm_val = 1,which_species = 0, rr=random_draw[j])
    C_post_trt_gam[j] = translate_f(pre_or_post = 1,arm_val = 1,which_species = 1, rr=random_draw[j])
  }
  cc1 = exp((C_post_trt_fun-A_pre_trt_fun)-(D_post_con_fun-B_pre_con_fun))
  did_funestus = 100*(1-cc1)## 100*(1-exp(stan4$coefficients[5]))
    
  cc2 = exp((C_post_trt_gam-A_pre_trt_gam)-(D_post_con_gam-B_pre_con_gam))
  did_gambiae = 100*(1-cc2) ## 100*(1-exp(stan4$coefficients[5]+stan4$coefficients[8]))
  
  return(data.frame(did_funestus = did_funestus,
                    did_gambiae = did_gambiae,
                    B_pre_con_fun = exp(B_pre_con_fun),
                    A_pre_trt_fun = exp(A_pre_trt_fun),
                    D_post_con_fun = exp(D_post_con_fun),
                    C_post_trt_fun = exp(C_post_trt_fun),
                    B_pre_con_gam = exp(B_pre_con_gam),
                    A_pre_trt_gam = exp(A_pre_trt_gam),
                    D_post_con_gam = exp(D_post_con_gam),
                    C_post_trt_gam = exp(C_post_trt_gam)))
}

rrr = translate_uncert_f(num_reps = 1000)
quantile(rrr[,1],c(0.05,0.5,0.95))
quantile(rrr[,2],c(0.05,0.5,0.95))

par(mar = c(3,4,1,1))
par(mfrow = c(1,3))
boxplot(sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 0 & dstat$arm == 0]),
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 0 & dstat$arm == 1]),
        NA,
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 1 & dstat$arm == 0]),
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 1 & dstat$arm == 1]), 
        yaxt = "n",border = "grey",
        ylim = sqrt(c(0,120)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")
text(4,sqrt(92),"An. funestus")
axis(2, las = 2, at = sqrt(seq(20,120,20)),labels = seq(20,120,20))
axis(2, las = 2, at = sqrt(c(0,5,12)),labels = c(0,5,12))

abline(v=3,lty=1,col="grey34")
mtext("Pre-intervention                  Post-intervention",
      cex = 0.8,
      side = 1, line = 1)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5))),
         y1 = sqrt(as.numeric(quantile(rrr[,5],0.5))),
         col="darkred")
diff_before = sqrt(as.numeric(quantile(rrr[,4],0.5))) - sqrt(as.numeric(quantile(rrr[,3],0.5)))
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before,
         y1 = sqrt(as.numeric(quantile(rrr[,5],0.5)))+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before) /1.5
segments(x0 = 3,
         x1 = 4.5,
         y0 = sqrt(y2)-0.2,
         y1 = sqrt(as.numeric(quantile(rrr[,6],0.5))),
         lty=2,col="darkblue")

points(sqrt(c(as.numeric(quantile(rrr[,3],0.5)),
           as.numeric(quantile(rrr[,4],0.5)),
           as.numeric(quantile(rrr[,5],0.5)),
           as.numeric(quantile(rrr[,6],0.5)))) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,3],0.05)),
                 as.numeric(quantile(rrr[,4],0.05)),
                 NA,
                 as.numeric(quantile(rrr[,5],0.05)),
                 as.numeric(quantile(rrr[,6],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,3],0.95)),
               as.numeric(quantile(rrr[,4],0.95)),
               NA,
               as.numeric(quantile(rrr[,5],0.95)),
               as.numeric(quantile(rrr[,6],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}
legend("right",legend = c("Control","LSM"),
       col = c("red","blue"),pch=15,bty="n")
text(1,sqrt(120),"Kenya")
##################################
##
## Repeat for gambiae
##
################################
boxplot(sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 0 & dstat$arm == 0]),
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 0 & dstat$arm == 1]),
        NA,
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 1 & dstat$arm == 0]),
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 1 & dstat$arm == 1]), 
        ylim = sqrt(c(0,120)),yaxt = "n",border = "grey",
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")
text(4,sqrt(92),"An. gambiae")
axis(2, las = 2, at = sqrt(seq(20,120,20)),labels = seq(20,120,20))
axis(2, las = 2, at = sqrt(c(0,5,12)),labels = c(0,5,12))

abline(v=3,lty=1,col="grey34")
mtext("Pre-intervention         Post-intervention",
      cex = 0.8,
      side = 1, line = 1)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,7],0.5))),
         y1 = sqrt(as.numeric(quantile(rrr[,9],0.5))),
         col="darkred")
diff_before = sqrt(as.numeric(quantile(rrr[,8],0.5))) - sqrt(as.numeric(quantile(rrr[,7],0.5)))
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,7],0.5)))+diff_before,
         y1 = sqrt(as.numeric(quantile(rrr[,9],0.5)))+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(as.numeric(quantile(rrr[,7],0.5)))+diff_before) /1.5
segments(x0 = 3,
         x1 = 4.5,
         y0 = 3.15,#sqrt(y2)+0.22,
         y1 = sqrt(as.numeric(quantile(rrr[,10],0.5))),
         lty=2,col="darkblue")

points(sqrt(c(as.numeric(quantile(rrr[,7],0.5)),
              as.numeric(quantile(rrr[,8],0.5)),
              as.numeric(quantile(rrr[,9],0.5)),
              as.numeric(quantile(rrr[,10],0.5)))) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,7],0.05)),
                    as.numeric(quantile(rrr[,8],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,9],0.05)),
                    as.numeric(quantile(rrr[,10],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,7],0.95)),
                    as.numeric(quantile(rrr[,8],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,9],0.95)),
                    as.numeric(quantile(rrr[,10],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}
text(1,sqrt(120),"Kenya")

##################################################################
##
## Next we examine the model to look at village level effects
##
#########################################
stan5 <- rstanarm::stan_glm.nb(
  count_mos ~ binary * site_score * species, ##
  data = dstat,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 2000 # for speed of example only
)

print(stan5)
plot(stan5)
plot(stan5, plotfun = "trace")

(trace <- plot(stan5, "trace", pars = "(Intercept)",ylab = "Intercept"))
(trace <- plot(stan5, "trace", pars = "binary"))
(trace <- plot(stan5, "trace", pars = c("site_scoreb",
                                        "site_scorec",
                                        "site_scored")))
(trace <- plot(stan5, "trace", pars = "speciesgam"))
(trace <- plot(stan5, "trace", pars = c("binary:site_scoreb",
                                        "binary:site_scorec",
                                        "binary:site_scored")) )
(trace <- plot(stan5, "trace", pars = "binary:speciesgam"))
(trace <- plot(stan5, "trace", pars = c("site_scoreb:speciesgam",
                                        "site_scorec:speciesgam",
                                        "site_scored:speciesgam")) )
(trace <- plot(stan5, "trace", pars = c("binary:site_scoreb:speciesgam",
                                        "binary:site_scorec:speciesgam",
                                        "binary:site_scored:speciesgam")) )
(trace <- plot(stan5, "trace", pars = c("reciprocal_dispersion")))


summary(stan5)
priors = prior_summary(stan5)
names(priors)
priors$prior$scale
priors$prior$adjusted_scale

pp_check(stan5, plotfun = "stat", stat = "mean")
pp_check(stan5, plotfun = "dens_overlay")
yy_pp = posterior_predict(stan5)
n_sims = nrow(yy_pp)
subset = sample(n_sims, 100)
ppc_dens_overlay(log10(dstat$count_mos+1), log10(yy_pp[subset,]+1))

set00 = which(dstat$arm == 0 & dstat$binary == 0 & dstat$species == "fun")
set01 = which(dstat$arm == 0 & dstat$binary == 1 & dstat$species == "fun")
set10 = which(dstat$arm == 1 & dstat$binary == 0 & dstat$species == "fun")
set11 = which(dstat$arm == 1 & dstat$binary == 1 & dstat$species == "fun")
n_sims = nrow(yy_pp)
subset = sample(n_sims, 100)
ppc_dens_overlay(log10(dstat$count_mos[set00]+1), log10(yy_pp[subset,set00]+1))
ppc_dens_overlay(log10(dstat$count_mos[set01]+1), log10(yy_pp[subset,set01]+1))
ppc_dens_overlay(log10(dstat$count_mos[set10]+1), log10(yy_pp[subset,set10]+1))
ppc_dens_overlay(log10(dstat$count_mos[set11]+1), log10(yy_pp[subset,set11]+1))
ppc_dens_overlay(log10(dstat$count_mos+1), log10(yy_pp[subset,]+1))

set00 = which(dstat$arm == 0 & dstat$binary == 0 & dstat$species == "gam")
set01 = which(dstat$arm == 0 & dstat$binary == 1 & dstat$species == "gam")
set10 = which(dstat$arm == 1 & dstat$binary == 0 & dstat$species == "gam")
set11 = which(dstat$arm == 1 & dstat$binary == 1 & dstat$species == "gam")
ppc_dens_overlay(log10(dstat$count_mos[set00]+1), log10(yy_pp[subset,set00]+1))
ppc_dens_overlay(log10(dstat$count_mos[set01]+1), log10(yy_pp[subset,set01]+1))
ppc_dens_overlay(log10(dstat$count_mos[set10]+1), log10(yy_pp[subset,set10]+1))
ppc_dens_overlay(log10(dstat$count_mos[set11]+1), log10(yy_pp[subset,set11]+1))

# library(shinystan)

## equation
pre_or_post = 0 ## 0 or 1
arm_val = 0 ## 0 controls and 1 LSM
which_species = 0 ## 0 or 1 where 0 is funestus, 1 is gambiae
translate2_f = function(pre_or_post,
                       vill1,vill3,vill4, ## either 0 (con), or 1 (villages with LSM)
                       which_species){
  cnt = (stan5$coefficients[1] + 
              stan5$coefficients[2] * pre_or_post +
              stan5$coefficients[3] * vill1 +
              stan5$coefficients[4] * vill3 +
              stan5$coefficients[5] * vill4 +
              stan5$coefficients[6] * which_species +
              stan5$coefficients[7] * pre_or_post * vill1 +
              stan5$coefficients[8] * pre_or_post * vill3 +
              stan5$coefficients[9] * pre_or_post * vill4 +
              stan5$coefficients[10] * pre_or_post * which_species +
              stan5$coefficients[11] * vill1 * which_species +
              stan5$coefficients[12] * vill3 * which_species +
              stan5$coefficients[13] * vill4 * which_species +
              stan5$coefficients[14] * pre_or_post * vill1 * which_species +
              stan5$coefficients[15] * pre_or_post * vill3 * which_species +
              stan5$coefficients[16] * pre_or_post * vill4 * which_species)
  return(as.numeric(cnt))
}
B_pre_con_fun = translate2_f(pre_or_post = 0,
                             vill1 = 0,vill3 = 0,vill4 = 0,
                             which_species = 0 )
B_pre_con_gam = translate2_f(pre_or_post = 0,
                             vill1 = 0,vill3 = 0,vill4 = 0,
                             which_species = 1 )
A1_pre_trt_fun = translate2_f(pre_or_post = 0,
                             vill1 = 1,vill3 = 0,vill4 = 0,
                             which_species = 0 )
A1_pre_trt_gam = translate2_f(pre_or_post = 0,
                             vill1 = 1,vill3 = 0,vill4 = 0,
                             which_species = 1 )
A3_pre_trt_fun = translate2_f(pre_or_post = 0,
                              vill1 = 0,vill3 = 1,vill4 = 0,
                              which_species = 0 )
A3_pre_trt_gam = translate2_f(pre_or_post = 0,
                              vill1 = 0,vill3 = 1,vill4 = 0,
                              which_species = 1 )
A4_pre_trt_fun = translate2_f(pre_or_post = 0,
                              vill1 = 0,vill3 = 0,vill4 = 1,
                              which_species = 0 )
A4_pre_trt_gam = translate2_f(pre_or_post = 0,
                              vill1 = 0,vill3 = 0,vill4 = 1,
                              which_species = 1 )

D_post_con_fun = translate2_f(pre_or_post = 1,
                             vill1 = 0,vill3 = 0,vill4 = 0,
                             which_species = 0 )
D_post_con_gam = translate2_f(pre_or_post = 1,
                             vill1 = 0,vill3 = 0,vill4 = 0,
                             which_species = 1 )

C1_post_trt_fun = translate2_f(pre_or_post = 1,
                             vill1 = 1,vill3 = 0,vill4 = 0,
                             which_species = 0 )
C1_post_trt_gam = translate2_f(pre_or_post = 1,
                             vill1 = 1,vill3 = 0,vill4 = 0,
                             which_species = 1 )
C3_post_trt_fun = translate2_f(pre_or_post = 1,
                              vill1 = 0,vill3 = 1,vill4 = 0,
                              which_species = 0 )
C3_post_trt_gam = translate2_f(pre_or_post = 1,
                              vill1 = 0,vill3 = 1,vill4 = 0,
                              which_species = 1 )
C4_post_trt_fun = translate2_f(pre_or_post = 1,
                              vill1 = 0,vill3 = 0,vill4 = 1,
                              which_species = 0 )
C4_post_trt_gam = translate2_f(pre_or_post = 1,
                              vill1 = 0,vill3 = 0,vill4 = 1,
                              which_species = 1 )

(exp(B_pre_con_fun) - exp(D_post_con_fun))/exp(B_pre_con_fun)  ## 100*(1-exp(stan5$coefficients[2]))
(exp(B_pre_con_fun) - exp(A1_pre_trt_fun))/exp(B_pre_con_fun)  ## 100*(1-exp(stan5$coefficients[3]))


cv1 = exp((C1_post_trt_fun-A1_pre_trt_fun)-(D_post_con_fun-B_pre_con_fun))
100*(1-cv1)
cv3 = exp((C3_post_trt_fun-A3_pre_trt_fun)-(D_post_con_fun-B_pre_con_fun))
100*(1-cv3)
cv4 = exp((C4_post_trt_fun-A4_pre_trt_fun)-(D_post_con_fun-B_pre_con_fun))
100*(1-cv4)

(B_pre_con_gam - D_post_con_gam)/B_pre_con_gam  ## 100*(1-exp(stan5$coefficients[2]+stan5$coefficients[6]))
(B_pre_con_gam - A_pre_trt_gam)/B_pre_con_gam  ## 100*(1-exp(stan5$coefficients[3]+stan5$coefficients[7]))

cgv1 = exp((C1_post_trt_gam-A1_pre_trt_gam)-(D_post_con_gam-B_pre_con_gam))
100*(1-cgv1)
cgv3 = exp((C3_post_trt_gam-A3_pre_trt_gam)-(D_post_con_gam-B_pre_con_gam))
100*(1-cgv3)
cgv4 = exp((C4_post_trt_gam-A4_pre_trt_gam)-(D_post_con_gam-B_pre_con_gam))
100*(1-cgv4)

### Generate uncertainty
translate2_uncert_f = function(num_reps){
  B_pre_con_fun = B_pre_con_gam = 
    A1_pre_trt_fun = A1_pre_trt_gam = 
    A3_pre_trt_fun = A3_pre_trt_gam = 
    A4_pre_trt_fun = A4_pre_trt_gam = 
    D_post_con_fun = D_post_con_gam = 
    C1_post_trt_fun = C1_post_trt_gam = 
    C3_post_trt_fun = C3_post_trt_gam = 
    C4_post_trt_fun = C4_post_trt_gam = numeric(num_reps)
  
  translate2_f = function(pre_or_post,
                         vill1,vill3,vill4, ## either 0 (con), or 1 (villages with LSM)
                         which_species,rr){
    cnt = (sims5[rr,1] + 
                sims5[rr,2] * pre_or_post +
                sims5[rr,3] * vill1 +
                sims5[rr,4] * vill3 +
                sims5[rr,5] * vill4 +
                sims5[rr,6] * which_species +
                sims5[rr,7] * pre_or_post * vill1 +
                sims5[rr,8] * pre_or_post * vill3 +
                sims5[rr,9] * pre_or_post * vill4 +
                sims5[rr,10] * pre_or_post * which_species +
                sims5[rr,11] * vill1 * which_species +
                sims5[rr,12] * vill3 * which_species +
                sims5[rr,13] * vill4 * which_species +
                sims5[rr,14] * pre_or_post * vill1 * which_species +
                sims5[rr,15] * pre_or_post * vill3 * which_species +
                sims5[rr,16] * pre_or_post * vill4 * which_species)
    
    return(as.numeric(cnt))
  }
  sims5 = as.matrix(stan5)
  n_sims = nrow(sims5)
  random_draw = sample(n_sims,num_reps)
  for(j in 1:num_reps){
    B_pre_con_fun[j] = translate2_f(pre_or_post = 0,vill1 = 0,vill3 = 0,vill4 = 0,which_species = 0, rr=random_draw[j])
    B_pre_con_gam[j] = translate2_f(pre_or_post = 0,vill1 = 0,vill3 = 0,vill4 = 0,which_species = 1, rr=random_draw[j])
    A1_pre_trt_fun[j] = translate2_f(pre_or_post = 0,vill1 = 1,vill3 = 0,vill4 = 0,which_species = 0, rr=random_draw[j])
    A1_pre_trt_gam[j] = translate2_f(pre_or_post = 0,vill1 = 1,vill3 = 0,vill4 = 0,which_species = 1, rr=random_draw[j])
    A3_pre_trt_fun[j] = translate2_f(pre_or_post = 0,vill1 = 0,vill3 = 1,vill4 = 0,which_species = 0, rr=random_draw[j])
    A3_pre_trt_gam[j] = translate2_f(pre_or_post = 0,vill1 = 0,vill3 = 1,vill4 = 0,which_species = 1, rr=random_draw[j])
    A4_pre_trt_fun[j] = translate2_f(pre_or_post = 0,vill1 = 0,vill3 = 0,vill4 = 1,which_species = 0, rr=random_draw[j])
    A4_pre_trt_gam[j] = translate2_f(pre_or_post = 0,vill1 = 0,vill3 = 0,vill4 = 1,which_species = 1, rr=random_draw[j])
    
    D_post_con_fun[j] = translate2_f(pre_or_post = 1,vill1 = 0,vill3 = 0,vill4 = 0,which_species = 0, rr=random_draw[j])
    D_post_con_gam[j] = translate2_f(pre_or_post = 1,vill1 = 0,vill3 = 0,vill4 = 0,which_species = 1, rr=random_draw[j])
    C1_post_trt_fun[j] = translate2_f(pre_or_post = 1,vill1 = 1,vill3 = 0,vill4 = 0,which_species = 0, rr=random_draw[j])
    C1_post_trt_gam[j] = translate2_f(pre_or_post = 1,vill1 = 1,vill3 = 0,vill4 = 0,which_species = 1, rr=random_draw[j])
    C3_post_trt_fun[j] = translate2_f(pre_or_post = 1,vill1 = 0,vill3 = 1,vill4 = 0,which_species = 0, rr=random_draw[j])
    C3_post_trt_gam[j] = translate2_f(pre_or_post = 1,vill1 = 0,vill3 = 1,vill4 = 0,which_species = 1, rr=random_draw[j])
    C4_post_trt_fun[j] = translate2_f(pre_or_post = 1,vill1 = 0,vill3 = 0,vill4 = 1,which_species = 0, rr=random_draw[j])
    C4_post_trt_gam[j] = translate2_f(pre_or_post = 1,vill1 = 0,vill3 = 0,vill4 = 1,which_species = 1, rr=random_draw[j])
  }
  
  cc1 = exp((C1_post_trt_fun-A1_pre_trt_fun)-(D_post_con_fun-B_pre_con_fun))
  did_funestusV1 = 100*(1-cc1)
  cc3 = exp((C3_post_trt_fun-A3_pre_trt_fun)-(D_post_con_fun-B_pre_con_fun))
  did_funestusV3 = 100*(1-cc3)
  cc4 = exp((C4_post_trt_fun-A4_pre_trt_fun)-(D_post_con_fun-B_pre_con_fun))
  did_funestusV1 = 100*(1-cc4)
  
  ccg1 = exp((C1_post_trt_gam-A1_pre_trt_gam)-(D_post_con_gam-B_pre_con_gam))
  did_gambiaeV1 = 100*(1-ccg1) ## 100*(1-exp(stan4$coefficients[5]+stan4$coefficients[8]))
  ccg3 = exp((C3_post_trt_gam-A3_pre_trt_gam)-(D_post_con_gam-B_pre_con_gam))
  did_gambiaeV3 = 100*(1-ccg3) ## 100*(1-exp(stan4$coefficients[5]+stan4$coefficients[8]))
  ccg4 = exp((C4_post_trt_gam-A4_pre_trt_gam)-(D_post_con_gam-B_pre_con_gam))
  did_gambiaeV4 = 100*(1-ccg4) ## 100*(1-exp(stan4$coefficients[5]+stan4$coefficients[8]))
  
  return(data.frame(did_funestusV1 = did_funestusV1,
                    did_funestusV3 = did_funestusV3,
                    did_funestusV4 = did_funestusV4,
                    did_gambiaeV1 = did_gambiaeV1,
                    did_gambiaeV3 = did_gambiaeV3,
                    did_gambiaeV4 = did_gambiaeV4,
                    
                    B_pre_con_fun = exp(B_pre_con_fun),
                    A1_pre_trt_fun = exp(A1_pre_trt_fun),
                    A3_pre_trt_fun = exp(A3_pre_trt_fun),
                    A4_pre_trt_fun = exp(A4_pre_trt_fun),
                    D_post_con_fun = exp(D_post_con_fun),
                    C1_post_trt_fun = exp(C1_post_trt_fun),
                    C3_post_trt_fun = exp(C3_post_trt_fun),
                    C4_post_trt_fun = exp(C4_post_trt_fun),
                    
                    B_pre_con_gam = exp(B_pre_con_gam),
                    A1_pre_trt_gam = exp(A1_pre_trt_gam),
                    A3_pre_trt_gam = exp(A3_pre_trt_gam),
                    A4_pre_trt_gam = exp(A4_pre_trt_gam),
                    D_post_con_gam = exp(D_post_con_gam),
                    C1_post_trt_gam = exp(C1_post_trt_gam),
                    C3_post_trt_gam = exp(C3_post_trt_gam),
                    C4_post_trt_gam = exp(C4_post_trt_gam)))
}

rrr = translate2_uncert_f(num_reps = 1000)
quantile(rrr[,1],c(0.05,0.5,0.95)) ## village 1, funestus
quantile(rrr[,2],c(0.05,0.5,0.95)) ## village 3, funestus
quantile(rrr[,3],c(0.05,0.5,0.95)) ## village 4, funestus

quantile(rrr[,4],c(0.05,0.5,0.95)) ## village 1, gam
quantile(rrr[,5],c(0.05,0.5,0.95)) ## village 3, gam
quantile(rrr[,6],c(0.05,0.5,0.95)) ## village 4, gam



##########################################
##
## Plotting the diff-in-diff for village level
par(mar = c(3,4,1,1))
par(mfrow = c(2,3))
boxplot(sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 0 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 0 & dstat$site_score == "b"]),
        NA,
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 1 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 1 & dstat$site_score == "b"]), 
        yaxt = "n",border = "grey",
        ylim = sqrt(c(0,120)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")
text(4,sqrt(92),"An. funestus")
text(4,sqrt(110),"Musilongo")
axis(2, las = 2, at = sqrt(seq(20,120,20)),labels = seq(20,120,20))
axis(2, las = 2, at = sqrt(c(0,5,12)),labels = c(0,5,12))

abline(v=3,lty=1,col="grey34")
mtext(c(        "Pre-intervention   Post-intervention" ),
      side = 1, line = 1,cex = 0.6)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_fun),
         y1 = sqrt(D_post_con_fun),
         col="darkred")
diff_before = sqrt(A1_pre_trt_fun) - sqrt(B_pre_con_fun)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_fun)+diff_before,
         y1 = sqrt(D_post_con_fun)+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(exp(B_pre_con_fun))+diff_before) /1.5
segments(x0 = 3,
         x1 = 4.5,
         y0 = sqrt(y2),
         y1 = sqrt(exp(C1_post_trt_fun)),
         lty=2,col="darkblue")

points(sqrt(c(B_pre_con_fun,
              A1_pre_trt_fun,
              D_post_con_fun,
              C1_post_trt_fun)) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,7],0.05)),
                    as.numeric(quantile(rrr[,8],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,11],0.05)),
                    as.numeric(quantile(rrr[,12],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,7],0.95)),
                    as.numeric(quantile(rrr[,8],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,11],0.95)),
                    as.numeric(quantile(rrr[,12],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}
legend("right",legend = c("Control","LSM"),
       col = c("red","blue"),pch=15,bty="n")


#################################
##
## Vill 3
boxplot(sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 0 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 0 & dstat$site_score == "c"]),
        NA,
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 1 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 1 & dstat$site_score == "c"]), 
        yaxt = "n",border = "grey",
        ylim = sqrt(c(0,120)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")
text(4,sqrt(92),"An. funestus")
text(4,sqrt(110),"Kezege")
axis(2, las = 2, at = sqrt(seq(20,120,20)),labels = seq(20,120,20))
axis(2, las = 2, at = sqrt(c(0,5,12)),labels = c(0,5,12))

abline(v=3,lty=1,col="grey34")
mtext(c(        "Pre-intervention   Post-intervention" ),
      side = 1, line = 1,cex = 0.6)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_fun),
         y1 = sqrt(D_post_con_fun),
         col="darkred")
diff_before = sqrt(A3_pre_trt_fun) - sqrt(B_pre_con_fun)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_fun)+diff_before,
         y1 = sqrt(D_post_con_fun)+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(B_pre_con_fun)+diff_before) /1.5
segments(x0 = 3,
         x1 = 4.5,
         y0 = sqrt(y2)-diff_before,
         y1 = sqrt(C3_post_trt_fun),
         lty=2,col="darkblue")

points(sqrt(c(B_pre_con_fun,
              A3_pre_trt_fun,
              D_post_con_fun,
              C3_post_trt_fun)) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,7],0.05)),
                    as.numeric(quantile(rrr[,9],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,11],0.05)),
                    as.numeric(quantile(rrr[,13],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,7],0.95)),
                    as.numeric(quantile(rrr[,9],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,11],0.95)),
                    as.numeric(quantile(rrr[,13],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}


#################################
##
## Vill 4
boxplot(sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 0 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 0 & dstat$site_score == "d"]),
        NA,
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 1 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "fun" & dstat$binary == 1 & dstat$site_score == "d"]), 
        yaxt = "n",border = "grey",
        ylim = sqrt(c(0,120)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")
text(4,sqrt(92),"An. funestus")
text(4,sqrt(110),"Wamondo")
axis(2, las = 2, at = sqrt(seq(20,120,20)),labels = seq(20,120,20))
axis(2, las = 2, at = sqrt(c(0,5,12)),labels = c(0,5,12))

abline(v=3,lty=1,col="grey34")
mtext(c(        "Pre-intervention   Post-intervention" ),
      side = 1, line = 1,cex = 0.6)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_fun),
         y1 = sqrt(D_post_con_fun),
         col="darkred")
diff_before = sqrt(A4_pre_trt_fun) - sqrt(B_pre_con_fun)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_fun)+diff_before,
         y1 = sqrt(D_post_con_fun)+diff_before,
         lty=1,col="darkblue")

segments(x0 = 3,
         x1 = 4.5,
         y0 = sqrt(0.6),
         y1 = sqrt(C4_post_trt_fun),
         lty=2,col="darkblue")

points(sqrt(c(B_pre_con_fun,
              A4_pre_trt_fun,
              D_post_con_fun,
              C4_post_trt_fun)) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,7],0.05)),
                    as.numeric(quantile(rrr[,10],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,11],0.05)),
                    as.numeric(quantile(rrr[,14],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,7],0.95)),
                    as.numeric(quantile(rrr[,10],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,11],0.95)),
                    as.numeric(quantile(rrr[,14],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}
#################################
##
## gambiae
boxplot(sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 0 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 0 & dstat$site_score == "b"]),
        NA,
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 1 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 1 & dstat$site_score == "b"]), 
        yaxt = "n",border = "grey",
        ylim = sqrt(c(0,120)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")
text(4,sqrt(92),"An. gambiae")
text(4,sqrt(110),"Musilongo")
axis(2, las = 2, at = sqrt(seq(20,120,20)),labels = seq(20,120,20))
axis(2, las = 2, at = sqrt(c(0,5,12)),labels = c(0,5,12))

abline(v=3,lty=1,col="grey34")
mtext(c(        "Pre-intervention   Post-intervention" ),
      side = 1, line = 1,cex = 0.6)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_gam),
         y1 = sqrt(D_post_con_gam),
         col="darkred")
diff_before = sqrt(A1_pre_trt_gam) - sqrt(B_pre_con_gam)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_gam)+diff_before,
         y1 = sqrt(D_post_con_gam)+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(B_pre_con_gam)+diff_before) /1.5
segments(x0 = 3,
         x1 = 4.5,
         y0 = sqrt(y2)+diff_before*1.28,
         y1 = sqrt(C1_post_trt_gam),
         lty=2,col="darkblue")

points(sqrt(c(B_pre_con_gam,
              A1_pre_trt_gam,
              D_post_con_gam,
              C1_post_trt_gam)) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,15],0.05)),
                    as.numeric(quantile(rrr[,16],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,19],0.05)),
                    as.numeric(quantile(rrr[,20],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,15],0.95)),
                    as.numeric(quantile(rrr[,16],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,19],0.95)),
                    as.numeric(quantile(rrr[,20],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}
#################################
##
## Vill 3
boxplot(sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 0 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 0 & dstat$site_score == "c"]),
        NA,
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 1 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 1 & dstat$site_score == "c"]), 
        yaxt = "n",border = "grey",
        ylim = sqrt(c(0,120)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")
text(4,sqrt(92),"An. gambiae")
text(4,sqrt(110),"Kezege")
axis(2, las = 2, at = sqrt(seq(20,120,20)),labels = seq(20,120,20))
axis(2, las = 2, at = sqrt(c(0,5,12)),labels = c(0,5,12))

abline(v=3,lty=1,col="grey34")
mtext(c(        "Pre-intervention     Post-intervention" ),
      side = 1, line = 1,cex=0.6)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_gam),
         y1 = sqrt(D_post_con_gam),
         col="darkred")
diff_before = sqrt(A3_pre_trt_gam) - sqrt(B_pre_con_gam)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_gam)+diff_before,
         y1 = sqrt(D_post_con_gam)+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(B_pre_con_gam)+diff_before) /1.5
segments(x0 = 3,
         x1 = 4.5,
         y0 = sqrt(y2)+diff_before+0.75,
         y1 = sqrt(C3_post_trt_gam),
         lty=2,col="darkblue")

points(sqrt(c(B_pre_con_gam,
              A3_pre_trt_gam,
              D_post_con_gam,
              C3_post_trt_gam)) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,15],0.05)),
                    as.numeric(quantile(rrr[,17],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,19],0.05)),
                    as.numeric(quantile(rrr[,21],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,15],0.95)),
                    as.numeric(quantile(rrr[,17],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,19],0.95)),
                    as.numeric(quantile(rrr[,21],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}



boxplot(sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 0 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 0 & dstat$site_score == "d"]),
        NA,
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 1 & dstat$site_score == "a"]),
        sqrt(dstat$count_mos[dstat$species == "gam" & dstat$binary == 1 & dstat$site_score == "d"]), 
        yaxt = "n",border = "grey",
        ylim = sqrt(c(0,120)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")
text(4,sqrt(92),"An. gambiae")
text(4,sqrt(110),"Wamondo")
axis(2, las = 2, at = sqrt(seq(20,120,20)),labels = seq(20,120,20))
axis(2, las = 2, at = sqrt(c(0,5,12)),labels = c(0,5,12))

abline(v=3,lty=1,col="grey34")
mtext(c(        "Pre-intervention   Post-intervention" ),
      side = 1, line = 1,cex = 0.6)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_gam),
         y1 = sqrt(D_post_con_gam),
         col="darkred")
diff_before = sqrt(A4_pre_trt_gam) - sqrt(B_pre_con_gam)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(B_pre_con_gam)+diff_before,
         y1 = sqrt(D_post_con_gam)+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(B_pre_con_gam)+diff_before) /1.5
segments(x0 = 3,
         x1 = 4.5,
         y0 = sqrt(y2)-2,
         y1 = sqrt(C4_post_trt_gam),
         lty=2,col="darkblue")

points(sqrt(c(B_pre_con_gam,
              A4_pre_trt_gam,
              D_post_con_gam,
              C4_post_trt_gam)) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,15],0.05)),
                    as.numeric(quantile(rrr[,18],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,19],0.05)),
                    as.numeric(quantile(rrr[,22],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,15],0.95)),
                    as.numeric(quantile(rrr[,18],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,19],0.95)),
                    as.numeric(quantile(rrr[,22],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}


