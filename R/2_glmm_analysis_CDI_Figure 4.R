##########################################
##
## Tia 2024
##
############################################
library(dplyr)
tia = read.csv("C:/Users/esherrar/Documents/RProjects/lsm/data/Tia2024/Data_larvae_Tia_2024.06.19.csv")

tia$date_form = lubridate::dmy(tia$Date)
tia$year = lubridate::year(tia$date_form)
tia$month = lubridate::month(tia$date_form)
tia$time_series = ifelse(tia$year == "2019",tia$month,tia$month+12)
tia$arm = ifelse(tia$Study.arm == "LLIN", 0, 1)
tia$larvalate = tia$L3.Anopheles.spp. + tia$L4.Anopheles.spp.
tia$site_score = ifelse(tia$Locality == "Kakologo","b",
                        ifelse(tia$Locality == "Nambatiourkaha", "c", "a"))


checks1 <- tia %>% group_by(Study.arm,time_series,Locality) %>%
  summarise(Learly = sum(c(L1.Anopheles.spp., L2.Anopheles.spp.)),
            Llate = sum(c(L3.Anopheles.spp., L4.Anopheles.spp.)),
            tot = sum(Total.anopheles),
            meanearly = mean(c(dAn..L1L2)),
            meanlate = mean(c(dAn..L3L4)),## these data for Fig 2A
            meanall = mean(c(dTotal.anopheles))
  )

checks1$time = rep(seq(0,22,2),2)

## This is Figure 2A from the paper

plot(checks1$meanlate[checks1$Study.arm == "LLIN"] ~ checks1$time[checks1$Study.arm == "LLIN"],
     lwd = 2, col = "aquamarine2",type = "l",ylim = c(0,3),
     ylab = "mean Anopheles #", 
     xlab = "Time in months since Bti trial starts")
lines(checks1$meanlate[checks1$Study.arm != "LLIN"] ~ checks1$time[checks1$Study.arm != "LLIN"],
      lwd = 2, col = "darkred")

Tab2Info <- ti2 %>% group_by(Legend,Sampling,Site) %>%
  summarise(peopleHHsTotal = sum(Sleeper),
            gam = sum(An..gambiae),
            nil = sum(An..nili),
            pha = sum(An..pharoensis),
            fun = sum(An..funestus),## these data for Fig 2A
            zei = sum(An..ziemanni),
            totalAnoph = sum(anoph_any)
  )
#############
##
##

ti2 = read.csv("C:/Users/esherrar/Documents/RProjects/lsm/data/Tia2024/Data_Culicidianfauna_Tia_2024.06.19_ellie.csv", header = TRUE)

ti2$date_form = lubridate::dmy(ti2$Date)
ti2$year = lubridate::year(ti2$date_form)
ti2$month = lubridate::month(ti2$date_form)
ti2$time_series = ifelse(ti2$year == "2019",ti2$month,ti2$month+12)
ti2$arm = ifelse(ti2$Legend == "LLIN-only arm", 0, 1)
ti2$Site = ifelse(ti2$Site == "Kol\xe9kaha","other",
                  ti2$Site)
ti2$site_score = ifelse(ti2$Site == "Kakologo","b",
                        ifelse(ti2$Site == "Nambatiourkaha", "c", "a"))
ti2$gam = ti2$An..gambiae

checks2 <- ti2 %>% group_by(Legend,Sampling,Site) %>%
  summarise(peopleHHsTotal = sum(Sleeper),
            gam = sum(An..gambiae),
            nil = sum(An..nili),
            pha = sum(An..pharoensis),
            fun = sum(An..funestus),## these data for Fig 2A
            zei = sum(An..ziemanni),
            totalAnoph = sum(anoph_any)
  )
checks2$perPeopleHHsTotalAnoph = checks2$totalAnoph/checks2$peopleHHsTotal

splitted <- t(sapply(checks2$Sampling, function(x) substring(x, first=c(1,8), last=c(7,9))))
df2 = cbind(checks2, splitted)
df2$time = as.numeric(df2$"...12")
df2 <- df2[order(df2$time),]

df2$perPeopleHHsTotalgam = df2$gam/df2$peopleHHsTotal

plot(df2$perPeopleHHsTotalAnoph[df2$Legend == "LLIN-only arm"] ~ 
       df2$time[df2$Legend == "LLIN-only arm"],
     lwd = 2, col = "aquamarine2",type = "l",ylim = c(0,10),
     ylab = "Mean Anopheles per person per HH", 
     xlab = "Time in months since Bti trial starts")
lines(df2$perPeopleHHsTotalAnoph[df2$Legend != "LLIN-only arm"] ~ 
        df2$time[df2$Legend != "LLIN-only arm"],
      lwd = 2, col = "darkred")

lines(df2$perPeopleHHsTotalgam[df2$Legend == "LLIN-only arm"] ~ 
        df2$time[df2$Legend == "LLIN-only arm"],
      lwd = 1, col = "aquamarine2")
lines(df2$perPeopleHHsTotalgam[df2$Legend != "LLIN-only arm"] ~ 
        df2$time[df2$Legend != "LLIN-only arm"],
      lwd = 1, col = "darkred")

################################################################
##
## Analysis using raw data

## 1 adults
## 2 late stage larvae (Fig 2A)
tia$binary = ifelse(tia$time_series == 3, 0, 1)

modL = data.frame(count_larva = tia$larvalate,
                 count_bites = tia$Total.anopheles,
                 arm = tia$arm,
                 site_score = tia$site_score,
                 month = tia$time_series,
                 binary = tia$binary,
                 sitea = tia$Locality)

library("rstanarm")
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan1 <- rstanarm::stan_glmer.nb(
  count_larva ~ binary * arm + (1|sitea), ## species not important so removed
  data = modL,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 2000 # for speed of example only
)  

print(stan1)
plot(stan1)
plot(stan1, plotfun = "trace")

(trace <- plot(stan1, "trace", pars = "(Intercept)"))
(trace <- plot(stan1, "trace", pars = "binary"))
(trace <- plot(stan1, "trace", pars = "arm"))
(trace <- plot(stan1, "trace", pars = "binary:arm"))

(trace <- plot(stan1, "trace", pars = c("b[(Intercept) sitea:Kakologo]",
                                        "b[(Intercept) sitea:Kolekaha]",
                                        "b[(Intercept) sitea:Lofinekaha]",
                                        "b[(Intercept) sitea:Nambatiourkaha]")) )
(trace <- plot(stan1, "trace", pars = c("reciprocal_dispersion",
                                        "Sigma[sitea:(Intercept),(Intercept)]")) )

(trace <- plot(stan1, "trace", pars = c("reciprocal_dispersion",
                                        "Sigma[sitea:(Intercept),(Intercept)]")) )

plot(stan1)
priors = prior_summary(stan1)
priors$prior$adjusted_scale

# plot(stan1, plotfun = "trace")
# summary(stan1)
# prior_summary(stan1)
# pp_check(stan1, plotfun = "stat", stat = "mean")
# pp_check(stan1, plotfun = "dens_overlay")
# yy_pp = posterior_predict(stan1)
# n_sims = nrow(yy_pp)
# subset = sample(n_sims, 100)
# ppc_dens_overlay(log10(modL$count_larva+1), log10(yy_pp[subset,]+1))
# library(shinystan)

## equation
translate_f = function(pre_or_post,
                       arm_val){
  cnt = (stan1$coefficients[1] + 
              stan1$coefficients[2] * pre_or_post +
              stan1$coefficients[3] * arm_val +
              stan1$coefficients[4] * pre_or_post * arm_val)
  return(as.numeric(cnt))
}
B_pre_con = translate_f(pre_or_post = 0,arm_val = 0)
A_pre_trt = translate_f(pre_or_post = 0,arm_val = 1)

D_post_con = translate_f(pre_or_post = 1,arm_val = 0)
C_post_trt = translate_f(pre_or_post = 1,arm_val = 1)

cc2 = exp((C_post_trt-A_pre_trt)-(D_post_con-B_pre_con))
100*(1-cc2)

ccc = exp(((stan1$coefficients[1] + 
  stan1$coefficients[2] * 1 +
  stan1$coefficients[3] * 1 +
  stan1$coefficients[4] * 1 * 1) - (stan1$coefficients[1] +
                                      stan1$coefficients[3] * 1)) -
  ((stan1$coefficients[1] + 
      stan1$coefficients[2] * 1) - (stan1$coefficients[1])))
1-ccc
## 100*(1-exp(stan1$coefficients[4]))


### Generate uncertainty
translate_uncert_f = function(num_reps){
  B_pre_con = A_pre_trt = 
    D_post_con = C_post_trt = numeric(num_reps)
  
  translate_f = function(pre_or_post,
                         arm_val,
                         which_species,rr){
    cnt = (sims1[rr,1] + 
                sims1[rr,2] * pre_or_post +
                sims1[rr,3] * arm_val +
                sims1[rr,4] * pre_or_post * arm_val )
    
    return(as.numeric(cnt))
  }
  sims1 = as.matrix(stan1)
  n_sims = nrow(sims1)
  random_draw = sample(n_sims,num_reps)
  for(j in 1:num_reps){
    B_pre_con[j] = translate_f(pre_or_post = 0,arm_val = 0, rr=random_draw[j])
    A_pre_trt[j] = translate_f(pre_or_post = 0,arm_val = 1,rr=random_draw[j])
    
    D_post_con[j] = translate_f(pre_or_post = 1,arm_val = 0,rr=random_draw[j])
    C_post_trt[j] = translate_f(pre_or_post = 1,arm_val = 1,rr=random_draw[j])
  }
  ccc = exp((C_post_trt - A_pre_trt) -
              (D_post_con - B_pre_con))
  did = 100*(1-ccc)
  return(data.frame(did = did,
                    B_pre_con = exp(B_pre_con),
                    A_pre_trt = exp(A_pre_trt),
                    D_post_con = exp(D_post_con),
                    C_post_trt = exp(C_post_trt)))
}

rrr = translate_uncert_f(num_reps = 1000)
quantile(rrr[,1],c(0.05,0.5,0.95))

# par(mar = c(3,4,1,1))
# par(mfrow = c(2,1))
boxplot(sqrt(modL$count_larva[modL$binary == 0 & modL$arm == 0]),
        sqrt(modL$count_larva[modL$binary == 0 & modL$arm == 1]),
        NA,
        sqrt(modL$count_larva[modL$binary == 1 & modL$arm == 0]),
        sqrt(modL$count_larva[modL$binary == 1 & modL$arm == 1]), 
        yaxt = "n",border = "grey",
        
        ylim = sqrt(c(0,160)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")

text(1,sqrt(160),"Cote D'Ivoire")
axis(2, las = 2, at = sqrt(seq(0,160,20)),labels = seq(0,160,20))
abline(v=3,lty=1,col="grey34")

mtext("Pre-intervention      Post-intervention",
      cex = 0.8,
      side = 1, line = 1)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,2],0.5))),
         y1 = sqrt(as.numeric(quantile(rrr[,4],0.5))),
         col="darkred")
diff_before = sqrt(as.numeric(quantile(rrr[,3],0.5))) - sqrt(as.numeric(quantile(rrr[,2],0.5)))
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,2],0.5)))+diff_before,
         y1 = sqrt(as.numeric(quantile(rrr[,4],0.5)))+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before) 
segments(x0 = 3,
         x1 = 4.5,
         y0 = sqrt(34),
         y1 = sqrt(as.numeric(quantile(rrr[,5],0.5))),
         lty=2,col="darkblue")

points(sqrt(c(as.numeric(quantile(rrr[,2],0.5)),
              as.numeric(quantile(rrr[,3],0.5)),
              as.numeric(quantile(rrr[,4],0.5)),
              as.numeric(quantile(rrr[,5],0.5)))) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,2],0.05)),
                    as.numeric(quantile(rrr[,3],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,4],0.05)),
                    as.numeric(quantile(rrr[,5],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,2],0.95)),
                    as.numeric(quantile(rrr[,3],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,4],0.95)),
                    as.numeric(quantile(rrr[,5],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}
legend("right",legend = c("Control","LSM"),
       col = c("red","blue"),pch=15,bty="n")

# boxplot(modL$count_larva ~ modL$arm + modL$binary,
#         ylab = "Estimated Anopheles late stage larvae per dip",
#         xlab = "",xaxt="n",
#         ylim=c(0,120),
#         yaxt = "n") ## 0 no LSM, 1 some, 0 before, 1 after
# axis(2, las = 2, seq(0,100,20))
# axis(1, las = 1, at=c(1,2,3,4), labels = c("No LSM", "LSM","No LSM", "LSM"))
# mtext("Before                                   After", side = 1, line = 2.5)
# points(modL$count_larva[modL$arm == 0 & modL$binary == 0] ~ 
#          sample(size = length(modL$count_larva[modL$arm == 0 & modL$binary == 0]),
#                 x = rnorm(n = 20,mean = 1, sd = 0.15),
#                 replace = TRUE),
#        pch=19,col = "grey40")
# points(modL$count_larva[modL$arm == 1 & modL$binary == 0] ~ sample(size = length(modL$count_larva[modL$arm == 1 & modL$binary == 0]),
#                                                                 x = rnorm(n = 20,mean = 2, sd = 0.15),
#                                                                 replace = TRUE),
#        pch=19,col = "darkgreen")
# 
# points(modL$count_larva[modL$arm == 0 & modL$binary == 1] ~ sample(size = length(modL$count_larva[modL$arm == 0 & modL$binary == 1]),
#                                                                 x = rnorm(n = 20,mean = 3, sd = 0.15),
#                                                                 replace = TRUE),
#        pch=19,col = "grey40")
# points(modL$count_larva[modL$arm == 1 & modL$binary == 1] ~ sample(size = length(modL$count_larva[modL$arm == 1 & modL$binary == 1]),
#                                                                 x = rnorm(n = 20,mean = 4, sd = 0.15),
#                                                                 replace = TRUE),
#        pch=19,col = "darkgreen")
# 
# 
# 
# 
# mid_bound1 = c(exp(posterior_interval(stan1, prob = 0.5)[1,1]),
#                exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[3,1]),
#                exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[2,1]),
#                exp(posterior_interval(stan1, prob = 0.5)[1,1] + posterior_interval(stan1, prob = 0.5)[4,1]))
# 
# lower_bound1 = c(exp(posterior_interval(stan1, prob = 0.8)[1,1]),
#                  exp(posterior_interval(stan1, prob = 0.8)[1,1] + posterior_interval(stan1, prob = 0.8)[3,1]),
#                  exp(posterior_interval(stan1, prob = 0.8)[1,1] + posterior_interval(stan1, prob = 0.8)[2,1]),
#                  exp(posterior_interval(stan1, prob = 0.8)[1,1] + posterior_interval(stan1, prob = 0.8)[4,1]))
# 
# upper_bound1 = c(exp(posterior_interval(stan1, prob = 0.8)[1,2]),
#                  exp(posterior_interval(stan1, prob = 0.8)[1,2] + posterior_interval(stan1, prob = 0.8)[3,2]),
#                  exp(posterior_interval(stan1, prob = 0.8)[1,2] + posterior_interval(stan1, prob = 0.8)[2,2]),
#                  exp(posterior_interval(stan1, prob = 0.8)[1,2] + posterior_interval(stan1, prob = 0.8)[4,2]))
# 
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
# 
# points(mid_bound1 ~ segs1,pch=19,cex = 1.2, col = "red")

####################################################
##
## Repeat for the two villages
stan2 <- rstanarm::stan_glm.nb(
  count_larva ~ binary * site_score, ## species not important so removed
  data = modL,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 2000 # for speed of example only
)  

print(stan2)
plot(stan2)
plot(stan2, plotfun = "trace")
(trace <- plot(stan2, "trace", pars = "(Intercept)"))
(trace <- plot(stan2, "trace", pars = "binary"))
(trace <- plot(stan2, "trace", pars = c("site_scoreb",
                                        "site_scorec")) )
(trace <- plot(stan2, "trace", pars = c("binary:site_scoreb",
                                        "binary:site_scorec")) )
(trace <- plot(stan2, "trace", pars = "reciprocal_dispersion"))

summary(stan2)
priors = prior_summary(stan2)
priors$prior$adjusted_scale
pp_check(stan2, plotfun = "stat", stat = "mean")
pp_check(stan2, plotfun = "dens_overlay")
yy_pp = posterior_predict(stan2)
n_sims = nrow(yy_pp)
subset = sample(n_sims, 100)
# rstanarm::ppc_dens_overlay(log10(modL$count_larva+1), log10(yy_pp[subset,]+1))
# library(shinystan)

## equation
translate_f = function(pre_or_post,
                       vill1,vill2){
  cnt = (stan2$coefficients[1] + 
           stan2$coefficients[2] * pre_or_post +
           stan2$coefficients[3] * vill1 +
           stan2$coefficients[4] * vill2 +
           stan2$coefficients[5] * pre_or_post * vill1 +
           stan2$coefficients[6] * pre_or_post * vill2)
  return(as.numeric(cnt))
}
B_pre_con = translate_f(pre_or_post = 0,vill1 = 0,vill2 = 0)
A1_pre_trt = translate_f(pre_or_post = 0,vill1 = 1,vill2 = 0)
A2_pre_trt = translate_f(pre_or_post = 0,vill1 = 0,vill2 = 1)

D_post_con = translate_f(pre_or_post = 1,vill1 = 0,vill2 = 0)
C1_post_trt = translate_f(pre_or_post = 1,vill1 = 1,vill2 = 0)
C2_post_trt = translate_f(pre_or_post = 1,vill1 = 0,vill2 = 1)

ccv1 = exp((C1_post_trt-A1_pre_trt)-(D_post_con-B_pre_con))
100*(1-ccv1)
ccv2 = exp((C2_post_trt-A2_pre_trt)-(D_post_con-B_pre_con))
100*(1-ccv2)

## 100*(1-exp(stan1$coefficients[5]))


### Generate uncertainty
translate_uncert_f = function(num_reps){
  B_pre_con = A_pre_trt = 
    D_post_con = C_post_trt = numeric(num_reps)
  
  translate_f = function(pre_or_post,
                         vill1,vill2,rr){
    cnt = (sims2[rr,1] + 
             sims2[rr,2] * pre_or_post +
             sims2[rr,3] * vill1 +
             sims2[rr,4] * vill2 +
             sims2[rr,5] * pre_or_post * vill1 + 
             sims2[rr,6] * pre_or_post * vill2)
    
    return(as.numeric(cnt))
  }
  sims2 = as.matrix(stan2)
  n_sims = nrow(sims2)
  random_draw = sample(n_sims,num_reps)
  for(j in 1:num_reps){
    B_pre_con[j] = translate_f(pre_or_post = 0,vill1 = 0,vill2 = 0, rr=random_draw[j])
    A1_pre_trt[j] = translate_f(pre_or_post = 0,vill1 = 1,vill2 = 0, rr=random_draw[j])
    A2_pre_trt[j] = translate_f(pre_or_post = 0,vill1 = 0,vill2 = 1, rr=random_draw[j])
    
    D_post_con[j] = translate_f(pre_or_post = 1,vill1 = 0,vill2 = 0, rr=random_draw[j])
    C1_post_trt[j] = translate_f(pre_or_post = 1,vill1 = 1,vill2 = 0, rr=random_draw[j])
    C2_post_trt[j] = translate_f(pre_or_post = 1,vill1 = 0,vill2 = 1, rr=random_draw[j])
  }
  ccc1 = exp((C1_post_trt - A1_pre_trt) -
              (D_post_con - B_pre_con))
  did1 = 100*(1-ccc1)
  ccc2 = exp((C2_post_trt - A2_pre_trt) -
               (D_post_con - B_pre_con))
  did2 = 100*(1-ccc2)
  return(data.frame(did1 = did1,
                    did2 = did2,
                    B_pre_con = exp(B_pre_con),
                    A1_pre_trt = exp(A1_pre_trt),
                    A2_pre_trt = exp(A2_pre_trt),
                    D_post_con = exp(D_post_con),
                    C1_post_trt = exp(C1_post_trt),
                    C2_post_trt = exp(C2_post_trt)))
}

rrr = translate_uncert_f(num_reps = 1000)
quantile(rrr[,1],c(0.05,0.5,0.95))
quantile(rrr[,2],c(0.05,0.5,0.95))
par(mfrow = c(1,2))
par(mar = c(4,4,1,1))
boxplot(sqrt(modL$count_larva[modL$binary == 0 & modL$arm == 0]),
        sqrt(modL$count_larva[modL$binary == 0 & modL$site_score == "b"]),
        NA,
        sqrt(modL$count_larva[modL$binary == 1 & modL$arm == 0]),
        sqrt(modL$count_larva[modL$binary == 1 & modL$site_score == "b"]), 
        yaxt = "n",border = "grey",
        
        ylim = sqrt(c(0,160)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")

text(1,sqrt(160),"Kakologo")
axis(2, las = 2, at = sqrt(seq(0,160,20)),labels = seq(0,160,20))
abline(v=3,lty=1,col="grey34")

mtext("Pre-intervention      Post-intervention",
      cex = 0.8,
      side = 1, line = 1)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5))),
         y1 = sqrt(as.numeric(quantile(rrr[,6],0.5))),
         col="darkred")
diff_before = sqrt(as.numeric(quantile(rrr[,3],0.5))) - sqrt(as.numeric(quantile(rrr[,4],0.5)))
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before,
         y1 = sqrt(as.numeric(quantile(rrr[,6],0.5)))+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before) 
segments(x0 = 3,
         x1 = 4.5,
         y0 = 2,
         y1 = sqrt(as.numeric(quantile(rrr[,7],0.5))),
         lty=2,col="darkblue")

points(sqrt(c(as.numeric(quantile(rrr[,3],0.5)),
              as.numeric(quantile(rrr[,4],0.5)),
              as.numeric(quantile(rrr[,6],0.5)),
              as.numeric(quantile(rrr[,7],0.5)))) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,3],0.05)),
                    as.numeric(quantile(rrr[,4],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,6],0.05)),
                    as.numeric(quantile(rrr[,7],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,3],0.95)),
                    as.numeric(quantile(rrr[,4],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,6],0.95)),
                    as.numeric(quantile(rrr[,7],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}
legend("topright",legend = c("Control","LSM"),
       col = c("red","blue"),pch=15,bty="n")


boxplot(sqrt(modL$count_larva[modL$binary == 0 & modL$arm == 0]),
        sqrt(modL$count_larva[modL$binary == 0 & modL$site_score == "c"]),
        NA,
        sqrt(modL$count_larva[modL$binary == 1 & modL$arm == 0]),
        sqrt(modL$count_larva[modL$binary == 1 & modL$site_score == "c"]), 
        yaxt = "n",border = "grey",
        
        ylim = sqrt(c(0,160)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")

text(1.5,sqrt(160),"Nambatiourkaha")
axis(2, las = 2, at = sqrt(seq(0,160,20)),labels = seq(0,160,20))
abline(v=3,lty=1,col="grey34")

mtext("Pre-intervention      Post-intervention",
      cex = 0.8,
      side = 1, line = 1)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5))),
         y1 = sqrt(as.numeric(quantile(rrr[,6],0.5))),
         col="darkred")
diff_before = sqrt(as.numeric(quantile(rrr[,5],0.5))) - sqrt(as.numeric(quantile(rrr[,3],0.5)))
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before,
         y1 = sqrt(as.numeric(quantile(rrr[,6],0.5)))+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before) 
segments(x0 = 3,
         x1 = 4.5,
         y0 = sqrt(12),
         y1 = sqrt(as.numeric(quantile(rrr[,8],0.5))),
         lty=2,col="darkblue")

points(sqrt(c(as.numeric(quantile(rrr[,3],0.5)),
              as.numeric(quantile(rrr[,5],0.5)),
              as.numeric(quantile(rrr[,6],0.5)),
              as.numeric(quantile(rrr[,8],0.5)))) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,3],0.05)),
                    as.numeric(quantile(rrr[,5],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,6],0.05)),
                    as.numeric(quantile(rrr[,8],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,3],0.95)),
                    as.numeric(quantile(rrr[,5],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,6],0.95)),
                    as.numeric(quantile(rrr[,8],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}
par(xpd=NA,cex = 1.11)

text(x = -8, y = 13.5,"(A)",cex=0.8)
text(x = -0.7, y = 13.5,"(B)",cex=0.8)
#################################
##
## And estimate the same for the adult biting
##
##################################
ti2 = read.csv("C:/Users/esherrar/Documents/RProjects/lsm/data/Tia2024/Data_Culicidianfauna_Tia_2024.06.19_ellie.csv", header = TRUE)

ti2$date_form = lubridate::dmy(ti2$Date)
ti2$year = lubridate::year(ti2$date_form)
ti2$month = lubridate::month(ti2$date_form)
ti2$time_series = ifelse(ti2$year == "2019",ti2$month,ti2$month+12)
ti2$arm = ifelse(ti2$Legend == "LLIN-only arm", 0, 1)
ti2$Site = ifelse(ti2$Site == "Kol\xe9kaha","other",
                  ti2$Site)
ti2$site_score = ifelse(ti2$Site == "Kakologo","b",
                        ifelse(ti2$Site == "Nambatiourkaha", "c", "a"))
ti2$gam = ti2$An..gambiae

ti2$binary = ifelse(ti2$time_series == 3, 0, 1)

modL2 = data.frame(count_bites = ti2$anoph_any,
                  arm = ti2$arm,
                  month = ti2$time_series,
                  binary = ti2$binary,
                  peopeHH = ti2$Sleeper,
                  sitea = ti2$Site,
                  site_score = ti2$site_score)

stan1b <- rstanarm::stan_glmer.nb(
  count_bites ~ binary * arm + (1|peopeHH) + (1|sitea), ## species not important so removed
  data = modL2,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 1000 # for speed of example only
)  

print(stan1b)
plot(stan1b)
# plot(stan1b, plotfun = "trace")
summary(stan1b)
prior_summary(stan1b)
pp_check(stan1b, plotfun = "stat", stat = "mean")
pp_check(stan1b, plotfun = "dens_overlay")
yy_pp = posterior_predict(stan1b)
n_sims = nrow(yy_pp)
subset = sample(n_sims, 100)
ppc_dens_overlay(log10(modL2$count_bites+1), log10(yy_pp[subset,]+1))
# library(shinystan)

## equation
translate_f = function(pre_or_post,
                       arm_val){
  cnt = (stan1b$coefficients[1] + 
           stan1b$coefficients[2] * pre_or_post +
           stan1b$coefficients[3] * arm_val +
           stan1b$coefficients[4] * pre_or_post * arm_val)
  return(as.numeric(cnt))
}
B_pre_con = translate_f(pre_or_post = 0,arm_val = 0)
A_pre_trt = translate_f(pre_or_post = 0,arm_val = 1)

D_post_con = translate_f(pre_or_post = 1,arm_val = 0)
C_post_trt = translate_f(pre_or_post = 1,arm_val = 1)

cc2 = exp((C_post_trt-A_pre_trt)-(D_post_con-B_pre_con))
100*(1-cc2)

ccc = exp(((stan1b$coefficients[1] + 
              stan1b$coefficients[2] * 1 +
              stan1b$coefficients[3] * 1 +
              stan1b$coefficients[4] * 1 * 1) - (stan1b$coefficients[1] +
                                                  stan1b$coefficients[3] * 1)) -
            ((stan1b$coefficients[1] + 
                stan1b$coefficients[2] * 1) - (stan1b$coefficients[1])))
1-ccc
## 100*(1-exp(stan1$coefficients[4]))


### Generate uncertainty
translate_uncert_f = function(num_reps){
  B_pre_con = A_pre_trt = 
    D_post_con = C_post_trt = numeric(num_reps)
  
  translate_f = function(pre_or_post,
                         arm_val,rr){
    cnt = (sims1b[rr,1] + 
             sims1b[rr,2] * pre_or_post +
             sims1b[rr,3] * arm_val +
             sims1b[rr,4] * pre_or_post * arm_val )
    
    return(as.numeric(cnt))
  }
  sims1b = as.matrix(stan1b)
  n_sims = nrow(sims1b)
  random_draw = sample(n_sims,num_reps)
  for(j in 1:num_reps){
    B_pre_con[j] = translate_f(pre_or_post = 0,arm_val = 0, rr=random_draw[j])
    A_pre_trt[j] = translate_f(pre_or_post = 0,arm_val = 1,rr=random_draw[j])
    
    D_post_con[j] = translate_f(pre_or_post = 1,arm_val = 0,rr=random_draw[j])
    C_post_trt[j] = translate_f(pre_or_post = 1,arm_val = 1,rr=random_draw[j])
  }
  ccc = exp((C_post_trt - A_pre_trt) -
              (D_post_con - B_pre_con))
  did = 100*(1-ccc)
  return(data.frame(did = did,
                    B_pre_con = exp(B_pre_con),
                    A_pre_trt = exp(A_pre_trt),
                    D_post_con = exp(D_post_con),
                    C_post_trt = exp(C_post_trt)))
}

rrr = translate_uncert_f(num_reps = 1000)
quantile(rrr[,1],c(0.05,0.5,0.95))


boxplot(modL2$count_bites ~ modL2$arm + modL2$binary,
        ylab = "Anopheles mosquitoes adults estimated counts",
        xlab = "",xaxt="n",
        ylim=c(0,120),
        yaxt = "n") ## 0 no LSM, 1 some, 0 before, 1 after
axis(2, las = 2, seq(0,100,20))
axis(1, las = 1, at=c(1,2,3,4), labels = c("No LSM", "LSM","No LSM", "LSM"))
mtext("Before                                   After", side = 1, line = 2.5)
points(modL2$count_bites[modL2$arm == 0 & modL2$binary == 0] ~ sample(size = length(modL2$count_bites[modL2$arm == 0 & modL2$binary == 0]),
                                                                   x = rnorm(n = 20,mean = 1, sd = 0.15),
                                                                   replace = TRUE),
       pch=19,col = "grey40")
points(modL2$count_bites[modL2$arm == 1 & modL2$binary == 0] ~ sample(size = length(modL2$count_bites[modL2$arm == 1 & modL2$binary == 0]),
                                                                   x = rnorm(n = 20,mean = 2, sd = 0.15),
                                                                   replace = TRUE),
       pch=19,col = "darkgreen")

points(modL2$count_bites[modL2$arm == 0 & modL2$binary == 1] ~ sample(size = length(modL2$count_bites[modL2$arm == 0 & modL2$binary == 1]),
                                                                   x = rnorm(n = 20,mean = 3, sd = 0.15),
                                                                   replace = TRUE),
       pch=19,col = "grey40")
points(modL2$count_bites[modL2$arm == 1 & modL2$binary == 1] ~ sample(size = length(modL2$count_bites[modL2$arm == 1 & modL2$binary == 1]),
                                                                   x = rnorm(n = 20,mean = 4, sd = 0.15),
                                                                   replace = TRUE),
       pch=19,col = "darkgreen")




mid_bound1 = c(exp(posterior_interval(stan2, prob = 0.5)[1,1]),
               exp(posterior_interval(stan2, prob = 0.5)[1,1] + posterior_interval(stan2, prob = 0.5)[3,1]),
               exp(posterior_interval(stan2, prob = 0.5)[1,1] + posterior_interval(stan2, prob = 0.5)[2,1]),
               exp(posterior_interval(stan2, prob = 0.5)[1,1] + posterior_interval(stan2, prob = 0.5)[4,1]))

lower_bound1 = c(exp(posterior_interval(stan2, prob = 0.8)[1,1]),
                 exp(posterior_interval(stan2, prob = 0.8)[1,1] + posterior_interval(stan2, prob = 0.8)[3,1]),
                 exp(posterior_interval(stan2, prob = 0.8)[1,1] + posterior_interval(stan2, prob = 0.8)[2,1]),
                 exp(posterior_interval(stan2, prob = 0.8)[1,1] + posterior_interval(stan2, prob = 0.8)[4,1]))

upper_bound1 = c(exp(posterior_interval(stan2, prob = 0.8)[1,2]),
                 exp(posterior_interval(stan2, prob = 0.8)[1,2] + posterior_interval(stan2, prob = 0.8)[3,2]),
                 exp(posterior_interval(stan2, prob = 0.8)[1,2] + posterior_interval(stan2, prob = 0.8)[2,2]),
                 exp(posterior_interval(stan2, prob = 0.8)[1,2] + posterior_interval(stan2, prob = 0.8)[4,2]))

segs1 = c(0.96,1.96,2.96,3.96)
for(i in 1:4){
  
  segments(x0 = segs1[i],
           x1 = segs1[i],
           y0 = lower_bound1[i],
           y1 = upper_bound1[i],
           lwd=4,col="red")
  
}

points(mid_bound1 ~ segs1,pch=19,cex = 1.2, col = "red")

#################
##
## Now the village level differences
stan2b <- rstanarm::stan_glm.nb(
  count_bites ~ binary * site_score, ## species not important so removed
  data = modL2,
  link = "log",
  # prior_aux = exponential(1.5, autoscale=TRUE),
  chains = 4, iter = 2000 # for speed of example only
)  

print(stan2b)
plot(stan2b)
plot(stan2b, plotfun = "trace")
summary(stan2b)
prior_summary(stan2b)
pp_check(stan2b, plotfun = "stat", stat = "mean")
pp_check(stan2b, plotfun = "dens_overlay")
yy_pp = posterior_predict(stan2b)
n_sims = nrow(yy_pp)
subset = sample(n_sims, 100)
# ppc_dens_overlay(log10(modL2$count_bites+1), log10(yy_pp[subset,]+1))
# library(shinystan)

## equation
translate_f = function(pre_or_post,
                       vill1,vill2){
  cnt = (stan2b$coefficients[1] + 
           stan2b$coefficients[2] * pre_or_post +
           stan2b$coefficients[3] * vill1 +
           stan2b$coefficients[4] * vill2 +
           stan2b$coefficients[5] * pre_or_post * vill1 +
           stan2b$coefficients[6] * pre_or_post * vill2)
  return(as.numeric(cnt))
}
B_pre_con = translate_f(pre_or_post = 0,vill1 = 0,vill2 = 0)
A1_pre_trt = translate_f(pre_or_post = 0,vill1 = 1,vill2 = 0)
A2_pre_trt = translate_f(pre_or_post = 0,vill1 = 0,vill2 = 1)

D_post_con = translate_f(pre_or_post = 1,vill1 = 0,vill2 = 0)
C1_post_trt = translate_f(pre_or_post = 1,vill1 = 1,vill2 = 0)
C2_post_trt = translate_f(pre_or_post = 1,vill1 = 0,vill2 = 1)

ccv1 = exp((C1_post_trt-A1_pre_trt)-(D_post_con-B_pre_con))
100*(1-ccv1)
ccv2 = exp((C2_post_trt-A2_pre_trt)-(D_post_con-B_pre_con))
100*(1-ccv2)

## 100*(1-exp(stan1$coefficients[5]))


### Generate uncertainty
translate_uncert_f = function(num_reps){
  B_pre_con = A_pre_trt = 
    D_post_con = C_post_trt = numeric(num_reps)
  
  translate_f = function(pre_or_post,
                         vill1,vill2,rr){
    cnt = (sims2[rr,1] + 
             sims2[rr,2] * pre_or_post +
             sims2[rr,3] * vill1 +
             sims2[rr,4] * vill2 +
             sims2[rr,5] * pre_or_post * vill1 + 
             sims2[rr,6] * pre_or_post * vill2)
    
    return(as.numeric(cnt))
  }
  sims2 = as.matrix(stan2b)
  n_sims = nrow(sims2)
  random_draw = sample(n_sims,num_reps)
  for(j in 1:num_reps){
    B_pre_con[j] = translate_f(pre_or_post = 0,vill1 = 0,vill2 = 0, rr=random_draw[j])
    A1_pre_trt[j] = translate_f(pre_or_post = 0,vill1 = 1,vill2 = 0, rr=random_draw[j])
    A2_pre_trt[j] = translate_f(pre_or_post = 0,vill1 = 0,vill2 = 1, rr=random_draw[j])
    
    D_post_con[j] = translate_f(pre_or_post = 1,vill1 = 0,vill2 = 0, rr=random_draw[j])
    C1_post_trt[j] = translate_f(pre_or_post = 1,vill1 = 1,vill2 = 0, rr=random_draw[j])
    C2_post_trt[j] = translate_f(pre_or_post = 1,vill1 = 0,vill2 = 1, rr=random_draw[j])
  }
  ccc1 = exp((C1_post_trt - A1_pre_trt) -
               (D_post_con - B_pre_con))
  did1 = 100*(1-ccc1)
  ccc2 = exp((C2_post_trt - A2_pre_trt) -
               (D_post_con - B_pre_con))
  did2 = 100*(1-ccc2)
  return(data.frame(did1 = did1,
                    did2 = did2,
                    B_pre_con = exp(B_pre_con),
                    A1_pre_trt = exp(A1_pre_trt),
                    A2_pre_trt = exp(A2_pre_trt),
                    D_post_con = exp(D_post_con),
                    C1_post_trt = exp(C1_post_trt),
                    C2_post_trt = exp(C2_post_trt)))
}

rrr = translate_uncert_f(num_reps = 1000)
quantile(rrr[,1],c(0.05,0.5,0.95))
quantile(rrr[,2],c(0.05,0.5,0.95))


par(mfrow = c(1,2))
par(mar = c(4,4,1,1))
boxplot(sqrt(modL2$count_bites[modL2$binary == 0 & modL2$arm == 0]),
        sqrt(modL2$count_bites[modL2$binary == 0 & modL2$site_score == "b"]),
        NA,
        sqrt(modL2$count_bites[modL2$binary == 1 & modL2$arm == 0]),
        sqrt(modL2$count_bites[modL2$binary == 1 & modL2$site_score == "b"]), 
        yaxt = "n",border = "grey",
        
        ylim = sqrt(c(0,160)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")

text(1,sqrt(160),"Kakologo")
axis(2, las = 2, at = sqrt(seq(0,160,20)),labels = seq(0,160,20))
abline(v=3,lty=1,col="grey34")

mtext("Pre-intervention      Post-intervention",
      cex = 0.8,
      side = 1, line = 1)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5))),
         y1 = sqrt(as.numeric(quantile(rrr[,6],0.5))),
         col="darkred")
diff_before = sqrt(as.numeric(quantile(rrr[,3],0.5))) - sqrt(as.numeric(quantile(rrr[,4],0.5)))
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before,
         y1 = sqrt(as.numeric(quantile(rrr[,6],0.5)))+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before) 
segments(x0 = 3,
         x1 = 4.5,
         y0 = 2.9,
         y1 = sqrt(as.numeric(quantile(rrr[,7],0.5))),
         lty=2,col="darkblue")

points(sqrt(c(as.numeric(quantile(rrr[,3],0.5)),
              as.numeric(quantile(rrr[,4],0.5)),
              as.numeric(quantile(rrr[,6],0.5)),
              as.numeric(quantile(rrr[,7],0.5)))) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,3],0.05)),
                    as.numeric(quantile(rrr[,4],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,6],0.05)),
                    as.numeric(quantile(rrr[,7],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,3],0.95)),
                    as.numeric(quantile(rrr[,4],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,6],0.95)),
                    as.numeric(quantile(rrr[,7],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}
legend("topright",legend = c("Control","LSM"),
       col = c("red","blue"),pch=15,bty="n")


boxplot(sqrt(modL2$count_bites[modL2$binary == 0 & modL2$arm == 0]),
        sqrt(modL2$count_bites[modL2$binary == 0 & modL2$site_score == "c"]),
        NA,
        sqrt(modL2$count_bites[modL2$binary == 1 & modL2$arm == 0]),
        sqrt(modL2$count_bites[modL2$binary == 1 & modL2$site_score == "c"]), 
        yaxt = "n",border = "grey",
        
        ylim = sqrt(c(0,160)),
        col = adegenet::transp(c("red","blue",NA,"red","blue"),0.4),
        ylab = "Mosquito counts",xlab = "Treatment")

text(1.5,sqrt(160),"Nambatiourkaha")
axis(2, las = 2, at = sqrt(seq(0,160,20)),labels = seq(0,160,20))
abline(v=3,lty=1,col="grey34")

mtext("Pre-intervention      Post-intervention",
      cex = 0.8,
      side = 1, line = 1)
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5))),
         y1 = sqrt(as.numeric(quantile(rrr[,6],0.5))),
         col="darkred")
diff_before = sqrt(as.numeric(quantile(rrr[,5],0.5))) - sqrt(as.numeric(quantile(rrr[,3],0.5)))
segments(x0 = 1.5,
         x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before,
         y1 = sqrt(as.numeric(quantile(rrr[,6],0.5)))+diff_before,
         lty=1,col="darkblue")
y2 = (3*sqrt(as.numeric(quantile(rrr[,3],0.5)))+diff_before) 
segments(x0 = 3,
         x1 = 4.5,
         y0 = 2.6,
         y1 = sqrt(as.numeric(quantile(rrr[,8],0.5))),
         lty=2,col="darkblue")

points(sqrt(c(as.numeric(quantile(rrr[,3],0.5)),
              as.numeric(quantile(rrr[,5],0.5)),
              as.numeric(quantile(rrr[,6],0.5)),
              as.numeric(quantile(rrr[,8],0.5)))) ~ c(1,2,4,5),col = c("red","blue","red","blue"),pch=19) 

#
mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,3],0.05)),
                    as.numeric(quantile(rrr[,5],0.05)),
                    NA,
                    as.numeric(quantile(rrr[,6],0.05)),
                    as.numeric(quantile(rrr[,8],0.05))) )#
maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,3],0.95)),
                    as.numeric(quantile(rrr[,5],0.95)),
                    NA,
                    as.numeric(quantile(rrr[,6],0.95)),
                    as.numeric(quantile(rrr[,8],0.95))) )
col_segs = c("red","blue",NA,"red","blue")
for(i in c(1,2,4,5)){
  segments(x0 = i,
           x1 = i,
           y0 = mins_stan4[i],
           y1 = maxs_stan4[i],col = col_segs[i],lwd=2)       
  
}
par(xpd=NA,cex = 1.11)

text(x = -8, y = 13.5,"(A)",cex=0.8)
text(x = -0.7, y = 13.5,"(B)",cex=0.8)

###
###
write.csv(modL,"C:/Users/esherrar/Documents/RProjects/lsm/data/Tia2024/larvalSummaryData.csv")#
write.csv(modL2,"C:/Users/esherrar/Documents/RProjects/lsm/data/Tia2024/adultSummaryData.csv")
