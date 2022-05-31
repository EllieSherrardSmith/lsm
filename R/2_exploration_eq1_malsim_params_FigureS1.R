#################
##
##
par(mfrow=c(2,2))
#################
phi_lsm = 0
rho_alpha = seq(5,-6,length=15)
rho_beta = 0.1
time = 1:365


roh_lsm = array(data=NA, dim=c(365,15))
for(i in 1:15){
  roh_lsm[,i] = phi_lsm + (((1-phi_lsm) * (1 + exp(-rho_alpha[i])))/
                             (1 + exp(-rho_alpha[i]) * exp(rho_beta * time)))
}

## plot the param with an increasing estimate for rho_alpha or rho_beta
cols = adegenet::transp("darkred",seq(0.1,1,length=15))

plot(roh_lsm[1:100,1] ~ time[1:100], ylim=c(0,1),pch="",
     ylab = "rho_alpha (ranging -6 to +5)",
     xlab = "Time in days",yaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
for(i in 1:15){
  lines(roh_lsm[1:100,i] ~ time[1:100], col = cols[i])
}

#################
##
##
#################
phi_lsm = 0
rho_alpha = -4
rho_beta = seq(0.01,0.2,length=15)
time = 1:365


roh_lsm2 = array(data=NA, dim=c(365,15))
for(i in 1:15){
  roh_lsm2[,i] = phi_lsm + (((1-phi_lsm) * (1 + exp(-rho_alpha)))/
                             (1 + exp(-rho_alpha) * exp(rho_beta[i] * time)))
  }

## plot the param with an increasing estimate for rho_alpha or rho_beta

plot(roh_lsm2[1:100,1] ~ time[1:100], ylim=c(0,1),pch="",
     ylab = "rho_beta (ranging 0.01 - 0.2)",
     xlab = "Time in days",yaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
for(i in 1:15){
  lines(roh_lsm2[1:100,i] ~ time[1:100], col = cols[i])
}


#################
##
##
#################
phi_lsm = seq(0,1,length=15)
rho_alpha = -4
rho_beta = 0.1
time = 1:365


roh_lsm3 = array(data=NA, dim=c(365,15))
for(i in 1:15){
  roh_lsm3[,i] = phi_lsm[i] + (((1-phi_lsm[i]) * (1 + exp(-rho_alpha)))/
                              (1 + exp(-rho_alpha) * exp(rho_beta * time)))
}

## plot the param with an increasing estimate for rho_alpha or rho_beta

plot(roh_lsm3[1:100,1] ~ time[1:100], ylim=c(0,1),pch="",
     ylab = "phi_lsm (ranging 0 to 1)",
     xlab = "Time in days",yaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
for(i in 1:15){
  lines(roh_lsm3[1:100,i] ~ time[1:100], col = cols[i])
}


#################
# 0 0.67	
# 1 1	
# 0.22 0.45	
# 0.88 0.11	
# 1 1	
# 1 1

phi_lsm = c(0, 0.22, 0.88,
            0.67,0.45,0.11)
rho_alpha = -4
rho_beta = 0.057
time = 1:365


roh_lsm3 = array(data=NA, dim=c(365,6))
for(i in 1:6){
  roh_lsm3[,i] = phi_lsm[i] + (((1-phi_lsm[i]) * (1 + exp(-rho_alpha)))/
                                 (1 + exp(-rho_alpha) * exp(rho_beta * time)))
}

## plot the param with an increasing estimate for rho_alpha or rho_beta

plot(roh_lsm3[1:100,1] ~ time[1:100], ylim=c(0,1),pch="",
     ylab = "phi_lsm (reflecting Table 1 estimates)",
     xlab = "Time in days",yaxt="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
for(i in 1:3){
  lines(roh_lsm3[1:100,i] ~ time[1:100], col = "red",lty=i)
  lines(roh_lsm3[1:100,i+3] ~ time[1:100], col = "blue",lty=i)
}




