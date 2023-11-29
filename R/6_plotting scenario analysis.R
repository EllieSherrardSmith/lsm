##########################################
##
## We will want to draw a matrix of the cases averted / reduction in prev due to lsm
## need to compare like-with-like, the baseline scenario will be no lsm with matched
## itn_cov, resistance, phiB/phiI, eir_target

set_up = read.csv("R/malariasimulation/scenario_simulations/inputs_set_up.csv",header=TRUE)
## X.1 is the unique ID for the sim runs

# Align the set_up row with the model results
# Calculate the total clinical cases over 1, 2 and 3 years given the scenario counterfactual
## choose the set up for the comparison
base_set_up = subset(set_up, set_up$cc_impact == 1 & set_up$phiB == 0.9 &
                       set_up$itn_cov == 0.8 & set_up$eir_target == 0.08)
## identify the code
vec_set = base_set_up$ID 

# set_up[1:378,] These are the counterfactual simulations
# par(mfrow=c(2,2))
# for(i in c(1)){
# 
#   test = read.csv(paste0("outputs/sim_",vec_set[i],"_run.csv"),header=TRUE)
#   plot(test$pv_182.5_1825 ~ test$timestep,xlim=c(2300,3650),ylim=c(0,1),pch="")
#   lines(test$pv_182.5_1825 ~ test$timestep)
# 
#   for(j in c(1,4,12,16)){
#     xx = j
#     treat = read.csv(paste0("outputs/sim_",vec_set[i]+378*j,"_run.csv"),header=TRUE)
#     lines(treat$pv_182.5_1825 ~ treat$timestep,col="blue")
#     
#   }
# }

# #Subset months, we only use -12 to 48
# library(dplyr)
# test2 = test %>% slice(which(row_number() %% 30 == 1))
# treat2 = treat %>% slice(which(row_number() %% 30 == 1))
# plot(test2$pv_182.5_1825 ~ test2$timestep,xlim=c(2300,3650),ylim=c(0,1),pch="")
# lines(test2$pv_182.5_1825 ~ test2$timestep,lty=1)
# lines(treat2$pv_182.5_1825 ~ treat2$timestep,col="red")
# abline(v=2555)

pull_data_f = function(xx){
  base = expand.grid(sims = c(1:21))
  trt = expand.grid(sims = c(1:21))
  xx = xx ## this is a marker for the 21 sets of data for cc levels (level 1, xx = 0, is the baseline)
  for(i in 1:21){
    
    
    test = read.csv(paste0("outputs/sim_",vec_set[i],"_run.csv"),header=TRUE)
    ## lsm switched on 7*365 run length is 10 years so we have up to 3 follow up years
    
    # the number of detectable cases through time using the `n_detect_[n]_[nn]`
    # the number of clinical cases through time using the   `n_inc_clinical_[n]_[nn]`
    # the number of severe cases through time using the     `n_inc_severe_[n]_[nn]`
    # the prevalence pv_[n]_[nn]
    
    base$detect_U5_y1[i] = sum(test$n_detect_182.5_1825[2555:2920])  ## year 1 
    base$detect_U5_y2[i] = sum(test$n_detect_182.5_1825[2921:3285])  ## year 2
    base$detect_U5_y3[i] = sum(test$n_detect_182.5_1825[3286:3650])  ## year 3
    base$detect_U5_over3years[i] = sum(test$n_detect_182.5_1825[2555:3650])  ## years 1 to 3
    
    base$clin_inc_U5_y1[i] = sum(test$n_inc_clinical_182.5_1825[2555:2920])  ## year 1 
    base$clin_inc_U5_y2[i] = sum(test$n_inc_clinical_182.5_1825[2921:3285])  ## year 2
    base$clin_inc_U5_y3[i] = sum(test$n_inc_clinical_182.5_1825[3286:3650])  ## year 3
    base$clin_inc_U5_over3years[i] = sum(test$n_inc_clinical_182.5_1825[2555:3650])  ## years 1 to 3
    
    base$prev_U5_y1[i] = mean(test$pv_182.5_1825[2555:2920])
    base$prev_U5_y2[i] = mean(test$pv_182.5_1825[2921:3285])
    base$prev_U5_y3[i] = mean(test$pv_182.5_1825[3286:3650])
    base$prev_U5_over3years[i] = mean(test$pv_182.5_1825[2555:3650])
    
    base$detect_all_y1[i] = sum(test$n_detect_0_36500[2555:2920])  ## year 1 
    base$detect_all_y2[i] = sum(test$n_detect_0_36500[2921:3285])  ## year 2
    base$detect_all_y3[i] = sum(test$n_detect_0_36500[3286:3650])  ## year 3
    base$detect_all_over3years[i] = sum(test$n_detect_0_36500[2555:3650])  ## years 1 to 3
    
    base$clin_inc_all_y1[i] = sum(test$n_inc_clinical_0_36500[2555:2920])  ## year 1 
    base$clin_inc_all_y2[i] = sum(test$n_inc_clinical_0_36500[2921:3285])  ## year 2
    base$clin_inc_all_y3[i] = sum(test$n_inc_clinical_0_36500[3286:3650])  ## year 3
    base$clin_inc_all_over3years[i] = sum(test$n_inc_clinical_0_36500[2555:3650])  ## years 1 to 3
    
    base$prev_all_y1[i] = mean(test$pv_0_36500[2555:2920])
    base$prev_all_y2[i] = mean(test$pv_0_36500[2921:3285])
    base$prev_all_y3[i] = mean(test$pv_0_36500[3286:3650])
    base$prev_all_over3years[i] = mean(test$pv_0_36500[2555:3650])
    
    # Calc
    # Calculate the total clinical cases over 1, 2 and 3 years given the scenario with lsm for n levels
    
    treat = read.csv(paste0("outputs/sim_",vec_set[i]+378*xx,"_run.csv"),header=TRUE)
    
    
    trt$detect_U5_y1[i] = sum(treat$n_detect_182.5_1825[2555:2920])  ## year 1 
    trt$detect_U5_y2[i] = sum(treat$n_detect_182.5_1825[2921:3285])  ## year 2
    trt$detect_U5_y3[i] = sum(treat$n_detect_182.5_1825[3286:3650])  ## year 3
    trt$detect_U5_over3years[i] = sum(treat$n_detect_182.5_1825[2555:3650])  ## years 1 to 3
    
    trt$clin_inc_U5_y1[i] = sum(treat$n_inc_clinical_182.5_1825[2555:2920])  ## year 1 
    trt$clin_inc_U5_y2[i] = sum(treat$n_inc_clinical_182.5_1825[2921:3285])  ## year 2
    trt$clin_inc_U5_y3[i] = sum(treat$n_inc_clinical_182.5_1825[3286:3650])  ## year 3
    trt$clin_inc_U5_over3years[i] = sum(treat$n_inc_clinical_182.5_1825[2555:3650])  ## years 1 to 3
    
    trt$prev_U5_y1[i] = mean(treat$pv_182.5_1825[2555:2920])
    trt$prev_U5_y2[i] = mean(treat$pv_182.5_1825[2921:3285])
    trt$prev_U5_y3[i] = mean(treat$pv_182.5_1825[3286:3650])
    trt$prev_U5_over3years[i] = mean(treat$pv_182.5_1825[2555:3650])
    
    trt$detect_all_y1[i] = sum(treat$n_detect_0_36500[2555:2920])  ## year 1 
    trt$detect_all_y2[i] = sum(treat$n_detect_0_36500[2921:3285])  ## year 2
    trt$detect_all_y3[i] = sum(treat$n_detect_0_36500[3286:3650])  ## year 3
    trt$detect_all_over3years[i] = sum(treat$n_detect_0_36500[2555:3650])  ## years 1 to 3
    
    trt$clin_inc_all_y1[i] = sum(treat$n_inc_clinical_0_36500[2555:2920])  ## year 1 
    trt$clin_inc_all_y2[i] = sum(treat$n_inc_clinical_0_36500[2921:3285])  ## year 2
    trt$clin_inc_all_y3[i] = sum(treat$n_inc_clinical_0_36500[3286:3650])  ## year 3
    trt$clin_inc_all_over3years[i] = sum(treat$n_inc_clinical_0_36500[2555:3650])  ## years 1 to 3
    
    trt$prev_all_y1[i] = mean(treat$pv_0_36500[2555:2920])
    trt$prev_all_y2[i] = mean(treat$pv_0_36500[2921:3285])
    trt$prev_all_y3[i] = mean(treat$pv_0_36500[3286:3650])
    trt$prev_all_over3years[i] = mean(treat$pv_0_36500[2555:3650])
  }
  
  ## first set just look at 5% LSM...
  
  # Calculate the absolute and relative difference in these estimates
  
  abs_diff = base - trt ## each
  rel_diff = (base - trt)/base
  
  return(list(base,
              trt,
              abs_diff,
              rel_diff))
}


cc0 = pull_data_f(xx = 0)
cc1 = pull_data_f(xx = 1)
cc2 = pull_data_f(xx = 2)
cc3 = pull_data_f(xx = 3)
cc4 = pull_data_f(xx = 4)
cc5 = pull_data_f(xx = 5)
cc6 = pull_data_f(xx = 6)
cc7 = pull_data_f(xx = 7)
cc8 = pull_data_f(xx = 8)
cc9 = pull_data_f(xx = 9)
cc10 = pull_data_f(xx = 10)
cc11 = pull_data_f(xx = 11)
cc12 = pull_data_f(xx = 12)
cc13 = pull_data_f(xx = 13)
cc14 = pull_data_f(xx = 14)
cc15 = pull_data_f(xx = 15)
cc16 = pull_data_f(xx = 16)
cc17 = pull_data_f(xx = 17)
cc18 = pull_data_f(xx = 18)
cc19 = pull_data_f(xx = 19)
cc20 = pull_data_f(xx = 20)


## plot the data as matrix of level of resistance versus lsm effect for: 
## each level of eir (low, med, high) with 40% ITN
## each level of eir with 60$ ITN
## each level of eir with 80% ITN
## for the seasonal and perennial settings


## Relative difference in xx
vec_min_selections = c(1,3,7,11,14,20,21)
zf1 = c(cc0[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc1[[4]]$clin_inc_all_over3years,
        # cc2[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc3[[4]]$clin_inc_all_over3years,
        cc4[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc5[[4]]$clin_inc_all_over3years,
        # cc6[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc7[[4]]$clin_inc_all_over3years,
        cc8[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc9[[4]]$clin_inc_all_over3years,
        # cc10[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc11[[4]]$clin_inc_all_over3years,
        cc12[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc13[[4]]$clin_inc_all_over3years,
        # cc14[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc15[[4]]$clin_inc_all_over3years,
        cc16[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc17[[4]]$clin_inc_all_over3years,
        # cc18[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc19[[4]]$clin_inc_all_over3years,
        cc20[[4]]$clin_inc_all_over3years[vec_min_selections])

zf2 = c(cc0[[3]]$clin_inc_all_over3years[vec_min_selections],
        # cc1[[4]]$clin_inc_all_over3years,
        # cc2[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc3[[4]]$clin_inc_all_over3years,
        cc4[[3]]$clin_inc_all_over3years[vec_min_selections],
        # cc5[[4]]$clin_inc_all_over3years,
        # cc6[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc7[[4]]$clin_inc_all_over3years,
        cc8[[3]]$clin_inc_all_over3years[vec_min_selections],
        # cc9[[4]]$clin_inc_all_over3years,
        # cc10[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc11[[4]]$clin_inc_all_over3years,
        cc12[[3]]$clin_inc_all_over3years[vec_min_selections],
        # cc13[[4]]$clin_inc_all_over3years,
        # cc14[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc15[[4]]$clin_inc_all_over3years,
        cc16[[3]]$clin_inc_all_over3years[vec_min_selections],
        # cc17[[4]]$clin_inc_all_over3years,
        # cc18[[4]]$clin_inc_all_over3years[vec_min_selections],
        # cc19[[4]]$clin_inc_all_over3years,
        cc20[[3]]$clin_inc_all_over3years[vec_min_selections])




x = rep(seq(0,1,length=21)[vec_min_selections],6)
y = rep(seq(0,0.6,length=6),each=length(seq(0,1,length=21)[vec_min_selections]))

xcoords = unique(x)
ycoords = unique(y)

surface.mat1 = matrix(zf1,ncol=length(ycoords),nrow=length(xcoords),byrow=F)

## Figure 
plot.new()
par(mar=c(4,5.5,2,2))
# major tick size and direction, < 0 means outside

library(gplots)
library(colorRamps)
library(adegenet)

## Functions
filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), 
            zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), 
            nlevels = 20, 
            color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), 
            plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...)   {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    # further modified by Carey McGilliard and Bridget Ferris
    # to allow multiple plots on one page
    
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    # on.exit(par(par.orig))
    # w <- (3 + mar.orig[2]) * par("csi") * 2.54
    # par(las = las)
    # mar <- mar.orig
    plot.new()
    # par(mar=mar)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
      stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
      storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                    col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      box()
    if (missing(plot.title)) 
      title(...)
    else plot.title
    invisible()
  }


filled.legend <-function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                                       length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                          ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                          levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                          col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                          key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                          axes = TRUE, frame.plot = axes, ...){
  # modification of filled.contour by Carey McGilliard and Bridget Ferris
  # designed to just plot the legend
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  #  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #  on.exit(par(par.orig))
  #  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  #  par(las = las)
  #  mar <- mar.orig
  #  mar[4L] <- mar[2L]
  #  mar[2L] <- 1
  #  par(mar = mar)
  # plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
}


heat.colors_Rev = function (n, alpha = 1) {
  if ((n <- as.integer(n[1L])) > 0) {
    j <- n%/%4
    i <- n - j
    rev(c(rainbow(i, start = 0, end = 1/6, alpha = alpha), if (j > 
                                                               0) hsv(h = 1/6, s = seq.int(from = 1 - 1/(2 * j), 
                                                                                           to = 1/(2 * j), length.out = j), v = 1, alpha = alpha)))
  }
  else character()
}
matrix_funs = function(surface.matrix,minimum_val,maximum_val,
                       # upps,
                       # uni,
                       # levs,
                       colschoice,
                       cols_conts){
  filled.contour3(xcoords,
                  ycoords,
                  surface.matrix,
                  color=colschoice,
                  plot.axes = { axis(2, at = seq(0, 0.6, by = 0.1), seq(0, 0.6, by = 0.1),las=2)
                    axis(1, at = seq(0, 1, by = 0.2), labels = seq(0, 100, by = 20)) },
                  xlim = c(min(xcoords),max(xcoords)),
                  ylim = c(min(ycoords),max(ycoords)),
                  zlim = c(minimum_val,maximum_val))
  
  # # the contour part - draw iso-lines
  # contour(xcoords,
  #         ycoords,
  #         surface.matrix,
  #         color=blue2green,
  #         xlab = "",
  #         ylab = "",
  #         nlevels = levs, 
  #         levels = seq(0,upps,by=uni),
  #         xlim = c(min(xcoords),max(xcoords)),
  #         ylim = c(min(ycoords),max(ycoords)),
  #         zlim = c(minimum_val,maximum_val),
  #         add=TRUE,                 # add the contour plot to filled-contour,
  #         #thus making an overlay
  #         col = cols_conts         # color of overlay-lines
  # )
}

##mid left plot: Explanation figure, Pervalence estimate given different phiI and phiB estimates
par(mfrow=c(1,1))
par(new = "TRUE",  
    plt = c(0.2,0.85,0.2,0.85)  )                 # major tick size and direction, < 0 means outside


matrix_funs(surface.mat1,min(surface.mat1),max(surface.mat1),
            # upps=100,uni = 10,
            # levs = 6,
            colschoice = heat.colors,
            cols_conts="")
# text(0.95,0.95,"D",cex=1.2,col="")

#Add some figure labels
par(xpd=NA,cex = 1.1)
text(x = -0.15,y = 0.3,"Proportional reduction in carrying capacity",srt = 90,cex = 0.9)
text(x = 0.5,y = -0.1,"Mosquito survival at bioassay (%)",cex = 0.9)

# text(x = 0.5,y = 0.65,"Absolute efficacy against all-age clinical incidence",cex = 0.9)
text(x = 0.5,y = 0.65,"Efficacy against all-age clinical incidence",cex = 0.9)


######################################################################
#Add a legend:
par(new = "TRUE",
    plt = c(0.88,0.9,0.2,0.6),   # define plot region for legend
    las = 1,
    cex.axis = 1)
#
filled.legend(
  xcoords,
  ycoords,
  surface.mat1,
  color = heat.colors,
  plot.title = "Efficacy (%)",
  xlab = "",
  ylab = "",
  xlim = c(0,100),
  ylim = c(0,0.6),
  zlim = c(min(surface.mat1),max(surface.mat1)))
