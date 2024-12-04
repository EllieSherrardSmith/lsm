# drat:::add("mrc-ide")
# install.packages("didehpc")
# setwd("Y:/Ellie/lsm")

options(didehpc.cluster = "fi--dideclusthn",
        didehpc.username = "esherrar")
didehpc::didehpc_config(cores = 5)

drat:::add("mrc-ide")

root <- "malaria_lsm"
packages <- c("malariasimulation", "malariaEquilibrium","tibble", "dplyr","tidyverse")
sources = c("cluster/source_files.R")
src <- conan::conan_sources(c("mrc-ide/malariasimulation", "mrc-ide/malariaEquilibrium"))
ctx <- context::context_save(root, packages = packages,
                             sources = sources, package_sources = src)

########
# options(
#   didehpc.username = "esherrar")#, ## ah1114
#   # didehpc.home = "Y:/Ellie")
# drat:::add("mrc-ide")
# 
# ## information on my computer
# malaria <- didehpc::path_mapping("Malaria", "Y:", "//fi--didenas1/Malaria", "Y:") ## 
# 
# config <- didehpc::didehpc_config(shares = malaria, use_rrq =  FALSE, 
#                                   cluster = "fi--didemrchnb", cores = 1,
#                                   parallel = FALSE)
# 
# packages <- c("malariasimulation", "malariaEquilibrium", "reshape2", "rio", "cali", "individual", "gridExtra", "ggplot2",
#               "tidyverse", "gghighlight", "ggpubr")
# 
# sources = c("cluster/source_files.R")
# 
# src <- conan::conan_sources(packages = packages, c("mrc-ide/cali", "mrc-ide/malariasimulation", "mrc-ide/malariaEquilibrium"))
# 
# ctx <- context::context_save(path = "new_contexts1", 
#                              packages = packages, 
#                              sources = sources,
#                              package_sources = src)

obj = didehpc::queue_didehpc(ctx)
# Don'tForget!

#Set up runs
all_set_up <- read.csv("R/malariasimulation/scenario_simulations/inputs_set_up.csv")

sims_data = data.frame(root_eir = all_set_up$root_eir,
                       phiB  = all_set_up$phiB,
                       phiI  = all_set_up$phiI,
                       itn_cov  = all_set_up$itn_cov,
                       dn0_med  = all_set_up$dn0_med,
                       rn0_med  = all_set_up$rn0_med,
                       gamman_med  = all_set_up$gamman_med,
                       cc_impact  = all_set_up$cc_impact)
for(i in 1:nrow(all_set_up)){
  sims_data$row_num[i] = i
} 

# ###
# test = lsm_sims_f(root_eir = sims_data[1,1],
#            phiB = sims_data[1,2],
#            phiI = sims_data[1,3],
#            itn_cov = sims_data[1,4],
#            dn0_med = sims_data[1,5],
#            rn0_med = sims_data[1,6],
#            gamman_med = sims_data[1,7],
#            cc_impact = sims_data[1,8],
#            row_num = sims_data[1,9])

head(all_set_up)
#Launch jobs
grp <- obj$enqueue_bulk(sims_data[1201:1300,],
                        lsm_sims_f, do_call = TRUE, progress = TRUE)