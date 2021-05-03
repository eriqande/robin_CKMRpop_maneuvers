library(tidyverse)
library(hexbin)
##library(rdist)
##CKMRpop::install_spip(Dir = system.file(package = "CKMRpop"))
##remotes::install_github("eriqande/CKMRpop", build_vignettes = TRUE)
library(CKMRpop)
#install_spip(Dir = system.file(package = "CKMRpop"))
#vignette("species_1_simulation", package = "CKMRpop")

SPD <- species_1_life_history

### define the scenario
survival = 0.7
alpha = 3  ## age at maturity
omega = 10 ## maximum age
adultlifespan = omega-alpha+1
phi = 1 ## ratio of Vk to kbar
##femalefecundity = c(0,0,alpha:omega)  ## fecundity proportional to age
femalefecundity = c(0,0,rep(1,adultlifespan))  ## constant fecundity
##femalefecundity = c(rep(0,9),1)  ## sweepstakes RS; only BOFFFs reproduce
cohort_size <- 20000
SampleSize = 500 ## fixed number to subsample from each cohort
samp_frac <- 0.1
SPD$`number-of-years` <- 56  # run the sim forward for 100 years
samp_start_year <- 51
samp_stop_year <- 55
range = paste(samp_start_year,"-",samp_stop_year,sep="")
NReps = 2
SPD[[4]] = c(0,0,1,1,1,1,1,1,1,1)  ## prob of reproducing at each age
##SPD[[4]] = c(0,0,1,0,1,0,1,0,1,0)
##SPD[[4]] = c(0,0,rep(0.1,8))
SPD[[6]] = femalefecundity
sampspace = paste(samp_frac," ")
samplingplan <- noquote(c(range, samp_frac, "0 0 0 0 0 0 0 0 0"))

  a=as.numeric(Sys.time())
  set.seed(a)

##Scenario B
SPD[[1]] = omega
SPD[[2]] = c(1,1,1,rep(survival,(SPD[[1]]-3)))
SPD[[3]] = SPD[[2]]

SPD[[5]] = SPD[[4]]

SPD[[7]] = SPD[[6]]
SPD[[8]] = "negbin"
SPD[[9]] = 1/phi
SPD[[10]] = SPD[[9]]
SPD[[11]] = -1
SPD[[12]] = 0.5

L <- leslie_from_spip(SPD, cohort_size)

# then we add those to the spip parameters
SPD$`initial-males` <- floor(L$stable_age_distro_fem)
SPD$`initial-females` <- floor(L$stable_age_distro_male)

# tell spip to use the cohort size
SPD$`cohort-size` <- paste("const", cohort_size, collapse = " ")
SPD$`fixed-cohort-size` <- ""  # define this flag and give it an empty string as an argument

sfspace = paste("0 ")
range = paste(samp_start_year,"-",samp_stop_year,sep="")
SPD$`discard-all` <- 0
#SPD$`gtyp-ppn-fem-pre` <- paste(range, "0 ", samp_frac, paste(rep(sfspace, SPD$'max-age' - 2), collapse = ""))
#SPD$`gtyp-ppn-male-pre` <- SPD$`gtyp-ppn-fem-pre`
SPD$`gtyp-ppn-fem-post` <- samplingplan
SPD$`gtyp-ppn-male-post` <- SPD$`gtyp-ppn-fem-post`

SPD$`cohort-size` <- paste("const", "200000", collapse = " ")

spip_dir <- run_spip(
  pars = SPD
)

