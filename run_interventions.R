
library(deSolve) #package for running differential equations
library(dplyr)
source('karamoja_setup.R') #load setup and model code/functions
source('karamoja_model.R')

# for server, specify which array, and number of samples per array
n_samples <- 100 #100 arrays of 100 samples each
part <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) #used if running on the server
if(is.na(part)) {
  part <- 1 #used if not running on the server
}

## read posterior parameters and corresponding population sizes from calibration
# choose which version: weights, in_bounds, or threshold_25
version <- "weights" #in_bounds, threshold_25, or weights
sim.pars <- read.csv(paste0('calibration/select_results_', version, '.csv')) 
load(paste0("calibration/select_pop_", version, ".rda")) 
if(version=="weights") {
  sim.pop <- select3.pop #select3 if weights, select1 if threshold25, select2 if in_bounds
  #resample params according to weights
  resamples <- sample(row.names(sim.pars), size=n_samples, replace=T,
                      prob=sim.pars$weight)
} else if(version=="threshold_25") {
  sim.pop <- select1.pop
  #no weights
  resamples <- sample(sim.pars$X, size=n_samples, replace=T)
} else if(version=="in_bounds") {
  sim.pop <- select2.pop
  #no weights
  resamples <- sample(sim.pars$X, size=n_samples, replace=T)
}
sim.pars <- sim.pars[resamples, ]
sim.pop <- bind_rows(sim.pop)

# set up
pars$tpt <- pars$tpt_increase #set to post-burn-in TPT rate
pars.baseline <- pars #set baseline parameter values from setup.R
results_all <- vector('list', n_samples) #create list to store raw results



#Intervention 1: CXR in active case finding#
#asymptomatic treatment
covg1 <- rep(-log(1-0.2), n_samples) 
sens_cxr1 <- rnorm(n=n_samples, mean=0.888, sd=(0.936-0.838)/(2*1.96)) 
p_return_sputum1 <- rbeta(n=n_samples, shape1=55.77, shape2=23.90)
initiate1 <- rbeta(n=n_samples, shape1=14.98, shape2=1.66) 
tx_sub1 <- cbind(covg1*sens_cxr1*p_return_sputum1*initiate1, 
                 covg1*sens_cxr1*p_return_sputum1*initiate1) #asymptomatic treatment initiation rate when CXR added to case finding - stratified by age
#higher coverage variation
covg1_hi <- rep(-log(1-0.5), n_samples)
tx_sub1_hi <- cbind(covg1_hi*sens_cxr1*p_return_sputum1*initiate1, 
                    covg1_hi*sens_cxr1*p_return_sputum1*initiate1) #asymptomatic treatment initiation rate when CXR added to case finding - stratified by age
#also create a high coverage variation of background HSS (intervention 0), based on these parameters
tx0_hi <- cbind(covg1_hi*p_return_sputum1*initiate1,
                covg1_hi*p_return_sputum1*initiate1) #for symptomatic TB only

#Intervention 2, part 1: to intervention 1, add HHC investigation for all TB cases (including non-PBC)#
#component 1: treatment of non-PBC HH contacts w/ TB
covg2 <- rep(0.9, n_samples)
p_tb_if_hhc <- rnorm(n=n_samples, mean=1.6, sd=(1.9-1.4)/(2*1.96))/100 #prevalence of active TB among HHCs 
p_tb <- rnorm(n=n_samples, mean=504, sd=(652-335)/(2*1.96))/100000 #urban TB prevalence in Uganda
p_non_pbc <- 1 - sim.pars$p_pbc
p_hhc_if_tb <-   p_tb_if_hhc*sim.pars$p_hhc*p_non_pbc/p_tb 
p_return_sputum2 <- rbeta(n=n_samples, shape1=438.5267, shape2=89.81873) 
omega2_add <- cbind(covg2*p_hhc_if_tb*p_return_sputum2, covg2*p_hhc_if_tb*p_return_sputum2)  
#component 2: TPT for non-PBC HH contacts
p_hhc_if_ltbi <- sim.pars$rr_ltbi_hhc*sim.pars$p_hhc*p_non_pbc/
  (1-sim.pars$p_hhc*p_non_pbc + sim.pars$rr_ltbi_hhc*sim.pars$p_hhc*p_non_pbc) 

tpt2_add <- cbind(covg2*p_hhc_if_ltbi*sim.pars$tpt_init* #convert initiation ratio to rate
                    sim.pars$tpt_comp*(1-sim.pars$tpt_eff), 
                  covg2*p_hhc_if_ltbi*sim.pars$tpt_init*
                    sim.pars$tpt_comp*(1-sim.pars$tpt_eff)) #additional TPT/LTBI clearance rate with enhanced HH contact investigation, stratified by age
#Intervention 2, part 2: to intervention 1, add nutritional support for people on TB treatment
mu_fail_mult2 <- rbeta(n=n_samples, shape1=16.89429, shape2=8.32107) #relative risk death/failure during tx for ppl w/ undernutrition when nut. support added 


#Intervention 3: To intervention 2, add MUAC screening & TBI test during ACF campaigns - positive on both get offered TPT
covg3 <- covg1 # same 20% coverage as in intervention 1
muac_sens <- rbeta(n=n_samples, shape1=28*0.9286, shape2=28*(1-0.9286)) #from https://link.springer.com/article/10.1186/s12889-020-09294-0
tst_ltfu <- 0.3 #assume 30% do not return to have their TST read
tpt3_add <- cbind(covg3*muac_sens*(1-tst_ltfu)*sim.pars$tpt_init*sim.pars$tpt_comp*(1-sim.pars$tpt_eff), 
                  covg3*muac_sens*(1-tst_ltfu)*sim.pars$tpt_init*sim.pars$tpt_comp*(1-sim.pars$tpt_eff)) #additional TPT/LTBI clearance rate, stratified by age, for undernourished only
#higher ACF coverage variation
covg3_hi <- covg1_hi #50% coverage
tpt3_hi_add <- cbind(covg3_hi*muac_sens*tst_ltfu*sim.pars$tpt_init*sim.pars$tpt_comp*(1-sim.pars$tpt_eff), 
                  covg3_hi*muac_sens*tst_ltfu*sim.pars$tpt_init*sim.pars$tpt_comp*(1-sim.pars$tpt_eff)) 

#Intervention 4: to intervention 3, add nutritional support for everyone flagged by MUAC during ACF
covg4 <- covg1 #same 20% coverage as other interventions
uw_baseline <- rbeta(n=n_samples, shape1=1275, shape2=3311-1275) #prevalence of UW at baseline from RATIONS, among adults in intervention group
uw_end <- rbeta(n=n_samples, shape1=892, shape2=3204-892) #prevalence of UW post-intervention from RATIONS, among adults in intervention group
nut_efficacy <- (uw_baseline-uw_end)/uw_baseline
nut_trans4 <- covg4*muac_sens*nut_efficacy #% of undernourished transitioning to adequately nourished
#higher coverage ACF variation
covg4_hi <- covg1_hi #50% coverage
nut_trans4_hi <- covg4_hi*muac_sens*nut_efficacy #% of undernourished transitioning to adequately nourished

tic <- Sys.time() #keep track of how long simulations take

for(js in 1:nrow(sim.pars)){
  print(js)
  id <- sim.pars[js, "id"]
  print(id)
  number <- sim.pars[js, "number"]
  
  result <- data.frame()
  
  ## set model parameters based on each row in sim.pars
  pars <- pars.baseline
  pars$birth.rate[1] <- pars$birth.rate[2] <- as.numeric(sim.pars[js,'birth.rate'])
  
  pars$mu[1] <- as.numeric(sim.pars[js,'mu.1'])
  pars$mu[2] <- as.numeric(sim.pars[js,'mu.2'])
  odds_death_AA_2 <- as.numeric(sim.pars[js, 'mu.OR.2'])*as.numeric(sim.pars[js, 'mu.AA'])/(1-as.numeric(sim.pars[js, 'mu.AA']))
  odds_death_SA_2 <- as.numeric(sim.pars[js, 'mu.OR.2'])*as.numeric(sim.pars[js, 'mu.SA'])/(1-as.numeric(sim.pars[js, 'mu.SA']))
  pars$mu.AA <- c(as.numeric(sim.pars[js, 'mu.AA']), odds_death_AA_2/(1+odds_death_AA_2))
  pars$mu.SA <- c(as.numeric(sim.pars[js, 'mu.SA']), odds_death_SA_2/(1+odds_death_SA_2))
  pars$mu.tx.AA <- c(as.numeric(sim.pars[js, "mu.tx.AA.risk.1"]), as.numeric(sim.pars[js, "mu.tx.AA.risk.2"]))
  pars$mu.tx.SA <- c(as.numeric(sim.pars[js, "mu.tx.risk.1"]), as.numeric(sim.pars[js, "mu.tx.risk.2"]))
  pars$beta <- as.numeric(sim.pars[js,'beta'])
  pars$beta.age[1] <- as.numeric(sim.pars[js,'beta.age.1'])
  pars$beta.age[2] <- as.numeric(sim.pars[js,'beta.age.2'])
  pars$beta.stage[1] <- as.numeric(sim.pars[js,'beta.stage.1'])
  pars$beta.stage[2] <- as.numeric(sim.pars[js,'beta.stage.2'])
  pars$xi <- as.numeric(sim.pars[js,'xi'])
  pars$p[1] <- as.numeric(sim.pars[js,'p'])*as.numeric(sim.pars[js,'prog.age.1'])
  pars$p[2] <- as.numeric(sim.pars[js,'p'])*as.numeric(sim.pars[js,'prog.age.2'])
  pars$s <- as.numeric(sim.pars[js,'s'])
  pars$phi[1] <- as.numeric(sim.pars[js,'phi'])*as.numeric(sim.pars[js,'prog.age.1'])
  pars$phi[2] <- as.numeric(sim.pars[js,'phi'])*as.numeric(sim.pars[js,'prog.age.2'])
  pars$prog.risk[1] <- as.numeric(sim.pars[js, "prog.risk.1"])
  pars$prog.risk[2] <- as.numeric(sim.pars[js, "prog.risk.2"])
  pars$sus.risk[1] <- as.numeric(sim.pars[js, "sus.risk.1"])
  pars$sus.risk[2] <- as.numeric(sim.pars[js, "sus.risk.2"])
  pars$r1 <- as.numeric(sim.pars[js,'r1'])
  pars$r2 <- pars$r1 - as.numeric(sim.pars[js,'rdiff'])
  pars$w <- as.numeric(sim.pars[js,'w'])
  pars$k[1] <- as.numeric(sim.pars[js,'k.risk.1'])
  pars$k[2] <- as.numeric(sim.pars[js,'k.risk.2'])
  pars$k.AA[1] <- as.numeric(sim.pars[js,'k.AA.risk.1'])
  pars$k.AA[2] <- as.numeric(sim.pars[js,'k.AA.risk.2'])
  pars$omega[1] <- as.numeric(sim.pars[js,'omega'])*(1+5*sim.pars[js, 'increase.omega'])
  pars$omega[2] <- as.numeric(sim.pars[js,'omega'])*(1+5*sim.pars[js, 'increase.omega'])
  #omega increases to what it was at end of 5 yrs of improvements-then no further improvements
  pars$tpt <- rep(as.numeric(sim.pars[js, 'tpt_increase']), 2) #no further increase, just starts out at what it increased to at the end of burn-in
  
  ## set y.start based on sim.pop - skip over ID and time columns and add 0s for extra states at the end
  ystart <- c(unlist(sim.pop[sim.pop$id==id & sim.pop$time==max(sim.pop$time),3:(2+max(matind))]), 
              rep(0, n_extra_states))
  
  ## run for 10 years, with a monthly timestep
  h <- 1/12 #timestep =  1/12 year
  tmax <- 20 #run for 10 years
  ll <- tmax/h + 1
  burnin <- 0
  
  ## 0a. NO INTERVENTION (STATUS QUO/BASELINE HSS) - lower coverage ACF
  y.trans0 <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars, ystart=ystart))
  y.trans0 <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"=0, "intervention_full"="0",  "covg_lab"="NA", y.trans0)
  
  #add to list of output
  result <- rbind(result, y.trans0)
  
 
 
  ## 1a. ADD CHEST XRAY SCREENING TO CASE FINDING - lower coverage
  pars1a <- pars
  pars1a$tx_sub <- tx_sub1[js,]
  y.trans1a <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars1a, ystart=ystart))
  y.trans1a <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"=1, "intervention_full"="1a", "covg_lab"="Low", y.trans1a)
  #add to list of output
  result <- rbind(result, y.trans1a)
  
  ## 1b. ADD CHEST XRAY SCREENING TO CASE FINDING - higher coverage
  pars1b <- pars
  pars1b$tx_sub <- tx_sub1_hi[js,]
  y.trans1b <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars1b, ystart=ystart))
  y.trans1b <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"=1, "intervention_full"="1b", "covg_lab"="High", y.trans1b)
  #add to list of output
  result <- rbind(result, y.trans1b)
  
  ## 2a. To 1a, add: ENHANCED HH CONTACT INVESTIGATION & NUTRITIONAL SUPPORT DURING TREATMENT - with lower coverage ACF
  pars2a <- pars1a
  pars2a$omega <- pars2a$omega + omega2_add[js,]
  pars2a$tpt <- pars2a$tpt + tpt2_add[js,]
  pars2a$mu.tx.AA[2] <- pars2a$mu.tx.AA[2]*mu_fail_mult2[js]
  pars2a$mu.tx.SA[2] <- pars2a$mu.tx.SA[2]*mu_fail_mult2[js]
  fail.tx.2a <- (1-(pars2a$mu.tx.SA[2] + pars2a$k[2]))*mu_fail_mult2[js]
  pars2a$k[2] <- 1-(pars2a$mu.tx.SA[2] + fail.tx.2a)
  y.trans2a <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars2a, ystart=ystart))
  y.trans2a <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"=2, "intervention_full"="2a", "covg_lab"="Low", y.trans2a)
  #add to list of output
  result <- rbind(result, y.trans2a)
  
  ## 2b. To 1b, add: ENHANCED HH CONTACT INVESTIGATION & NUTRITIONAL SUPPORT DURING TREATMENT - with higher coverage ACF
  pars2b <- pars1b
  pars2b$omega <- pars2b$omega + omega2_add[js,]
  pars2b$tpt <- pars2b$tpt + tpt2_add[js,]
  pars2b$mu.tx.AA[2] <- pars2b$mu.tx.AA[2]*mu_fail_mult2[js]
  pars2b$mu.tx.SA[2] <- pars2b$mu.tx.SA[2]*mu_fail_mult2[js]
  fail.tx.2b <- (1-(pars2b$mu.tx.SA[2] + pars2b$k[2]))*mu_fail_mult2[js]
  pars2b$k[2] <- 1-(pars2b$mu.tx.SA[2] + fail.tx.2b)
  y.trans2b <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars2b, ystart=ystart))
  y.trans2b <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"=2, "intervention_full"="2b", "covg_lab"="High", y.trans2b)
  #add to list of output
  result <- rbind(result, y.trans2b)
  
  ## 3a. To 2a, add MUAC SCREEN, TBI TEST, TPT DURING CASE FINDING WAVES (lower coverage ACF)
  pars3a <- pars2a
  pars3a$tpt <- pars3a$tpt + tpt3_add[js,]
  y.trans3a <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars3a, ystart=ystart))
  y.trans3a <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="3", "intervention_full"="3a", "covg_lab"="Low", y.trans3a)
  result <- rbind(result, y.trans3a)
  
  ## 3b. To 2b, add MUAC SCREEN, TBI TEST, TPT DURING CASE FINDING WAVES (higher coverage ACF)
  pars3b <- pars2b
  pars3b$tpt <- pars3b$tpt + tpt3_hi_add[js,]
  y.trans3b <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars3b, ystart=ystart))
  y.trans3b <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="3", "intervention_full"="3b", "covg_lab"="High", y.trans3b)
  result <- rbind(result, y.trans3b)
  
  ## 4a. To 3a, add NUTRITIONAL SUPPORT FOR THOSE WITH LOW MUAC (lower coverage ACF)
  pars4a <- pars3a
  ystart4a <- ystart
  ystart4a[matind[,"Hi",]] <- ystart[matind[,"Hi",]]*(1-nut_trans4[js])
  ystart4a[matind[,"Lo",]] <- ystart[matind[,"Lo",]] + ystart[matind[,"Hi",]]*nut_trans4[js]
  y.trans4a <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars4a, ystart=ystart4a))
  y.trans4a <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="4", "intervention_full"="4a", "covg_lab"="Low", y.trans4a)
  result <- rbind(result, y.trans4a)
  
  ## 4b. To 3b, add NUTRITIONAL SUPPORT FOR THOSE WITH LOW MUAC (higher coverage ACF)
  pars4b <- pars3b
  ystart4b <- ystart
  ystart4b[matind[,"Hi",]] <- ystart[matind[,"Hi",]]*(1-nut_trans4_hi[js])
  ystart4b[matind[,"Lo",]] <- ystart[matind[,"Lo",]] + ystart[matind[,"Hi",]]*nut_trans4_hi[js]
  y.trans4b <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars4b, ystart=ystart4b))
  y.trans4b <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="4", "intervention_full"="4b", "covg_lab"="High", y.trans4b)
  result <- rbind(result, y.trans4b)
  
  #NEW: intervention "5" runs nutritional support (intervention 4) without the TPT first (intervention 3)
  #5a. To 2a, add NUTRITIONAL SUPPORT FOR THOSE WITH LOW MUAC (lower coverage ACF)
  pars5a <- pars2a
  ystart5a <- ystart
  ystart5a[matind[,"Hi",]] <- ystart[matind[,"Hi",]]*(1-nut_trans4[js])
  ystart5a[matind[,"Lo",]] <- ystart[matind[,"Lo",]] + ystart[matind[,"Hi",]]*nut_trans4[js]
  y.trans5a <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars5a, ystart=ystart5a))
  y.trans5a <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="5", "intervention_full"="5a", "covg_lab"="Low", y.trans5a)
  result <- rbind(result, y.trans5a)
  
  #5b. To 2b, add NUTRITIONAL SUPPORT FOR THOSE WITH LOW MUAC (lower coverage ACF)
  pars5b <- pars2b
  ystart5b <- ystart
  ystart5b[matind[,"Hi",]] <- ystart[matind[,"Hi",]]*(1-nut_trans4_hi[js])
  ystart5b[matind[,"Lo",]] <- ystart[matind[,"Lo",]] + ystart[matind[,"Hi",]]*nut_trans4_hi[js]
  y.trans5b <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars5b, ystart=ystart5b))
  y.trans5b <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="5", "intervention_full"="5b", "covg_lab"="High", y.trans5b)
  result <- rbind(result, y.trans5b)
  
  #add to summary results
  results_all <- rbind(results_all, result)
}

toc <- Sys.time() #keep track of how long simulations take
print(toc-tic)

write.csv(results_all, paste0("interventions/int_out_", version, "_", part, ".csv"), 
          row.names=F)

#also finalize and save sim.pars so we can link parameter sets to output
sim.pars <- sim.pars %>% mutate(run=1:n())
#add intervention parameters (only include those that vary)
sim.pars <- cbind(sim.pars, part, sens_cxr1, p_return_sputum1, initiate1,
                  p_hhc_if_tb, p_return_sputum2, p_hhc_if_ltbi, mu_fail_mult2,
                  muac_sens, nut_efficacy)
#save to file
write.csv(sim.pars, paste0("interventions/pars_out_", version, "_", part, ".csv"),
          row.names=F)

