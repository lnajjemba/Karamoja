
library(deSolve) 
library(dplyr)
source('karamoja_setup.R') 
source('karamoja_model.R')


n_samples <- 100 
part <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) 
if(is.na(part)) {
  part <- 1 
}


version <- "weights" 
sim.pars <- read.csv(paste0('calibration/select_results_', version, '.csv')) 
load(paste0("calibration/select_pop_", version, ".rda")) 
if(version=="weights") {
  sim.pop <- select3.pop 
 
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
pars$tpt <- pars$tpt_increase 
pars.baseline <- pars 
results_all <- vector('list', n_samples) 


covg1 <- rep(-log(1-0.2), n_samples) 
sens_cxr1 <- rnorm(n=n_samples, mean=0.888, sd=(0.936-0.838)/(2*1.96)) 
p_return_sputum1 <- rbeta(n=n_samples, shape1=55.77, shape2=23.90)
initiate1 <- rbeta(n=n_samples, shape1=14.98, shape2=1.66) 
tx_sub1 <- cbind(covg1*sens_cxr1*p_return_sputum1*initiate1, 
                 covg1*sens_cxr1*p_return_sputum1*initiate1) 
covg1_hi <- rep(-log(1-0.5), n_samples)
tx_sub1_hi <- cbind(covg1_hi*sens_cxr1*p_return_sputum1*initiate1, 
                    covg1_hi*sens_cxr1*p_return_sputum1*initiate1) 
tx0_hi <- cbind(covg1_hi*p_return_sputum1*initiate1,
                covg1_hi*p_return_sputum1*initiate1) 


covg2 <- rep(0.9, n_samples)
p_tb_if_hhc <- rnorm(n=n_samples, mean=1.6, sd=(1.9-1.4)/(2*1.96))/100 
p_tb <- rnorm(n=n_samples, mean=504, sd=(652-335)/(2*1.96))/100000 
p_non_pbc <- 1 - sim.pars$p_pbc
p_hhc_if_tb <-   p_tb_if_hhc*sim.pars$p_hhc*p_non_pbc/p_tb 
p_return_sputum2 <- rbeta(n=n_samples, shape1=438.5267, shape2=89.81873) 
omega2_add <- cbind(covg2*p_hhc_if_tb*p_return_sputum2, covg2*p_hhc_if_tb*p_return_sputum2)  
p_hhc_if_ltbi <- sim.pars$rr_ltbi_hhc*sim.pars$p_hhc*p_non_pbc/
  (1-sim.pars$p_hhc*p_non_pbc + sim.pars$rr_ltbi_hhc*sim.pars$p_hhc*p_non_pbc) 

tpt2_add <- cbind(covg2*p_hhc_if_ltbi*sim.pars$tpt_init* 
                    sim.pars$tpt_comp*(1-sim.pars$tpt_eff), 
                  covg2*p_hhc_if_ltbi*sim.pars$tpt_init*
                    sim.pars$tpt_comp*(1-sim.pars$tpt_eff)) 
mu_fail_mult2 <- rbeta(n=n_samples, shape1=16.89429, shape2=8.32107) 

covg3 <- covg1 
muac_sens <- rbeta(n=n_samples, shape1=28*0.9286, shape2=28*(1-0.9286)) 
tst_ltfu <- 0.3 
tpt3_add <- cbind(covg3*muac_sens*(1-tst_ltfu)*sim.pars$tpt_init*sim.pars$tpt_comp*(1-sim.pars$tpt_eff), 
                  covg3*muac_sens*(1-tst_ltfu)*sim.pars$tpt_init*sim.pars$tpt_comp*(1-sim.pars$tpt_eff)) 

covg3_hi <- covg1_hi
tpt3_hi_add <- cbind(covg3_hi*muac_sens*tst_ltfu*sim.pars$tpt_init*sim.pars$tpt_comp*(1-sim.pars$tpt_eff), 
                  covg3_hi*muac_sens*tst_ltfu*sim.pars$tpt_init*sim.pars$tpt_comp*(1-sim.pars$tpt_eff)) #additional TPT/LTBI clearance rate, stratified by age, for undernourished only


covg4 <- covg1 
uw_baseline <- rbeta(n=n_samples, shape1=1275, shape2=3311-1275) 
uw_end <- rbeta(n=n_samples, shape1=892, shape2=3204-892) 
nut_efficacy <- (uw_baseline-uw_end)/uw_baseline
nut_trans4 <- covg4*muac_sens*nut_efficacy 
covg4_hi <- covg1_hi 
nut_trans4_hi <- covg4_hi*muac_sens*nut_efficacy 
tic <- Sys.time() 

for(js in 1:nrow(sim.pars)){
  print(js)
  id <- sim.pars[js, "id"]
  print(id)
  number <- sim.pars[js, "number"]
  
  result <- data.frame()
  
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
  pars$tpt <- rep(as.numeric(sim.pars[js, 'tpt_increase']), 2) 
  
  ystart <- c(unlist(sim.pop[sim.pop$id==id & sim.pop$time==max(sim.pop$time),3:(2+max(matind))]), 
              rep(0, n_extra_states))
  
  h <- 1/12
  tmax <- 20
  ll <- tmax/h + 1
  burnin <- 0
  
  y.trans0 <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars, ystart=ystart))
  y.trans0 <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"=0, "intervention_full"="0",  "covg_lab"="NA", y.trans0)
  
  result <- rbind(result, y.trans0)
  
  pars1a <- pars
  pars1a$tx_sub <- tx_sub1[js,]
  y.trans1a <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars1a, ystart=ystart))
  y.trans1a <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"=1, "intervention_full"="1a", "covg_lab"="Low", y.trans1a)
  result <- rbind(result, y.trans1a)
  
  pars1b <- pars
  pars1b$tx_sub <- tx_sub1_hi[js,]
  y.trans1b <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars1b, ystart=ystart))
  y.trans1b <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"=1, "intervention_full"="1b", "covg_lab"="High", y.trans1b)
  result <- rbind(result, y.trans1b)
  
  pars2a <- pars1a
  pars2a$omega <- pars2a$omega + omega2_add[js,]
  pars2a$tpt <- pars2a$tpt + tpt2_add[js,]
  pars2a$mu.tx.AA[2] <- pars2a$mu.tx.AA[2]*mu_fail_mult2[js]
  pars2a$mu.tx.SA[2] <- pars2a$mu.tx.SA[2]*mu_fail_mult2[js]
  fail.tx.2a <- (1-(pars2a$mu.tx.SA[2] + pars2a$k[2]))*mu_fail_mult2[js]
  pars2a$k[2] <- 1-(pars2a$mu.tx.SA[2] + fail.tx.2a)
  y.trans2a <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars2a, ystart=ystart))
  y.trans2a <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"=2, "intervention_full"="2a", "covg_lab"="Low", y.trans2a)
  result <- rbind(result, y.trans2a)
  
  pars2b <- pars1b
  pars2b$omega <- pars2b$omega + omega2_add[js,]
  pars2b$tpt <- pars2b$tpt + tpt2_add[js,]
  pars2b$mu.tx.AA[2] <- pars2b$mu.tx.AA[2]*mu_fail_mult2[js]
  pars2b$mu.tx.SA[2] <- pars2b$mu.tx.SA[2]*mu_fail_mult2[js]
  fail.tx.2b <- (1-(pars2b$mu.tx.SA[2] + pars2b$k[2]))*mu_fail_mult2[js]
  pars2b$k[2] <- 1-(pars2b$mu.tx.SA[2] + fail.tx.2b)
  y.trans2b <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars2b, ystart=ystart))
  y.trans2b <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"=2, "intervention_full"="2b", "covg_lab"="High", y.trans2b)
  result <- rbind(result, y.trans2b)
  
  pars3a <- pars2a
  pars3a$tpt <- pars3a$tpt + tpt3_add[js,]
  y.trans3a <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars3a, ystart=ystart))
  y.trans3a <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="3", "intervention_full"="3a", "covg_lab"="Low", y.trans3a)
  result <- rbind(result, y.trans3a)
  
  pars3b <- pars2b
  pars3b$tpt <- pars3b$tpt + tpt3_hi_add[js,]
  y.trans3b <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars3b, ystart=ystart))
  y.trans3b <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="3", "intervention_full"="3b", "covg_lab"="High", y.trans3b)
  result <- rbind(result, y.trans3b)
  
  pars4a <- pars3a
  ystart4a <- ystart
  ystart4a[matind[,"Hi",]] <- ystart[matind[,"Hi",]]*(1-nut_trans4[js])
  ystart4a[matind[,"Lo",]] <- ystart[matind[,"Lo",]] + ystart[matind[,"Hi",]]*nut_trans4[js]
  y.trans4a <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars4a, ystart=ystart4a))
  y.trans4a <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="4", "intervention_full"="4a", "covg_lab"="Low", y.trans4a)
  result <- rbind(result, y.trans4a)
  
  pars4b <- pars3b
  ystart4b <- ystart
  ystart4b[matind[,"Hi",]] <- ystart[matind[,"Hi",]]*(1-nut_trans4_hi[js])
  ystart4b[matind[,"Lo",]] <- ystart[matind[,"Lo",]] + ystart[matind[,"Hi",]]*nut_trans4_hi[js]
  y.trans4b <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars4b, ystart=ystart4b))
  y.trans4b <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="4", "intervention_full"="4b", "covg_lab"="High", y.trans4b)
  result <- rbind(result, y.trans4b)
  
  pars5a <- pars2a
  ystart5a <- ystart
  ystart5a[matind[,"Hi",]] <- ystart[matind[,"Hi",]]*(1-nut_trans4[js])
  ystart5a[matind[,"Lo",]] <- ystart[matind[,"Lo",]] + ystart[matind[,"Hi",]]*nut_trans4[js]
  y.trans5a <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars5a, ystart=ystart5a))
  y.trans5a <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="5", "intervention_full"="5a", "covg_lab"="Low", y.trans5a)
  result <- rbind(result, y.trans5a)
  
  pars5b <- pars2b
  ystart5b <- ystart
  ystart5b[matind[,"Hi",]] <- ystart[matind[,"Hi",]]*(1-nut_trans4_hi[js])
  ystart5b[matind[,"Lo",]] <- ystart[matind[,"Lo",]] + ystart[matind[,"Hi",]]*nut_trans4_hi[js]
  y.trans5b <- as.data.frame(ode_model(h=h, tmax=tmax, params=pars5b, ystart=ystart5b))
  y.trans5b <- cbind("id"=id, "number"=number, "part"=part, "run"=js, "intervention"="5", "intervention_full"="5b", "covg_lab"="High", y.trans5b)
  result <- rbind(result, y.trans5b)
  
  results_all <- rbind(results_all, result)
}

toc <- Sys.time()
print(toc-tic)

write.csv(results_all, paste0("interventions/int_out_", version, "_", part, ".csv"), 
          row.names=F)

sim.pars <- sim.pars %>% mutate(run=1:n())
sim.pars <- cbind(sim.pars, part, sens_cxr1, p_return_sputum1, initiate1,
                  p_hhc_if_tb, p_return_sputum2, p_hhc_if_ltbi, mu_fail_mult2,
                  muac_sens, nut_efficacy)
write.csv(sim.pars, paste0("interventions/pars_out_", version, "_", part, ".csv"),
          row.names=F)

