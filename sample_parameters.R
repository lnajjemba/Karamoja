require('lhs') #package used for more efficient random sampling (Latin Hypercube Sampling)

nsmpl <- 1000000 #10,000 (originally 500,000) parameter set samples
N.p.hi <- 0.35  # proportion in high risk group (undernourished)
par.range <- read.csv('par_range.csv', row.names=1) #lower and upper bounds on parameters
au.pars <- randomLHS(nsmpl, ncol(par.range)) #500,000 x 28 samples, each uniform in [0, 1]

#convert from 0-1 to parameter lower and upper bounds
for(j in 1:ncol(par.range)){
  #mostly using uniform priors, but a few parameters are treated differently because we have more evidence on them
  if(names(par.range)[j]=="mu.OR.2") { #Odds ratio of death during treatment if undernourished
    samples <- qlnorm(p=au.pars[,j], meanlog=0.56, sdlog=0.22) #fits mean 1.8 [CI 1.1-2.7]
    samples[samples<1] <- 1 #truncate at 1
    au.pars[,j] <- samples 
  } else if(names(par.range)[j]=="unsuccess.RR.2") { #Relative risk of unfavorable treatment outcome if undernourished
    samples <- qlnorm(p=au.pars[,j], meanlog=0.40, sdlog=0.12) #fits mean 1.5 [CI 1.2-1.9]
    samples[samples<1] <- 1 #truncate at 1
    au.pars[,j] <- samples 
  } else if(names(par.range)[j]=="p_pbc") { #% of notifications that pulmonary bac-confirmed
    au.pars[,j] <- rbeta(n=nsmpl, shape1=88.23943, shape2=45.45667) 
  } else if(names(par.range)[j]=="p_hhc") { #% of people that are a HH contact of a notified TB case
    au.pars[,j] <- qnorm(p=au.pars[,j], mean=0.013, sd=(0.017-0.008)/(2*1.96))
  } else if(names(par.range)[j]=="rr_ltbi_hhc") { #relative prevalence of LTBI in HHs with TB
    au.pars[,j] <- qnorm(p=au.pars[,j], mean=3, sd=(3.3-2.7)/(2*1.96)) #fits mean 3.0 [2.7-3.3]
  } else if(names(par.range)[j]=="tpt_init") { #TPT initiation ratio
    au.pars[,j] <- qnorm(p=au.pars[,j], mean=0.84, sd=(0.847-0.838)/(2*1.96))
  } else if(names(par.range)[j]=="tpt_comp") { #TPT completion ratio
    au.pars[,j] <- qnorm(p=au.pars[,j], mean=0.929, sd=(0.949-0.902)/(2*1.96))
  } else if(names(par.range)[j]=="tpt_eff") { #TPT efficacy
    au.pars[,j] <- qbeta(p=au.pars[,j], shape1=3.08, shape2=1.66) #fits (to the extent possible) mean 65% [12-90%]
  }
  else { #rest are normally distributed (mostly with wide priors)
    au.pars[,j] <- par.range['low',j] + (par.range['high',j] - par.range['low',j])*au.pars[,j]
  }
}
colnames(au.pars) <- colnames(par.range)
#for treatment outcomes, override using multinomial distribution (so we never get failures < 0%)
multinom_samples <- rmultinom(n=nsmpl, size=250, prob=c(0.9, 0.035, 0.065))/250
au.pars[,"k"] <- multinom_samples[1,]
au.pars[,"mu.tx"] <- multinom_samples[2,]

#make sure mortality risk never increases while on treatment
flags <- which(au.pars[,"mu.SA"] < au.pars[,"mu.tx"])
while(length(flags)>0) {
  au.pars[flags, "mu.SA"] <- runif(n=length(flags), min=par.range["low","mu.SA"], max=par.range["high","mu.SA"])
  au.pars[flags, "mu.tx"] <- runif(n=length(flags), min=par.range["low", "mu.tx"], max=par.range["high","mu.tx"])
  flags <- which(au.pars[,"mu.SA"] < au.pars[,"mu.tx"])
  print(length(flags))
}

#calculate treatment outcomes - all stay as proportions (not rates)
#deaths
odds_die_all <- au.pars[,"mu.tx"]/(1-au.pars[,"mu.tx"])
odds_die_risk.1 <- odds_die_all/(N.p.hi*au.pars[,"mu.OR.2"] + (1-N.p.hi))
odds_die_risk.2 <- odds_die_risk.1*au.pars[,"mu.OR.2"]
p_die_risk.1 <- odds_die_risk.1/(1+odds_die_risk.1)
p_die_risk.2 <- odds_die_risk.2/(1+odds_die_risk.2)
#failures
p_unsuccess_all <- 1-au.pars[,"k"]
p_unsuccess_risk.1 <- p_unsuccess_all/(N.p.hi*(1*au.pars[,"unsuccess.RR.2"]) + (1-N.p.hi))
p_unsuccess_risk.2 <- p_unsuccess_risk.1*au.pars[,"unsuccess.RR.2"]
p_fail_risk.1 <- pmax(0, p_unsuccess_risk.1 - p_die_risk.1) #to ensure never negative
p_fail_risk.2 <- pmax(p_fail_risk.1, p_unsuccess_risk.2 - p_die_risk.2) #to ensure never negative
#successes
p_success_risk.1 <- 1 - (p_fail_risk.1 + p_die_risk.1)
p_success_risk.2 <- 1 - (p_fail_risk.2 + p_die_risk.2)
#for asymptomatic TB, deaths=0, failures=same as symptomatic, rest=successes
p_success_AA_risk.1 <- 1 - p_fail_risk.1
p_success_AA_risk.2 <- 1 - p_fail_risk.2
#add to au.pars
au.pars <- cbind(au.pars, "mu.tx.risk.1"=p_die_risk.1, "mu.tx.risk.2"=p_die_risk.2,
                 "k.risk.1"=p_success_risk.1, "k.risk.2"=p_success_risk.2,
                 "mu.tx.AA.risk.1"=0, "mu.tx.AA.risk.2"=0,
                 "k.AA.risk.1"=p_success_AA_risk.1, "k.AA.risk.2"=p_success_AA_risk.2)

#calculate TPT increase in last 2 years of burn in
p_hhc_if_ltbi <- au.pars[,"rr_ltbi_hhc"]*au.pars[,"p_hhc"]*au.pars[,"p_pbc"]/
  (1-au.pars[,"p_hhc"]*au.pars[,"p_pbc"] + au.pars[,"rr_ltbi_hhc"]*au.pars[,"p_hhc"]*au.pars[,"p_pbc"])
tpt_increase <- p_hhc_if_ltbi*(1-exp(-1*au.pars[,"tpt_init"]))* #convert initiation ratio to rate
  au.pars[,"tpt_comp"]*(1-au.pars[,"tpt_eff"])
#add to au.pars
au.pars <- cbind(au.pars, "tpt_increase"=tpt_increase)

#save as 100 separate parameter set files
npart <- 100 #divide this into 100 parts to be run in parallel (100 parameter sets each). 
#Set npart=1 if running calibration on your own computer (not the server)
for(jp in 1:npart){
	s.from <- (jp-1)*nsmpl/npart + 1
	s.to <- (jp)*nsmpl/npart
	write.csv(au.pars[s.from:s.to,], file=paste0('calibration/pars_',jp,'.csv'), row.names=F)
	}
