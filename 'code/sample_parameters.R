require('lhs') 

nsmpl <- 1000000 
N.p.hi <- 0.35  
par.range <- read.csv('par_range.csv', row.names=1)
au.pars <- randomLHS(nsmpl, ncol(par.range)) 

for(j in 1:ncol(par.range)){
  if(names(par.range)[j]=="mu.OR.2") { 
    samples <- qlnorm(p=au.pars[,j], meanlog=0.56, sdlog=0.22) 
    samples[samples<1] <- 1 
    au.pars[,j] <- samples 
  } else if(names(par.range)[j]=="unsuccess.RR.2") { 
    samples <- qlnorm(p=au.pars[,j], meanlog=0.40, sdlog=0.12) 
    samples[samples<1] <- 1 
    au.pars[,j] <- samples 
  } else if(names(par.range)[j]=="p_pbc") { 
    au.pars[,j] <- rbeta(n=nsmpl, shape1=88.23943, shape2=45.45667) 
  } else if(names(par.range)[j]=="p_hhc") { 
    au.pars[,j] <- qnorm(p=au.pars[,j], mean=0.013, sd=(0.017-0.008)/(2*1.96))
  } else if(names(par.range)[j]=="rr_ltbi_hhc") { 
    au.pars[,j] <- qnorm(p=au.pars[,j], mean=3, sd=(3.3-2.7)/(2*1.96)) 
  } else if(names(par.range)[j]=="tpt_init") { 
    au.pars[,j] <- qnorm(p=au.pars[,j], mean=0.84, sd=(0.847-0.838)/(2*1.96)) 
  } else if(names(par.range)[j]=="tpt_comp") { 
    au.pars[,j] <- qnorm(p=au.pars[,j], mean=0.929, sd=(0.949-0.902)/(2*1.96))
  } else if(names(par.range)[j]=="tpt_eff") {
    au.pars[,j] <- qbeta(p=au.pars[,j], shape1=3.08, shape2=1.66) 
  }
  else { 
    au.pars[,j] <- par.range['low',j] + (par.range['high',j] - par.range['low',j])*au.pars[,j]
  }
}
colnames(au.pars) <- colnames(par.range)
multinom_samples <- rmultinom(n=nsmpl, size=250, prob=c(0.9, 0.035, 0.065))/250
au.pars[,"k"] <- multinom_samples[1,]
au.pars[,"mu.tx"] <- multinom_samples[2,]

flags <- which(au.pars[,"mu.SA"] < au.pars[,"mu.tx"])
while(length(flags)>0) {
  au.pars[flags, "mu.SA"] <- runif(n=length(flags), min=par.range["low","mu.SA"], max=par.range["high","mu.SA"])
  au.pars[flags, "mu.tx"] <- runif(n=length(flags), min=par.range["low", "mu.tx"], max=par.range["high","mu.tx"])
  flags <- which(au.pars[,"mu.SA"] < au.pars[,"mu.tx"])
  print(length(flags))
}

odds_die_all <- au.pars[,"mu.tx"]/(1-au.pars[,"mu.tx"])
odds_die_risk.1 <- odds_die_all/(N.p.hi*au.pars[,"mu.OR.2"] + (1-N.p.hi))
odds_die_risk.2 <- odds_die_risk.1*au.pars[,"mu.OR.2"]
p_die_risk.1 <- odds_die_risk.1/(1+odds_die_risk.1)
p_die_risk.2 <- odds_die_risk.2/(1+odds_die_risk.2)
#failures
p_unsuccess_all <- 1-au.pars[,"k"]
p_unsuccess_risk.1 <- p_unsuccess_all/(N.p.hi*(1*au.pars[,"unsuccess.RR.2"]) + (1-N.p.hi))
p_unsuccess_risk.2 <- p_unsuccess_risk.1*au.pars[,"unsuccess.RR.2"]
p_fail_risk.1 <- pmax(0, p_unsuccess_risk.1 - p_die_risk.1) 
p_fail_risk.2 <- pmax(p_fail_risk.1, p_unsuccess_risk.2 - p_die_risk.2)
p_success_risk.1 <- 1 - (p_fail_risk.1 + p_die_risk.1)
p_success_risk.2 <- 1 - (p_fail_risk.2 + p_die_risk.2)
p_success_AA_risk.1 <- 1 - p_fail_risk.1
p_success_AA_risk.2 <- 1 - p_fail_risk.2
au.pars <- cbind(au.pars, "mu.tx.risk.1"=p_die_risk.1, "mu.tx.risk.2"=p_die_risk.2,
                 "k.risk.1"=p_success_risk.1, "k.risk.2"=p_success_risk.2,
                 "mu.tx.AA.risk.1"=0, "mu.tx.AA.risk.2"=0,
                 "k.AA.risk.1"=p_success_AA_risk.1, "k.AA.risk.2"=p_success_AA_risk.2)

p_hhc_if_ltbi <- au.pars[,"rr_ltbi_hhc"]*au.pars[,"p_hhc"]*au.pars[,"p_pbc"]/
  (1-au.pars[,"p_hhc"]*au.pars[,"p_pbc"] + au.pars[,"rr_ltbi_hhc"]*au.pars[,"p_hhc"]*au.pars[,"p_pbc"])
tpt_increase <- p_hhc_if_ltbi*(1-exp(-1*au.pars[,"tpt_init"]))*
  au.pars[,"tpt_comp"]*(1-au.pars[,"tpt_eff"])
au.pars <- cbind(au.pars, "tpt_increase"=tpt_increase)

npart <- 100  
for(jp in 1:npart){
	s.from <- (jp-1)*nsmpl/npart + 1
	s.to <- (jp)*nsmpl/npart
	write.csv(au.pars[s.from:s.to,], file=paste0('calibration/pars_',jp,'.csv'), row.names=F)
	}
