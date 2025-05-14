library(deSolve) #package for running differential equations
source('karamoja_setup.R') #load setup and model code/functions
source('karamoja_model.R')

#specify which parameter set samples file to use
part <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) #used if running on the server
if(is.na(part)) {
  part <- 1 #used if not running on the server
}

## read parameters for simulation
sim.pars <- read.csv(paste0('calibration/pars_',part,'.csv'))

pars.baseline <- pars #set baseline parameter values from setup.R
result <- vector('list', nrow(sim.pars)) #create list to store calibration results
#also track model outcomes that will be compared with calibration targets
res.sum <- array(NA,c(14,length(result)))
row.names(res.sum) <- c('prev', 'mort', 'notif_2022', 'notif_2017', 'prop_notif_u15', 
                        'rel_inc_undernut', 'prop_prev_sub', 'inc',
                        'ltbi_o15', 'ltbi_u15', 'prev_u15_low', 'prev_u15_high',
                        'prev_o15_low', 'prev_o15_high')

tic <- Sys.time() #keep track of how long calibration takes

for(js in 1:nrow(sim.pars)){
  print(js)
  
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
	pars$omega[1] <- as.numeric(sim.pars[js,'omega'])
	pars$omega[2] <- as.numeric(sim.pars[js,'omega'])
	pars$k[1] <- as.numeric(sim.pars[js,'k.risk.1'])
	pars$k[2] <- as.numeric(sim.pars[js,'k.risk.2'])
	pars$k.AA[1] <- as.numeric(sim.pars[js,'k.AA.risk.1'])
	pars$k.AA[2] <- as.numeric(sim.pars[js,'k.AA.risk.2'])
	pars$tpt_increase <- rep(as.numeric(sim.pars[js, 'tpt_increase']), 2)
	##print(pars)
	
	## run to steady state
	h <- 1 #timestep =  1 year
	tmax <- 500 #run for 500 years
	ll <- tmax/h + 1
	burnin <- 1
	y.trans <- ode_model(h=h, tmax=tmax, params=pars, ystart=c(as.vector(ystart),rep(0, n_extra_states)))
	#print(head(y.trans, 1))
	#print(tail(y.trans, 1))
	
	#for first few random parameter sets, check for burn-in - plot prevalence over time
	if(js %in% 1:5) {
	  pdf(width=10, height=10, file=paste0("calibration/burnin_", part, "_", js, ".pdf"))
	  plot(100000*rowSums(y.trans[c(1+as.vector(matind["AA",,]), 1+as.vector(matind["SA",,]))])/
	         rowSums(y.trans[1+as.vector(matind)]), ann=FALSE) 
	  dev.off()
	}

	## extract population size at the end of transience and adjust ystart accordingly
	N.trans <- rowSums(y.trans[ll,-c(1, 1+max(matind)+(1:n_extra_states))])
	ystart.cont <- c(as.numeric(y.trans[ll,-c(1, 1+max(matind)+(1:n_extra_states))])*N/N.trans, 
	                 rep(0, n_extra_states))
	
  ## run final 5 years, adjusting a couple of parameters (demographics and tx improvements)
	h <- 1/12 #timestep = monthly
	tmax <- 5 #run for 5 years to include tx improvements - may adjust this approach later
	burnin <- 0
	pars$increase.omega <- as.numeric(sim.pars[js,'increase.omega'])
	pars$t_tpt_increase <- 4 #increase in years 4 and 5 only. set this to INF to never increase the TPT rate

	y <- ode_model(h=h, tmax=tmax, params=pars, ystart=ystart.cont)
	
	## save results - raw model output
	result[[js]] <- list(y=y,ystart=ystart.cont,pars=pars)
	
	## also calculate output that will be compared to targets (used in calibration)
	nt <- 1/h #number of timesteps per year
	ll <- (tmax-1)*nt+1 #first time step of the last modeled year
	select <- (nt+1):(tmax*nt+1) #all time steps except the 1st modeled year
	z <- y[select,-c(1, 1+max(matind)+(1:n_extra_states))]
	#prevalence among adults, both risk groups 
	res.sum[1,js] <- as.numeric(rowSums(z[ll,as.vector(matind[c('AA', 'SA'),,2])])/rowSums(z[ll,as.vector(matind[,,2])]))*100000
	#annual mortality all ages and risk groups
	res.sum[2,js] <- as.numeric((diff(y[,30],nt)+ diff(y[,31],nt))/rowSums(z[,as.vector(matind[,,])]))[ll]*100000
	#annual notifications, both ages and risk groups - at end of 5 years of tx improvements 
	res.sum[3,js] <- as.numeric((diff(y[,32],nt)+ diff(y[,33],nt))/rowSums(z[,as.vector(matind[,,])]))[ll]*100000
	#annual notifications, both ages and risk groups - beginning of 5 years of tx improvements 
	res.sum[4,js] <- as.numeric((diff(y[,32],nt)+ diff(y[,33],nt))/rowSums(z[,as.vector(matind[,,])]))[1]*100000
	#% of notifications that are among children, both risk groups
	res.sum[5,js] <- as.numeric((diff(y[,38],nt)+ diff(y[,39],nt))/(diff(y[,32],nt)+ diff(y[,33],nt)))[ll]
	#relative incidence among people w/ under-nutrition, all ages (commented out version is relative prevalence)
	res.sum[6,js] <- as.numeric((diff(y[,41],nt))/(diff(y[,40],nt)))[ll]/
	  (rowSums(z[ll,as.vector(matind[,2,])])/rowSums(z[ll,as.vector(matind[,1,])]))
	#res.sum[6,js] <- as.numeric(rowSums(z[ll,as.vector(matind[c('AA', 'SA'),2,])])/rowSums(z[ll,as.vector(matind[,2,])]))/
	  #as.numeric(rowSums(z[ll,as.vector(matind[c('AA', 'SA'),1,])])/rowSums(z[ll,as.vector(matind[,1,])]))
	#proportion of prevalence that is subclinical, adults
	res.sum[7,js] <- as.numeric(rowSums(z[ll,as.vector(matind['AA',,2])])/rowSums(z[ll,as.vector(matind[c('AA', 'SA'),,2])]))
	#incidence among adults, both risk groups
	res.sum[8,js] <- as.numeric((diff(y[,40],nt)+ diff(y[,41],nt))/rowSums(z[,as.vector(matind[,,2])]))[ll]*100000

	#other outputs of interest to track
	#LTBI prevalence among adults, both risk groups
	res.sum[9,js] <- as.numeric(rowSums(z[ll,as.vector(matind[c('LE', 'LL'),,2])])/rowSums(z[ll,as.vector(matind[,,2])]))
	#LTBI prevalence among < 15s, both risk groups
	res.sum[10,js] <- as.numeric(rowSums(z[ll,as.vector(matind[c('LE', 'LL'),,1])])/rowSums(z[ll,as.vector(matind[,,1])]))
	#TB prevalence among < 15s, low-risk group
	res.sum[11,js] <- as.numeric(rowSums(z[ll,as.vector(matind[c('AA', 'SA'),1,1])])/rowSums(z[ll,as.vector(matind[,1,1])]))*100000
	#TB prevalence among < 15s, high-risk group (undernutrition)
	res.sum[12,js] <- as.numeric(rowSums(z[ll,as.vector(matind[c('AA', 'SA'),2,1])])/rowSums(z[ll,as.vector(matind[,2,1])]))*100000
	#TB prevalence among adults, low-risk group
	res.sum[13,js] <- as.numeric(rowSums(z[ll,as.vector(matind[c('AA', 'SA'),1,2])])/rowSums(z[ll,as.vector(matind[,1,2])]))*100000
	#TB prevalence among adults, high-risk group (undernutrition)
	res.sum[14,js] <- as.numeric(rowSums(z[ll,as.vector(matind[c('AA', 'SA'),2,2])])/rowSums(z[ll,as.vector(matind[,2,2])]))*100000
}

toc <- Sys.time()
print(toc-tic)

#save raw output
save(result, file=paste0('calibration/model_output_', part, '.rda'))

#combine summary output with params and save 
out_comb <- cbind(sim.pars, t(res.sum))
write.csv(out_comb, file=paste0('calibration/summary_output_', part, '.csv'))
