## DEFINE "STARTER" PARAMETER VALUES
n_extra_states <- 16 #number of extra outcomes to track in the model.R function

# birth and mortality rates
birth.rate <- c(0.045, 0.045) #birth rate (low risk, high risk)
mu <- c(0.006, 0.025)  # PLACEHOLDER VALUES; background (non-TB-specific) mortality rate (< 15 yrs, > 15 yrs)
mu.AA <- c(0, 0) #TB-specific mortality rate, asymptomatic untreated TB, by risk group
mu.SA <- 0.2*c(1, 1.5) #TB-specific mortality rate, symptomatic untreated TB, by risk group (low, high)
mu.tx.AA <- c(0, 0) #TB-specific mortality rate, asymptomatic TB on treatment
mu.tx.SA <- 0.1*c(1, 1.5) #TB-specific mortality rate, symptomatic TB on treatment, by risk group (low, high)

# population size, 2017 (before tx initiation improvements kick in)
N <- 1059100
N.p.hi <- 0.35  # proportion in high risk group
N.p.lo <- 1-N.p.hi 

# transmission rate (aka effective contact rate - rate at which an infected person infects a susceptible person )
beta <- 26
# beta gets multiplied by the % of the population infectious and the % susceptible, so the value of 26 is not meaningful by itself

# age-specific transmission scalar (relative infectiousness)
beta.age <- c(0, 1) #< 15, > 15
# stage-specific transmission scalar (relative infectiousness)
beta.stage <- c(0.2,1) #asymptomatic, symptomatic

# relative susceptibility,and progression risk of ppl w/ and w/out undernutrition
# these are all set to 1 for the low-risk (non-undernourished) group
sus.risk <- c(1, 1.2) #may set to 1 and only raise progression risk, but using placeholder value of 2 for now
prog.risk <- c(1, 1.2) #progression; 2 is a placeholder value

# protection against reinfection (relative risk of reinfection if already latently infected)
xi <- 0.4

# rapid progression rate, by age
p <- 0.035*c(0.725,1)

# stabilization rate (1 divided by duration of early latent period)
s <- 1/5

# reactivation rate, by age
phi <- 0.001*c(0.725,1)

# symptom progression
r1 <- 2.05

# symptom regression
r2 <- 2

# self-resolution
w <- 0.85

# treatment success proportions, by risk group
k <- 0.9*c(1, 1)   
k.AA <- k # treatment success among asymptomatic active

# treatment rates, by risk group
omega <- 0.8*c(1, 1)  #don't vary this by risk group

# increase in treatment rate at the end of burn-in (to replicate CAST)
increase.omega <- 0 #start at 0 during burn-in, adjust after burn-in in the other scripts

# tpt rate (LTBI to recovered)
tpt <- c(0, 0) #varied by age. value during burn-in (before increase at end of burn-in period)
tpt_increase <- c(0.05, 0.05) #varied by age - increase in years 4 and 5 at the end of the burn-in period
t_tpt_increase <- Inf #time at which to increase from TPT to TPT_increase

# subclinical treatment initiation rate
tx_sub <- c(0, 0) #varied by age. No subclinical treatment in calibration and most interventions

pars <- list(birth.rate=birth.rate,
			mu=mu,
			mu.AA=mu.AA,
			mu.SA=mu.SA,
			mu.tx.AA=mu.tx.AA,
			mu.tx.SA=mu.tx.SA,
			beta=beta,
			beta.stage=beta.stage,
			beta.age=beta.age,
			sus.risk=sus.risk,
			prog.risk=prog.risk,
			xi=xi,
			p=p,
			s=s,
			phi=phi,
			r1=r1,
			r2=r2,
			w=w,
			k=k,
			k.AA=k.AA,
			omega=omega,
			increase.omega=increase.omega,
			tpt=tpt,
			tpt_increase=tpt_increase,
			t_tpt_increase=t_tpt_increase,
			tx_sub=tx_sub
			)


## DEFINE A "HELPER ARRAY" OF INDICES TO BE USED IN THE MODEL

matind <- array(1:(6*2*2),c(6,2,2))
#dimensions are TB states (6 types), risk groups (low and high), ages (< 15 and > 15)
dimnames(matind)[[1]] <- c('U','LE', 'LL', 'AA', 'SA', 'R') #uninfected, early latent, late latent, asymptomatic active, symptomatic active, recovered
dimnames(matind)[[2]] <- c('Lo','Hi') #low-risk vs. high-risk
dimnames(matind)[[3]] <- c('0-14','15+') #< 15 vs. > 15



ystart <- array(1:(6*2*2),c(6,2,2))

#percent of low-risk and high-risk population by TB state (U, LE, LL, AA, SA, R)
ystart.epi <- matrix(c(c(0.3, 0.1, 0.5, 0.001, 0.001, 0.1-0.002),
                       c(0.3, 0.1, 0.5, 0.001, 0.001, 0.1-0.002)), nrow=2, byrow=T)

#percent of population low-risk vs. high-risk
ystart.space <- c(N.p.lo,N.p.hi)

#percent of population < 15 vs. > 15
ystart.age <- c(0.5,0.5)
			
#combine these percents together to estimate percent in each TB state-age-setting combination
for(tb_state in 1:6){
	for(setting in 1:2){
		for(age in 1:2){
			ystart[tb_state,setting,age] <- N*ystart.epi[setting,tb_state]*ystart.space[setting]*ystart.age[age]
		}
	}
}

