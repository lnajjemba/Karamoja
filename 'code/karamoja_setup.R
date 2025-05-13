
n_extra_states <- 16 


birth.rate <- c(0.045, 0.045) 
mu <- c(0.006, 0.025)  
mu.AA <- c(0, 0)
mu.SA <- 0.2*c(1, 1.5) 
mu.tx.AA <- c(0, 0) 
mu.tx.SA <- 0.1*c(1, 1.5) 

N <- 1059100
N.p.hi <- 0.35  
N.p.lo <- 1-N.p.hi 


beta <- 26

beta.age <- c(0, 1) #< 15, > 15

beta.stage <- c(0.2,1) 


sus.risk <- c(1, 1.2) 
prog.risk <- c(1, 1.2) 


xi <- 0.4


p <- 0.035*c(0.725,1)


s <- 1/5


phi <- 0.001*c(0.725,1)


r1 <- 2.05


r2 <- 2


w <- 0.85


k <- 0.9*c(1, 1)   
k.AA <- k 
omega <- 0.8*c(1, 1)  
increase.omega <- 0 
tpt <- c(0, 0) 
tpt_increase <- c(0.05, 0.05) 
t_tpt_increase <- Inf 


tx_sub <- c(0, 0) 

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




matind <- array(1:(6*2*2),c(6,2,2))

dimnames(matind)[[1]] <- c('U','LE', 'LL', 'AA', 'SA', 'R') 
dimnames(matind)[[2]] <- c('Lo','Hi') 
dimnames(matind)[[3]] <- c('0-14','15+') 



ystart <- array(1:(6*2*2),c(6,2,2))


ystart.epi <- matrix(c(c(0.3, 0.1, 0.5, 0.001, 0.001, 0.1-0.002),
                       c(0.3, 0.1, 0.5, 0.001, 0.001, 0.1-0.002)), nrow=2, byrow=T)


ystart.space <- c(N.p.lo,N.p.hi)


ystart.age <- c(0.5,0.5)
			

for(tb_state in 1:6){
	for(setting in 1:2){
		for(age in 1:2){
			ystart[tb_state,setting,age] <- N*ystart.epi[setting,tb_state]*ystart.space[setting]*ystart.age[age]
		}
	}
}

