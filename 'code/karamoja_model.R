
model <- function(t, X, params){
	with(as.list(params),{ 
		
		DX <- array(0, length(as.vector(matind)) + n_extra_states) 

		names.nat.hist <- c('U','LE', 'LL', 'AA', 'SA', 'R') 
		names.setting <- c('Lo','Hi')
		names.age <- c('0-14','15+')
		
		
		if(burnin==0) {
		  if(t>=t_tpt_increase) {
		    tpt <- tpt_increase
		  }
		}

		
		for(jn in 1:length(names.nat.hist)){
			for(js in 1:length(names.setting)){
				for(ja in 1:length(names.age)){
					DX[matind[jn,js,ja]] <- DX[matind[jn,js,ja]] - mu[ja]*X[matind[jn,js,ja]]
					
					if(names.nat.hist[jn] == 'AA'){ 
						DX[matind[jn,js,ja]] <- DX[matind[jn,js,ja]] - mu.AA[js]*X[matind[jn,js,ja]] -
						  tx_sub[ja]*mu.tx.AA[js]*X[matind[jn,js,ja]]
					}
					else if	(names.nat.hist[jn] == 'SA'){ 
						DX[matind[jn,js,ja]] <- DX[matind[jn,js,ja]] - mu.SA[js]*X[matind[jn,js,ja]] -
						  (1+t*increase.omega)*omega[js]*mu.tx.SA[js]*X[matind[jn,js,ja]] 
					}		
				}
			}
		}
		
		
		if(burnin==1) { 
		  size.pop <- c(sum(X[as.vector(matind[,1,])]) , sum(X[as.vector(matind[,2,])]))
		  deaths.rate <- -1*c(sum(DX[as.vector(matind[,1,])]), sum(DX[as.vector(matind[,2,])]))/size.pop
		  births.adj.2 <- birth.rate[1] + size.pop[1]*(deaths.rate[2] - deaths.rate[1])/(size.pop[1] + size.pop[2])
		  births.adj.1 <- births.adj.2 + deaths.rate[1] - deaths.rate[2]
		  births.adj <- c(births.adj.1, births.adj.2)
		  for(js in 1:length(names.setting)){
		    DX[matind['U',js,'0-14']] <- DX[matind['U',js,'0-14']] + births.adj[js]*sum(X[as.vector(matind[,js,])])
		  }
		} else { 
		  for(js in 1:length(names.setting)){
		    DX[matind['U',js,'0-14']] <- DX[matind['U',js,'0-14']] + birth.rate[js]*sum(X[as.vector(matind[,js,])]) 
		  }
		}

		
		for(jn in 1:length(names.nat.hist)){
			for(js in 1:length(names.setting)){
				DX[matind[jn,js,'0-14']] <- DX[matind[jn,js,'0-14']] - 1/15*X[matind[jn,js,'0-14']] 
				DX[matind[jn,js,'15+']] <- DX[matind[jn,js,'15+']] + 1/15*X[matind[jn,js,'0-14']] 
			}
		}

		
		foi <- 0 
		for(jn in c(1:2)){ 
			for(ja in 1:length(names.age)){ 
				tb <- c(4,5)
				foi <- foi + beta.stage[jn]*beta.age[ja]*sum(X[as.vector(matind[tb[jn],,ja])]) 				
			}
		}
		
		lambda <- c(0,0)
		lambda[1] <- sus.risk[1]*beta*foi/(sum(X[as.vector(matind)])) 
		lambda[2] <- sus.risk[2]*beta*foi/(sum(X[as.vector(matind)]))
		
	
		for(js in 1:length(names.setting)){
			for(ja in 1:length(names.age)){
		
				 
				DX[matind['U',js,ja]] <- DX[matind['U',js,ja]] - lambda[js]*X[matind['U',js,ja]]
				DX[matind['LE',js,ja]] <- DX[matind['LE',js,ja]] + lambda[js]*X[matind['U',js,ja]]
		
				
				DX[matind['R',js,ja]] <- DX[matind['R',js,ja]] - xi*lambda[js]*X[matind['R',js,ja]]
				DX[matind['LE',js,ja]] <- DX[matind['LE',js,ja]] + xi*lambda[js]*X[matind['R',js,ja]]
		
				DX[matind['LL',js,ja]] <- DX[matind['LL',js,ja]] - xi*lambda[js]*X[matind['LL',js,ja]]
				DX[matind['LE',js,ja]] <- DX[matind['LE',js,ja]] + xi*lambda[js]*X[matind['LL',js,ja]]	
		
				
				DX[matind['LE',js,ja]] <- DX[matind['LE',js,ja]] - prog.risk[js]*p[ja]*X[matind['LE',js,ja]]
				DX[matind['AA',js,ja]] <- DX[matind['AA',js,ja]] + prog.risk[js]*p[ja]*X[matind['LE',js,ja]]

				
				DX[matind['LE',js,ja]] <- DX[matind['LE',js,ja]] - s*X[matind['LE',js,ja]]
				DX[matind['LL',js,ja]] <- DX[matind['LL',js,ja]] + s*X[matind['LE',js,ja]]
				
				
				DX[matind['LE',js,ja]] <- DX[matind['LE',js,ja]] - tpt[ja]*X[matind['LE',js,ja]]
				DX[matind['LL',js,ja]] <- DX[matind['LL',js,ja]] - tpt[ja]*X[matind['LL',js,ja]]
				DX[matind['R',js,ja]] <- DX[matind['R',js,ja]] + tpt[ja]*(X[matind['LE',js,ja]] + X[matind['LL',js,ja]])

				
				DX[matind['LL',js,ja]] <- DX[matind['LL',js,ja]] - prog.risk[js]*phi[ja]*X[matind['LL',js,ja]]
				DX[matind['AA',js,ja]] <- DX[matind['AA',js,ja]] + prog.risk[js]*phi[ja]*X[matind['LL',js,ja]]

				
				DX[matind['AA',js,ja]] <- DX[matind['AA',js,ja]] - r1*X[matind['AA',js,ja]]
				DX[matind['SA',js,ja]] <- DX[matind['SA',js,ja]] + r1*X[matind['AA',js,ja]]

				
				DX[matind['SA',js,ja]] <- DX[matind['SA',js,ja]] - r2*X[matind['SA',js,ja]]
				DX[matind['AA',js,ja]] <- DX[matind['AA',js,ja]] + r2*X[matind['SA',js,ja]]
		
				
				DX[matind['AA',js,ja]] <- DX[matind['AA',js,ja]] - w*X[matind['AA',js,ja]]
				DX[matind['R',js,ja]] <- DX[matind['R',js,ja]] + w*X[matind['AA',js,ja]]

				
				DX[matind['SA',js,ja]] <- DX[matind['SA',js,ja]] - (1+t*increase.omega)*omega[js]*(1-mu.tx.SA[js])*k[js]*X[matind['SA',js,ja]]
				DX[matind['R',js,ja]] <- DX[matind['R',js,ja]] + (1+t*increase.omega)*omega[js]*(1-mu.tx.SA[js])*k[js]*X[matind['SA',js,ja]]
				
				
				DX[matind['AA',js,ja]] <- DX[matind['AA',js,ja]] - tx_sub[ja]*(1-mu.tx.AA[js])*k.AA[js]*X[matind['AA',js,ja]]
				DX[matind['R',js,ja]] <- DX[matind['R',js,ja]] + tx_sub[ja]*(1-mu.tx.AA[js])*k.AA[js]*X[matind['AA',js,ja]]
				
			}
		}
		
		
		for(ja in 1:length(names.age)){
		
			
			DX[length(as.vector(matind))+1] <- DX[length(as.vector(matind))+1] + lambda[1]*X[matind['U',1,ja]] 
			DX[length(as.vector(matind))+2] <- DX[length(as.vector(matind))+2] + lambda[2]*X[matind['U',2,ja]] 
		
			
			DX[length(as.vector(matind))+3] <- DX[length(as.vector(matind))+3] + 
			  prog.risk[1]*p[ja]*X[matind['LE',1,ja]] + prog.risk[1]*phi[ja]*X[matind['LL',1,ja]]
			DX[length(as.vector(matind))+4] <- DX[length(as.vector(matind))+4] + 
			  prog.risk[2]*p[ja]*X[matind['LE',2,ja]] + prog.risk[2]*phi[ja]*X[matind['LL',2,ja]]
		
			
			DX[length(as.vector(matind))+5] <- DX[length(as.vector(matind))+5] + mu.AA[js]*X[matind['AA',1,ja]] +  mu.SA[js]*X[matind['SA',1,ja]] + 
			  tx_sub[ja]*mu.tx.AA[js]*X[matind['AA',1,ja]] + (1+t*increase.omega)*omega[js]*mu.tx.SA[js]*X[matind['SA',1,ja]]
			DX[length(as.vector(matind))+6] <- DX[length(as.vector(matind))+6] + mu.AA[js]*X[matind['AA',2,ja]] +  mu.SA[js]*X[matind['SA',2,ja]] + 
			  tx_sub[ja]*mu.tx.AA[js]*X[matind['AA',2,ja]] + (1+t*increase.omega)*omega[js]*mu.tx.SA[js]*X[matind['SA',2,ja]]
			
			
			DX[length(as.vector(matind))+7] <- DX[length(as.vector(matind))+7] + (1+t*increase.omega)*omega[1]*X[matind['SA',1,ja]] +
			  tx_sub[1]*X[matind['AA',1,ja]]
			DX[length(as.vector(matind))+8] <- DX[length(as.vector(matind))+8] + (1+t*increase.omega)*omega[2]*X[matind['SA',2,ja]] + 
			  tx_sub[2]*X[matind['AA',2,ja]]
			
			
			DX[length(as.vector(matind))+9] <- DX[length(as.vector(matind))+9] + prog.risk[1]*p[ja]*X[matind['LE',1,ja]] 
			DX[length(as.vector(matind))+10] <- DX[length(as.vector(matind))+10] + prog.risk[2]*p[ja]*X[matind['LE',2,ja]] 
			
			
			DX[length(as.vector(matind))+11] <- DX[length(as.vector(matind))+11] + prog.risk[1]*phi[ja]*X[matind['LL',1,ja]] 
			DX[length(as.vector(matind))+12] <- DX[length(as.vector(matind))+12] + prog.risk[2]*phi[ja]*X[matind['LL',2,ja]]
		}
		
		
		DX[length(as.vector(matind))+13] <- DX[length(as.vector(matind))+13] + (1+t*increase.omega)*omega[1]*X[matind['SA',1,1]] +
		  tx_sub[1]*X[matind['AA',1,1]]
		DX[length(as.vector(matind))+14] <- DX[length(as.vector(matind))+14] + (1+t*increase.omega)*omega[2]*X[matind['SA',2,1]] +
		  tx_sub[2]*X[matind['AA',2,1]]
		
		
		
		DX[length(as.vector(matind))+15] <- DX[length(as.vector(matind))+15] + 
		  prog.risk[1]*p[2]*X[matind['LE',1,2]] + prog.risk[1]*phi[2]*X[matind['LL',1,2]]
		DX[length(as.vector(matind))+16] <- DX[length(as.vector(matind))+16] + 
		  prog.risk[2]*p[2]*X[matind['LE',2,2]] + prog.risk[2]*phi[2]*X[matind['LL',2,2]]
		
	
		list(DX) 
	})
}



ode_model <- function (h, tmax, params, ystart){
  sol <- as.data.frame(
                     lsoda(
                           ystart,
                           times=seq(0,tmax,by=h),
                           func=model,
                           parms=params,
                           rtol=1e-8,
                           atol=1e-8
                           )
                     )  
}
