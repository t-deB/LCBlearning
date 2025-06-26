#First Theta is constant
set.seed(230205)
theta1 <- 1
omega <- c(NaN,10^-9)
a <- c(NaN,1+10 ^-9) 
jobs <- c(2,2)
runs <- 1000
tableOPT <- matrix(data = c(1,1.172,1.297,1.396,1.480,1.552,1.616,1.673,1.725,1.773,2,2.196,2.343,2.461,2.560,2.646,2.723,2.792,2.854,2.912,3,3.206,3.363,3.491,3.599,3.694,3.778,3.854,3.923,3.987,4,4.212,4.375,4.509,4.623,4.723,4.813,4.893,4.967,5.035,5,5.216,5.383,5.521,5.639,5.744,5.837, 5.921, 5.999, 6.07,6,6.218,6.388,6.530,6.651,6.759,6.855,6.942,7.023,7.09,7,7.220,7.392,7.536,7.66,7.77,7.869,7.959,8.041,8.118,8,8.221,8.396,8.541,8.667,8.779,8.88,8.972,9.057,9.135,9,9.223,9.398,9.545,9.673,9.787,9.889,9.983,10.069,10.149),byrow = TRUE, nrow = 9,ncol = 10)

delta <- c(0.95,0.9,0.75,0.5,0.25)
for(pa in 1:5){
priortheta <- c((a/omega))
priortheta[1] = theta1
pexp <- pgamma((a[2]-1)/omega[2],a[2],rate = omega[2])
priortheta[1] = qgamma(delta[pa]*(1-pexp) + pexp,a[2],rate = omega[2]) - 1
theta1 = priortheta[1]
print(theta1)
maxjobs <- max(jobs)
nclasses <- length(priortheta)



lseptscore <- 0
lcbscore <- c(0,0,0,0,0)
slcbscore <- c(0,0,0,0,0)
optscore <- 0
ocTheta <- matrix(data = NA, nrow = nclasses, ncol = max(jobs))

for (i in 1:runs){
  
  realtheta <- rep(0,nclasses)
  realtheta[1] = priortheta[1]
  for(j in 2:nclasses){
  realtheta[j] <- rgamma(1, a[j], rate = omega[j])
  }
  for(j in 1:nclasses){
    for(l in 1:jobs[j]){
  ocTheta[j,l] <- rexp(1, realtheta[j])
    }
  }
  
## SEPT
  septjobs <- rep(1,nclasses)
  septtime <- 0
  septproc <- 0
  septtheta <- priortheta
  
  for(j in 1:sum(jobs)){
    bestclass <- -1
    besttheta <- -1 * Inf
    for(l in 1:nclasses){
      if(septjobs[l] <= jobs[l]){
        if(septtheta[l] > besttheta){
          bestclass <- l
          besttheta <- septtheta[l]
        }
      }
    }
    septproc <- septproc + ocTheta[bestclass,septjobs[bestclass]]
    septtime <- septtime + septproc
    septjobs[bestclass] <- septjobs[bestclass]+1
  }
  septscore <- septscore + septtime
  
  
  ## l-sept
  lseptjobs <- rep(1,nclasses)
  lsepttime <- 0
  lseptproc <- 0
  lsepttheta <- priortheta
  la <- a
  lomega <- omega
  
  for(j in 1:sum(jobs)){
    bestclass <- -1
    besttheta <- Inf
    for(l in 1:nclasses){
      
      if(lseptjobs[l] <= jobs[l]){
        if( l == 1){
          bestclass <- 1
          besttheta <- 1/lsepttheta[1]
        }
        else if(lomega[l]/(la[l]-1) < besttheta){
          bestclass <- l
          besttheta <- lomega[l]/(la[l]-1)
        }
      }
    }
    lseptproc <- lseptproc + ocTheta[bestclass,lseptjobs[bestclass]]
    lsepttime <- lsepttime + lseptproc
    
    if(bestclass > 1){
      la[bestclass] <- la[bestclass] + 1
      lomega[bestclass] <- lomega[bestclass] + ocTheta[bestclass,lseptjobs[bestclass]]
      lsepttheta[bestclass] <- la[bestclass]/lomega[bestclass]
    }
    lseptjobs[bestclass] <- lseptjobs[bestclass]+1
  }
  lseptscore <- lseptscore + lsepttime


  # LCB 
  
  delta <- c(0.95,0.9,0.75,0.5,0.25)
  
  for(k in 1:length(delta)){
    lcbjobs <- rep(1,nclasses)
    lcbtime <- 0
    lcbproc <- 0
    lcbtheta <- priortheta
    lcba <- a
    lcbomega <- omega
    lcbMeasure <- rep(0,nclasses)
    for(l in 2:nclasses){
    pexp <- pgamma((lcba[l]-1)/lcbomega[l],lcba[l],rate = lcbomega[l])
    lcbMeasure[l] <- qgamma((1-pexp)*delta[k] + pexp,lcba[l],rate = lcbomega[l]) 
    }
    lcbMeasure[1] <- lcbtheta[1]
    lcbMeasure <- round(lcbMeasure,15)
    for(j in 1:sum(jobs)){
      bestclass <- -1
      bestMeasure <- -1 * Inf
      for(l in 1:nclasses){
        if(lcbjobs[l] <= jobs[l]){
          if(lcbMeasure[l] > bestMeasure){
            bestclass <- l
            bestMeasure <- lcbMeasure[l]
          }
        }
      }
      lcbproc <- lcbproc + ocTheta[bestclass,lcbjobs[bestclass]]
      lcbtime <- lcbtime + lcbproc
      
      if(bestclass > 1){
        lcba[bestclass] <- lcba[bestclass] + 1
        lcbomega[bestclass] <- lcbomega[bestclass] + ocTheta[bestclass,lcbjobs[bestclass]]
        lcbtheta[bestclass] <- lcba[bestclass]/lcbomega[bestclass]
        pexp <- pgamma((lcba[bestclass]-1)/lcbomega[bestclass],lcba[bestclass],rate = lcbomega[bestclass])
        lcbMeasure[bestclass] <- qgamma((1-pexp)*delta[k] + pexp,lcba[bestclass],rate = lcbomega[bestclass])
        lcbMeasure <- round(lcbMeasure,15)
      }
      lcbjobs[bestclass] <- lcbjobs[bestclass]+1
    }
    lcbscore[k] <- lcbscore[k] + lcbtime
  }
  
  ## Smart-LCB
  delta <- c(0.95,0.9,0.75,0.5,0.25)
  
  for(k in 1:length(delta)){
    slcbjobs <- rep(1,nclasses)
    slcbtime <- 0
    slcbproc <- 0
    slcbtheta <- priortheta
    slcba <- a
    slcbomega <- omega
    slcbMeasure <- rep(0,nclasses)
    z <- log(jobs-slcbjobs+1)/log(maxjobs)
    if(maxjobs == 1){
      z <- rep(0,nclasses)
    }
    
    for(l in 2:nclasses){
      
      pexp <- pgamma((slcba[l]-1)/slcbomega[l],slcba[l],rate = slcbomega[l])
      slcbMeasure[l] <- qgamma(z[l]*(1-pexp)*delta[k] + pexp,slcba[l],rate = slcbomega[l])
      
    }
    slcbMeasure[1] <- slcbtheta[1]
    slcbMeasure <- round(slcbMeasure,15)
    for(j in 1:sum(jobs)){
      
      bestclass <- -1
      bestMeasure <- -1 * Inf
      for(l in 1:nclasses){
        if(slcbjobs[l] <= jobs[l]){
          if(bestMeasure < slcbMeasure[l]){
            bestclass <- l
            bestMeasure <- slcbMeasure[l]
          }
        }
      }
      slcbproc <- slcbproc + ocTheta[bestclass,slcbjobs[bestclass]]
      
      slcbtime <- slcbtime + slcbproc
      
      if(bestclass > 1){
        slcba[bestclass] <- slcba[bestclass] + 1
        slcbomega[bestclass] <- slcbomega[bestclass] + ocTheta[bestclass,slcbjobs[bestclass]]
        slcbtheta[bestclass] <- slcba[bestclass]/slcbomega[bestclass]
        if(maxjobs > 1){
        z[bestclass] <- log(jobs[bestclass]-slcbjobs[bestclass])/log(maxjobs)
        }
        pexp <- pgamma((slcba[bestclass]-1)/lcbomega[bestclass],lcba[bestclass],rate = lcbomega[bestclass])
        slcbMeasure[bestclass] <- qgamma(z[bestclass]*(1-pexp)*delta[k] +pexp,slcba[bestclass],rate = slcbomega[bestclass])
        slcbMeasure <- round(slcbMeasure,15)
      }
      slcbjobs[bestclass] <- slcbjobs[bestclass]+1
    }
    slcbscore[k] <- slcbscore[k] + slcbtime
  }
  
  
  ## OPT
  opttime <- 0
  optproc <- 0
  if(nclasses == 2){
    if(jobs[2] <= 10){
      if(a[2] >= 2){
        if(a[2] + jobs[2] - 1 <= 10){
          opttime <- 0
          optproc <- 0
          optjobs <- rep(1,nclasses)
          opttheta <- priortheta
          opta <- a
          optomega <- omega
          
          
          for(j in 1:sum(jobs)){
            bestclass <- -1
            if(jobs[1] >= optjobs[1]){
              bestclass <- 1
            }
            else if(jobs[2] >= optjobs[2]){
              bestclass <- 2
            }
            if(jobs[2] >= optjobs[2]){
            if((optomega[2]*opttheta[1]) < tableOPT[(opta[2]-1),(jobs[2]+1-optjobs[2])]){
              
                bestclass <- 2
              }
            }
            
            optproc <-  optproc + ocTheta[bestclass,optjobs[bestclass]]
            opttime <- opttime + optproc
            if(bestclass == 2){
              opta[2] <- opta[2]+1
              optomega[2] <- optomega[2] + ocTheta[bestclass,optjobs[bestclass]]
              opttheta[2] <- opta[2]/optomega[2]
            }
            optjobs[bestclass] <- optjobs[bestclass] + 1
          }
          
        }
      }
    }
  }
  optscore <- optscore + opttime
}

#print(lseptscore)
#print(lcbscore)
#print(slcbscore)
#print(optscore)


print(lcbscore/lseptscore)
print(slcbscore/lseptscore)
#print(lcbscore/optscore)
#print(slcbscore/optscore)

}
