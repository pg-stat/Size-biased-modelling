# MODEL FITTING 
#CONSIDERING THE CONSTRAINT ON pij
#CONSIDERING Sij ~ Neg-binomial(lambdai, varS)

# IN THE ISRO MODEL OF DEY ET AL, 
# IT IS ASSUMED THAT, THE DETECTED BUGS ARE BORN IN A PARTICULAR MISSION AND THEY ARE NOT LEAKED FROM EARLIER MISSIONS
# IT IS ALSO ASSUMED THAT, THE UNDETECTED BUGS CAN BE ANYWHERE IN ANY TRIPLET OF (MISSION, PHASE, SOFTWARE),
# AND HENCE THE NUMBER OF TEST CASES THAT CAN DETECT THE UNDETECTED BUGS ARE A SUM OF ALL TEST CASES
# FROM ALL MISSIONS AND PHASES.


# IDEALLY, SOME NEW BUGS ARE BORN IN DIFFERENT MISSIONS. THIS REPRODUCTIVE RATE IS VERY HARD TO ESTIMATE.
# HENCE IT IS BEST TO ASSUME A CLOSED POPULATION MODEL IN SIMULATION STUDY (LIKE WHAT WE DID FOR OUT DEY ET AL.).

# # VERSION INFO:  risk <- 1-exp(-S/TCmax) # M
# pmat <- 1 - exp(-TC)

# J = 20, K = 8, N = 400, lam0 = 50
# J = 20, K = 8, N = 400, lam0 = 80
# J = 20, K = 8, N = 200, lam0 = 80
# J = 30, K = 8, N = 200, lam0 = 80


{##--DO ALL
  
  {##--DO ALL
    
    
    rm(list=ls())
    gc()
    # cat("\014") # CLEAR CONSOLE
    
    #----- PACKAGES -----#
    #====================#
    # install.packages("mcmcr")
    # install.packages("coda")
    # install.packages("actuar")
    # install.packages("extraDistr")
    # install.packages("abind")
    # install.packages("xlsx")
    # install.packages("goeveg")
    
    #library(xlsx)
    library(actuar)  #rztpois
    library(extraDistr) #rtpois
    library(coda)       #convergence diagnostics
    library(mcmcr)
    library(abind)
    library(basicMCMCplots)
    library(mcmcse)
    library(goeveg)
    
    # list.of.packages <- c('actuar', 'extraDistr', 'coda', 'mcmcr', 'abind', 'basicMCMCplots')
    # new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    # if(length(new.packages) > 0) install.packages(new.packages, dependencies = T)
    
    #================#
    #----- DATA -----#
    #================#
    # path1 <- "D:/myCodes/RovquantSD_Codes/Shared by Pallabi/Pallabi230210"
    # WD <- "D:/Results/Pallabi_ISRO"
    
    path1 <- "D:/1. PALLABI/PHD/ISRO/Program/Model fitting/PG_new/Final/sim"
    WD <- "D:/1. PALLABI/PHD/ISRO/Program/Model fitting/PG_new/Final/Results"
    
    
    
    
    # NOTE: to run the following code, you must bring in the functions
    # in PART 2 and the dataset in PART 3 (see below)
    
    
    
    # Simulate a dataset and estimate parameters using algorithm 1
    # The spNmix function *does not* update the latent z[i,j,k] variables,
    # ie, it is the unconditional-on-z formulation.
     set.seed(343892)
    
    
    ## R functions required to fit the model using the code in PART 1
    
    
    ilogit <- function(inValues) {
      a <- 1.0 / (1.0 + exp(-inValues))
      return(a)
    }
    
    logit = function(x) {
      a <- log(x / (1-x))
      return(a)
    }
    
    
    
    
    niters <- 10000
    J <- 30 # No. of missions
    K <- 8 # No. of phases
    TC <- matrix(sample(0:50, J*K, replace = T), J, K ) # rpois(J*K, lambda = 20) # Test cases
    N <- 100 # True no. of bugs
    M <- 400 # max possible no. of bugs
    # initial values
    TCmax = max(TC)
    lam0 <- rep(80,M) # Mean size of each bug
    nu=1.5
    
    ##------------------------------
    ## Negative binomial
    ##------------------------------
    # size <- 200
    varS <- 200 ### varS <- lam0+lam0^2/size
    size <- 80^2 / (varS - 80)
    size
    S0 = rnbinom(N, size = 80^2 / (varS - 80), mu = 80) # N x 1
    #S0 = rpois(N, lambda = 80) # N x 1
    var(S0)
    
    # S=matrix(rtpois((J*K),lam0,a=TC,b=500),J,K)
    S0[is.infinite(S0)] <- max(TC)*2        
    S0[S0<=min(TC)] <- max(TC)*2
    # S0 = max(TC)*2
    
    r <- 2e-1 # 1/(J*K) # 
    a <- 2e-2 # 1-1/(J*K+1)
    tune=c(abs(logit(a))/10, 0.1, 2,5,abs(logit(r))/10)
    
    z <- rep(0, M)
    
    
    z[1:N] <- 1
    
    psi <- N/M
    
    ##=================================================================
    ##-- The detection prob. function for multinomial distribution
    ##=================================================================
    
    
    p.multi <- function(pmat,S){ # N x (J*K + 1)
      probs <- c(pmat)
      risk <- 1-exp(-(S^nu)/TCmax) # M x 1
      p <- do.call(rbind, lapply(1:length(S), function(i){ # M x (J*K + 1)
        c(risk[i] * probs / sum(probs), 1 - risk[i]) 
        # The last entry is for the event where the bug is not detected
        # in any of the mission and in any of the phases
      }))
      # p <- c(a * probs / sum(probs), 1 - a)
      return(p)
    }
    
    b0 <- 1
    b1 <- -1
    b2 <- -3.2
    
    ##=========================
    ##-- The detection prob.
    ##=========================
    pmat <- 1 - exp(-TC)
    p <- p.multi(pmat,S0) # N x (J*K + 1)
    
    p[,J*K+1]
    
    ##=========================
    ##-- The detection history
    ##=========================
    
    data.full0 <- NULL
    w.arr0 <- array(NA, dim= c(N,J,K),
                    dimnames = list(paste0("Bugs",1:N), paste0("Missions",1:J), paste0("Phases",1:K)))
    u0 <- rep(0, N)
    ii <- 1
    
    for(ii in 1:N){
      
      this.bug <- c(rmultinom(1, size = 1, prob = p[ii,])) # (J*K + 1)
      
      data.full0 <- rbind(data.full0, this.bug)
      u0[ii] <- this.bug[length(this.bug)] # u_i = 1 if i-th bug is not detected
      
      w.arr0[ii,,] <- matrix(this.bug[1:(length(this.bug)-1)],# # detection history in missions and phases
                             nrow=J, ncol=K)
    }##ii
    
    detected0 <- apply(w.arr0, 1, sum) > 0
    
    whichNotDetected0 <- which(detected0==0)
    n <- sum(detected0)
    
    all.equal(1-as.numeric(detected0), u0) # TRUE
    
    
    ##=========================
    ##-- The detected data set
    ##=========================
    data.full <- data.full0[detected0,]
    u <- u0[detected0] # all zeros since u_1 = 1 if i-th bug is not detected
    w.arr <- w.arr0[detected0,,]
    S <- S0[detected0]
    
    ##=========================
    ##-- Data augmentation
    ##=========================
    data.full.aug <- rbind(data.full, matrix(c(rep(0, J*K), 1), M-n, J*K+1, byrow = T))
    u.aug <- c(u, rep(1, M-n))
    w.arr.aug <- abind(w.arr, array(0, c(M-n,J,K)), along = 1)
    S.aug <- c(S0[detected0], S0[!detected0], lam0[(N+1):M])
    
    y.mat <- apply(w.arr, c(2,3), sum)
    
    
    ##============================
    ##-- The likelihood functions
    ##============================
    
    llfn <- function(data0, p0, z0){
      ll <- sum(unlist(lapply(1:nrow(data0),function(i){
        out <- 0
        if(z0[i]==1){
          out <- dmultinom(data0[i,],size=1, prob=p0[i,], log=TRUE)
        }
        return(out)
      })))
      return(ll)
    }
    
    
    
    # llfn2.slow <- function(data0, p0, z0){
    #   ll <- sum(unlist(lapply(1:nrow(data0),function(i){
    #     z0[i]*dmultinom(data0[i,],size=1, prob=p0[i,], log=TRUE)
    #   })))
    #   return(ll)
    # }
    # 
    # llfn3.slow <- function(w.arr,u, p, z){
    #   ll <- 0
    #   for(ii in 1:dim(w.arr)[1]){
    #     if(z[ii]==1){
    #       data <- c( c(w.arr[ii,,]), u[ii])
    #       ll <- ll + dmultinom(data, size = 1, prob = p[ii,], log=TRUE)
    #     }#else{
    #     #   ll <- ll + 0
    #     # }
    #   }#for
    #   return(ll)
    # }
    
    ##================================================
    ##-- The likelihood function for a single bug 
    ##================================================
    llfn.ind <- function(w.i,u.i, p.i, z.i){
      ll <- 0
      if(z.i==1){
        data <- c(c(w.i), u.i)
        ll <- dmultinom(data, size = 1, prob = p.i, log=TRUE)
      }
      return(ll)
    }
    
    ##================================================
    ##-- Some checks
    ##================================================
    # checks
    sum(w.arr)
    sum(y.mat)
    
    sum(u.aug)
    sum(u.aug) == M-sum(y.mat)
    
    sum(z)
    
    
    probs <- c(pmat)
    
    risk <- 1-exp(-(S.aug^nu)/TCmax) # M x 1
    
    p <- do.call(rbind, lapply(1:length(S.aug), function(i){ # M x (J*K + 1)
      c(risk[i] * probs / sum(probs), 1 - risk[i])
    }))
    all.equal(p, p.multi(pmat,S.aug)) # FALSE because the S is rearranged now because of augmentation
    
    apply(p, 1, sum) # p is row stochastic
    sum(apply(p, 1, sum) == 1)
    # data <- datafn(w.arr,u, z)
    data <- data.full.aug
    ll <- llfn(data, p, z)
    ll
    #llfn3.slow(w.arr.aug,u.aug, p, z)
    
    #llfn2.slow(data, p, z)
    
    ##================================================
    ##-- Some useful stats from observed data
    ##================================================
    
    detected <- apply(w.arr.aug, 1, sum) > 0
    n <- sum(detected)
    print(n)
    whichNotDetected <- which(detected == 0)
    whichDetected <- which(detected > 0)
    
    
    ## R code to fit the model to simulated data using the 2 MCMC algorithms
    
  }##--DO ALL
  
  ##==========================================================###
  ###=========== Draw posterior samples using MCMC ===========###
  ###=========== This version *updates* data w.arr ===========### 
  ##==========================================================###
  
  sizeBiased <- function(w.arr,y.mat,u,data,TC,
                         r,
                         lam0,S,psi,z,
                         J, K, M,varS, niters, tune=c(0.1, 0.1, 2,5,5e-6)){  #, xlims, ylims, #
    # monitorS=FALSE)
    
    # initial values
    
    
    
    ww.arr <- w.arr
    yy.mat <- y.mat
    uu <- u
    
    detected <- apply(ww.arr, 1, sum) > 0
    
    n <- sum(detected)
    n
    whichNotDetected <- which(detected == 0)
    whichDetected <- which(detected > 0)
    
    whichCellNotDetected <- which(yy.mat == 0, arr.ind = T)
    whichCellNotDetected
    whichCellNotDetected.key <- apply(whichCellNotDetected, 1, function(x){paste0(x[1],",",x[2])})
    whichCellNotDetected.key
    
    whichMissionNotDetected <- which(rowSums(yy.mat) == 0, arr.ind = T)
    whichMissionNotDetected
    
    # "2,1" %in% whichCellNotDetected.key
    #-------
    # checks
    detected1 <- uu == 0
    names(detected1) <- names(detected)
    all.equal(detected, detected1)
    #-------
    
    N <- sum(z)
    
    logit.a <- logit(a)
    logit.r <- logit(r)
    pmat <- 1 - exp(-TC)
    
    p <- p.multi(pmat, S)
    
    ll <-  llfn(data, p, z)
    
    S.names <- paste0("S[,", 1:M, "]")
    
    paramnames <- c("r",
                    "psi", "N", "ll")#
    out <- matrix(NA, nrow=niters, ncol=length(paramnames))
    colnames(out) <- paramnames
    
    # if(monitorS)
    Sout <- array(NA, c(M, niters))
    lam_out <- array(NA, c(M, niters))
    # total_size <- array(NA, c(M, niters))
    # detected_size <- array(NA, c(M, niters))
    SZ=array(NA, c(M, niters))
    SU=array(NA, c(M, niters))
    remaining_size <- array(0)
    
    
    
    batchsize <- 100 #500
    batch = 0
    min.eta = 0.01
    eta = min.eta
    
    # Setting counters
    naccept.b0 = naccept.b1 = naccept.b2 = naccept.lam0 = 0
    naccept.S = 0
    naccept.r = 0
    
    Sups <- 0
    zups <- 0
    
    
    
    cat("\ninitial values =", c(r,
                                psi, sum(z), ll), "\n\n")
    
    cat("Current time ", format(Sys.time(), "%H:%M:%S"), "\n\n")
    
    # c(S)
    iter = 1
    
    # niters <- 500
    
    for(iter in 1:niters) {
      
      if(iter %% 500 ==0) {
        cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
        cat("current =", out[iter-1,!colnames(out)%in%S.names], "\n")
        cat("  Acceptance rates\n")
        cat("    S =", naccept.S/batchsize, "\n")
      }
      
      if (iter%%batchsize == 0){
        if (1/sqrt(batch) < min.eta)  eta = 1/sqrt(batch)
        batch = batch + 1
      }
      
      ##--------------##
      # update lam0
      ##--------------##
      #lam0.cand <- rnorm(1, lam0, tune[4])
      # lam0.cand=rgamma(1,shape = 50,rate=0.5)
      # if(lam0.cand>0) {
        
        #priorS <- sum(dnbinom(S, size = lam0^2 / (varS - lam0), mu = lam0, log = TRUE))
        #priorS <- sum(dnbinom(S, size = lam0^2 / (varS - lam0), mu = lam0, log = TRUE))+dgamma(lam0,shape = 50,rate=0.5,log=TRUE)+
        #  dexp(lam0,rate =1/lam0.cand,log=TRUE)
        #priorScand <- sum(dnbinom(S, size = lam0.cand^2 / (varS - lam0.cand), mu = lam0.cand, log = TRUE))
        #priorScand <- sum(dnbinom(S, size = lam0.cand^2 / (varS - lam0.cand), mu = lam0.cand, log = TRUE))+dgamma(lam0.cand,shape = 50,rate=0.5,log=TRUE)+
        #  dexp(lam0.cand,rate =1/lam0,log=TRUE)
        
        #if(runif(1) < exp( priorScand - priorS ) ) {
        #  priorS <- priorScand
        #  lam0<-lam0.cand
        #  naccept.lam0 <- naccept.lam0 + 1
        #}
      #}#lam0.cand>0
      
      ##--------------##
      # update z
      ##--------------##
      
      ncprob <- p[,J*K+1] # M x 1
      fc<- ncprob*psi / (ncprob*psi + 1 - psi) # M x 1
      z[whichNotDetected] <- rbinom(M-n, 1, fc[whichNotDetected])
      z[whichDetected] <- 1
      N <- sum(z)
      
      
      ##--------------##
      # update psi
      ##--------------##
      psi <- rbeta(1, 1+N, 1+M-N)
      
      ##--------------##
      # update S
      ##--------------##
      
      for(ii in 1:M) {
        jump <- 3
        this.Scand <- ifelse(runif(1) < 0.5, S[ii]+jump, S[ii]-jump)
        condition <- this.Scand > 0
        
        if(!condition)
          next
        
        Scand <- S
        
        Scand[ii] <- this.Scand
        
        probs <- c(pmat)
        
        risk.cand.i <- 1-exp(-(Scand[ii]^nu)/TCmax)
        
        p.cand.i <- c(risk.cand.i * probs / sum(probs), 1 - risk.cand.i)
        llcand.i <- llfn.ind(ww.arr[ii,,],uu[ii], p.cand.i, z[ii])
        priorScand <- dnbinom(Scand[ii], size = ((lam0[ii])^2) / (varS - lam0[ii]), mu = lam0[ii], log = TRUE)
        #priorScand <- dpois(Scand[ii], lambda = lam0[ii], log = TRUE)
        
        risk.i <- 1-exp(-(S[ii]^nu)/TCmax)
        p.i <- c(risk.i * probs / sum(probs), 1 - risk.i)
        ll.i <- llfn.ind(ww.arr[ii,,],uu[ii], p.i, z[ii])
        priorS <- dnbinom(S[ii], size = ((lam0[ii])^2) / (varS - lam0[ii]), mu = lam0[ii], log = TRUE)
        #priorS <- dpois(S[ii], lambda = lam0[ii], log = TRUE)
        
        
        if(runif(1) < exp( (llcand.i+priorScand)  - (ll.i+priorS) ) ){
          
          S <- Scand
          naccept.S <- naccept.S + 1
        }
        
        this.lamcand=rgamma(1,shape = 50,rate=0.5)
        condition1 <- this.lamcand < varS
        
        if(!condition1)
          next
        lam0.cand=lam0
        lam0.cand[ii]=this.lamcand
        priorLamcand=dgamma(lam0.cand[ii],shape = 50,rate=0.5,log=TRUE)
        priorLam=dgamma(lam0[ii],shape = 50,rate=0.5,log=TRUE)
        priorSupdatecand=dnbinom(S[ii], size = ((lam0.cand[ii])^2) / (varS - lam0.cand[ii]), mu = lam0.cand[ii], log = TRUE)
        #priorSupdatecand=dpois(S[ii], lambda = lam0.cand[ii], log = TRUE)
        
        priorSupdate=dnbinom(S[ii], size = ((lam0[ii])^2) / (varS - lam0[ii]), mu = lam0[ii], log = TRUE)
        #priorSupdate=dpois(S[ii], lambda = lam0[ii], log = TRUE)
        
        if(runif(1) < exp( (llcand.i+priorSupdatecand+priorLamcand)  - (ll.i+priorSupdate+priorLam) ) ){
          
          lam0 <- lam0.cand
          
          naccept.lam0 <- naccept.lam0 + 1
        }
        
      }#ii
      
      
      # ##--------------##
      # # Log-likelihood
      # ##--------------##
      p <- p.multi(pmat, S)
      ll <- llfn(data, p, z)
      
      
      if (iter%%batchsize == 0){
        
        SigmaDiff = ifelse(naccept.b1 > 0.44*batchsize, exp(2*eta), exp(-2*eta))
        if(iter <= niters){ tune[2] = tune[2] * SigmaDiff}
        
        SigmaDiff = ifelse(naccept.b2 > 0.44*batchsize, exp(2*eta), exp(-2*eta))
        if(iter <= niters){ tune[3] = tune[3] * SigmaDiff}
        
        SigmaDiff = ifelse(naccept.lam0 > 0.44*batchsize, exp(2*eta), exp(-2*eta))
        if(iter <= niters){ tune[4] = tune[4] * SigmaDiff}
        
        SigmaDiff = ifelse(naccept.r > 0.44*batchsize, exp(2*eta), exp(-2*eta))
        if(iter <= niters){ tune[5] = tune[5] * SigmaDiff}
        
        naccept.b0 = 0   # reset counter for next batch
        naccept.b1 = 0   # reset counter for next batch
        naccept.b2 = 0   # reset counter for next batch
        naccept.lam0 = 0   # reset counter for next batch
        naccept.r = 0   # reset counter for next batch
        # naccept.a = 0   # reset counter for next batch
        naccept.S = 0   # reset counter for next batch
      }##if (iter%%batchsize == 0){
      
      out[iter,] <- c(r,
                      psi, N, ll)
      # if(monitorS)
      Sout[,iter] <- S
      lam_out[,iter] = lam0
      SZ[,iter]=Sout[,iter]*z
      SU[,iter]=Sout[,iter]*(1-u)
      
    }##iter
    
    total_size=apply(SZ,2,sum)
    print(length(total_size))
    detected_size=apply(SU,2,sum)
    remaining_size=total_size-detected_size
    
    list(out=out, Sout=Sout, lam_out=lam_out, total_size=total_size, detected_size=detected_size, remaining_size=remaining_size)
  }#sizeBiased
  
  ##==========================================================###
  ##==========================================================###
  ##==========================================================###
  
  {##--DO ALL
    
    start.time <- Sys.time()
    stm <- proc.time()
    
    niters <- 50000
    
    MCMCout1 <- sizeBiased(w.arr=w.arr.aug,y.mat=y.mat,u=u.aug,data=data.full.aug,TC=TC,
                           r=r,
                           lam0=lam0,S=S.aug,psi=psi,z=z,
                           J=J, K=K, M=M,varS=varS,
                           niters=niters,
                           tune=tune)
    
    MCMCout2 <- sizeBiased(w.arr=w.arr.aug,y.mat=y.mat,u=u.aug,data=data.full.aug,TC=TC,
                           r=r,
                           lam0=lam0,S=S.aug,psi=psi,z=z,
                           J=J, K=K, M=M,varS=varS,
                           niters=niters,
                           tune=tune)
    
    MCMCout3 <- sizeBiased(w.arr=w.arr.aug,y.mat=y.mat,u=u.aug,data=data.full.aug,TC=TC,
                           r=r,
                           lam0=lam0,S=S.aug,psi=psi,z=z,
                           J=J, K=K, M=M,varS=varS,
                           niters=niters,
                           tune=tune)
    
    RunTime <- proc.time()-stm 
    
    end.time <- Sys.time()
    print(end.time-start.time)
    
    
  }##--DO ALL
  
  {##--DO ALL
    
    
    burnin <- niters/2
    
    summary(MCMCout1$out)
    
    j <- 4
    PAR.ESTIMATES1 <- do.call(rbind, lapply(1:ncol(MCMCout1$out), function(j){
      x <- MCMCout1$out[,j]
      simvalue<-get(colnames(MCMCout1$out)[j])
      out <- c((mean(x)-simvalue)/simvalue*100, sd(x)/mean(x)*100, mean(x), sd(x), quantile(x, 0.025), quantile(x, 0.5), quantile(x, 0.975))
    }))
    dimnames(PAR.ESTIMATES1) <- list(colnames(MCMCout1$out), c("RB%", "CV%", "Mean", "SD", "2.5%", "50%", "97.5%"))
    print(PAR.ESTIMATES1)
    # Histories
    
    summary(MCMCout2$out)
    
    j <- 4
    PAR.ESTIMATES2 <- do.call(rbind, lapply(1:ncol(MCMCout2$out), function(j){
      x <- MCMCout2$out[,j]
      simvalue<-get(colnames(MCMCout2$out)[j])
      out <- c((mean(x)-simvalue)/simvalue*100, sd(x)/mean(x)*100, mean(x), sd(x), quantile(x, 0.025), quantile(x, 0.5), quantile(x, 0.975))
    }))
    dimnames(PAR.ESTIMATES2) <- list(colnames(MCMCout2$out), c("RB%", "CV%", "Mean", "SD", "2.5%", "50%", "97.5%"))
    print(PAR.ESTIMATES2)
    
    summary(MCMCout3$out)
    
    j <- 4
    PAR.ESTIMATES3 <- do.call(rbind, lapply(1:ncol(MCMCout3$out), function(j){
      x <- MCMCout3$out[,j]
      simvalue<-get(colnames(MCMCout3$out)[j])
      out <- c((mean(x)-simvalue)/simvalue*100, sd(x)/mean(x)*100, mean(x), sd(x), quantile(x, 0.025), quantile(x, 0.5), quantile(x, 0.975))
    }))
    dimnames(PAR.ESTIMATES3) <- list(colnames(MCMCout3$out), c("RB%", "CV%", "Mean", "SD", "2.5%", "50%", "97.5%"))
    print(PAR.ESTIMATES3)
    
    
    plot.mcmc <- function(samples1,samples2,samples3,
                          par,
                          sim.values = NULL,
                          niters=50000,
                          burnin=25000,
                          if.hist=F){
      # par(mfrow=c(1,2))
      for(this.par in par){
        sim.value <- get(this.par)
        # Traceplots
        #plot(samples[1:niters,this.par], type="l", ylab=this.par)
        plot(samples1[1:niters,this.par], type="l", ylab=this.par,col="red",main="Traceplot")
        lines(samples2[1:niters,this.par],col=rgb(0, 1, 0, 0.5))
        lines(samples3[1:niters,this.par],col=rgb(0, 0, 1, 0.5))
        abline(h = sim.value, lwd = 2, col = 'red')
        
        # Densities
        
        if(if.hist==F) 
        { 
          plot(density(c(samples1[(burnin+1):niters,this.par],samples2[(burnin+1):niters,this.par],samples3[(burnin+1):niters,this.par])),
                       xlab=this.par,main="density plot")
          
        } else if(if.hist==T)
        # else if(if.hist==T)
        # {
        #   hist(samples1[(burnin+1):niters,this.par], xlab=this.par, freq=FALSE, main="",col="red")
        #   hist(samples2[(burnin+1):niters,this.par],freq=FALSE, col=rgb(0, 1, 0, 0.5),add=TRUE)
        #   hist(samples3[(burnin+1):niters,this.par],freq=FALSE, col=rgb(0, 0, 1, 0.5),add=TRUE)
        # }
        
        
        {
          plot(table(samples1[(burnin+1):niters,this.par])/sum(table(samples1[(burnin+1):niters,this.par])), type="h", xlab=this.par, ylab="relative frequency", main="Relative frequency distribution",col="red")
          lines(table(samples2[(burnin+1):niters,this.par])/sum(table(samples2[(burnin+1):niters,this.par])),type="h", col=rgb(0, 1, 0, 0.5))
          lines(table(samples3[(burnin+1):niters,this.par])/sum(table(samples3[(burnin+1):niters,this.par])),type="h", col=rgb(0, 0, 1, 0.5))
        }
        
        abline(v = sim.value, lwd = 2, col = 'red')
      }#for
    }#plot.mcmc
    
    ts = format(Sys.time(), "%y%m%d_%H%M%S")
    modelName <- paste0("Results_Multinomial_Pallabi")
    if(!dir.exists(file.path(WD))){dir.create(file.path(WD))}
    if(!dir.exists(file.path(WD,modelName))){dir.create(file.path(WD, modelName))}
    
    outfolder <- file.path(WD,modelName,paste0("Results_Multinomial_Snegbinom_SperBug_", ts))
    if(!dir.exists(outfolder)){dir.create(outfolder,recursive=TRUE)}
    
    save(MCMCout1, MCMCout2,MCMCout3,PAR.ESTIMATES1,PAR.ESTIMATES2,PAR.ESTIMATES3,
         w.arr,y.mat,u,TC,
         r,
         lam0,S,psi,z,varS, 
         J, K, M, ll,
         niters, burnin,RunTime,
         tune, file = file.path(outfolder, paste0("mcmcsamples_",ts, ".RData"))
    )
    
    graphics.off()
    # ts = format(Sys.time(), '%d%m%y_%H%M%S')
    path <- file.path(outfolder, paste("Traceplots_", ts, ".pdf", sep = ""))
    pdf(file=path, width = 8, height = 8)
    
    par(mfrow= c(1, 1))
    parnames <- c('psi') 
    sim.values <- c(psi) 
    names(sim.values) <- parnames
    plot.mcmc(samples1=MCMCout1$out,samples2=MCMCout2$out,samples3=MCMCout3$out,par=parnames,niters=niters,burnin=burnin,if.hist=F)
    plot.mcmc(samples1=MCMCout1$out,samples2=MCMCout2$out,samples3=MCMCout3$out,par=c('N'),niters=niters,burnin=burnin,if.hist=T)
    
    for(i in c(1:3, (M-2):M)){
      
      plot(MCMCout1$lam_out[i,1:niters], type="l", ylab=paste0('lam0[', i, "]"),main="Traceplot", col="red")
      lines(MCMCout2$lam_out[i,1:niters], type="l",col=rgb(0, 1, 0, 0.5))
      lines(MCMCout3$lam_out[i,1:niters], type="l",col=rgb(0, 0, 1, 0.5))
      #abline(h = S.aug[i], lwd = 2, col = 'red')
      hist(MCMCout1$lam_out[i,1:niters], xlab=paste0('lam0[', i, "]"), freq=FALSE, main="Relative frequency distribution",col="red")
      hist(MCMCout2$lam_out[i,1:niters], freq=FALSE,add=TRUE,col=rgb(0, 1, 0, 0.5))
      hist(MCMCout3$lam_out[i,1:niters], freq=FALSE,add=TRUE,col=rgb(0, 0, 1, 0.5))
      #abline(v = S.aug[i], lwd = 2, col = 'red')
    }
    
    
    for(i in c(1:3, (M-2):M)){
      
      plot(MCMCout1$Sout[i,1:niters], type="l", ylab=paste0('S[', i, "]"),main="Traceplot",col="red")
      lines(MCMCout2$Sout[i,1:niters], type="l",col=rgb(0, 1, 0, 0.5))
      lines(MCMCout3$Sout[i,1:niters], type="l",col=rgb(0, 0, 1, 0.5))
      #abline(h = S.aug[i], lwd = 2, col = 'red')
      # hist(MCMCout1$Sout[i,1:niters], xlab=paste0('S[', i, "]"), freq=FALSE, main="",col="red")
      # hist(MCMCout2$Sout[i,1:niters], freq=FALSE,add=TRUE,col=rgb(0, 1, 0, 0.5))
      # hist(MCMCout3$Sout[i,1:niters], freq=FALSE,add=TRUE,col=rgb(0, 0, 1, 0.5))
      plot(table(MCMCout1$Sout[i,(burnin+1):niters])/sum(table(MCMCout1$Sout[i,(burnin+1):niters])), type="h", xlab=paste0('S[', i, "]"), ylab="relative frequency", main="Relative frequency distribution",col="red")
      lines(table(MCMCout2$Sout[i,(burnin+1):niters])/sum(table(MCMCout2$Sout[i,(burnin+1):niters])),type="h",col=rgb(0, 1, 0, 0.5))
      lines(table(MCMCout3$Sout[i,(burnin+1):niters])/sum(table(MCMCout3$Sout[i,(burnin+1):niters])),type="h", col=rgb(0, 0, 1, 0.5))
      #abline(v = S.aug[i], lwd = 2, col = 'red')
      
      # plot(MCMCout1$Sout[i,1:niters], type="l", ylab=paste0('S[', i, "]"))
      # abline(h = S.aug[i], lwd = 2, col = 'red')
      # hist(MCMCout1$Sout[i,1:niters], xlab=paste0('S[', i, "]"), freq=FALSE, main="")
      # abline(v = S.aug[i], lwd = 2, col = 'red')
    }
    # }
    
    
    graphics.off()
    
    
  }##--DO ALL
  
}##--DO ALL



g1=apply(MCMCout1$out[(burnin+1):niters,],2,geweke.diag,frac1=0.25,frac2=0.25)
g2=apply(MCMCout2$out[(burnin+1):niters,],2,geweke.diag,frac1=0.25,frac2=0.25)
g3=apply(MCMCout3$out[(burnin+1):niters,],2,geweke.diag,frac1=0.25,frac2=0.25)

gel_diag=list()
ess=0
mc_se_ch1=list()
mc_se_ch2=list()
mc_se_ch3=list()
for(par in 1:4)
{
  MCMCout_out=list(MCMCout1$out[(burnin+1):niters,par],MCMCout2$out[(burnin+1):niters,par],MCMCout3$out[(burnin+1):niters,par])
  MCMC_ch=as.mcmc.list(lapply(MCMCout_out,mcmc))
  
  mc_se_ch1[[par]]=mcse(x = MCMC_ch[[1]], method = "bm", g = NULL)
  mc_se_ch2[[par]]=mcse(x = MCMC_ch[[2]], method = "bm", g = NULL)
  mc_se_ch3[[par]]=mcse(x = MCMC_ch[[3]], method = "bm", g = NULL)
  
  gel_diag[[par]]=gelman.diag(MCMC_ch)
  ess[par]=effectiveSize(MCMC_ch)
}

g1_lam=apply(MCMCout1$lam_out[(M-2):M,(burnin+1):niters],1,geweke.diag,frac1=0.25,frac2=0.25)
g2_lam=apply(MCMCout2$lam_out[(M-2):M,(burnin+1):niters],1,geweke.diag,frac1=0.25,frac2=0.25)
g3_lam=apply(MCMCout3$lam_out[(M-2):M,(burnin+1):niters],1,geweke.diag,frac1=0.25,frac2=0.25)

gel_diag_lam=list()
ess_lam=0
for(i in c((M-2):M)){
MCMCout_lamout=list(MCMCout1$lam_out[i,(burnin+1):niters],MCMCout2$lam_out[i,(burnin+1):niters],MCMCout3$lam_out[i,(burnin+1):niters])
MCMC_ch_lam=as.mcmc.list(lapply(MCMCout_lamout,mcmc))
gel_diag_lam[[i]]=gelman.diag(MCMC_ch_lam)
ess_lam[i]=effectiveSize(MCMC_ch_lam)
}


g1_S=apply(MCMCout1$Sout[(M-2):M,(burnin+1):niters],1,geweke.diag,frac1=0.25,frac2=0.25)
g2_S=apply(MCMCout2$Sout[(M-2):M,(burnin+1):niters],1,geweke.diag,frac1=0.25,frac2=0.25)
g3_S=apply(MCMCout3$Sout[(M-2):M,(burnin+1):niters],1,geweke.diag,frac1=0.25,frac2=0.25)

gel_diag_S=list()
ess_S=0
for(i in c((M-2):M)){
  MCMCout_Sout=list(MCMCout1$Sout[i,(burnin+1):niters],MCMCout2$Sout[i,(burnin+1):niters],MCMCout3$Sout[i,(burnin+1):niters])
  MCMC_ch_S=as.mcmc.list(lapply(MCMCout_Sout,mcmc))
  gel_diag_S[[i]]=gelman.diag(MCMC_ch_S)
  ess_S[i]=effectiveSize(MCMC_ch_S)
}



save(g1,g2,g3,g1_lam,g2_lam,g3_lam,g1_S,g2_S,g3_S,gel_diag,gel_diag_S,ess,ess_S,mc_se_ch1,mc_se_ch2,mc_se_ch3,
     file = file.path(outfolder, paste0("diag",ts,".RData")))


eps=100
reliability1=mean(MCMCout1$remaining_size<eps);reliability1
reliability2=mean(MCMCout2$remaining_size<eps);reliability2
reliability3=mean(MCMCout3$remaining_size<eps);reliability3

mean_chains=data.frame("psi"=c(mean(MCMCout1$out[(burnin+1):niters,2]),
  mean(MCMCout2$out[(burnin+1):niters,2]),mean(MCMCout3$out[(burnin+1):niters,2])),"N"=c(mean(MCMCout1$out[(burnin+1):niters,3]),
  mean(MCMCout2$out[(burnin+1):niters,3]),mean(MCMCout3$out[(burnin+1):niters,3])),"lam1"=c(mean(MCMCout1$lam_out[1,(burnin+1):niters]),
  mean(MCMCout2$lam_out[1,(burnin+1):niters]),mean(MCMCout3$lam_out[1,(burnin+1):niters])),"lam2"=
  c(mean(MCMCout1$lam_out[2,(burnin+1):niters]),mean(MCMCout2$lam_out[1,(burnin+1):niters]),
  mean(MCMCout3$lam_out[1,(burnin+1):niters])),"lam399"=c(mean(MCMCout1$lam_out[399,(burnin+1):niters]),
  mean(MCMCout2$lam_out[399,(burnin+1):niters]),mean(MCMCout3$lam_out[399,(burnin+1):niters])),
  "lam400"=c(mean(MCMCout1$lam_out[400,(burnin+1):niters]),mean(MCMCout2$lam_out[400,(burnin+1):niters]),
  mean(MCMCout3$lam_out[400,(burnin+1):niters])),"S1"=c(mean(MCMCout1$Sout[1,(burnin+1):niters]),
  mean(MCMCout2$Sout[1,(burnin+1):niters]),mean(MCMCout3$Sout[1,(burnin+1):niters])), 
  "S2"=c(mean(MCMCout1$Sout[2,(burnin+1):niters]), mean(MCMCout2$Sout[2,(burnin+1):niters]),
  mean(MCMCout3$Sout[2,(burnin+1):niters])),"S399"=c(mean(MCMCout1$Sout[399,(burnin+1):niters]),
  mean(MCMCout2$Sout[399,(burnin+1):niters]),mean(MCMCout3$Sout[399,(burnin+1):niters])), "S400"=
  c(mean(MCMCout1$Sout[400,(burnin+1):niters]), mean(MCMCout2$Sout[400,(burnin+1):niters]),mean(MCMCout3$Sout[400,(burnin+1):niters])))

sd_chains=data.frame("psi"=c(sd(MCMCout1$out[(burnin+1):niters,2]),
    sd(MCMCout2$out[(burnin+1):niters,2]),sd(MCMCout3$out[(burnin+1):niters,2])),
    "N"=c(sd(MCMCout1$out[(burnin+1):niters,3]), sd(MCMCout2$out[(burnin+1):niters,3]),
    sd(MCMCout3$out[(burnin+1):niters,3])),"lam1"=c(sd(MCMCout1$lam_out[1,(burnin+1):niters]),
    sd(MCMCout2$lam_out[1,(burnin+1):niters]),sd(MCMCout3$lam_out[1,(burnin+1):niters])),"lam2"=
    c(sd(MCMCout1$lam_out[2,(burnin+1):niters]),sd(MCMCout2$lam_out[1,(burnin+1):niters]),
    sd(MCMCout3$lam_out[1,(burnin+1):niters])),"lam399"=c(sd(MCMCout1$lam_out[399,(burnin+1):niters]),
    sd(MCMCout2$lam_out[399,(burnin+1):niters]),sd(MCMCout3$lam_out[399,(burnin+1):niters])),
    "lam400"=c(sd(MCMCout1$lam_out[400,(burnin+1):niters]),sd(MCMCout2$lam_out[400,(burnin+1):niters]),
    sd(MCMCout3$lam_out[400,(burnin+1):niters])),"S1"=c(sd(MCMCout1$Sout[1,(burnin+1):niters]),
    sd(MCMCout2$Sout[1,(burnin+1):niters]),sd(MCMCout3$Sout[1,(burnin+1):niters])), 
    "S2"=c(sd(MCMCout1$Sout[2,(burnin+1):niters]), sd(MCMCout2$Sout[2,(burnin+1):niters]),
    sd(MCMCout3$Sout[2,(burnin+1):niters])),"S399"=c(sd(MCMCout1$Sout[399,(burnin+1):niters]),
    sd(MCMCout2$Sout[399,(burnin+1):niters]),sd(MCMCout3$Sout[399,(burnin+1):niters])), "S400"=
    c(sd(MCMCout1$Sout[400,(burnin+1):niters]), sd(MCMCout2$Sout[400,(burnin+1):niters]),sd(MCMCout3$Sout[400,(burnin+1):niters])))

cv_chains=data.frame("psi"=c(cv(MCMCout1$out[(burnin+1):niters,2]),
                             cv(MCMCout2$out[(burnin+1):niters,2]),cv(MCMCout3$out[(burnin+1):niters,2])),
                     "N"=c(cv(MCMCout1$out[(burnin+1):niters,3]), cv(MCMCout2$out[(burnin+1):niters,3]),
                           cv(MCMCout3$out[(burnin+1):niters,3])),"lam1"=c(cv(MCMCout1$lam_out[1,(burnin+1):niters]),
                          cv(MCMCout2$lam_out[1,(burnin+1):niters]),cv(MCMCout3$lam_out[1,(burnin+1):niters])),"lam2"=
                       c(cv(MCMCout1$lam_out[2,(burnin+1):niters]),cv(MCMCout2$lam_out[1,(burnin+1):niters]),
                         cv(MCMCout3$lam_out[1,(burnin+1):niters])),"lam399"=c(cv(MCMCout1$lam_out[399,(burnin+1):niters]),
                        cv(MCMCout2$lam_out[399,(burnin+1):niters]),cv(MCMCout3$lam_out[399,(burnin+1):niters])),
                     "lam400"=c(cv(MCMCout1$lam_out[400,(burnin+1):niters]),cv(MCMCout2$lam_out[400,(burnin+1):niters]),
                        cv(MCMCout3$lam_out[400,(burnin+1):niters])),"S1"=c(cv(MCMCout1$Sout[1,(burnin+1):niters]),
                      cv(MCMCout2$Sout[1,(burnin+1):niters]),cv(MCMCout3$Sout[1,(burnin+1):niters])), 
                     "S2"=c(cv(MCMCout1$Sout[2,(burnin+1):niters]), cv(MCMCout2$Sout[2,(burnin+1):niters]),
                            cv(MCMCout3$Sout[2,(burnin+1):niters])),"S399"=c(cv(MCMCout1$Sout[399,(burnin+1):niters]),
                            cv(MCMCout2$Sout[399,(burnin+1):niters]),cv(MCMCout3$Sout[399,(burnin+1):niters])), "S400"=
                       c(cv(MCMCout1$Sout[400,(burnin+1):niters]), cv(MCMCout2$Sout[400,(burnin+1):niters]),cv(MCMCout3$Sout[400,(burnin+1):niters])))

t(round(mean_chains,4))
t(round(sd_chains,4))
100*t(round(cv_chains,4)) 

#d1=load("D:/1. PALLABI/PHD/ISRO/Program/Model fitting/PG_new/Final/Results/Results_Multinomial_Pallabi/nu_1_S_NB/mcmcsamples_231129_082528.RData")

