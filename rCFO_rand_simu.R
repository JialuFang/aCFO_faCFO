if (!require(CFO, quietly = TRUE)) {
  install.packages("CFO")
} else {
  cat("CFO package is already installed.\n")
}
library(CFO)
library(parallel)
source("utilities.R")
source("./setting/Sce_6dose_phi20_L.R")

CFO.rand.next <- function(target, cys, cns, currdose, prior.para=list(alp.prior=target, bet.prior=1-target),
                          cutoff.eli=0.95, early.stop=0.95){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  overdose.fn <- function(phi, threshold, prior.para=list()){
    y <- prior.para$y
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    # print(data.frame("prob of overdose" = pp))
    if ((pp >= threshold) & (prior.para$n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- bet.prior + n1 - y1
    bet2 <- bet.prior + n2 - y2
    fn.min <- function(x){
      dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2)) 
    }
    fn.max <- function(x){
      pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
    }
    const.min <- integrate(fn.min, lower=0, upper=0.99, subdivisions=1000, rel.tol = 1e-10)$value
    const.max <- integrate(fn.max, lower=0, upper=1, rel.tol = 1e-10)$value
    p1 <- integrate(fn.min, lower=0, upper=phi)$value/const.min
    p2 <- integrate(fn.max, lower=0, upper=phi)$value/const.max
    
    list(p1=p1, p2=p2)
  }
  
  OR.values <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior, type){
    ps <- prob.int(phi, y1, n1, y2, n2, alp.prior, bet.prior)
    if (type=="L"){
      pC <- 1 - ps$p2
      pL <- 1 - ps$p1
      oddsC <- pC/(1-pC)
      oddsL <- pL/(1-pL)
      OR <- oddsC*oddsL
      
    }else if (type=="R"){
      pC <- 1 - ps$p1
      pR <- 1 - ps$p2
      oddsC <- pC/(1-pC)
      oddsR <- pR/(1-pR)
      OR <- (1/oddsC)/oddsR
    }
    return(OR)
  }
  
  All.OR.table <- function(phi, n1, n2, type, alp.prior, bet.prior){
    ret.mat <- matrix(rep(0, (n1+1)*(n2+1)), nrow=n1+1)
    for (y1cur in 0:n1){
      for (y2cur in 0:n2){
        ret.mat[y1cur+1, y2cur+1] <- OR.values(phi, y1cur, n1, y2cur, n2, alp.prior, bet.prior, type)
      }
    }
    ret.mat
  }
  
  # compute the marginal prob when lower < phiL/phiC/phiR < upper
  # i.e., Pr(Y=y|lower<phi<upper)
  margin.phi <- function(y, n, lower, upper){
    C <- 1/(upper-lower)
    fn <- function(phi) {
      dbinom(y, n, phi)*C
    }
    integrate(fn, lower=lower, upper=upper)$value
  }
  
  # Obtain the table of marginal distribution of (y1, y2) 
  # after intergrate out (phi1, phi2)
  # under H0 and H1
  # H0: phi1=phi, phi < phi2 < 2phi
  # H1: phi2=phi, 0   < phi1 < phi
  margin.ys.table <- function(n1, n2, phi, hyperthesis){
    if (hyperthesis=="H0"){
      p.y1s <- dbinom(0:n1, n1, phi)
      p.y2s <- sapply(0:n2, margin.phi, n=n2, lower=phi, upper=2*phi)
    }else if (hyperthesis=="H1"){
      p.y1s <- sapply(0:n1, margin.phi, n=n1, lower=0, upper=phi)
      p.y2s <- dbinom(0:n2, n2, phi)
    }
    p.y1s.mat <- matrix(rep(p.y1s, n2+1), nrow=n1+1)
    p.y2s.mat <- matrix(rep(p.y2s, n1+1), nrow=n1+1, byrow=TRUE)
    margin.ys <- p.y1s.mat * p.y2s.mat
    margin.ys
  }
  
  # Obtain the optimal gamma for the hypothesis test
  optim.gamma.fn <- function(n1, n2, phi, type, alp.prior, bet.prior){
    OR.table <- All.OR.table(phi, n1, n2, type, alp.prior, bet.prior) 
    ys.table.H0 <- margin.ys.table(n1, n2, phi, "H0")
    ys.table.H1 <- margin.ys.table(n1, n2, phi, "H1")
    
    argidx <- order(OR.table)
    sort.OR.table <- OR.table[argidx]
    sort.ys.table.H0 <- ys.table.H0[argidx]
    sort.ys.table.H1 <- ys.table.H1[argidx]
    n.tol <- length(sort.OR.table)
    
    if (type=="L"){
      errs <- rep(0, n.tol-1)
      for (i in 1:(n.tol-1)){
        err1 <- sum(sort.ys.table.H0[1:i])
        err2 <- sum(sort.ys.table.H1[(i+1):n.tol])
        err <- err1 + err2
        errs[i] <- err
      }
      min.err <- min(errs)
      if (min.err > 1){
        gam <- 0
        min.err <- 1
      }else {
        minidx <- which.min(errs)
        gam <- sort.OR.table[minidx]
      }
    }else if (type=='R'){
      errs <- rep(0, n.tol-1)
      for (i in 1:(n.tol-1)){
        err1 <- sum(sort.ys.table.H1[1:i])
        err2 <- sum(sort.ys.table.H0[(i+1):n.tol])
        err <- err1 + err2
        errs[i] <- err
      }
      min.err <- min(errs)
      if (min.err > 1){
        gam <- 0
        min.err <- 1
      }else {
        minidx <- which.min(errs)
        gam <- sort.OR.table[minidx]
      }
      
    }
    list(gamma=gam, min.err=min.err)
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  cover.doses <- c(0,0,0)
  
  for (i in 1:3){
    cy <- cys[i]
    cn <- cns[i]
    if (is.na(cn)){
      cover.doses[i] <- NA
    }else{
      prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior))
      if (overdose.fn(target, cutoff.eli, prior.para)){
        cover.doses[i:3] <- 1
        break()
      }
    }
  }
  
  if (cutoff.eli != early.stop) {
    cy <- cys[1]
    cn <- cns[1]
    if (is.na(cn)){
      cover.doses[1] <- NA
    }else{
      prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior))
      if (overdose.fn(target, early.stop, prior.para)){
        cover.doses[1:3] <- 1
      }
    }
  }
  
  cover.doses <- ifelse(is.na(cys), NA, cover.doses)
  
  position <- which(cover.doses == 1)[1]
  overtox <- c(-1, 0, 1)[position] + currdose
  prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior))
  if ((cover.doses[2] == 1)&(currdose == 1)){
    index <- NA
    decision <- "stop"
  } else {
    if (cover.doses[2] == 1){
      index <- -1
      decision <- "de-escalation"
    }
    else{
      if (is.na(cys[1]) & (cover.doses[3]==1)){
        index <- 0
        decision <- "stay"
      }
      else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
        gam2 <- optim.gamma.fn(cns[2], cns[3], target, "R", alp.prior, bet.prior)$gamma 
        OR.v2 <- OR.values(target, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        if (OR.v2>gam2){
          index <- 1
          decision <- "escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (is.na(cys[3]) | (cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], target, "L", alp.prior, bet.prior)$gamma 
        OR.v1 <- OR.values(target, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        if (OR.v1>gam1){
          index <- -1
          decision <- "de-escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], target, "L", alp.prior, bet.prior)$gamma 
        gam2 <- optim.gamma.fn(cns[2], cns[3], target, "R", alp.prior, bet.prior)$gamma 
        OR.v1 <- OR.values(target, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        OR.v2 <- OR.values(target, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        pL <- OR.v1/(OR.v1+OR.v2)
        v1 <- OR.v1 > gam1
        v2 <- OR.v2 > gam2
        if (v1 & !v2){
          if (runif(1) < pL){
            index <- -1
            decision <- "de-escalation"
          } else {
            index <- 0
            decision <- "stay"
          }
        }else if (!v1 & v2){
          if (runif(1) < pL){
            index <- 0
            decision <- "stay"
          } else {
            index <- 1
            decision <- "escalation"
          }
        }else if (v1 & v2){
          if (OR.v1 == OR.v2){
            print('equal occur')
            index <- 0
            decision <- "stay"
          } else {
            if (runif(1) < pL){
              index <- -1
              decision <- "de-escalation"
            } else {
              index <- 1
              decision <- "escalation"
            }
          }
        }else{
          index <- 0
          decision <- "stay"
        }
      }
    }
  }
  
  if (decision=='stop'){
    nextdose <- 99
  }else{
    nextdose <- currdose+index
  }
  
  out <- list(target=target, cys=cys, cns=cns, decision=decision, currdose = currdose, 
              nextdose=nextdose, overtox=overtox)
  class(out) <- c("cfo_decision", "cfo")
  return(out)
}

CFO.rand.simu <- function(design, target, p.true, init.level=1, ncohort, cohortsize,
                          prior.para=list(alp.prior=target, bet.prior=1-target),
                          cutoff.eli=0.95, early.stop=0.95, seed=NULL){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.9){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  overdose.fn <- function(phi, threshold, prior.para=list()){
    y <- prior.para$y
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    # print(data.frame("prob of overdose" = pp))
    if ((pp >= threshold) & (prior.para$n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  set.seed(seed)
  earlystop <- 0
  ndose <- length(p.true)
  doselist <- rep(0, ncohort)
  currdose <- init.level
  
  ays <- rep(0, ndose) # number of responses for different doses.
  ans <- rep(0, ndose) # number of subject for different doses.
  tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
  DLTlist <- c()
  
  for (i in 1:ncohort){
    pc <- p.true[currdose]
    doselist[i] <- currdose
    
    # sample from current dose
    cres <- rbinom(cohortsize, 1, pc)
    DLTlist <- c(DLTlist, cres)
    
    # update results
    ays[currdose] <- ays[currdose] + sum(cres)
    ans[currdose] <- ans[currdose] + cohortsize
    
    cy <- ays[currdose]
    cn <- ans[currdose]
    
    prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
    
    if (overdose.fn(target, cutoff.eli, prior.para)){
      tover.doses[currdose:ndose] <- 1
    }
    
    if (currdose == 1){
      if (cutoff.eli != early.stop) {
        cy <- ays[1]
        cn <- ans[1]
        prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
        if (overdose.fn(target, early.stop, prior.para)){
          tover.doses[1:ndose] <- 1
        }
      }
    }
    
    
    if (tover.doses[1] == 1){
      earlystop <- 1
      break()
    }
    
    prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior))
    if (design == 'CFO'){
      # the results for current 3 dose levels
      if (currdose!=1){
        cys <- ays[(currdose-1):(currdose+1)]
        cns <- ans[(currdose-1):(currdose+1)]
      }else{
        cys <- c(NA, ays[1:(currdose+1)])
        cns <- c(NA, ans[1:(currdose+1)])
      }
      currdose <- CFO.rand.next(target, cys, cns, currdose, prior.para, cutoff.eli, early.stop)$nextdose
    }else if (design == 'aCFO'){
      currdose <- aCFO.next(target, ays, ans, currdose, prior.para, cutoff.eli, early.stop)$nextdose
    }else{
      stop("The input design is invalid; it can only be set as 'CFO' or 'aCFO'.")
    }
  }
  
  if (earlystop==0){
    MTD <- CFO.selectmtd(target, ans, ays, prior.para, cutoff.eli, early.stop, verbose=FALSE)$MTD
  }else{
    MTD <- 99
  }
  
  correct <- 0
  if (MTD == target){
    correct <- 1
  }
  
  npercent <- ans[which(p.true == target)]/(ncohort*cohortsize)
  ptoxic <- sum(ans[which(p.true > target)])/(ncohort*cohortsize)
  
  out<-list(target=target, MTD=MTD, correct=correct, npatients=ans, ntox=ays, npercent=npercent, 
            over.doses=tover.doses, cohortdose=doselist, ptoxic=ptoxic, patientDLT=DLTlist,
            sumDLT=sum(DLTlist), earlystop=earlystop)
  class(out) <- c("cfo_trial", "cfo")
  return(out)
}

run.fn <- function(i){
  if (i%%100 == 0){
    print(i)
  }
  
  p.true.all <- gen.rand.doses(ndose, target, mu1=mu, mu2=mu)
  p.true <- p.true.all$p.true
  tmtd <- p.true.all$mtd.level
  
  randCFO <- CFO.rand.simu('CFO', target, p.true, init.level, ncohort, cohortsize, prior.para, 
                           cutoff.eli=0.95, early.stop=0.95, seed = seeds[i])
  randCFO.res <- list(MTD=randCFO$MTD, dose.ns=randCFO$npatients, DLT.ns=randCFO$ntox, 
                      p.true=p.true, target=target, over.doses=randCFO$over.doses)
  CFO <- CFO.simu('CFO', target, p.true, init.level, ncohort, cohortsize, prior.para, 
                  cutoff.eli=0.95, early.stop=0.95, seed = seeds[i])
  CFO.res <- list(MTD=CFO$MTD, dose.ns=CFO$npatients, DLT.ns=CFO$ntox, 
                  p.true=p.true, target=target, over.doses=CFO$over.doses)
  ress <- list(
    cfo=CFO.res,
    randcfo=randCFO.res,
    paras=list(p.true=p.true, 
               mtd=tmtd, 
               prior.para=prior.para,
               target=target,
               ncohort=ncohort,
               cohortsize=cohortsize)
  )
  ress
}

diff_list <- c(0.05, 0.07, 0.1, 0.15)
mu_list <- c(0.344, 0.484, 0.637, 0.852)
ndose <- length(p.trues[[1]])
init.level <- 1

for (i in 2:4){
  mu <- mu_list[i]
  diff <- diff_list[i]
  prior.para <- list(alp.prior=target, bet.prior=1-target)
  
  nsimu <- 5000
  seeds <- 1:nsimu
  
  file.name <- paste0("./test_data/","rand_MTDSimu_",ndose, "dose_phi_", target*100, "_random_",  diff*100, "_S.Rdata")
  
  t <- system.time({
    results <- mclapply(1:nsimu, run.fn, mc.cores=40)
  })
  print(t)
  
  print(post.process.random(results))
  save(results, file=file.name)
  
  cfo.ress <- lapply(1:nsimu, function(i)results[[i]]$cfo)
  randcfo.ress <- lapply(1:nsimu, function(i)results[[i]]$randcfo)
  sum.all <- list(
    CFO = phase1.post.fn(cfo.ress),
    randCFO = phase1.post.fn(randcfo.ress)
  )
  print(phase.I.pretty.tb(sum.all))
}