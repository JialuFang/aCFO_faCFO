library(magrittr)
library(parallel)
library(BOIN)

source("CRM_tox_utils.R")
source("BOIN_tox_utils.R")
source("CFO_tox_utils.R")
source("ACFO_tox_utils.R")
source("utilities.R")
source("./setting/Sce_7dose_phi20.R")

run.fn <- function(i){
  set.seed(seeds[i])
  if (i%%100 == 0){
    print(i)
  }
  cfo.res <- CFO.simu.fn(target, p.true, ncohort=ncohort, init.level=init.level,  
                     cohortsize=cohortsize, add.args=add.args)
  crm.res <- crm.simu.fn(target=target, p.true=p.true, 
                         init.level=init.level, cohortsize=cohortsize, ncohort=ncohort)
  boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, 
                           init.level=init.level, cohortsize=cohortsize)
  acfo.res <- ACFO.simu.fn(target, p.true, ncohort=ncohort, init.level=init.level,
                          cohortsize=cohortsize, add.args=add.args)
  ress <- list(
    cfo=cfo.res,
    crm=crm.res,
    boin=boin.res,
    acfo=acfo.res, 
    paras=list(p.true=p.true, 
               mtd=tmtd, 
               add.args=add.args,
               target=target,
               ncohort=ncohort,
               cohortsize=cohortsize)
  )
  ress
  
}

for (i in 1:length(p.trues)){
  init.level <- 1
  nsimu <- 5000
  
  add.args <- list(alp.prior=target, bet.prior=1-target)
  
  idx <- i
  p.true <- p.trues[[idx]]
  ndose <- length(p.true)
  tmtd <- MTD.level(target, p.true)
  
  seeds <- 1:nsimu
  
  t <- system.time({
    results <- mclapply(1:nsimu, run.fn, mc.cores=40)
  })
  print(t)
  file.name <- paste0("./data_revision1/","MTDSimu_",ndose, "dose_phi_", target*100, "_fix_",  idx, ".Rdata")
  save(results, file=file.name)
  
  print(post.process.random(results))
  
  cfo.ress <- lapply(1:nsimu, function(i)results[[i]]$cfo)
  acfo.ress <- lapply(1:nsimu, function(i)results[[i]]$acfo)
  crm.ress <- lapply(1:nsimu, function(i)results[[i]]$crm)
  boin.ress <- lapply(1:nsimu, function(i)results[[i]]$boin)
  sum.all <- list(
    CFO = phase1.post.fn(cfo.ress),
    ACFO = phase1.post.fn(acfo.ress),
    CRM = phase1.post.fn(crm.ress),
    BOIN = phase1.post.fn(boin.ress)
  )
  print(tmtd)
  print(phase.I.pretty.tb(sum.all))
}

