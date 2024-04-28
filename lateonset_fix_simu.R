library(magrittr)
library(parallel)
library(survival)
library(BOIN)
library(dfcrm)

source("./setting/Sce_7dose_phi20.R")
source("utilities.R")
source("CFO_tox_utils.R")
source("ACFO_tox_utils.R")
source("lateonset_utils.R")


# 2813
run.fn <- function(j){
  set.seed(seeds[j])
  if (j %% 1 == 0){
    print(j)
  }
  fcfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                      accrual, tite.dist, accrual.dist, design=1, impute.method="frac", add.args=add.args)
  titecfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                         accrual, tite.dist, accrual.dist, design=1, impute.method="TITE", add.args=add.args)
  facfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                       accrual, tite.dist, accrual.dist, design=2, impute.method="frac", add.args=add.args)
  titeacfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                          accrual, tite.dist, accrual.dist, design=2, impute.method="TITE", add.args=add.args)
  titecrm.res <- Simu.Fn (target, p.true, tau, cohortsize, ncohort, 
                      accrual, tite.dist, accrual.dist, design=4, impute.method="TITE", add.args=list())
  titeboin.res <- Simu.Fn (target, p.true, tau, cohortsize, ncohort, 
                       accrual, tite.dist, accrual.dist, design=5, impute.method="TITE", add.args=list())
  benchacfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                           accrual, tite.dist, accrual.dist, design=3, impute.method="No", add.args=add.args)
  ress <- list(
    fcfo=fcfo.res,
    titecfo=titecfo.res,
    facfo=facfo.res,
    titeacfo=titeacfo.res,
    titecrm=titecrm.res,
    titeboin=titeboin.res,
    benchacfo=benchacfo.res,
    paras=list(p.true=p.true, 
               mtd=tmtd, 
               add.args=add.args,
               target=target,
               ncohort=ncohort,
               cohortsize=cohortsize)
  )
  ress
}

safe.run.fn <- function(i) {
  tryCatch(
    run.fn(i),
    error = function(e) {
      cat("Error in run.fn with arguments", list(...), ": ", e$message, "\n")
      i
    }
  )
}

for (i in 2:2)
  {
  tau <- 3
  accrual <- 6
  tite.dist <- 2
  accrual.dist <- 1
  init.dose <- 1
  
  add.args <- list()
  
  nsimu <- 5000
  
  idx <- i
  p.true <- p.trues[[idx]]
  print(p.true)
  ndose <- length(p.true)
  tmtd <- MTD.level(target, p.true)
  
  seeds <- c(1:5000)
  
  t <- system.time({
    results <- mclapply(1:nsimu, run.fn, mc.cores=40)
  })
  print(t)
  file.name <- paste0("./data_revision1/","TITE_MTDSimu_",ndose, "dose_phi_", target*100, "_fix_",  idx, ".Rdata")
  save(results, file=file.name)
  
  print(post.process.random(results))
}
