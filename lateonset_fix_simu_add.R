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


run.fn <- function(j){
  set.seed(seeds[j])
  if (j %% 100 == 0){
    print(j)
  }
  titecrm.res <- Simu.Fn (target, p.true, tau, cohortsize, ncohort, 
                          accrual, tite.dist, accrual.dist, design=4, impute.method="TITE", add.args=list())
  titeboin.res <- Simu.Fn (target, p.true, tau, cohortsize, ncohort, 
                           accrual, tite.dist, accrual.dist, design=5, impute.method="TITE", add.args=list())
  ress <- list(
    titecrm=titecrm.res,
    titeboin=titeboin.res,
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
  file.name <- paste0("./data_revision1/","add_TITE_MTDSimu_",ndose, "dose_phi_", target*100, "_fix_",  idx, ".Rdata")
  save(results, file=file.name)
  
  print(post.process.random(results))
}
