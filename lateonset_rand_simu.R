library(magrittr)
library(parallel)
library(survival)
library(BOIN)
library(dfcrm)

source("utilities.R")
source("CFO_tox_utils.R")
source("ACFO_tox_utils.R")
source("lateonset_utils.R")

tau <- 3
accrual <- 6
tite.dist <- 2
accrual.dist <- 1
init.dose <- 1

target <- 0.2
ndose <- 7
ncohort <- 12
cohortsize <- 3

diff_list <- c(0.05, 0.07, 0.1, 0.15)
mu_list <- c(0.36, 0.50, 0.65, 0.87)
for (i in 1:4){
  mu <- mu_list[i]
  diff <- diff_list[i]
  run.fn <- function(k){
    set.seed(seeds[k])
    if (k %% 1 == 0){
      print(k)
    }
    p.true.all <- gen.rand.doses(ndose, target, mu1=mu, mu2=mu)
    p.true <- p.true.all$p.true
    tmtd <- p.true.all$mtd.level
    
    fcfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                        accrual, tite.dist, accrual.dist, design=1, impute.method="frac", add.args=list())
    titecfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                           accrual, tite.dist, accrual.dist, design=1, impute.method="TITE", add.args=list())
    facfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                         accrual, tite.dist, accrual.dist, design=2, impute.method="frac", add.args=list())
    titeacfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                            accrual, tite.dist, accrual.dist, design=2, impute.method="TITE", add.args=list())
    titecrm.res <- Simu.Fn (target, p.true, tau, cohortsize, ncohort, 
                            accrual, tite.dist, accrual.dist, design=4, impute.method="TITE", add.args=list())
    titeboin.res <- Simu.Fn (target, p.true, tau, cohortsize, ncohort, 
                             accrual, tite.dist, accrual.dist, design=5, impute.method="TITE", add.args=list())
    benchacfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                             accrual, tite.dist, accrual.dist, design=3, impute.method="No", add.args=list())
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
  
  nsimu <- 5000
  seeds <- 1:nsimu
  
  file.name <- paste0("./data_revision1/","TITE_MTDSimu_",ndose, "dose_phi_", target*100, "_random_",  diff*100, ".Rdata")
  save(results, file=file.name)
  
  t <- system.time({
    results <- mclapply(1:nsimu, run.fn, mc.cores=1)
  })
  print(t)
  
  post.process.random(results)

}

