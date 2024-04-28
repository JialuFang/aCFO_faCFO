library(BOIN)

boin.simu.fn <- function(target, p.true, ncohort, cohortsize, init.level=1){
    if (cohortsize > 1) {
        temp = get.boundary(target, ncohort, cohortsize, n.earlystop = ncohort*cohortsize)$full_boundary_tab
    } else {
        temp = get.boundary(target, ncohort, cohortsize, n.earlystop = ncohort*cohortsize)$boundary_tab
    }
    
    b.e = temp[2, ]
    b.d = temp[3, ]
    b.elim = temp[4, ]
    ndose <- length(p.true)
    y <- rep(0, ndose)
    n <- rep(0, ndose)
    earlystop = 0
    d = init.level
    npts <- cohortsize * ncohort
    elimi = rep(0, ndose)
    for (i in 1:ncohort) {
        newcohort = runif(cohortsize) < p.true[d]
        if ((sum(n) + cohortsize) >= npts) {
            nremain = npts - sum(n)
            y[d] = y[d] + sum(newcohort[1:nremain])
            n[d] = n[d] + nremain
            break
        }
        else {
            y[d] = y[d] + sum(newcohort)
            n[d] = n[d] + cohortsize
        }
        if (!is.na(b.elim[n[d]])) {
            if (y[d] >= b.elim[n[d]]) {
                elimi[d:ndose] = 1
                if (d == 1) {
                    earlystop = 1
                    break
                }
            }
        }
        if (y[d] <= b.e[n[d]] && d != ndose) {
            if (elimi[d + 1] == 0)  # elimination rule
                d = d + 1
        }
        else if (y[d] >= b.d[n[d]] && d != 1) {
            d = d - 1
        }
        else {
            d = d
        }
    }
    if (earlystop == 1) {
        MTD <- 99
    }else {
        MTD <- select.mtd(target, n, y)$MTD
    }
    list(MTD=MTD, dose.ns=n, DLT.ns=y, p.true=p.true, target=target)
}
