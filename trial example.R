library(CFO)

cfo <- CFO.simu(design = 'CFO', target, p.true, init.level, ncohort, cohortsize, seed = 13)
acfo <- CFO.simu(design = 'aCFO', target, p.true, init.level, ncohort, cohortsize, seed = 13)
cfo$MTD
acfo$MTD

x1_no <- c(0.75,1,1.25, 1.75,2.25, 2.75,3,3.25, 3.75,4,4.25, 4.75,5,5.25, 5.75,6,6.25, 7,7.25, 7.75,8,8.25, 8.75,9,9.25,
           9.75,10,10.25, 10.75,11, 12,12.25)
y1_no <- c(1,1,1, 2,2, 2,2,2,  2, 2, 2,  3, 3, 3,  4, 4, 4,  5, 5,  5, 5, 5,  5, 5, 5,  6, 6, 6,  7, 7,  7, 7) +0.1
x1_yes <-c(2, 6.75, 11.25, 11.75)
y1_yes <-c(2,    5,    7,      7) + 0.1

x2_no <- c(0.75,1,1.25, 1.75,2.25, 2.75,3,3.25, 3.75,4,4.25, 4.75,5,5.25, 5.75,6,6.25, 7,7.25, 7.75,8,8.25, 8.75,9,9.25, 
           9.75,10,10.25, 10.75,11, 12,12.25)
y2_no <- c(1,1,1, 2,2, 2,2,2,  2, 2, 2,  3, 3, 3,  4, 4, 4,  5, 5,  4, 4, 4,  5, 5, 5,  5, 5, 5,  6, 6,  5, 5) - 0.1
x2_yes <-c(2, 6.75, 11.25, 11.75)
y2_yes <-c(2,    5,    6,  5) - 0.1


fig.name <- paste0("./use_revision1/", "Trial Example", ".jpeg")
png(filename=fig.name, unit="in", height=7, width=12, res=300)
par(cex.axis = 1.8, mar = c(5,5,2,2))
plot(x1_no, y1_no, pch = 1, cex=1.8, xlab = "cohort index", ylab = "Dose level", 
     cex.lab = 1.8, ylim = c(0.4,7.5), xaxt="n")
axis(side=1, at=1:14)
points(x1_yes, y1_yes, pch = 16, cex=1.8)
points(x2_no, y2_no, pch = 0, cex=1.8)
points(x2_yes, y2_yes, pch = 15, cex=1.8)

legend("bottomright", 
       pch = c(1, 16, 0, 15), 
       c("no DLT in aCFO ", "DLT in aCFO", "no DLT in CFO ", "DLT in CFO"),
       cex = 1.5)
dev.off()
