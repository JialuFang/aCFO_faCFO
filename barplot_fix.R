library(ggplot2)
library(cowplot)
library(gridExtra)
source("utilities.R")
source("./setting/Sce_5dose_phi33.R")
ndose <- length(p.trues[[1]])
res.ls <- list()

for (flag in 1:length(p.trues)){
  print(flag)
  file.name <- paste0("./data_revision1/","TITE_MTDSimu_",ndose, "dose_phi_", target*100, "_fix_",  flag, ".Rdata")
  load(file.name)
  res.ls[[flag]] <-  post.process.random(results)[c("fcfo", "facfo", "titecfo", "titeacfo", "titecrm","titeboin", "benchacfo"), ]
  #res.ls[[flag]] <-  post.process.random(results)[c("cfo","acfo","crm","boin"), ]
  print(res.ls[[flag]])
  flag <- flag + 1
}

res.ls <- lapply(res.ls, function(df) {
  df[is.na(df)] <- 0
  return(df)
})

res.ls[[length(res.ls)+1]] <- Reduce(`+`, res.ls)/ length(res.ls)

res.ls <- lapply(res.ls, function(df) {
  df[, !colnames(df) %in% "Risk.of.HT", drop = FALSE]
})

dim <- dim(res.ls[[1]])


#title<- c("MTD selection", "MTD allocation", "Overdose selection", "Overdose allocation", "Average DLT rate")
title<- c("MTD selection", "MTD allocation", "Overdose selection", "Overdose allocation", 
          "Average DLT rate", "Average trial duration")


plots_list <- list()
for (i in 1:(length(title)-1)){
  grp.names <- c(1:length(p.trues), 'avg')
  #m.names <- c("CFO", "aCFO", "CRM", "BOIN")
  m.names <-   c("fCFO", "f-aCFO", "TITE-CFO", "TITE-aCFO", "TITE-CRM", "TITE-BOIN", "Benchmark aCFO")
  g.var <- rep(grp.names, each=length(m.names))
  m.var <- rep(m.names, times=length(p.trues)+1)
  v.var <- unlist(lapply(res.ls, function(x) x[,i]))*100
  data <- data.frame(g=factor(g.var, levels=grp.names), m=factor(m.var, levels=m.names), v=v.var)
  
  p <- ggplot(data = data, mapping = aes(x = g, y = v, fill = m)) + geom_bar(stat = 'identity', position = 'dodge') +
    theme(legend.key.size = unit(2, 'mm'), 
          legend.text = element_text(size = 6),
          legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
          legend.position = 'right',
          plot.title = element_text(size = 12,hjust=0.5)) +
    xlab("") + ylab("Percentage (%)") + 
    guides(fill=guide_legend(title = "")) + 
    ggtitle(title[i])
  
  plots_list[[i]] <- p
}

i <- length(title)
grp.names <- c(1:length(p.trues), 'avg')
#m.names <- c("CFO", "aCFO", "CRM", "BOIN")
m.names <-   c("fCFO", "f-aCFO", "TITE-CFO", "TITE-aCFO", "TITE-CRM", "TITE-BOIN", "Benchmark aCFO")
g.var <- rep(grp.names, each=length(m.names))
m.var <- rep(m.names, times=length(p.trues)+1)
v.var <- unlist(lapply(res.ls, function(x) x[,i]))
data <- data.frame(g=factor(g.var, levels=grp.names), m=factor(m.var, levels=m.names), v=v.var)

p <- ggplot(data = data, mapping = aes(x = g, y = v, fill = m)) + geom_bar(stat = 'identity', position = 'dodge') +
  theme(legend.key.size = unit(2, 'mm'), 
        legend.text = element_text(size = 6),
        legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
        legend.position = 'right',
        plot.title = element_text(size = 12,hjust=0.5)) +
  xlab("") + ylab("Time (in months)") + 
  guides(fill=guide_legend(title = "")) + 
  ggtitle(title[i])

plots_list[[i]] <- p

grid.arrange(grobs = plots_list, nrow = 3, ncol = 2)
