#source("file:///Home/bccs/sofija/enolase/clustering/Rw/fig5.R")
#Active clusters & max MI
library(RColorBrewer)
cols <- brewer.pal(4,"Set1")

ps.options(
    width=6.83,
    height=3,
    horizontal = FALSE,
    onefile = FALSE,
    paper = "special",
    pointsize=9)
postscript("/Home/bccs/sofija/enolase/clustering/figures/fig5.eps",
    family=c("/Home/bccs/sofija/arial.afm",
        "/Home/bccs/sofija/arialbd.afm",
        "/Home/bccs/sofija/ariali.afm",
        "/Home/bccs/sofija/arialbi.afm"))
par(mfrow=c(1,2))
#par(mai = c(0.3, 0.3, 0.3, 0.3)) # removes axis labels

d <- read.table("/Home/bccs/sofija/enolase/clustering/max_mi.csv",sep="\t")
d3 <- d[2]
d3 <- sapply(d3, as.numeric)
d3 <- na.omit(d3) 
plot(d3, type="l",  col=cols[2], ylim=range(d3), xlab="number of clusters",ylab="max MI")

d <- read.table("/Home/bccs/sofija/enolase/clustering/active_clusters.csv",header=T,sep="\t")
d1 <- na.omit(sapply(d[2], as.numeric))
plot(d1, type="l", col=cols[1],xlab="number of clusters",ylab="clusters-classificators")

dev.off()

