#source("file:///Home/bccs/sofija/enolase/clustering/Rw/fig6.R")
#FDR2

ps.options(
    width=3.27,
    height=3,
    horizontal = FALSE,
    onefile = FALSE,
    paper = "special",
    pointsize=9)
postscript("/Home/bccs/sofija/enolase/clustering/figures/fig6.eps",
    family=c("/Home/bccs/sofija/arial.afm",
        "/Home/bccs/sofija/arialbd.afm",
        "/Home/bccs/sofija/ariali.afm",
        "/Home/bccs/sofija/arialbi.afm"))

d <- read.table("/Home/bccs/sofija/enolase/clustering/FDR2(1;260).csv",header=T,sep="\t")

d1 <- na.omit(sapply(d[2], as.numeric))

library(RColorBrewer)
cols <- brewer.pal(4,"Set1")

plot(d1, type="l", col=cols[3],xlab="number of clusters",ylab="FDR")
dev.off()

#png(filename = "/Home/bccs/sofija/enolase/clustering/figures/FDR2.png",height=220, width=235.53 , bg="white")
plot(d1, type="l", col=cols[3],xlab="number of clusters",ylab="FDR")

dev.off()
