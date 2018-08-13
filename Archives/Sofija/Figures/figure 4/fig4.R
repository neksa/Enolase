#source("file:///Home/bccs/sofija/enolase/clustering/Rw/fig4.R")
ps.options(pointsize = 9, height=1.0, width=6.83)
postscript("/Home/bccs/sofija/enolase/clustering/figures/procedure_figure/part1.eps",
family=c("/Home/bccs/sofija/arial.afm",
"/Home/bccs/sofija/arialbd.afm",
"/Home/bccs/sofija/ariali.afm",
"/Home/bccs/sofija/arialbi.afm"))
par(mfrow=c(1,4))
par(mai = c(0.3, 0.3, 0.3, 0.3)) 
x <- c(247,1,126,186)
for (i in x) { 
    tmp = paste("/Home/bccs/sofija/enolase/clustering/", "_MI_matrix.csv", sep = as.character(i))
    print(tmp)    
    d <- read.table(tmp,sep=":" )
    d3 <- d[3]
    d3 <- sapply(d3, as.numeric)
    d3 <- na.omit(d3)  
    hist(d3,ann=FALSE)
    abline(v=mean(d3),col="blue")
}
dev.off()

ps.options(pointsize = 9, height=1.0, width=6.83)
postscript("/Home/bccs/sofija/enolase/clustering/figures/procedure_figure/part2.eps",
family=c("/Home/bccs/sofija/arial.afm",
"/Home/bccs/sofija/arialbd.afm",
"/Home/bccs/sofija/ariali.afm",
"/Home/bccs/sofija/arialbi.afm"))
par(mfrow=c(1,4))
par(mai = c(0.3, 0.3, 0.3, 0.3)) 
#x <- c(1,126,186,218)
x <- c(218,233,241,245)
for (i in x) { 
    tmp = paste("/Home/bccs/sofija/enolase/clustering/", "_MI_matrix.csv", sep = as.character(i))
    print(tmp)    
    d <- read.table(tmp,sep=":" )
    d3 <- d[3]
    d3 <- sapply(d3, as.numeric)
    d3 <- na.omit(d3)
    hist(d3,ann=FALSE)
    abline(v=mean(d3),col="blue")
}
dev.off()
