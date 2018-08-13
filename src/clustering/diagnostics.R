#!/usr/local/bin/Rscript --vanilla

DATA <- "/Users/agoncear/projects/Enolase/data/"

library(data.table)


plot_freqs <- function(w) {
    f <- fread(paste(DATA, "F_", w, ".tab", sep=""))
    png(paste("F_", w, ".png", sep=""))
    plot(f$V1, f$V2, type="s", xlab="profile", ylab="freq")
    dev.off()

    g <- fread(paste(DATA, "G_", w, ".tab", sep=""))
    # plot(f$V1, f$V2, type="s", xlab="profile", ylab="freq")
    png(paste("G_", w, ".png", sep=""))
    plot(g$V3, type="s", xlab="profile pair", ylab="freq")
    dev.off()

    png(paste("G_hist_", w, ".png", sep=""))
    hist(g$V3, xlab="G", main="pair freq histogram")
    dev.off()    
}

weights = c('one', 'zero_one', 'nprot', 'log_nprot')
sapply(weights, plot_freqs)

################
m <- fread(paste(DATA, "1/MI_matrix.csv", sep=""), sep=":")
png("MI_distribution.png")
hist(m$V3, xlab="pairwise MI between profile combinations", main="")
dev.off()
