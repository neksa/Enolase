
a <- read.csv("profile_stats.csv", sep="\t", header=T)

ps.options(width=3.42, height=3.42, horizontal = FALSE, onefile = FALSE, paper = "special", pointsize=8)

postscript("profile_entropy_vs_K.eps")
plot(a$matches, a$entropy, log="x", pch=" ", ylim=c(0, 31*log2(20)), xlim=c(1, max(a$matches) * 10), xlab="Matches (K), log10 scale", ylab="Entropy, bits")

# plot("", ylim=c(0,120), xlim=c(0, 12200), log="x", xlab="Matches (K), log10 scale", ylab="Entropy, bits")
# text(log10(a$matches), a$entropy, labels=a$profile) 

text(a$matches, a$entropy, labels=a$profile) 
dev.off()
