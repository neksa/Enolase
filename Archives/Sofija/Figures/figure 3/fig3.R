#source("file:///Home/bccs/sofija/enolase/fig3.R")
library(ggplot2)
library(RColorBrewer)
library("grid")
library(extrafont)
library(gridExtra)

d <- read.table("/Home/bccs/sofija/enolase/fig5_cumulative_coverage/cumulative_coverage_annotated.csv",sep="\t", header=T)

postscript(
    file = "/Home/bccs/sofija/enolase/fig3.eps",
    horizontal=F,
    onefile=F,
    width=6.83 , height=4,
    family=c("/Home/bccs/sofija/arial.afm",
    "/Home/bccs/sofija/arialbd.afm",
    "/Home/bccs/sofija/ariali.afm",
    "/Home/bccs/sofija/arialbi.afm"),
    pointsize=9
)

d1 <- na.omit(sapply(d[1], as.numeric))
d2 <- na.omit(sapply(d[2], as.numeric))
d3 <- na.omit(sapply(d[3], as.character))
d3 <- c(d3)
DF <- data.frame(x=d1, y=d2, z=d3)


cols <- c("#1B9E77", "#E7298A", "#CCCCCC", "#7570B3", "#E6AB02")


a <- c(cols)[DF$z]
u_c <- data.frame(a, d3, stringsAsFactors=FALSE)
u_c <- unique(u_c)

d$profile = factor(d$profile,levels=d$profile,ordered=TRUE)
d$pr_num = c(1:34)

m <- ggplot(d, aes(pr_num, coverage)) 
m <- m + geom_line(aes(pr_num, coverage))
m <- m + geom_point(aes(colour=factor(EF),fill = factor(EF)), shape=21, size = 3) 
m <- m + scale_fill_manual(name="Profiles' EFs",values=cols) 
m <- m + scale_colour_manual(name="Profiles' EFs",values=rep(c("gray30"), each=5))

m <- m + scale_x_continuous('number of profiles',limits = c(1,34),breaks = round(seq(1, 34, by = 5),1))
m <- m + scale_y_continuous('cumulative coverage',limits = c(0.4,1.0),breaks = round(seq(0.4, 1.0, by =0.1),1))
label.these <- c("83","81","82","64","46","84")
m <- m + geom_text(aes(label = profile), color = "black",size = 3,vjust = 1.1,hjust = -0.5, data = d[d$profile %in% label.these, ])
m <- m + theme(
    legend.key = element_rect(fill='white'),
    legend.background = element_rect(fill='white'),
    legend.position = "bottom",
    legend.text = element_text(
        family = "Arial",
        colour="black",
        size=9),
    legend.title = element_text(
        family = "Arial",
        face = "plain",
        colour="black",
        size=9),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
    panel.background = element_rect(fill='white', colour='black'),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.ticks=element_line(colour = 'black'),
    axis.title.x = element_text(size=9,vjust = 0.12, family = "Arial"),
    axis.title.y = element_text(size=9,vjust = 0.12, family = "Arial"),
    axis.text.x = element_text(
        family = "Arial",
        colour="black",
        size=9,
        vjust=1),

    axis.text.y = element_text(
        family = "Arial",
        colour="black",
        size=9,
        hjust=1)
)
m <- m + guides(fill = guide_legend(ncol = 3))



#__________________________________________________________
library(ggplot2)
library(RColorBrewer)
library("grid")
library(extrafont)

d <- read.csv("/Home/bccs/sofija/enolase/fig4_matches_entropy/profile_stats_modified.csv", sep="\t", header=T)

d1 <- na.omit(sapply(d$profile, as.numeric))
d2 <- na.omit(sapply(d$matches, as.numeric))
d3 <- na.omit(sapply(d$entropy, as.numeric))
d4 <- na.omit(sapply(d$EF, as.character))
d4 <- c(d4)
DF <- data.frame(w=d1 ,x=d2, y=d3, z=d4)

cols <- c("#1B9E77", "#E7298A", "#CCCCCC", "#7570B3", "#E6AB02")

a <- c(cols)[DF$z]
u_c <- data.frame(a, d4, stringsAsFactors=FALSE)
u_c <- unique(u_c)

m1 <- ggplot(d, aes(matches, entropy)) 
m1 <- m1 + geom_point(aes(colour=factor(EF),fill = factor(EF)), shape=21, size = 3) 
m1 <- m1 + scale_fill_manual(values=cols) 
m1 <- m1 + scale_colour_manual(values=rep(c("gray30"), each=5))
m1 <- m1 + guides(fill=FALSE)


m1 <- m1 + scale_x_log10('matches, log10 scale', limits = c(1,10000), breaks = c(1,10,100,1000,10000))
m1 <- m1 + scale_y_continuous('entropy, bits',limits = c(0,120),breaks = round(seq(0, 120, by = 20),1))
label.these <- c("83","81","82","64","84")
m1 <- m1 + geom_text(aes(label = profile), color = "black",size = 3,vjust = 0.5,hjust = -0.5, data = d[d$profile %in% label.these, ])
m1 <- m1 + theme(
    legend.position="none",
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
    panel.background = element_rect(fill='white', colour='black'),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.ticks=element_line(colour = 'black'),
    axis.title.x = element_text(size=9,vjust = 0.12, family = "Arial"),
    axis.title.y = element_text(size=9,vjust = 0.12, family = "Arial"),
    axis.text.x = element_text(
        family = "Arial",
        colour="black",
        size=9,
        vjust=1),

    axis.text.y = element_text(
        family = "Arial",
        colour="black",
        size=9,
        hjust=1)
)


g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}
legend <- g_legend(m)

lwidth <- sum(legend$width)

p3 <- grid.arrange(arrangeGrob(m1 + theme(legend.position="none"),
                                        m + theme(legend.position="none"),
                                        widths=c(2/5, 3/5), nrow=1), legend, 
                     heights=c(4/5,1/5), nrow=2)


print(p3)
dev.off()

