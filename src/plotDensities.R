args <- commandArgs(trailingOnly = TRUE)
fileName <- args[1]
numStates <- args[2]
out_data <- read.table(fileName, T, sep="\t")

plot.multidens <- function(s, mycolors)
{
        junk.x = NULL
        junk.y = NULL
        for( i in 1:length(s)){
                junk.x = c(junk.x, density(s[[i]])$x)
                junk.y = c(junk.y, density(s[[i]])$y)
        }
        xr <- range(junk.x)
        yr <- range(junk.y)
        xrbreaks <- seq(from=xr[1], to=xr[2], length.out=20)
	d1 <- hist(s[[1]], breaks=xrbreaks, plot=F)
        #d1$y <- length(s[[1]])/sum(d1$y) * d1$y 
        plot(d1, 
             xlim = xr,
             ylim = c(0,100000),
             main = NA, col = mycolors[1], xlab="Average JSD score", ylab="Density",lwd=1)
        for(i in 2:length(s)){
                d1 <- hist(s[[i]], breaks=xrbreaks, plot=F)
		#d1$y <- length(s[[i]])/sum(d1$y) * d1$y
                lines(d1, col = mycolors[i])
        }
}

svg("segments_density.svg", h=5,w=10,pointsize=12)
par(mfrow=c(1,2))
d1 = density(out_data$stat)
d1$y = length(out_data$stat)/sum(d1$y) * d1$y
plot(d1, main="Distribution of segment stats", xlab="Stats", ylab="Segment density") #, xlim=c(0,.1))
if(numStates == 5){
        my_colors <- rainbow(numStates, alpha=0.5)
        plot.multidens(list(as.numeric(out_data$stat[which(out_data$state==5)]),
                        as.numeric(out_data$stat[which(out_data$state==4)]),
                        as.numeric(out_data$stat[which(out_data$state==3)]),
                        as.numeric(out_data$stat[which(out_data$state==2)]),
                        as.numeric(out_data$stat[which(out_data$state==1)])), my_colors)
        legend("topright", c("S5", "S4", "S3", "S2", "S1"), col=my_colors, lwd=2, lty=1, bty="n",cex=1.2)
}
if(numStates == 4){
        my_colors <- rainbow(numStates,alpha=0.5)
        plot.multidens(list(as.numeric(out_data$stat[which(out_data$state==4)]),
                        as.numeric(out_data$stat[which(out_data$state==3)]),
                        as.numeric(out_data$stat[which(out_data$state==2)]),
                        as.numeric(out_data$stat[which(out_data$state==1)])), my_colors)

        legend("topright", c("S4", "S3", "S2", "S1"), col=my_colors, lwd=2, lty=1, bty="n",cex=1.2)
}
if(numStates == 3){
        my_colors <- rainbow(numStates,alpha=0.5)
        plot.multidens(list(as.numeric(out_data$stat[which(out_data$state==3)]),
                        as.numeric(out_data$stat[which(out_data$state==2)]),
                        as.numeric(out_data$stat[which(out_data$state==1)])), my_colors)
        legend("topright", c("S3", "S2", "S1"), col=my_colors, lwd=2, lty=1, bty="n",cex=1.2)
}
dev.off()

