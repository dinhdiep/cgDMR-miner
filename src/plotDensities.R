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
        plot(density(s[[1]]), xlim = xr, ylim = yr, main = NA, col = mycolors[1], xlab="Average Dispersion Score", ylab="Density",lwd=1)
        for(i in 2:length(s)){
                lines(density(s[[i]]), xlim = xr, ylim = yr, col = mycolors[i], lwd=1)
        }
}
svg("segments_density.svg", h=5,w=10,pointsize=12)
par(mfrow=c(1,2))
if(numStates == 5){
        my_colors <- rainbow(numStates)
        plot.multidens(list(as.numeric(out_data$stat[which(out_data$stat_states==5)]),
                        as.numeric(out_data$stat[which(out_data$stat_states==4)]),
                        as.numeric(out_data$stat[which(out_data$stat_states==3)]),
                        as.numeric(out_data$stat[which(out_data$stat_states==2)]),
                        as.numeric(out_data$stat[which(out_data$stat_states==1)])), my_colors)
        legend("topright", c("S5", "S4", "S3", "S2", "S1"), col=my_colors, lwd=2, lty=1, bty="n",cex=1.2)
        plot.multidens(list(as.numeric(out_data$stat[which(out_data$stat_random_states==5)]),
                        as.numeric(out_data$stat[which(out_data$stat_random_states==4)]),
                        as.numeric(out_data$stat[which(out_data$stat_random_states==3)]),
                        as.numeric(out_data$stat[which(out_data$stat_random_states==2)]),
                        as.numeric(out_data$stat[which(out_data$stat_random_states==1)])), my_colors)
        legend("topright", c("rS5", "rS4", "rS3", "rS2", "rS1"), col=my_colors, lwd=2, lty=1, bty="n",cex=1.2)
}
if(numStates == 4){
        my_colors <- rainbow(numStates)
        plot.multidens(list(as.numeric(out_data$stat[which(out_data$stat_states==4)]),
                        as.numeric(out_data$stat[which(out_data$stat_states==3)]),
                        as.numeric(out_data$stat[which(out_data$stat_states==2)]),
                        as.numeric(out_data$stat[which(out_data$stat_states==1)])), my_colors)

        legend("topright", c("S4", "S3", "S2", "S1"), col=my_colors, lwd=2, lty=1, bty="n",cex=1.2)
        plot.multidens(list(as.numeric(out_data$stat[which(out_data$stat_random_states==4)]),
                        as.numeric(out_data$stat[which(out_data$stat_random_states==3)]),
                        as.numeric(out_data$stat[which(out_data$stat_random_states==2)]),
                        as.numeric(out_data$stat[which(out_data$stat_random_states==1)])), my_colors)
        legend("topright", c("rS4", "rS3", "rS2", "rS1"), col=my_colors, lwd=2, lty=1, bty="n",cex=1.2)
}
if(numStates == 3){
        my_colors <- rainbow(numStates)
        plot.multidens(list(as.numeric(out_data$stat[which(out_data$stat_states==3)]),
                        as.numeric(out_data$stat[which(out_data$stat_states==2)]),
                        as.numeric(out_data$stat[which(out_data$stat_states==1)])), my_colors)
        legend("topright", c("S3", "S2", "S1"), col=my_colors, lwd=2, lty=1, bty="n",cex=1.2)
        plot.multidens(list(as.numeric(out_data$stat[which(out_data$stat_random_states==3)]),
                        as.numeric(out_data$stat[which(out_data$stat_random_states==2)]),
                        as.numeric(out_data$stat[which(out_data$stat_random_states==1)])), my_colors)
        legend("topright", c("rS3", "rS2", "rS1"), col=my_colors, lwd=2, lty=1, bty="n",cex=1.2)
}
dev.off()

