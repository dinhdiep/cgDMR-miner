msg.trap <- capture.output( suppressMessages( library( mhsmm ) ) )
msg.trap <- capture.output( suppressMessages( library( methods ) ) )

args <- commandArgs(trailingOnly = TRUE)
fileName <- args[1]
outName <- args[2]
statCol <- args[3]
max.gap <- 10000 # 10 kbp

# Read in data table
rawdata <- read.table(fileName, T)
data <- na.omit(rawdata)
head(data)

# positions 
pos <- as.integer(data[,"middle_pos"]) 

# stats
stats_0 <- data[, statCol]
stats_0[which(stats_0 < 0.0001)] <- 0.0001 # make the minimum 0.0001

if(statCol=="jsd"){
  print("values are jsd, performs logit transformation")
  meandiff <- log(stats_0/(1 - stats_0),10)
  #meandiff <- stats_0
}else{
  print("values are chsq, performs log transformation")
  meandiff <- log(stats_0, 10)
}

print(range(meandiff))

# use the max.gap to create segments first
data.gap <- diff(pos) > max.gap
gap.idx <- which(data.gap)
starts<-c(1, gap.idx + 1)
ends<-c(gap.idx, length(pos))

if(length(starts) == 1){
  meandiff.data <- list(x=meandiff, N=length(pos))
  pos.lst <- list(x=pos)
  stats.lst <- list(x=stats_0)
}else{
  meandiff.lst <- mapply(function(s,e, data) data[s:e], starts, ends, MoreArgs=list(meandiff))
  pos.lst <- mapply(function(s,e, data) data[s:e], starts, ends, MoreArgs=list(pos))
  stats.lst <-mapply(function(s,e, data) data[s:e], starts, ends, MoreArgs=list(stats_0))

  # Find segments with five or fewer sites to remove
  list.lens <- as.numeric(summary(pos.lst)[,1])
  keep.rows <- which(list.lens >= 5) 
 
  meandiff.data <- list(x=c(meandiff.lst[keep.rows], recursive=TRUE), N=list.lens[keep.rows])
  pos.data <- list(x=c(pos.lst[keep.rows], recursive=TRUE))
  stats.data <- list(x=c(stats.lst[keep.rows], recursive=TRUE))
}

# tracking which position are kept
keep.pos <- unlist(pos.data)

# setting up the initial probabilities for HMM

# 4 HMM states
N <- 4
initial<-rep(1/N, N)

# find the peaks which are the starting mu
bins <- hist(meandiff, plot=FALSE)
stat_counts <- data.frame(counts=bins$counts, mids=bins$mids)
peaks <- sort(stat_counts$mids[order(-stat_counts$counts)][1:N])

# transition and emission probabilities
A<- matrix(c(0.6,0.3,0.05,0.05, 0.175,0.6,0.175,0.05, 0.05,0.175,0.6,0.175, 0.05,0.05,0.3,0.6), N, byrow = TRUE)
B<- list(mu = peaks, sigma=rep(sd(meandiff),N))
model0 <- hmmspec(init=initial, trans=A, parms.emis=B, dens.emis=dnorm.hsmm)
print(model0)

H <- hmmfit(meandiff.data, model0, mstep=mstep.norm)
print(H$model)

meandiff.s <- predict.hmmspec(H$model, meandiff.data, method="viterbi")	
	
mydf <- data.frame(pos=as.numeric(keep.pos), stat=as.numeric(unlist(stats.data)),
		state = as.numeric(meandiff.s$s))

write.table(mydf, file=paste0(outName, ".HMM4.sites.txt"), quote=F, row=F, sep="\t")

