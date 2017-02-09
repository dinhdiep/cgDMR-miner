msg.trap <- capture.output( suppressMessages( library( DNAcopy ) ) )
msg.trap <- capture.output( suppressMessages( library( fastseg ) ) )
msg.trap <- capture.output( suppressMessages( library( mhsmm ) ) )
msg.trap <- capture.output( suppressMessages( library( data.table ) ) )
msg.trap <- capture.output( suppressMessages( library( GenomicRanges ) ) )
msg.trap <- capture.output( suppressMessages( library( RColorBrewer ) ) )
msg.trap <- capture.output( suppressMessages( library( methods ) ) )
msg.trap <- capture.output( suppressMessages( library( mgcv ) ) )
msg.trap <- capture.output( suppressMessages( library( graphics ) ) )
msg.trap <- capture.output( suppressMessages( library( entropy ) ) )
msg.trap <- capture.output( suppressMessages( library( MASS ) ) )

args <- commandArgs(trailingOnly = TRUE)
fileName <- args[1]
outName <- args[2]
mode <- args[3]

my_cpu <- 2

df <- read.table(fileName, T, sep="\t", stringsAsFactor=FALSE)
summary(df)

get_hmm_seg <- function(df){
        max.gap = 10000
	segments <- data.frame()
	o1 <- order(df$pos)
	meandiff <- df$stat[o1]
	pos <- df$pos[o1]
	data.gap <- diff(pos) > max.gap
	gap.idx <- which(data.gap)
	starts<-c(1, gap.idx + 1)
	ends<-c(gap.idx, length(pos))
	keep.pos <- pos
	if(length(starts) == 1){
        	meandiff.data <- list(x=meandiff, N=length(pos))
	}else{
	        meandiff.lst <- mapply(function(s,e, data) data[s:e], starts, ends, MoreArgs=list(meandiff))
        	####
	        pos.lst <- mapply(function(s,e, data) data[s:e], starts, ends, MoreArgs=list(pos))
        	list.lens <- as.numeric(summary(pos.lst)[,1])
	        keep.rows <- which(list.lens >= 5) #need at least 5 data points per segment and have variabity in methylation means
        	####
        	meandiff.data <- list(x=c(meandiff.lst[keep.rows], recursive=TRUE), N=list.lens[keep.rows])
        	keep.pos <- unlist(pos.lst[keep.rows])
	}


	# setting up the initial probabilities
	N=5
	bins <- hist(meandiff)
	stat_counts <- data.frame(counts=bins$counts, mids=bins$mids)
	mypeaks <- sort(stat_counts$mids[order(-stat_counts$counts)][1:N])
	initial<-rep(1/N, N)
	#A<- matrix(c(0.5,0.5,0, 0.4,0.2,0.4, 0,0.5,0.5), N, byrow = TRUE)
	A<- matrix(c(0.5,0.5,0,0,0, 0.3,0.4,0.3,0,0, 0,0.3,0.4,0.3,0, 0,0,0.3,0.4,0.3, 0,0,0,0.5,0.5), N, byrow = TRUE)
	B<- list(mu = mypeaks, sigma=rep(sd(meandiff),N))
	model0 <- hmmspec(init=initial, trans=A, parms.emis=B, dens.emis=dnorm.hsmm)
	print(model0)
	H <- hmmfit(meandiff.data, model0, mstep=mstep.norm)
	print(H$model)
	meandiff.s <- predict.hmmspec(H$model, meandiff.data, method="viterbi")

	mydf <- data.frame(pos=as.numeric(keep.pos), stat=as.numeric(unlist(meandiff.data$x)),
                        	state = as.numeric(meandiff.s$s), id=rep(1, length(keep.pos)))
	
        cur_id = 1
        mydf$id[1] = cur_id
	cur_state = mydf$state[1]
        for( j in 2:length(keep.pos) ){
        	if(cur_state != mydf$state[j] || mydf$pos[j] - mydf$pos[j-1] > max.gap){
        		cur_id = cur_id+1
                	cur_state = mydf$state[j]
        	}
		mydf$id[j] = cur_id
        }
       	collapsed <- data.table(mydf, key="id")
	segments <- data.frame(start=c(segments$start, collapsed[,min(pos),by=id]$V1),
			end=c(segments$end, collapsed[,max(pos)+1,by=id]$V1), 
			stat=c(segments$stat, collapsed[,mean(stat),by=id]$V1),
			state=c(segments$state, collapsed[,mean(state),by=id]$V1))
        write.table(data.frame(chrom=rep(df$chr[1], nrow(segments)), segments), file=paste0(outName, ".seg.txt"), quote=F, row=F, sep="\t")
        write.table(mydf, file=paste0(outName, ".HMM.sites.txt"), quote=F, row=F, sep="\t")
}

get_fastseg <- function(df){
        max.gap = 10000
	cur_chrom = df$chr[1]
        o1 <- order(df$pos)
        mystat <- df$stat[o1]
        mypos <- df$pos[o1]
        data.gap <- diff(df$pos[o1]) > max.gap
        gap.idx <- which(data.gap)
        starts <- c(1, gap.idx+1)
        ends <- c(gap.idx, length(o1))
        segments <- data.frame()
        cur_id = 1;
        for( i in 1:length(starts) ){
                # only process if segment is >= 10
                e = ends[i]
                s = starts[i]
                if(e - s >= 10){
                        seg.data <- data.frame(stat = mystat[s:e], pos=mypos[s:e], id=c(s:e))
                      	res <- fastseg(seg.data$stat, minSeg=1, cyberWeight=20, alpha=0.25) 
			N=length(res$ID)
                        for( j in 1:N ){
                                for( k in res$startRow[j]:res$endRow[j] ){
					seg.data$id[k] = cur_id
				}
				cur_id = cur_id+1
                        }
                        cur_id = cur_id+1
                        collapsed <- data.table(seg.data, key="id")
                        segments <- data.frame(start=c(segments$start, collapsed[,min(pos),by=id]$V1),
                                        end=c(segments$end, collapsed[,max(pos)+1,by=id]$V1),
                                        stat=c(segments$stat, collapsed[,mean(stat),by=id]$V1))
                }
        }
	segments <- data.frame(chrom=rep(cur_chrom, nrow(segments)), segments)
        write.table(segments, file=paste0(outName, ".seg.txt"), quote=F, row=F, sep="\t")
}

get_segment <- function(df){
        o1 <- order(df$pos)
        CNA.object <- CNA(df$stat[o1], df$chr, df$pos[o1], data.type="logratio", sampleid="test")
        smoothed.CNA.object <- smooth.CNA(CNA.object)
        segment.smoothed.CNA.object <- segment(smoothed.CNA.object,verbose=1)
        segment.smoothed.CNA.object$output$loc.end = segment.smoothed.CNA.object$output$loc.end + 1
        write.table(segment.smoothed.CNA.object$output[,2:6], file=paste0(outName, ".seg.txt"), quote=F, row=F, sep="\t")
}

plotHeatScatter <- function(x,y,main,xlab,ylab,numBins){
	h1 <- hist(x, plot=F)
	h2 <- hist(y, plot=F)
	top <- max(h1$counts, h2$counts)
	k <- kde2d(x, y, n=numBins)
	rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))	
	r <- rf(32)
	ColorLevels <- round(seq(min(k$z)*length(k$x)/sum(k$x), max(k$z)*length(k$x)/sum(k$x), length=5),2)
	# margins
	oldpar <- par()
	par(mar=c(3,3,1,1))
	layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
	image(k, col=r) #plot the image
	legend("topright", legend=ColorLevels, pch=19, col=rf(5),bg="white")
	par(mar=c(0,2,1,0))
	barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
	par(mar=c(2,0,0.5,1))
	barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)
}


my_meandiff_function <- function(x){
        controls_mean <- mean(as.numeric(x[controls_pos]), na.rm=T)
        cases_mean <- mean(as.numeric(x[cases_pos]), na.rm=T)
        if(is.na(controls_mean) || is.na(cases_mean)){
                return (NA)
        }else{
                return (cases_mean - controls_mean)
        }
}

my_log_coefvsqr_function <- function(x){
        a <- as.numeric(x[keep_col])
        cur_mean <- mean(a, na.rm=TRUE)
        cur_sd <- sd(a,na.rm=TRUE)
        if(is.na(cur_sd) || is.na(cur_mean)){
                return(NA)
        }
        if(cur_mean > 0 & cur_sd > 0){
                y <- cur_sd/cur_mean
                return (log(y*y))
        }else{
                return(NA)
        }
}

my_entropy_function <- function(x){
	dat <- as.numeric(x[keep_col])
	if(sd(dat, na.rm=T) == 0 || is.na(sd(dat,na.rm=T))){
		return(0)
	}
	if(is.na(mean(dat, na.rm=T))){
		return(NA)
	}else{
		dat[is.na(dat)] = median(dat, na.rm=T)
		return(entropy(dat))
	}
}

my_mean_function <- function(x){
	a <- as.numeric(x[keep_col])
	return(mean(a, na.rm=TRUE))
}

my_jsD_function <- function(x){
	dat <- as.numeric(x[keep_col])
	# interpolate the missing values
	if(sd(dat, na.rm=T) == 0 || is.na(sd(dat,na.rm=T))){
		return(0)
	}
	if(is.na(mean(dat, na.rm=T))){
		return(NA)
	}else{
		dat[is.na(dat)] = median(dat, na.rm=T)
		n <- length(dat)
		myQ <- rep(1/n,n)
		scaled_dat <- scale(dat, center=F)
		myP <- scaled_dat/sum(scaled_dat)
		myMid <- 0.5*(myP+myQ)
		myJSD <- sqrt(0.5*KL.plugin(myP, myMid)+0.5*KL.plugin(myQ, myMid))
		return(myJSD)
	}
}

#get_fastseg(df)

if(grepl("CBS", mode)){
	get_segment(df)
}else{
	get_hmm_seg(df)
}

