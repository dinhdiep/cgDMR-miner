msg.trap <- capture.output( suppressMessages( library( DNAcopy ) ) )
msg.trap <- capture.output( suppressMessages( library( fastseg ) ) )
msg.trap <- capture.output( suppressMessages( library( depmixS4 ) ) )
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
analysis_type <- as.numeric(args[2])
label_groups <- args[3]
outName <- args[4]

my_cpu <- 2

data <- read.table(fileName, T, sep="\t", stringsAsFactor=FALSE)
summary(data)
my_train_chr <- data$chr[1]
using_samples = read.table(label_groups, F, sep="\t", stringsAsFactor=FALSE)
keep_col = unlist(sapply(using_samples$V1, grep, colnames(data)))

# for meandiff function
controls <- using_samples[which(grepl("ctrl", using_samples$V2)),]
cases <- using_samples[which(grepl("case", using_samples$V2)),]
controls_pos = unlist(sapply(controls$V1, grep, colnames(data)))
cases_pos = unlist(sapply(cases$V1, grep, colnames(data)))
print(controls_pos)
print(cases_pos)
####################

get_hmm_seg <- function(df){
        max.gap = 10000
        o1 <- order(df$pos)
        mystat <- df$stat[o1]
        mypos <- df$pos[o1]
        data.gap <- diff(df$pos[o1]) > max.gap
        gap.idx <- which(data.gap)
        starts <- c(1, gap.idx+1)
        ends <- c(gap.idx, length(o1))
        mydf <- data.frame()
	segments <- data.frame()
        cur_id = 1;
	#initTrans = matrix(c(0.5,0.5,0,0,0, 0.3,0.4,0.3,0,0, 0,0.3,0.4,0.3,0, 0,0,0.3,0.4,0.3, 0,0,0,0.5,0.5), 5, byrow = TRUE)
	initTrans = matrix(rep(0.2,5*5), 5, byrow = TRUE)
        for( i in 1:length(starts) ){
                # only process if segment is >= 10
                e = ends[i]
                s = starts[i]
                if(e - s >= 10){
                        hmm.data <- data.frame(stat = mystat[s:e], pos=mypos[s:e], id=c(s:e))
                        mod <- depmix(response = stat ~ 1, data=hmm.data, nstates=5, trstart=initTrans, family=gaussian())
                        fm <- fit(mod, emcontrol=em.control(tol=1e-6, crit="relative"))
                        fitted.state <- posterior(fm)$state
                        N=length(fitted.state)
                        cur_state = fitted.state[1]
                        hmm.data$id[1] = cur_id
                        for( j in 2:N ){
                                if(cur_state != fitted.state[j]){
                                        cur_id = cur_id+1
                                        cur_state = fitted.state[j]
                                }
                                hmm.data$id[j] = cur_id
                        }
                        cur_id = cur_id+1
			hmm.data$state = fitted.state
                        mydf <- data.frame(pos=c(mydf$pos, hmm.data$pos), stat=c(mydf$stat, hmm.data$stat),
                                state = c(mydf$state, fitted.state), id=c(mydf$id, hmm.data$id))
        		collapsed <- data.table(hmm.data, key="id")
			segments <- data.frame(start=c(segments$start, collapsed[,min(pos),by=id]$V1),
					end=c(segments$end, collapsed[,max(pos)+1,by=id]$V1), 
					stat=c(segments$stat, collapsed[,mean(stat),by=id]$V1),
					state=c(segments$state, collapsed[,mean(state),by=id]$V1))
                }
        }
        write.table(data.frame(chrom=rep(df$chr[1], nrow(segments)), segments), file=paste0(outName, ".HMM.seg.txt"), quote=F, row=F, sep="\t")
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
        write.table(segments, file=paste0(outName, ".Cyber.seg.txt"), quote=F, row=F, sep="\t")
}

get_segment <- function(df){
        o1 <- order(df$pos)
        CNA.object <- CNA(df$stat[o1], df$chr, df$pos[o1], data.type="logratio", sampleid="test")
        smoothed.CNA.object <- smooth.CNA(CNA.object)
        segment.smoothed.CNA.object <- segment(smoothed.CNA.object,verbose=1)
        segment.smoothed.CNA.object$output$loc.end = segment.smoothed.CNA.object$output$loc.end + 1
        write.table(segment.smoothed.CNA.object$output[,2:6], file=paste0(outName, ".CNA.seg.txt"), quote=F, row=F, sep="\t")
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

if(analysis_type==6){
	# calculating
        print("Begin calculating summary statistics...")
        myMeanDiff<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_meandiff_function(x)), mc.cores=my_cpu))
        myMean<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_mean_function(x)), mc.cores=my_cpu))
	print(summary(myMeanDiff))
	print(summary(myMean))
        keep_rows <- which(!is.na(myMeanDiff))
	# plot
	pdf(paste0(outName, ".meandiff.pdf"))
        plotHeatScatter(x=myMean[keep_rows], y=myMeanDiff[keep_rows], main="Groups mean difference", ylab="Groups mean difference", xlab="Mean methylation", numBins=100)
        dev.off()
        # write table
	df <- data.frame(chr=data$chr, pos=data$start, stat=myMeanDiff)
        write.table(df,file=paste0(outName,".meandiff.txt"), quote=F,sep="\t",row=F,col=T)
	#get_fastseg(df)
	#get_segment(df)
	#get_hmm_segment(df)
}

if(analysis_type == 5){
	# calculating
        print("Begin calculating summary statistics...")
        myEntropy<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_entropy_function(x)), mc.cores=my_cpu))
        myMean<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_mean_function(x)), mc.cores=my_cpu))
        print(summary(myEntropy))
        print(summary(myMean))
        keep_rows <- which(!is.na(myEntropy) & !is.na(myMean))
        # plot
        pdf(paste0(outName, ".entropy.pdf"))
        plotHeatScatter(x=myMean[keep_rows], y=myEntropy[keep_rows], main="Shannon Entropy", ylab="Shannon Entropy", xlab="Mean methylation", numBins=100)
        dev.off()
        # write table
        df <- data.frame(chr=data$chr, pos=data$start, stat=myEntropy)
        write.table(df,file=paste0(outName,".entropy.txt"), quote=F,sep="\t",row=F,col=T)
        #get_fastseg(df)
        #get_segment(df)
	#get_hmm_segment(df)
}

if(analysis_type == 1){
	# calculating 
	print("Begin calculating summary statistics...")
	my_jsD_hypo<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_jsD_function(x)), mc.cores=my_cpu))
	my_jsD_hyper<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) {y=1-as.numeric(x); my_jsD_function(y);}), mc.cores=my_cpu))
	myMean<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_mean_function(x)), mc.cores=my_cpu))
	myMinJSD <- pmin(my_jsD_hypo, my_jsD_hyper, na.rm=T)
	#hypo_rows <- which(my_jsD_hypo <= my_jsD_hyper & !is.na(my_jsD_hypo))
	#print(myMinJSD)
	#if(!is.na(hypo_rows)){
	#	myMinJSD[hypo_rows] = -1 * myMinJSD[hypo_rows]
	#}
	# plot
	pdf(paste0(outName, ".JSD.pdf"))
	keep_rows = which(!is.na(my_jsD_hypo) & !is.na(myMean))
	plotHeatScatter(x=myMean[keep_rows], y=my_jsD_hypo[keep_rows], 
		main="Jensen-Shannon Distance", ylab="Jensen-Shannon Distance", xlab="Mean methylation", numBins=100)
	keep_rows = which(!is.na(my_jsD_hyper) & !is.na(myMean))
	plotHeatScatter(x=myMean[keep_rows], y=my_jsD_hyper[keep_rows], 
		main="Jensen-Shannon Distance", ylab="Jensen-Shannon Distance", xlab="Mean methylation", numBins=100)
	keep_rows = which(!is.na(myMinJSD) & !is.na(myMean))
	plotHeatScatter(x=myMean[keep_rows], y=myMinJSD[keep_rows], 
		main="Jensen-Shannon Distance", ylab="Jensen-Shannon Distance", xlab="Mean methylation", numBins=100)
	dev.off()
	# write table
	df <- data.frame(chr=data$chr[keep_rows], pos=data$start[keep_rows], stat=myMinJSD[keep_rows])
	write.table(df, file=paste0(outName,".jsd.txt"), quote=F,sep="\t",row=F,col=T)
	#get_fastseg(df)
	#get_segment(df)
	#get_hmm_seg(df)
}

if(analysis_type == 2){
        # calculating
        print("Begin calculating summary statistics...")
        my_jsD_hypo<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_jsD_function(x)), mc.cores=my_cpu))
        my_jsD_hyper<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) {y=1-as.numeric(x); my_jsD_function(y);}), mc.cores=my_cpu))
        myMean<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_mean_function(x)), mc.cores=my_cpu))
        myMinJSD <- pmin(my_jsD_hypo, my_jsD_hyper, na.rm=T)
        hypo_rows <- which(my_jsD_hypo <= my_jsD_hyper & !is.na(my_jsD_hypo))
        #print(myMinJSD)
        if(!is.na(hypo_rows)){
               myMinJSD[hypo_rows] = -1 * myMinJSD[hypo_rows]
        }
        # plot
        pdf(paste0(outName, ".JSD.pdf"))
        keep_rows = which(!is.na(my_jsD_hypo) & !is.na(myMean))
        plotHeatScatter(x=myMean[keep_rows], y=my_jsD_hypo[keep_rows],
                main="Jensen-Shannon Distance", ylab="Jensen-Shannon Distance", xlab="Mean methylation", numBins=100)
        keep_rows = which(!is.na(my_jsD_hyper) & !is.na(myMean))
        plotHeatScatter(x=myMean[keep_rows], y=my_jsD_hyper[keep_rows],
                main="Jensen-Shannon Distance", ylab="Jensen-Shannon Distance", xlab="Mean methylation", numBins=100)
        keep_rows = which(!is.na(myMinJSD) & !is.na(myMean))
        plotHeatScatter(x=myMean[keep_rows], y=myMinJSD[keep_rows],
                main="Jensen-Shannon Distance", ylab="Jensen-Shannon Distance", xlab="Mean methylation", numBins=100)
        dev.off()
        # write table
        df <- data.frame(chr=data$chr[keep_rows], pos=data$start[keep_rows], stat=myMinJSD[keep_rows])
        write.table(df, file=paste0(outName,".jsd.txt"), quote=F,sep="\t",row=F,col=T)
        #get_fastseg(df)
        #get_segment(df)
        #get_hmm_seg(df)
}

if(analysis_type == 3){
	# calculating
	print("Begin calculating summary statistics...")
	myMean<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_mean_function(x)), mc.cores=my_cpu))
	myLogCVsq<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_log_coefvsqr_function(x)), mc.cores=my_cpu))
	keep_rows <- which(!is.na(myLogCVsq))
	myMean <- myMean[keep_rows]
	myLogCVsq<-myLogCVsq[keep_rows]
	model.data <- data.frame(x=myMean, y=myLogCVsq)
	gam.obj <- gam(y ~ s(x), data=model.data)
	mean.seq <- data.frame(x=seq(0,1,length=100))
	preds <- predict(gam.obj, newdata=mean.seq)
	predicted<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) predict.gam(gam.obj,newdata=data.frame(x=myMean[i]), type="response", se=TRUE)$fit, mc.cores=my_cpu))
	##model_se<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) predict.gam(gam.obj,newdata=data.frame(x=myMean[i]), type="response", se=TRUE)$se.fit, mc.cores=my_cpu))
	myResiduals <- residuals(gam.obj, type="response")
	keep_rows <- which(!is.na(myResiduals))
	# plotting
	pdf(paste0(outName, ".GAM.pdf"))
	plotHeatScatter(x=myMean, y=myLogCVsq, main="General additive model", ylab="Log(CV^2)", xlab="Mean methylation",numBins=100)
	plotHeatScatter(x=myMean[keep_rows], y=myResiduals[keep_rows], main="Residuals", ylab="Residuals", xlab="Predicted", numBins=100)
	dev.off()
	# write table
	df <- data.frame(chr=data$chr[keep_rows], pos=data$start[keep_rows], stat=exp(myResiduals[keep_rows]+predicted[keep_rows]))
	#df <- data.frame(chr=data$chr[keep_rows], pos=data$start[keep_rows], stat=myLogCVsq[keep_rows])
	write.table(df, file=paste0(outName,".coefv.txt"), quote=F,sep="\t",row=F,col=T)
	#get_fastseg(df)
	#get_segment(df)
	#get_hmm_seg(df)
}

if(analysis_type == 4){
        # calculating
        print("Begin calculating summary statistics...")
        myMean<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_mean_function(1-as.numeric(x))), mc.cores=my_cpu))
        myLogCVsq<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) apply(data[i,], 1, function(x) my_log_coefvsqr_function(1-as.numeric(x))), mc.cores=my_cpu))
        keep_rows <- which(!is.na(myLogCVsq))
        myMean <- myMean[keep_rows]
        myLogCVsq<-myLogCVsq[keep_rows]
        model.data <- data.frame(x=myMean, y=myLogCVsq)
        gam.obj <- gam(y ~ s(x), data=model.data)
        mean.seq <- data.frame(x=seq(0,1,length=100))
        preds <- predict(gam.obj, newdata=mean.seq)
        predicted<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) predict.gam(gam.obj,newdata=data.frame(x=myMean[i]), type="response", se=TRUE)$fit, mc.cores=my_cpu))
        ##model_se<- unlist(mclapply(splitIndices(nrow(data), my_cpu), function(i) predict.gam(gam.obj,newdata=data.frame(x=myMean[i]), type="response", se=TRUE)$se.fit, mc.cores=my_cpu))
        myResiduals <- residuals(gam.obj, type="deviance")
        keep_rows <- which(!is.na(myResiduals))
        # plotting
        pdf(paste0(outName, ".GAM.pdf"))
        plotHeatScatter(x=myMean, y=myLogCVsq, main="General additive model", ylab="Log(CV^2)", xlab="Mean methylation",numBins=100)
        plotHeatScatter(x=myMean[keep_rows], y=myResiduals[keep_rows], main="Residuals", ylab="Residuals", xlab="Predicted", numBins=100)
        dev.off()
        # write table
        df <- data.frame(chr=data$chr[keep_rows], pos=data$start[keep_rows], stat=exp(myResiduals[keep_rows]+predicted[keep_rows]))
        #df <- data.frame(chr=data$chr[keep_rows], pos=data$start[keep_rows], stat=myLogCVsq[keep_rows])
        write.table(df, file=paste0(outName,".coefv.txt"), quote=F,sep="\t",row=F,col=T)
        #get_fastseg(df)
        #get_segment(df)
        #get_hmm_seg(df)
}

